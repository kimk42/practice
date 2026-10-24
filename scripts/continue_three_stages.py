"""Continue existing calculations, retry time limits once, then run gated MD.

Timeout retries preserve the method/thresholds and use fresh directories. Other
errors stop the workflow. This process does not claim unrun stages succeeded.
"""
import json
import os
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs'
STATE=RUN/'three-stage-status.json'
ACTIVE=RUN/'workflow-active-parameters.json'
SOURCES={
 'COP10':'models/copolymer/cbaa1-bma-bpaa-dp10-seed2026.mol',
 'PC4':'models/pc-graft/bpa-pc-dp4-phenolic-ends.mol',
 'GRAFT':'models/pc-graft/cop10-graft-bpa-pc4-one-site.mol',
}
SUB={'COP10':'cop10','PC4':'pc4','GRAFT':'graft'}

def save(path,data):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data,indent=2)+'\n');temporary.replace(path)

def quantum_progress(folder,name):
    file=folder/(name+'.acpype')/'sqm.out'
    lines=file.read_text().splitlines() if file.exists() else []
    steps=[x for x in lines if x.startswith('xmin')]
    return steps[-1].strip() if steps else 'No reported optimization iteration yet; not converged'

def main():
    RUN.mkdir(exist_ok=True)
    active={n:f'runs/parameterization-optimized/{SUB[n]}' for n in SOURCES}
    if ACTIVE.exists():active=json.loads(ACTIVE.read_text())
    save(ACTIVE,active)
    attempts={n:int('parameterization-long' in active[n]) for n in SOURCES}
    processes={}
    packing_process=None
    packing_attempt=int((RUN/'pc-packing-long').exists())
    deadline=time.monotonic()+14400
    while time.monotonic()<deadline:
        summary={}
        for name in SOURCES:
            directory=ROOT/active[name]
            log=directory/(name+'.log')
            text=log.read_text() if log.exists() else ''
            ready=(directory/(name+'.acpype')/(name+'_GMX.itp')).is_file()
            if not ready and 'Timed out!' in text and not attempts[name]:
                destination=f'runs/parameterization-long/{SUB[name]}'
                if (ROOT/destination/(name+'.acpype')).exists():
                    raise RuntimeError('Refusing to overwrite a pre-existing long attempt: '+name)
                environment=os.environ.copy();environment['MD_SQM_MAX_SECONDS']='10800'
                with (RUN/(name+'-long-launch.log')).open('w') as output:
                    process=subprocess.Popen(['bash','scripts/parameterize_candidate.sh',SOURCES[name],name,destination],
                                             cwd=ROOT,env=environment,stdout=output,stderr=subprocess.STDOUT)
                processes[name]=process;attempts[name]=1;active[name]=destination;save(ACTIVE,active)
                directory=ROOT/destination;text='';ready=False
            if not ready and ('ERROR:' in text or 'ACPYPE FAILED' in text or 'Timed out!' in text):
                raise RuntimeError(name+': parameterization failed; inspect '+active[name]+'/'+name+'.log')
            if name in processes and processes[name].poll() not in (None,0) and not ready:
                raise RuntimeError(name+': long attempt exited with '+str(processes[name].returncode))
            summary[name]={'topology_file_present':ready,'attempt':attempts[name]+1,
                           'path':active[name],'quantum_progress':quantum_progress(directory,name)}
        packing_ready=(ROOT/'models/pc-patch/patch.json').is_file()
        packing_log=RUN/('pc-patch-long.log' if packing_attempt else 'pc-patch-resumed.log')
        packing_text=packing_log.read_text() if packing_log.exists() else ''
        if not packing_ready and 'Traceback' in packing_text:
            checkpoint=RUN/'pc-packing/packing.restart'
            if not packing_attempt and 'TimeoutExpired' in packing_text and checkpoint.exists():
                with (RUN/'pc-patch-long.log').open('w') as output:
                    packing_process=subprocess.Popen(['python','scripts/pack_pc_patch.py','--density','0.6',
                        '--resume',str(checkpoint),'--timeout','10800','--run-dir',str(RUN/'pc-packing-long')],
                        cwd=ROOT,stdout=output,stderr=subprocess.STDOUT)
                packing_attempt=1
            else:
                raise RuntimeError('Corrected packing did not complete; inspect '+str(packing_log.relative_to(ROOT)))
        if packing_process and packing_process.poll() not in (None,0) and not packing_ready:
            raise RuntimeError('Long packing attempt failed; inspect runs/pc-patch-long.log')
        current={'state':'stage1_quantum_parameters_pending','parameterization':summary,
                 'stage2_packing_complete':packing_ready,'packing_attempt':packing_attempt+1,'stage3_MD_complete':False,
                 'note':'File presence is not parameter approval; all audits run before MD.'}
        save(STATE,current)
        if all(x['topology_file_present'] for x in summary.values()) and packing_ready:
            from surface_validation import parameter_gate
            for name in SOURCES:
                if parameter_gate(name) is None:raise RuntimeError(name+': parameter audit is still pending')
            current['state']='running_stage1_audit_and_stages2_3';save(STATE,current)
            with (RUN/'three-stage-validation.log').open('w') as output:
                subprocess.run(['python','scripts/surface_validation.py','--run'],cwd=ROOT,stdout=output,
                               stderr=subprocess.STDOUT,check=True)
            result=json.loads((RUN/'surface-validation-status.json').read_text())
            if result.get('state')!='passed_initialization_check':
                raise RuntimeError('Surface validation did not pass: '+str(result))
            current['state']='three_stages_completed';current['stage3_MD_complete']=True
            current['result']=result;save(STATE,current);return
        time.sleep(15)
    raise RuntimeError('Four-hour workflow wait limit reached; inspect active job logs without overwriting them')

if __name__=='__main__':
    try:main()
    except Exception as error:
        save(STATE,{'state':'blocked','reason':str(error),'stage3_MD_complete':False})
        raise
