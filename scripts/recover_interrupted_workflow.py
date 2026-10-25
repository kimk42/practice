"""Recover interrupted quantum jobs in fresh directories and launch gated controls.

No SQM optimization restart coordinates were saved in the interrupted runs, so
this repeats AM1-BCC from the same input and records that loss explicitly.
"""
import json,os,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main():
    # Refuse concurrent recovery rather than duplicating expensive SQM jobs.
    result=subprocess.run(['pgrep','-x','sqm'],capture_output=True,text=True)
    if result.returncode==0:raise RuntimeError('SQM already running; inspect before recovery')
    manifest=ROOT/'runs/workflow-active-parameters.json'
    active=json.loads(manifest.read_text());previous=dict(active);jobs=[]
    sources={'COP10':'models/copolymer/cbaa1-bma-bpaa-dp10-seed2026.mol','GRAFT':'models/pc-graft/cop10-graft-bpa-pc4-one-site.mol'}
    for name,source in sources.items():
        if (ROOT/active[name]/(name+'.acpype')/(name+'_GMX.itp')).exists():continue
        path=Path(tempfile.mkdtemp(prefix=name.lower()+'-recovery-',dir=ROOT/'runs'))
        active[name]=str(path.relative_to(ROOT))
        env=os.environ.copy();env['MD_SQM_MAX_SECONDS']=os.environ.get('MD_SQM_MAX_SECONDS','86400')
        with (path/'launch.log').open('w') as log:
            proc=subprocess.Popen(['bash','scripts/parameterize_candidate.sh',source,name,str(path)],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        jobs.append({'name':name,'pid':proc.pid,'directory':active[name]})
    manifest.write_text(json.dumps(active,indent=2)+'\n')
    (ROOT/'runs/recovery-manifest.json').write_text(json.dumps({'previous':previous,'active':active,'jobs':jobs,'quantum_seconds_limit':int(os.environ.get('MD_SQM_MAX_SECONDS','86400')),'restart':'same chemical inputs/method, optimization starts afresh because no usable saved optimizer coordinates exist'},indent=2)+'\n')
    with (ROOT/'runs/fxii-recovery-launch.log').open('w') as log:
        subprocess.Popen(['python','scripts/recover_fxii_control.py'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    # Give the quantum jobs their full configured limit, then perform gated surface MD.
    env=os.environ.copy();env['MD_WORKFLOW_WAIT_SECONDS']=str(int(env.get('MD_SQM_MAX_SECONDS','86400'))+3600)
    subprocess.run(['python','scripts/continue_three_stages.py'],cwd=ROOT,env=env,check=True)
if __name__=='__main__':main()
