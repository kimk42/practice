"""Persist the authorized NPT→check→interface→large-solvation→report workflow.
This supervisor adopts the already running NPT job and stops on execution failures.
"""
import argparse,json,os,subprocess,sys,time,traceback
from pathlib import Path
try:
    from .extend_pe_npt import atomic_json
    from .report_pe_followup import publish
except ImportError:
    from extend_pe_npt import atomic_json
    from report_pe_followup import publish

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/pe-followup01';NPT=RUN/'01-npt-100ps';OUT=ROOT/'artifacts/pe-followup';LARGE=ROOT/'data/fibrinogen-core-pe-solvated01'

def state(status,**extra):
    payload={'status':status,'updated_unix':time.time(),'authorized_scope':['100 ps NPT extension','Structural diagnostics','10 ps supported-interface preparation and protein restraint release','Six-chain core/PE solvation and neutralization','Merge all results and data into portable report and archive'],'scientific_equilibrium_established':False,**extra}
    atomic_json(RUN/'pipeline.json',payload);print(status,flush=True);return payload


def command(script,args,log):
    with (RUN/log).open('w') as f:
        p=subprocess.Popen([sys.executable,str(ROOT/'md'/script),*map(str,args)],cwd=ROOT,env={**os.environ,'OPENMM_CPU_THREADS':'4','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},stdout=f,stderr=subprocess.STDOUT)
        atomic_json(RUN/'active-process.json',{'pid':p.pid,'script':script,'started_unix':time.time()})
        while p.poll() is None:time.sleep(10)
        if p.returncode:raise RuntimeError(f'{script} exited {p.returncode}; see {RUN/log}')


def main():
    RUN.mkdir(exist_ok=True);lock=RUN/'supervisor.lock'
    # Exclusive creation prevents two supervisors from launching the same stages.
    fd=os.open(lock,os.O_WRONLY|os.O_CREAT|os.O_EXCL);os.write(fd,str(os.getpid()).encode());os.close(fd)
    try:
        state('NPT_running');publish(RUN,OUT,LARGE);last_report_ps=-10
        while not (NPT/'result.json').exists():
            if (NPT/'failure.json').exists():raise RuntimeError((NPT/'failure.json').read_text())
            progress=json.loads((NPT/'progress.json').read_text())
            if time.time()-progress['updated_unix']>600:raise RuntimeError('No NPT checkpoint update for 10 minutes; inspect worker before restarting')
            if progress['extension_time_ps']>=last_report_ps+5:
                publish(RUN,OUT,LARGE);last_report_ps=progress['extension_time_ps']
            time.sleep(10)
        state('structure_check_running');publish(RUN,OUT,LARGE)
        command('analyze_pe_followup.py',[NPT,ROOT/'data/pe-protein-model04',RUN/'02-structure-check'],'structure.log')
        analysis=json.loads((RUN/'02-structure-check/analysis.json').read_text())
        if not analysis['allow_next_preparation']:raise RuntimeError('Structure checks require review: '+str(analysis['failures']))
        state('interface_running');publish(RUN,OUT,LARGE)
        command('relax_pe_interface.py',[NPT,ROOT/'data/pe-protein-model04',RUN/'02-structure-check',RUN/'03-interface-release'],'interface.log')
        state('solvation_running');publish(RUN,OUT,LARGE)
        command('solvate_core_pe.py',[ROOT/'data/fibrinogen-core-pe02',LARGE],'solvation.log')
        state('completed',note='All requested execution/preparation stages completed; equilibrium and production readiness not established.');publish(RUN,OUT,LARGE)
        command('collect_results.py',[],'collection.log')
        # Validate every archived file against the embedded manifest, including trajectories.
        import zipfile,hashlib
        archive=ROOT/'artifacts/collected-results/results-complete.zip'
        with zipfile.ZipFile(archive) as z:
            manifest=json.loads(z.read('FILES.json'))
            for entry in manifest['files']:
                data=z.read(entry['path'])
                if len(data)!=entry['bytes'] or hashlib.sha256(data).hexdigest()!=entry['sha256']:raise ValueError('Archive integrity failure: '+entry['path'])
        atomic_json(RUN/'archive-validation.json',{'status':'all_archive_sha256_checks_passed','files':len(manifest['files']),'bytes':archive.stat().st_size,'validated_unix':time.time()})
        print('PIPELINE AND COMPLETE REPORT FINISHED',flush=True)
    except BaseException as e:
        state('failed',error=f'{type(e).__name__}: {e}',traceback=traceback.format_exc())
        try:publish(RUN,OUT,LARGE)
        except Exception:traceback.print_exc()
        raise

if __name__=='__main__':main()
