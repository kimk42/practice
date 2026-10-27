"""Checkpointed10ps baseline follow-ups after completed pilot/control QA.
No full-buffer or adsorption/equilibrium claim is made.
"""
import json,subprocess,fcntl,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];RUN=ROOT/'runs'
def main():
    lock=(RUN/'paired-interface-followup.lock').open('w')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print('Paired follow-up already active');return
    qa=json.loads((ROOT/'models/analysis/fxii-water-100ps-QA.json').read_text())
    if qa['disulfide_geometry_flags']:raise RuntimeError('FXII QA review required')
    index=json.loads((ROOT/'models/validation/velocity-pilots-audit.json').read_text())
    if len(index['records'])!=18:raise RuntimeError('Pilot gate incomplete')
    for name in ['coating-control-100ps','fxii-water-control-100ps']:
        if json.loads((RUN/name/'result.json').read_text())['duration_ps']!=100:raise RuntimeError('Control incomplete')
    if subprocess.run(['pgrep','-x','gmx'],capture_output=True).returncode==0:raise RuntimeError('Another GROMACS job is active;preserve it and wait')
    completed=[];status=RUN/'paired-interface-followup-status.json'
    for source,target in [('fxii-coating-near-pilot','fxii-coating-near-10ps'),('fxii-bare-near-pilot','fxii-bare-near-10ps')]:
        status.write_text(json.dumps({'state':'running','pid':os.getpid(),'job':target,'completed':completed,'duration_ps':10,'conditions':'legacy300K,neutralizing ions only,unglycosylated predicted FXII;not cPBSzI'},indent=2)+'\n')
        with (RUN/(target+'-workflow.log')).open('w') as log:
            subprocess.run(['python','scripts/run_control_followup.py','--source-dir',str(RUN/source),'--basename','md','--output-dir',str(RUN/target),'--duration-ps','10'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
        result=json.loads((RUN/target/'result.json').read_text());completed.append({'case':target,'result':result})
    status.write_text(json.dumps({'state':'paired10ps numerical follow-ups completed','completed':completed,'scope':'baseline initial relaxation only;not matched-buffer validation or adsorption comparison'},indent=2)+'\n')
if __name__=='__main__':
    try:main()
    except Exception as error:
        (RUN/'paired-interface-followup-error.json').write_text(json.dumps({'state':'blocked','reason':str(error)},indent=2)+'\n');raise
