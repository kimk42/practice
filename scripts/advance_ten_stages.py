"""Sequential CPU workflow for authorized pilots/controls; never label gates complete prematurely."""
import json,subprocess,time,os,fcntl
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];RUN=ROOT/'runs'
def active_recovery():return subprocess.run(['pgrep','-f','^python scripts/recover_expanded_em.py$'],capture_output=True).returncode==0
def main():
    lock=(RUN/'ten-stage-execution.lock').open('w')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print('Ten-stage worker already active');return
    error=RUN/'ten-stage-execution-error.json'
    if error.exists():error.rename(RUN/f'ten-stage-previous-error-{time.time_ns()}.json')
    status=RUN/'ten-stage-execution.json';completed=[]
    def save(state,**extra):status.write_text(json.dumps({'state':state,'pid':os.getpid(),'completed_jobs':completed,**extra},indent=2)+'\n')
    while active_recovery() and not (RUN/'surface-expanded-18x12/result.json').exists():
        save('waiting_for_existing_expanded_EM');time.sleep(15)
    jobs=[
      ('expanded_surface_EM',['python','scripts/recover_expanded_em.py'],'surface-expanded-18x12/result.json'),
      ('coated_FXII_pose2',['python','scripts/run_fxii_surface_pilot.py'],'fxii-surface-pilot/result.json'),
      ('matched_bare_PC_EM',['python','scripts/minimize_bare_pc.py'],'bare-pc-matched/result.json'),
      ('bare_PC_FXII_pose2',['python','scripts/run_fxii_surface_pilot.py','--surface-dir',str(RUN/'bare-pc-matched'),'--output-dir',str(RUN/'fxii-bare-pc-pilot')],'fxii-bare-pc-pilot/result.json'),
      ('coated_FXII_pose3',['python','scripts/run_fxii_surface_pilot.py','--output-dir',str(RUN/'fxii-coating-pose3-pilot'),'--pose-file','pose-3.gro'],'fxii-coating-pose3-pilot/result.json'),
      ('bare_PC_FXII_pose3',['python','scripts/run_fxii_surface_pilot.py','--surface-dir',str(RUN/'bare-pc-matched'),'--output-dir',str(RUN/'fxii-bare-pc-pose3-pilot'),'--pose-file','pose-3.gro'],'fxii-bare-pc-pose3-pilot/result.json'),
      ('velocity_seed_pilots',['python','scripts/run_velocity_pilots.py'],'velocity-pilots/result.json'),
      ('coating_100ps_control',['python','scripts/run_control_followup.py','--source-dir',str(RUN/'surface-nvt-10ps'),'--basename','nvt','--output-dir',str(RUN/'coating-control-100ps')],'coating-control-100ps/result.json'),
      ('FXII_water_100ps_control',['python','scripts/run_control_followup.py','--source-dir',str(RUN/'fxii-unglycosylated-topology'),'--basename','md-recovery','--output-dir',str(RUN/'fxii-water-control-100ps')],'fxii-water-control-100ps/result.json'),
    ]
    for name,command,result in jobs:
        if not (RUN/result).exists():
            save('running',job=name,result=result)
            with (RUN/(name+'-workflow.log')).open('w') as log:subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
        if not (RUN/result).exists():raise RuntimeError(name+': no verified result')
        completed.append(name);save('job_completed',job=name)
        subprocess.run(['python','scripts/report_ten_stages.py'],cwd=ROOT,check=True)
    save('bounded_jobs_completed_further_research_gates_pending',note='Short controls do not establish equilibration/production; glycoform, fibrinogen and sampling/resource gates remain explicit')
    subprocess.run(['python','scripts/report_ten_stages.py'],cwd=ROOT,check=True)
if __name__=='__main__':
    try:main()
    except Exception as error:
        (RUN/'ten-stage-execution-error.json').write_text(json.dumps({'state':'blocked','reason':str(error)},indent=2)+'\n')
        subprocess.run(['python','scripts/report_ten_stages.py'],cwd=ROOT)
        raise
