"""Run short velocity sensitivity pilots after the active serial control worker."""
import json,subprocess,time,fcntl
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];RUN=ROOT/'runs'
def main():
    lock=(RUN/'velocity-queue.lock').open('w')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print('Velocity queue already active');return
    status=RUN/'velocity-queue-status.json'
    while True:
        f=RUN/'ten-stage-execution.json';state=json.loads(f.read_text()) if f.exists() else {}
        if (RUN/'ten-stage-execution-error.json').exists():raise RuntimeError('Main queue reports an error; inspect it before further pilots')
        if state.get('state')=='bounded_jobs_completed_further_research_gates_pending':break
        status.write_text(json.dumps({'state':'waiting_for_serial_controls','scope':'short velocity sensitivity, not production'},indent=2)+'\n');time.sleep(15)
    status.write_text(json.dumps({'state':'running_velocity_pilots'},indent=2)+'\n')
    subprocess.run(['python','scripts/run_velocity_pilots.py'],cwd=ROOT,check=True)
    subprocess.run(['python','scripts/report_ten_stages.py'],cwd=ROOT,check=True)
    status.write_text(json.dumps({'state':'velocity_pilots_completed_not_production'},indent=2)+'\n')
if __name__=='__main__':
    try:main()
    except Exception as e:
        (RUN/'velocity-queue-status.json').write_text(json.dumps({'state':'blocked','reason':str(e)},indent=2)+'\n');raise
