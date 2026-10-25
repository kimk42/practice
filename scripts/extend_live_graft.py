"""Preserve live SQM while replacing ACPYPE's non-adjustable alarm watchdog.

Pause only ACPYPE (not its antechamber/SQM children). Once the existing
antechamber finishes, reuse its completed BCC mol2 with a new ACPYPE process.
The old pending alarm must never be resumed after successful handoff.
"""
import hashlib,json,os,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];RUN=ROOT/'runs'
def alive(pid):
    try:
        state=Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0]
        return state!='Z'
    except FileNotFoundError:return False

def main():
    raise RuntimeError('Deprecated: pausing the wrapper may block its stdout pipe; recover completed SQM with finalize_completed_sqm.py instead')
    parent=21924;child=21996;sqm=22018;controller=21882
    active=json.loads((RUN/'workflow-active-parameters.json').read_text())
    directory=ROOT/active['GRAFT'];work=directory/'GRAFT.acpype'
    status=RUN/'live-extension-status.json'
    def save(state,**extra):
        status.write_text(json.dumps({'state':state,'parent_pid':parent,'sqm_pid':sqm,'directory':active['GRAFT'],**extra},indent=2)+'\n')
    if status.exists():raise RuntimeError('Extension already recorded; inspect before any repeat')
    # Verify process identity and input before signaling any PID.
    for pid,expected in [(parent,'acpype'),(child,'antechamber'),(sqm,'sqm')]:
        if Path(f'/proc/{pid}/comm').read_text().strip()!=expected:raise RuntimeError('PID identity mismatch')
    if not alive(sqm):raise RuntimeError('Original SQM is no longer alive')
    elapsed=int(subprocess.check_output(['ps','-p',str(parent),'-o','etimes='],text=True).strip())
    remaining=86400-elapsed
    if remaining<=0:raise RuntimeError('Already beyond requested 24-hour total limit')
    deadline=time.monotonic()+remaining
    os.kill(parent,signal.SIGSTOP)
    # Original controller has a separate 13-hour deadline; replace it too.
    os.kill(controller,signal.SIGTERM)
    save('live_SQM_preserved_alarm_parent_paused',total_limit_hours=24,seconds_remaining_at_handoff=remaining)
    success=False
    try:
        while alive(child):
            if time.monotonic()>deadline:raise TimeoutError('Extended 24-hour total limit reached')
            time.sleep(15)
        mol2=work/'GRAFT_bcc_gaff2.mol2'
        output=(work/'sqm.out').read_text()
        if 'Calculation Completed' not in output:raise RuntimeError('SQM did not complete successfully')
        if not mol2.exists() or '@<TRIPOS>ATOM' not in mol2.read_text():raise RuntimeError('Completed BCC mol2 is absent')
        digest=hashlib.sha256(mol2.read_bytes()).hexdigest()
        # Copy the complete intermediate before topology conversion/cleanup.
        (directory/'completed-bcc-before-handoff.mol2').write_bytes(mol2.read_bytes())
        os.kill(parent,signal.SIGKILL)
        save('completed_BCC_reused_for_topology',bcc_sha256=digest)
        # Use the same basename/input/method and existing directory: ACPYPE's
        # documented file-presence branch skips Antechamber for completed mol2.
        command=['acpype','-i',str(ROOT/'models/pc-graft/cop10-graft-bpa-pc4-one-site.mol'),'-b','GRAFT','-a','gaff2','-c','bcc','-n','0','-m','1','-s','86400']
        with (directory/'GRAFT-handoff.log').open('w') as log:
            subprocess.run(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,check=True)
        if hashlib.sha256(mol2.read_bytes()).hexdigest()!=digest:raise RuntimeError('Completed BCC charges changed during handoff')
        save('topology_handoff_completed',bcc_sha256=digest)
        success=True
    finally:
        if not success and alive(parent):
            # Restore the original watchdog; its pending alarm terminates the
            # calculation on an extension failure rather than leaving it parked.
            os.kill(parent,signal.SIGCONT)
    env=os.environ.copy();env['MD_WORKFLOW_WAIT_SECONDS']='14400'
    with (RUN/'extended-surface-controller.log').open('w') as log:
        subprocess.run(['python','scripts/continue_three_stages.py'],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    save('surface_controller_finished')
if __name__=='__main__':
    try:main()
    except Exception as error:
        (RUN/'live-extension-error.json').write_text(json.dumps({'error':str(error)},indent=2)+'\n')
        raise
