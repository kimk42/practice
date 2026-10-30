"""Run a bounded postprocessing process tree without leaving wrapper children."""
import os
import signal
import subprocess


def run_group(command, *, cwd, env, log, timeout, pass_fds=()):
    process=subprocess.Popen(command,cwd=cwd,env=env,stdin=subprocess.DEVNULL,
                             stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
                             pass_fds=pass_fds)
    try:
        code=process.wait(timeout=timeout)
    except BaseException:
        try:
            os.killpg(process.pid,signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=2)
        except subprocess.TimeoutExpired:
            pass
        # A shell may exit before its child, including a child ignoring TERM.
        try:
            os.killpg(process.pid,signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()
        raise
    if code:
        raise subprocess.CalledProcessError(code,command)
    return code
