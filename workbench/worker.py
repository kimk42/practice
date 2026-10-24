"""Independent process: survives browser/server closure, honors checkpointed controls."""
import argparse
import json
from pathlib import Path
import traceback
from .manager import atomic
from .config import ROOT


def main(job):
 from .engine import Engine,Cancelled,ReviewRequired
 from .report import verify,publish
 job=Path(job);engine=Engine(job);status='failed'
 try:
  result=engine.run();verify(engine);status=result['status'];engine.state(status,message='Executed preparation evaluated; physical equilibrium is not established',result=result,last=engine.rows[-1])
 except Cancelled as e:
  status='cancelled';engine.state(status,message=str(e),last=engine.rows[-1] if engine.rows else None)
 except ReviewRequired as e:
  status='needs_review';engine.state(status,message=str(e),last=engine.rows[-1] if engine.rows else None)
 except BaseException as e:
  engine.state('failed',message=f'{type(e).__name__}: {e}',traceback=traceback.format_exc(),last=engine.rows[-1] if engine.rows else None);traceback.print_exc()
 finally:
  try:publish(job)
  except Exception:traceback.print_exc()
  lock=job.parent/'.active'
  if lock.exists():
   try:
    if json.loads(lock.read_text())['job_id']==job.name:lock.unlink()
   except (OSError,ValueError):pass
 if status=='failed':raise SystemExit(1)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--job',required=True);a=p.parse_args();main(a.job)
