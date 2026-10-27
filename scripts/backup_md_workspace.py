"""Create a verified local recovery archive;this is not off-environment backup."""
import argparse,datetime,hashlib,json,subprocess,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--destination',type=Path,default=Path('/workspace/md-backups'));a=p.parse_args()
    for name in ['gmx','sqm','antechamber']:
        if subprocess.run(['pgrep','-x',name],capture_output=True).returncode==0:raise RuntimeError('Active calculation detected;preserve it and create backup at a consistent checkpoint later')
    a.destination=a.destination.resolve()
    if a.destination==ROOT or ROOT in a.destination.parents:raise ValueError('Backup must be outside the source worktree')
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');out=a.destination/('cbaa-surface-'+stamp);out.mkdir(parents=True)
    print(str(out),flush=True)
    (out/'status.json').write_text(json.dumps({'state':'creating_local_archive','external_backup_complete':False},indent=2)+'\n')
    snapshots={str(f.relative_to(ROOT)):(f.stat().st_size,f.stat().st_mtime_ns) for f in (ROOT/'runs').rglob('*.cpt')}
    (out/'worktree.patch').write_bytes(git('diff','--binary','HEAD'))
    (out/'worktree-status.txt').write_bytes(git('status','--short'))
    bundle=out/'repository.bundle';subprocess.run(['git','bundle','create',str(bundle),'md/cbaa-surface'],cwd=ROOT,check=True)
    with (out/'bundle-verify.log').open('w') as log:subprocess.run(['git','bundle','verify',str(bundle)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    archive=out/'project-and-runs.tar.zst'
    with (out/'archive.log').open('w') as log:subprocess.run(['tar','--use-compress-program=zstd -T2 -3','-cf',str(archive),'--exclude=./.git','--exclude=*/__pycache__','--exclude=*.pyc','-C',str(ROOT),'.'],stdout=log,stderr=subprocess.STDOUT,check=True)
    with (out/'archive-verify.log').open('w') as log:subprocess.run(['zstd','-t',str(archive)],stdout=log,stderr=subprocess.STDOUT,check=True)
    with (out/'archive-contents.txt').open('w') as listing:subprocess.run(['tar','--zstd','-tf',str(archive)],stdout=listing,check=True)
    now={str(f.relative_to(ROOT)):(f.stat().st_size,f.stat().st_mtime_ns) for f in (ROOT/'runs').rglob('*.cpt')}
    if now!=snapshots:raise RuntimeError('Checkpoint files changed during backup;do not mark this archive consistent')
    checkpoint_files=[{'path':name,'size_bytes':info[0],'sha256':sha(ROOT/name)} for name,info in sorted(snapshots.items())]
    files=[{'file':f.name,'size_bytes':f.stat().st_size,'sha256':sha(f)} for f in [archive,bundle,out/'worktree.patch']]
    report={'state':'verified_local_recovery_backup','source_worktree':str(ROOT),'branch':git('branch','--show-current').decode().strip(),'head':git('rev-parse','HEAD').decode().strip(),'files':files,'checkpoints':checkpoint_files,'checkpoint_snapshot_unchanged':True,'contains_private_raw_documents':True,'external_backup_complete':False,'limitation':'same cloud environment;does not protect against deletion/loss of the environment;do not upload this archive to the public repository'}
    (out/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');(out/'status.json').write_text(json.dumps({'state':'verified_local_recovery_backup','external_backup_complete':False},indent=2)+'\n')
    (out/'RESTORE.txt').write_text('Private recovery archive. Do not publish.\n1. Verify SHA256 values in manifest.json.\n2. git clone --branch md/cbaa-surface repository.bundle restored-md\n3. In restored-md,apply worktree.patch with git apply using its absolute path.\n4. Extract project-and-runs.tar.zst into restored-md (tar --zstd -xf ARCHIVE -C restored-md).\n5. Reinstall shared tools with the project setup;the installed /workspace/md-tools is not bundled.\n6. Inspect job/checkpoint status before resuming;archiving does not start any calculation.\n')
    (a.destination/'latest-backup.json').write_text(json.dumps({'directory':str(out),'head':report['head'],'state':report['state'],'external_backup_complete':False},indent=2)+'\n')
    print(json.dumps({'state':report['state'],'directory':str(out),'archive_bytes':archive.stat().st_size,'checkpoints':len(checkpoint_files),'external_backup_complete':False}),flush=True)
if __name__=='__main__':main()
