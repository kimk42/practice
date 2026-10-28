"""Extend the accepted2ps azide control by100ps with its unchanged TPR/checkpoint."""
import hashlib
import fcntl
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'runs/azide-rigid-water-control'
OUT = ROOT/'runs/azide-rigid-water-102ps'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    lock=(ROOT/'runs/azide-solution-continuation.lock').open('w')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (OUT/'completed.json').exists():
        print((OUT/'completed.json').read_text()); return
    audit = json.loads((ROOT/'models/buffer-parameters/azide-rigid-water-control-audit.json').read_text())
    spec = {name: sha(SOURCE/name) for name in ['md.tpr','md.cpt']}
    assert spec['md.tpr'] == audit['input_and_output_sha256']['md.tpr']
    assert spec['md.cpt'] == audit['input_and_output_sha256']['md.cpt']
    if (OUT/'spec.json').exists() and json.loads((OUT/'spec.json').read_text()) != spec:
        raise RuntimeError('Source changed;existing continuation preserved')
    OUT.mkdir(exist_ok=True)
    (OUT/'spec.json').write_text(json.dumps(spec,indent=2)+'\n')
    def run(args, log):
        with (OUT/log).open('a') as stream:
            subprocess.run(['gmx',*args],cwd=OUT,stdout=stream,stderr=subprocess.STDOUT,check=True)
    if not (OUT/'nvt.tpr').exists():
        shutil.copy2(SOURCE/'md.cpt',OUT/'source.cpt')
        run(['convert-tpr','-s',str(SOURCE/'md.tpr'),'-extend','100','-o','nvt.tpr'],'convert.log')
    args=['mdrun','-s','nvt.tpr','-deffnm','nvt','-ntmpi','1','-ntomp','2','-cpt','1','-maxh','1']
    if (OUT/'nvt.cpt').exists():
        args += ['-cpi','nvt.cpt','-append']
    else:
        if list(OUT.glob('nvt.part*.log')):
            raise RuntimeError('Partial run without a checkpoint requires inspection')
        args += ['-cpi','source.cpt','-noappend']
    run(args,'launch.log')
    assert all(sha(SOURCE/name)==digest for name,digest in spec.items())
    energy_files=list(OUT.glob('nvt.part*.edr'))
    if len(energy_files)!=1:raise RuntimeError('Inspect continuation energy parts before accepting completion')
    with (OUT/'duration-check.log').open('a') as stream:
        subprocess.run(['gmx','energy','-f',str(energy_files[0]),'-o','duration-check.xvg'],
                       cwd=OUT,input='Potential\n0\n',text=True,stdout=stream,stderr=subprocess.STDOUT,check=True)
    times=[float(line.split()[0]) for line in (OUT/'duration-check.xvg').read_text().splitlines()
           if line.strip() and not line.startswith(('#','@'))]
    if not times or abs(times[-1]-102)>0.001:
        raise RuntimeError('Continuation not finished;checkpoint preserved,resume this same directory')
    (OUT/'completed.json').write_text(json.dumps({
        'state':'mdrun returned;duration/energy/geometry and solution diagnostics require analysis',
        'source_duration_ps':2,'requested_total_ps':102,
        'method':'convert-tpr extension plus original checkpoint;no velocity or model reinitialization',
        'physical_force_field_approved':False},indent=2)+'\n')
    print(OUT)


if __name__ == '__main__': main()
