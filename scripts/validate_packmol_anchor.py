"""Executable probes of PBC distance and rigid-body boundary crossing.

The artificial two-point rod is a geometry regression, not a chemical molecule.
"""
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from pc_length_parameters import ROOT,sha
from run_pc_sqm_checkpointed import write_json


def probe(out,binary,points,poses):
    out.mkdir()
    lines=[f'HETATM{i+1:5d}  C   MOL     1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  ' for i,(x,y,z) in enumerate(points)]
    (out/'molecule.pdb').write_text('\n'.join(lines)+'\nEND\n')
    (out/'source.restart').write_text('\n'.join(' '.join(map(str,row)) for row in poses)+'\n')
    (out/'packmol.inp').write_text('tolerance 2.004\nprecision 0.01\nfiletype pdb\noutput packed.pdb\nseed 2026\n'
        'pbc 0 0 0 12 12 12\nnloop 10\nmaxit 200\nrestart_from source.restart\nrestart_to packing.restart\n'
        f'structure molecule.pdb\nnumber {len(poses)}\nend structure\n')
    with (out/'packmol.inp').open() as inp,(out/'packmol.log').open('x') as log:
        result=subprocess.run([str(binary)],cwd=out,stdin=inp,stdout=log,stderr=subprocess.STDOUT,timeout=30)
    text=(out/'packmol.log').read_text();xyz=None
    if (out/'packed.pdb').exists():
        rows=[r for r in (out/'packed.pdb').read_text().splitlines() if r.startswith(('ATOM  ','HETATM'))]
        xyz=np.array([[float(r[j:j+8]) for j in (30,38,46)] for r in rows])
    return result.returncode,'Success!' in text,xyz


def main():
    build=json.loads((ROOT/'models/surface-revision-v2/packmol-anchor-build.json').read_text())
    binary=Path(build['binary']);original=Path('/workspace/md-tools/root/usr/bin/packmol')
    if sha(binary)!=build['binary_sha256']:raise ValueError('Variant binary changed')
    out=Path(tempfile.mkdtemp(prefix='packmol-anchor-validation-',dir=ROOT/'runs'))
    rod=[[0.,0.,0.],[30.,0.,0.]];pose=[[18.,6.,6.,0.,0.,0.]]
    old_code,old_success,_=probe(out/'original-rod',original,rod,pose)
    reference=Path(build['unmodified_reference_binary'])
    if sha(reference)!=build['unmodified_reference_sha256']:raise ValueError('Reference binary changed')
    ref_code,ref_success,_=probe(out/'same-version-reference-rod',reference,rod,pose)
    code,success,xyz=probe(out/'variant-rod',binary,rod,pose)
    if old_success or ref_success or not success or code!=0:raise ValueError('Boundary-crossing executable comparison failed')
    if not np.any((xyz<0)|(xyz>12)) or abs(np.linalg.norm(xyz[1]-xyz[0])-30)>.002:
        raise ValueError('Rod did not cross the primary cell with rigid geometry preserved')
    code2,success2,xy2=probe(out/'periodic-pair',binary,[[0.,0.,0.]],[[.5,6,6,0,0,0],[11.5,6,6,0,0,0]])
    delta=xy2[1]-xy2[0];delta-=12*np.rint(delta/12);distance=float(np.linalg.norm(delta))
    if code2 or not success2 or distance<1.998:raise ValueError('Periodic close pair was not resolved')
    report={'state':'variant executable geometry probes passed','original_rod_exit':old_code,'original_rod_success':old_success,
            'unmodified_same_version_rod_exit':ref_code,'unmodified_same_version_rod_success':ref_success,
            'source_version':build['source_version'],'reference_binary_sha256':sha(reference),
            'variant_rod_exit':code,'variant_rod_distance_A':float(np.linalg.norm(xyz[1]-xyz[0])),
            'variant_rod_crosses_primary_cell':True,'initial_periodic_pair_distance_A':1.,'final_periodic_pair_distance_A':distance,
            'binary_sha256':sha(binary),'original_binary_unchanged':sha(original)==build['original_binary_sha256'],
            'probe_directory':str(out.relative_to(ROOT)),'scope':'geometry implementation tests;no physical substrate validation',
            'remaining':'every real PC output must independently pass periodic interchain,self-image and bond checks'}
    write_json(ROOT/'models/surface-revision-v2/packmol-anchor-validation.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':main()
