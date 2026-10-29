"""Build a separate Packmol variant with one PBC box anchor per molecule.

The installed Packmol21.0.2 and inspected source21.2.3 add an automatic inside-box
constraint to every atom even when the input selects atom1 only. This variant
is built from source21.2.3 and changes that automatic assignment to
iatom==1; periodic interatomic distances are unchanged and independently audited.
The original source and installed executable are never modified.
"""
import json
from pathlib import Path
import shutil
import subprocess
from pc_length_parameters import ROOT,sha
from run_pc_sqm_checkpointed import write_json

TOOLS=Path('/workspace/md-tools')
SOURCE=TOOLS/'packmol-source'
DEST=TOOLS/'packmol-anchor-pbc-v1'


def build():
    commit=subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip()
    if commit!='14e50c65fa9120b58e1ba33ad482c7c7260f72b2':
        raise ValueError('Packmol source revision differs from reviewed21.2.3 source')
    original=SOURCE/'app/packmol.f90';text=original.read_text()
    old='        if (using_pbc) then\n          nratom(icart) = nratom(icart) + 1'
    new='        if (using_pbc .and. iatom == 1) then\n          nratom(icart) = nratom(icart) + 1'
    if text.count(old)!=1:raise ValueError('Source differs from reviewed automatic restriction block')
    header=' PBC on: Constraint automatically added to non-fixed atoms:'
    if text.count(header)!=1:raise ValueError('Unexpected source diagnostic')
    changed=text.replace(old,new).replace(header,' PBC anchor variant: automatic box constraint applies to atom1 of each molecule:')
    if DEST.exists():
        if (DEST/'packmol').exists() or (DEST/'app/packmol.f90').read_text()!=changed:
            raise RuntimeError('Existing binary or differing source retained; no overwrite')
    else:
        shutil.copytree(SOURCE,DEST,ignore=shutil.ignore_patterns('.git','*.o','*.mod','packmol'))
        (DEST/'app/packmol.f90').write_text(changed)
    compiler=TOOLS/'root/usr/bin/x86_64-linux-gnu-gfortran-14'
    # The compiler is locally extracted; use the matching system GCC runtime
    # search directory for crt objects and its LTO linker plugin.
    command=['make','-j2',f'FORTRAN={compiler}','FLAGS=-O2 -B/usr/lib/gcc/x86_64-linux-gnu/14/']
    attempt=len(list(DEST.glob('build*.log')))+1
    with (DEST/f'build-attempt-{attempt:02d}.log').open('x') as log:
        subprocess.run(command,cwd=DEST,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    reference=TOOLS/'packmol-reference-21.2.3'
    if not reference.exists():
        shutil.copytree(SOURCE,reference,ignore=shutil.ignore_patterns('.git','*.o','*.mod','packmol'))
        with (reference/'build.log').open('x') as log:
            subprocess.run(command,cwd=reference,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    if sha(reference/'app/packmol.f90')!=sha(original):raise ValueError('Reference source differs')
    report={'source_directory':str(SOURCE),'source_main_sha256':sha(original),
            'source_version':'21.2.3','source_git_commit':subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip(),
            'modified_main_sha256':sha(DEST/'app/packmol.f90'),'binary':str(DEST/'packmol'),
            'binary_sha256':sha(DEST/'packmol'),'original_binary_sha256':sha(TOOLS/'root/usr/bin/packmol'),
            'unmodified_reference_binary':str(reference/'packmol'),'unmodified_reference_sha256':sha(reference/'packmol'),
            'compiler':subprocess.check_output([str(compiler),'--version'],text=True).splitlines()[0],
            'build_command':command,'patch':'automatic PBC restraint assignment changes using_pbc to using_pbc .and. iatom==1;diagnostic label updated',
            'unchanged_relative_to_source21.2.3':'PBC cell/distance kernels,geometry,atom radii,target tolerance and optimizer',
            'physical_approval':False,'requires':'boundary-crossing and periodic-distance regression probes before PC use'}
    write_json(ROOT/'models/surface-revision-v2/packmol-anchor-build.json',report)
    print(json.dumps(report,indent=2))


if __name__=='__main__':build()
