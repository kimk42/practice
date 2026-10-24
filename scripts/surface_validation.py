"""Gate and run exploratory surface validation after parameterization finishes.

Never fill missing force-field constants or override GROMACS warnings. Outputs
are initialization checks, not evidence of an equilibrated PC/coating interface.
"""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
import numpy as np
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'runs/parameterization-optimized'
STATUS=ROOT/'runs/surface-validation-status.json'
INPUTS={
 'COP10': ('cop10',ROOT/'models/copolymer/cbaa1-bma-bpaa-dp10-seed2026.mol',313),
 'PC4': ('pc4', ROOT/'models/pc-graft/bpa-pc-dp4-phenolic-ends.mol',132),
 'GRAFT': ('graft',ROOT/'models/pc-graft/cop10-graft-bpa-pc4-one-site.mol',445),
}

def section_rows(text, section):
    active=False
    for line in text.splitlines():
        value=line.split(';',1)[0].strip()
        if value.startswith('['):active=value.strip('[] ').lower()==section
        elif active and value and not value.startswith('#'):yield value.split()

def parameter_gate(name):
    sub,input_path,natoms=INPUTS[name]
    active_file=ROOT/'runs/workflow-active-parameters.json'
    active=json.loads(active_file.read_text()) if active_file.exists() else {}
    base=ROOT/active[name] if name in active else BASE/sub
    folder=base/(name+'.acpype')
    itp=folder/(name+'_GMX.itp');gro=folder/(name+'_GMX.gro')
    if not itp.exists() or not gro.exists():
        log=base/(name+'.log')
        if log.exists() and re.search(r'^.*(?:ERROR:|ACPYPE FAILED)',log.read_text(),re.M):
            raise ValueError(name+': parameterization reported an error; inspect '+str(log.relative_to(ROOT)))
        return None
    frcmods=list(folder.glob('*.frcmod'))
    if not frcmods:raise ValueError(name+': no frcmod evidence')
    for p in frcmods:
        if re.search(r'ATTN|need\s+revision',p.read_text(),re.I):
            raise ValueError(name+': force-field terms require review in '+p.name)
    text=itp.read_text();atoms=list(section_rows(text,'atoms'))
    if len(atoms)!=natoms or int(gro.read_text().splitlines()[1])!=natoms:
        raise ValueError(name+': atom-count mismatch')
    total=sum(float(x[6]) for x in atoms)
    if abs(total)>1e-4:raise ValueError(name+': non-neutral charge '+str(total))
    for row in atoms:
        if row[1].lower() in {'du','xx','dummy'}:raise ValueError(name+': untyped atom')
    mol2_files=list(folder.glob('*bcc*gaff2*.mol2'))
    if len(mol2_files)!=1:raise ValueError(name+': expected one final BCC/GAFF2 mol2')
    source=Chem.MolFromMolFile(str(input_path),removeHs=False)
    target=Chem.MolFromMol2File(str(mol2_files[0]),sanitize=False,removeHs=False)
    if target is None:raise ValueError(name+': cannot inspect atom order')
    def graph(m):
        return [a.GetAtomicNum() for a in m.GetAtoms()],{tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in m.GetBonds()}
    if graph(source)!=graph(target):raise ValueError(name+': input/topology atom order or connectivity changed')
    molecule=list(section_rows(text,'moleculetype'))[0][0]
    types={r[0]:r for r in section_rows(text,'atomtypes')}
    match=re.search(r'^\s*\[\s*moleculetype\s*\]',text,re.M|re.I)
    if not match:raise ValueError(name+': missing moleculetype')
    return {'name':name,'folder':folder,'itp':itp,'gro':gro,'mol2':mol2_files[0],
            'charge':total,'molecule':molecule,'types':types,'atoms':atoms,'body':text[match.start():], 'frcmods':frcmods}

def status(value):
    STATUS.parent.mkdir(parents=True,exist_ok=True)
    STATUS.write_text(json.dumps(value,indent=2)+'\n')
    print(json.dumps(value),flush=True)

def run(args,cwd,name):
    with (cwd/(name+'.log')).open('w') as log:
        subprocess.run(args,cwd=cwd,stdout=log,stderr=subprocess.STDOUT,check=True)

def execute(pc,graft,patch):
    out=ROOT/'runs/surface-validation'
    if out.exists():raise ValueError('Existing surface-validation run preserved; use a fresh run location before rerunning')
    out.mkdir()
    types={}
    for item in [pc,graft]:
        for name,row in item['types'].items():
            if name in types and row!=types[name]:raise ValueError('Conflicting atom type: '+name)
            types[name]=row
        (out/(item['name']+'.itp')).write_text(item['body'])
    (out/'gaff2-types.itp').write_text('[ atomtypes ]\n'+'\n'.join(' '.join(v) for v in types.values())+'\n')
    (out/'topol.top').write_text('#include "amber99sb-ildn.ff/forcefield.itp"\n#include "gaff2-types.itp"\n#include "PC4.itp"\n#include "GRAFT.itp"\n#include "amber99sb-ildn.ff/tip3p.itp"\n[ system ]\nExploratory single-graft short-PC patch\n[ molecules ]\n'+graft['molecule']+' 1\n'+pc['molecule']+' '+str(patch['free_pc_oligomers'])+'\n')
    coords=ROOT/'models/pc-patch/pc4-grafted-cop10-patch.pdb'
    records=[x for x in coords.read_text().splitlines() if x.startswith(('ATOM  ','HETATM'))]
    expected=445+132*patch['free_pc_oligomers']
    if len(records)!=expected:raise ValueError('Patch molecule count does not match topology')
    xyz=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in records])
    lateral=patch['patch_region_A'][0]/10
    height=max(8.,float(np.ptp(xyz[:,2]))/10+3.)
    prepared=ROOT/'models/solvated-patch/preparation.json'
    reused=False
    if prepared.exists():
        evidence=json.loads(prepared.read_text())
        wet=prepared.parent/'solvated.gro'
        if (hashlib.sha256(coords.read_bytes()).hexdigest()==evidence.get('source_patch_sha256')
            and wet.exists() and hashlib.sha256(wet.read_bytes()).hexdigest()==evidence.get('solvated_gro_sha256')
            and evidence.get('dry_patch_atoms')==expected):
            shutil.copy2(wet,out/'solvated.gro')
            with (out/'topol.top').open('a') as topology:topology.write('SOL '+str(evidence['water_molecules'])+'\n')
            reused=True
    if not reused:
        run(['gmx','editconf','-f',str(coords),'-o','dry.gro','-box',str(lateral),str(lateral),str(height),'-c'],out,'box')
        run(['gmx','solvate','-cp','dry.gro','-cs','spc216.gro','-p','topol.top','-o','solvated.gro'],out,'solvate')
    # Graph/element ordering was audited above. Normalize only names/residue IDs
    # to the actual generated topology; no coordinate, element or charge changes.
    wet_path=out/'solvated.gro';gro_lines=wet_path.read_text().splitlines()
    index=2;residue=1
    for item in [graft]+[pc]*patch['free_pc_oligomers']:
        for atom in item['atoms']:
            line=gro_lines[index]
            gro_lines[index]=f'{residue%100000:5d}{atom[3][:5]:<5}{atom[4][:5]:>5}'+line[15:]
            index+=1
        residue+=1
    if index!=2+expected:raise ValueError('Coordinate naming map length mismatch')
    wet_path.write_text('\n'.join(gro_lines)+'\n')
    (out/'coordinate-preparation.json').write_text(json.dumps({'reused_hash_verified_solvation':reused,'solute_atoms_name_mapped_after_graph_audit':expected},indent=2)+'\n')
    (out/'em.mdp').write_text('integrator = steep\nnsteps = 5000\nemtol = 1000\nemstep = 0.005\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n')
    run(['gmx','grompp','-f','em.mdp','-c','solvated.gro','-p','topol.top','-o','em.tpr'],out,'em-prep')
    run(['gmx','mdrun','-deffnm','em','-ntmpi','1','-ntomp','2'],out,'em-run')
    emlog=(out/'em.log').read_text()
    if 'converged to Fmax <' not in emlog:raise ValueError('Energy minimization did not meet the force threshold; no MD started')
    (out/'check.mdp').write_text('integrator = md\ndt = 0.0005\nnsteps = 400\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\nconstraints = h-bonds\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\npcoupl = no\ngen-vel = yes\ngen-temp = 300\ngen-seed = 2026\nnstenergy = 20\npbc = xyz\n')
    run(['gmx','grompp','-f','check.mdp','-c','em.gro','-p','topol.top','-o','check.tpr'],out,'md-prep')
    run(['gmx','mdrun','-deffnm','check','-ntmpi','1','-ntomp','2'],out,'md-run')
    run(['gmx','check','-e','check.edr'],out,'energy-check')
    if not (out/'check.gro').is_file():raise ValueError('No final MD coordinates')
    dest=ROOT/'models/forcefield';dest.mkdir(exist_ok=True)
    manifest={}
    for item in [pc,graft]:
        for file in [item['itp'],item['gro'],item['mol2'],*item['frcmods']]:
            target=dest/(item['name']+'-'+file.name);shutil.copy2(file,target)
            manifest[target.name]={'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'source':str(file.relative_to(ROOT)),'net_charge':item['charge']}
    (dest/'provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return {'state':'passed_initialization_check','output':str(out.relative_to(ROOT)),
            'em':'force threshold met','md':'400 steps x0.0005ps =0.2ps at300K','protein_present':False,
            'limitation':'initialization only: low-density short-OH-ended BPA-PC patch, not an equilibrated experimental coating or validated adsorption model'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--watch',action='store_true');p.add_argument('--run',action='store_true');p.add_argument('--max-wait',type=int,default=1800);a=p.parse_args()
    deadline=time.monotonic()+a.max_wait
    previous=None
    lock=None
    if a.watch or a.run:
        STATUS.parent.mkdir(parents=True,exist_ok=True)
        lock=(STATUS.parent/'surface-validation.lock').open('w')
        try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            print('Another surface-validation worker is active; no second run started',flush=True);return
    while True:
        try:
            pc=parameter_gate('PC4');graft=parameter_gate('GRAFT')
            patch_file=ROOT/'models/pc-patch/patch.json';patch=json.loads(patch_file.read_text()) if patch_file.exists() else None
            if not a.watch and not a.run:
                print(json.dumps({'state':'inspection_only','PC4_ready':pc is not None,'GRAFT_ready':graft is not None,'packed_patch_ready':patch is not None}),flush=True);return
            if pc and graft and patch:
                status({'state':'running_parameter_and_surface_checks'});status(execute(pc,graft,patch));return
            if not a.watch or time.monotonic()>=deadline:
                status({'state':'prerequisites_pending','PC4_ready':pc is not None,'GRAFT_ready':graft is not None,'packed_patch_ready':patch is not None});return
            current={'state':'waiting_for_quantum_parameters_and_packing','PC4_ready':pc is not None,'GRAFT_ready':graft is not None,'packed_patch_ready':patch is not None}
            if current!=previous:status(current);previous=current
        except Exception as e:
            status({'state':'blocked','reason':str(e),'md_success':False});raise
        time.sleep(10)
if __name__=='__main__':main()
