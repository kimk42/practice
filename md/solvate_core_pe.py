"""Solvate and neutralize the large dry core/finite-PE model, without claiming equilibration."""
import argparse,hashlib,json,os,shutil,subprocess,time
from pathlib import Path
import numpy as np
from openmm import unit
from openmm.app import AmberPrmtopFile,AmberInpcrdFile,PDBxFile
try:
    from .extend_pe_npt import atomic_json
except ImportError:
    from extend_pe_npt import atomic_json

try:
    from .project_paths import amber_home, resolve_project_path
except ImportError:
    from project_paths import amber_home, resolve_project_path


def solvate(dry,output,amberhome=None):
    amberhome=amberhome or amber_home()
    dry,out,amber=map(lambda x:Path(x).resolve(),(dry,output,amberhome));out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((dry/'manifest.json').read_text())
    for name in ['dry.prmtop','dry.inpcrd']:
        if hashlib.sha256((dry/name).read_bytes()).hexdigest()!=manifest['sha256'][name]:raise ValueError('Dry input checksum differs')
    started=time.monotonic()
    # Recreate the reviewed dry unit using its exact LEaP recipe, including 25 disulfides.
    # LEaP does not load prmtop/inpcrd as a unit; avoid an unsupported loadAmberParm command.
    recipe=(dry/'build.leap').read_text().split('saveamberparm model dry.prmtop dry.inpcrd')[0]
    portable_lines=[]
    for line in recipe.splitlines():
        if line.startswith('loadamberparams '):line='loadamberparams '+str(resolve_project_path(line.split(' ',1)[1]))
        elif line.startswith('PEO = loadmol2 '):line='PEO = loadmol2 '+str(resolve_project_path(line.removeprefix('PEO = loadmol2 ')))
        portable_lines.append(line)
    recipe='\n'.join(portable_lines)+'\n'
    if recipe.count('bond model.')!=25:raise ValueError('Unexpected disulfide recipe')
    shutil.copyfile(dry/'dry-complex.pdb',out/'dry-complex.pdb')
    script=recipe+'solvatebox model TIP3PBOX 10.0\naddionsrand model Na+ 0\nsaveamberparm model solvated.prmtop solvated.inpcrd\nquit\n'
    (out/'solvate.leap').write_text(script)
    atomic_json(out/'progress.json',{'status':'running_solvation','source':str(dry),'padding_A':10,'water':'TIP3P','ions':'neutralization only; no prescribed 150 mM salt'})
    try:
        with (out/'tleap.log').open('w') as f:subprocess.run([str(amber/'bin/tleap'),'-f','solvate.leap'],cwd=out,env={**os.environ,'AMBERHOME':str(amber)},stdout=f,stderr=subprocess.STDOUT,check=True)
        prm=AmberPrmtopFile(str(out/'solvated.prmtop'));coord=AmberInpcrdFile(str(out/'solvated.inpcrd'));original=AmberPrmtopFile(str(dry/'dry.prmtop'));origcoord=AmberInpcrdFile(str(dry/'dry.inpcrd'))
        n=original.topology.getNumAtoms();atoms=list(prm.topology.atoms());old=list(original.topology.atoms())
        if [(a.name,a.residue.name) for a in atoms[:n]]!=[(a.name,a.residue.name) for a in old]:raise ValueError('Dry atom order/identity changed')
        # LEaP translates the full solute to its box; subtract the common translation.
        p=np.asarray(coord.positions.value_in_unit(unit.nanometer));q=np.asarray(origcoord.positions.value_in_unit(unit.nanometer));delta=p[:n]-q;maxerr=float(np.linalg.norm(delta-delta.mean(axis=0),axis=1).max())
        if maxerr>1e-5:raise ValueError('Dry solute geometry changed during solvation')
        counts={}
        for r in prm.topology.residues():counts[r.name]=counts.get(r.name,0)+1
        ss={tuple(sorted((a.index,b.index))) for a,b in prm.topology.bonds() if a.name==b.name=='SG'};oldss={tuple(sorted((a.index,b.index))) for a,b in original.topology.bonds() if a.name==b.name=='SG'}
        if ss!=oldss or len(ss)!=25:raise ValueError('Source disulfides changed')
        charge=float(sum(prm._prmtop.getCharges()))
        if abs(charge)>1e-4:raise ValueError(f'Non-neutral solvated model: {charge}')
        if counts.get('PEO')!=690 or counts.get('ACE')!=6 or counts.get('NME')!=6:raise ValueError('Dry-model identity counts changed')
        if coord.boxVectors is None or len(p)!=len(atoms):raise ValueError('Missing periodic cell or positions')
        box=np.array([v.value_in_unit(unit.nanometer) for v in coord.boxVectors]);volume=float(np.linalg.det(box));mass=sum(prm._prmtop.getMasses())
        result={'status':'solvated_and_neutralized_input_not_equilibrated','atoms':len(atoms),'dry_atoms':n,'water_molecules':counts.get('HOH',counts.get('WAT',0)),'residue_counts':counts,'box_vectors_nm':box.tolist(),'volume_nm3':volume,'initial_total_density_g_cm3':mass*1.66053906660e-3/volume,'net_charge_e':charge,'disulfides':len(ss),'max_dry_geometry_change_after_translation_nm':maxerr,'elapsed_seconds':time.monotonic()-started,'source_dry_manifest':manifest,'equilibrium_established':False,'production_ready':False,'limitations':['Observed six-chain core, not full-length fibrinogen; alphaC and other unresolved tails absent','Deglycosylated approximation; caps do not restore missing segments','Only neutralizing sodium ions, not 150 mM salt','Finite oligomer PE patch, not validated crystalline bulk surface','Solvation and neutralization are input preparation, not energy minimization, MD, adsorption or binding-energy results'],'sha256':{name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['solvated.prmtop','solvated.inpcrd','solvate.leap']}}
        # Full solvated coordinates use mmCIF rather than overflowing PDB atom serial numbers.
        prm.topology.setPeriodicBoxVectors(coord.boxVectors)
        with (out/'solvated.cif').open('w') as f:PDBxFile.writeFile(prm.topology,coord.positions,f)
        result['sha256']['solvated.cif']=hashlib.sha256((out/'solvated.cif').read_bytes()).hexdigest()
        from scipy.spatial import cKDTree
        heavy=np.array([a.index for a in atoms[:n] if a.element.symbol!='H'])
        water_oxygen=np.array([a.index for a in atoms[n:] if a.residue.name in ('HOH','WAT') and a.element.symbol=='O'])
        lengths=np.diag(box)
        distances=cKDTree(np.mod(p[heavy],lengths),boxsize=lengths).query(np.mod(p[water_oxygen],lengths))[0]
        result['min_water_oxygen_dry_solute_heavy_distance_nm']=float(distances.min())
        if distances.min()<.12:raise ValueError('Severe water/solute overlap after solvation')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        art=Path(__file__).resolve().parents[1]/'artifacts/core-solvated';art.mkdir(parents=True,exist_ok=True)
        ca=np.array([a.index for a in atoms[:n] if a.name=='CA' and a.element.symbol=='C'])
        carbon=np.array([a.index for a in atoms[:n] if a.residue.name=='PEO' and a.element.symbol=='C'])
        sodium=np.array([a.index for a in atoms if a.element.symbol=='Na'])
        subset=water_oxygen[::max(1,len(water_oxygen)//2500)]
        fig,axes=plt.subplots(2,1,figsize=(14,7))
        for ax,(i,j,title) in zip(axes,[(0,2,'Side view (X/Z)'),(0,1,'Top view (X/Y)')]):
            ax.scatter(p[subset,i],p[subset,j],s=1,c='#85bde0',alpha=.35,label='Water O (subsample)')
            ax.scatter(p[carbon,i],p[carbon,j],s=1,c='#cc934c',alpha=.55,label='PE carbons')
            ax.scatter(p[ca,i],p[ca,j],s=4,c='#225c93',label='Protein CA')
            if len(sodium):ax.scatter(p[sodium,i],p[sodium,j],s=18,c='#883ba0',label='Na+')
            ax.set(xlabel='X (nm)',ylabel=('Z' if j==2 else 'Y')+' (nm)',title=title);ax.set_aspect('equal',adjustable='datalim');ax.grid(alpha=.15)
        axes[0].legend(loc='upper right',ncol=4,fontsize=8);fig.suptitle(f'Actual solvated six-chain observed core + finite PE patch: {len(atoms):,} atoms')
        fig.text(.5,.015,'Input preparation only — no large-system minimization, equilibration or adsorption MD; water dots are a visual subsample.',ha='center',fontsize=10)
        fig.tight_layout(rect=(0,.04,1,.96));fig.savefig(art/'model.png',dpi=150);fig.savefig(art/'model.svg');plt.close(fig)
        result['visualization']='artifacts/core-solvated/model.png'
        atomic_json(out/'manifest.json',result);atomic_json(out/'progress.json',{'status':'completed','atoms':len(atoms),'water_molecules':result['water_molecules']});print(json.dumps({k:v for k,v in result.items() if k!='source_dry_manifest'},indent=2),flush=True);return result
    except BaseException as e:
        atomic_json(out/'failure.json',{'status':'failed','error':f'{type(e).__name__}: {e}'});raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('dry');p.add_argument('output');a=p.parse_args();solvate(a.dry,a.output)
