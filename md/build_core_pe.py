"""Prepare a size-matched dry PE patch + capped six-chain observed fibrinogen core.
This is a parameterized preparation artifact, not solvated/equilibrated adsorption MD.
"""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
import numpy as np
from openmm import unit
from openmm.app import AmberPrmtopFile,AmberInpcrdFile,NoCutoff,element

try:
    from .project_paths import amber_home, resolve_project_path
except ImportError:
    from project_paths import amber_home, resolve_project_path


def build(core,polymer,output,amberhome):
    core=Path(core).resolve();polymer=Path(polymer).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    metadata=json.loads((core/'manifest.json').read_text());relax=core/'local-relaxation-restrained'
    report=json.loads((relax/'report.json').read_text())
    if report['status']!='local_preparation_completed_not_equilibrated' or report['nonbonded_pairs_under_1_2A_after']!=0:
        raise ValueError('Core local-relaxation diagnostics require resolution')
    protein=AmberPrmtopFile(str(core/'core.prmtop'));coords=AmberInpcrdFile(str(relax/'core-relaxed.inpcrd'))
    pe=AmberPrmtopFile(str(polymer/'polymer.prmtop'));px=AmberInpcrdFile(str(polymer/'polymer.inpcrd'))
    xyz=np.array(coords.positions.value_in_unit(unit.nanometer));atoms=list(protein.topology.atoms());maps=metadata['residue_mapping']
    ca=[a.index for a in atoms if a.name=='CA' and maps[a.residue.index]['source_label_seq'] is not None]
    center=xyz[ca].mean(axis=0);_,_,vt=np.linalg.svd(xyz[ca]-center,full_matrices=False);basis=vt.T
    if np.linalg.det(basis)<0:basis[:,-1]*=-1
    xyz=(xyz-center)@basis
    pxyz=np.array(px.positions.value_in_unit(unit.nanometer));pc=[a.index for a in pe.topology.atoms() if a.element==element.carbon]
    pcenter=pxyz[pc].mean(axis=0);_,_,vt=np.linalg.svd(pxyz[pc]-pcenter,full_matrices=False);pbasis=vt.T
    if np.linalg.det(pbasis)<0:pbasis[:,-1]*=-1
    pxyz=(pxyz-pcenter)@pbasis
    dx=np.ptp(pxyz[:,0])+.30;dy=.48;dz=.44
    ext=np.ptp(xyz,axis=0);nx=int(np.ceil((ext[0]+4)/dx));ny=int(np.ceil((ext[1]+4)/dy));layers=2
    pe_carbon_top=pxyz[pc,2].max()+dz
    xyz[:,2]+=pe_carbon_top+.5-xyz[:,2].min()
    native_protein_names=protein._prmtop.getAtomNames();native_pe_names=pe._prmtop.getAtomNames()
    serial=0;lines=[]
    def atomline(name,resname,chain,resid,pos,elm):
        nonlocal serial
        serial+=1
        if serial>99999:raise ValueError('Dry model exceeds PDB serial limit; use mmCIF writer before scaling further')
        x,y,z=pos*10;field=name.ljust(4) if len(elm)==2 or len(name)==4 else (' '+name).ljust(4)
        lines.append(f'ATOM  {serial:5d} {field} {resname:3s} {chain}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elm:>2s}')
    for r in protein.topology.residues():
        m=maps[r.index]
        for a in r.atoms():atomline(native_protein_names[a.index],m['name'],m['source_chain'],r.index+1,xyz[a.index],a.element.symbol)
        if m['name'] in ('NME','CA'):lines.append('TER')
    count=0;lower_res=[]
    for layer in range(layers):
        for row in range(ny):
            for column in range(nx):
                resid=len(maps)+count+1;offset=np.array([(column-(nx-1)/2)*dx,(row-(ny-1)/2)*dy,layer*dz])
                if layer==0:lower_res.append(resid)
                for a in pe.topology.atoms():atomline(native_pe_names[a.index],'PEO','P',resid,pxyz[a.index]+offset,a.element.symbol)
                lines.append('TER');count+=1
    allowed_names={m['name'] for m in maps}|{'PEO'}
    if any(l[17:20].strip() not in allowed_names for l in lines if l.startswith('ATOM')):raise ValueError('PDB column alignment invalid')
    lines.append('END');(out/'dry-complex.pdb').write_text('\n'.join(lines)+'\n')
    ss=[(a.residue.index+1,b.residue.index+1) for a,b in protein.topology.bonds() if a.name==b.name=='SG']
    leap=f'source leaprc.protein.ff14SB\nsource leaprc.gaff2\nsource leaprc.water.tip3p\nloadamberparams {polymer/"PE.frcmod"}\nPEO = loadmol2 {polymer/"PE.mol2"}\nmodel = loadpdb dry-complex.pdb\n'
    leap+=''.join(f'bond model.{a}.SG model.{b}.SG\n' for a,b in ss)
    leap+='saveamberparm model dry.prmtop dry.inpcrd\nquit\n';(out/'build.leap').write_text(leap)
    with (out/'tleap.log').open('w') as log:subprocess.run([str(Path(amberhome)/'bin/tleap'),'-f','build.leap'],cwd=out,env={**os.environ,'AMBERHOME':str(amberhome)},stdout=log,stderr=subprocess.STDOUT,check=True)
    prm=AmberPrmtopFile(str(out/'dry.prmtop'));system=prm.createSystem(nonbondedMethod=NoCutoff,constraints=None)
    actual_ss={(min(a.residue.index+1,b.residue.index+1),max(a.residue.index+1,b.residue.index+1)) for a,b in prm.topology.bonds() if a.name==b.name=='SG'}
    if actual_ss!={(min(a,b),max(a,b)) for a,b in ss}:raise ValueError('Disulfides changed in mixed model')
    if sum(r.name=='PEO' for r in prm.topology.residues())!=count:raise ValueError('Polymer identity count mismatch')
    box=np.array([nx*dx+2,ny*dy+2,float(xyz[:,2].max()-pxyz[:,2].min())+2])
    expected_water=float(np.prod(box)*33.4)
    summary={'status':'parameterized_dry_model_not_MD_ready','protein_model':'six-chain observed 3GHG core, neutral caps, calcium, no glycans','source_core':str(core),'source_core_prmtop_sha256':metadata['sha256']['core.prmtop'],
             'source_relaxed_coordinates_sha256':report['relaxed_inpcrd_sha256'],'observed_protein_residues':1913,'protein_caps':12,'disulfides':25,'calcium':4,
             'polymer':'PE finite C24H50 oligomer patch; not validated bulk crystal','polymer_chains':count,'patch_grid':[nx,ny,layers],'patch_pitch_nm':[float(dx),dy,dz],
             'dry_atoms':system.getNumParticles(),'rough_solvation_box_nm':box.tolist(),'rough_upper_box_water_molecules_before_solute_exclusion':expected_water,
             'rough_total_atoms_before_solute_exclusion':float(3*expected_water+system.getNumParticles()),'lower_layer_residue_ids_for_future_restraints':lower_res,
             'next_required':['Solvation and neutralization','Diagnose remaining local geometry and restrained minimization in explicit solvent','Pressure/density/interface equilibration','Repeated orientations and seeds before adsorption interpretation'],
             'limitations':metadata['limitations']+['Dry coordinates and parameter files are not adsorption dynamics or equilibrium results'],
             'sha256':{name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['dry-complex.pdb','dry.prmtop','dry.inpcrd','build.leap']}}
    (out/'manifest.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k!='lower_layer_residue_ids_for_future_restraints'},indent=2));return summary
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('core');p.add_argument('polymer');p.add_argument('output');p.add_argument('--amberhome',default=str(amber_home()));a=p.parse_args();build(a.core,a.polymer,a.output,a.amberhome)
