"""Build an exploratory finite PE oligomer patch + real fibrinogen fragment + explicit TIP3P water.
No bulk-PE or adsorption validation is implied by successful model construction.
"""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
import numpy as np
from openmm import unit,Vec3,NonbondedForce
from openmm.app import PDBFile,AmberPrmtopFile,AmberInpcrdFile,Topology,element,ForceField,NoCutoff

try:
    from .project_paths import amber_home, resolve_project_path
except ImportError:
    from project_paths import amber_home, resolve_project_path


def principal_orientation(xyz):
    xyz=np.asarray(xyz,dtype=float)
    centered=xyz-xyz.mean(axis=0)
    _,_,vt=np.linalg.svd(centered,full_matrices=False)
    basis=vt.T
    if np.linalg.det(basis)<0:basis[:,-1]*=-1
    return centered@basis


def build(protein_path,polymer_dir,output,amberhome):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    protein_path=Path(protein_path).resolve();polymer_dir=Path(polymer_dir).resolve();amber=Path(amberhome).resolve()
    p=PDBFile(str(protein_path));pe=AmberPrmtopFile(str(polymer_dir/'polymer.prmtop'));x=AmberInpcrdFile(str(polymer_dir/'polymer.inpcrd'))
    protein=principal_orientation(p.positions.value_in_unit(unit.nanometer))
    poly=np.array(x.positions.value_in_unit(unit.nanometer));carbons=[a.index for a in pe.topology.atoms() if a.element==element.carbon]
    # Align the extended oligomer carbon backbone along X.
    center=poly[carbons].mean(axis=0);_,_,vt=np.linalg.svd(poly[carbons]-center,full_matrices=False);basis=vt.T
    if np.linalg.det(basis)<0:basis[:,-1]*=-1
    poly=(poly-center)@basis
    # Finite patch, not an infinite periodic crystal. Two columns, twelve rows, two layers.
    nx,ny,nz=2,12,2;dx=np.ptp(poly[:,0])+.30;dy=.48;dz=.44
    slab=np.vstack([poly+[(ix-(nx-1)/2)*dx,(iy-(ny-1)/2)*dy,iz*dz] for iz in range(nz) for iy in range(ny) for ix in range(nx)])
    heavy_protein=np.array([a.index for a in p.topology.atoms() if a.element!=element.hydrogen])
    heavy_poly=np.array([a.index for a in pe.topology.atoms() if a.element!=element.hydrogen])
    upper_heavy=max((poly[heavy_poly,2]+iz*dz).max() for iz in range(nz))
    initial_gap_nm=.5
    protein[:,2]+=upper_heavy+initial_gap_nm-protein[heavy_protein,2].min()
    top=Topology();positions=[];mapping={};chain=top.addChain('A')
    disulfide_res={a.residue for a,b in p.topology.bonds() if a.name==b.name=='SG'}|{b.residue for a,b in p.topology.bonds() if a.name==b.name=='SG'}
    residue_map=[]
    for residue in p.topology.residues():
        names={a.name for a in residue.atoms()};name=residue.name
        if residue in disulfide_res:name='CYX'
        if name=='HIS':name='HIP' if {'HD1','HE2'}<=names else 'HID' if 'HD1' in names else 'HIE'
        r=top.addResidue(name,chain,id=residue.id);residue_map.append({'output_residue':r.index+1,'source_fragment_residue':residue.id,'template':name})
        for atom in residue.atoms():
            mapping[atom]=top.addAtom(atom.name,atom.element,r);positions.append(protein[atom.index])
    for a,b in p.topology.bonds():top.addBond(mapping[a],mapping[b])
    disulfides=[(a.residue.index+1,b.residue.index+1) for a,b in p.topology.bonds() if a.name==b.name=='SG']
    # OpenMM PDB exporter normalizes HIS/CYS names; restore reviewed Amber template names below.
    for copy in range(nx*ny*nz):
        chain=top.addChain('P');r=top.addResidue('PEO',chain);offset=copy*len(poly);m={}
        for atom in pe.topology.atoms():m[atom]=top.addAtom(atom.name,atom.element,r);positions.append(slab[offset+atom.index])
        for a,b in pe.topology.bonds():top.addBond(m[a],m[b])
    with (out/'complex.pdb').open('w') as f:PDBFile.writeFile(top,np.array(positions)*unit.nanometer,f)
    lines=[line for line in (out/'complex.pdb').read_text().splitlines() if not line.startswith('CONECT')];names={i+1:r['template'] for i,r in enumerate(residue_map)}
    for i,line in enumerate(lines):
        if line.startswith('ATOM  '):
            number=int(line[22:26])
            if number in names:
                line=line[:17]+names[number].rjust(3)+line[20:]
                if number==1 and line[12:16].strip()=='H':line=line[:12]+' H1 '+line[16:]
                lines[i]=line
    (out/'complex.pdb').write_text('\n'.join(lines)+'\n')
    checked=PDBFile(str(out/'complex.pdb'))
    if sum(r.name=='PEO' for r in checked.topology.residues())!=nx*ny*nz:raise ValueError('PDB export lost PE identities')
    # Bind source parameter hashes to the prepared combined system.
    env={**os.environ,'AMBERHOME':str(amber),'OMP_NUM_THREADS':'4'}
    leap=f'''source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams {polymer_dir/'PE.frcmod'}
PEO = loadmol2 {polymer_dir/'PE.mol2'}
complex = loadpdb complex.pdb
'''
    leap+=''.join(f'bond complex.{a}.SG complex.{b}.SG\n' for a,b in disulfides)
    protein_system=ForceField('amber14-all.xml').createSystem(p.topology,nonbondedMethod=NoCutoff,constraints=None)
    protein_nb=next(f for f in protein_system.getForces() if isinstance(f,NonbondedForce))
    protein_charge=sum(protein_nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(protein_system.getNumParticles()))
    neutralizer='Na+' if protein_charge<-.1 else 'Cl-' if protein_charge>.1 else None
    leap+='check complex\nsolvatebox complex TIP3PBOX 10.0\n'
    if neutralizer:leap+=f'addionsrand complex {neutralizer} 0\n'
    leap+='saveamberparm complex complex.prmtop complex.inpcrd\nsavepdb complex solvated.pdb\nquit\n'
    (out/'build.leap').write_text(leap)
    with (out/'tleap.log').open('w') as f:subprocess.run([str(amber/'bin/tleap'),'-f','build.leap'],cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
    parameterized=AmberPrmtopFile(str(out/'complex.prmtop'));coords=AmberInpcrdFile(str(out/'complex.inpcrd'))
    if coords.boxVectors is None:raise ValueError('Missing periodic box')
    counts={}
    for r in parameterized.topology.residues():counts[r.name]=counts.get(r.name,0)+1
    if counts.get('PEO')!=nx*ny*nz:raise ValueError('Lost polymer residues')
    manifest={'status':'prepared','protein':'actual 3GHG gamma-chain fragment, source C label sequence 143–394; artificial termini',
              'protein_source_sha256':hashlib.sha256(protein_path.read_bytes()).hexdigest(),'polymer':'PE','polymer_model':'finite C24H50 oligomer patch','polymer_chains':nx*ny*nz,
              'patch_layout':{'nx':nx,'ny':ny,'layers':nz,'dx_nm':float(dx),'dy_nm':dy,'dz_nm':dz},'initial_heavy_atom_plane_gap_nm':initial_gap_nm,
              'protein_forcefield':'Amber ff14SB','polymer_forcefield':'GAFF2 / AM1-BCC','water_model':'TIP3P explicit water','ions':'neutralizing Na+/Cl- only; no prescribed bulk salt concentration',
              'particles':parameterized.topology.getNumAtoms(),'residue_counts':counts,'protein_residue_mapping':residue_map,
              'sha256':{n:hashlib.sha256((out/n).read_bytes()).hexdigest() for n in ['complex.prmtop','complex.inpcrd','complex.pdb','build.leap']},
              'polymer_parameters_sha256':json.loads((polymer_dir/'manifest.json').read_text())['sha256'],
              'limitations':['Finite oligomer patch, not validated crystalline bulk PE','Finite patch edges and methyl chain ends are present','No polymer phase or macroscopic crystallinity claim','Protein is an isolated gamma fragment, not intact fibrinogen','No adsorption/affinity conclusion from a short pilot']}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps({k:v for k,v in manifest.items() if k!='protein_residue_mapping'},indent=2))
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('protein');p.add_argument('polymer');p.add_argument('output');p.add_argument('--amberhome',default=str(amber_home()));a=p.parse_args();build(a.protein,a.polymer,a.output,a.amberhome)
