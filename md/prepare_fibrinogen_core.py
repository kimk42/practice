"""Prepare the observed six-chain 3GHG core with neutral caps and source-defined disulfides.
Missing alphaC/tails are NOT reconstructed. Glycans/knob peptides are omitted explicitly.
"""
import argparse,hashlib,json,os,subprocess
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from openmm.app import AmberPrmtopFile,AmberInpcrdFile,NoCutoff
from openmm import NonbondedForce,unit

try:
    from .project_paths import amber_home, resolve_project_path
except ImportError:
    from project_paths import amber_home, resolve_project_path


def unit_vector(x):return x/np.linalg.norm(x)

def internal_atom(a,b,c,length,angle_deg,dihedral_deg):
    # Coordinates in Å; a,b,c are the previous three atoms.
    bc=unit_vector(c-b);normal=unit_vector(np.cross(b-a,bc));perp=np.cross(normal,bc)
    theta,phi=np.deg2rad([angle_deg,dihedral_deg])
    return c+length*(-np.cos(theta)*bc+np.sin(theta)*(np.cos(phi)*perp+np.sin(phi)*normal))

def acetyl_cap(residue,other_positions):
    n,ca,c=[residue[k]['xyz'] for k in ('N','CA','C')]
    best=None
    tree=cKDTree(other_positions)
    for phi in range(-180,180,15):
        carbon=internal_atom(c,ca,n,1.329,121.7,phi)
        methyl=internal_atom(ca,n,carbon,1.522,116.2,180)
        oxygen=internal_atom(ca,n,carbon,1.229,123.0,0)
        coords=np.array([carbon,methyl,oxygen]);clearance=float(tree.query(coords)[0].min())
        if best is None or clearance>best[0]:best=(clearance,{'C':carbon,'CH3':methyl,'O':oxygen},phi)
    return best

def methylamide_cap(residue):
    ca,c,o=[residue[k]['xyz'] for k in ('CA','C','O')]
    along=unit_vector(ca-c);perp=unit_vector((o-c)-np.dot(o-c,along)*along)
    n=c+1.335*(np.cos(np.deg2rad(116.2))*along-np.sin(np.deg2rad(116.2))*perp)
    methyl=internal_atom(ca,c,n,1.458,121.7,180)
    return {'N':n,'C':methyl}


def prepare(source,output,amberhome):
    source=Path(source).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    c=MMCIF2Dict(str(source))
    if c['_entry.id'][0]!='3GHG':raise ValueError('This protocol is reviewed specifically for 3GHG')
    j=c['_pdbx_struct_assembly_gen.assembly_id'].index('1')
    if c['_pdbx_struct_assembly_gen.oper_expression'][j]!='1':raise ValueError('Assembly transform must be implemented explicitly')
    selected=set(c['_pdbx_struct_assembly_gen.asym_id_list'][j].split(','))
    if not set('ABCDEF')<=selected:raise ValueError('Unexpected first assembly')
    residues={ch:defaultdict(dict) for ch in 'ABCDEF'};calcium=[]
    for i,ch in enumerate(c['_atom_site.label_asym_id']):
        if ch not in selected or c['_atom_site.pdbx_PDB_model_num'][i]!='1':continue
        atom=c['_atom_site.label_atom_id'][i];element=c['_atom_site.type_symbol'][i]
        xyz=np.array([float(c[f'_atom_site.Cartn_{k}'][i]) for k in 'xyz'])
        if element.upper()=='CA':calcium.append({'source_asym':ch,'xyz':xyz});continue
        if ch not in residues or element in ('H','D'):continue
        seq=int(c['_atom_site.label_seq_id'][i]);occ=float(c['_atom_site.occupancy'][i])
        r=residues[ch][seq]
        if atom not in r or r[atom]['occupancy']<occ:
            r[atom]={'xyz':xyz,'occupancy':occ,'name':c['_atom_site.label_comp_id'][i],'element':element,'auth_seq':c['_atom_site.auth_seq_id'][i]}
    if len(calcium)!=4:raise ValueError('Expected four source calcium ions in assembly 1')
    expected=[174,401,381,174,401,382]
    for ch,count in zip('ABCDEF',expected):
        seq=sorted(residues[ch])
        if len(seq)!=count or seq!=list(range(seq[0],seq[-1]+1)):raise ValueError('Unexpected gaps/count in '+ch)
        if any(not {'N','CA','C','O'}<=r.keys() for r in residues[ch].values()):raise ValueError('Incomplete observed backbone')
    disulfides=[]
    for i,t in enumerate(c['_struct_conn.conn_type_id']):
        a,b=[c[f'_struct_conn.ptnr{k}_label_asym_id'][i] for k in (1,2)]
        if t=='disulf' and a in residues and b in residues:
            ra,rb=[int(c[f'_struct_conn.ptnr{k}_label_seq_id'][i]) for k in (1,2)]
            if 'SG' not in residues[a][ra] or 'SG' not in residues[b][rb]:raise ValueError('Missing disulfide atom')
            disulfides.append((a,ra,b,rb))
    if len(disulfides)!=25:raise ValueError('Unexpected source disulfide count')
    ss_res={(a,ra) for a,ra,b,rb in disulfides}|{(b,rb) for a,ra,b,rb in disulfides}
    all_positions=np.array([a['xyz'] for rs in residues.values() for r in rs.values() for a in r.values()])
    lines=[];mapping=[];serial=0;global_residue=0;lookup={};cap_geometry=[]
    def write_res(chain,name,atoms,source_seq=None,source_auth=None,cap=False):
        nonlocal serial,global_residue
        global_residue+=1
        mapping.append({'output_residue':global_residue,'source_chain':chain,'source_label_seq':source_seq,'source_auth_seq':source_auth,'name':name,'artificial_cap':cap})
        if source_seq is not None:lookup[(chain,source_seq)]=global_residue
        for atom,xyz,elm in atoms:
            serial+=1;x,y,z=xyz
            record='HETATM' if name=='CA' else 'ATOM  '
            atom_field=atom.ljust(4) if len(elm)==2 else (' '+atom).ljust(4)
            lines.append(f'{record}{serial:5d} {atom_field} {name:3s} {chain}{global_residue:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elm:>2s}')
    for ch in 'ABCDEF':
        seq=sorted(residues[ch]);first,last=seq[0],seq[-1]
        # Exclude the terminal residue itself from ACE clearance scoring; retain other chains.
        others=np.array([a['xyz'] for chain,rs in residues.items() for n,r in rs.items() if (chain,n)!=(ch,first) for a in r.values()])
        clearance,ace,phi=acetyl_cap(residues[ch][first],others);nme=methylamide_cap(residues[ch][last])
        cap_geometry.append({'chain':ch,'ACE_phi_deg':phi,'ACE_min_other_residue_distance_A':clearance,'N_cap_C_N_bond_A':float(np.linalg.norm(ace['C']-residues[ch][first]['N']['xyz'])),'C_cap_C_N_bond_A':float(np.linalg.norm(nme['N']-residues[ch][last]['C']['xyz']))})
        write_res(ch,'ACE',[(name,pos,'O' if name=='O' else 'C') for name,pos in ace.items()],cap=True)
        for n in seq:
            r=residues[ch][n];name=next(iter(r.values()))['name'];name='CYX' if (ch,n) in ss_res else 'HIE' if name=='HIS' else name
            write_res(ch,name,[(a,v['xyz'],v['element']) for a,v in r.items() if a!='OXT'],n,next(iter(r.values()))['auth_seq'])
        write_res(ch,'NME',[(name,pos,'N' if name=='N' else 'C') for name,pos in nme.items()],cap=True)
        lines.append('TER')
    for ca in calcium:write_res('I','CA',[('CA',ca['xyz'],'Ca')]);lines.append('TER')
    lines.append('END');(out/'core-heavy-capped.pdb').write_text('\n'.join(lines)+'\n')
    leap='source leaprc.protein.ff14SB\nsource leaprc.water.tip3p\ncore = loadpdb core-heavy-capped.pdb\n'
    leap+=''.join(f'bond core.{lookup[(a,ra)]}.SG core.{lookup[(b,rb)]}.SG\n' for a,ra,b,rb in disulfides)
    leap+='check core\nsaveamberparm core core.prmtop core.inpcrd\nsavepdb core core.pdb\nquit\n'
    (out/'prepare.leap').write_text(leap)
    env={**os.environ,'AMBERHOME':str(Path(amberhome).resolve()),'OMP_NUM_THREADS':'4'}
    with (out/'tleap.log').open('w') as f:subprocess.run([str(Path(amberhome)/'bin/tleap'),'-f','prepare.leap'],cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
    prm=AmberPrmtopFile(str(out/'core.prmtop'));crd=AmberInpcrdFile(str(out/'core.inpcrd'))
    top_res=list(prm.topology.residues())
    ss=[(a.residue.index+1,b.residue.index+1) for a,b in prm.topology.bonds() if a.name==b.name=='SG']
    intended={frozenset((lookup[(a,ra)],lookup[(b,rb)])) for a,ra,b,rb in disulfides}
    if {frozenset(x) for x in ss}!=intended:raise ValueError('Source disulfides not preserved exactly')
    if len(top_res)!=len(mapping) or sum(r.name=='ACE' for r in top_res)!=6 or sum(r.name=='NME' for r in top_res)!=6:raise ValueError('Residue/cap identities changed')
    xyz=np.array(crd.positions.value_in_unit(unit.angstrom));errors=[]
    for atom in prm.topology.atoms():
        m=mapping[atom.residue.index]
        if m['source_label_seq'] is not None and atom.name in residues[m['source_chain']][m['source_label_seq']]:
            original=residues[m['source_chain']][m['source_label_seq']][atom.name]['xyz'];errors.append(float(np.linalg.norm(xyz[atom.index]-original)))
    if not errors or max(errors)>.002:raise ValueError('Observed heavy atoms moved unexpectedly')
    system=prm.createSystem(nonbondedMethod=NoCutoff,constraints=None)
    nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce))
    charge=sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles()))
    backbone=[xyz[a.index] for a in prm.topology.atoms() if a.name=='CA' and mapping[a.residue.index]['source_label_seq'] is not None]
    _,_,vt=np.linalg.svd(np.array(backbone)-np.mean(backbone,axis=0),full_matrices=False)
    extent=np.ptp((np.array(backbone)-np.mean(backbone,axis=0))@vt.T,axis=0)/10
    # A bounding-box water estimate, not an executed solvation or atom count claim.
    box=extent+np.array([4.,6.,5.]);estimated_waters=float(np.prod(box)*33.4)
    manifest={'status':'parameterized_observed_core_not_equilibrated','entry':'3GHG','assembly':'1','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
              'observed_residues':sum(expected),'atoms':system.getNumParticles(),'caps':{'ACE':6,'NME':6},'calcium_ions':len(calcium),'source_disulfides':25,'source_interchain_disulfides':15,'net_charge_e':charge,
              'protein_forcefield':'Amber ff14SB','calcium_parameters':'Amber TIP3P ion parameters (12-6), no explicit waters in this core bundle',
              'residue_mapping':mapping,'cap_geometry':cap_geometry,'observed_heavy_atoms_checked':len(errors),'max_source_heavy_atom_displacement_A':max(errors),
              'CA_principal_extent_nm':extent.tolist(),'rough_explicit_water_box_nm':box.tolist(),'rough_box_water_molecules_before_solute_exclusion':estimated_waters,
              'missing_policy':'Retain only observed continuous core; all absent terminal regions including long alphaC are omitted, not reconstructed. Neutral ACE/NME caps replace artificial terminal charges.',
              'protonation':'Standard ff14SB residue states; HIS uses HIE; no site-specific pKa assessment',
              'glycan_policy':'Deglycosylated approximation: glycans removed from working core only, original retained','knob_policy':'A/B knob peptides excluded',
              'limitations':['Not full-length fibrinogen; long unresolved alphaC and other tails remain absent','Caps reduce artificial terminal charge but do not reproduce missing segments','No solvent equilibration, adsorption or binding-energy result','Glycans and their influence are omitted; calcium treatment is approximate','New surface dimensions must accommodate entire core, not reuse the small gamma-fragment patch'],
              'sha256':{name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['core-heavy-capped.pdb','core.prmtop','core.inpcrd','core.pdb','prepare.leap']}}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps({k:v for k,v in manifest.items() if k not in ['residue_mapping','sha256']},indent=2))
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source');p.add_argument('output');p.add_argument('--amberhome',default=str(amber_home()));a=p.parse_args();prepare(a.source,a.output,a.amberhome)
