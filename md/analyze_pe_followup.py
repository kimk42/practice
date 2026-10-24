"""Geometric and descriptive checks for a finite PE/protein/water preparation.
No equilibrium significance test or adsorption/affinity inference is made.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from openmm import XmlSerializer,unit
from openmm.app import AmberPrmtopFile,element
try:
    from .extend_pe_npt import atomic_json
except ImportError:
    from extend_pe_npt import atomic_json

try:
    from .project_paths import amber_home, resolve_project_path
except ImportError:
    from project_paths import amber_home, resolve_project_path


def window_summary(rows,window_ps=10):
    end=float(rows[-1]['extension_time_ps']); windows=[]
    for left in np.arange(0,end,window_ps):
        subset=[x for x in rows if left<float(x['extension_time_ps'])<=left+window_ps]
        if len(subset)<2:continue
        t=np.array([float(x['extension_time_ps']) for x in subset]);d={'start_ps':float(left),'end_ps':float(min(end,left+window_ps)),'samples':len(subset)}
        for key in ['total_density_g_cm3','volume_nm3','potential_kJ_mol','temperature_K']:
            v=np.array([float(x[key]) for x in subset]);d[key]={'mean':float(v.mean()),'standard_deviation':float(v.std()),'linear_slope_per_ps':float(np.polyfit(t-t[0],v,1)[0])}
        windows.append(d)
    return windows


def structural_snapshot(prm,state,system):
    xyz=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    vectors=np.array([v.value_in_unit(unit.nanometer) for v in state.getPeriodicBoxVectors()]);box=np.diag(vectors)
    if not np.allclose(vectors,np.diag(box)):raise ValueError('Nonorthorhombic box')
    atoms=list(prm.topology.atoms());ca=[a.index for a in atoms if a.residue.index<252 and a.name=='CA'];water=[a.index for a in atoms if a.residue.name in ('HOH','WAT') and a.element==element.oxygen];pe=[a.index for a in atoms if a.residue.name=='PEO' and a.element==element.carbon];ph=[a.index for a in atoms if a.residue.index<252 and a.element!=element.hydrogen]
    def distances(a,b):return cKDTree(np.mod(xyz[b],box),boxsize=box).query(np.mod(xyz[a],box))[0]
    # All metrics here use minimum-image distances, except intramolecular CA Rg.
    water_pe=distances(water,pe);water_protein=distances(water,ph)
    dirs=[];lengths=[];centers=[]
    for residue in prm.topology.residues():
        if residue.name!='PEO':continue
        ids=[a.index for a in residue.atoms() if a.element==element.carbon]
        p=xyz[ids];d=p[-1]-p[0];d-=box*np.round(d/box);norm=np.linalg.norm(d)
        dirs.append(d/norm);lengths.append(norm);centers.append(p.mean(axis=0))
    centers=np.mod(np.array(centers),box)
    nearest_centers=cKDTree(centers,boxsize=box).query(centers,k=2)[0][:,1]
    tree=cKDTree(np.mod(xyz[pe],box),boxsize=box)
    carbon_dist,carbon_index=tree.query(np.mod(xyz[pe],box),k=25)
    carbon_chain=np.array([atoms[i].residue.index for i in pe])
    other_chain=carbon_chain[carbon_index]!=carbon_chain[:,None]
    min_interchain=float(np.where(other_chain,carbon_dist,np.inf).min())
    dirs=np.array(dirs);tensor=1.5*np.einsum('ni,nj->ij',dirs,dirs)/len(dirs)-.5*np.eye(3)
    # A surface-normal profile is referenced to the initial principal slab Z axis, not a validated instantaneous surface.
    zpe=xyz[pe,2];ang=2*np.pi*zpe/box[2];center=np.angle(np.mean(np.exp(1j*ang)))*box[2]/(2*np.pi)
    z=(xyz[water,2]-center+box[2]/2)%box[2]-box[2]/2
    hist,edges=np.histogram(z,bins=np.linspace(-box[2]/2,box[2]/2,41));number_density=hist/(box[0]*box[1]*np.diff(edges))
    constraint_error=0
    for i in range(system.getNumConstraints()):
        a,b,d=system.getConstraintParameters(i);target=d.value_in_unit(unit.nanometer);delta=xyz[a]-xyz[b];delta-=box*np.round(delta/box)
        constraint_error=max(constraint_error,abs(np.linalg.norm(delta)-target)/target)
    caxyz=xyz[ca];rg=np.sqrt(np.mean(np.sum((caxyz-caxyz.mean(axis=0))**2,axis=1)))
    return {'CA_radius_of_gyration_nm':float(rg),'PE_chain_end_to_end_mean_nm':float(np.mean(lengths)),'PE_nearest_chain_COM_distance_mean_nm':float(np.mean(nearest_centers)),'PE_nearest_chain_COM_distance_max_nm':float(np.max(nearest_centers)),'PE_min_interchain_carbon_distance_nm':min_interchain,'PE_chain_end_to_end_min_nm':float(np.min(lengths)),'PE_chain_end_to_end_max_nm':float(np.max(lengths)),'PE_backbone_nematic_order_max_eigenvalue':float(np.linalg.eigvalsh(tensor).max()),'waters_oxygen_within_0_35nm_of_PE_carbon':int(np.sum(water_pe<.35)),'waters_oxygen_within_0_35nm_of_protein_heavy':int(np.sum(water_protein<.35)),'min_water_oxygen_PE_carbon_distance_nm':float(water_pe.min()),'min_water_oxygen_protein_heavy_distance_nm':float(water_protein.min()),'max_relative_constraint_error':float(constraint_error),'water_O_z_profile':{'coordinate':'Z relative to circular mean of PE carbon coordinates, initial slab axis; not a fitted dynamic interface','bin_center_nm':((edges[:-1]+edges[1:])/2).tolist(),'oxygen_number_density_nm_minus3':number_density.tolist()},'counts':{'water_molecules':len(water),'PE_chains':len(lengths),'protein_CA':len(ca)}}


def analyze(run,bundle,output):
    run,bundle,out=map(Path,(run,bundle,output));out.mkdir(parents=True,exist_ok=True)
    r=json.loads((run/'result.json').read_text())
    if r['status']!='completed_NPT_extension_equilibrium_not_established':raise ValueError('NPT extension not completed')
    cfg=json.loads((run/'config.json').read_text());prm=AmberPrmtopFile(str(bundle/'complex.prmtop'));system=XmlSerializer.deserialize((run/'system.xml').read_text())
    initial=XmlSerializer.deserialize((resolve_project_path(cfg['source_stage'])/'state.xml').read_text());final=XmlSerializer.deserialize((run/'state.xml').read_text())
    first=structural_snapshot(prm,initial,system);last=structural_snapshot(prm,final,system)
    with (run/'observables.csv').open() as f:rows=list(csv.DictReader(f))
    windows=window_summary(rows)
    failures=[]
    if last['PE_min_interchain_carbon_distance_nm']<.15:failures.append('Severe interchain PE carbon overlap')
    if last['max_relative_constraint_error']>1e-3:failures.append('Excessive final bond-constraint error')
    if min(last['min_water_oxygen_PE_carbon_distance_nm'],last['min_water_oxygen_protein_heavy_distance_nm'])<.12:failures.append('Severe water/solute heavy-atom overlap')
    if not .5<float(rows[-1]['total_density_g_cm3'])<1.5:failures.append('Unexpected total density requires review')
    result={'status':'preliminary_geometric_checks_passed' if not failures else 'requires_review','allow_next_preparation':not failures,'equilibrium_established':False,'production_ready':False,'initial_structure':first,'final_structure':last,'descriptive_10ps_windows':windows,'failures':failures,'interpretation':['Window means, standard deviations and slopes are descriptive, with correlated samples; no statistical equilibrium conclusion.','Rg and chain alignment can reveal gross changes but cannot validate native protein conformation or crystalline PE.','Water profile includes solute-excluded space and averages over the finite patch and its edges; it is not a bulk-water density or adsorption measure.','A successful geometric check allows more preparation, not production adsorption or free-energy claims.']}
    try:
        from .dcd_diagnostics import trajectory_checks
    except ImportError:
        from dcd_diagnostics import trajectory_checks
    result['trajectory_diagnostics']=trajectory_checks(run/'trajectory.dcd',prm,initial.getPositions(asNumpy=True).value_in_unit(unit.nanometer),out)
    # Independently cross-check the saved restart state with the last trajectory frame.
    try:
        from .dcd_diagnostics import frames
    except ImportError:
        from dcd_diagnostics import frames
    for frame_step,frame_xyz,frame_box in frames(run/'trajectory.dcd'):
        pass
    state_xyz=final.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    delta=frame_xyz-state_xyz;delta-=frame_box*np.round(delta/frame_box)
    mismatch=float(np.linalg.norm(delta,axis=1).max())
    if mismatch>1e-5:raise ValueError('Last DCD frame disagrees with saved restart coordinates')
    state_box=np.array([v.value_in_unit(unit.nanometer) for v in final.getPeriodicBoxVectors()])
    if not np.allclose(frame_box,np.diag(state_box),atol=1e-7):raise ValueError('Last DCD cell disagrees with restart cell')
    if abs(float(rows[-1]['volume_nm3'])-np.linalg.det(state_box))>1e-6:raise ValueError('CSV and restart volumes disagree')
    result['trajectory_diagnostics']['max_DCD_restart_minimum_image_difference_nm']=mismatch
    result['trajectory_diagnostics']['DCD_restart_and_CSV_box_consistent']=True

    if result['trajectory_diagnostics']['frames_read']!=100 or result['trajectory_diagnostics']['last_step']!=50000:raise ValueError('Incomplete 100 ps trajectory')
    atomic_json(out/'analysis.json',result)
    with (out/'water-profile.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['snapshot','relative_z_nm','water_oxygen_number_density_nm_minus3'])
        for label,snap in [('initial',first),('final',last)]:
            profile=snap['water_O_z_profile']
            w.writerows((label,z,d) for z,d in zip(profile['bin_center_nm'],profile['oxygen_number_density_nm_minus3']))
    print(json.dumps({'status':result['status'],'failures':failures},ensure_ascii=False),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run');p.add_argument('bundle');p.add_argument('output');a=p.parse_args();analyze(a.run,a.bundle,a.output)
