"""Endpoint QA only; no adsorption or equilibration inference from0.2ps."""
import argparse,json,hashlib,re
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def gro(path):
 l=path.read_text().splitlines();return l[2:-1],np.array([[float(x[20:28]),float(x[28:36]),float(x[36:44])] for x in l[2:-1]])
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--directory',type=Path,default=ROOT/'runs/fxii-surface-pilot');parser.add_argument('--prefix',default='initial-fxii-coating');args=parser.parse_args()
 source=args.directory.resolve();result=json.loads((source/'result.json').read_text());raw,initial=gro(source/'dry.gro');_,final=gro(source/'md.gro')
 surface=Path(result['surface_directory']);surface=surface if surface.is_absolute() else ROOT/surface
 start=json.loads((surface/'result.json').read_text())['dry_atoms'];stop=start+8995
 heavy=np.array([not r[10:15].strip().startswith('H') for r in raw[start:stop]])
 a=initial[start:stop][heavy];b=final[start:stop][heavy];x=a-a.mean(0);y=b-b.mean(0);u,_,vt=np.linalg.svd(x.T@y);rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
 rmsd=float(np.sqrt(np.mean(np.sum((x@rotation-y)**2,axis=1))))
 protein=final[start:stop];substrate=final[:start];z_lower_bound=float(min(protein[:,2].min()-substrate[:,2].max(),18.-protein[:,2].max()+substrate[:,2].min()))
 template=(ROOT/'runs/fxii-unglycosylated-topology/processed.gro').read_text().splitlines()[2:-1];sg={int(r[:5]):i for i,r in enumerate(template) if r[10:15].strip()=='SG'}
 annotations=json.loads((ROOT/'models/proteins/fxii-audit.json').read_text())['disulfide_annotations'];distances=[]
 for row in annotations:
  i,j=row['residues'];distances.append({'precursor_residues':[i,j],'endpoint_SG_distance_nm':float(np.linalg.norm(protein[sg[i]]-protein[sg[j]]))})
 performance=re.search(r'Performance:\s+([\d.]+)',(source/'md.log').read_text())
 report={'source':str(source.relative_to(ROOT)),'endpoint_sha256':hashlib.sha256((source/'md.gro').read_bytes()).hexdigest(),'duration_ps':result['duration_ps'],'protein_heavy_atom_fitted_RMSD_nm':rmsd,'RMSD_reference':'placed pre-EM coordinates versus post-EM/MD endpoint; includes minimization','endpoint_vertical_separation_lower_bound_including_z_periodicity_nm':z_lower_bound,'endpoint_contact_exclusion':'if vertical bound exceeds0.45nm, protein/surface contacts below0.45nm are excluded at this endpoint; not evidence of antifouling','annotated_disulfides':distances,'short_run_performance_ns_day':float(performance.group(1)) if performance else None,'scope':'Endpoint numerical/geometry QA only; not full trajectory convergence, adsorption, native structure validation or activation'}
 out=ROOT/'models/analysis';out.mkdir(exist_ok=True);(out/(args.prefix+'-QA.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='annotated_disulfides'}))
if __name__=='__main__':main()
