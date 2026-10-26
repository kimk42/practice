"""Single-frame protein periodic-image clearance; short-reference reconstruction only."""
import json,itertools,argparse
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from openmm import XmlSerializer,unit
from openmm.app import AmberPrmtopFile,element
from workbench.config import ROOT,MODELS
from workbench.manager import atomic
parser=argparse.ArgumentParser();parser.add_argument('--job',required=True);args=parser.parse_args()
if len(args.job)!=32 or any(c not in '0123456789abcdef' for c in args.job): raise ValueError('Invalid job id')
job=ROOT/'runs/workbench'/args.job
s=XmlSerializer.deserialize((job/'final-state.xml').read_text());xyz=s.getPositions(asNumpy=True).value_in_unit(unit.nanometer);box=np.diag([v.value_in_unit(unit.nanometer) for v in s.getPeriodicBoxVectors()]);ref=np.load(job/'reference.npy')
m=MODELS['pe-core'];src=XmlSerializer.deserialize((ROOT/m['source']/'state.xml').read_text());b0=np.diag([v.value_in_unit(unit.nanometer) for v in src.getPeriodicBoxVectors()]);pred=ref*(box/b0);delta=xyz-pred;delta-=box*np.round(delta/box);unwrapped=pred+delta
prm=AmberPrmtopFile(str(ROOT/m['prmtop']));atoms=list(prm.topology.atoms());maps=json.loads((ROOT/'data/fibrinogen-core03/manifest.json').read_text())['residue_mapping'];ph=[a.index for a in atoms if a.residue.index<len(maps) and maps[a.residue.index]['source_label_seq'] is not None and a.element!=element.hydrogen];p=unwrapped[ph];tree=cKDTree(p)
image_dist=[]
for shift in itertools.product([-1,0,1],repeat=3):
 if shift==(0,0,0):continue
 image_dist.append((float(tree.query(p+np.array(shift)*box)[0].min()),shift))
minimum,shift=min(image_dist)
bond_errors=[]
for a,b in prm.topology.bonds():
 if a.index<len(ref) and b.index<len(ref):
  before=ref[a.index]-ref[b.index];before-=b0*np.round(before/b0)
  after=unwrapped[a.index]-unwrapped[b.index]
  if np.linalg.norm(after)>.4:bond_errors.append([a.index,b.index,float(np.linalg.norm(after))])
r={'source_job':job.name,'context_time_ps':float(s.getTime().value_in_unit(unit.picosecond)),'box_lengths_nm':box.tolist(),'protein_extent_nm':np.ptp(p,axis=0).tolist(),'protein_nearest_periodic_self_image_heavy_nm':minimum,'nearest_image_shift':shift,'short_reference_reconstruction_bonds_over_0_4nm':bond_errors[:20],'long_reconstructed_bond_count':len(bond_errors),'minimum_image_clearance_gate_nm':2.0,'geometry_gate_passed':minimum>2 and not bond_errors,'interpretation':'Single final-frame finite protein self-image separation using short-preparation reference reconstruction; not general unwrapping, phase validation or equilibrium.'}
out=job/'analysis/boundary-audit.json';atomic(out,r);print(json.dumps(r,indent=2))
