"""Geometric research diagnostics, never an estimator of affinity or retention."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from openmm import XmlSerializer,unit
from openmm.app import AmberPrmtopFile,element
from .dcd_diagnostics import frames
from workbench.config import ROOT,MODELS
from workbench.manager import atomic


def diagnostics(xyz,box,reference,reference_box,ca,protein_heavy,polymer_heavy,water_O,residue_labels):
 xyz=np.asarray(xyz);box=np.asarray(box);reference=np.asarray(reference)
 if box.shape!=(3,) or np.any(box<=0) or not np.isfinite(xyz).all():raise ValueError('Finite coordinates and positive orthorhombic box required')
 # Recover a finite patch near the scaled initial image. This is not valid for
 # long diffusion/rearrangement beyond the reference image and must not be
 # mistaken for general long-trajectory molecule unwrapping.
 predicted=reference*(box/np.asarray(reference_box));delta=xyz-predicted;delta-=box*np.round(delta/box);unwrapped=predicted+delta
 p=unwrapped[np.asarray(polymer_heavy)];_,singular,axes=np.linalg.svd(p-p.mean(axis=0),full_matrices=False);normal=axes[-1]
 if normal[2]<0:normal=-normal
 separation=float(singular[-1]/max(singular[-2],1e-12));pca_valid=bool(singular[-2]>1e-8 and separation<.5)
 ca_xyz=unwrapped[np.asarray(ca)];_,ca_s,ca_axes=np.linalg.svd(ca_xyz-ca_xyz.mean(axis=0),full_matrices=False);axis=ca_axes[0]
 protein_axis_valid=bool(ca_s[0]>1e-8 and ca_s[1]/ca_s[0]<.5)
 tree=cKDTree(np.mod(xyz[np.asarray(polymer_heavy)],box),boxsize=box);d=tree.query(np.mod(xyz[np.asarray(protein_heavy)],box))[0];contacts=sorted({str(residue_labels[i]) for i in np.flatnonzero(d<.45)});waters=np.asarray(water_O,dtype=int);bridging=0
 if len(waters):
  w=np.mod(xyz[waters],box);dp=cKDTree(np.mod(xyz[np.asarray(protein_heavy)],box),boxsize=box).query(w)[0];ds=tree.query(w)[0];bridging=int(np.count_nonzero((dp<.35)&(ds<.35)))
 return {'minimum_protein_polymer_heavy_nm':float(d.min()),'contact_residues':contacts,'contact_residue_count':len(contacts),'protein_axis_angle_to_fitted_patch_plane_deg':float(np.rad2deg(np.arcsin(np.clip(abs(np.dot(axis,normal)),0,1)))) if pca_valid and protein_axis_valid else None,'patch_PCA_normal':normal.tolist(),'patch_plane_fit_RMS_height_nm':float(singular[-1]/np.sqrt(len(p))),'patch_plane_PCA_separation_ratio':separation,'patch_plane_fit_usable':pca_valid,'protein_axis_fit_usable':protein_axis_valid,'water_O_within_0_35nm_of_both_groups':bridging,'interpretation':'Distance contacts and double-proximity water count only, not bonds, binding sites, affinity, residence/off-rate or experimental retention. Finite-patch PCA and reference-image reconstruction are short-preparation diagnostics only.'}


def analyze(job,output,root=ROOT):
 job=Path(job);root=Path(root);output=Path(output);output.mkdir(parents=True,exist_ok=True);cfg=json.loads((job/'config.json').read_text());m=MODELS[cfg['model_id']];prm=AmberPrmtopFile(str(root/m['prmtop']));atoms=list(prm.topology.atoms());source=XmlSerializer.deserialize((root/m['source']/'state.xml').read_text());reference_box=np.diag([v.value_in_unit(unit.nanometer) for v in source.getPeriodicBoxVectors()]);reference=np.load(job/'reference.npy');core=cfg['model_id']=='pe-core'
 maps=json.loads((root/'data/fibrinogen-core03/manifest.json').read_text())['residue_mapping'] if core else None
 native=lambda a:(a.residue.index<len(maps) and maps[a.residue.index]['source_label_seq'] is not None) if core else a.residue.index<252
 ph=np.array([a.index for a in atoms if native(a) and a.element!=element.hydrogen]);ca=np.array([a.index for a in atoms if native(a) and a.name=='CA']);poly=np.array([a.index for a in atoms if a.residue.name=='PEO' and a.element!=element.hydrogen]);water=np.array([a.index for a in atoms if a.residue.name in ('HOH','WAT') and a.element==element.oxygen]);labels=[str(maps[atoms[i].residue.index]['source_chain'])+':'+str(maps[atoms[i].residue.index]['source_label_seq']) if core else 'fragment:'+str(atoms[i].residue.id) for i in ph]
 by_step={r['context_step']:r for r in map(json.loads,(job/'observations.jsonl').read_text().splitlines())};rows=[];hashes={}
 for segment in json.loads((job/'trajectory-index.json').read_text()):
  path=job/segment['path']
  if not path.exists() or path.stat().st_size==0:continue
  h=hashlib.sha256()
  with path.open('rb') as f:
   while block:=f.read(8*1024**2):h.update(block)
  hashes[segment['path']]=h.hexdigest()
  for step,xyz,box in frames(path):
   context=step+segment['context_step_offset'];obs=by_step[context];r=diagnostics(xyz,box,reference,reference_box,ca,ph,poly,water,labels);r.update(context_step=context,context_time_ps=obs['context_time_ps'],additional_time_ps=obs['additional_time_ps'],stage=obs['stage'],restraint_scale=obs['restraint_scale']);rows.append(r)
 if not rows:raise ValueError('No saved trajectory frames available')
 with (output/'geometric-observables.csv').open('w',newline='') as f:
  keys=['context_step','context_time_ps','additional_time_ps','stage','restraint_scale','minimum_protein_polymer_heavy_nm','contact_residue_count','protein_axis_angle_to_fitted_patch_plane_deg','patch_plane_fit_RMS_height_nm','patch_plane_fit_usable','water_O_within_0_35nm_of_both_groups'];writer=csv.DictWriter(f,fieldnames=keys,lineterminator='\n');writer.writeheader();writer.writerows({k:r[k] for k in keys} for r in rows)
 summary={'status':'executed_preparation_trajectory_geometry_analyzed','job_id':job.name,'model_id':cfg['model_id'],'frames':len(rows),'saved_frame_context_time_range_ps':[rows[0]['context_time_ps'],rows[-1]['context_time_ps']],'contacts_ever_observed':sorted(set(x for r in rows for x in r['contact_residues'])),'last_frame':rows[-1],'trajectory_sha256':hashes,'geometry_rows':rows,'source_lineage_note':'Only this job’s saved DCD segments are analyzed. Parent trajectories are preserved separately; carried parent observations are not counted as new frames.','limitations':['Preparation trajectory, not production or an adsorption experiment.','Small-model fragment:* labels are internal prepared-residue numbers, not full-protein site assignments.','Contact threshold 0.45 nm; water double-proximity threshold 0.35 nm to each group; no hydrogen-bond angular criterion.','No contact lifetime, off-rate, adsorption free energy or material ranking inferred.','Plane-normal fitting uses short-preparation reference-image reconstruction; not a general long-time unwrapping method.']};atomic(output/'geometry-analysis.json',summary);return summary

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--job',required=True);p.add_argument('--output',required=True);a=p.parse_args();print(json.dumps({k:v for k,v in analyze(ROOT/'runs/workbench'/a.job,a.output).items() if k!='geometry_rows'},indent=2))
