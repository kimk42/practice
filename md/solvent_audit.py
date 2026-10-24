"""Single-frame periodic solvent-region/void diagnostics, not equilibrium validation."""
import argparse
import hashlib
import html
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from openmm import XmlSerializer,unit
from openmm.app import AmberPrmtopFile,element
from workbench.config import ROOT,MODELS
from workbench.manager import atomic


def region_statistics(xyz,lengths,solute_heavy,water_O,water_mass_da,samples=100000,seed=1729,cutoff_nm=.8):
 xyz=np.asarray(xyz);lengths=np.asarray(lengths);tree=cKDTree(np.mod(xyz[solute_heavy],lengths),boxsize=lengths);water=np.mod(xyz[water_O],lengths);fraction_water=tree.query(water)[0]>=cutoff_nm;points=np.random.default_rng(seed).random((samples,3))*lengths;eligible=tree.query(points)[0]>=cutoff_nm;f=float(eligible.mean())
 if not 0<f<1:raise ValueError('No usable mixed solvent/solute exclusion domain')
 volume=float(np.prod(lengths)*f);se=float(np.prod(lengths)*np.sqrt(f*(1-f)/samples));distance=cKDTree(water,boxsize=lengths).query(points[eligible])[0];density=int(fraction_water.sum())*water_mass_da*.00166053906660/volume
 return {'samples':samples,'seed':seed,'solute_heavy_exclusion_nm':cutoff_nm,'solvent_region_volume_nm3':volume,'volume_MC_standard_error_nm3':se,'water_O_in_region':int(fraction_water.sum()),'water_molecular_mass_Da':water_mass_da,'solvent_region_water_density_g_cm3':float(density),'density_geometry_MC_standard_error_g_cm3':float(density*se/volume),'region_volume_fraction_with_nearest_water_O_over_0_4nm':float(np.mean(distance>.4)),'region_volume_fraction_with_nearest_water_O_over_0_6nm':float(np.mean(distance>.6)),'interpretation':'One frame. Monte Carlo errors describe geometric volume integration only, not MD sampling or thermodynamic uncertainty. Low-water-proximity regions are geometric diagnostics, not identified vapor phases.'}


def audit(job,output,root=ROOT):
 job=Path(job);output=Path(output);output.mkdir(parents=True,exist_ok=True);cfg=json.loads((job/'config.json').read_text());model=MODELS[cfg['model_id']];statefile=job/'final-state.xml';state=XmlSerializer.deserialize(statefile.read_text());xyz=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer);box=np.array([v.value_in_unit(unit.nanometer) for v in state.getPeriodicBoxVectors()]);lengths=np.diag(box)
 if not np.allclose(box,np.diag(lengths)):raise ValueError('Orthorhombic state required')
 prm=AmberPrmtopFile(str(Path(root)/model['prmtop']));atoms=list(prm.topology.atoms());waters=[r for r in prm.topology.residues() if r.name in ('HOH','WAT')];water_O=np.array([a.index for a in atoms if a.residue.name in ('HOH','WAT') and a.element==element.oxygen]);solute=np.array([a.index for a in atoms if a.residue.name not in ('HOH','WAT') and a.element!=element.hydrogen]);mass=np.array(prm._prmtop.getMasses());water_mass=float(sum(mass[a.index] for a in waters[0].atoms()));stats=region_statistics(xyz,lengths,solute,water_O,water_mass)
 h=hashlib.sha256(statefile.read_bytes()).hexdigest();result={'status':'single_frame_solvent_region_geometry_audited','source_job':job.name,'source_final_state_sha256':h,'context_time_ps':state.getTime().value_in_unit(unit.picosecond),'whole_mixture_density_g_cm3':float(mass.sum()*.00166053906660/np.prod(lengths)),'whole_cell_nm3':float(np.prod(lengths)),'water_molecules':len(water_O),'domain_definition':'Points at least 0.8 nm from every non-water heavy atom, including ions; minimum-image periodic Euclidean distances. Water counted using its oxygen in the same domain.','statistics':stats,'equilibrium_established':False,'limitations':['Not a calibrated liquid-density or surface-validation test.','Single frame may include cavities and strong spatial heterogeneity; this is not necessarily thermodynamic bulk water.','No adsorption, phase or retention inference.','Whole-mixture density includes protein, PE and ions and must not be compared directly with pure-water density.']};atomic(output/'solvent-audit.json',result)
 body='<h1>PE 준비 계산의 용매 공간 점검</h1><p>누적 '+str(round(result['context_time_ps'],3))+' ps의 저장된 한 프레임을 검사했습니다. 평형 검증이 아닙니다.</p><p>전체 혼합계 밀도와, 용질에서 떨어진 공간의 물 밀도를 구분합니다. 물이 드문 공간의 비율도 기하 지표로 기록합니다. 공간 적분의 Monte Carlo 오차를 MD의 통계 오차로 해석하지 않습니다.</p><pre>'+html.escape(json.dumps(result,ensure_ascii=False,indent=2))+'</pre>'
 (output/'report.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><title>PE 용매 공간 진단</title><style>body{max-width:1000px;margin:auto;padding:24px;font:16px/1.7 system-ui}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style><body>'+body+'</body></html>');return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--job',required=True);p.add_argument('--output',required=True);a=p.parse_args();print(json.dumps(audit(ROOT/'runs/workbench'/a.job,a.output),indent=2))
