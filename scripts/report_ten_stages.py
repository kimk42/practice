"""Regenerate an evidence-based interim report without promoting pilots to production."""
import json,csv
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];RUN=ROOT/'runs';OUT=ROOT/'models/analysis'
def read(path):
    p=ROOT/path
    return json.loads(p.read_text()) if p.exists() else None

def main():
    OUT.mkdir(exist_ok=True)
    initial=read('runs/fxii-surface-pilot/result.json');surface=read('runs/coating-control-100ps/result.json');protein=read('runs/fxii-water-control-100ps/result.json');bare=read('runs/fxii-bare-pc-pilot/result.json')
    velocity=read('runs/velocity-pilots/result.json')
    stages=[
      {'step':1,'work':'initial FXII/coating numerical analysis','state':'passed numerical pilot only' if initial else 'pending FXII initialization','evidence':'runs/fxii-surface-pilot/result.json'},
      {'step':2,'work':'coating relaxation/equilibration','state':'100ps preliminary control available; equilibrium unproven' if surface else '100ps preliminary control queued','gate':'stable density/hydration/structure and sufficient sampling; temperature alone is insufficient'},
      {'step':3,'work':'surface model sensitivity','state':'0/3/6 graft coordinate candidates prepared and preprocessed; representative model unresolved','gate':'experimental DP, density, film thickness and crosslink distribution; sensitivity trajectories'},
      {'step':4,'work':'FXII protein-only control','state':'100ps preliminary control available; equilibrium unproven' if protein else '0.2ps passed;100ps queued','gate':'conformational assessment, adequate relaxation and matched conditions'},
      {'step':5,'work':'bare-PC versus coated-PC controls','state':'both short initial pilots available' if bare and initial else 'bare-PC preparation/preflight passed; matched simulations queued','gate':'matched physical models and converged sampling before efficacy comparison'},
      {'step':6,'work':'orientation/starting-state sampling','state':'18 short velocity/distance/orientation pilots passed;not production replicas' if velocity else 'two accepted orientations;bare/coated PC;0.5/2nm gap subset;three velocity seeds queued;short numerical sensitivity only','gate':'rigid orientations of one predicted conformation are not an ensemble'},
      {'step':7,'work':'replicated production MD','state':'blocked: no validated production model or converged pilot budget','gate':'relaxation/sampling criteria, independent conformations and suitable GPU/HPC resources; short velocity pilots do not satisfy this stage'},
      {'step':8,'work':'glycosylated FXII comparison','state':'nine annotated sites and five citation records retrieved;glycoform/parameters unresolved','gate':'site-specific glycan structures/occupancy, reproducible glycan building and validated compatible force-field/1-4 treatment'},
      {'step':9,'work':'full fibrinogen simulation','state':'Alpha/Gamma-A isoform sequence candidates aligned;955 construct residues unmodeled,1081 relative to this full-length candidate','gate':'select Alpha versus Alpha-E/gamma isoforms, full-length disordered-region ensembles and glycans; no arbitrary loop filling'},
      {'step':10,'work':'integrated mechanism/uncertainty report','state':'interim evidence report generated;mechanistic comparison not yet supported','gate':'matched repeated adsorption trajectories and uncertainty; no activation or antifouling conclusion from initialization'},
    ]
    metrics={}
    qa=read('models/analysis/initial-fxii-coating-QA.json')
    if qa:metrics['initial_FXII_coating_QA']=qa
    bare_qa=read('models/analysis/initial-fxii-bare-pc-QA.json')
    if bare_qa:metrics['initial_FXII_bare_PC_QA']=bare_qa
    water_qa=read('models/analysis/fxii-water-100ps-QA.json')
    if water_qa:metrics['FXII100ps_structural_QA']=water_qa
    paired=read('runs/paired-interface-followup-status.json')
    if paired:metrics['paired_interface_followup_status_snapshot']=paired
    for key,path in [('coating10ps','runs/surface-nvt-10ps/result.json'),('coating100ps','runs/coating-control-100ps/result.json'),('protein100ps','runs/fxii-water-control-100ps/result.json')]:
        value=read(path)
        if value:metrics[key]=value
    density=RUN/'surface-nvt-10ps/density-profile.xvg'
    if density.exists():
        values=np.array([[float(x) for x in line.split()] for line in density.read_text().splitlines() if line and not line.startswith(('#','@'))])
        if not np.isfinite(values).all():raise ValueError('Nonfinite density profile')
        with (OUT/'initial-coating-density.csv').open('w') as f:
            w=csv.writer(f,lineterminator="\n");w.writerow(['z_nm','coating_polymer_kg_m3','PC_kg_m3','water_kg_m3']);w.writerows(values)
        metrics['density_scope']='Short10ps descriptive profile of the provisional6x6nm cell, not equilibrated bulk density or validated brush height'
    report={'stages':stages,'completed_outputs':metrics,'all_ten_stages_completed':False,'scope':'Interim numerical/model-preparation evidence only. No adsorption free energy, coating efficacy, native glycoform or FXII activation result.'}
    (OUT/'ten-stage-status.json').write_text(json.dumps(report,indent=2)+'\n')
    text='# Ten-stage execution: interim evidence\n\nThis report follows the latest ten steps requested by the user. Prepared inputs, short numerical checks and research conclusions are distinct. No complete ten-stage claim is made.\n\n| Step | Work | Actual state |\n|---|---|---|\n'
    for s in stages:text+=f"| {s['step']} | {s['work']} | {s['state']} |\n"
    text+='\n## Measured evidence\n\nThe small coating10ps run used27826 atoms, returned21 finite energy records and mean sampled temperature297.87K. Its initial z-density profile is in models/analysis/initial-coating-density.csv. The FXII water0.2ps pilot passed, with20 annotated disulfides retained. Neither fact establishes equilibration. Expanded/protein-surface result files are accepted only when their strict preparation/force/energy gates pass.\n\n## Research gates\n\nFull-length protein modeling and a representative coating remain essential. FGA canonical P02671 is866-residue Alpha-E; the alternative Alpha isoform and Gamma-A sequences were retrieved via EBI and aligned to3GHG constructs. The provisional six-chain sequence candidate has2994 residues versus1913 observed residues, leaving1081 absent coordinates relative to that candidate (955 is the construct-relative count). Fibrinopeptides are retained; only signal peptides are removed, because the target is fibrinogen rather than fibrin. Thus using the canonical sequence to fill3GHG gaps could select the wrong isoform. The glycosylation papers establish site evidence, but the retrieved annotations do not uniquely determine all nine site-specific native glycan structures/occupancies. No placeholder glycans or missing force-field terms are introduced.\n\nThe current host has5 CPU cores and no exposed GPU. The380158-atom FXII/coating0.2ps pilot measured0.366ns/day on4 threads. This short-run benchmark corresponds to roughly2.7days per1ns per run if sustained, not a validated long-run forecast. Long, replicated sampling requires model/relaxation validation and suitable compute allocation before production lengths are selected.\n\n## Interpretation limits\n\nNo protein adsorption efficacy, adsorption free energy, physiological glycoform effect or proteolytic FXII activation has been demonstrated. Failure to observe adsorption in a short trajectory would not establish antifouling. The final mechanistic report remains gated on matched converged repeated observations.\n'
    (ROOT/'research/TEN_STAGE_PROGRESS.md').write_text(text)
if __name__=='__main__':main()
