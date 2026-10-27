"""Structural/numerical QA for the completed100ps water control,not native validation."""
import json,subprocess,hashlib,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from surface_validation import section_rows
def coords(path):
    rows=path.read_text().splitlines()[2:2+8995]
    return np.array([[float(r[20:28]),float(r[28:36]),float(r[36:44])] for r in rows])
def main():
    source=ROOT/'runs/fxii-water-control-100ps';result=json.loads((source/'result.json').read_text())
    if result['state']!='bounded numerical follow-up completed':raise RuntimeError('Control incomplete')
    out=source/'structural-QA'
    if out.exists():raise RuntimeError('Existing QA preserved')
    out.mkdir()
    def run(args,name,stdin='0\n'):
        with (out/(name+'.log')).open('w') as log:subprocess.run(['gmx',*args],cwd=source,input=stdin,text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
    run(['trjconv','-s','nvt.tpr','-f','nvt.xtc','-o',str(out/'whole.xtc'),'-pbc','mol'],'whole-trajectory')
    for name in ['source','nvt']:run(['trjconv','-s','nvt.tpr','-f',name+'.gro','-o',str(out/(name+'-whole.gro')),'-pbc','mol'],name+'-whole')
    run(['rms','-s','nvt.tpr','-f',str(out/'whole.xtc'),'-o',str(out/'backbone-rmsd.xvg'),'-tu','ps'],'rmsd','Backbone\nBackbone\n')
    v=np.array([[float(x) for x in r.split()] for r in (out/'backbone-rmsd.xvg').read_text().splitlines() if r and not r.startswith(('#','@'))])
    if not np.isfinite(v).all() or v[-1,0]-v[0,0]<99.99:raise RuntimeError('Incomplete/nonfinite structural trace')
    mass=np.array([float(r[7]) for r in section_rows((source/'topol.top').read_text(),'atoms')]);assert len(mass)==8995
    initial=coords(out/'source-whole.gro');final=coords(out/'nvt-whole.gro')
    def rg(x):return float(np.sqrt(np.average(np.sum((x-np.average(x,axis=0,weights=mass))**2,axis=1),weights=mass)))
    template=(ROOT/'runs/fxii-unglycosylated-topology/processed.gro').read_text().splitlines()[2:-1];sg={int(r[:5]):i for i,r in enumerate(template) if r[10:15].strip()=='SG'}
    links=[]
    for row in json.loads((ROOT/'models/proteins/fxii-audit.json').read_text())['disulfide_annotations']:
        a,b=row['residues'];distance=float(np.linalg.norm(final[sg[a]]-final[sg[b]]));links.append({'residues':[a,b],'distance_nm':distance})
    flagged=[x for x in links if not .15<=x['distance_nm']<=.30]
    report={'state':'100ps numerical structural QA completed','duration_ps':100,'trajectory_sha256':hashlib.sha256((source/'nvt.xtc').read_bytes()).hexdigest(),'backbone_RMSD_endpoint_nm':float(v[-1,1]),'backbone_RMSD_max_nm':float(v[:,1].max()),'radius_gyration_initial_nm':rg(initial),'radius_gyration_endpoint_nm':rg(final),'annotated_disulfides':links,'disulfide_geometry_flags':flagged,'scope':'single predicted unglycosylated conformation;these short structural metrics are not equilibrium or native-conformation validation'}
    target=ROOT/'models/analysis';(target/'fxii-water-100ps-QA.json').write_text(json.dumps(report,indent=2)+'\n')
    (target/'fxii-water-100ps-backbone-rmsd.xvg').write_bytes((out/'backbone-rmsd.xvg').read_bytes())
    if flagged:raise RuntimeError('Disulfide geometry review required')
    print(json.dumps({k:v for k,v in report.items() if k!='annotated_disulfides'}))
if __name__=='__main__':main()
