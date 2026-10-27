"""Sparse-frame geometry/contact descriptions,not adsorption statistics or free energy."""
import json,subprocess,hashlib,csv
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def load_xvg(path):return np.array([[float(x) for x in row.split()] for row in path.read_text().splitlines() if row and not row.startswith(('#','@'))])
def main():
    reports=[];out=ROOT/'models/analysis';out.mkdir(exist_ok=True)
    for case in ['fxii-coating-near-10ps','fxii-bare-near-10ps']:
        run=ROOT/'runs'/case;result=json.loads((run/'result.json').read_text());spec=json.loads((run/'spec.json').read_text());parent=Path(spec['source_gro']).parent;pilot=json.loads((parent/'result.json').read_text());surface=Path(pilot['surface_directory']);n=json.loads((surface/'result.json').read_text())['dry_atoms']
        rows=(run/'nvt.gro').read_text().splitlines()[2:-1];groups={'Protein_heavy':[i+1 for i in range(n,n+8995) if not rows[i][10:15].strip().startswith('H')],'Surface_heavy':[i+1 for i in range(n) if not rows[i][10:15].strip().startswith('H')]}
        analysis=run/'contact-QA'
        if analysis.exists():raise RuntimeError('Existing contact analysis preserved')
        analysis.mkdir();text=''
        for name,ids in groups.items():text+='[ '+name+' ]\n'+'\n'.join(' '.join(map(str,ids[i:i+15])) for i in range(0,len(ids),15))+'\n'
        (analysis/'heavy.ndx').write_text(text)
        with (analysis/'mindist.log').open('w') as log:subprocess.run(['gmx','mindist','-s','nvt.tpr','-f','nvt.xtc','-n',str(analysis/'heavy.ndx'),'-od',str(analysis/'min-distance.xvg'),'-on',str(analysis/'contacts.xvg'),'-d','0.45','-ng','1'],cwd=run,input='0\n1\n',text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
        distance=load_xvg(analysis/'min-distance.xvg');contacts=load_xvg(analysis/'contacts.xvg')
        if not np.isfinite(distance).all() or not np.isfinite(contacts).all() or not np.array_equal(distance[:,0],contacts[:,0]):raise RuntimeError('Invalid trace')
        with (out/(case+'-sampled-contacts.csv')).open('w') as f:w=csv.writer(f,lineterminator='\n');w.writerow(['time_ps','minimum_heavy_atom_distance_nm','gmx_mindist_contact_count_cutoff_0_45nm']);w.writerows(np.column_stack([distance[:,0],distance[:,1],contacts[:,1]]))
        report={'case':case,'state':'10ps baseline geometry/energy follow-up completed','duration_ps':10,'structure_frames':len(distance),'sample_times_ps':distance[:,0].tolist(),'sampled_minimum_distance_nm':float(distance[:,1].min()),'endpoint_minimum_distance_nm':float(distance[-1,1]),'sampled_contact_counts':contacts[:,1].tolist(),'temperature_second_half_mean_K':result['temperature_second_half_mean_K'],'trajectory_sha256':hashlib.sha256((run/'nvt.xtc').read_bytes()).hexdigest(),'conditions':'legacy300K,neutralizing ions only;not cPBSzI','limits':'sparse structural frames from one geometry/seed;not contact probabilities,adsorption kinetics,equilibrium or efficacy;raw total energies of differently sized systems are not adsorption free energies'}
        reports.append(report)
    (out/'paired-interface-10ps-summary.json').write_text(json.dumps({'cases':reports,'scope':'descriptive numerical follow-ups only;no comparative adsorption conclusion'},indent=2)+'\n')
    print(json.dumps(reports))
if __name__=='__main__':main()
