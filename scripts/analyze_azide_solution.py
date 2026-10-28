"""Descriptive hydration/ion-distance diagnostics;no equilibrium or accuracy approval."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/azide-rigid-water-102ps'


def table(path):return np.array([[float(v) for v in line.split()] for line in path.read_text().splitlines() if line.strip() and not line.startswith(('#','@'))])


def main():
    out=RUN/'solution-analysis'
    if out.exists():raise RuntimeError('Existing analysis preserved')
    assert 'Finished mdrun' in (RUN/'nvt.part0002.log').read_text()
    assert 'LINCS WARNING' not in (RUN/'nvt.part0002.log').read_text()
    out.mkdir()
    groups={'TERMINAL_N':[3,5],'CENTRAL_N':[4],'SODIUM':[6],'WATER_O':list(range(7,2632,3)),'AZIDE':[1,2,3,4,5]}
    (out/'groups.ndx').write_text(''.join('[ '+key+' ]\n'+' '.join(map(str,value))+'\n' for key,value in groups.items()))
    def run(args,name,stdin=None):
        with (out/(name+'.log')).open('w') as stream:
            subprocess.run(['gmx',*args],cwd=out,input=stdin,text=True,stdout=stream,stderr=subprocess.STDOUT,check=True)
    run(['energy','-f',str(RUN/'nvt.part0002.edr'),'-o','energy.xvg'],'energy','Potential\nTemperature\n0\n')
    energy=table(out/'energy.xvg');assert energy[0,0]==2 and energy[-1,0]==102 and len(energy)==2001 and np.isfinite(energy).all()
    traj=['-s',str(RUN/'nvt.tpr'),'-f',str(RUN/'nvt.part0002.trr'),'-n','groups.ndx']
    rdfs={}
    for reference in ['TERMINAL_N','CENTRAL_N']:
        for label,start,end in [('first',2,52),('second',52,102)]:
            stem=reference.lower()+'-'+label
            run(['rdf',*traj,'-ref',f'group "{reference}"','-sel','group "WATER_O"',
                 '-b',str(start),'-e',str(end),'-bin','0.01','-rmax','1.2','-o',stem+'.xvg','-cn',stem+'-cn.xvg'],stem)
            data=table(out/(stem+'.xvg'));cn=table(out/(stem+'-cn.xvg'))
            assert np.isfinite(data).all() and np.isfinite(cn).all()
            mask=(data[:,0]>=.2)&(data[:,0]<=.45);peak=data[mask][np.argmax(data[mask,1])]
            rdfs[stem]={'time_window_ps':[start,end],'peak_search_interval_nm':[.2,.45],
                        'peak_r_nm':float(peak[0]),'peak_g_r':float(peak[1]),
                        'coordination_at_0_35nm':float(np.interp(.35,cn[:,0],cn[:,1]))}
    run(['mindist',*traj,'-od','terminal-na-distance.xvg','-on','terminal-na-contact.xvg','-d','0.35'],
        'terminal-na','TERMINAL_N\nSODIUM\n')
    distance=table(out/'terminal-na-distance.xvg');assert len(distance)==len(energy) and np.isfinite(distance).all()
    run(['trjconv',*traj,'-o','sites.g96','-ndec','9'],'geometry','AZIDE\n')
    lines=(out/'sites.g96').read_text().splitlines()
    coords=np.array([[[float(v) for v in row.split()] for row in lines[i+1:i+6]]
                     for i,line in enumerate(lines) if line=='POSITIONRED'])
    assert len(coords)==2001
    disp=coords-coords[:,:1];disp-=3*np.rint(disp/3);coords=coords[:,:1]+disp
    left=coords[:,2]-coords[:,3];right=coords[:,4]-coords[:,3]
    lengths=np.stack([np.linalg.norm(left,axis=1),np.linalg.norm(right,axis=1)],axis=1)
    bond_error=float(np.abs(lengths-.11626).max())
    collinear=float(np.linalg.norm(left+right,axis=1).max())
    assert max(bond_error,collinear)<5e-6
    result={'state':'checkpoint continuation and descriptive solution diagnostics completed',
            'duration_total_ps':102,'continuation_duration_ps':100,'sampled_frames':2001,
            'temperature_mean_K':float(energy[:,2].mean()),
            'temperature_second_half_mean_K':float(energy[energy[:,0]>=52,2].mean()),
            'maximum_NN_error_nm':bond_error,'maximum_terminal_symmetry_error_nm':collinear,
            'water_RDF_windows':rdfs,
            'terminal_N_to_Na_min_distance_range_nm':[float(distance[:,1].min()),float(distance[:,1].max())],
            'fraction_saved_frames_below_0_35nm':float((distance[:,1]<.35).mean()),
            'contact_threshold_scope':'0.35nm is a screening definition;not a validated ion-pair boundary or equilibrium probability',
            'concentration_scope':'one ion pair in27nm3,about61.5mM each;not experimental3.08mM azide or full cPBSzI',
            'sampling_limits':'one ion pair/trajectory;two50ps windows share their boundary frame;RDFs are descriptive,not independent replicates or a convergence proof',
            'physical_force_field_approved':False,'full_cpbszi_ready':False,
            'source_sha256':{name:hashlib.sha256((RUN/name).read_bytes()).hexdigest() for name in ['nvt.tpr','nvt.cpt','nvt.part0002.trr','source.cpt']}}
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'models/buffer-parameters/azide-solution-102ps-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    for name in rdfs:
        data=table(out/(name+'.xvg'))
        np.savetxt(ROOT/f'models/analysis/azide-{name}-rdf.csv',data,delimiter=',',header='r_nm,g_r',comments='')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
