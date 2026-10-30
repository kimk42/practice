"""Continue accepted dry PC controls from20ps to100ps, preserving all state.

This is a bounded relaxation diagnostic, not physical substrate approval.
One worker per DP; lock descriptors remain inherited by GROMACS children.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

from pc_length_parameters import ROOT, sha
from run_pc_sqm_checkpointed import write_json

BASE=ROOT/'runs/pc-dry-npt100-v1'


def run(command,out,label,lock_fd,stdin=None):
    with (out/(label+'.log')).open('a') as stream:
        subprocess.run(['gmx',*command],cwd=out,input=stdin,text=True,
                       stdout=stream,stderr=subprocess.STDOUT,check=True,
                       pass_fds=(lock_fd,))


def check_tpr_comparison(text):
    lines=[s.strip() for s in text.split('comparing inputrec',1)[-1].splitlines() if s.strip()]
    differences=[s for s in lines if not s.startswith('comparing ')]
    if 'comparing inputrec' not in text or differences!=['inputrec->nsteps (40000 - 200000)']:
        raise ValueError('TPR comparison is not limited to the requested step extension')


def resume_command(out):
    command=['mdrun','-s','npt.tpr','-deffnm','npt','-ntmpi','1','-ntomp','1','-cpt','1','-maxh','1']
    if (out/'npt.cpt').exists():
        return command+['-cpi','npt.cpt','-append']
    if list(out.glob('npt.part*')) or (out/'npt.log').exists():
        raise RuntimeError('Partial MD without checkpoint;preserve and review instead of resetting')
    if not (out/'source.cpt').exists():
        raise ValueError('Source checkpoint missing')
    return command+['-cpi','source.cpt','-noappend']


def table(path):
    return np.array([[float(v) for v in s.split()] for s in path.read_text().splitlines()
                     if s and not s.startswith(('#','@'))])


def analyze(out,source,lock_fd):
    energies=list(out.glob('npt.part*.edr'));logs=list(out.glob('npt.part*.log'));trajectories=list(out.glob('npt.part*.xtc'))
    if not (len(energies)==len(logs)==len(trajectories)==1):
        raise ValueError('Unexpected continuation parts require review')
    text=logs[0].read_text()
    if 'Finished mdrun' not in text or re.search(r'LINCS WARNING|constraint.*warning|Fatal error',text,re.I):
        raise ValueError('Run completion/numerical-warning check failed')
    run(['energy','-f',energies[0].name,'-o','energy.xvg'],out,'energy-extract',lock_fd,
        'Potential\nTemperature\nPressure\nVolume\nDensity\n0\n')
    data=table(out/'energy.xvg')
    if data.shape!=(1601,6) or not np.isfinite(data).all() or not np.allclose(data[[0,-1],0],[20.,100.]):
        raise ValueError('Continuation missing finite frames or correct20-100ps interval')
    if np.any(data[:,[2,4,5]]<=0):
        raise ValueError('Nonpositive temperature/volume/density')
    run(['check','-f',trajectories[0].name],out,'trajectory-check',lock_fd)
    if not re.search(r'^Coords\s+801\s+0\.1',(out/'trajectory-check.log').read_text(),re.M):
        raise ValueError('Unexpected trajectory cadence/frame count')
    windows=[]
    for start in (20,40,60,80):
        block=data[(data[:,0]>=start)&((data[:,0]<start+20) if start<80 else (data[:,0]<=100))]
        windows.append({'start_ps':start,'end_ps':start+20,'upper_endpoint_included':start==80,
                        'frames':len(block),'temperature_mean_K':float(block[:,2].mean()),
                        'pressure_mean_bar':float(block[:,3].mean()),
                        'density_mean_kg_m3':float(block[:,5].mean()),
                        'density_std_kg_m3':float(block[:,5].std())})
    density_change=100*(windows[-1]['density_mean_kg_m3']/windows[-2]['density_mean_kg_m3']-1)
    late=data[data[:,0]>=80]
    return {'case':out.name,'state':'20-to100ps NPT continuation and density diagnostics completed',
            'continuation_duration_ps':80.,'total_NPT_ps':100.,'finite_energy_frames':len(data),
            'coordinate_frames':801,'temperature_mean_K':float(data[:,2].mean()),
            'density_final_kg_m3':float(data[-1,5]),'density_windows':windows,
            'last_two_windows_density_change_percent':density_change,
            'last20ps_density_slope_kg_m3_ps':float(np.polyfit(late[:,0],late[:,5],1)[0]),
            'source':str(source.relative_to(ROOT)),
            'outputs_sha256':{p.name:sha(p) for p in [out/'npt.tpr',out/'npt.cpt',out/'source.cpt',energies[0],trajectories[0]]},
            'physical_approval':False,'equilibrium_claim':False,
            'scope':'unchanged300K/1bar finite OH-ended dry oligomer diagnostic;no coating,FXII,cPBSzI',
            'limitations':'short correlated windows and one common intrachain conformer;PC4 and PC8 initial packing protocols differ;neither molecular weight nor equilibrium density is experimentally established'}


def execute(case,lock_fd):
    source=ROOT/'runs/pc-bulk-controls-v2'/case
    previous=json.loads((source/'result.json').read_text())
    if previous['state']!='EM/2ps NVT/20ps NPT numerical control completed':
        raise ValueError('Source control is not complete')
    for name,digest in previous['checkpoints_sha256'].items():
        if sha(source/name)!=digest:raise ValueError('Source checkpoint changed')
    dp=previous['PC_DP'];parameter=json.loads((ROOT/f'runs/pc-length-forcefield-v2/pc{dp}/parameter-audit.json').read_text())
    if sha(source/'pc.itp')!=parameter['outputs_sha256'][f'PC{dp}_GMX.itp']:
        raise ValueError('Source forcefield differs from the audited candidate')
    spec={'source_sha256':{name:sha(source/name) for name in ['npt.tpr','npt.cpt','pc.itp','topol.top','result.json']},
          'original_NPT_end_ps':20.,'new_NPT_end_ps':100.,'extension_ps':80.,
          'protocol_change':'nsteps only;retain checkpoint velocities,thermostat/barostat state and original forcefield',
          'reason':'observed density drift in the initial20ps;bounded follow-up before substrate selection',
          'physical_approval':False}
    out=BASE/case;out.mkdir(exist_ok=True)
    if (out/'spec.json').exists() and json.loads((out/'spec.json').read_text())!=spec:
        raise ValueError('Continuation source/spec changed')
    write_json(out/'spec.json',spec)
    if (out/'completed.json').exists():
        completed=json.loads((out/'completed.json').read_text())
        if any(sha(out/name)!=digest for name,digest in completed['outputs_sha256'].items()):
            raise ValueError('Completed continuation evidence changed')
        return completed
    if not (out/'npt.tpr').exists():
        shutil.copy2(source/'npt.cpt',out/'source.cpt')
        run(['convert-tpr','-s',str(source/'npt.tpr'),'-until','100','-o','npt.tpr'],out,'convert',lock_fd)
        run(['check','-s1',str(source/'npt.tpr'),'-s2','npt.tpr'],out,'tpr-comparison',lock_fd)
        check_tpr_comparison((out/'tpr-comparison.log').read_text())
        write_json(out/'tpr-audit.json',{'sha256':sha(out/'npt.tpr'),'only_nsteps_changed':True})
    tpr_audit=json.loads((out/'tpr-audit.json').read_text())
    if sha(out/'npt.tpr')!=tpr_audit['sha256'] or sha(out/'source.cpt')!=spec['source_sha256']['npt.cpt']:
        raise ValueError('Continuation TPR/source state changed')
    write_json(BASE/f'pc{dp}-status.json',{'state':'running','case':case,'worker_pid':os.getpid(),
                                         'updated_utc':datetime.now(timezone.utc).isoformat()})
    # A final GRO indicates a clean mdrun exit, which can also be a -maxh stop.
    # Its energy endpoint decides whether another checkpoint continuation is needed.
    finals=list(out.glob('npt.part*.gro'))
    done=False
    if finals:
        energies=list(out.glob('npt.part*.edr'))
        if len(energies)!=1:raise ValueError('Unexpected energy parts')
        run(['energy','-f',energies[0].name,'-o','endpoint.xvg'],out,'endpoint',lock_fd,'Temperature\n0\n')
        done=abs(table(out/'endpoint.xvg')[-1,0]-100)<1e-6
    if not done:run(resume_command(out),out,'launch',lock_fd)
    result=analyze(out,source,lock_fd)
    if any(sha(source/name)!=digest for name,digest in spec['source_sha256'].items()):
        raise ValueError('Original source changed during continuation')
    write_json(out/'completed.json',result)
    write_json(ROOT/f'models/surface-revision-v2/{case}-npt100-audit.json',result)
    return result


def main(dp):
    BASE.mkdir(exist_ok=True)
    with (BASE/f'pc{dp}.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        status=BASE/f'pc{dp}-status.json'
        if status.exists() and json.loads(status.read_text()).get('state')=='held_preserved':
            raise RuntimeError('Prior failure is held;review preserved evidence before any retry')
        try:
            for seed in (2026,2027):
                result=execute(f'pc{dp}-rho0.6-seed{seed}',lock.fileno())
                print(json.dumps({'case':result['case'],'state':result['state']}),flush=True)
            write_json(BASE/f'pc{dp}-status.json',{'state':'completed','cases':2,'physical_approval':False})
        except Exception as error:
            write_json(BASE/f'pc{dp}-status.json',{'state':'held_preserved','reason':str(error)})
            raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dp',type=int,choices=[4,8],required=True)
    main(parser.parse_args().dp)
