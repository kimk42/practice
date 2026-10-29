"""Full-precision coordinate checkpoints using SQM AM1 analytic gradients.

SciPy L-BFGS-B controls geometry; every energy/gradient comes from installed SQM.
Accepted coordinates are saved as double-precision JSON, with raw SQM evidence.
After an interruption coordinates are retained, while optimizer memory restarts.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from scipy.optimize import minimize
import scipy

from pc_length_parameters import ROOT, preflight, sha
from run_pc_sqm_checkpointed import BASE as NATIVE, amber_environment, read_initial, write_json

BASE=ROOT/'runs/pc-length-sqm-lbfgs-v3'
TARGET=.0005
NUMBER=r'[-+]?(?:\d+\.\d*|\.\d+)(?:[eEdD][-+]?\d+)?'


def validate_coordinates(path, xyz, elements):
    """Match SQM's rounded final table to the exact evaluated coordinates."""
    text = path.read_text().rsplit('Final Structure', 1)[-1]
    rows = re.findall(r'^\s*QMMM:\s+(\d+)\s+(\d+)\s+([A-Z][a-z]?)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', text, re.M)
    symbols = {1:'H',6:'C',8:'O',15:'P'}
    if [int(r[0]) for r in rows] != list(range(1,len(elements)+1)) or [r[2] for r in rows] != [symbols[z] for z in elements]:
        raise ValueError('SQM final atom order/elements differ')
    rounded = np.array([[float(v) for v in r[3:6]] for r in rows])
    if not np.allclose(rounded,xyz,atol=.000051,rtol=0):
        raise ValueError('SQM coordinates differ from accepted input beyond output rounding')


def parse_singlepoint(path,natoms):
    text=path.read_text()
    if 'Calculation Completed' not in text or 'SCF Converged' not in text:
        raise ValueError('SQM single-point completion/SCF convergence missing')
    if re.search(r'SCF.*(?:failed|not converg)|FATAL|XMIN ERROR',text,re.I):
        raise ValueError('SQM error in preserved output')
    energies=re.findall(r'QMMM: SCF Energy =\s*('+NUMBER+r')\s+KCal/mol',text)
    if not energies:raise ValueError('Missing full-precision energy')
    block=text.rsplit('Forces on QM atoms from SCF calculation',1)[-1]
    rows=re.findall(r'^QMMM: Atm\s+(\d+):\s*(.*)$',block,re.M)
    if [int(r[0]) for r in rows]!=list(range(1,natoms+1)):
        raise ValueError('Missing or unordered final analytic gradient')
    values=[re.findall(NUMBER,r[1]) for r in rows]
    if any('*' in row[1] or len(v)!=3 for row,v in zip(rows,values)):
        raise ValueError('SQM printed gradient overflow or incomplete vector; trial is not accepted')
    gradient=np.array([[float(v.replace('D','E').replace('d','e')) for v in row] for row in values])
    energy=float(energies[-1]);rms=float(np.sqrt(np.mean(gradient**2)))
    if gradient.shape!=(natoms,3) or not np.isfinite(gradient).all() or not np.isfinite(energy):
        raise ValueError('Nonfinite energy/gradient')
    return energy,gradient,{'rms_gradient_kcal_mol_A':rms,'gradient_gate_passed':rms<=TARGET,
                             'output_sha256':sha(path),'method':'SQM analytic derivative;sign and units checked by finite differences'}


class Evaluator:
    def __init__(self,out,elements,names,lock_fd=None,qmcharge=0):
        self.out=out;self.elements=elements;self.names=names;self.env=amber_environment();self.cache={};self.lock_fd=lock_fd
        if not isinstance(qmcharge,int):raise ValueError('QM charge must be an explicit integer')
        self.qmcharge=qmcharge
        self.count=max([int(p.name.split('-')[1]) for p in out.glob('evaluation-*')]+[0])

    @staticmethod
    def key(x):return hashlib.sha256(np.asarray(x,dtype=np.float64).tobytes()).hexdigest()

    def evaluate(self,x):
        xyz=np.asarray(x,dtype=np.float64).reshape(-1,3);key=self.key(xyz)
        if key in self.cache:
            r=self.cache[key];return r['energy'],r['gradient'].ravel()
        if len(xyz)!=len(self.elements) or not np.isfinite(xyz).all():raise ValueError('Invalid trial coordinates')
        if shutil.disk_usage(self.out).free < 2*1024**3:
            raise RuntimeError('Less than2GiB free; evaluated coordinates and outputs preserved')
        self.count+=1;folder=self.out/f'evaluation-{self.count:06d}';folder.mkdir()
        text=("AM1 analytic-gradient evaluation\n &qmmm\n"
              f" qm_theory='AM1', qmcharge={self.qmcharge}, scfconv=1.d-10, ndiis_attempts=700,\n"
              " maxcyc=0, verbosity=4, diag_routine=6, pseudo_diag=0,\n /\n")
        text+=''.join(f'{z:3d} {name:5s} '+ ' '.join(f'{v:.17g}' for v in row)+'\n'
                      for z,name,row in zip(self.elements,self.names,xyz))
        (folder/'sqm.in').write_text(text)
        # The exact optimizer coordinates are also saved independently of the
        # formatted QM input/output. This is not reconstructed output precision.
        np.save(folder/'input-coordinates.npy',xyz)
        kwargs={'pass_fds':(self.lock_fd,)} if self.lock_fd is not None else {}
        with (folder/'launch.log').open('x') as log:
            subprocess.run([self.env['AMBERHOME']+'/bin/sqm','-O','-i','sqm.in','-o','sqm.out'],
                           cwd=folder,env=self.env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=7200,**kwargs)
        energy,gradient,info=parse_singlepoint(folder/'sqm.out',len(xyz))
        validate_coordinates(folder/'sqm.out',xyz,self.elements)
        totals=re.findall(r'Total Mulliken Charge\s*=\s*('+NUMBER+')',(folder/'sqm.out').read_text())
        if not totals or abs(float(totals[-1])-self.qmcharge)>.0005:
            raise ValueError('Reported total QM charge differs from requested integer')
        record={'coordinates_A':xyz.tolist(),'energy_kcal_mol':energy,'gradient_kcal_mol_A':gradient.tolist(),
                'evaluation':folder.name,'input_coordinate_sha256':sha(folder/'input-coordinates.npy'),
                'qmcharge':self.qmcharge,**info}
        write_json(folder/'evaluation.json',record)
        if len(self.cache) >= 100:
            self.cache.pop(next(iter(self.cache)))
        self.cache[key]={'energy':energy,'gradient':gradient,'record':record,'path':folder}
        return energy,gradient.ravel()


def optimize(out,elements,names,xyz,mol=None,lock_fd=None,coordinate_scale=1.,qmcharge=0):
    if not 0 < coordinate_scale <= 1:
        raise ValueError('Invalid optimizer coordinate scale')
    evaluator=Evaluator(out,elements,names,lock_fd,qmcharge)
    accepted=0
    def save_accepted(x):
        nonlocal accepted
        evaluator.evaluate(x);value=evaluator.cache[evaluator.key(np.asarray(x).reshape(-1,3))]
        record=dict(value['record']);accepted+=1;record['accepted_iteration_this_session']=accepted
        if mol is not None:
            p=np.asarray(record['coordinates_A']);lengths=[float(np.linalg.norm(p[b.GetBeginAtomIdx()]-p[b.GetEndAtomIdx()])) for b in mol.GetBonds()]
            if min(lengths)<.7 or max(lengths)>2.0:raise ValueError('Accepted geometry fails bond-length envelope')
            record['bond_lengths_A']=[min(lengths),max(lengths)]
        write_json(out/'accepted-checkpoint.json',record)
        write_json(out/'status.json',{'state':'running','worker_pid':os.getpid(),'accepted_iteration_this_session':accepted,
                   'evaluations_total':evaluator.count,'rms_gradient_kcal_mol_A':record['rms_gradient_kcal_mol_A'],
                   'updated_utc':datetime.now(timezone.utc).isoformat()})
        print('accepted',accepted,'evaluations',evaluator.count,'RMS',record['rms_gradient_kcal_mol_A'],flush=True)
        return record
    def callback(intermediate_result):
        record=save_accepted(intermediate_result.x*coordinate_scale)
        if record['gradient_gate_passed']:raise StopIteration
    def scaled_evaluate(y):
        energy,gradient=evaluator.evaluate(y*coordinate_scale)
        return energy,gradient*coordinate_scale
    initial=save_accepted(xyz.ravel())
    if not initial['gradient_gate_passed']:
        result=minimize(scaled_evaluate,xyz.ravel()/coordinate_scale,jac=True,method='L-BFGS-B',callback=callback,
                        options={'gtol':1e-8*coordinate_scale,'ftol':0.,'maxiter':4000,'maxfun':20000,'maxls':40,'maxcor':10})
        final=json.loads((out/'accepted-checkpoint.json').read_text())
        if not final['gradient_gate_passed']:
            write_json(out/'optimizer-stop.json',{'message':str(result.message),'status':int(result.status),'success':bool(result.success)})
            raise RuntimeError('Optimizer stopped before the explicit AM1 gradient criterion')
    final=json.loads((out/'accepted-checkpoint.json').read_text())
    final.update(state='AM1_gradient_gate_passed',optimizer='SQM analytic gradient with SciPy L-BFGS-B',
                 final_output=str((out/final['evaluation']/'sqm.out').relative_to(ROOT)) if ROOT in out.parents else str(out/final['evaluation']/'sqm.out'))
    write_json(out/'completed.json',final);write_json(out/'status.json',{'state':'AM1_gradient_gate_passed','rms_gradient_kcal_mol_A':final['rms_gradient_kcal_mol_A']})
    return final


def run(dp):
    mol,info=preflight(dp);NATIVE.mkdir(exist_ok=True);BASE.mkdir(exist_ok=True)
    # Shared with native-SQM worker, so the same chain cannot run twice.
    with (NATIVE/f'pc{dp}.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        out=BASE/f'pc{dp}';out.mkdir(exist_ok=True)
        if (out/'completed.json').exists():return
        native=NATIVE/f'pc{dp}'
        elements,names,original=read_initial(native/'original-sqm.in')
        if elements != [a.GetAtomicNum() for a in mol.GetAtoms()]:
            raise ValueError('Original QM atom order differs from molecule graph')
        spec={'molecule_source':info,'engine_sha256':sha(Path(amber_environment()['AMBERHOME'])/'bin/sqm'),
              'AM1_SCF_threshold':1e-10,'AM1_gradient_RMS_target':TARGET,'scipy_version':scipy.__version__,
              'diagonalization':'DSYEVD,full diagonalization;Hamiltonian unchanged',
              'restart':'accepted coordinates kept at double precision;L-BFGS history restarts after interruption'}
        if (out/'spec.json').exists() and json.loads((out/'spec.json').read_text())!=spec:raise ValueError('Existing spec differs')
        write_json(out/'spec.json',spec)
        settings_path=ROOT/'config/pc-qm-continuation.json'
        settings=json.loads(settings_path.read_text()).get('optimizer_coordinate_scale',{})
        coordinate_scale=float(settings.get(str(dp),1.))
        policy={'coordinate_scale':coordinate_scale,'physical_gradient_gate_unchanged':TARGET,
                'meaning':'optimizer y=x/scale; energy E(scale*y), gradient scale*dE/dx; Hamiltonian unchanged'}
        policy_path=out/'optimizer-settings.json'
        if policy_path.exists() and json.loads(policy_path.read_text())!=policy:
            raise ValueError('Optimizer settings changed without a recorded recovery review')
        write_json(policy_path,policy)
        if (out/'accepted-checkpoint.json').exists():
            saved=json.loads((out/'accepted-checkpoint.json').read_text());xyz=np.array(saved['coordinates_A'])
            if sha(out/saved['evaluation']/'sqm.out')!=saved['output_sha256']:raise ValueError('Checkpoint evidence changed')
        else:
            checkpoint=native/'coordinate-checkpoint.json'
            if not checkpoint.exists():raise RuntimeError('Wait for native coordinate checkpoint;do not reset the chain')
            saved=json.loads(checkpoint.read_text());xyz=np.array(saved['coordinates_A'])
            if saved['atoms'] != elements or saved['names'] != names:
                raise ValueError('Native coordinate checkpoint atom order differs')
            candidates=list(native.glob(f'segment-{saved["segment"]:04d}-attempt-*/sqm.out'))
            if not any(sha(p)==saved['output_sha256'] for p in candidates):
                raise ValueError('Native checkpoint output evidence missing or changed')
            write_json(out/'initial-coordinate-provenance.json',{'source':str(checkpoint.relative_to(ROOT)),'sha256':sha(checkpoint),
                                                              'native_logged_iterations':saved['total_logged_iterations']})
        try:optimize(out,elements,names,xyz,mol,lock.fileno(),coordinate_scale)
        except Exception as error:
            write_json(out/'status.json',{'state':'blocked_preserved','reason':str(error)});raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dp',type=int,choices=[4,8,16],required=True)
    a=p.parse_args();run(a.dp)
