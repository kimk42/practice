"""Bounded FIRE diagnostic using unchanged SQM AM1 analytic gradients.

Energy line searches stalled near the PC8 numerical energy resolution. This
separate force-based relaxation is not an automatic parameter-export route.
Final derivatives still require the original RMS gate and further review.
FIRE: Bitzek et al., Physical Review Letters97,170201 (2006).
"""
from datetime import datetime, timezone
import fcntl
import json
import os

import numpy as np

from pc_length_parameters import ROOT, preflight, sha
from pc_sqm_lbfgs import Evaluator, TARGET
from run_pc_sqm_checkpointed import read_initial, write_json

OUT = ROOT/'runs/pc8-force-polish-v4'


def fire_step(gradient, velocity, dt, alpha, positive_steps):
    force = -gradient
    power = float(np.sum(velocity*force))
    if power > 0:
        velocity = (1-alpha)*velocity + alpha*force*np.linalg.norm(velocity)/np.linalg.norm(force)
        if positive_steps > 5:
            dt = min(dt*1.1, .02)
            alpha *= .99
        positive_steps += 1
    else:
        velocity = np.zeros_like(force)
        alpha = .1
        dt *= .5
        positive_steps = 0
    velocity = velocity + dt*force
    displacement = dt*velocity
    largest = float(np.linalg.norm(displacement,axis=1).max())
    if largest > .001:
        displacement *= .001/largest
    return displacement, velocity, dt, alpha, positive_steps


def main():
    with (ROOT/'runs/pc-length-sqm-checkpointed-v2/pc8.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        source = ROOT/'runs/pc-length-sqm-lbfgs-v3/pc8'
        original = json.loads((source/'accepted-checkpoint.json').read_text())
        if sha(source/original['evaluation']/'sqm.out') != original['output_sha256']:
            raise ValueError('Original accepted evidence changed')
        mol, info = preflight(8)
        elements, names, _ = read_initial(ROOT/'runs/pc-length-sqm-checkpointed-v2/pc8/original-sqm.in')
        spec = {'source_checkpoint_sha256':sha(source/'accepted-checkpoint.json'),
                'source_evaluation':original['evaluation'], 'molecule':info,
                'energy_model':'unchanged AM1,SCF1e-10,DSYEVD,full diagonalization',
                'optimizer':'FIRE with artificial unit masses;not physical dynamics',
                'initial_dt':.005,'max_dt':.02,'max_atom_displacement_A':.001,
                'max_evaluations':1000,'required_RMS_gradient_kcal_mol_A':TARGET,
                'scope':'separate diagnostic;not an automatic forcefield/MD handoff'}
        OUT.mkdir(exist_ok=True)
        if (OUT/'completed.json').exists():
            return
        if (OUT/'spec.json').exists() and json.loads((OUT/'spec.json').read_text()) != spec:
            raise ValueError('Existing diagnostic specification differs')
        write_json(OUT/'spec.json',spec)
        if (OUT/'checkpoint.json').exists():
            state = json.loads((OUT/'checkpoint.json').read_text())
            if sha(OUT/state['evaluation']/'sqm.out') != state['output_sha256']:
                raise ValueError('Diagnostic checkpoint evidence changed')
            xyz=np.array(state['coordinates_A']);velocity=np.array(state['velocity'])
            dt=state['dt'];alpha=state['alpha'];positive=state['positive_steps']
        else:
            xyz=np.array(original['coordinates_A']);velocity=np.zeros_like(xyz);dt=.005;alpha=.1;positive=0
        ev=Evaluator(OUT,elements,names,lock.fileno())
        try:
            while ev.count < 1000:
                energy, flat_gradient = ev.evaluate(xyz.ravel())
                gradient=flat_gradient.reshape(-1,3)
                result=dict(ev.cache[ev.key(xyz)]['record'])
                lengths=[float(np.linalg.norm(xyz[b.GetBeginAtomIdx()]-xyz[b.GetEndAtomIdx()])) for b in mol.GetBonds()]
                if min(lengths)<.7 or max(lengths)>2.0 or energy-original['energy_kcal_mol']>.05:
                    raise ValueError('Diagnostic geometry/energy excursion rejected')
                result.update(velocity=velocity.tolist(),dt=dt,alpha=alpha,positive_steps=positive,
                              bond_lengths_A=[min(lengths),max(lengths)])
                write_json(OUT/'checkpoint.json',result)
                write_json(OUT/'status.json',{'state':'diagnostic_running','worker_pid':os.getpid(),
                           'evaluation':result['evaluation'],'rms_gradient_kcal_mol_A':result['rms_gradient_kcal_mol_A'],
                           'updated_utc':datetime.now(timezone.utc).isoformat()})
                print(result['evaluation'],energy,result['rms_gradient_kcal_mol_A'],flush=True)
                if result['gradient_gate_passed']:
                    result.update(state='gradient_gate_passed_pending_independent_numerical_review',
                                  physical_approval=False,automatic_export=False)
                    write_json(OUT/'completed.json',result)
                    write_json(OUT/'status.json',{'state':result['state'],'automatic_export':False})
                    return
                displacement,velocity,dt,alpha,positive=fire_step(gradient,velocity,dt,alpha,positive)
                if dt < 1e-12:
                    raise ValueError('FIRE timestep stalled')
                xyz=xyz+displacement
            raise RuntimeError('Bounded diagnostic exhausted without gradient convergence')
        except Exception as error:
            write_json(OUT/'status.json',{'state':'held_preserved','reason':str(error)})
            raise


if __name__ == '__main__':
    main()
