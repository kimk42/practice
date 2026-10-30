"""Independently check PC8's force-polished endpoint before numerical use."""
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import subprocess

import numpy as np

from pc_length_parameters import ROOT, preflight, sha
from pc_sqm_lbfgs import parse_singlepoint, validate_coordinates, TARGET
from run_pc_sqm_checkpointed import amber_environment, read_initial, write_json

OUT = ROOT/'runs/pc8-force-polish-review-v4'


def main():
    with (ROOT/'runs/pc-length-sqm-checkpointed-v2/pc8.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        base=ROOT/'runs/pc8-force-polish-v4'
        accepted=json.loads((base/'completed.json').read_text())
        raw=base/accepted['evaluation']/'sqm.out'
        if sha(raw)!=accepted['output_sha256']:
            raise ValueError('Accepted endpoint changed')
        mol, info=preflight(8)
        elements,names,_=read_initial(ROOT/'runs/pc-length-sqm-checkpointed-v2/pc8/original-sqm.in')
        xyz=np.array(accepted['coordinates_A'])
        energy,gradient,gate=parse_singlepoint(raw,len(elements))
        validate_coordinates(raw,xyz,elements)
        if not gate['gradient_gate_passed']:
            raise ValueError('Original gradient gate failed')
        # Proper rotation, determinant+1, followed by translation. Chemistry and
        # pair distances are unchanged; compare gradients in the original frame.
        angle=.37
        rotation=np.array([[np.cos(angle),-np.sin(angle),0.],[np.sin(angle),np.cos(angle),0.],[0.,0.,1.]])
        direction=gradient/np.linalg.norm(gradient)
        cases={'tight-density':xyz,'rotated-tight-density':xyz@rotation.T+np.array([1.25,-.75,.5])}
        for h in [.0005,-.0005,.001,-.001]:
            cases[f'displaced{h:+.4f}']=xyz+h*direction
        OUT.mkdir()  # Review results are immutable; interrupted reviews require inspection.
        write_json(OUT/'spec.json',{'source_sha256':sha(base/'completed.json'),'raw_sha256':sha(raw),
                   'SCF_energy_threshold':1e-10,'tight_p_conv':1,'original_RMS_gate':TARGET,
                   'rotation_energy_tolerance_kcal_mol':1e-4,'gradient_rotation_RMS_tolerance':5e-5,
                   'meaning':'numerical robustness only;not a Hessian/minimum or physical forcefield validation'})
        env=amber_environment();results={}
        for name,trial in cases.items():
            folder=OUT/name;folder.mkdir()
            text=('PC8 endpoint independent numerical review\n &qmmm\n'
                  " qm_theory='AM1', qmcharge=0, scfconv=1.d-10, ndiis_attempts=700,\n"
                  ' maxcyc=0,verbosity=4,diag_routine=6,pseudo_diag=0,tight_p_conv=1,\n /\n')
            text+=''.join(f'{z:3d} {label:5s} '+' '.join(f'{v:.17g}' for v in row)+'\n'
                          for z,label,row in zip(elements,names,trial))
            (folder/'sqm.in').write_text(text);np.save(folder/'input-coordinates.npy',trial)
            with (folder/'launch.log').open('x') as log:
                subprocess.run([env['AMBERHOME']+'/bin/sqm','-O','-i','sqm.in','-o','sqm.out'],cwd=folder,env=env,
                               stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300,pass_fds=(lock.fileno(),))
            e,g,audit=parse_singlepoint(folder/'sqm.out',len(elements))
            validate_coordinates(folder/'sqm.out',trial,elements)
            result={'energy_kcal_mol':e,'gradient_kcal_mol_A':g.tolist(),**audit}
            write_json(folder/'result.json',result);results[name]=result
            print(name,e,result['rms_gradient_kcal_mol_A'],flush=True)
        tight=results['tight-density'];rotated=results['rotated-tight-density']
        rotation_error=float(np.sqrt(np.mean((np.array(rotated['gradient_kcal_mol_A'])@rotation-np.array(tight['gradient_kcal_mol_A']))**2)))
        finite=[]
        for h in [.0005,.001]:
            numeric=(results[f'displaced{h:+.4f}']['energy_kcal_mol']-results[f'displaced{-h:+.4f}']['energy_kcal_mol'])/(2*h)
            finite.append({'step_norm_A':h,'central_difference':numeric,
                           'projected_tight_gradient':float(np.sum(np.array(tight['gradient_kcal_mol_A'])*direction))})
        passed=(tight['gradient_gate_passed'] and rotated['gradient_gate_passed'] and
                abs(rotated['energy_kcal_mol']-tight['energy_kcal_mol'])<1e-4 and rotation_error<5e-5)
        report={'checked_utc':datetime.now(timezone.utc).isoformat(),'numerical_endpoint_robustness_passed':bool(passed),
                'original_RMS_gradient':gate['rms_gradient_kcal_mol_A'],'tight_density_RMS_gradient':tight['rms_gradient_kcal_mol_A'],
                'rotated_RMS_gradient':rotated['rms_gradient_kcal_mol_A'],'gradient_rotation_RMS_difference':rotation_error,
                'rotation_energy_difference_kcal_mol':rotated['energy_kcal_mol']-tight['energy_kcal_mol'],
                'finite_difference_diagnostics':finite,'accepted_output':str(raw.relative_to(ROOT)),
                'accepted_output_sha256':sha(raw),'molecule_source':info,
                'limitations':'Finite-difference precision at tiny displacements remains limited;no Hessian/minimum,physical-PC or equilibrium approval. Only numerical parameter candidates may be considered after explicit review.',
                'physical_approval':False}
        write_json(OUT/'result.json',report)
        write_json(ROOT/'models/surface-revision-v2/pc8-force-polish-review-20261030.json',report)
        if not passed:
            raise ValueError('Independent endpoint robustness failed;no export')


if __name__=='__main__':
    main()
