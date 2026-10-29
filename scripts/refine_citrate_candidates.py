"""Refine saved citrate geometries with an explicit final AM1 gradient gate.

Existing charges/controls remain untouched. New BCC outputs are candidates and
are not exported to MD by this script; solution validation remains separate.
"""
import fcntl
import json
from pathlib import Path
import re
import subprocess
from rdkit import Chem
import numpy as np
from pc_length_parameters import ROOT,sha
from pc_sqm_lbfgs import optimize
from run_pc_sqm_checkpointed import read_initial,parse_completed,amber_environment,write_json
from sqm_antechamber_view import write_view

BASE=ROOT/'runs/citrate-gradient-refinement-v2'


def run():
    BASE.mkdir(exist_ok=True)
    with (BASE/'worker.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        reports=[]
        for name,charge in [('citrate3-',-3),('hydrogen_citrate2-',-2)]:
            old=ROOT/'runs/buffer-charge-candidates-g5tzm8pz'/name
            source=ROOT/f'models/buffer-parameters/{name}.mol'
            mol=Chem.MolFromMolFile(str(source),removeHs=False)
            if Chem.GetFormalCharge(mol)!=charge:raise ValueError('Source charge differs')
            elements,names,_=read_initial(old/'sqm.in')
            if elements!=[a.GetAtomicNum() for a in mol.GetAtoms()]:raise ValueError('Source atom order differs')
            xyz,_=parse_completed(old/'sqm.out',elements)
            out=BASE/name;out.mkdir(exist_ok=True)
            spec={'species':name,'charge':charge,'source_mol_sha256':sha(source),
                  'source_sqm_sha256':sha(old/'sqm.out'),'original_bcc_sha256':sha(old/'candidate.mol2'),
                  'method':'SQM AM1 analytic gradient with L-BFGS;SCF1e-10;RMS gradient<=0.0005kcal/mol/A',
                  'initial_coordinates':'legacy final SQM table rounded to0.0001A;not a new random conformer',
                  'scope':'separate refined gas-phase candidate;no aqueous physical approval or populations'}
            if (out/'spec.json').exists() and json.loads((out/'spec.json').read_text())!=spec:
                raise ValueError('Existing refinement spec differs')
            write_json(out/'spec.json',spec)
            if not (out/'completed.json').exists():
                if (out/'accepted-checkpoint.json').exists():
                    cp=json.loads((out/'accepted-checkpoint.json').read_text())
                    if sha(out/cp['evaluation']/'sqm.out')!=cp['output_sha256']:raise ValueError('Saved QM evidence changed')
                    xyz=np.array(cp['coordinates_A'])
                optimize(out,elements,names,xyz,mol,lock.fileno(),qmcharge=charge)
            done=json.loads((out/'completed.json').read_text())
            sqm=ROOT/done['final_output']
            if sha(sqm)!=done['output_sha256'] or not done['gradient_gate_passed']:
                raise ValueError('Refined gradient evidence fails')
            param=out/'bcc-compatible-v2'
            if not (param/'candidate.mol2').exists():
                if param.exists():raise RuntimeError('Partial BCC output preserved;inspect before recovery')
                param.mkdir()
                write_view(sqm,param/'antechamber-sqm.out',mol.GetNumAtoms())
                cmd=['antechamber','-i','antechamber-sqm.out','-fi','sqmout','-o','candidate.mol2','-fo','mol2',
                     '-c','bcc','-at','gaff2','-j','4','-nc',str(charge),'-m','1','-eq','1','-dr','no','-pf','n','-s','2']
                with (param/'antechamber.log').open('x') as log:
                    subprocess.run(cmd,cwd=param,env=amber_environment(),stdout=log,stderr=subprocess.STDOUT,
                                   check=True,timeout=600,pass_fds=(lock.fileno(),))
            text=(param/'candidate.mol2').read_text()
            atoms=[s.split() for s in text.split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>')[0].splitlines() if s.strip()]
            bonds=[s.split() for s in text.split('@<TRIPOS>BOND')[1].split('@<TRIPOS>')[0].splitlines() if s.strip()]
            if len(atoms)!=mol.GetNumAtoms():raise ValueError('BCC atom count differs')
            if [a[5][0] for a in atoms]!=[a.GetSymbol().lower() for a in mol.GetAtoms()]:raise ValueError('BCC element order differs')
            graph={tuple(sorted((b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1))) for b in mol.GetBonds()}
            if graph!={tuple(sorted(map(int,r[1:3]))) for r in bonds}:raise ValueError('BCC graph differs')
            q=np.array([float(a[8]) for a in atoms])
            if not np.isfinite(q).all() or abs(q.sum()-charge)>len(q)*.0005:
                raise ValueError('Charge residual exceeds printed-atomic-charge rounding envelope')
            warnings=[r for r in (param/'antechamber.log').read_text().splitlines() if 'warning' in r.lower()]
            if warnings:raise ValueError('BCC warning requires review: '+str(warnings))
            report={'species':name,'state':'explicit AM1 gradient accepted;new BCC candidate retained for review',
                    'qmcharge':charge,'RMS_gradient_kcal_mol_A':done['rms_gradient_kcal_mol_A'],
                    'energy_kcal_mol':done['energy_kcal_mol'],'evaluations':len(list(out.glob('evaluation-*'))),
                    'raw_BCC_net_charge_e':float(q.sum()),'raw_charge_residual_e':float(q.sum()-charge),
                    'candidate_mol2':str((param/'candidate.mol2').relative_to(ROOT)),
                    'candidate_mol2_sha256':sha(param/'candidate.mol2'),'accepted_sqm_sha256':sha(sqm),
                    'original_bcc_unchanged':sha(old/'candidate.mol2')==spec['original_bcc_sha256'],
                    'MD_topology_exported':False,'physical_approval':False,'pH_population_assigned':False,
                    'remaining':['new candidate bond-order/charge symmetry review and integer-charge serialization',
                                 'bonded topology and isolated numerical water control','aqueous interactions and protonation microstate/speciation assessment']}
            write_json(out/'candidate-audit.json',report);reports.append(report)
        write_json(ROOT/'models/buffer-parameters/citrate-gradient-refinement-audit.json',{'species':reports,'full_buffer_ready':False})
        print(json.dumps(reports,indent=2),flush=True)


if __name__=='__main__':run()
