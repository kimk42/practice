"""Inspect preserved candidate calculations,including method overrides;never approve MD here."""
import json,hashlib,re,shutil
from pathlib import Path
import numpy as np
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1]
def main():
    manifest=json.loads((ROOT/'runs/buffer-charge-active.json').read_text());base=ROOT/manifest['directory'];audit=ROOT/'models/buffer-parameters/charge-candidate-audit.json'
    if not (base/'first-pass-audit.json').exists():shutil.copy2(audit,base/'first-pass-audit.json')
    items=[s for s in json.loads((ROOT/'models/buffer-parameters/polyatomic-species-audit.json').read_text())['candidates'] if s['name']!='azide-'];reports=[]
    for item in items:
        directory=base/item['name'];row={'species':item['name'],'expected_charge':item['formal_charge'],'directory':str(directory.relative_to(ROOT)),'approved_for_MD':False}
        try:
            log=(directory/'calculation.log').read_text();warnings=[s for s in log.splitlines() if 'warning' in s.lower()]
            override=[s for s in warnings if "automatically changes the prediction type from '1' to '4' for bcc/abcg2" in s]
            other=[s for s in warnings if s not in override]
            row.update({'all_warnings_preserved':warnings,'requested_prediction_mode':1 if override else 4,'actual_prediction_mode':4,'method':'AM1-BCC/GAFF2 with full bond typing;isolated gas-phase candidate','method_override_review':'tool explicitly overrode requested mode1 to4;original mode-preservation assumption was not met;actual method is recorded,not suppressed' if override else 'explicit mode4','other_warnings':other})
            if other:raise ValueError('Unresolved calculation warning')
            if 'Calculation Completed' not in (directory/'sqm.out').read_text():raise ValueError('No quantum completion marker')
            mol=Chem.MolFromMolFile(str(directory/'input.mol'),removeHs=False)
            if Chem.GetFormalCharge(mol)!=item['formal_charge']:raise ValueError('Input charge differs')
            sections={};section=None
            for line in (directory/'candidate.mol2').read_text().splitlines():
                if line.startswith('@<TRIPOS>'):section=line[9:];sections[section]=[]
                elif section and line.strip():sections[section].append(line.split())
            atoms=sections['ATOM'];numbers={'c':6,'h':1,'o':8,'p':15};xyz=np.array([[float(v) for v in r[2:5]] for r in atoms]);charges=np.array([float(r[8]) for r in atoms])
            if [a.GetAtomicNum() for a in mol.GetAtoms()]!=[numbers[r[5][0]] for r in atoms]:raise ValueError('Element/order mismatch')
            target={tuple(sorted((int(r[1])-1,int(r[2])-1))):r[3] for r in sections['BOND']};edges={tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))):str(int(b.GetBondTypeAsDouble())) for b in mol.GetBonds()}
            if set(target)!=set(edges):raise ValueError('Connectivity changed')
            row.update({'observed_net_partial_charge':float(charges.sum()),'net_charge_residual_e':float(charges.sum()-item['formal_charge']),'coordinates_finite':bool(np.isfinite(xyz).all())})
            if abs(charges.sum()-item['formal_charge'])>1e-4 or not np.isfinite(xyz).all():raise ValueError('Charge/finite mismatch')
            lengths=[float(np.linalg.norm(xyz[i]-xyz[j])) for i,j in edges]
            if min(lengths)<.65 or max(lengths)>2.1:raise ValueError('Bonded geometry requires review')
            changes=[{'atoms_1based':[i+1 for i in edge],'source_order':value,'output_order':target[edge]} for edge,value in edges.items() if value!=target[edge]]
            row.update({'state':'calculated candidate numerical checks passed;physical/representation validation pending','atoms':len(atoms),'net_partial_charge':float(charges.sum()),'partial_charge_range_e':[float(charges.min()),float(charges.max())],'bond_length_range_A':[min(lengths),max(lengths)],'bond_order_representation_changes':changes,'candidate_mol2_sha256':hashlib.sha256((directory/'candidate.mol2').read_bytes()).hexdigest(),'qm_output_sha256':hashlib.sha256((directory/'sqm.out').read_bytes()).hexdigest(),'remaining':['source/output resonance representation and charge symmetry review','multi-anion gas-phase versus aqueous electrostatics validation','compatible nonbonded and ion-pairing behavior;no populations assigned']})
        except Exception as e:row.update({'state':'candidate blocked','reason':str(e)})
        reports.append(row)
    audit.write_text(json.dumps({'state':'preserved quantum candidates reviewed without rerunning','species':reports,'azide':'excluded;typing warning unresolved','full_buffer_ready':False,'MD_topologies_exported':False},indent=2)+'\n')
    print(json.dumps([{'species':r['species'],'state':r['state'],'net_charge':r.get('net_partial_charge')} for r in reports]))
if __name__=='__main__':main()
