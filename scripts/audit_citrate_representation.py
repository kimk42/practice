"""Read-only citrate graph/resonance audit; no charges or populations assigned."""
import json
from pathlib import Path
import re
from rdkit import Chem
from pc_length_parameters import ROOT, sha
from run_pc_sqm_checkpointed import write_json


def review(name, candidate=None):
    source=ROOT/f'models/buffer-parameters/{name}.mol'
    folder=ROOT/'runs/buffer-charge-candidates-g5tzm8pz'/name
    candidate=Path(candidate) if candidate else folder/'candidate.mol2'
    mol=Chem.MolFromMolFile(str(source),removeHs=False)
    expected=-3 if name=='citrate3-' else -2
    if mol is None or Chem.GetFormalCharge(mol)!=expected:
        raise ValueError('Formal source charge differs')
    sections={};section=None
    for line in candidate.read_text().splitlines():
        if line.startswith('@<TRIPOS>'):
            section=line[9:];sections[section]=[]
        elif section and line.strip():sections[section].append(line.split())
    atoms=sections['ATOM'];n=mol.GetNumAtoms()
    if len(atoms)!=n or [int(a[0]) for a in atoms]!=list(range(1,n+1)):
        raise ValueError('Atom numbering/count differs')
    elements={'o':8,'c':6,'h':1}
    if [elements.get(a[5][0]) for a in atoms]!=[a.GetAtomicNum() for a in mol.GetAtoms()]:
        raise ValueError('Atom elements/order differ')
    charges=[float(a[8]) for a in atoms]
    if abs(sum(charges)-expected)>1e-6:raise ValueError('Incorrect net charge')
    bonds={tuple(sorted((int(a[1])-1,int(a[2])-1))):a[3] for a in sections['BOND']}
    original={tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))):str(int(b.GetBondTypeAsDouble())) for b in mol.GetBonds()}
    if bonds.keys()!=original.keys():raise ValueError('Bond graph differs')
    groups=mol.GetSubstructMatches(Chem.MolFromSmarts('[O-]-[C](=[O])'))
    allowed={tuple(sorted((c,o))) for _,c,o in groups}
    changes=[]
    for edge,order in original.items():
        if bonds[edge]!=order:
            if edge not in allowed or order!='2' or bonds[edge]!='1':
                raise ValueError('Bond-order change outside carboxylate serialization')
            changes.append([i+1 for i in edge])
    checked=[]
    for o1,c,o2 in groups:
        if [atoms[o1][5],atoms[c][5],atoms[o2][5]]!=['o','c','o']:
            raise ValueError('Carboxylate atom types differ')
        if abs(charges[o1]-charges[o2])>1e-8:raise ValueError('Carboxylate oxygen charges differ')
        if any(bonds[tuple(sorted((c,o)))]!='1' for o in (o1,o2)):
            raise ValueError('Inconsistent equivalent carboxylate bond representation')
        checked.append({'atoms_1based':[o1+1,c+1,o2+1], 'equal_oxygen_charge_e':charges[o1]})
    acids=mol.GetSubstructMatches(Chem.MolFromSmarts('[O;H1]-[C](=[O])'))
    locations=[]
    for oh,c,o in acids:
        if atoms[oh][5]!='oh' or atoms[o][5]!='o' or bonds[tuple(sorted((c,o)))]!='2':
            raise ValueError('Protonated carboxylic acid representation changed')
        attached=[a for a in mol.GetAtomWithIdx(c).GetNeighbors() if a.GetAtomicNum()==6]
        assert len(attached)==1
        central=any(a.GetAtomicNum()==8 for a in attached[0].GetNeighbors())
        locations.append({'carboxyl_carbon_1based':c+1,'site':'central' if central else 'terminal'})
    if len(groups)!=(-expected) or len(acids)!=(3+expected):
        raise ValueError('Citrate protonation counts differ')
    sqm=folder/'sqm.out';text=sqm.read_text()
    progress=re.findall(r'^xmin\s+(\d+)\s+([-\d.]+) kcal/mol\s+([-\d.]+) kcal/\(mol\*A\)',text,re.M)
    return {'species':name,'atoms':n,'formal_charge':expected,'partial_charge_sum_e':sum(charges),
            'carboxylates':checked,'protonated_carboxyl_sites':locations,
            'changed_bond_orders_1based':changes,
            'representation_decision':'consistent with equivalent GAFF carboxylate serialization;formal MOL graph remains chemistry authority',
            'microstate_limit':'hydrogen citrate candidate protonates one terminal carboxyl;central/terminal protonation populations not established' if acids else None,
            'QM_log_review':{'last_printed_iteration':int(progress[-1][0]),'last_printed_RMS_kcal_mol_A':float(progress[-1][2]),
                'final_gradient_numerically_certified':False,
                'interpretation':'legacy trace prints every10iterations;last printed gradient need not be the final gradient. This does not establish failed convergence or satisfy a new explicit final-gradient gate.'},
            'source_sha256':sha(source),'candidate_mol2_sha256':sha(candidate),'sqm_sha256':sha(sqm),
            'charges_modified':False,'MD_started':False,'approved_for_full_buffer_MD':False}


if __name__=='__main__':
    reports=[review(name) for name in ('citrate3-','hydrogen_citrate2-')]
    write_json(ROOT/'models/buffer-parameters/citrate-representation-review.json',{
        'species':reports,'method':'read existing formal graphs,MOL2,QM logs;no new QM or MD',
        'resolved':'carboxylate bond-order serialization and within-group charge symmetry for the existing candidates',
        'remaining':['explicit final-gradient evidence before any revised charge assignment',
                     'protonation microstate/speciation choices','aqueous electrostatics,hydration and ion-pairing validation'],
        'full_buffer_ready':False})
    print(json.dumps(reports,indent=2))
