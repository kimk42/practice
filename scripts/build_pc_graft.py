"""Build explicitly provisional BPA-PC oligomer and one post-UV graft geometry.

The user reports PC grade unknown; BPA-PC and a single methyl-site attachment
are illustrative choices, not measured substrate chemistry or conversion.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

ROOT=Path(__file__).resolve().parents[1]

def pc_graph(dp):
    graph=Chem.RWMol(); right=None
    for _ in range(dp):
        unit=Chem.MolFromSmiles('Oc1ccc(C(C)(C)c2ccc(O)cc2)cc1')
        oxygens=[a.GetIdx() for a in unit.GetAtoms() if a.GetAtomicNum()==8]
        offset=graph.GetNumAtoms(); graph.InsertMol(unit)
        left,current_right=[offset+i for i in oxygens]
        if right is not None:
            carbon=graph.AddAtom(Chem.Atom(6)); oxygen=graph.AddAtom(Chem.Atom(8))
            graph.AddBond(right,carbon,Chem.BondType.SINGLE)
            graph.AddBond(carbon,oxygen,Chem.BondType.DOUBLE)
            graph.AddBond(carbon,left,Chem.BondType.SINGLE)
        right=current_right
    mol=graph.GetMol();Chem.SanitizeMol(mol)
    return mol

def relax(mol, seed, name, out):
    mol=Chem.AddHs(mol)
    params=AllChem.ETKDGv3();params.randomSeed=seed;params.useRandomCoords=True;params.maxIterations=2000
    print('Embedding',name,mol.GetNumAtoms(),'atoms',flush=True)
    if AllChem.EmbedMolecule(mol,params)!=0:raise RuntimeError('Embedding failed: '+name)
    if not AllChem.MMFFHasAllMoleculeParams(mol):raise RuntimeError('Missing generic geometry params')
    ff=AllChem.MMFFGetMoleculeForceField(mol,AllChem.MMFFGetMoleculeProperties(mol))
    for batch in range(10):
        print(name,'geometry batch',batch+1,flush=True)
        if ff.Minimize(maxIts=1000)==0:break
    else:raise RuntimeError('Geometry minimization not converged: '+name)
    mol.SetProp('_Name',name);mol.SetProp('STATUS','Provisional geometry only; no production parameters or film equilibration')
    with Chem.SDWriter(str(out/(name+'.sdf'))) as w:w.write(mol)
    p=out/(name+'.sdf');lines=p.read_text().splitlines();lines[1]='     RDKit          3D';p.write_text('\n'.join(lines)+'\n')
    Chem.MolToPDBFile(mol,str(out/(name+'.pdb')))
    Chem.MolToMolFile(mol,str(out/(name+'.mol')))
    return mol

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pc-dp',type=int,default=4);a=p.parse_args()
    if a.pc_dp<2:p.error('Use at least two BPA units')
    out=ROOT/'models'/'pc-graft';out.mkdir(parents=True,exist_ok=True)
    pc=pc_graph(a.pc_dp)
    expected=f'C{16*a.pc_dp-1}H{14*a.pc_dp+2}O{3*a.pc_dp-1}'
    assert rdMolDescriptors.CalcMolFormula(pc)==expected
    cop=Chem.SDMolSupplier(str(ROOT/'models/copolymer/cbaa1-bma-bpaa-dp10-seed2026.sdf'),removeHs=True)[0]
    ketones=[]
    for bond in cop.GetBonds():
        if bond.GetBondType()!=Chem.BondType.DOUBLE:continue
        c,o=bond.GetBeginAtom(),bond.GetEndAtom()
        if c.GetAtomicNum()==8:c,o=o,c
        if c.GetAtomicNum()==6 and o.GetAtomicNum()==8 and sum(n.GetIsAromatic() for n in c.GetNeighbors())==2:
            ketones.append((c.GetIdx(),o.GetIdx()))
    assert len(ketones)==1
    methyls=[x.GetIdx() for x in pc.GetAtoms() if x.GetAtomicNum()==6 and x.GetTotalNumHs()==3
             and any(sum(n.GetIsAromatic() for n in center.GetNeighbors())==2 for center in x.GetNeighbors())]
    assert methyls
    c,o=ketones[0];offset=cop.GetNumAtoms();partner=offset+methyls[0]
    graph=Chem.RWMol(Chem.CombineMols(cop,pc));graph.GetBondBetweenAtoms(c,o).SetBondType(Chem.BondType.SINGLE)
    graph.GetAtomWithIdx(o).SetNoImplicit(False)
    graph.AddBond(c,partner,Chem.BondType.SINGLE)
    graft=graph.GetMol();Chem.SanitizeMol(graft)
    # Hydrogen abstraction/recombination transfers H from PC methyl to ketone O.
    assert graft.GetAtomWithIdx(o).GetTotalNumHs()==1
    assert graft.GetAtomWithIdx(partner).GetTotalNumHs()==2
    assert len(Chem.GetMolFrags(graft))==1 and Chem.GetFormalCharge(graft)==0
    def elements(m):return Counter(x.GetSymbol() for x in Chem.AddHs(m).GetAtoms())
    assert elements(graft)==elements(cop)+elements(pc)
    pc_name=f'bpa-pc-dp{a.pc_dp}-phenolic-ends'
    graft_name=f'cop10-graft-bpa-pc{a.pc_dp}-one-site'
    relax(pc,2026,pc_name,out)
    relax(graft,2027,graft_name,out)
    report={'substrate_grade':'experimentally unknown; BPA-PC is a user-accepted provisional assumption',
            'pc_dp':a.pc_dp,'pc_formula':expected,'pc_end_groups':'phenolic OH; finite-oligomer prototype',
            'copolymer_counts':{'CBAA1':7,'BMA':2,'BPAA':1},'graft_formula':rdMolDescriptors.CalcMolFormula(graft),
            'selected_attachment':'one BPA methyl C-H site; ketone C becomes alcohol-bearing C with new C-C bond',
            'mapped_reaction_atoms_1based':{'BPAA_C':c+1,'BPAA_O':o+1,'PC_C':partner+1},
            'crosslink_assumption':'one BPAA reacted in DP10 (1/1) solely to demonstrate one possible graft; not measured conversion',
            'not_included':['bulk PC surface/film','interchain crosslinks','production force-field parameters','equilibration'],
            'status':'chemical-graph/geometry candidates, not an equilibrated surface',
            'smiles':{'pc':Chem.MolToSmiles(pc),'graft':Chem.MolToSmiles(graft)}}
    (out/'assumptions.json').write_text(json.dumps(report,indent=2)+'\n');print('PC/graft candidates built',flush=True)
if __name__=='__main__':main()
