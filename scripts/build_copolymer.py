"""Build a real connected 7:2:1 terpolymer graph and unreacted 3D prototype.

No production force-field assignment or surface/UV reaction is implied.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import random
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

MONOMERS = {
    'CBAA1': 'C=CC(=O)NCCC[N+](C)(C)CC(=O)[O-]',
    'BMA': 'C=C(C)C(=O)OCCCC',
    'BPAA': 'C=CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1',
}
RESIDUES = {'CBAA1': 'CBA', 'BMA': 'BMA', 'BPAA': 'BPA'}
ROOT = Path(__file__).resolve().parents[1]


def build_graph(sequence, seed):
    graph = Chem.RWMol()
    units = []
    for unit_id, component in enumerate(sequence, 1):
        monomer = Chem.MolFromSmiles(MONOMERS[component])
        vinyl = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in monomer.GetBonds()
                 if b.GetBondType() == Chem.BondType.DOUBLE
                 and all(a.GetAtomicNum() == 6 for a in (b.GetBeginAtom(), b.GetEndAtom()))]
        if len(vinyl) != 1:
            raise ValueError('Ambiguous vinyl group')
        left, right = vinyl[0]
        if monomer.GetAtomWithIdx(left).GetTotalNumHs() != 2:
            left, right = right, left
        offset = graph.GetNumAtoms()
        graph.InsertMol(monomer)
        left += offset
        right += offset
        graph.GetBondBetweenAtoms(left, right).SetBondType(Chem.BondType.SINGLE)
        for idx in range(offset, graph.GetNumAtoms()):
            graph.GetAtomWithIdx(idx).SetIntProp('unit_id', unit_id)
        if units:
            graph.AddBond(units[-1]['alpha'], left, Chem.BondType.SINGLE)
        units.append({'unit_id': unit_id, 'component': component, 'head': left, 'alpha': right})
    mol = graph.GetMol()
    Chem.SanitizeMol(mol)
    randomizer = random.Random(seed + 10000)
    alpha_indices = {u['alpha'] for u in units}
    for idx, _ in Chem.FindMolChiralCenters(mol, includeUnassigned=True):
        if idx in alpha_indices:
            mol.GetAtomWithIdx(idx).SetChiralTag(randomizer.choice([
                Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW]))
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return Chem.AddHs(mol), units


def validate_graph(mol, units, counts):
    if len(Chem.GetMolFrags(mol)) != 1 or Chem.GetFormalCharge(mol) != 0:
        raise ValueError('Disconnected or non-neutral chain')
    charged = Counter(a.GetFormalCharge() for a in mol.GetAtoms())
    if charged[1] != counts['CBAA1'] or charged[-1] != counts['CBAA1']:
        raise ValueError('Zwitterion count mismatch')
    # Independent elemental balance: addition polymerization preserves monomer
    # atoms; hydrogen-terminated open chain adds two terminal H atoms.
    expected = Counter({'H': 2})
    for name, count in counts.items():
        for atom in Chem.AddHs(Chem.MolFromSmiles(MONOMERS[name])).GetAtoms():
            expected[atom.GetSymbol()] += count
    actual = Counter(a.GetSymbol() for a in mol.GetAtoms())
    if actual != expected:
        raise ValueError(f'Elemental balance mismatch: {actual} != {expected}')
    cc_double = [b for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE
                 and b.GetBeginAtom().GetAtomicNum() == b.GetEndAtom().GetAtomicNum() == 6]
    if cc_double:
        raise ValueError('Unreacted vinyl carbon-carbon double bond remains')
    if mol.GetRingInfo().NumRings() != 2 * counts['BPAA']:
        raise ValueError('BPAA aromatic ring count mismatch')
    for prev, nxt in zip(units, units[1:]):
        bond = mol.GetBondBetweenAtoms(prev['alpha'], nxt['head'])
        if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
            raise ValueError('Incorrect backbone link')
    for unit in units:
        if unit['component'] == 'BMA':
            atom = mol.GetAtomWithIdx(unit['alpha'])
            if any(n.GetAtomicNum() == 1 for n in atom.GetNeighbors()) and unit != units[-1]:
                raise ValueError('Interior methacrylate alpha carbon wrongly hydrogenated')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dp', type=int, default=20)
    parser.add_argument('--seed', type=int, default=2026)
    parser.add_argument('--output', type=Path, default=ROOT / 'models' / 'copolymer')
    args = parser.parse_args()
    if args.dp < 10 or args.dp % 10:
        parser.error('DP must be a positive multiple of ten (at least ten)')
    counts = {'CBAA1': 7 * args.dp // 10, 'BMA': 2 * args.dp // 10, 'BPAA': args.dp // 10}
    sequence = [name for name, n in counts.items() for _ in range(n)]
    random.Random(args.seed).shuffle(sequence)
    mol, units = build_graph(sequence, args.seed)
    validate_graph(mol, units, counts)
    print('Connectivity verified; embedding', mol.GetNumAtoms(), 'atoms', flush=True)
    params = AllChem.ETKDGv3()
    params.randomSeed = args.seed
    params.useRandomCoords = True
    params.numThreads = 2
    params.maxIterations = 2000
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError('3D embedding failed')
    print('3D embedding complete; relaxing geometry', flush=True)
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise RuntimeError('MMFF geometry preparation unsupported')
    properties = AllChem.MMFFGetMoleculeProperties(mol)
    forcefield = AllChem.MMFFGetMoleculeForceField(mol, properties)
    initial_energy = forcefield.CalcEnergy()
    # Reuse the optimizer state across bounded batches; fail rather than silently
    # presenting a nonconverged geometry as minimized.
    converged = False
    for batch in range(10):
        print('MMFF batch', batch + 1, flush=True)
        if forcefield.Minimize(maxIts=1000) == 0:
            converged = True
            break
    if not converged:
        raise RuntimeError('MMFF geometry relaxation did not converge')
    coordinates = np.asarray(mol.GetConformer().GetPositions())
    distances = [float(np.linalg.norm(coordinates[b.GetBeginAtomIdx()] - coordinates[b.GetEndAtomIdx()]))
                 for b in mol.GetBonds()]
    if not np.isfinite(coordinates).all() or min(distances) < 0.8 or max(distances) > 2.0:
        raise ValueError('Invalid coordinates or bond lengths')
    # Populate chemically distinct residues for downstream atom mapping.
    counters = Counter()
    atom_mapping = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            atom.SetIntProp('unit_id', atom.GetNeighbors()[0].GetIntProp('unit_id'))
        unit_id = atom.GetIntProp('unit_id')
        component = sequence[unit_id - 1]
        counters[unit_id] += 1
        name = f'{atom.GetSymbol()}{counters[unit_id]:03}'
        info = Chem.AtomPDBResidueInfo()
        info.SetName(name)
        info.SetResidueName(RESIDUES[component])
        info.SetResidueNumber(unit_id)
        info.SetChainId('A')
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
        atom_mapping.append({'atom_index_1based': atom.GetIdx() + 1, 'unit_id': unit_id,
                             'component': component, 'pdb_atom_name': name,
                             'element': atom.GetSymbol(), 'formal_charge': atom.GetFormalCharge()})
    reactive_sites = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        carbon, oxygen = (a, b) if a.GetAtomicNum() == 6 else (b, a)
        if carbon.GetAtomicNum() != 6 or oxygen.GetAtomicNum() != 8:
            continue
        aromatic_neighbors = [n for n in carbon.GetNeighbors() if n.GetIsAromatic()]
        if len(aromatic_neighbors) == 2:
            reactive_sites.append({'unit_id': carbon.GetIntProp('unit_id'),
                                   'ketone_C_1based': carbon.GetIdx() + 1,
                                   'ketone_O_1based': oxygen.GetIdx() + 1,
                                   'state': 'unreacted; attachment partner not assigned'})
    if len(reactive_sites) != counts['BPAA']:
        raise ValueError('Benzophenone ketone mapping mismatch')
    stem = f'cbaa1-bma-bpaa-dp{args.dp}-seed{args.seed}'
    args.output.mkdir(parents=True, exist_ok=True)
    mol.SetProp('_Name', stem)
    mol.SetProp('SEQUENCE', ','.join(sequence))
    mol.SetProp('STATUS', 'Unreacted terpolymer geometry only; not a production force-field topology')
    with Chem.SDWriter(str(args.output / f'{stem}.sdf')) as writer:
        writer.write(mol)
    sdf = args.output / f'{stem}.sdf'
    lines = sdf.read_text().splitlines()
    lines[1] = '     RDKit          3D'
    sdf.write_text('\n'.join(lines) + '\n')
    Chem.MolToPDBFile(mol, str(args.output / f'{stem}.pdb'))
    mol_path = args.output / f'{stem}.mol'
    Chem.MolToMolFile(mol, str(mol_path))
    mdl_lines = mol_path.read_text().splitlines()
    mdl_lines[1] = '     RDKit          3D'
    mol_path.write_text('\n'.join(mdl_lines) + '\n')
    reloaded = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
    if reloaded is None or Chem.MolToSmiles(reloaded) != Chem.MolToSmiles(mol):
        raise ValueError('SDF round-trip changed molecular graph/stereochemistry')
    report = {
        'name': stem, 'dp': args.dp, 'seed': args.seed, 'counts': counts, 'sequence': sequence,
        'atom_count': mol.GetNumAtoms(), 'formula': rdMolDescriptors.CalcMolFormula(mol),
        'molar_mass_g_mol': Descriptors.MolWt(mol), 'formal_charge': Chem.GetFormalCharge(mol),
        'monomer_smiles': MONOMERS, 'polymer_smiles': Chem.MolToSmiles(Chem.RemoveHs(mol)),
        'assumptions': ['7:2:1 is the user-approved modeling repeat-count ratio; experimental feed/incorporation basis unresolved',
                        'DP is a computational prototype, not measured chain length',
                        'Exact-count uniformly shuffled sequence; no reactivity-ratio inference',
                        'Random backbone stereochemical tags; real tacticity unmeasured',
                        'Hydrogen-terminated backbone; real initiator/end groups unmodeled',
                        'Isolated unreacted chain; no PC substrate or assigned UV crosslinks'],
        'geometry': {'method': 'ETKDGv3 then generic MMFF geometry minimization only',
                     'converged': converged, 'initial_energy_kcal_mol': initial_energy,
                     'final_energy_kcal_mol': forcefield.CalcEnergy(),
                     'bond_length_range_A': [min(distances), max(distances)]},
        'validated': ['exact repeat counts', 'single connected component', 'elemental balance',
                      'zwitterion formal charges', 'no remaining vinyl C=C', 'BPAA aromatic ring count',
                      'head-to-tail backbone links', 'BMA alpha-methyl connectivity',
                      'finite coordinates and bond lengths', 'SDF graph/stereochemistry round-trip'],
        'bpaa_potential_reaction_sites': reactive_sites, 'atom_mapping': atom_mapping,
        'status': 'not ready for MD: partial charges, bonded/nonbonded parameters and post-UV network unresolved',
    }
    (args.output / f'{stem}.json').write_text(json.dumps(report, indent=2) + '\n')
    print('Built', stem, counts, report['formula'], 'net charge', report['formal_charge'], flush=True)


if __name__ == '__main__':
    main()
