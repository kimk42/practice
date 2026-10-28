"""Independent PC length candidates and quantitative review of existing profiles."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from build_pc_graft import pc_graph, relax

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'models/surface-revision-v2'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def profile(name):
    path=ROOT/f'models/analysis/{name}-density.csv'
    rows=list(csv.DictReader(path.open()))
    z=np.array([float(x['z_nm']) for x in rows]); dz=float(np.median(np.diff(z)))
    pc=np.array([float(x['PC_kg_m3']) for x in rows]); water=np.array([float(x['water_kg_m3']) for x in rows])
    coating=np.array([float(x['coating_polymer_kg_m3']) for x in rows])
    assert np.allclose(np.diff(z),dz) and pc.sum()>0 and coating.sum()>0
    # Whole bins containing the5th through95th percentiles of PC mass.
    cumulative=np.cumsum(pc)/pc.sum()
    lo,hi=[int(np.searchsorted(cumulative,q)) for q in [.05,.95]]
    return {'source':str(path.relative_to(ROOT)),'sha256':sha(path),
            'PC_mass_5_to_95_percent_interval_bin_edges_nm':[float(z[lo]-dz/2),float(z[hi]+dz/2)],
            'PC_mean_in_interval_kg_m3':float(pc[lo:hi+1].mean()),
            'water_mean_in_same_interval_kg_m3':float(water[lo:hi+1].mean()),
            'coating_mass_centroid_z_nm':float(np.dot(z,coating)/coating.sum()),
            'PC_mass_centroid_z_nm':float(np.dot(z,pc)/pc.sum()),
            'scope':'xy/time-averaged profile descriptor;not local porosity,bulk density or an equilibrium/convergence test'}


def main():
    if OUT.exists(): raise RuntimeError('Existing revision candidates preserved')
    OUT.mkdir()
    profiles=[profile(n) for n in ['initial-coating','coating100ps']]
    (OUT/'existing-profile-assessment.json').write_text(json.dumps(profiles,indent=2)+'\n')
    candidates=[]
    for dp in [8,16]:
        folder=OUT/f'pc{dp}';folder.mkdir()
        mol=relax(pc_graph(dp),2026+dp,f'bpa-pc-dp{dp}-phenolic-ends',folder)
        assert Chem.GetFormalCharge(mol)==0 and len(Chem.GetMolFrags(mol))==1
        assert rdMolDescriptors.CalcMolFormula(mol)==f'C{16*dp-1}H{14*dp+2}O{3*dp-1}'
        assert len(mol.GetSubstructMatches(Chem.MolFromSmarts('[O;H1]-c')))==2
        xyz=mol.GetConformer().GetPositions()
        lengths=[float(np.linalg.norm(xyz[b.GetBeginAtomIdx()]-xyz[b.GetEndAtomIdx()])) for b in mol.GetBonds()]
        assert np.isfinite(xyz).all() and min(lengths)>.07*10 and max(lengths)<.20*10
        candidates.append({'PC_BPA_units':dp,'formula':rdMolDescriptors.CalcMolFormula(mol),
                           'atoms':mol.GetNumAtoms(),'phenolic_ends':2,'OH_ends_per_BPA_unit':2/dp,
                           'bond_lengths_A':[min(lengths),max(lengths)],
                           'coordinate_files':{str(p.relative_to(ROOT)):sha(p) for p in folder.glob('*.mol')},
                           'state':'connected MMFF-relaxed geometry;no new MD forcefield assigned'})
    report={'state':'surface-revision inputs and baseline diagnostics prepared',
            'purpose':'PC4/PC8/PC16 finite-chain-end sensitivity;not a measured PC molecular-weight distribution',
            'legacy_PC4_OH_ends_per_BPA_unit':.5,'new_candidates':candidates,
            'existing_profile_assessment':profiles,
            'required_before_surface_MD':['derive/review consistent longer-PC parameters without copying terminal charges onto interior repeats',
                 'pack independent substrate realizations and assess density/chain-end/water response',
                 'choose boundary/support and coating scenarios explicitly'],
            'full_FXII_runs_started':False}
    (OUT/'candidate-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(candidates,indent=2))


if __name__=='__main__':main()
