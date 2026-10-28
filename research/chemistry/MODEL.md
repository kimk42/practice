# CBAA-1 chemical model and parameterization gate

## Connectivity reconstructed from primary figures
Source: Sarker et al., Chemical Science (2023), DOI 10.1039/d3sc01977b, Figure 1 (polymer) and Scheme 1 (small hydration probes). Local source images: ../d3sc01977b-f1.jpg and ../d3sc01977b-s1.jpg.

CBAA-1 acrylamide monomer, zwitterionic carboxylate state:
`C=CC(=O)NCCC[N+](C)(C)CC(=O)[O-]`
Formula C10H18N2O3, molecular weight 214.265 g/mol, net formal charge zero.

Polymer repeat: `–CH2–CH[–C(=O)–NH–(CH2)3–N+(CH3)2–CH2–C(=O)O−]–`.
Here “1” describes the single methylene between quaternary nitrogen and the carboxylate carbon, not the propyl linker to the amide. Polymerization consumes the acrylamide C=C bond. Tacticity and end groups remain modeling choices.

Scheme 1 instead represents CBAA1 hydration with trimethylglycine: `C[N+](C)(C)CC(=O)[O-]`, C5H11NO2. It lacks the acrylamide, amide-propyl linker and polymer backbone. Do not use its small-molecule model as the polymer repeat.

## Generated candidates
- cbaa1_monomer.sdf/.pdb: isolated monomer.
- cbaa1_capped_trimer_example.sdf/.pdb: saturated 3-repeat backbone with hydrocarbon end caps, 3 explicit example stereocenters; not a selected tacticity, ATRP end chemistry or grafted surface.
- paper_small_molecule_betaine.sdf/.pdb: hydration-probe comparator only.
SDF is the connectivity/charge authority; PDB alone does not retain all chemical information. RDKit 2025.9.6 ETKDG coordinates and generic MMFF minimization were used solely to generate geometries. No MMFF parameters are proposed for production GROMACS. Validation: one connected component, neutral formal charge, expected 1/3 zwitterion pairs, monomer formula, embedding success and geometry optimizer convergence. Formal charge validation does not establish atomistic partial-charge accuracy. See validation.json and build_candidates.py.

## Force-field decision
Current implementation uses GAFF2/AM1-BCC for polymer/PC with Amber99sb-ildn/TIP3P for protein/water. The earlier CHARMM36/CGenFF direction was superseded;it is not an instruction for new inputs. COP10,PC4 and one graft have numerical assignment checks but still need physical hydration/conformation validation. New chain lengths require consistent interior/end parameters and cannot inherit whole-PC4 terminal charges by simple concatenation.

The cited2023 paper reports CHARMM36/TIP3P with separately derived RESP charges;that literature protocol does not justify combining its parameters with the current Amber implementation. Maintain documented water/ion/mixing conventions and source provenance. MMFF geometry relaxation is not an MD forcefield.

## Surface and protein gates
The confirmed substrate is polycarbonate and the coating is UV-immobilized CBAA-1/BMA/BPAA at user-reported 7:2:1 composition. PC grade/repeat chemistry, DP, film loading/thickness, attachment chemistry and crosslink conversion remain unresolved. User pH7.4,cPBSzI and room temperature are known;exact temperature and unmeasured morphology may be explicit sensitivity assumptions. Begin with coating-only minimization/equilibration and assess density/height, water distribution and ion partitioning before protein addition. Silica and end-grafted homopolymer assumptions have been superseded.

3GHG requires biological-assembly selection, missing-region/glycan/disulfide review; 8OS5 is N-terminal FXII domains, not full FXII. Use separate fibrinogen and FXII systems initially. Full-protein surface MD requires box/atom-count and compute estimates. Classical MD can examine adsorption/conformation; it cannot demonstrate proteolytic FXII activation.

## Superseding user correction
The actual target is a UV-fixed 7:2:1 CBAA-1/BMA/BPAA coating on polycarbonate, not silica-supported CBAA homopolymer. See ../COATING.md for latest evidence and pending protocol details. Earlier silica and end-grafted brush assumptions must not drive model generation.

## Added component candidates
See ../COATING.md for evidence. BMA: C=C(C)C(=O)OCCCC, n-butyl methacrylate, supported by Lin 2018 primary SI. BPAA: C=CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1, N-(4-benzoylphenyl)acrylamide, now directly confirmed by uploaded Lin2020 Figure1 and section4.2, and Lin2025 section2.1. Existing candidate connectivity matches; target batch parameters remain unvalidated. Updated build_candidates.py generates both with charge/connectivity/geometry checks. Sequence candidates in sequences-dp100.json are only composition examples; they contain no force-field topology or crosslinks.


## User-approved three-component prototype
The user approved directly building CBAA-1/BMA/BPAA with a7:2:1 repeat-count ratio. See ../../models/README.md and ../../scripts/build_copolymer.py for the connected pre-UV chain. The older homotrimer/sequence-only artifacts are comparators, not the current target. DP20, seeded random sequence/stereochemistry and H termination are explicit prototype assumptions. BPAA identity is established by uploaded primary papers; a matching published terpolymer formulation is not a prerequisite for building the authorized model. Production parameterization and surface/post-UV network validation remain necessary.
