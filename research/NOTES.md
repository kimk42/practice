# CBAA-1 coating and blood-protein MD: literature notes

Current target and execution instructions: see ../README.md and COATING.md. The substrate is polycarbonate with UV-fixed CBAA-1/BMA/BPAA coating (7:2:1), not silica.

## Retrieved primary literature
- Sarker et al., including Shaoyi Jiang, Hydration behaviors of nonfouling zwitterionic materials (2023), DOI https://doi.org/10.1039/d3sc01977b. Full text: PMC10337769.xml. Compares TMAO and CBAA1/2/3. Describes preparatory GROMACS 2019.6 classical MD using CHARMM36/TIP3P with separately estimated partial charges and subsequent AIMD. This is not a ready-made validated polymer-brush topology. Exact monomer/repeat chemistry must be checked against Scheme 1, Figure 1 and SI before implementation; small zwitterion hydration models must not be mistaken for polymer repeats.
- Crystal structure of coagulation factor XII N-terminal domains 1–5 (2025), DOI https://doi.org/10.1107/S2059798325005297. Full text: PMC12216678.xml. PDB 8OS5. Covers N-terminal domains, not full-length FXII or activated FXIIa. Reports residues 18–278 in two crystallographic copies and 77–278 in another; numbering must be mapped to the selected construct/reference sequence. Crystal dimer must not automatically be used as physiological oligomer.

## Downloaded structures, not prepared simulation inputs
- 3GHG.pdb, 3GHG.cif: human fibrinogen; DOI https://doi.org/10.1021/bi802205g. Entry metadata lists assemblies 1 and 2. PDB contains extensive REMARK 465 missing-residue annotations across the deposited copies. Select biological assembly and evaluate missing regions, glycans, disulfides and ligands before topology generation. Deposited asymmetric unit is not automatically one biological molecule.
- 8OS5.cif: candidate domain-level FXII surface-binding model, explicitly a fragment. Full-length FXII modeling requires separate evaluation; isolated protease domains cannot stand in for surface-binding heavy chain.

## Execution and remaining work
CPU-only GROMACS 2025.2 available via /workspace/md-tools/bin/gmx. Water engine smoke test passed previously. No brush/protein trajectory has been generated. The smoke test uses Amber water setup solely for installation validation; it does not choose the research force field.
Crossref, Europe PMC full text, RCSB metadata and coordinate downloads now work. Official Lin 2018/2020 SI downloads via Figshare/S3 now work. CBAA-1 structure figures have been checked and component candidates generated; see chemistry/MODEL.md. Next: obtain compatible polymer/PC parameters and resolve BPAA primary synthesis, DP, film thickness/loading, crosslink state, salt/pH/temperature and protein scope; validate coating hydration before adding proteins. Whole fibrinogen plus a surface will need a resource estimate and likely GPU infrastructure for meaningful sampling. Adsorption and classical MD alone do not demonstrate FXII proteolytic activation.

Update: Europe PMC supplementaryFiles download succeeded. Archived ZIP, extracted SC-014-D3SC01977B-s001.pdf and CBAA-SI.txt are retained here. SI includes the partial-charge derivation method (RESP fitting / Gaussian16); this does not by itself provide a validated polymer repeat topology. The direct PMC PDF URL returned 404, but the Europe PMC ZIP route worked.
