# CBAA-1/BMA/BPAA coating on polycarbonate: GROMACS research

User target: UV-immobilized CBAA-1-co-BMA-co-BPAA coating on polycarbonate, reported composition **7:2:1**, followed by separate fibrinogen and Factor XII interaction studies. The user confirmed a7:2:1 synthesis molar ratio;measured incorporation remains unreported. Current scope is coating and FXII/cPBSzI. This is a photocrosslinked/immobilized coating; an end-grafted homopolymer brush is not assumed.

## Working directory and setup

Use `/workspace/practice-md`, branch `md/cbaa-surface`. The original `/workspace/practice` checkout is a separate worktree and must not be switched, reset or modified for this project. This cloud task is already isolated; do not create more worktrees unless requested. The linked worktree's Git metadata lives inside the original checkout; preserve both directories in this cloud instance.

```bash
cd /workspace/practice-md
bash scripts/install.sh
source scripts/activate.sh
gmx --version
```

The tested setup is Debian 13, GROMACS 2025.2 CPU-only and RDKit 2025.9.6. Shared binaries, package cache and Python environment remain outside Git in `/workspace/md-tools`; these are installed tools, not research data. `MD_TOOLS_DIR` can select another installation location. `config/requirements-chem.lock` records exact tested Python dependency versions. The apt setup verifies Debian signatures/checksums; no system package installation or privileged writes are used. If the pinned Debian package disappears from the mirror, diagnose and update the version deliberately rather than bypassing verification.

## Reproduce and validate

```bash
source scripts/activate.sh
python scripts/fetch_inputs.py
python research/chemistry/build_candidates.py
python scripts/build_copolymer.py --dp 20 --seed 2026
bash scripts/smoke.sh
```

`fetch_inputs.py` verifies official Lin 2018/2020 supplementary PDF sizes and publisher checksums, and retains fibrinogen/FXII coordinate candidates. It does not retrieve the paywalled main article or all historical searches. Existing source files are reused; downloaded coordinate SHA-256 hashes are recorded in `research/input-inventory.json`. Public input downloads need Figshare/RCSB network access. Papers and raw inputs remain locally present but are ignored by Git; their sources are recorded in research notes/inventory. Do not distribute publisher PDFs as repository source.

The chemistry script generates monomer/capped-oligomer geometry candidates and validates connectivity, expected formal charges and geometry minimization. Its generic MMFF geometry optimization is **not** a production force field. Earlier DP100 sequence examples are retained as illustrative artifacts;their deleted generator is not part of the current reproduction commands. DP100 is not a measured chain length. Small generated SDF/PDB/JSON candidates are versioned. Source structure downloads and simulation output are excluded.

The water smoke check performs solvation, PME preprocessing, minimization, 100 steps of MD and energy-file validation. Fresh results go in ignored `runs/water-smoke.*`. This only validates GROMACS execution; it does not validate the coating, proteins, adsorption or photochemistry. It uses an Amber water setup solely as an engine check.

## New three-component model

The user approved directly constructing CBAA-1/BMA/BPAA at a 7:2:1 repeat-count ratio. `scripts/build_copolymer.py` builds an actual connected pre-UV chain; the default prototype uses DP20 (14/4/2), hydrogen end groups and seeded random backbone stereochemistry. These are explicit computational assumptions. See models/README.md for outputs and validation. The chemistry uses BPAA, not Nile blue. Identical formulation in a paper is not a prerequisite for building this model.

## PC/graft preparation and completed numerical pilots

The user reports PC grade unknown and accepts BPA-PC as the initial assumption. `scripts/build_pc_graft.py` builds finite OH-ended PC4 and one DP10-to-PC4 covalent-graft geometry; this is not a bulk slab or measured UV conversion. See models/pc-graft-README.md. Exploratory force-field assignment now uses GAFF2/AM1-BCC through ACPYPE2026.10.15 with optimized OpenBLAS, replacing the earlier tentative CHARMM route; future proteins must use compatible Amber conventions. PC4/COP10/GRAFT parameter-assignment checks and short surface/protein numerical controls completed. See CURRENT_MD_STATUS.md and research/DESIGN_REVIEW_2026-10-28.md for actual results and limits;inspect process state before reporting a running job.

## Project map and current limits

- `models/`: actual three-component graph/3D prototypes and assumptions.
- `config/system.json`: confirmed user target versus unresolved experimental inputs.
- `research/COATING.md`: current coating evidence and primary/secondary source distinctions.
- `research/chemistry/MODEL.md`: chemistry definitions, provenance and parameterization gates.
- `research/chemistry/`: small candidate structures, scripts and validation results.
- `research/NOTES.md`: earlier research findings; newer COATING.md supersedes silica assumptions.
- `scripts/`: installation, activation, source fetching and water smoke validation.
- `migration-sha256.json`: original research file hashes verified during migration, before documented edits.
- `runs/`: ignored outputs including preserved earlier setup tests/logs.

BMA is confirmed as n-butyl methacrylate in Lin 2018 SI. Uploaded Lin2020/2025 primary papers directly confirm BPAA as N-(4-benzoylphenyl)acrylamide; candidate connectivity matches. The target batch and its force-field parameters remain unvalidated. The main papers also contain fluorescent ternary copolymers (CBAA/BPAA/Nile blue in Lin2020 Figure1b; CBAA/BMA/Nile blue in Lin2025 Figure7a). An earlier binary-only description was incorrect; see research/PAPER-READTHROUGH.md. Those fluorescent systems do not establish the user's CBAA/BMA/BPAA7:2:1 formulation or polycarbonate substrate. Main texts are now available from user-uploaded PDFs (locally retained, Git-ignored); website access is no longer a blocker for reading these articles. Reference UV312nm/600mJ/cm² and primary synthesis data are recorded in research/COATING.md and config/system.json, separately from target experimental conditions.

Before production inference: evaluate explicit scenarios for PC grade,molecular weight/DP,film loading/thickness and crosslinks;the user UV protocol is recorded in research/USER_EXPERIMENT.md. retain the user-approved 7:2:1 repeat-count modeling assumption while distinguishing it from measured incorporation; derive compatible unreacted/reacted polymer and PC parameters; validate brush/coating hydration; then prepare each protein assembly including missing residues/glycans/disulfides. 8OS5 is only an N-terminal FXII fragment, not full-length FXII or FXIIa. Fixed-topology classical MD represents a post-UV network, not photochemistry or proteolytic activation. Coating/protein trajectories exist as numerical pilots;none establishes adsorption efficacy or FXII activation.
