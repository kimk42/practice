# CBAA-1/BMA/BPAA coating on polycarbonate: GROMACS research

User target: UV-immobilized CBAA-1-co-BMA-co-BPAA coating on polycarbonate, reported composition **7:2:1**, followed by separate fibrinogen and Factor XII interaction studies. The ratio basis (molar/mass; feed/incorporated) remains unconfirmed. This is a photocrosslinked/immobilized coating; an end-grafted homopolymer brush is not assumed.

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
python research/chemistry/sequence_candidates.py
bash scripts/smoke.sh
```

`fetch_inputs.py` verifies official Lin 2018/2020 supplementary PDF sizes and publisher checksums, and retains fibrinogen/FXII coordinate candidates. It does not retrieve the paywalled main article or all historical searches. Existing source files are reused; downloaded coordinate SHA-256 hashes are recorded in `research/input-inventory.json`. Public input downloads need Figshare/RCSB network access. Papers and raw inputs remain locally present but are ignored by Git; their sources are recorded in research notes/inventory. Do not distribute publisher PDFs as repository source.

The chemistry script generates monomer/capped-oligomer geometry candidates and validates connectivity, expected formal charges and geometry minimization. Its generic MMFF geometry optimization is **not** a production force field. The sequence script generates five reproducible DP100 examples with exact 70/20/10 counts under an explicitly provisional incorporated-molar-ratio assumption; DP100 is not a measured chain length. Small generated SDF/PDB/JSON candidates are versioned. Source structure downloads and simulation output are excluded.

The water smoke check performs solvation, PME preprocessing, minimization, 100 steps of MD and energy-file validation. Fresh results go in ignored `runs/water-smoke.*`. This only validates GROMACS execution; it does not validate the coating, proteins, adsorption or photochemistry. It uses an Amber water setup solely as an engine check.

## Project map and current limits

- `config/system.json`: confirmed user target versus unresolved experimental inputs.
- `research/COATING.md`: current coating evidence and primary/secondary source distinctions.
- `research/chemistry/MODEL.md`: chemistry definitions, provenance and parameterization gates.
- `research/chemistry/`: small candidate structures, scripts and validation results.
- `research/NOTES.md`: earlier research findings; newer COATING.md supersedes silica assumptions.
- `scripts/`: installation, activation, source fetching and water smoke validation.
- `migration-sha256.json`: original research file hashes verified during migration, before documented edits.
- `runs/`: ignored outputs including preserved earlier setup tests/logs.

BMA is confirmed as n-butyl methacrylate in Lin 2018 SI. A 2022 review citing Lin 2020 identifies BPAA as N-(4-benzoylphenyl)acrylamide; exact target-batch identity still requires its primary synthesis record. Lin 2020 SI describes a CBAA/BPAA binary copolymer on PVC and does not establish the user's 7:2:1 terpolymer or polycarbonate substrate. ACS main text remains blocked by HTTP 403; Figshare SI downloads work.

Before production: confirm BPAA identity, PC grade, composition basis, molecular weight/DP, film loading/thickness, UV protocol and crosslink constraints; derive compatible unreacted/reacted polymer and PC parameters; validate brush/coating hydration; then prepare each protein assembly including missing residues/glycans/disulfides. 8OS5 is only an N-terminal FXII fragment, not full-length FXII or FXIIa. Fixed-topology classical MD represents a post-UV network, not photochemistry or proteolytic activation. No coating/protein trajectory has been generated.
