# Provisional BPA-PC and post-UV graft

The user reports **PC grade unknown** and accepts BPA-PC as the initial modeling assumption. `scripts/build_pc_graft.py` constructs a finite four-BPA-unit carbonate oligomer with phenolic OH ends, plus one candidate covalent graft to the DP10 7/2/1 copolymer. This is a local chemical model, not an equilibrated bulk PC film.

The phenolic terminal chemistry, oligomer DP, selected BPA methyl C-H attachment site and one reacted BPAA out of one in DP10 are explicit construction choices. They are not measured end-group chemistry, surface density or UV conversion. There is no inference that all experimental BPAA groups react.

The reaction candidate follows the source mechanism: reduce benzophenone C=O to C-OH, remove one H from a PC methyl group and form a new bond from that carbon to the former ketone carbon. Total elemental composition and neutral formal charge are conserved and the graft is one connected molecule. No arbitrary attachment force or unparameterized bond is added to a running GROMACS system.

```bash
source scripts/activate.sh
python scripts/build_pc_graft.py --pc-dp 4
```

Coordinates, connectivity and explicit assumptions are in `pc-graft/`. Force-field assignment must be done **after** the reaction, since changing a ketone and methyl to an alcohol and bonded methylene changes charges and bonded terms. Do not reuse the unreacted chain charges for this graft. A bulk slab requires packing and equilibration of longer PC chains and separate validation of density, roughness and end-group effects.

## Force-field route

GAFF2 / AM1-BCC via ACPYPE2026.10.15 (bundled AmberTools) is the executable **exploratory** route selected here. It is distinct from the previously considered CHARMM/CGenFF route. Future protein parameters and water/ion parameters must use compatible Amber conventions; do not combine GAFF topologies with CHARMM defaults. Automated analog terms and charge distributions require review before quantitative adsorption claims. Tool/topology generation or short stable trajectories alone do not validate zwitterionic-polymer hydration or protein-surface energetics.

## Parameterization execution

```bash
bash scripts/parameterize_candidate.sh models/copolymer/cbaa1-bma-bpaa-dp10-seed2026.mol COP10 runs/parameterization-optimized/cop10
bash scripts/parameterize_candidate.sh models/pc-graft/bpa-pc-dp4-phenolic-ends.mol PC4 runs/parameterization-optimized/pc4
bash scripts/parameterize_candidate.sh models/pc-graft/cop10-graft-bpa-pc4-one-site.mol GRAFT runs/parameterization-optimized/graft
```

New calculations now have a3-hour SQM stage limit (the original30-minute attempts timed out and were preserved) and must run in a new output directory; the script preserves existing outputs. It explicitly initializes the bundled AmberTools PATH so ACPYPE does not overwrite the optimized OpenBLAS library path, verifies dynamic-library resolution, and limits each job to one BLAS/OMP thread. Prior reference-BLAS attempts are retained in runs/parameterization and runs/parameterization-openblas as interrupted diagnostics, not completed charge assignments.

The AM1-BCC calculations are not equivalent to the earlier generic MMFF geometry preparation. Retain the SQM convergence output, ACPYPE logs and generated frcmod files. Check parameter-review warnings, total charges and atom mapping before topology assembly or GROMACS startup. No missing/zero force constants may be filled by hand just to make the simulation run.

`python scripts/pack_pc_patch.py` prepares a small finite PC patch around the trial graft, with a provisional0.6g/cm³ **low-density starting** PC mass/region-volume packing target. Earlier1.2g/cm³ and0.6g/cm³ attempts used an incorrect fixed-graft translation and are preserved as diagnostics. Their termination does not establish that a particular material density is infeasible. The lower-density initializer is not a replacement for subsequent density/film equilibration. The measured experimental density and film thickness remain unknown; this is an initialization choice, not evidence of equilibration. Packing must report success before the patch is used. Its output is not a production trajectory.

## Gated follow-up validation

```bash
source scripts/activate.sh
python scripts/surface_validation.py             # read-only inspection of prerequisites
python scripts/surface_validation.py --watch --max-wait 1800
```

Use --run for a one-shot execution after prerequisites are ready. A lock prevents duplicate watcher jobs. The watcher waits for PC4/GRAFT topologies and a successfully packed patch, then audits atom counts, charge totals, atom-index connectivity against the molecular inputs, required frcmod files and review flags. It stops on unresolved parameter warnings or conflicting atom types. Only after those checks does it prepare an Amber-compatible TIP3P-water system, demand minimization toFmax<1000kJ/mol/nm and attempt a400-step0.2ps initialization test at300K. No warning override is passed to grompp. This is a prerequisite check, not a physical validation of polymer hydration, bulk-PC density, crosslink fraction or adsorption. The MD branch of the watcher is not considered verified until actual result files report its completion.

Read `runs/surface-validation-status.json` and the source logs under `runs/parameterization-optimized/`. Long-running processes exist only in the current runtime; a snapshot does not preserve them. A failed or pending calculation is not a completed force-field assignment. Check outputs before restarting in a new directory; do not overwrite existing calculation evidence.


## Fixed-placement correction and continuation
The original input double-translated an already positioned graft: the output displaced it by(30,30,21.685)Å. Packmol's fixed-coordinate handling was checked against its official source. The corrected input is `fixed 0. 0. 0. 0. 0. 0.` because the input PDB is already in the desired absolute position. A separate actual Packmol regression run preserved all graft coordinates exactly (maximum displacement0Å).

The timed-out0.6g/cm³ attempt is retained in runs/pc-packing-shifted-diagnostic. Its38 free-PC poses were recovered by rigid fits, with maximum reconstruction error0.000597Å, and converted to Packmol's documented restart Euler format, with independent reconstruction checks. The corrected packing was resumed without loosening the2Å distance requirement:

```bash
source scripts/activate.sh
python scripts/pack_pc_patch.py --resume runs/pc-packing-shifted-diagnostic/recovered.restart --timeout 1800
```

New runs write `runs/pc-packing/packing.restart` at coordinate checkpoints. `--timeout` controls only one Packmol subprocess, not the cloud lifetime or MD duration. The initial600s cap was a script choice; default is now10800s (3hours). Do not rerun into an existing active packing directory. The actual full-patch success and later MD remain pending until their result files say otherwise.


## Three-stage continuation requested by the user
`python scripts/continue_three_stages.py` supervises the three stages: complete/audit the COP10,PC4,GRAFT force-field assignments; accept the corrected packing and solvate; then run force-converged minimization and the short initialization MD. Progress is in `runs/three-stage-status.json`; active parameterization directories are atomically recorded in `runs/workflow-active-parameters.json`. The script performs one fresh-directory retry only for a timeout, with the same scientific method/thresholds and a3-hour SQM limit; other calculation/parameter errors stop the workflow. A packing timeout can resume once from its saved checkpoint in a fresh run directory with a3-hour limit. The supervisor itself has a4-hour wait limit. These are resource guards, not simulation durations or evidence of completion. Do not start duplicate supervisors or overwrite active directories.
