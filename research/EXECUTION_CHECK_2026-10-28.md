# Execution check and packing recovery

At2026-10-28 19:15UTC no SQM,Antechamber,Packmol or GROMACS process remained.
The container init process had started about100seconds earlier. This supports
an environment restart, but does not prove the exact termination cause or time.
The earlier `running` JSON records were stale; original copies were retained
before their states were changed to `interrupted_preserved`.

PC8's final recorded optimization iteration is440, RMS gradient0.0128kcal/mol/A;
PC16's is60, gradient0.3534kcal/mol/A. Neither output contains geometry convergence
or calculation completion. No intermediate-coordinate/restart file was found in
either parameter directory. The energy/gradient log cannot reconstruct atomic
coordinates. No charge assignment passed, and neither quantum calculation was
silently restarted from its initial geometry. Continuing that route requires a
new optimization and a restart-capable execution design; it is not a checkpoint
continuation of these440/60 iterations.

Packing revision2 reached four PC4 cases. Two0.6g/cm3 cases reported Packmol
success but failed the independent output-distance criterion (minimum distances
1.997411825 and1.997387043A; required2.0A with0.002A PDB serialization allowance).
One1.0g/cm3 case reached its one-hour limit; the other was interrupted. All four
have complete48-line rigid-body restart files. Eight PC8/16 cases had not started.

The optimizer-success/output-audit discrepancy is explained by Packmol21 source:
`getinp.f90` defaults `precision` to0.01; `fparc.f90` measures maximum violation as
target_distance_squared minus actual_distance_squared. Thus a2A target can accept
about1.9975A before PDB rounding. The independent audit correctly rejected these
outputs. It has not been relaxed.

`scripts/recover_pc_packing.py` prepares revision3 with precision0.0001 and the
unchanged2.0A target/audit. Four cases use their saved poses; eight remain new,
previously unstarted cases. Near-converged0.6g/cm3 cases use local refinement with
`disable_movebad` and1000 inner iterations. All previous files are retained.
The queue was launched in `runs/pc-dry-packing-v3`; inspect its live process/logs
before reporting current activity. Quantum calculations and MD are not launched
by this worker. A packing pass would still establish only starting geometry.

All43 pre-existing MD checkpoint hashes match the recorded values. A new full
local archive completed and verified at
`/workspace/md-backups/cbaa-surface-20261028T191757Z`:1,551,059,854 bytes,43checkpoints,
repository HEAD9ff4e18. It includes the preserved interrupted outputs and corrected
runtime states, and predates revision3 recovery. This remains an in-environment
backup, not an external copy. No new substrate or FXII/cPBSzI MD has started.
