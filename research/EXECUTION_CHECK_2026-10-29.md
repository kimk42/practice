# Runtime check: 2026-10-29 15:37 UTC, updated 15:39 UTC

Live processes, recent logs and immutable evidence were checked together. The
PC16 SQM worker and existing supervisor remained active; neither was restarted.
PC8 had stopped with L-BFGS-B status2 (`ABNORMAL`), not a timeout or successful
gradient convergence. Its last accepted evaluation008140 has RMS gradient
0.0020228355 kcal/mol/A, above the unchanged0.0005 criterion.

The accepted coordinates exactly match their saved NumPy array; raw-output and
coordinate hashes, atom ordering and bond geometry passed reinspection. The old
stop/status/settings/checkpoint records were copied into a timestamped recovery
directory. One bounded continuation was authorized from these same coordinates
with fresh L-BFGS memory and unchanged AM1, SCF, optimizer and gradient settings.
The existing supervisor launched it under the existing per-chain lock. A further
failure remains held; no automatic repeated retry or criterion relaxation was
introduced. Late failed trials differ from the accepted coordinates by less than
1e-13 A but exhibit electronic/numerical energy variability; the underlying cause
has not been established.

Follow-up: this single recovery attempt also stopped with status2 before accepting
a new optimization step. Evaluation008178 exactly reproduces the old accepted
coordinates, raw SQM-output hash and gradient. PC8 is therefore held again with
no further retry scheduled. Energy/gradient numerical diagnostics are needed
before selecting another optimization route. PC16 and the reviewed packing
continuation remain live; they are independent of this PC8 stop.

Five of twelve packing cases pass independent geometry checks. The new
PC8/rho0.6/seed2027 result was independently recomputed: minimum interchain
distance2.000938 A, nearest own-image distance24.914203 A and maximum bond-length
serialization difference0.002122 A. This is a6336-atom starting cell, not an
equilibrated substrate. The two old PC4 low-density starts retain the previously
documented whole-molecule confinement bias.

PC8/rho0.6/seed2026 timed out while an inherited `disable_movebad` option held
its placement in a local packing minimum. Saved poses were retained and one
reviewed continuation restores standard Packmol molecule relocation, as already
validated for the other low-density PC8 seed. The box, seed, validated generator
binary and independent distance/bond criteria are unchanged. The six high-density
cases failed with normal relocation already enabled and remain held for further
geometry-method review. No failed packing was passed to MD.

The same two PC4 dry EM/2ps NVT/20ps NPT controls and two isolated citrate
EM/2ps controls remain complete. No new dry MD control has completed since the
earlier report. All43 legacy, six additional PC4 and two citrate checkpoint
hashes match their archived references. Neither numerical completion nor a
packing pass establishes density equilibrium, aqueous accuracy or adsorption.

The runtime snapshot is
`models/validation/execution-check-20261029T1537Z.json`; the two explicit recovery
reviews are `models/surface-revision-v2/pc8-linesearch-recovery-20261029.json`
and `models/surface-revision-v2/pc8-seed2026-relocation-review-20261029.json`.
Live process/lock state and `runs/pc-completion-status.json` supersede this dated
snapshot. Full coating/FXII/cPBSzI remains unstarted while substrate and buffer
physical/speciation gates remain. No additional user input is required for the
currently resumed preparation work.
