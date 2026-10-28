# PC length preparation: active calculations and revised packing

PC8 and PC16 whole-chain GAFF2/AM1-BCC assignments started independently at
2026-10-28 16:26 UTC. Each uses one BLAS thread and the existing 86400-second
ACPYPE limit. These are new finite-chain candidates, not an extension or reset
of any coating, buffer, or FXII trajectory. No new substrate MD has started.
Live process state and log growth must be checked before reporting these jobs
as still running; a saved `running` status alone is insufficient.

The input checks confirm one neutral connected molecule, two phenolic ends,
7/15 carbonate links, and 264/528 atoms for PC8/PC16. No PC4 terminal or repeat
charges were copied. Geometry seeds and assumed OH ends remain model choices,
not measured PC molecular weight or end chemistry. One charge conformation per
length is a limited starting point, not a conformationally averaged fit.

## Charge-conversion finding

The stricter audit of the preserved PC4 parameters identified an earlier
unreported ACPYPE conversion: its 132 printed AM1 charges sum to -0.004e, as
does the BCC MOL2, although SQM reports total charge -0.000 and uses qmcharge=0.
Installed ACPYPE2026.10.15 `balanceCharges` distributes this residual over atoms
tied for the largest absolute raw charge. The three carbonate carbon atoms
35,54,73 change from +0.772200 to +0.773533e; the serialized ITP total is
-0.000001e. This is converter charge balancing, not an intentionally charged
PC molecule. The installed converter source hash and exact differences are in
`models/surface-revision-v2/pc4-charge-conversion-audit.json`.

No old charge, topology, or trajectory was changed. Earlier numerical stability
results still describe the original ITP actually simulated. The new audit
reconstructs the converter operation rather than incorrectly requiring the raw
MOL2 and ITP charges to be identical. Nontrivial balancing is explicitly flagged
for review before new MD. PC8/16 raw residuals and redistribution must be compared
with PC4 after completion; this document does not preapprove them. Four regression
checks verify read-only behavior and rejection of charge corruption (even at
unchanged total charge), bond-order corruption, and missing-term flags.

## Dry packing design

Twelve cases cover PC4/PC8/PC16, starting densities 0.6/1.0g/cm3, and packing seeds
2026/2027. Every cell has 192 BPA units and 6336 atoms. These densities are
algorithmic starts, not measured or equilibrated PC densities.

| DP | Chains | OH ends | Carbonate links | Molecular mass sum, Da | Cubic side at 0.6/1.0g/cm3, nm |
| --- | ---: | ---: | ---: | ---: | --- |
| 4 | 48 | 96 | 144 | 47575.008 | 5.087357 / 4.290843 |
| 8 | 24 | 48 | 168 | 48198.864 | 5.109498 / 4.309517 |
| 16 | 12 | 24 | 180 | 48510.792 | 5.120496 / 4.318794 |

Initial volumes use RDKit standard atomic masses; density must be recomputed
from topology masses before interpreting MD. Dry preparation is cubic and periodic
in all three directions. The planned 36nm2 lateral area applies to the later slab
comparison, not to these bulk-preparation boxes. This avoids imposing a thin film
before assessing the dry substrate. Different packing seeds share the same
intrachain geometry at each DP and are not independent conformational ensembles.

The first preparation incorrectly constrained every atom inside the primary
cell, preventing chains from crossing periodic faces. That new packing worker
was stopped while preserving its restart and all intermediate files; neither
quantum job nor any prior MD was stopped. Revision2 restricts only atom1 of each
chain to the primary cell and retains periodic intermolecular distances. Its
first PC4 case reads the 48 saved rigid-body poses, with the restart hash recorded;
it does not reset that case randomly. The other eleven cases had not started.
Revision1 remains under `runs/pc-dry-packing-v1`; its plan is marked superseded.

Revision2 is `models/surface-revision-v2/dry-packing-plan-v2.json` and
`runs/pc-dry-packing-v2`. Each case must report Packmol success and pass an
independent element/count, periodic bond-geometry, interchain-distance and
chain/self-image check. The requested separation is 2.0A; a 0.002A allowance
accounts only for output PDB rounding. Failures retain their outputs and are
not retried automatically. No topology, density, or physical validation is
implied by a packing pass.

## Inspection and next gates

```bash
source scripts/activate.sh
ps -C sqm,packmol -o pid,etimes,pcpu,rss,stat
cat runs/pc-length-parameters-v1/pc8/status.json
cat runs/pc-length-parameters-v1/pc16/status.json
tail -20 runs/pc-dry-packing-v2/queue.log
python scripts/test_pc_parameter_audit.py -v
```

Do not relaunch the parameter commands into existing directories. On completion,
review SQM geometry/SCF convergence, converter charge redistribution, atom and
bond order, masses, LJ and bonded terms. Then perform strict GROMACS preprocessing,
minimization and bounded numerical controls in new versioned directories. Bulk
density/structure and finite-chain effects need separate time-window and
independent-conformation assessment; short room-temperature MD cannot establish
equilibration of a glassy polymer. Hydrated slab/coating work follows those checks.
Full FXII/cPBSzI remains pending the surface and buffer physical/speciation gates.

All 43 pre-existing checkpoint hashes were recorded and verified unchanged.
The last complete local recovery archive remains
`/workspace/md-backups/cbaa-surface-20261028T160117Z` (HEAD60af023); it predates
these active calculations. A new full archive is deferred until a consistent
idle point, without stopping calculations merely to back them up. Small scripts,
plans and audits are saved on `md/cbaa-surface`; raw outputs remain ignored.
