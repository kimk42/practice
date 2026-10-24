# Solvated coordinate initialization

Corrected, accepted patch: one COP10–PC4 covalent graft plus38 free PC4 oligomers,5461 solute atoms. The PC grade is unknown; BPA-PC and the short OH-ended substrate are explicit initial assumptions. The PC packing target is0.6g/cm³, not an equilibrated bulk density.

GROMACS editconf and solvate prepared a6×6×8nm periodic water box with7455 waters and27826 total atoms. GROMACS check read the coordinates successfully; atom counts, solvent count, box dimensions and finite coordinates were independently checked. Hashes link the solvated coordinates to the accepted dry patch in preparation.json.

No salt or protein has been added. This is not a minimized, equilibrated or validated interaction model. AM1-BCC charge calculations and production topology assembly remain pending. Water coordinates use the standard spc216 starting geometry; the gated Amber-compatible setup uses TIP3P parameters and must minimize/validate geometry before MD.

Reproduce coordinate generation from the project root after successful packing:

```bash
source scripts/activate.sh
gmx editconf -f models/pc-patch/pc4-grafted-cop10-patch.pdb \
  -o models/solvated-patch/dry-box.gro -box 6 6 8 -c
gmx solvate -cp models/solvated-patch/dry-box.gro -cs spc216.gro \
  -o models/solvated-patch/solvated.gro
gmx check -f models/solvated-patch/solvated.gro
```

The full gated workflow checks hashes before reusing the saved coordinates. Do not interpret this water-box file as a trajectory or a successful energy minimization.
