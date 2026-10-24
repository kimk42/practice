# CBAA-1/BMA/BPAA terpolymer prototypes

The user authorized constructing the three-component model directly, using a **7:2:1 repeat-count ratio**. This does not require an identical formulation in the reference papers. BPAA is N-(4-benzoylphenyl)acrylamide; Nile blue is not included.

```bash
source scripts/activate.sh
python scripts/build_copolymer.py --dp 20 --seed 2026
# Smaller exact-ratio prototype:
python scripts/build_copolymer.py --dp 10 --seed 2026
# Optional GROMACS coordinates; dry box is a format check, not a solvated system:
gmx editconf -f models/copolymer/cbaa1-bma-bpaa-dp20-seed2026.pdb \
  -o models/copolymer/cbaa1-bma-bpaa-dp20-seed2026.gro -c -d 1.2 -bt cubic
```

Default prototype: **20 repeats = 14 CBAA-1 + 4 n-butyl methacrylate + 2 BPAA**, in a reproducibly shuffled head-to-tail chain. The polymerizable C=C bonds become single backbone bonds; BMA retains its alpha-methyl group; BPAA retains its unreacted benzophenone ketone. The model has hydrogen-terminated backbone ends and randomly sampled backbone stereochemical tags. DP, end chemistry, sequence statistics and stereochemistry are modeling choices, not measured properties.

Outputs in `copolymer/`:
- `.sdf`: molecular connectivity, bond orders, formal charges, stereochemistry and coordinates; preferred chemistry input.
- `.pdb`: coordinates with repeat-specific CBA/BMA/BPA residues and atom names.
- `.gro`, when converted with GROMACS: coordinate-format compatibility only.
- `.json`: sequence, atom-to-repeat mapping, assumptions, validation results and atom indices of potential BPAA reaction sites.

The builder independently checks elemental balance, neutral total formal charge, expected zwitterions, remaining vinyl bonds, BPAA rings, backbone links and BMA substitution. Geometry uses ETKDG plus generic MMFF minimization and is round-tripped through SDF. These are chemistry/geometry validations, not production MD force-field validation. Generic MMFF energies describe only this preparation procedure and are not adsorption energies or free energies.

**This is a pre-UV isolated chain, not a finished PC-supported coating.** No partial charges, production atom types/torsions, surface bonds or crosslink fraction are assigned. The two BPAA sites are potential reaction sites; they are not two automatically reacted crosslinks. Next parameterize the unreacted and reacted components consistently, select an explicit PC repeat/substrate model, and construct candidate post-UV networks with documented attachment/crosslink assumptions. Larger chains, film loading and crosslink fractions require separate resource and sampling choices.

Earlier CBAA homotrimer and DP100 sequence-only files remain historical comparators in `research/chemistry/`; the actual three-component 3D prototype is here.
