# Protein preparation evidence

Run `source scripts/activate.sh`, then `python scripts/fetch_fxii.py --include-computed`, `python scripts/fetch_fibrinogen.py`, and `python scripts/audit_protein_inputs.py`. Downloads are public; source hashes are in the download manifests and audits. Raw coordinates/API files remain ignored. Install uses pinned gemmi0.7.3.

## FXII

The derived PDB is the AlphaFoldDB v6 predicted human mature zymogen, residues20–615; it is not an experimental conformation. The official615-residue sequence and checksum were verified before removing signal peptide1–19. No cleavage at the activation site was introduced. No missing amino-acid coordinates were invented. Glycans are absent despite nine UniProt glycosylation annotations. Approximately28% of mature residues have pLDDT below70; interdomain orientation and flexible regions require assessment using the downloaded PAE and suitable ensembles.

`python scripts/validate_fxii_control.py` creates an explicitly unglycosylated protein-only numerical pilot in ignored runs/fxii-unglycosylated-topology. It checks that pdb2gmx creates exactly the20 annotated disulfide pairs, then solvates, neutralizes, minimizes and attempts0.2ps at300K. Amber99sb-ildn/TIP3P and automatic standard protonation are provisional pilot assumptions. The cutoff electrostatics TPR is used only for ion placement; EM/MD use PME. Successful preparation/initialization does not validate a physiological glycoform, pH, conformation, equilibrium or adsorption. Existing dynamics are preserved rather than overwritten.

## Fibrinogen

RCSB3GHG assembly1 contains the six protein chains, with1913 observed amino-acid residues out of2868 declared construct residues. Therefore955 construct residues have no observed coordinates; full native sequence coverage would require additional canonical-sequence checking. The large missing alpha-chain regions cannot be repaired as a few ordinary loops. Bound peptides, carbohydrate and ions remain in the raw assembly for explicit handling. No truncated core was passed off as full fibrinogen and no arbitrary completion was used to initiate MD.

The observed coordinate extents are approximately43.5×7.1×16.6nm. Even this incomplete core cannot be simulated on the current6×6nm patch without major periodic self-interaction/geometry issues. A rectangular water-only bounding-box estimate with1.2nm padding each side gives roughly8300nm³, or830000 solvent atoms before solute, excluded-volume, substrate and interface geometry corrections. This is a rough resource estimate, not a validated box. Protein orientation can change box requirements. FXII's current candidate extents are about10.1×8.7×9.0nm, also larger than the patch.

## Required before physiological comparison

Select an explicit native glycoform/ensemble and protonation model for FXII, and a justified full-length or explicitly scoped fibrinogen construct. Build representative coating area, density, chain-length/crosslink distributions and hydration before comparing proteins. Protein-only equilibration, matched bare-PC controls, orientation sampling, independent replicas and convergence checks precede adsorption conclusions. The present CPU-only environment and short pilots cannot establish adsorption free energies or FXII activation.
