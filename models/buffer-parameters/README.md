# cPBSzI parameter preparation

The assay recipe is config/experimental-buffer.json. This directory does not yet constitute a complete cPBSzI force field or matched-buffer simulation.

`scripts/prepare_buffer_ions.py` reads the installed Amber `frcmod.ionsjc_tip3p` (Joung–Cheatham2008 Ewald/TIP3P), records its SHA256/numerical values and converts Amber Rmin/2 and kcal/mol to GROMACS sigma/nm and kJ/mol. The three ions use distinct JC atom/molecule names, avoiding replacement of existing baseline Na/Cl definitions. Na/Cl/I form one consistent proposal. The original default ions.itp lacks an iodide molecule definition; an atom type named I alone is insufficient evidence of a usable iodide model. A neutral2Na/1Cl/1I water control passed GROMACS preprocessing and Fmax<1000 EM. This is import/numerical evidence,not new physical validation or the assay's concentration mixture.

`scripts/prepare_buffer_species.py` creates formal chemical graphs for citrate3−,one hydrogen-citrate2− microstate,H2PO4−,HPO4²− and azide−. No force-field partial charges,validated topologies,protonation fractions or production coordinates are supplied for those polyatomic species. Speciation,phosphate hydration/counterions and parameters remain explicit gates. The assay's approximate NaOH amount must not be modeled as the same amount of free hydroxide atpH7.4.

Ongoing neutralizing-ion numerical controls retain their original settings. A separately built matched-buffer model requires all gates above; partial ion support must not be advertised as complete cPBSzI MD.

## Latest bonded diagnostics

`scripts/audit_buffer_bonded.py` performs only typing/bonded lookup in ignored runs/. `--inspect-existing` rereads preserved outputs without rerunning jobs. The absence of parmchk2 missing-term flags is insufficient: azide emits an explicit assigned-bond-types warning and is blocked; citrate and phosphate show changed bond-order representations, which may involve resonance but require review before a charge/force-field model is approved. No diagnostic charge fields are exported for MD. The azide literature assessment supports a conditional pKa/speciation estimate only, not transferable atomistic parameters.

config/cpbszi-system-preparation.json records analytical integer-count alternatives for the existing18×12×18nm box. Phosphate hydrate cases remain separate, species populations and the macromolecular counterion policy are unselected, and no final Na count or runnable matched-buffer topology is asserted.
