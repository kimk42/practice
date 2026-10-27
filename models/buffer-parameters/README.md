# cPBSzI parameter preparation

The assay recipe is config/experimental-buffer.json. This directory does not yet constitute a complete cPBSzI force field or matched-buffer simulation.

`scripts/prepare_buffer_ions.py` reads the installed Amber `frcmod.ionsjc_tip3p` (Joung–Cheatham2008 Ewald/TIP3P), records its SHA256/numerical values and converts Amber Rmin/2 and kcal/mol to GROMACS sigma/nm and kJ/mol. The three ions use distinct JC atom/molecule names, avoiding replacement of existing baseline Na/Cl definitions. Na/Cl/I form one consistent proposal. The original default ions.itp lacks an iodide molecule definition; an atom type named I alone is insufficient evidence of a usable iodide model. A neutral2Na/1Cl/1I water control passed GROMACS preprocessing and Fmax<1000 EM. This is import/numerical evidence,not new physical validation or the assay's concentration mixture.

`scripts/prepare_buffer_species.py` creates formal chemical graphs for citrate3−,one hydrogen-citrate2− microstate,H2PO4−,HPO4²− and azide−. No force-field partial charges,validated topologies,protonation fractions or production coordinates are supplied for those polyatomic species. Speciation,phosphate hydration/counterions and parameters remain explicit gates. The assay's approximate NaOH amount must not be modeled as the same amount of free hydroxide atpH7.4.

Ongoing neutralizing-ion numerical controls retain their original settings. A separately built matched-buffer model requires all gates above; partial ion support must not be advertised as complete cPBSzI MD.

## Latest bonded diagnostics

`scripts/audit_buffer_bonded.py` performs only typing/bonded lookup in ignored runs/. `--inspect-existing` rereads preserved outputs without rerunning jobs. The absence of parmchk2 missing-term flags is insufficient: azide emits an explicit assigned-bond-types warning and is blocked; citrate and phosphate show changed bond-order representations, which may involve resonance but require review before a charge/force-field model is approved. No diagnostic charge fields are exported for MD. The azide literature assessment supports a conditional pKa/speciation estimate only, not transferable atomistic parameters.

config/cpbszi-system-preparation.json records analytical integer-count alternatives for the existing18×12×18nm box. Phosphate hydrate cases remain separate, species populations and the macromolecular counterion policy are unselected, and no final Na count or runnable matched-buffer topology is asserted.

## Published azide route

The open SI of DOI10.1021/acs.jpcb.5c02025 was inspected,including section2 and references7–9. It uses LAMMPS/OPLS-AA in DMSO and refers separately to original azide bonded,charge andLJ models. The source chain is recorded in azide-published-model-route.json. The original publisher endpoints returned proxy403 on this runtime;this is an access failure,not proof of a paywall. No numeric parameters were guessed. Water-model/mixing-rule compatibility and aqueous validation remain gates before any Amber/TIP3P adoption. A different article using covalent azidoadamantane was excluded as an ion-parameter source.

## Computed charge candidates

The new candidate calculations exclude azide and are stored only under runs/. Initial requests for-j1 were automatically changed by Antechamber to-j4 for AM1-BCC;this explicit method warning and actual method are preserved in charge-candidate-audit.json. `scripts/review_buffer_charge_candidates.py` rereads existing outputs without recomputing. Citrate candidates pass numerical conservation/connectivity checks but need aqueous/charge-model validation. Both phosphate candidate outputs have a-0.001e net-charge residual;no arbitrary renormalization was applied and they remain blocked. The candidate generator now requests the actual supported-j4 method explicitly for any future fresh calculation. No MD topology was exported and no candidate is approved as an experimental-buffer model.

## Charge-output precision tracing

`scripts/audit_sqm_charge_precision.py` reads existing SQM/BCC outputs only. Both phosphate outputs sum to their3-decimal printed atomic-charge sum (-1.001/-2.001),while SQM reports the configured integer total (-1/-2). The pre-BCC sums match the printed sums and BCC preserves them. The inspected public upstream charge-print routine usesF14.3 for atomic values andF12.3 for totals,consistent with this loss of precision. This source is corroboration,not proven exact build provenance. No true extra digits were recovered,no charges renormalized and no executable modified. Higher-precision supported output or an independently validated charge route remains necessary.

The public-source network additions are saved only as an environment configuration draft;config/source-access-requirements.json describes the pending review/application. Do not treat a saved draft orunknown runtime policy as current connectivity.
