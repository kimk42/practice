# Bonded coverage audit, not a ready force field

`graft-bonded-coverage.frcmod` was generated with the installed AmberTools atom typing and parmchk2, without quantum charge assignment. It checks bonded-term availability for the445-atom COP10–PC4 graft while the full AM1-BCC calculations continue. It must **not** be used as a complete MD topology or as proof of reliable adsorption energetics.

No `ATTN`/`need revision` unassigned-term markers were found. Several torsions/impropers are analog/general terms; reported penalty scores reach49.6, including carbonate/aromatic impropers. These are recorded rather than hidden. Zero torsional amplitudes supplied by the force field are not manually inserted missing parameters. Physical torsional/hydration validation is still required before production claims.

The source run is retained at `runs/bonded-coverage/graft/`; generated typing-only charges are not assigned AM1-BCC charges and are not supplied to the surface-validation workflow.
