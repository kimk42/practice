# Ten-stage execution: interim evidence

This report follows the latest ten steps requested by the user. Prepared inputs, short numerical checks and research conclusions are distinct. No complete ten-stage claim is made.

| Step | Work | Actual state |
|---|---|---|
| 1 | initial FXII/coating numerical analysis | passed numerical pilot only |
| 2 | coating relaxation/equilibration | 100ps preliminary control available; equilibrium unproven |
| 3 | surface model sensitivity | 0/3/6 graft coordinate candidates prepared and preprocessed; representative model unresolved |
| 4 | FXII protein-only control | 0.2ps passed;100ps queued |
| 5 | bare-PC versus coated-PC controls | both short initial pilots available |
| 6 | orientation/starting-state sampling | two accepted orientations;bare/coated PC;0.5/2nm gap subset;three velocity seeds queued;short numerical sensitivity only |
| 7 | replicated production MD | blocked: no validated production model or converged pilot budget |
| 8 | glycosylated FXII comparison | nine annotated sites and five citation records retrieved;glycoform/parameters unresolved |
| 9 | full fibrinogen simulation | Alpha/Gamma-A isoform sequence candidates aligned;955 construct residues unmodeled,1081 relative to this full-length candidate |
| 10 | integrated mechanism/uncertainty report | interim evidence report generated;mechanistic comparison not yet supported |

## Measured evidence

The small coating10ps run used27826 atoms, returned21 finite energy records and mean sampled temperature297.87K. Its initial z-density profile is in models/analysis/initial-coating-density.csv. The FXII water0.2ps pilot passed, with20 annotated disulfides retained. Neither fact establishes equilibration. Expanded/protein-surface result files are accepted only when their strict preparation/force/energy gates pass.

## Research gates

Full-length protein modeling and a representative coating remain essential. FGA canonical P02671 is866-residue Alpha-E; the alternative Alpha isoform and Gamma-A sequences were retrieved via EBI and aligned to3GHG constructs. The provisional six-chain sequence candidate has2994 residues versus1913 observed residues, leaving1081 absent coordinates relative to that candidate (955 is the construct-relative count). Fibrinopeptides are retained; only signal peptides are removed, because the target is fibrinogen rather than fibrin. Thus using the canonical sequence to fill3GHG gaps could select the wrong isoform. The glycosylation papers establish site evidence, but the retrieved annotations do not uniquely determine all nine site-specific native glycan structures/occupancies. No placeholder glycans or missing force-field terms are introduced.

The current host has5 CPU cores and no exposed GPU. The380158-atom FXII/coating0.2ps pilot measured0.366ns/day on4 threads. This short-run benchmark corresponds to roughly2.7days per1ns per run if sustained, not a validated long-run forecast. Long, replicated sampling requires model/relaxation validation and suitable compute allocation before production lengths are selected.

## Interpretation limits

No protein adsorption efficacy, adsorption free energy, physiological glycoform effect or proteolytic FXII activation has been demonstrated. Failure to observe adsorption in a short trajectory would not establish antifouling. The final mechanistic report remains gated on matched converged repeated observations.
