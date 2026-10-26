# Current MD status

## Completed

- Connected7:2:1 copolymer, explicit single BPAA–PC graft, corrected5461-atom low-density PC/coating patch, and27826-atom water-solvated prototype.
- COP10, PC4 and GRAFT GAFF2/AM1-BCC parameters all passed graph/order, atom count, charge, mass and missing-term-flag gates. See models/validation/parameter-gates.json. These checks do not establish physical accuracy of analog torsions or the coating model.
- Combined coating EM converged to Fmax below1000kJ/mol/nm in734 steps, followed by400 steps of0.0005ps MD at300K. See models/validation/surface-initialization.json and ignored runs/surface-validation.
- Separate PC4 and predicted unglycosylated mature FXII water pilots each passed initial EM/0.2ps MD. See models/validation/. FXII has20 verified annotated disulfides; no native glycoform or equilibrium validation is claimed.

## Follow-up completed

The10ps fixed-volume surface-only follow-up completed in runs/surface-nvt-10ps:27826 atoms,21 finite energy records, mean sampled temperature297.87K and no LINCS/constraint warnings detected. See models/validation/surface-nvt-10ps.json. It continued the initial state/velocities; no protein or pressure coupling was added. This remains a short numerical check, not equilibrium validation.

## Recovery details

GRAFT SQM completed after its original12-hour limit was extended to24hours without restarting it. Pausing the original ACPYPE wrapper caused its unread stdout pipe to block Antechamber postprocessing. The wrapper/watchdog were stopped after preserving the completed SQM output. Standard Antechamber `-fi sqmout -c bcc -at gaff2` applied BCC from the converged AM1 output without another quantum optimization; all445 original atoms and graph edges matched, and the topology gate passed. The active manifest points to runs/graft-completed-sqm-08tfya1h, which contains source hashes and logs. scripts/finalize_completed_sqm.py supports this recovery method. The incident-specific live-extension script is disabled; ordinary new jobs default to24hours.

## Before protein–surface production

The present0.6g/cm³ starting PC patch, DP10 polymer, short OH-ended PC4 and single graft are provisional. Coating density, film loading/thickness, DP and graft/crosslink distributions need representative modeling. Both protein candidates exceed the6×6nm patch. Fibrinogen3GHG lacks955 declared-construct residues; FXII requires explicit glycoform, protonation and conformational-ensemble assessment. See models/proteins/README.md and PLAN_10_STEPS.md. No protein–surface adsorption MD, replicated production or adsorption free-energy result has been completed.

Execution processes require the environment to remain alive; saved files alone do not resume interrupted jobs.

## Expanded surface and authorized FXII pilot

The user authorized an unglycosylated FXII–surface placement/numerical pilot. An18×12×18nm surface candidate with380103 atoms is now minimizing. It comprises6 correlated translated copies of the small coating candidate; it does not establish representative morphology. Three FXII orientations were prepared; two satisfy lateral-padding criteria. The pilot worker waits for expanded EM, then places the first accepted pose, neutralizes and attempts strict EM plus0.2ps MD. See research/FXII_SURFACE_PILOT.md and runs/fxii-surface-pilot-status.json. No protein–surface run success is claimed before that worker completes.


## Ten-stage autonomous workflow

The user authorized continuation through the latest ten-stage plan. Expanded EM passed in2862 steps and the380158-atom unglycosylated mature FXII/coating0.2ps pilot passed. Endpoint QA is in models/analysis/initial-fxii-coating-QA.json; the initial2nm gap remains, so this is not adsorption evidence. Matched bare-PC EM, accepted orientation pilots, short0.5/2nm gap/velocity sensitivity tests and100ps preliminary coating/protein controls are queued sequentially. Runtime status is runs/ten-stage-execution.json and runs/velocity-queue-status.json. See PLAN_10_STEPS.md and research/TEN_STAGE_PROGRESS.md for actual scientific gates; stages7–9 and the final mechanistic interpretation are not claimed complete.

The current5-core host exposes no GPU; the short coupled pilot measured0.366ns/day on4 threads. Production lengths are deferred until model/relaxation and sampling/resource evidence support them. Glycosylation citations and explicit fibrinogen Alpha/Gamma-A sequence candidates were retrieved. The candidate six-chain sequence has2994 residues;1081 lack3GHG observed coordinates, compared with955 missing within the smaller declared construct. No missing loops/glycans/force-field terms were fabricated.

Matched bare-PC EM and the corresponding FXII0.2ps pilot also passed. The worker has advanced to the second accepted orientation; longer preliminary controls and short velocity/gap sensitivity remain queued. These two initial successes do not establish a comparative adsorption effect.
