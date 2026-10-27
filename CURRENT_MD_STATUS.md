# Current MD status

Latest verified state: both100ps controls and18 short pilots completed. FXII100ps backbone/SG geometry QA passed the stated numerical checks,without claiming native/equilibrium validation. Separate checkpointed10ps near-interface continuations for coated/bare PC are running sequentially;see runs/paired-interface-followup-status.json. Full cPBSzI is not approved: azide numeric source/compatibility and polyatomic-charge/speciation gates remain. Earlier in-progress descriptions below are retained as historical milestones.

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

## User experiment conditions received

Synthesis ratio7:2:1 is confirmed molar. The user reports0.5wt% aqueous coating solution,30min contact,drain/dry,9000mJ/cm² with three UV exposures using the cited-paper wavelength, and protein exposure for2h at room temperature,pH7.4. Buffer spelling `cPBSzI`,dose basis and meaning of `1% plasma concentration` await clarification. See research/USER_EXPERIMENT.md. Current300K/neutralizing-ion numerical controls continue as a provisional baseline; they are not relabeled as matched-buffer/2h adsorption simulations.

## Uploaded SOP resolved assay details

The uploaded SOP identifies cPBSzI1x containing NaCl,NaI,citrate,phosphate and azide;see config/experimental-buffer.json. User confirms9000mJ/cm² per exposure×3 and purified-protein1% plasma-reference concentration. Main SOP Fbg maps to0.03mg/mL;late0.3mg/mL Fbg-only entry is flagged as inconsistent. FXII concentration is absent. User room temperature overrides the SOP37°C reference. Existing neutralizing-ion numerical controls are not relabeled as matched-buffer simulations.

## Matched-buffer preparation progress

Published Joung–Cheatham/TIP3P Na/Cl/I parameters were converted with source hashes; a small neutral ion/water GROMACS import/EM control passed. Formal citrate/phosphate/azide graphs are prepared without invented partial charges or speciation fractions. Full cPBSzI MD remains gated on polyatomic parameters and charge/speciation accounting. A primary study's375nM physiological human-FXII reference gives a provisional1%=3.75nM assay mapping;this is not user-measured concentration. Existing numerical controls continue unchanged.

## Latest checkpoint-safe continuation

Coating100ps follow-up completed:201 energy frames,mean sampled temperature300.0525K;density profile and trajectory/topology hashes are saved in models/analysis/coating100ps-density*. FXII100ps water follow-up is alive with saved checkpoints;no restart or duplicate worker was launched. Buffer-only diagnostic typing completed without charge assignment: azide produced an explicit bond-typing warning; citrate/phosphate changed bond-order representation. These findings block full-buffer topology approval. Analytical molecule-count alternatives and a charge-balance constraint are prepared in config/cpbszi-system-preparation.json without choosing unmeasured hydrate/speciation fractions. No new mandatory user input is needed for the running baseline/continued parameter research.


## Completed controls and next paired baseline continuation

Both100ps controls and all18 short orientation/gap/velocity pilots completed. The18 energy traces/seeds/final files were independently reread without restarting them. FXII water control mean temperature299.568K(second-half300.012K);backbone endpoint RMSD0.2415nm(max0.2594nm),Rg3.1749→3.2262nm,and no annotated disulfide geometry flags. These are short numerical structural metrics,not equilibrium/native validation. See models/analysis/fxii-water-100ps-QA.json.

The accepted0.5nm-gap coated/bare interface states now have checkpointed10ps NVT continuations,one at a time,without reinitializing velocities or changing their source checkpoints. Status is runs/paired-interface-followup-status.json. These remain the legacy300K/neutralizing-ion baseline,not matched cPBSzI or adsorption evidence.

For azide,the reviewed2025 SI identifies three original bonded/charge/LJ sources. Original numeric-value URLs returned proxy403;no purchase/access bypass was attempted. Published-model details and water/mixing-rule compatibility remain unverified,so no replacement azide parameters were exported. Source route is models/buffer-parameters/azide-published-model-route.json. No new user input is mandatory for ongoing baseline continuations or further parameter research.


## Latest completed paired baseline and buffer charge-candidate review

Both10ps near-interface follow-ups completed and were analyzed without restarting them. Only three structural frames were saved per case(0,5,10ps);none show heavy-atom contacts below0.45nm. Sampled minimum distances are1.078nm(coated) and1.322nm(bare). Sparse-frame noncontact is not evidence of antifouling or an adsorption probability. See models/analysis/paired-interface-10ps-summary.json.

Four buffer AM1-BCC/GAFF2 charge candidates were computed in separate ignored directories;azide was excluded. Antechamber explicitly changed the requested prediction mode1 to4. That method override was retained and reviewed in the audit;the assumption of mode1 bond preservation was not fulfilled. The same preserved outputs were inspected without rerunning:two citrate candidates pass numerical charge/connectivity checks but are not physically approved;phosphate candidates have a-0.001e net-charge residual and remain blocked without ad-hoc correction. No candidate topology was exported or used as full cPBSzI MD. Source access and water/speciation/ion-behavior validation remain unresolved.

All repository changes remain local on md/cbaa-surface;no GitHub push,PR or new publication was performed.
