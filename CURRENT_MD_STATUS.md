# Current MD status

Latest continuation 2026-10-30: PC16 resumed from its verified evaluation002345 after the environment lost the workers. PC8's separate bounded force-polish passed the unchanged RMS gate and independent tighter-SCF/rigid-rotation checks; its reviewed endpoint is undergoing strict AM1-BCC/GAFF2 postprocessing. Six of twelve packing cases pass geometry; all six high-density cases retain over1000 serious interchain clashes each and remain held. Both citrate and both phosphate water controls completed102ps and descriptive hydration/Na-distance checks;no physical approval follows. Full coating/FXII/cPBSzI remains unstarted. See [the report](research/CONTINUATION_2026-10-30.md), active processes, `runs/pc8-force-polish-review-v4/result.json`, and the citrate/phosphate continuation status files. The supervisor's last completed-loop snapshot may predate live parameter processing.

Runtime check 2026-10-29 15:39 UTC: PC16 remains live. PC8 stopped with an abnormal L-BFGS line search above the unchanged gradient gate; one reviewed continuation from verified accepted coordinates failed the same way and is held without further automatic retries. Five of twelve packing cases pass geometry checks, including the newly rechecked PC8 low-density seed2027 case. The other low-density PC8 case resumes saved poses with standard molecule relocation; six high-density cases remain held after packing failures. Two of twelve dry MD controls and both citrate 2ps controls are complete. All51 archived MD checkpoint hashes match. Full FXII/cPBSzI remains unstarted. See [this runtime check](research/EXECUTION_CHECK_2026-10-29.md).

Latest continuation2026-10-29: PC8/16 QM resumed from verified accepted coordinates after an environment interruption. Two isolated citrate candidates completed explicit final-gradient checks and neutral-water EM/2ps controls;the exported charges match their original candidates. Full coating/FXII/cPBSzI MD remains unstarted.

Packing correction: Packmol's automatic PBC constraint had still confined every atom to the primary box,despite the input selecting atom1. Prior v2-v4 claims that chains could cross the boundary were incorrect. A separately built21.2.3 anchor-only variant passed boundary-crossing and minimum-image tests against unmodified executables. Unfinished cases now use this variant in new attempt folders with saved poses and unchanged independent distance/bond checks. See [the latest report](research/CONTINUATION_2026-10-29.md) and `runs/pc-completion-status.json`;physical and equilibrium gates remain.

## Earlier continuation snapshots

Latest continuation2026-10-28: the environment became usable again. PC4 native AM1 optimization passed the explicit gradient gate; two6336-atom dry PC4 controls each completed EM,2ps NVT and20ps NPT. Independent energy/checkpoint checks confirm completion, but density increases by5.57–6.48% between the first/last5ps windows, so neither case is an equilibrated substrate. See [the continuation report](research/PC_CONTINUATION_2026-10-28.md) and `models/surface-revision-v2/pc4-dry-controls-interim-audit.json`.

PC8/16 retained completed native coordinate checkpoints after another environment interruption. Continuation now uses SQM analytic gradients with L-BFGS and exact accepted-coordinate storage; the first PC16 unscaled trial was rejected for gradient-print overflow, and a reviewed smaller optimizer step is selected. The unchanged physical RMS threshold is0.0005kcal/mol/A. Interrupted packing resumes in new attempt directories from saved poses with the unchanged independent distance gate. `runs/pc-completion-status.json`, per-chain status files and live processes are the runtime authorities. The supervisor advances accepted cases through strict parameter review and the same short dry controls; scientific failures are held. Full coating/FXII/cPBSzI MD is still not started.

All43 prior MD checkpoints and six additional PC4 checkpoint files were verified before continuation. The complete local backup remains `cbaa-surface-20261028T191757Z`; it predates these new results. Do not stop live work to create a full backup.

A verified selective snapshot at `/workspace/md-backups/pc-continuation-20261028T214651Z` adds the completed PC4 results and selected accepted QM/packing checkpoints without stopping live work; it is not a full simultaneous backup. Citrate's existing carboxylate serialization and equal oxygen-pair charges were reviewed without changing parameters. Its protonation microstates and aqueous accuracy remain unresolved; see research/CITRATE_REPRESENTATION_REVIEW_2026-10-28.md.

## Historical milestones

At19:15UTC, the original PC8/16 jobs had stopped without restart coordinates. Revision3 packing recovery was launched then; its optimizer stopping-floor failure was subsequently superseded by revision4's conservative generator margin. Those original outputs remain preserved. See research/EXECUTION_CHECK_2026-10-28.md for that earlier snapshot.

Previous launch: PC8/PC16 whole-chain AM1-BCC assignments started at2026-10-28 16:26UTC. The stronger PC4 audit records ACPYPE's existing -0.004e raw-charge residual redistribution;old parameters are unchanged. See research/PC_LENGTH_PREPARATION_2026-10-28.md for launch details and scientific gates.

Previously completed: design-review follow-up through five bounded stages. PC8/PC16 geometries prepared;water-referenced cPBSzI count alternatives replace full-box planning;two separately charge-projected phosphate candidates passed EM/2ps numerical controls;unchanged rigid azide continued from2ps to102ps with hydration/distance diagnostics. These are not physical forcefield or adsorption approvals. See research/FIVE_STEP_REVISION_2026-10-28.md and config/model-revision-policy.json. Earlier descriptions below are historical milestones.

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

## Next-step precision/network diagnosis

The phosphate0.001e deviation was traced read-only through printed3-decimal SQM atomic charges,pre-BCC input and BCC output. Integer QM totals are correct;no missing digits were invented or charges changed. The public upstream Fortran print routine corroborates fixed3-decimal output but is not claimed as the exact bundled build. Existing baseline MD remains complete with no active calculations. Original-parameter access is blocked at currently unlisted domains;an additive4-domain network draft and corrected startup guidance were saved for user review,not applied or published. Existing custom hosts,presets,repositories,credentials and install script were preserved.

## Formal buffer bookkeeping review

Enumerated3552 candidate population/macrocharge cases across three phosphate hydration hypotheses and formal macromolecule charges0/+6. Every case satisfies integer global charge balance with the conditional sodium count. This enumerates only the specified Cit3-/HCit2- and HPO4²-/H2PO4- basis, not pH7.4 speciation or a selected counterion policy. No partial charges were changed and no MD topology was exported. See models/buffer-parameters/formal-charge-balance-audit.json. Runtime revision6 still lacks the four requested source-access domains;the saved environment draft is not applied,so blocked original-source requests were not repeated. Existing checkpoints and completed runs remain untouched.

## Source access after environment revision8

The four requested hosts are now in the enforced runtime policy. OpenAlex and Semantic Scholar returned200. They identify a green-OA repository resolver for the1993 azide LJ paper:http://hdl.handle.net/11380/741371. Numeric parameters remain unread:the AIP PDFs returned HTTP403 and Elsevier returned tunnel403. Indexes disagree about OA status of the2011 bonded paper;neither index label proves retrievable full text. Saved an additive hdl.handle.net network draft;not runtime-applied. No MD started and no parameters changed. Request details:research/proteins/modeling-evidence/azide-primary-access-revision8.json.

## Repository access verified at revision17

The iris.unimore.it host is now allowed in the enforced runtime policy. A browser-user-agent request returned200 and the actual repository metadata was read. It explicitly states “Non ci sono file associati a questo prodotto.” (no attached files). The green-OA index record therefore did not yield an original PDF or numeric LJ values. No further domain addition is presently justified. Azide source-value/water-model compatibility and phosphate charge-precision/physical-validation gates remain open. No MD process was active and no calculation/checkpoint was reset. See research/proteins/modeling-evidence/azide-repository-access-revision17.json.

## Uploaded azide papers: numeric source gate resolved

User supplied1998/2011/2025 PDFs;originals and complete text extraction are preserved privately in ignored raw/. Inspected parameter/method passages and visually verified1998 p5515 and2011 p214504-7. Published terminal/central charges are -0.9215/+0.8430e (sum -1),equal NN distance0.11626nm,and all-N LJ sigma0.332nm/epsilon0.3550124kJ/mol (1998 conversion;2011 rounds to0.355). Sigma is explicitly identified in1998,so2011 word“radius” is not interpreted as half-sigma. Units and LJ landmarks passed arithmetic checks. Values from1993 are independently repeated by the two uploaded primary studies;that original PDF is not required to read these values.

Important correction:2011 is a geometry/model reference,not an established source of transferable harmonic azide force constants. The1998 bath simulation constrains linear solute geometry and includes field-dependent charges. The2025 DMSO SI says no bond/angle constraints but does not list a complete flexible azide bonded table in section2. No arbitrary harmonic constants or polarizable-to-fixed-charge equivalence was assumed. Next gate:choose/verify a nonsingular linear-ion GROMACS representation and assess aqueous Amber/TIP3P/JC compatibility. No topology exported or new MD started;phosphate and speciation gates also remain. See models/buffer-parameters/published-azide-candidate.json.

## Rigid azide numerical control completed

Built a separate fixed-reference-charge azide candidate using two noninteracting mass carriers and three interacting virtual N sites. Carrier positions/masses analytically preserve the three equal-N nuclei total mass,COM,and transverse inertia;one carrier-distance constraint avoids a singular collinear constraint triangle. GROMACS2025.2 local manual equation284 supports the linear virtual-site construction. No harmonic azide constants were invented. This implements an exploratory rigid fixed-charge approximation,not the1998 field-dependent-charge model.

A separate3nm cubic system with875 TIP3P waters,one JC Na andone azide (2631 coordinate sites including carriers) passed strict grompp without maxwarn,EM in13steps (Fmax818.001),and2ps NVT at300K. Forty-one saved frames had finite energy,NN distances within1.964e-7nm and NNN angles179.999899–179.999987degrees. No LINCS warning was detected. Mean sampled temperature306.946K,endpoint300.079K;this short initial relaxation is not an equilibrium claim. Full-precision TRR-derived G96 coordinates were used for geometry checks after six-significant-digit XVG proved insufficient;no tolerance was relaxed.

Preparation initially failed because atomtype directives followed a molecule definition;split type/molecule includes corrected the order without maxwarn. A recovery command used a wrong cwd and accidentally replaced the initial failed-grompp stdout file before archival;that failure is recorded in representation provenance and conversation. No MD/checkpoint existed then. Successful preparation/run logs,trajectory andcheckpoint are retained under runs/azide-rigid-water-control. All pre-existing coating/protein trajectories/checkpoints remained untouched.

Candidate ITP files and source/numerical audits are under models/buffer-parameters/azide-rigid-*. The physical model is NOT approved for full cPBSzI:hydration/ion-pair transfer to TIP3P/JC sodium and phosphate charge precision/speciation remain open. No full-buffer or new FXII calculation started.

## Five-step design revision completed

See research/FIVE_STEP_REVISION_2026-10-28.md for exact results,water-reference assumptions,charge projection,checkpoint continuation and remaining scientific limitations. New phosphate simulations and azide continuation have completed;no full cPBSzI/FXII or new surface MD has started.
