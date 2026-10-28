# Authorized next ten stages

This is the earlier ten-stage framework accepted by the user. The2026-10-28 design review and five-step revision supersede its execution priorities; Actual evidence/status is generated in research/TEN_STAGE_PROGRESS.md and models/analysis/ten-stage-status.json. Short numerical pilots are not relabeled as equilibration or production.

| Step | Work | Execution/gate |
|---|---|---|
| 1 | Analyze the FXII/coating initial check | Expanded EM and0.2ps coupled pilot passed; endpoint RMSD/disulfide/separation QA in models/analysis. |
| 2 | Equilibrate coating alone |100ps numerical control completed; density/hydration/structure and longer relaxation determine whether equilibrium is established. |
| 3 | Improve/test surface models |0/3/6 graft input matrix preprocessed; DP, film density/thickness and crosslink distribution remain unresolved physical-model gates. |
| 4 | Equilibrate FXII alone |0.2ps water check passed;100ps numerical control completed. Single predicted conformation/protonation remains provisional. |
| 5 | Compare bare and coated PC | Restore PC chemistry using original whole-PC4 parameters; matched bare-PC EM and corresponding FXII pilot completed. |
| 6 | Sample starting orientations/distances |Two accepted orientations across both surfaces;0.5/2nm gap subset and three velocity seeds,18 short numerical runs in total. These are not independent conformational production replicas. |
| 7 | Run replicated production |Gated on physically defined/equilibrated models and sampling/compute assessment. No unvalidated months-long CPU production run is launched blindly. |
| 8 | Compare glycosylated FXII |Nine annotated sites/five citation records retrieved; site-specific structures/occupancy and compatible glycan parameters/building must be resolved. |
| 9 | Extend to full fibrinogen |Alpha and Gamma-A isoform sequence candidates retrieved/aligned;2994 candidate residues versus1913 observed. No arbitrary coordinates supplied for1081 absent residues; a six-chain AF3 sequence-candidate input is prepared but not run. |
| 10 | Integrate mechanism/uncertainty |Interim evidence report generated now and after completed jobs; final adsorption/efficacy/free-energy/activation conclusions remain gated on appropriate repeated data. |

Runtime: scripts/advance_ten_stages.py serializes4-thread CPU jobs. scripts/continue_velocity_after_controls.py handles the queued short velocity tests after the active control worker; completed results are reused. scripts/report_ten_stages.py generates the interim report. Failures stop dependent work. EM coordinate snapshots are saved every100 steps; long MD checkpoints every minute. Environment shutdown can still stop execution, so process status must be checked before reporting that a job is active.

Current pilot conditions are300K, fixed volume and neutralizing ions only; user pH7.4 and cPBSzI recipe are documented;their complete MD mapping remains unapproved. The initial2nm separation and0.2ps pilots do not sample adsorption. Glycosylated/native fibrinogen models are not silently replaced by the unglycosylated FXII or incomplete3GHG core.
