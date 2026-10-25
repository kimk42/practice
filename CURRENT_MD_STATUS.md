# Current MD workflow status

Completed: connected7:2:1 copolymer models, a chemically explicit single-PC/BPAA graft candidate, corrected5461-atom low-density PC/coating packing and27826-atom water-solvated prototype. These are initial assumptions, not an experimentally calibrated coating.

PC4 GAFF2/AM1-BCC parameterization and topology/connectivity/charge audit passed. A PC4-only water control passed EM and400 steps of0.0005ps MD at300K;21 energy frames cover0–0.2ps. The small durable result is models/validation/pc4-water-control.json. It establishes numerical initialization only.

The larger COP10 and GRAFT quantum parameterizations remain active in directories recorded by runs/workflow-active-parameters.json. The controller writes runs/three-stage-status.json and will audit all parameters, then execute combined-surface EM/initialization only after they pass. It will stop on parameter or force failures; it does not invent fallback charges or bypass warnings. Running processes depend on the current environment remaining alive. Restarting the environment does not resume them automatically.

Both requested protein references are downloaded and audited. An unglycosylated predicted mature FXII protein-only pilot is separately running; its authoritative eventual outcome is runs/fxii-unglycosylated-topology/result.json or failure logs. Fibrinogen3GHG has extensive unmodeled regions and both proteins exceed the current patch dimensions. See models/proteins/README.md. Neither protein–surface production simulation nor full equilibration has been completed.

The ten-stage research plan remains PLAN_10_STEPS.md. Quantitative adsorption and mechanistic conclusions require resolved protein/coating models, larger systems, suitable compute resources, replicas and convergence. A short initialization success is not substituted for those stages.

## Interrupted-run recovery

The previously reported jobs were no longer alive on inspection. `scripts/recover_interrupted_workflow.py` restarted COP10/GRAFT in fresh directories, retaining original logs and the completed PC4 result. No usable optimizer-coordinate restart was present, so these are repeated AM1-BCC optimizations, not continuation from the last reported iteration. Per-job timeout is now12hours; the controller waits up to13hours before stopping. The manifest is runs/recovery-manifest.json. The FXII minimization was also restarted from its original TPR under em-recovery, preserving interrupted logs. Neither increasing a timeout nor launching recovery establishes convergence. These jobs still depend on the current execution environment staying alive.

## Live GRAFT deadline extension

The live GRAFT SQM process was preserved with the same PID and advanced to iteration290. Total execution allowance is now24hours from its original start. ACPYPE uses an in-process SIGALRM which cannot be changed by editing configuration, so scripts/extend_live_graft.py paused only the ACPYPE wrapper and replaced its watchdog; Antechamber and SQM continue. runs/live-extension-status.json and live-extension-error.json (if present) are authoritative for this handoff; the original three-stage status is no longer refreshed until topology handoff finishes. On successful completion, the new wrapper reuses the original finished AM1-BCC mol2, verifies its SHA256 is unchanged, and launches the normal parameter audits/surface initialization. Original wrapper/output logs remain preserved. New ordinary jobs also default to24hours. The live-extension script records incident-specific PIDs and must not be reused blindly. This extension is conditional on the execution environment remaining alive and does not guarantee convergence.
