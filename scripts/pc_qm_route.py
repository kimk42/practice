"""Explicit per-chain optimizer selection; shared lock prevents duplicate QM."""
import json
from pc_length_parameters import ROOT

CONFIG = ROOT/'config/pc-qm-continuation.json'


def route(dp):
    choice = json.loads(CONFIG.read_text())['routes'][str(dp)]
    if choice == 'native':
        return ROOT/'runs/pc-length-sqm-checkpointed-v2', 'run_pc_sqm_checkpointed.py'
    if choice == 'lbfgs':
        return ROOT/'runs/pc-length-sqm-lbfgs-v3', 'pc_sqm_lbfgs.py'
    if choice == 'reviewed':
        return ROOT/'runs/pc-length-sqm-reviewed-v4', 'pc_sqm_reviewed.py'
    raise ValueError('Unknown QM route; no automatic fallback')
