"""Centralised directory and file paths, all derived from the repository root."""

import os
from pathlib import Path


CONFIG_FILE_DIR = Path(__file__).parent.absolute()
ROOT_DIR = CONFIG_FILE_DIR.parent.parent

FIGURE_PATH = os.path.join(ROOT_DIR, 'Figs')
DATA_DIR = os.path.join(ROOT_DIR, 'data')
OBSERVATION_DIR = os.path.join(ROOT_DIR, 'Observations')
SIMULATION_DIR = os.path.join(ROOT_DIR, 'Simulations')
MODEL_RUNS_DIR = os.path.join(SIMULATION_DIR, 'Model_runs')
REFERENCES_SIMULATION_DIR = os.path.join(SIMULATION_DIR, 'Reference_simulations')
OPTIM_DIR = os.path.join(SIMULATION_DIR, 'Optimizations')
LOG_FILE = os.path.join(SIMULATION_DIR, 'Simulations_log.csv')

# Optimization log. It is kept OUTSIDE the repository when the companion log directory
# `_private_LOG/` is there, so that the same log can be shared across machines without
# going through this repository; otherwise it sits next to the simulation log, so that the
# public repository runs on its own. Resolved once, at import.
PRIVATE_LOG_PATH = os.path.join(ROOT_DIR, '_private_LOG')
OPT_LOG_FILE = (os.path.join(PRIVATE_LOG_PATH, 'Optimizations_log.csv')
                if os.path.isdir(PRIVATE_LOG_PATH)
                else os.path.join(SIMULATION_DIR, 'Optimizations_log.csv'))
