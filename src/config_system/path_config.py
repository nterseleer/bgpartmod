#TODO: Centralize file path configuration

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
LOG_FILE = os.path.join(SIMULATION_DIR, 'Simulations_log.csv')
OPT_LOG_FILE = os.path.join(SIMULATION_DIR, 'Optimizations_log.csv')
STRAINS_INIT_COND_FILE = os.path.join(OBSERVATION_DIR, 'init_condition.CSV')


PRIVATE_PATH = os.path.join(ROOT_DIR, '_private_LOG')
PRIVATE_OPT_LOG_FILE = os.path.join(PRIVATE_PATH, 'Optimizations_log.csv')


# =============================================================================
# OPTIMIZATION FOLDER STRUCTURE CONSTANTS
# =============================================================================

OPTIMIZATIONS_ROOT_DIR = os.path.join(SIMULATION_DIR, 'Optimizations')

OPTIMIZATIONS_README_SUB_FILE = 'README'
OPTIMIZATIONS_SUMMARY_SUB_FILE = 'summary'
OPTIMIZATIONS_LOG_SUB_FILE = 'log'
OPTIMIZATIONS_CALIBRATION_SUB_FILE = 'calibration'


OPTIMIZATIONS_SRC_SUB_DIR = "src"
OPTIMIZATIONS_OBS_SUB_DIR = os.path.join(OPTIMIZATIONS_SRC_SUB_DIR, "observations")
OPTIMIZATIONS_CONFIG_SUB_DIR = os.path.join(OPTIMIZATIONS_SRC_SUB_DIR, "config")
OPTIMIZATIONS_SETUP_SUB_DIR = os.path.join(OPTIMIZATIONS_SRC_SUB_DIR, "setup")
OPTIMIZATIONS_MOD_KWARGS_SUB_FILE = os.path.join(OPTIMIZATIONS_SRC_SUB_DIR, "modkwargs")

OPTIMIZATIONS_OUT_SUB_DIR = "out"
OPTIMIZATIONS_RESULTS_SUB_FILE = os.path.join(OPTIMIZATIONS_OUT_SUB_DIR,'result')
OPTIMIZATIONS_BEST_MODEL_SUB_DIR = os.path.join(OPTIMIZATIONS_OUT_SUB_DIR, "best_model")
OPTIMIZATIONS_WINNING_CONFIG_SUB_DIR = os.path.join(OPTIMIZATIONS_OUT_SUB_DIR, "winning_config")