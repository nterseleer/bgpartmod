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
OPTIM_DIR = os.path.join(SIMULATION_DIR, 'Optimizations')
LOG_FILE = os.path.join(SIMULATION_DIR, 'Simulations_log.csv')
OPT_LOG_FILE = os.path.join(SIMULATION_DIR, 'Optimizations_log.csv')
STRAINS_INIT_COND_FILE = os.path.join(OBSERVATION_DIR, 'init_condition.CSV')

# =============================================================================
# OPTIMIZATION FOLDER STRUCTURE CONSTANTS
# =============================================================================

# Standard subdirectories within each optimization folder
OPTIM_SUBDIRS = {
    'metadata': 'metadata',
    'results': 'results'
}

# Standard filenames for optimization files
OPTIM_FILENAMES = {
    'descriptions': 'descriptions.json',
    'readme': 'README.txt',
    'calibration_info': 'calibration_info.json',
    'config': 'config.pkl',
    'optimization_results': 'optimization_results.txt',
    'best_model': 'best_model.pkl',
    'summary': 'SUMMARY.pkl',
    'solver_results': 'solver_results.pkl'
}

# Description-specific subdirectories
DESC_SUBDIRS = {
    'config': 'config',
    'modkwargs': 'modkwargs',
    'observations': 'observations'
}

# Description-specific filenames
DESC_FILENAMES = {
    'base_config_json': 'base_config.json',
    'base_config_pkl': 'base_config.pkl',
    'winning_config_json': 'winning_config.json',
    'winning_config_pkl': 'winning_config.pkl',
    'modkwargs_json': 'modkwargs.json',
    'modkwargs_pkl': 'modkwargs.pkl'
}

PRIVATE_PATH = os.path.join(ROOT_DIR, '_private_LOG')
PRIVATE_OPT_LOG_FILE = os.path.join(PRIVATE_PATH, 'Optimizations_log.csv')
