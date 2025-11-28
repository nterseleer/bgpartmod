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

PRIVATE_PATH = os.path.join(ROOT_DIR, '_private')
PRIVATE_LOG_PATH = os.path.join(ROOT_DIR, '_private_LOG')
PRIVATE_OPT_LOG_FILE = os.path.join(PRIVATE_LOG_PATH, 'Optimizations_log.csv')
PRIVATE_FIGS = os.path.join(PRIVATE_PATH, 'figs')
