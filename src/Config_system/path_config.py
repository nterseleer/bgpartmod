#TODO: Centralize file path configuration

import os


ROOT_DIR = os.getcwd() #bgpartmod/
DATA_DIR = os.path.join(ROOT_DIR, 'data')
OBSERVATION_DIR = os.path.join(ROOT_DIR, 'Observations')
SIMULATION_DIR = os.path.join(ROOT_DIR, 'Simulations')
MODEL_RUNS_DIR = os.path.join(SIMULATION_DIR, 'Model_runs')
REFERENCES_SIMULATION_DIR = os.path.join(SIMULATION_DIR, 'References_simulations')
LOG_FILE = os.path.join(SIMULATION_DIR, 'Simulations_log.csv')
