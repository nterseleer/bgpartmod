"""
Manages model optimization including running, saving, and analyzing results.
"""

import os
import time
import pickle
import json
import dill
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any

import pandas as pd
import numpy as np

from . import functions as fns
from . import desolver
from src.core import model
from src.config_system import path_config as path_cfg

# Constants
BASE_DIR = '../Simulations'
OPTIMIZATIONS_DIR = path_cfg.OPTIM_DIR
LOG_FILE = path_cfg.OPT_LOG_FILE



class Optimization:
    """Manages model optimization including running, saving, and analyzing results."""

    def __init__(self, verbose=True):
        """Initialize optimization container."""
        self.name = None
        self.optdir = None
        self.files = None  # Will store all file paths
        self.config = None
        self.setup = None
        self.obs = None
        self.calibrated_vars = None
        self.verbose = verbose

    @classmethod
    def run_new(cls,
                dconf: Dict,
                modkwargs: Dict,
                obs: Any,
                optimized_parameters: List[str],
                bounds: Tuple[List[float], List[float]],
                calibrated_vars: Optional[List[str]] = None,
                name: Optional[str] = None,
                population_size: int = 90,
                num_cpus: int = 30,
                num_generations: int = 100,
                badlnl: float = -100000.,
                **solver_kwargs) -> 'Optimization':
        """
        Run new optimization.

        Args:
            dconf: Model configuration
            modkwargs: Model keyword arguments
            obs: Observation data
            optimized_parameters: Parameters to optimize
            bounds: (min_bounds, max_bounds) for parameters
            calibrated_vars: Variables to use in likelihood calculation
            name: Optional name (generated if None)
            population_size: DE population size
            num_generations: Number of generations
            badlnl: Score for failed evaluations
            **solver_kwargs: Additional solver parameters
        """
        instance = cls()
        instance.name = name or instance._get_next_name()
        instance.obs = obs
        instance._setup_directories()

        # Store calibrated variables
        instance.calibrated_vars = calibrated_vars or [
            'Phy_Chl', 'NH4_concentration', 'NO3_concentration',
            'DIP_concentration', 'DSi_concentration', 'Phy_C', 'TEPC_C'
        ]

        # Store configuration
        instance.config = {
            'dconf': dconf,
            'modkwargs': modkwargs,
            'optimized_parameters': optimized_parameters,
            'bounds': bounds,
            'population_size': population_size,
            'num_cpus': num_cpus,
            'num_generations': num_generations,
            'badlnl': badlnl,
            'solver_kwargs': solver_kwargs,
            'calibrated_vars': instance.calibrated_vars,
            'observation_info': {
                'station': obs.station,
                'name': obs.name
            }
        }

        instance.obs = obs  # Store observation data

        # Save initial state
        instance._save_configuration()
        instance._save_observations()
        instance._save_calibration_info()

        # Run optimization
        instance._run_optimization(obs)

        return instance

    @classmethod
    def load_existing(cls, name: str) -> 'Optimization':
        """Load existing optimization results."""
        instance = cls()
        instance.name = name

        # Set up directories and file paths
        instance._setup_directories()

        if not os.path.exists(instance.optdir):
            raise FileNotFoundError(f"Optimization {name} not found in {OPTIMIZATIONS_DIR}")

        # Load configuration
        if not os.path.exists(instance.files['config']):
            raise FileNotFoundError(f"Configuration file not found: {instance.files['config']}")

        with open(instance.files['config'], 'rb') as f:
            instance.config = pickle.load(f)
            instance.calibrated_vars = instance.config['calibrated_vars']

        # Load Observations
        obs_path = os.path.join(instance.optdir,
                                f"{instance.name}_observations",
                                f"{instance.name}_observations.pkl")
        if os.path.exists(obs_path):
            with open(obs_path, 'rb') as f:
                instance.obs = pickle.load(f)
        else:
            print(f"Warning: Could not load original observations from {obs_path}")

        return instance

    def _get_next_name(self) -> str:
        """Generate next available optimization name."""
        if not os.path.exists(OPTIMIZATIONS_DIR):
            os.makedirs(OPTIMIZATIONS_DIR)
            return 'OPT000'

        existing = [d for d in os.listdir(OPTIMIZATIONS_DIR)
                    if os.path.isdir(os.path.join(OPTIMIZATIONS_DIR, d))]
        numbers = [int(d[3:]) for d in existing if d.startswith('OPT')]
        next_num = max(numbers + [-1]) + 1
        return f'OPT{next_num:03d}'

    def _setup_directories(self):
        """Create optimization directory structure and define file paths."""
        self.optdir = os.path.join(OPTIMIZATIONS_DIR, self.name)
        os.makedirs(self.optdir, exist_ok=True)

        # Define all file paths with correct extensions
        self.files = {
            'results': os.path.join(self.optdir, f"{self.name}_optimization_results.txt"),  # Changed from .csv
            'config': os.path.join(self.optdir, f"{self.name}_config.pkl"),
            'best_model': os.path.join(self.optdir, f"{self.name}_best_model_dill.pkl"),
            # Added _dill to match original
            'winning_config': os.path.join(self.optdir, f"{self.name}_WINNING_CONFIG.pkl"),
            'summary': os.path.join(self.optdir, f"{self.name}_SUMMARY.pkl")
        }

    def _save_configuration(self):
        """Save all configuration information for reproducibility."""
        # Save main configuration
        with open(os.path.join(self.optdir, f"{self.name}_config.pkl"), 'wb') as f:
            pickle.dump(self.config, f)

        # Save human-readable versions
        fns.write_dict_to_file(self.config['dconf'],
                               f"{self.name}_CONFIG",
                               fdir=self.optdir)
        fns.write_dict_to_file(self.config['modkwargs'],
                               f"{self.name}_modkwargs",
                               fdir=self.optdir,
                               serialize_objects=False)

        # Write readme
        readme = (f"Optimization: {self.name}\n"
                  f"Created: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                  f"Parameters: {', '.join(self.config['optimized_parameters'])}\n"
                  f"Generations: {self.config['num_generations']}\n"
                  f"Population: {self.config['population_size']}")

        with open(os.path.join(self.optdir, f"{self.name}_README.txt"), 'w') as f:
            f.write(readme)

    def _load_configuration(self):
        """Load optimization configuration."""
        config_file = os.path.join(self.optdir, f"{self.name}_config.pkl")
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Configuration file not found: {config_file}")

        with open(config_file, 'rb') as f:
            self.config = pickle.load(f)

    def _save_observations(self):
        """Save observation data and metadata."""
        # Create Observations directory with optimization name prefix
        obs_dir = os.path.join(self.optdir, f"{self.name}_observations")
        os.makedirs(obs_dir, exist_ok=True)

        # Save full observation object with prefix
        obs_pkl_path = os.path.join(obs_dir, f"{self.name}_observations.pkl")
        with open(obs_pkl_path, 'wb') as f:
            pickle.dump(self.obs, f)

        # Save metadata summary with prefix
        summary = self.obs.create_summary()
        summary_path = os.path.join(obs_dir, f"{self.name}_observations_summary.json")
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)

        # Try to create symbolic link to original data with prefix
        try:
            source_path = os.path.abspath(os.path.join(self.obs.datadir, f"{self.obs.station}.feather"))
            link_path = os.path.join(obs_dir, f"{self.name}_original_data.feather")
            if os.path.exists(source_path) and not os.path.exists(link_path):
                os.symlink(source_path, link_path)
        except Exception as e:
            print(f"Could not create symbolic link to original data: {e}")

    def _save_calibration_info(self):
        """Save calibration variables information."""
        calib_info = {
            "calibrated_variables": self.calibrated_vars,
            "datetime": datetime.now().isoformat()
        }

        # Save in optimization directory with optimization name prefix
        calib_path = os.path.join(self.optdir, f"{self.name}_calibration_info.json")
        with open(calib_path, 'w') as f:
            json.dump(calib_info, f, indent=2)

    def _run_optimization(self, obs):
        """Run the optimization using DESolver."""
        self.obs = obs  # Store observation data for evaluation
        min_bounds, max_bounds = self.config['bounds']

        # Configure solver
        solver = desolver.DESolver(
            job=self,  # The solver will access self.files through this
            populationSize=self.config['population_size'],
            maxGenerations=self.config['num_generations'],
            minInitialValue=min_bounds,
            maxInitialValue=max_bounds,
            num_cpus=self.config['num_cpus'],
            **self.config['solver_kwargs']
        )

        # Run optimization
        start_time = time.time()
        self.solver_results = solver.Solve()
        self.runtime = time.time() - start_time
        # Save final results including runtime
        self._save_results()

    def evaluate_model(self, parameters):
        """Evaluate a parameter set (called by solver)."""
        try:
            # Create new configuration with updated parameters
            newconfig = fns.update_config(
                self.config['dconf'],
                self.config['optimized_parameters'],
                parameters
            )

            # Run the model
            trial = model.Model(newconfig, **self.config['modkwargs'])

            # Calculate likelihood using calibrated_vars
            lnl = trial.get_likelihood(
                self.obs,
                calibrated_vars=self.calibrated_vars,
                verbose=False,
                _cached_obs=self.obs.df
            )

            # Handle invalid results
            if lnl is None or np.isinf(lnl) or np.isnan(lnl):
                return self.config['badlnl']

            return lnl

        except (RuntimeWarning, FloatingPointError, ValueError, ZeroDivisionError) as e:
            # If warnings are treated as errors, they'll be caught here
            if self.verbose:
                print(f"Model evaluation failed with error: {e}")
            return self.config['badlnl']


    def _save_results(self):
        """Save optimization results and winning configuration."""
        # Save solver results
        results_dict = {
            'solver_results': self.solver_results,
            'runtime': self.runtime,
            'timestamp': datetime.now()
        }
        with open(os.path.join(self.optdir, 'solver_results.pkl'), 'wb') as f:
            pickle.dump(results_dict, f)

        # Process and save winning configuration
        if hasattr(self, 'summary'):
            winner = self.summary['best_parameters']
            winning_config = fns.update_config(
                self.config['dconf'],
                self.config['optimized_parameters'],
                winner
            )
            fns.write_dict_to_file(
                winning_config,
                f"{self.name}_WINNING_CONFIG",
                fdir=self.optdir
            )

    def process_results(self, rerun_best: bool = False):
        """Process optimization results to get best parameters and statistics."""
        # Load raw results (now expecting .txt file)
        if not os.path.exists(self.files['results']):
            raise FileNotFoundError(
                f"Results file not found: {self.files['results']}\n"
                f"Note: Expected format is 'name_optimization_results.txt'"
            )

        # Read as txt file with comma separator
        self.df = pd.read_csv(self.files['results'], sep=',')
        self.df['cost_raw'] = self.df['cost']
        self.df.loc[self.df['cost'] == self.config['badlnl'], 'cost'] = np.nan

        # Find best parameters
        best_idx = self.df['cost'].idxmax()
        self.winner = self.df.loc[best_idx]

        # Calculate summary statistics
        self.summary = {
            'best_parameters': {
                param: self.winner[param]
                for param in self.config['optimized_parameters']
            },
            'best_score': float(self.winner['cost']),
            'convergence': self._analyze_convergence(),
            'parameter_ranges': self._analyze_parameter_ranges()
        }

        # Save winning configuration
        winning_config = fns.update_config(
            self.config['dconf'],
            self.config['optimized_parameters'],
            self.summary['best_parameters'].values()
        )
        fns.write_dict_to_file(
            winning_config,
            f"{self.name}_WINNING_CONFIG",
            fdir=self.optdir
        )

        # Save summary stats
        fns.write_dict_to_file(
            self.summary,
            f"{self.name}_SUMMARY",
            fdir=self.optdir
        )

        return self.summary

    def _analyze_convergence(self) -> Dict:
        """Analyze optimization convergence."""
        generations = self.df['generation'].unique()

        return {
            'generations_completed': len(generations),
            'generations_planned': self.config['num_generations'],
            'score_progression': {
                'initial': float(self.df.loc[self.df['generation'] == 1, 'cost'].max()),
                'final': float(self.df.loc[self.df['generation'] == len(generations), 'cost'].max()),
                'improvement': float(self.df['cost'].max() - self.df.loc[self.df['generation'] == 1, 'cost'].max())
            },
            'last_improvement': int(self.df.loc[self.df['cost'] == self.df['cost'].max(), 'generation'].iloc[0])
        }

    def _analyze_parameter_ranges(self) -> Dict:
        """Analyze explored parameter ranges."""
        param_stats = {}
        for param in self.config['optimized_parameters']:
            values = self.df[param]
            param_stats[param] = {
                'min': float(values.min()),
                'max': float(values.max()),
                'mean': float(values.mean()),
                'std': float(values.std()),
                'best': float(self.winner[param]),
                'bounds': {
                    'min': float(self.config['bounds'][0][
                                     self.config['optimized_parameters'].index(param)
                                 ]),
                    'max': float(self.config['bounds'][1][
                                     self.config['optimized_parameters'].index(param)
                                 ])
                }
            }
        return param_stats

    def get_best_model(self, force_rerun: bool = False) -> Any:
        """
        Get best model from optimization.

        Args:
            force_rerun: Whether to rerun model even if saved version exists

        Returns:
            Best model instance
        """
        if not hasattr(self, 'summary'):
            self.process_results()

        model_file = os.path.join(self.optdir, f"{self.name}_best_model.pkl")

        # Return existing unless forced rerun
        if os.path.exists(model_file) and not force_rerun:
            with open(model_file, 'rb') as f:
                return dill.load(f)

        # Create and run best model
        best_config = fns.update_config(
            self.config['dconf'],
            self.config['optimized_parameters'],
            self.summary['best_parameters'].values()
        )

        # Add additional diagnostics for best model
        model_kwargs = {
            **self.config['modkwargs'],
            'verbose': True,
            'do_diagnostics': True,
            'full_diagnostics': False
        }

        best_model = model.Model(best_config, **model_kwargs)

        # Save model
        with open(model_file, 'wb') as f:
            dill.dump(best_model, f)

        return best_model

    def compare_model_vs_optim(self, obs, rerun: bool = False):
        """Compare rerun model likelihood with optimization results."""
        best_model = self.get_best_model(force_rerun=rerun)
        model_lnl = best_model.get_likelihood(obs, verbose=True)

        comparison = {
            'optimization_score': self.summary['best_score'],
            'model_score': model_lnl,
            'difference': model_lnl - self.summary['best_score']
        }

        print(f"Optimization score: {comparison['optimization_score']:.6f}")
        print(f"Model score: {comparison['model_score']:.6f}")
        print(f"Difference: {comparison['difference']:.6f}")

        return comparison