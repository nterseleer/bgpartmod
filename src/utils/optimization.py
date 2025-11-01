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

from src.utils import functions as fns
from src.utils import desolver
from src.core import model
from src.core import phys  # Import nécessaire pour la désérialisation pickle des objets Setup
from src.config_system import path_config as path_cfg
from src.config_model import config
from src.utils import simulation_manager as sim_manager

# Constants
OPTIMIZATIONS_DIR = path_cfg.OPTIM_DIR
LOG_FILE = path_cfg.OPT_LOG_FILE




class Optimization:
    """Manages model optimization including running, saves, and analyzing results."""

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
    def _get_next_id(cls) -> str:
        """Get next available optimization ID from the shared log."""
        return sim_manager.get_next_optimization_id()
    
    @classmethod 
    def _prompt_user_note(cls) -> str:
        """Prompt user for optimization note."""
        print("\n" + "="*60)
        print("OPTIMIZATION SETUP - User Note")
        print("="*60)
        note = input("Please describe this optimization (purpose, approach, expectations):\n> ")
        return note.strip()

    @classmethod
    def run_new(cls,
                dconf: Dict,
                modkwargs: Dict,
                obs: Any,
                optimized_parameters: List[str],
                bounds: Tuple[List[float], List[float]],
                setup: Optional[Any] = None,
                calibrated_vars: Optional[List[str]] = None,
                name: Optional[str] = None,
                user_note: Optional[str] = None,
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
            user_note: Optional user note (prompts if None)
            population_size: DE population size
            num_generations: Number of generations
            badlnl: Score for failed evaluations
            **solver_kwargs: Additional solver parameters
        """
        instance = cls()

        # Setup optimization log (prompt user if note not provided)
        instance.name = name or cls._get_next_id()
        if user_note is None:
            user_note = cls._prompt_user_note()

        # Store calibrated variables first to get the count
        instance.calibrated_vars = calibrated_vars or [
            'Phy_Chl', 'NH4_concentration', 'NO3_concentration',
            'DIP_concentration', 'DSi_concentration', 'Phy_C', 'TEPC_C'
        ]

        sim_manager.add_optimization_to_log(
            instance.name,
            len(optimized_parameters),
            user_note,
            calibrated_vars_count=len(instance.calibrated_vars)
        )

        instance.obs = obs
        instance._setup_directories()

        # Handle setup - extract from modkwargs if not provided explicitly
        if setup is not None:
            instance.setup = setup
            modkwargs = {k: v for k, v in modkwargs.items() if k != 'setup'}
        else:
            instance.setup = modkwargs.pop('setup', None)
            if instance.setup is None:
                raise ValueError("setup must be provided either explicitly or in modkwargs")

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
        instance._save_setup()
        instance._save_observations()

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

        # Load setup
        instance._load_setup()

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

    def _save_overview(self):
        """Save optimization overview file (replaces README and calibration_info)."""
        # Detect environment
        environment = "ecmwf" if "ecmwf" in os.getcwd().lower() or "copernicus" in os.getcwd().lower() else "local"

        # Build parameter bounds dictionary (compact format)
        parameter_bounds = {
            param: {
                "min": self.config['bounds'][0][i],
                "max": self.config['bounds'][1][i]
            }
            for i, param in enumerate(self.config['optimized_parameters'])
        }

        # Extract essential setup info
        setup_info = {
            "tmin": self.setup.tmin,
            "tmax": self.setup.tmax,
            "dt": self.setup.dt,
            "start_date": self.setup.start_date.strftime('%Y-%m-%d'),
            "station": self.obs.station if hasattr(self.obs, 'station') else "unknown"
        }

        overview = {
            "optimization_id": self.name,
            "created": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            "environment": environment,

            "parameters": parameter_bounds,

            "calibrated_variables": self.calibrated_vars,

            "optimization_settings": {
                "population_size": self.config['population_size'],
                "generations_planned": self.config['num_generations'],
                "num_cpus": self.config['num_cpus'],
                "badlnl": self.config['badlnl']
            },

            "setup": setup_info,

            "observations": self.config['observation_info'],

            "user_notes": {
                "pre": "",
                "post": ""
            }
        }

        # Save as JSON
        overview_path = os.path.join(self.optdir, f"{self.name}_OVERVIEW.json")
        with open(overview_path, 'w') as f:
            json.dump(overview, f, indent=2)

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

        # Save overview file
        self._save_overview()

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


    def _save_setup(self):
        """Save setup separately from modkwargs."""
        setup_path = os.path.join(self.optdir, f"{self.name}_setup.pkl")
        with open(setup_path, 'wb') as f:
            pickle.dump(self.setup, f)

        # Save human-readable version using Setup's to_dict method
        try:
            setup_dict = self.setup.to_dict()
            fns.write_dict_to_file(setup_dict, f"{self.name}_setup_readable", fdir=self.optdir)
        except Exception as e:
            print(f"Could not save human-readable setup: {e}")

        if self.verbose:
            print(f'Setup saved to {setup_path}')

    def _load_setup(self):
        """Load setup separately with backward compatibility."""
        setup_path = os.path.join(self.optdir, f"{self.name}_setup.pkl")
        if os.path.exists(setup_path):
            with open(setup_path, 'rb') as f:
                self.setup = pickle.load(f)
        else:
            # Backward compatibility: extract from modkwargs
            self.setup = self.config['modkwargs'].get('setup', None)
            if self.setup is None:
                print(f"Warning: Could not load setup for optimization {self.name}")
            else:
                print(f"Loaded setup from legacy modkwargs for optimization {self.name}")

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

    def evaluate_model(self, parameters):
        """Evaluate a parameter set (called by solver)."""
        try:
            # Create parameter dictionary from names and values
            param_dict = dict(zip(self.config['optimized_parameters'], parameters))

            # Create new configuration with updated parameters
            newconfig = fns.update_config(self.config['dconf'], param_dict)

            # Run the model
            trial = model.Model(newconfig, setup=self.setup, **self.config['modkwargs'])

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

        # Calculate runtime from timestamps in the file (robust to SCP transfers and manual edits)
        if not hasattr(self, 'runtime') and 'timestamp' in self.df.columns:
            try:
                # Parse first and last timestamps
                first_time = pd.to_datetime(self.df['timestamp'].iloc[0], format='%Y/%m/%d %H:%M:%S')
                last_time = pd.to_datetime(self.df['timestamp'].iloc[-1], format='%Y/%m/%d %H:%M:%S')
                # Calculate runtime in seconds
                self.runtime = (last_time - first_time).total_seconds()
            except Exception as e:
                if self.verbose:
                    print(f"Warning: Could not calculate runtime from timestamps: {e}")
                # Fallback: set runtime to None so it won't be logged
                self.runtime = None

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

        # Save winning configuration (only if it doesn't exist)
        winning_config_path = os.path.join(self.optdir, f"{self.name}_WINNING_CONFIG.json")
        if not os.path.exists(winning_config_path):
            # self.summary['best_parameters'] is already a dict
            winning_config = fns.update_config(
                self.config['dconf'],
                self.summary['best_parameters']
            )
            fns.write_dict_to_file(
                winning_config,
                f"{self.name}_WINNING_CONFIG",
                fdir=self.optdir
            )

        # Save summary stats (only if it doesn't exist)
        summary_path = os.path.join(self.optdir, f"{self.name}_SUMMARY.json")
        if not os.path.exists(summary_path):
            fns.write_dict_to_file(
                self.summary,
                f"{self.name}_SUMMARY",
                fdir=self.optdir
            )
        
        # Check boundary hits
        boundary_hits = self._check_boundary_hits(threshold_percent=5.0)

        # Calculate runtime_minutes (None if not available - will appear as empty cell in CSV)
        runtime_minutes = (self.runtime / 60.0) if (hasattr(self, 'runtime') and self.runtime is not None) else None

        # Update optimization log with full convergence information
        sim_manager.update_optimization_log_convergence(
            opt_id=self.name,
            likelihood_score=self.summary['best_score'],
            runtime_minutes=runtime_minutes,
            generations_completed=self.summary['convergence']['generations_completed'],
            generations_planned=self.summary['convergence']['generations_planned'],
            first_valid_score=self.summary['convergence']['first_valid_score'],
            first_valid_gen=self.summary['convergence']['first_valid_gen'],
            score_improvement=self.summary['convergence']['score_progression']['improvement'],
            last_improvement_gen=self.summary['convergence']['last_improvement'],
            params_at_lower_bound=boundary_hits['params_at_lower_bound'],
            params_at_upper_bound=boundary_hits['params_at_upper_bound']
        )

        return self.summary

    def _analyze_convergence(self) -> Dict:
        """Analyze optimization convergence with first valid score handling."""
        generations = self.df['generation'].unique()

        # Find first valid (non-badlnl) score using cost_raw (before NaN replacement)
        valid_scores = self.df[self.df['cost_raw'] != self.config['badlnl']]
        if len(valid_scores) > 0:
            first_valid_score = float(valid_scores.iloc[0]['cost'])
            first_valid_gen = int(valid_scores.iloc[0]['generation'])
            score_improvement = float(self.df['cost'].max() - first_valid_score)
        else:
            # All scores are badlnl - should not happen in successful optimization
            first_valid_score = self.config['badlnl']
            first_valid_gen = None
            score_improvement = 0.0

        return {
            'generations_completed': len(generations),
            'generations_planned': self.config['num_generations'],
            'score_progression': {
                'initial': float(self.df.loc[self.df['generation'] == 1, 'cost'].max()),
                'final': float(self.df.loc[self.df['generation'] == len(generations), 'cost'].max()),
                'improvement': score_improvement
            },
            'first_valid_score': first_valid_score,
            'first_valid_gen': first_valid_gen,
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

    def _check_boundary_hits(self, threshold_percent: float = 5.0) -> Dict[str, str]:
        """
        Check which parameters hit their boundaries.

        Args:
            threshold_percent: Percentage of range to consider as "at boundary" (default: 5%)

        Returns:
            Dict with 'params_at_lower_bound' and 'params_at_upper_bound' as comma-separated strings
        """
        lower_hits = []
        upper_hits = []

        for param in self.config['optimized_parameters']:
            idx = self.config['optimized_parameters'].index(param)
            min_bound = self.config['bounds'][0][idx]
            max_bound = self.config['bounds'][1][idx]
            best_val = self.summary['best_parameters'][param]

            range_width = max_bound - min_bound
            threshold = threshold_percent / 100.0 * range_width

            if best_val < min_bound + threshold:
                lower_hits.append(param)
            if best_val > max_bound - threshold:
                upper_hits.append(param)

        return {
            'params_at_lower_bound': ', '.join(lower_hits) if lower_hits else '',
            'params_at_upper_bound': ', '.join(upper_hits) if upper_hits else ''
        }

    def get_best_model(self, force_rerun: bool = False, savemodel: bool = False) -> Any:
        """
        Get best model from optimization.

        Args:
            force_rerun: Whether to rerun model even if saved version exists
            savemodel: Whether to save (dill) the model to pkl object

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
            self.summary['best_parameters']
        )

        # add all diagnostics for subsequent plots and post-processing
        best_config = fns.deep_update(best_config, config.plotting_diagnostics)

        # Add additional diagnostics for best model
        model_kwargs = {
            **self.config['modkwargs'],
            'verbose': True,
            'do_diagnostics': True,
            'full_diagnostics': False
        }

        best_model = model.Model(best_config, setup=self.setup, name=self.name, **model_kwargs)

        # Save model
        if savemodel:
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