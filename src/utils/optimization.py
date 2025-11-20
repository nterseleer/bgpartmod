"""
Manages model optimization including running, saving, and analyzing results.
"""

import os
import time
import pickle
import json
import dill
from dataclasses import dataclass
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


@dataclass
class ModelCase:
    """
    Configuration for a single evaluation case in optimization.

    A case represents one model configuration to evaluate: different station,
    species, time period, formulation, or any combination.

    Attributes:
        case_id: Unique identifier
        dconf: Model configuration dictionary
        setup: Physical setup instance
        obs: Observation data instance
        calibrated_vars: Variables for likelihood calculation
        weight: Weight in aggregated score
        case_label: Optional human-readable description
    """
    case_id: str
    dconf: Dict[str, Any]
    setup: Any
    obs: Any
    calibrated_vars: List[str]
    weight: float = 1.0
    case_label: Optional[str] = None

    def __post_init__(self):
        if self.case_label is None:
            self.case_label = self.case_id


class Optimization:
    """Manages model optimization including running, saves, and analyzing results."""

    def __init__(self, verbose=True):
        """Initialize optimization container."""
        self.name = None
        self.optdir = None
        self.files = None
        self.config = None
        self.setup = None
        self.obs = None
        self.calibrated_vars = None
        self.verbose = verbose
        self.cases = None
        self.is_multi_case = False
        self.param_metadata = None

    @staticmethod
    def _parse_parameter_name(param_name: str) -> Tuple[str, Optional[str]]:
        """
        Parse parameter name to extract base name and case_id.

        Returns:
            (base_param_name, case_id or None)
        """
        if '@' in param_name:
            base, case_id = param_name.rsplit('@', 1)
            return base, case_id
        return param_name, None

    @classmethod
    def _filter_parameters_for_case(cls, param_dict: Dict[str, float],
                                     case_id: str) -> Dict[str, float]:
        """
        Filter and clean parameters for a specific case.

        For multi-case optimizations, parameters can be:
        - Global (no @): apply to all cases
        - Case-specific (@case_id): apply only to that case

        This method:
        1. Keeps all global parameters (no @)
        2. Keeps only case-specific parameters matching the given case_id
        3. Removes the @case_id suffix from parameter names

        Args:
            param_dict: Dictionary of parameters, potentially with @case_id suffixes
            case_id: The case ID to filter for

        Returns:
            Filtered dictionary with cleaned parameter names

        Example:
            >>> params = {
            ...     'Phy+mu_max': 1.5,                          # global
            ...     'Phy+divide_water_depth_ratio@MOW1': 0.8,  # case-specific
            ...     'Phy+divide_water_depth_ratio@W05': 0.6    # other case
            ... }
            >>> _filter_parameters_for_case(params, 'MOW1')
            {'Phy+mu_max': 1.5, 'Phy+divide_water_depth_ratio': 0.8}
        """
        filtered_params = {}

        for param_name, value in param_dict.items():
            base_name, param_case_id = cls._parse_parameter_name(param_name)

            if param_case_id is None:
                # Global parameter - include for all cases
                filtered_params[param_name] = value
            elif param_case_id == case_id:
                # Case-specific parameter matching this case - include with cleaned name
                filtered_params[base_name] = value
            # else: parameter for another case - skip

        return filtered_params

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
    def _validate_case_specific_params(cls, optimized_parameters: List[str],
                                       cases: List[ModelCase]):
        """Validate that all @case_id references exist."""
        case_ids = {c.case_id for c in cases}
        for param_name in optimized_parameters:
            _, case_id = cls._parse_parameter_name(param_name)
            if case_id is not None and case_id not in case_ids:
                raise ValueError(
                    f"Parameter '{param_name}' references case_id '{case_id}', "
                    f"but no case with that ID exists. Available: {case_ids}"
                )

    def _categorize_parameters(self, optimized_parameters: List[str],
                               bounds: Tuple[List[float], List[float]]):
        """Categorize parameters into shared vs case-specific with metadata."""
        self.param_metadata = []
        for i, param_name in enumerate(optimized_parameters):
            base_name, case_id = self._parse_parameter_name(param_name)
            self.param_metadata.append({
                'full_name': param_name,
                'base_name': base_name,
                'case_id': case_id,
                'index': i,
                'bounds': (bounds[0][i], bounds[1][i])
            })

    def _build_case_param_dict(self, parameters: np.ndarray, case_id: str) -> Dict[str, float]:
        """
        Build parameter dictionary for a specific case.

        Includes shared parameters and case-specific parameters matching case_id.
        """
        param_dict = {}
        for meta in self.param_metadata:
            if meta['case_id'] is None:
                # Shared parameter
                param_dict[meta['base_name']] = parameters[meta['index']]
            elif meta['case_id'] == case_id:
                # Case-specific parameter for this case
                param_dict[meta['base_name']] = parameters[meta['index']]
        return param_dict

    def _group_parameters_by_case(self) -> Dict[str, Dict[str, float]]:
        """Group best parameters by case for readability."""
        grouped = {'shared': {}}
        for meta in self.param_metadata:
            value = self.summary['best_parameters'][meta['full_name']]
            if meta['case_id'] is None:
                grouped['shared'][meta['base_name']] = value
            else:
                if meta['case_id'] not in grouped:
                    grouped[meta['case_id']] = {}
                grouped[meta['case_id']][meta['base_name']] = value
        return grouped

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
            calibrated_vars_count=len(instance.calibrated_vars),
            case_mode="Single"
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
    def run_new_multi_case(cls,
                           cases: List[ModelCase],
                           modkwargs: Dict,
                           optimized_parameters: List[str],
                           bounds: Tuple[List[float], List[float]],
                           name: Optional[str] = None,
                           user_note: Optional[str] = None,
                           population_size: int = 90,
                           num_cpus: int = 30,
                           num_generations: int = 100,
                           badlnl: float = -100000.,
                           **solver_kwargs) -> 'Optimization':
        """
        Run multi-case optimization.

        Args:
            cases: List of ModelCase instances
            modkwargs: Model kwargs (shared across cases)
            optimized_parameters: Parameters to optimize (shared)
            bounds: Parameter bounds (shared)
            name: Optional name
            user_note: Optional note
            population_size: DE population size
            num_generations: Number of generations
            badlnl: Score for failed evaluations
            **solver_kwargs: Additional solver parameters
        """
        instance = cls()
        instance.is_multi_case = True
        instance.cases = cases

        # Validate and categorize parameters
        cls._validate_case_specific_params(optimized_parameters, cases)
        instance._categorize_parameters(optimized_parameters, bounds)

        # Setup optimization log
        instance.name = name or cls._get_next_id()
        if user_note is None:
            user_note = cls._prompt_user_note()

        # Count total calibrated vars (sum across cases)
        total_calibrated_vars_count = sum(len(c.calibrated_vars) for c in cases)

        sim_manager.add_optimization_to_log(
            instance.name,
            len(optimized_parameters),
            f"**MULTI X{len(cases)} ** | {user_note}",
            calibrated_vars_count=total_calibrated_vars_count,
            case_mode="Multi"
        )

        instance._setup_directories()

        # Store configuration
        instance.config = {
            'modkwargs': modkwargs,
            'optimized_parameters': optimized_parameters,
            'bounds': bounds,
            'population_size': population_size,
            'num_cpus': num_cpus,
            'num_generations': num_generations,
            'badlnl': badlnl,
            'solver_kwargs': solver_kwargs,
            'cases_info': [
                {
                    'case_id': c.case_id,
                    'case_label': c.case_label,
                    'calibrated_vars': c.calibrated_vars,
                    'observation_info': {'station': c.obs.station, 'name': c.obs.name},
                    'weight': c.weight
                }
                for c in cases
            ]
        }

        # Save initial state
        instance._save_configuration_multi_case()
        instance._save_cases()

        # Run optimization
        instance._run_optimization_multi_case()

        return instance

    @classmethod
    def load_existing(cls, name: str) -> 'Optimization':
        """Load existing optimization (single-case or multi-case)."""
        instance = cls()
        instance.name = name
        instance._setup_directories()

        if not os.path.exists(instance.optdir):
            raise FileNotFoundError(f"Optimization {name} not found in {OPTIMIZATIONS_DIR}")

        if not os.path.exists(instance.files['config']):
            raise FileNotFoundError(f"Configuration file not found: {instance.files['config']}")

        with open(instance.files['config'], 'rb') as f:
            instance.config = pickle.load(f)

        # Detect mode
        if 'cases_info' in instance.config:
            # Multi-case mode
            instance.is_multi_case = True
            instance._load_cases()
            # Rebuild param_metadata for loaded optimization
            instance._categorize_parameters(
                instance.config['optimized_parameters'],
                instance.config['bounds']
            )
        else:
            # Single-case mode
            instance.is_multi_case = False
            instance.calibrated_vars = instance.config['calibrated_vars']
            instance._load_setup()
            obs_path = os.path.join(instance.optdir, f"{instance.name}_observations",
                                    f"{instance.name}_observations.pkl")
            if os.path.exists(obs_path):
                with open(obs_path, 'rb') as f:
                    instance.obs = pickle.load(f)
            else:
                print(f"Warning: Could not load observations from {obs_path}")

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
        """Save single-case optimization overview."""
        environment = "ecmwf" if "ecmwf" in os.getcwd().lower() or "copernicus" in os.getcwd().lower() else "local"

        parameter_bounds = {
            param: {"min": self.config['bounds'][0][i], "max": self.config['bounds'][1][i]}
            for i, param in enumerate(self.config['optimized_parameters'])
        }

        setup_info = {
            "tmin": self.setup.tmin,
            "tmax": self.setup.tmax,
            "dt": self.setup.dt,
            "start_date": self.setup.start_date.strftime('%Y-%m-%d'),
            "station": self.obs.station if hasattr(self.obs, 'station') else "unknown"
        }

        overview = {
            "optimization_id": self.name,
            "mode": "single-case",
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
            "user_notes": {"pre": "", "post": ""}
        }

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

    def _save_cases(self):
        """Save all case configurations."""
        cases_dir = os.path.join(self.optdir, f"{self.name}_cases")
        os.makedirs(cases_dir, exist_ok=True)

        for case in self.cases:
            case_subdir = os.path.join(cases_dir, case.case_id)
            os.makedirs(case_subdir, exist_ok=True)

            with open(os.path.join(case_subdir, "dconf.pkl"), 'wb') as f:
                pickle.dump(case.dconf, f)
            fns.write_dict_to_file(case.dconf, f"{case.case_id}_dconf", fdir=case_subdir)

            with open(os.path.join(case_subdir, "setup.pkl"), 'wb') as f:
                pickle.dump(case.setup, f)

            with open(os.path.join(case_subdir, "observations.pkl"), 'wb') as f:
                pickle.dump(case.obs, f)

            with open(os.path.join(case_subdir, "calibrated_vars.json"), 'w') as f:
                json.dump(case.calibrated_vars, f, indent=2)

            metadata = {
                'case_id': case.case_id,
                'case_label': case.case_label,
                'weight': case.weight,
                'obs_station': case.obs.station,
                'obs_name': case.obs.name
            }
            with open(os.path.join(case_subdir, "metadata.json"), 'w') as f:
                json.dump(metadata, f, indent=2)

    def _load_cases(self):
        """Load all case configurations."""
        cases_dir = os.path.join(self.optdir, f"{self.name}_cases")
        self.cases = []

        for case_info in self.config['cases_info']:
            case_id = case_info['case_id']
            case_subdir = os.path.join(cases_dir, case_id)

            with open(os.path.join(case_subdir, "dconf.pkl"), 'rb') as f:
                dconf = pickle.load(f)

            with open(os.path.join(case_subdir, "setup.pkl"), 'rb') as f:
                setup = pickle.load(f)

            with open(os.path.join(case_subdir, "observations.pkl"), 'rb') as f:
                obs = pickle.load(f)

            case = ModelCase(
                case_id=case_id,
                dconf=dconf,
                setup=setup,
                obs=obs,
                calibrated_vars=case_info['calibrated_vars'],
                weight=case_info.get('weight', 1.0),
                case_label=case_info.get('case_label', case_id)
            )
            self.cases.append(case)

    def _save_configuration_multi_case(self):
        """Save multi-case configuration."""
        with open(os.path.join(self.optdir, f"{self.name}_config.pkl"), 'wb') as f:
            pickle.dump(self.config, f)

        fns.write_dict_to_file(self.config['modkwargs'], f"{self.name}_modkwargs",
                               fdir=self.optdir, serialize_objects=False)

        self._save_overview_multi_case()

    def _save_overview_multi_case(self):
        """Save multi-case optimization overview."""
        environment = "ecmwf" if "ecmwf" in os.getcwd().lower() or "copernicus" in os.getcwd().lower() else "local"

        parameter_bounds = {
            param: {"min": self.config['bounds'][0][i], "max": self.config['bounds'][1][i]}
            for i, param in enumerate(self.config['optimized_parameters'])
        }

        overview = {
            "optimization_id": self.name,
            "mode": "multi-case",
            "created": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            "environment": environment,
            "cases": [
                {
                    "case_id": c.case_id,
                    "case_label": c.case_label,
                    "weight": c.weight,
                    "obs_file": c.obs.station,
                    "calibrated_vars": c.calibrated_vars
                }
                for c in self.cases
            ],
            "parameters": parameter_bounds,
            "optimization_settings": {
                "population_size": self.config['population_size'],
                "generations_planned": self.config['num_generations'],
                "num_cpus": self.config['num_cpus'],
                "badlnl": self.config['badlnl']
            }
        }

        with open(os.path.join(self.optdir, f"{self.name}_OVERVIEW.json"), 'w') as f:
            json.dump(overview, f, indent=2)

    def _run_optimization_multi_case(self):
        """Run multi-case optimization using DESolver."""
        min_bounds, max_bounds = self.config['bounds']

        solver = desolver.DESolver(
            job=self,
            populationSize=self.config['population_size'],
            maxGenerations=self.config['num_generations'],
            minInitialValue=min_bounds,
            maxInitialValue=max_bounds,
            num_cpus=self.config['num_cpus'],
            **self.config['solver_kwargs']
        )

        start_time = time.time()
        self.solver_results = solver.Solve()
        self.runtime = time.time() - start_time

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
        """Evaluate parameter set (single-case or multi-case mode)."""
        try:
            if not self.is_multi_case:
                # Single-case mode
                param_dict = dict(zip(self.config['optimized_parameters'], parameters))
                newconfig = fns.update_config(self.config['dconf'], param_dict)
                trial = model.Model(newconfig, setup=self.setup, **self.config['modkwargs'])
                lnl = trial.get_likelihood(
                    self.obs,
                    calibrated_vars=self.calibrated_vars,
                    verbose=False,
                    _cached_obs=self.obs.df
                )
                if lnl is None or np.isinf(lnl) or np.isnan(lnl):
                    return self.config['badlnl']
                return lnl
            else:
                # Multi-case mode: aggregate scores across all cases
                total_lnl = 0.0
                for case in self.cases:
                    # Build case-specific parameter dict (handles @case_id suffixes)
                    case_param_dict = self._build_case_param_dict(parameters, case.case_id)
                    case_dconf = fns.update_config(case.dconf, case_param_dict)
                    trial = model.Model(case_dconf, setup=case.setup, **self.config['modkwargs'])
                    lnl = trial.get_likelihood(
                        case.obs,
                        calibrated_vars=case.calibrated_vars,
                        verbose=False,
                        _cached_obs=case.obs.df
                    )
                    if lnl is None or np.isinf(lnl) or np.isnan(lnl):
                        return self.config['badlnl']
                    total_lnl += case.weight * lnl
                return total_lnl

        except (RuntimeWarning, FloatingPointError, ValueError, ZeroDivisionError) as e:
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

        # Add grouped parameters for multi-case
        if self.is_multi_case and self.param_metadata is not None:
            self.summary['best_parameters_grouped'] = self._group_parameters_by_case()

        # Save winning configuration (only if it doesn't exist)
        winning_config_path = os.path.join(self.optdir, f"{self.name}_WINNING_CONFIG.json")
        if not os.path.exists(winning_config_path):
            if not self.is_multi_case:
                # Single-case: save one winning config
                winning_config = fns.update_config(
                    self.config['dconf'],
                    self.summary['best_parameters']
                )
                fns.write_dict_to_file(winning_config, f"{self.name}_WINNING_CONFIG", fdir=self.optdir)
            else:
                # Multi-case: save winning config per case
                for case in self.cases:
                    # Filter parameters for this specific case
                    filtered_params = self._filter_parameters_for_case(
                        self.summary['best_parameters'], case.case_id
                    )
                    winning_config = fns.update_config(case.dconf, filtered_params)
                    fns.write_dict_to_file(
                        winning_config,
                        f"{self.name}_WINNING_CONFIG_{case.case_id}",
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

    def get_best_model(self, case_id: Optional[str] = None,
                       force_rerun: bool = False, savemodel: bool = False) -> Any:
        """
        Get best model from optimization.

        Args:
            case_id: For multi-case, which case to reconstruct (None = all)
            force_rerun: Whether to rerun model
            savemodel: Whether to save model

        Returns:
            Best model instance (or dict of models if multi-case and case_id=None)
        """
        if not hasattr(self, 'summary'):
            self.process_results()

        if not self.is_multi_case:
            # Single-case mode
            model_file = os.path.join(self.optdir, f"{self.name}_best_model.pkl")
            if os.path.exists(model_file) and not force_rerun:
                with open(model_file, 'rb') as f:
                    return dill.load(f)

            best_config = fns.update_config(self.config['dconf'], self.summary['best_parameters'])
            best_config = fns.deep_update(best_config, config.plotting_diagnostics)
            model_kwargs = {
                **self.config['modkwargs'],
                'verbose': True,
                'do_diagnostics': True,
                'full_diagnostics': False
            }
            best_model = model.Model(best_config, setup=self.setup, name=self.name, **model_kwargs)

            if savemodel:
                with open(model_file, 'wb') as f:
                    dill.dump(best_model, f)
            return best_model

        else:
            # Multi-case mode
            if case_id is not None:
                # Reconstruct for one specific case
                case = next(c for c in self.cases if c.case_id == case_id)
                model_file = os.path.join(self.optdir, f"{self.name}_best_model_{case_id}.pkl")

                if os.path.exists(model_file) and not force_rerun:
                    with open(model_file, 'rb') as f:
                        return dill.load(f)

                # Filter parameters for this specific case (removes @case_id suffixes)
                filtered_params = self._filter_parameters_for_case(
                    self.summary['best_parameters'], case_id
                )
                best_config = fns.update_config(case.dconf, filtered_params)
                best_config = fns.deep_update(best_config, config.plotting_diagnostics)
                model_kwargs = {
                    **self.config['modkwargs'],
                    'verbose': True,
                    'do_diagnostics': True,
                    'full_diagnostics': False
                }
                best_model = model.Model(best_config, setup=case.setup,
                                         name=f"{self.name}_{case_id}", **model_kwargs)

                if savemodel:
                    with open(model_file, 'wb') as f:
                        dill.dump(best_model, f)
                return best_model

            else:
                # Reconstruct for all cases
                return {c.case_id: self.get_best_model(case_id=c.case_id,
                                                       force_rerun=force_rerun,
                                                       savemodel=savemodel)
                        for c in self.cases}

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