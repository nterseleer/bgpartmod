"""
Manages model optimization including running, saving, and analyzing results.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from datetime import datetime

import numpy as np
import pandas as pd
from typing import Optional, Any

from src.core import model
from src.optimizations import desolver
from src.optimizations.description import Description
from src.optimizations.optimization_metadata import OptimizationMetadata as MetaData
from src.optimizations.optimization_metastructure import OptimizationMetaStructure as MetaStructure
from src.utils import functions as fns
from src.utils import simulation_manager as sim_manager
from src.utils.functions import load_pkl, save_to_json, write_dict_to_file


@dataclass
class OptimizationConfig:
    optimized_parameters: tuple[str, ...]
    bounds: tuple[list[str | float], list[str | float]]
    population_size: int
    num_cpus: int
    num_generations: int
    badlnl: float
    solver_kwargs: dict[str, Any] | None
    calibrated_vars: tuple[str, ...]


class Optimization:
    """Manages model optimization including running, saves, and analyzing results."""

    def __init__(self, verbose=True):
        """Initialize optimization container."""
        self.name = None
        self.creation_date: datetime | None = None
        self.meta_data: MetaData | None = None  # Will store all file paths
        self.descriptions: list[Description] = []
        self.optimization_config: OptimizationConfig | None = None
        self.modkwargs: dict | None = None
        self.is_optimization_run: bool = False
        self.is_optimization_processed: bool = False
        self.verbose = verbose
        self.user_note = ""

    @classmethod
    def _get_next_id(cls) -> str:
        """Get next available optimization ID from the shared log."""
        return sim_manager.get_next_optimization_id()

    @classmethod
    def _prompt_user_note(cls) -> str:
        """Prompt user for optimization note."""
        print("\n" + "=" * 60)
        print("OPTIMIZATION SETUP - User Note")
        print("=" * 60)
        note = input("Please describe this optimization (purpose, approach, expectations):\n> ")
        return note.strip()

    @classmethod
    def build_from_metadata(cls, meta_data: MetaData, is_force_run: bool = False):
        instance = cls()
        instance.meta_data = meta_data
        instance.name = meta_data.optimization_name
        instance.creation_date = instance.meta_data.creation_date
        user_note = meta_data.user_note

        for meta_description in meta_data.descriptions:
            instance.descriptions.append(Description.from_meta_description(meta_description))
        instance.modkwargs = load_pkl(meta_data.modkwargs)

        instance.optimization_config = OptimizationConfig(**meta_data.optimization_config)
        instance.is_optimization_run = meta_data.is_optimization_run
        instance.is_optimization_processed = meta_data.is_optimization_processed

        # Run optimization
        sim_manager.add_optimization_to_log(
            instance.name, len(instance.optimization_config.optimized_parameters), user_note)

        if is_force_run:
            instance.meta_data.is_optimization_run = True
            instance._run_optimization()

        return instance

    @classmethod
    def run_new(cls, description: Description,
                modkwargs: dict,
                optimization_config: OptimizationConfig
                , name: Optional[str] = None
                ) -> 'Optimization':
        """Legacy single-configuration optimization method."""
        return cls.run_multi_config(
            descriptions=[description],
            modkwargs=modkwargs,
            optimization_config=optimization_config,
            name=name
        )

    @classmethod
    def run_multi_config(cls,
                         descriptions: list[Description],
                         modkwargs: dict,
                         optimization_config: OptimizationConfig,
                         name: Optional[str] = None) -> 'Optimization':

        instance = cls()
        instance.creation_date = datetime.now()
        instance.name = name or instance._get_next_id()
        instance.is_optimization_run = False
        instance.is_optimization_processed = False
        instance.modkwargs = modkwargs
        instance.descriptions = descriptions
        instance.user_note = instance._prompt_user_note()

        optimization_config.calibrated_vars = optimization_config.calibrated_vars or ['Phy_Chl', 'NH4_concentration',
                                                                                      'NO3_concentration',
                                                                                      'DIP_concentration',
                                                                                      'DSi_concentration', 'Phy_C',
                                                                                      'TEPC_C']
        instance.optimization_config = optimization_config
        instance.meta_data = MetaData.from_optimization(instance)

        sim_manager.add_optimization_to_log(instance.name, len(optimization_config.optimized_parameters),
                                            instance.user_note)

        # Run optimization
        instance.meta_data.is_optimization_run = True
        instance._save()
        instance._run_optimization()

        return instance

    def _save(self):
        meta_structure = self.meta_data.meta_structure

        self._save_readme(meta_structure)
        self._save_descriptions(meta_structure)
        self._save_meta_data()
        self._save_modkwargs(meta_structure)
        self._save_calibration_info(meta_structure)

    def _save_readme(self, meta_structure: MetaStructure):
        readme = (f"Optimization: {self.name}\n"
                  f"Created: {self.creation_date}\n"
                  f"Descriptions: {len(self.descriptions)}\n"
                  f"Parameters: {', '.join(self.optimization_config.optimized_parameters)}\n"
                  f"Generations: {self.optimization_config.num_generations}\n"
                  f"Population: {self.optimization_config.population_size}")

        with open(os.path.join(meta_structure.optimization_path, f"{self.name}_README.txt"), 'w') as f:
            f.write(readme)

    def _save_descriptions(self, meta_structure: MetaStructure):
        if self.verbose:
            print("saving descriptions")

        for description in self.descriptions:
            description.save(meta_structure)

        if self.verbose:
            print("Descriptions saved")

    def _save_meta_data(self):
        if self.verbose:
            print("saving metadata")

        self.meta_data.save()

        if self.verbose:
            print("Metadata saved")

    def _save_calibration_info(self, meta_structure: MetaStructure):
        if not self.is_optimization_processed:
            """Save calibration variables information."""
            calib_info = {
                "calibrated_variables": self.optimization_config.calibrated_vars,
                "datetime": datetime.now().isoformat()
            }
            # Save in optimization directory with optimization name prefix
            calib_path = meta_structure.calibration_file
            save_to_json(calib_info, calib_path)

    def _save_modkwargs(self, meta_structure):
        if not self.is_optimization_processed:
            if self.verbose:
                print("saving modkwargs")

            write_dict_to_file(self.modkwargs, meta_structure.modkwargs_file)

            if self.verbose:
                print("modkwargs saved")

    def _run_optimization(self):
        """Run the optimization using DESolver."""
        min_bounds, max_bounds = self.optimization_config.bounds
        # Configure solver
        solver_kwargs = self.optimization_config.solver_kwargs or {}
        solver = desolver.DESolver(
            job=self,  # The solver will access self.files through this
            populationSize=self.optimization_config.population_size,
            maxGenerations=self.optimization_config.num_generations,
            minInitialValue=min_bounds,
            maxInitialValue=max_bounds,
            num_cpus=self.optimization_config.num_cpus,
            **solver_kwargs
        )

        solver.Solve()

    def _if_result_can_be_processed(self, force_process: bool) -> bool:
        if self.is_optimization_processed and not force_process:
            warnings.warn("Optimization already processed. To reprocess optimization: set force_process to True")
            return False
        if not self.meta_data.is_optimization_run:
            warnings.warn("You can't process results until optimization is run")
            return False
        return True

    def _generate_best_params(self, metastructure: MetaStructure) -> None:
        if not os.path.exists(metastructure.results_file):
            raise FileNotFoundError(
                f"Results file not found: {metastructure}\n"
                f"Note: Expected format is 'name_optimization_results.txt'"
            )

        # Read as txt file with comma separator
        self.df = pd.read_csv(metastructure.results_file, sep=',')
        self.df['cost_raw'] = self.df['cost']
        self.df.loc[self.df['cost'] == self.optimization_config.badlnl, 'cost'] = np.nan

        # Find best parameters
        best_idx = self.df['cost'].idxmax()

        self.winner = self.df.loc[best_idx]

    def _generate_process_summary(self):
        self.summary = {
            'best_parameters': {param: self.winner[param] for param in self.optimization_config.optimized_parameters},
            'best_score': float(self.winner['cost']),
            'convergence': self._analyze_convergence(),
            'parameter_ranges': self._analyze_parameter_ranges()
        }

    def _generate_winning_configs(self) -> None:
        for description in self.descriptions:
            description.winning_config = fns.update_config(
                description.config,
                self.optimization_config.optimized_parameters,
                self.summary['best_parameters'].values()
            )

    def _generate_best_models(self) -> None:
        if not hasattr(self, 'summary'):
            self.process_results()

        model_kwargs = {
            **self.modkwargs,
            'verbose': True,
            'do_diagnostics': True,
            'full_diagnostics': False
        }

        for description in self.descriptions:
            description.best_model = model.Model(description.winning_config,
                                                 setup=description.setup,
                                                 name=description.name,
                                                 **model_kwargs)

    def process_results(self, force_process: bool = False) -> None:
        if not self._if_result_can_be_processed(force_process): return
        metastructure = self.meta_data.meta_structure

        self._generate_best_params(metastructure)
        self._generate_process_summary()
        self._generate_winning_configs()
        self._generate_best_models()
        self.is_optimization_processed = True
        self.meta_data.update_optimization_process_status(self)
        self._save()

        # Save summary stats (only if it doesn't exist)
        if not os.path.exists(metastructure.summary_file):
            save_to_json(self.summary, metastructure.summary_file)

        # Update optimization log with results
        if hasattr(self, 'runtime'):
            runtime_minutes = self.runtime / 60.0
            sim_manager.update_optimization_log_results(self.name, self.summary['best_score'], runtime_minutes)

    def evaluate_model(self, parameters) -> float:
        total_score = 0.0
        try:
            for description in self.descriptions:
                new_config = fns.update_config(
                    description.config,
                    self.optimization_config.optimized_parameters,
                    parameters
                )

                trial = model.Model(new_config, description.setup, **self.modkwargs)
                lnl = trial.get_likelihood(
                    description.observation,
                    calibrated_vars=self.optimization_config.calibrated_vars,
                    verbose=False,
                    _cached_obs=description.observation.df
                )

                if lnl is None or np.isinf(lnl) or np.isnan(lnl):
                    if self.verbose:
                        print(f"Description failed with invalid likelihood: {lnl}")
                    return self.optimization_config.badlnl

                total_score += lnl
        except (RuntimeWarning, FloatingPointError, ValueError, ZeroDivisionError) as e:
            # If warnings are treated as errors, they'll be caught here
            if self.verbose:
                print(f"Model evaluation failed with error: {e}")
            return self.optimization_config.badlnl

        return total_score

    def _analyze_convergence(self) -> dict:
        """Analyze optimization convergence."""
        generations = self.df['generation'].unique()

        return {
            'generations_completed': len(generations),
            'generations_planned': self.optimization_config.num_generations,
            'score_progression': {
                'initial': float(self.df.loc[self.df['generation'] == 1, 'cost'].max()),
                'final': float(self.df.loc[self.df['generation'] == len(generations), 'cost'].max()),
                'improvement': float(self.df['cost'].max() - self.df.loc[self.df['generation'] == 1, 'cost'].max())
            },
            'last_improvement': int(self.df.loc[self.df['cost'] == self.df['cost'].max(), 'generation'].iloc[0])
        }

    def _analyze_parameter_ranges(self) -> dict:
        """Analyze explored parameter ranges."""
        param_stats = {}
        for param in self.optimization_config.optimized_parameters:
            values = self.df[param]
            param_stats[param] = {
                'min': float(values.min()),
                'max': float(values.max()),
                'mean': float(values.mean()),
                'std': float(values.std()),
                'best': float(self.winner[param]),
                'bounds': {
                    'min': float(self.optimization_config.bounds[0][
                                     self.optimization_config.optimized_parameters.index(param)
                                 ]),
                    'max': float(self.optimization_config.bounds[1][
                                     self.optimization_config.optimized_parameters.index(param)
                                 ])
                }
            }
        return param_stats

    def compare_model_vs_optim(self, obs, rerun: bool = False) -> dict:
        """Compare rerun model likelihood with optimization results."""

        if not self.is_optimization_processed: self._generate_best_models()

        comparaisons = {'optimization_score': self.summary['best_score']}

        print(f"Optimization score: {comparaisons['optimization_score']:.6f}")

        for description in self.descriptions:
            model_lnl = description.best_model.get_likelihood(description.observation, verbose=True)
            comparaisons[description.name] = {'model_score': model_lnl,
                                              'difference': model_lnl - self.summary['best_score']}
            print(f"Model score for description \"{description.name}\":"
                  f" {comparaisons[description.name]['model_score']:.6f}")

        return comparaisons