import time
from functools import wraps
from collections import defaultdict
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Any
from scipy import integrate

from src.utils import functions as fns
from core import phys
from src.utils import evaluation
from src.components import phytoplankton as phyto
from src.config_model import varinfos


@staticmethod
def _track_time(func):
    """Decorator to track function execution time"""
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        start = time.time()
        result = func(self, *args, **kwargs)
        duration = time.time() - start
        self.perf_stats[func.__name__].append(duration)
        return result
    return wrapper

class Model:
    def __init__(
            self,
            config_dict: Dict[str, Any],
            setup: phys.Setup = phys.Setup(dt=0.001, tmax=20),
            euler: bool = True,
            output_time_offset: float = 0,
            output_config: Optional[Dict] = None,
            output_vars: Optional[List[str]] = None,
            name: str = 'Model',
            verbose: bool = True,
            debug_budgets: bool = False,
            aggregate_vars: Optional[List[str]] = None,
            do_diagnostics: bool = True,
            full_diagnostics: bool = False,
            keep_model_units: bool = False,
            aggressive_cleanup: bool = True,
    ):
        # Basic attributes
        self.setup = setup
        self.name = name
        self.verbose = verbose
        self.debug_budgets = debug_budgets
        self.do_diagnostics = do_diagnostics
        self.full_diagnostics = full_diagnostics
        self.keep_model_units = keep_model_units
        self.aggressive_cleanup = aggressive_cleanup
        self.error = False
        self.euler = euler

        # Default calibrated variables (typically used in optimization)
        self.calibrated_vars = [
            'Phy_C',
            'Phy_Chl',
            'DIN_concentration',
            'DIP_concentration',
            'DSi_concentration',
            'SPMC'
        ]
        self.output_config = output_config or varinfos.doutput
        self.output_vars = output_vars or list(self.output_config.keys())
        self.toutputoffset = output_time_offset
        self.config = config_dict.copy()

        # Store aggregate_vars for initialization
        self.aggregate_vars_init = aggregate_vars

        # Performance tracking
        self.perf_stats = defaultdict(list)
        self.start_time = time.time()

        # Initialize all components
        self._initialize_all()

        # Optional spin-up phase before main simulation
        if self._has_spinup_components():
            self._run_spinup_phase()

        # Run model and process results
        self._run_model()
        self._process_results()

        # Print performance statistics if verbose
        if verbose:
            self._report_performance()

    def _initialize_all(self):
        """Initialize all model components and parameters"""
        self._initialize_time_parameters()
        self._initialize_tracking_variables()
        self._initialize_components()
        self._setup_component_couplings()
        self._precompute_performance_optimizations()
        self.initial_state = self._create_initial_state_vector()

    def _initialize_time_parameters(self):
        """Initialize time-related parameters"""
        self.dates = self.setup.dates
        if self.setup.dt2 is not None:
            self.used_dt = self.setup.dt2
            self.dt_ratio = round(self.setup.dt / self.setup.dt2)
            self.dates1 = self.dates[::self.dt_ratio]
            self.two_dt = True
            self.dates1_set = set(self.dates1)
        else:
            self.dates1 = None
            self.used_dt = self.setup.dt
            self.dt_ratio = None
            self.two_dt = False
        self.t_span = self.setup.t_span

    def _initialize_tracking_variables(self):
        """Initialize variables for tracking model state"""
        self.nsubcomponents = 0
        self.pool_names = []
        self.generic_pool_names = []
        self.component_ids = []
        self.pool_ids = []
        self.ipools = {}
        self.pool_indices = {}

        # Track special components
        self.phyto_components = set()
        self.phyto_instances = []

        # Diagnostic tracking
        self.diag_pool_names = []
        self.diag_indices = {}
        self.n_diagnostics = None
        self.diag_component_ids = []

        # Initialize aggregate variables
        default_aggs = ['C_tot', 'N_tot', 'P_tot', 'Si_tot', 'Chl_tot',
                       'Cphy_tot', 'POC', 'PON', 'POP', 'DOC']
        self.aggregate_vars = {
            k: [] for k in (self.aggregate_vars_init or default_aggs)
        }

    @_track_time
    def _initialize_components(self):
        """Initialize each model component"""
        self.components = {}
        current_idx = 0

        for key, cfg in self.config.items():
            if key == 'formulation':  # Skip the formulation key
                continue

            instance = cfg['class'](name=key, **cfg.get('parameters', {}))
            instance.formulation = self.config.get('formulation', 'default')
            instance.set_ICs(**cfg.get('initialization', {}))

            # Get pools and update tracking
            pools = list(cfg.get('initialization', {}).keys())
            npools = len(pools)

            # Update tracking variables
            self._update_tracking_variables(key, instance, pools, npools, current_idx)

            # Handle diagnostics
            if self.do_diagnostics:
                diags = cfg.get('diagnostics', [])
                self._setup_diagnostics(instance, diags, key)

            # Track phytoplankton components
            if isinstance(instance, phyto.Phyto):
                self.phyto_components.add(key)
                self.phyto_instances.append(instance)

            self.components[key] = instance
            current_idx += npools

    def _update_tracking_variables(self, key, instance, pools, npools, current_idx):
        """Update component tracking variables"""
        self.nsubcomponents += npools
        self.pool_names.extend([f"{key}_{pool}" for pool in pools])
        self.generic_pool_names.extend([
            f"{instance.classname}_{key * ('concentration' in pool)}{pool}"
            for pool in pools
        ])
        self.component_ids.extend([key] * npools)
        self.pool_ids.extend(pools)
        self.ipools[key] = np.arange(current_idx, current_idx + npools)

        for idx, name in enumerate(self.pool_names[-npools:]):
            self.pool_indices[name] = current_idx + idx

    def _recursive_diagnostic_attrs(self, obj, prefix=''):
        """Recursively find diagnostic attributes."""
        diagnostic_attrs = []

        # Handle nested objects like Elms
        if hasattr(obj, '__dict__'):
            for attr, value in vars(obj).items():
                # Skip private attributes and setup
                if attr.startswith('_') or attr.startswith('coupled') or attr in ('setup', 'diagnostics'):
                    continue

                # Handle nested objects
                if hasattr(value, '__dict__'):
                    nested_attrs = self._recursive_diagnostic_attrs(value, prefix=f"{prefix}{attr}.")
                    diagnostic_attrs.extend(nested_attrs)

                # Check for numeric or None values
                elif isinstance(value, (int, float, np.number, type(None))):
                    diagnostic_attrs.append(f"{prefix}{attr}")

        return diagnostic_attrs

    def _setup_diagnostics(self, instance, diagnostics, key):
        """Setup diagnostic variables for a component"""
        # Only set diagnostics if explicitly defined or full_diagnostics is True
        if diagnostics and not self.full_diagnostics:
            # Use explicitly defined diagnostics
            instance.diagnostics = diagnostics
        elif self.full_diagnostics:
            # Auto-discover diagnostic attributes only in full_diagnostics mode
            instance.diagnostics = self._recursive_diagnostic_attrs(instance)
        else:
            # If no diagnostics explicitly defined and not in full mode, set empty list
            instance.diagnostics = []

        if instance.diagnostics:
            self.diag_component_ids.extend([key] * len(instance.diagnostics))
            self.diag_indices[key] = np.where(np.array(self.diag_component_ids) == key)[0]
            self.diag_pool_names.extend([f"{key}_{diag}" for diag in instance.diagnostics])

    @_track_time
    def _setup_component_couplings(self):
        """Setup couplings between components"""
        for key, cfg in self.config.items():
            if key == 'formulation':
                continue

            component = self.components[key]
            couplings = self._process_couplings(cfg.get('coupling', {}))
            component.setup = self.setup
            component.set_coupling(**couplings)

            # Register aggregates
            if 'aggregate' in cfg:
                for agg_key, agg_value in cfg['aggregate'].items():
                    self.aggregate_vars[agg_key].append(f"{key}_{agg_value}")

    def _process_couplings(self, coupling_dict):
        """Process coupling configurations"""
        couplings = {}
        for coupling_key, coupled_component in coupling_dict.items():
            if coupled_component is None:
                # Skip explicitly disabled couplings
                continue
            if isinstance(coupled_component, list):
                couplings[coupling_key] = [
                    self.components[name] for name in coupled_component
                ]
            else:
                couplings[coupling_key] = self.components[coupled_component]
        return couplings

    def _precompute_performance_optimizations(self):
        """Precompute arrays and mappings for better performance"""
        # Convert poolID to array for faster indexing
        self.pool_id_array = np.array(self.pool_ids)

        # Pre-compute indices for totals
        self.chl_tot_indices = [
            self.pool_indices[name]
            for name in self.aggregate_vars['Chl_tot']
        ]
        self.cphy_tot_indices = [
            self.pool_indices[name]
            for name in self.aggregate_vars['Cphy_tot']
        ]

        # Pre-compute component update mappings
        self.component_update_maps = {}
        for key, component in self.components.items():
            indices = self.ipools[key]
            pool_names = self.pool_id_array[indices]
            self.component_update_maps[key] = {
                str(name): idx for name, idx in zip(pool_names, indices)
            }

        # Pre-compute fast component indices for two_dt case
        if self.two_dt:
            self.fast_components = [comp for comp in self.components.values()
                                    if hasattr(comp, 'dt2') and comp.dt2]
            self.slow_components = [comp for comp in self.components.values()
                                    if not hasattr(comp, 'dt2') or not comp.dt2]

            # Pre-compute dt factors
            self.dt_factors = []
            for comp in self.components.values():
                # Get basic timestep factor
                factor = (1 if hasattr(comp, 'dt2') and comp.dt2
                          else self.dt_ratio)
                # Multiply by component's time conversion factor
                factor *= comp.time_conversion_factor
                self.dt_factors.extend([factor] * len(self.ipools[comp.name]))
            self.dt_factors = np.array(self.dt_factors)

            # Pre-compute zero arrays for slow components
            self.slow_zeros = {
                comp.name: np.zeros(len(self.ipools[comp.name]))
                for comp in self.slow_components
            }

    def _create_initial_state_vector(self):
        """Create initial state vector from all component ICs"""
        return np.concatenate([
            comp.ICs for comp in self.components.values()
        ])

    @_track_time
    def _compute_derivatives(self, t: float, y: np.ndarray, t_idx: int = None) -> np.ndarray:
        """Highly vectorized derivative computation."""
        # Update component states
        for key, component in self.components.items():
            update_map = self.component_update_maps[key]
            component.update_val(
                t=t,
                t_idx=t_idx,
                debugverbose=self.debug_budgets,
                **{name: y[idx] for name, idx in update_map.items()}
            )

        # Update totals
        self.setup.Chl_tot = y[self.chl_tot_indices].sum()
        self.setup.Cphy_tot = y[self.cphy_tot_indices].sum()

        if self.two_dt:
            return self._compute_two_timestep_derivatives(t, t_idx)
        return self._compute_single_timestep_derivatives(t, t_idx)

    @_track_time
    def _compute_single_timestep_derivatives(self, t: float, t_idx: int = None) -> np.ndarray:
        """Compute derivatives for single timestep case"""
        for component in self.components.values():
            try:
                component.get_coupled_processes_indepent_sinks_sources(t, t_idx=t_idx)
            except AttributeError:
                pass

        sources = np.hstack([
            comp.get_sources(t, t_idx=t_idx)
            for comp in self.components.values()
        ])

        sinks = np.hstack([
            comp.get_sinks(t, t_idx=t_idx)
            for comp in self.components.values()
        ])

        return sources - sinks

    @_track_time
    def _compute_two_timestep_derivatives(self, t: float, t_idx: int = None) -> np.ndarray:
        """Compute derivatives for two timestep case"""
        if t in self.dates1_set:
            active_components = self.components.values()
        else:
            active_components = self.fast_components

        for component in active_components:
            try:
                component.get_coupled_processes_indepent_sinks_sources(t, t_idx=t_idx)
            except AttributeError:
                continue

        if t in self.dates1_set:
            sources = np.hstack([
                comp.get_sources(t, t_idx=t_idx)
                for comp in self.components.values()
            ])
            sinks = np.hstack([
                comp.get_sinks(t, t_idx=t_idx)
                for comp in self.components.values()
            ])
        else:
            sources = np.hstack([
                comp.get_sources(t, t_idx=t_idx) if comp in self.fast_components
                else self.slow_zeros[comp.name]
                for comp in self.components.values()
            ])
            sinks = np.hstack([
                comp.get_sinks(t, t_idx=t_idx) if comp in self.fast_components
                else self.slow_zeros[comp.name]
                for comp in self.components.values()
            ])

            if np.all(sources == 0) and np.all(sinks == 0):
                return np.zeros_like(sources)


        return (sources - sinks) * self.dt_factors

    def _has_spinup_components(self) -> bool:
        """Check if any components require spin-up"""
        spinup_found = any(getattr(comp, 'spinup_days', 0) > 0 for comp in self.components.values())
        if self.verbose and spinup_found:
            spinup_components = [comp.name for comp in self.components.values()
                               if getattr(comp, 'spinup_days', 0) > 0]
            print(f'Found components requiring spin-up: {spinup_components}')
        return spinup_found

    @_track_time
    def _run_spinup_phase(self) -> None:
        """Run spin-up phase for components that require it"""
        spinup_components = [comp for comp in self.components.values()
                            if getattr(comp, 'spinup_days', 0) > 0]

        if not spinup_components:
            return

        max_spinup_days = max(comp.spinup_days for comp in spinup_components)

        if self.verbose:
            print(f'Running spin-up phase for {max_spinup_days} days for components: {[comp.name for comp in spinup_components]}')

        # Set spin-up flag
        self.setup.in_spinup_phase = True

        # Create spin-up time array
        spinup_steps = int(max_spinup_days / self.used_dt)
        spinup_dates = np.linspace(0, max_spinup_days, spinup_steps + 1)

        # Run spin-up using existing integration machinery
        y_spinup = self.initial_state.copy()

        for t_spinup in spinup_dates[1:]:
            # Map spin-up time to setup time index (use first setup time for all spin-up)
            t_setup = self.setup.dates[0]

            # Update ONLY spin-up component states
            for component in spinup_components:
                key = component.name
                update_map = self.component_update_maps[key]
                component.update_val(
                    t=t_setup,  # Use setup time instead of spinup time
                    debugverbose=False,
                    **{name: y_spinup[idx] for name, idx in update_map.items()}
                )

            # Update totals (only if needed by spinup components)
            self.setup.Chl_tot = y_spinup[self.chl_tot_indices].sum()
            self.setup.Cphy_tot = y_spinup[self.cphy_tot_indices].sum()

            # Compute derivatives ONLY for spin-up components
            derivatives = np.zeros_like(y_spinup)
            for component in spinup_components:
                key = component.name
                component_indices = self.ipools[key]
                component_sources = component.get_sources(t_setup)
                component_sinks = component.get_sinks(t_setup)

                # Apply time conversion factors for this component
                if self.two_dt:
                    factor = (1 if hasattr(component, 'dt2') and component.dt2 else self.dt_ratio)
                    factor *= component.time_conversion_factor
                    derivatives[component_indices] = (component_sources - component_sinks) * factor
                else:
                    derivatives[component_indices] = component_sources - component_sinks

            y_spinup = y_spinup + self.used_dt * derivatives

            if self.verbose and t_spinup in spinup_dates[::int(1 / self.used_dt)]:
                print(f'Spin-up integration for t = {t_spinup:.1f} days')

        self.initial_state = y_spinup

        # Reset spin-up flag
        self.setup.in_spinup_phase = False

        if self.verbose:
            print(f'Spin-up phase completed. Starting main simulation.')

    @_track_time
    def _compute_diagnostics(self, t: float, t_idx: int = None) -> np.ndarray:
        """Compute diagnostic variables"""
        if self.two_dt:
            if t in self.dates1_set:
                return np.hstack([
                    comp.get_diagnostic_variables()
                    for comp in self.components.values()
                    if hasattr(comp, 'get_diagnostic_variables')
                       and comp.diagnostics
                ])
            else:
                return np.hstack([
                    comp.get_diagnostic_variables()
                    if comp in self.fast_components
                    else np.full(len(self.diag_indices[comp.name]), np.nan)
                    for comp in self.components.values()
                    if hasattr(comp, 'get_diagnostic_variables')
                       and comp.diagnostics
                ])
        else:
            return np.hstack([
                comp.get_diagnostic_variables()
                for comp in self.components.values()
                if hasattr(comp, 'get_diagnostic_variables')
            ])

    @_track_time
    def _run_model(self) -> None:
        """Run the model using either Euler or ODE solver integration."""
        if self.verbose:
            print(f'Starting simulation {self.name}')

        if self.euler:
            self._run_euler_integration()
        else:
            self._run_ode_integration()

    @_track_time
    def _run_euler_integration(self) -> None:
        """Run model using Euler integration with pre-allocated arrays"""
        n_steps = len(self.dates)
        n_vars = len(self.initial_state)

        # Pre-allocate result arrays
        states = np.empty((n_steps, n_vars), dtype=np.float64)
        states[0] = self.initial_state

        if self.do_diagnostics:
            n_diags = len(self.diag_pool_names)
            diagnostics = np.empty((n_steps, n_diags), dtype=np.float64)
            diagnostics[0] = self._compute_diagnostics(self.setup.dates[0], t_idx=0)

        y = self.initial_state.copy()

        for t_idx, t in enumerate(self.dates[1:], start=1):
            derivatives = self._compute_derivatives(t, y, t_idx=t_idx)

            if np.isnan(derivatives).any():
                if self.verbose:
                    print(f'STOP MODEL: NaN values in derivatives: {derivatives}')
                states[t_idx:] = np.nan
                if self.do_diagnostics:
                    diagnostics[t_idx:] = np.nan
                self.error = True
                self.name += '-ERROR'
                break

            y = y + self.used_dt * derivatives
            states[t_idx] = y

            if self.do_diagnostics:
                diagnostics[t_idx] = self._compute_diagnostics(t, t_idx=t_idx)

            if self.verbose and t in self.dates[::int(1 / self.used_dt)]:
                print(f'Eulerian integration for t = {t}')

        self.t = self.dates
        self.y = states.T

        if self.do_diagnostics:
            self.diagnostics = diagnostics.T

    @_track_time
    def _run_ode_integration(self) -> None:
        """Run model using ODE solver"""
        try:
            results = integrate.solve_ivp(
                self._compute_derivatives,
                self.t_span,
                self.initial_state,
                method='DOP853',
                t_eval=self.dates
            )
            self.t = self.dates  # Use DatetimeIndex instead of results.t
            self.y = results.y

        except Exception as e:
            print(f'Error with {self.name}: {e}')
            self.error = True
            self.name += '-ERROR'
            self.t = self.dates
            self.y = np.nan

    @_track_time
    def _process_results(self) -> None:
        """Process model results into a pandas DataFrame with proper units and aggregated variables."""
        # Create base DataFrame with padding if needed
        yvals = self._pad_results(self.y.T)

        # Create initial DataFrame
        self.df = pd.DataFrame(yvals, index=self.t, columns=self.pool_names)

        # Store model units (temporarily needed for derived variable calculations)
        model_data = self.df.values
        self.df = pd.concat([
            self.df,
            pd.DataFrame(model_data, index=self.t, columns=[f'm{col}' for col in self.df.columns])
        ], axis=1)

        # Compute aggregate variables
        self._compute_aggregate_variables()
        # Add diagnostics if available
        if self.do_diagnostics:
            self._add_diagnostics()
        # Compute derived variables
        self._compute_derived_variables()

        # Batch process transformations
        transform_factors = np.array([
            self.output_config.get(pool, {}).get('trsfrm', 1)
            for pool in self.df.columns
        ])
        # self.df.iloc[:, :len(self.pool_names)] *= transform_factors
        self.df.iloc[:, :] *= transform_factors

        # Remove model-unit columns to save memory if not needed
        self._cleanup_dataframe()




    @_track_time
    def _pad_results(self, yvals: np.ndarray) -> np.ndarray:
        """Pad results with NaN values if model was interrupted."""
        if self.error:
            return np.pad(
                yvals,
                ((0, len(self.t) - yvals.shape[0]), (0, 0)),
                'constant',
                constant_values=np.nan
            )
        return yvals

    @_track_time
    def _compute_aggregate_variables(self) -> None:
        """Compute aggregate variables using optimized numpy operations."""
        # Pre-compute column indices for each aggregate variable
        agg_indices = {
            agg_var: [self.df.columns.get_loc(comp) for comp in components]
            for agg_var, components in self.aggregate_vars.items()
        }

        # Compute all aggregates at once using numpy operations
        df_values = self.df.values
        for agg_var, indices in agg_indices.items():
            self.df[agg_var] = df_values[:, indices].sum(axis=1)

    @_track_time
    def _compute_derived_variables(self) -> None:
        """Compute variables defined by expressions in output configuration."""
        for var in set(self.output_vars) - set(self.df.columns):
            if expression := self.output_config.get(var, {}).get('oprt'):
                try:
                    self.df[var] = fns.eval_expr(
                        expression,
                        subdf=self.df,
                        fulldf=self.df,
                        setup=self.setup,
                        model=self
                    )
                except KeyError:
                    if self.verbose:
                        print(f'KeyError in output preparation for {var}')
                    self.df[var] = np.nan

    @_track_time
    def _add_diagnostics(self) -> None:
        """Add diagnostic variables to the results DataFrame efficiently."""
        if not hasattr(self, 'diagnostics'):
            return
        ydiags = self._pad_results(self.diagnostics.T)
        diag_df = pd.DataFrame(ydiags, index=self.t, columns=self.diag_pool_names)
        with pd.option_context('future.no_silent_downcasting', True):
            diag_df = diag_df.bfill()
        # Skip diagnostic columns that already exist in main df to avoid duplicates
        diag_df = diag_df[[col for col in diag_df.columns if col not in self.df.columns]]
        self.df = pd.concat([self.df, diag_df], axis=1)

    def _cleanup_dataframe(self) -> None:
        """
        Remove model-unit columns (m-prefixed) and optionally raw solver results
        to aggressively reduce memory footprint.
        """
        memory_saved = 0

        # Step 1: Remove m-prefixed columns (model units)
        if not self.keep_model_units:
            m_cols = [col for col in self.df.columns
                     if col.startswith('m') and col[1:] in self.df.columns]
            if m_cols:
                self.df.drop(columns=m_cols, inplace=True)
                memory_saved += len(m_cols) * self.df.shape[0] * 8  # 8 bytes per float64
                if self.verbose:
                    print(f'Removed {len(m_cols)} model-unit columns')

        # Step 2: Aggressive cleanup - remove raw solver results
        if self.aggressive_cleanup:
            # Delete raw solver output arrays (self.y, self.diagnostics)
            # These are no longer needed after _process_results() completes
            if hasattr(self, 'y'):
                y_size = self.y.nbytes if hasattr(self.y, 'nbytes') else 0
                memory_saved += y_size
                del self.y

            if hasattr(self, 'diagnostics'):
                diag_size = self.diagnostics.nbytes if hasattr(self.diagnostics, 'nbytes') else 0
                memory_saved += diag_size
                del self.diagnostics

            if self.verbose and memory_saved > 0:
                memory_mb = memory_saved / (1024 ** 2)
                print(f'Aggressive cleanup: freed ~{memory_mb:.1f} MB of memory')

    def get_model_summary(self) -> Dict[str, Any]:
        """Get comprehensive model summary for logging"""
        return {
            'runtime': time.time() - self.start_time,
            'component_count': len(self.components),
            'variable_count': len(self.pool_names),
            'diagnostic_count': self.n_diagnostics,
            'performance': dict(self.perf_stats),
            'error_status': self.error
        }

    def _report_performance(self):
        """Report performance statistics"""
        print("\nPerformance Statistics:")
        print("-" * 80)
        for func_name, times in self.perf_stats.items():
            avg_time = np.mean(times)
            total_time = np.sum(times)
            calls = len(times)
            print(f"{func_name:30s}: {total_time:8.3f}s total, {avg_time * 1000:8.3f}ms/call ({calls} calls)")


    def get_likelihood(self, obs, verbose=False, calibrated_vars=None, **kwargs):
        if self.error:
            return None
        lnl = evaluation.calculate_likelihood(self.df, obs, calibrated_vars=calibrated_vars, verbose=verbose, **kwargs)
        return lnl

