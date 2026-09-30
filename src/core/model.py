import time
import warnings
from functools import wraps
from collections import defaultdict
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Any
from scipy import integrate

from src.utils import functions as fns
from src.core import phys
from src.utils import evaluation
from src.core import legacy
from src.config_model import varinfos

# Precision of the computation. The `dtype` of a Model only sets the precision in which
# its outputs are stored.
COMPUTE_DTYPE = np.float64


def _track_time(func):
    """Decorator accumulating [total seconds, call count] per method.

    Accumulated rather than appended: the derivative methods are called once per
    timestep, i.e. hundreds of thousands of times per run.
    """
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        start = time.time()
        result = func(self, *args, **kwargs)
        stats = self.perf_stats[func.__name__]
        stats[0] += time.time() - start
        stats[1] += 1
        return result
    return wrapper

class _ModelUnitNames:
    """Column lookup for derived expressions during _process_results: 'm<pool>' (model
    units) resolves to the pool column itself, which still holds model units at that
    point, since the unit transform comes after. Replaces a physical duplicate of every
    pool. Only valid before _transform_units: compute_derived_variables, run after the
    fact on converted columns, keeps plain self.df lookup."""

    def __init__(self, df: pd.DataFrame, pool_names: List[str]):
        self.df = df
        self.pools = set(pool_names)

    def __getitem__(self, name):
        if name not in self.df.columns and name[:1] == 'm' and name[1:] in self.pools:
            return self.df[name[1:]]
        return self.df[name]


class Model:
    def __init__(
            self,
            config_dict: Dict[str, Any],
            setup: phys.Setup = phys.Setup(dt=0.001, tmax=20),
            dtype: type = np.float64,
            euler: bool = True,
            output_config: Optional[Dict] = None,
            output_vars: Optional[List[str]] = None,
            name: str = 'Model',
            verbose: bool = True,
            verbose_print_tstep: float = 15.0,
            force_performance_report: bool = False,
            debug_budgets: bool = False,
            aggregate_vars: Optional[List[str]] = None,
            do_diagnostics: bool = True,
            full_diagnostics: bool = False,
            debug_mode: bool = False,
            debug_mode_check_interval: float = 5.0,
            keep_model_units: bool = False,
            aggressive_cleanup: bool = True,
            output_fast_steps: bool = False,
    ):
        """
        dtype: precision in which the outputs are stored (float32 halves their memory and
            disk footprint). The computation itself is always done in float64.
        output_fast_steps: two-timestep scheme only. Record the outputs at every fast step
            (dt2) instead of every slow step (dt): ten times more rows, for a time
            resolution that the slow step (14.4 min at dt = 0.01 d) already covers.
        """
        # Basic attributes
        self.setup = setup
        self.dtype = dtype
        self.name = name
        self.verbose = verbose
        self.verbose_print_tstep = verbose_print_tstep
        self.debug_budgets = debug_budgets
        self.do_diagnostics = do_diagnostics
        self.full_diagnostics = full_diagnostics
        self.debug_mode = debug_mode
        self.debug_mode_check_interval = debug_mode_check_interval
        self.keep_model_units = keep_model_units
        self.aggressive_cleanup = aggressive_cleanup
        self.error = False
        self.euler = euler
        self.output_fast_steps = output_fast_steps

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
        self.config = config_dict.copy()

        # Store aggregate_vars for initialization
        self.aggregate_vars_init = aggregate_vars

        # Performance tracking
        self.perf_stats = defaultdict(lambda: [0.0, 0])
        self.start_time = time.time()

        # Initialize all components
        self._initialize_all()

        # Run model and process results
        self._run_model()
        self._process_results()

        # Print performance statistics if verbose
        if verbose or force_performance_report:
            self._report_performance()

    def _initialize_all(self):
        """Initialize all model components and parameters"""
        self._initialize_time_parameters()
        self._initialize_tracking_variables()
        self._initialize_components()
        self._setup_component_couplings()
        self._initialize_vertical_coupling_ratios()
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

            # Pre-compute boolean mask for slow timesteps (optimization)
            self.is_slow_dt = np.zeros(len(self.dates), dtype=bool)
            self.is_slow_dt[::self.dt_ratio] = True
        else:
            self.dates1 = None
            self.used_dt = self.setup.dt
            self.dt_ratio = None
            self.two_dt = False
            self.is_slow_dt = None
        self.t_span = self.setup.t_span

    def _initialize_tracking_variables(self):
        """Initialize variables for tracking model state"""
        self.pool_names = []
        self.pool_ids = []
        self.ipools = {}
        self.pool_indices = {}

        # Diagnostic tracking
        self.diag_pool_names = []
        self.diag_indices = {}
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
            if key in legacy.NON_COMPONENT_KEYS:
                continue

            parameters = legacy.translate_parameters(key, cfg.get('parameters', {}))
            instance = cfg['class'](name=key, dtype=COMPUTE_DTYPE, **parameters)
            instance.set_ICs(**cfg.get('initialization', {}))

            # Get pools and update tracking
            pools = list(cfg.get('initialization', {}).keys())
            npools = len(pools)

            # Update tracking variables
            self._update_tracking_variables(key, pools, npools, current_idx)

            # Handle diagnostics
            if self.do_diagnostics:
                diags = cfg.get('diagnostics', [])
                self._setup_diagnostics(instance, diags, key)

            self.components[key] = instance
            current_idx += npools

    def _update_tracking_variables(self, key, pools, npools, current_idx):
        """Update component tracking variables"""
        self.pool_names.extend([f"{key}_{pool}" for pool in pools])
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
            if key in legacy.NON_COMPONENT_KEYS:
                continue

            component = self.components[key]
            couplings = self._process_couplings(cfg.get('coupling', {}))
            component.setup = self.setup
            component.debug_mode = self.debug_mode
            component.set_coupling(**couplings)

            # Register aggregates
            if 'aggregate' in cfg:
                for agg_key, agg_value in cfg['aggregate'].items():
                    self.aggregate_vars[agg_key].append(f"{key}_{agg_value}")

    def _initialize_vertical_coupling_ratios(self):
        """Initialize BGC/floc ratios for vertical coupling after all components are set up.

        Applies vertical_coupling_state from config (if present) before initialization,
        allowing restart from extracted simulation state.
        """
        for key, cfg in self.config.items():
            if key in legacy.NON_COMPONENT_KEYS or key not in self.components:
                continue
            component = self.components[key]
            # Apply pre-set state from config (e.g., extracted from previous simulation)
            if 'vertical_coupling_state' in cfg and hasattr(component, 'set_vertical_coupling_state'):
                component.set_vertical_coupling_state(cfg['vertical_coupling_state'])
            # Initialize (respects pre-set values)
            if hasattr(component, 'initialize_vertical_coupling_ratios'):
                component.initialize_vertical_coupling_ratios()

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
        """Precompute, once per run, everything the timestep loop would otherwise rebuild
        at each of its ~10^6 iterations."""
        # Convert poolID to array for faster indexing
        self.pool_id_array = np.array(self.pool_ids)

        # Pre-compute component update mappings
        self.component_update_maps = {}
        for key, component in self.components.items():
            indices = self.ipools[key]
            pool_names = self.pool_id_array[indices]
            self.component_update_maps[key] = {
                str(name): idx for name, idx in zip(pool_names, indices)
            }

        # Layout of the state vector, in component order: (component, [(pool name, index)],
        # slice of the component's pools). Drives both the state update and the assembly
        # of the derivatives.
        self.state_layout = [
            (comp, list(self.component_update_maps[key].items()),
             slice(self.ipools[key][0], self.ipools[key][-1] + 1))
            for key, comp in self.components.items()
        ]
        # Buffers into which the components write their sources and sinks
        n_pools = len(self.pool_names)
        self._sources = np.zeros(n_pools, dtype=COMPUTE_DTYPE)
        self._sinks = np.zeros(n_pools, dtype=COMPUTE_DTYPE)

        # Components exposing a pre-coupling step, and those carrying diagnostics:
        # resolved once instead of being probed at every derivative evaluation.
        self.precoupled_components = [
            comp for comp in self.components.values()
            if hasattr(comp, 'get_coupled_processes_indepent_sinks_sources')
        ]
        self.diag_components = [comp for comp in self.components.values() if comp.diagnostics]

        # Two-timestep scheme: what a fast step (dt2) touches
        if self.two_dt:
            self.fast_components = [comp for comp in self.components.values()
                                    if hasattr(comp, 'dt2') and comp.dt2]
            self._fast_ids = fast = {id(comp) for comp in self.fast_components}
            self.precoupled_fast_components = [comp for comp in self.precoupled_components
                                               if id(comp) in fast]

            # Pre-compute dt factors
            self.dt_factors = []
            for comp in self.components.values():
                # Get basic timestep factor
                factor = (1 if hasattr(comp, 'dt2') and comp.dt2
                          else self.dt_ratio)
                # Multiply by component's time conversion factor
                factor *= comp.time_conversion_factor
                self.dt_factors.extend([factor] * len(self.ipools[comp.name]))
            self.dt_factors = np.array(self.dt_factors, dtype=COMPUTE_DTYPE)

            # The fast pools, gathered into a sub-vector: its layout (same as state_layout,
            # slices now pointing into the sub-vector), indices in the full state vector,
            # dt factors, and source/sink buffers.
            self.fast_layout = []
            start = 0
            for comp, items, pools in self.state_layout:
                if id(comp) in fast:
                    n = pools.stop - pools.start
                    self.fast_layout.append((comp, items, slice(start, start + n)))
                    start += n
            self.fast_indices = np.array([idx for comp in self.fast_components
                                          for idx in self.ipools[comp.name]], dtype=int)
            self.fast_dt_factors = self.dt_factors[self.fast_indices]
            # Where the fast pools sit in the state vector: a slice when they are
            # contiguous (the usual case), which numpy reads and writes faster
            contiguous = len(self.fast_indices) and np.all(np.diff(self.fast_indices) == 1)
            self._fast_slots = (slice(self.fast_indices[0], self.fast_indices[-1] + 1)
                                if contiguous else self.fast_indices)
            self._fast_sources = np.zeros(start, dtype=COMPUTE_DTYPE)
            self._fast_sinks = np.zeros(start, dtype=COMPUTE_DTYPE)
            self.fast_diag_components = [comp for comp in self.diag_components if id(comp) in fast]

    def _create_initial_state_vector(self):
        """Create initial state vector from all component ICs"""
        return np.concatenate([
            comp.ICs for comp in self.components.values()
        ])

    def _update_components(self, layout, t, t_idx, y):
        """Hand each component of `layout` its pools, read from the state vector y.

        As Python floats: the components compute on scalars, and arithmetic on numpy
        scalars is several times slower."""
        values = y.tolist()
        for component, pools, _ in layout:
            component.update_val(t=t, t_idx=t_idx, debugverbose=self.debug_budgets,
                                 **{name: values[idx] for name, idx in pools})

    def _compute_derivatives(self, t: float, y: np.ndarray, t_idx: int = None) -> np.ndarray:
        """Derivatives of the whole state vector: every component is updated and evaluated.

        This is a full step. In the two-timestep scheme it is the slow step, and the dt
        factors then scale each pool to its own timestep. The fast steps in between, where
        only the dt2 components move, are done by _fast_step.
        """
        self._update_components(self.state_layout, t, t_idx, y)
        for component in self.precoupled_components:
            component.get_coupled_processes_indepent_sinks_sources(t, t_idx=t_idx)

        sources, sinks = self._sources, self._sinks
        for component, _, pools in self.state_layout:
            sources[pools] = component.get_sources(t, t_idx=t_idx)
        for component, _, pools in self.state_layout:
            sinks[pools] = component.get_sinks(t, t_idx=t_idx)

        if self.two_dt:
            return (sources - sinks) * self.dt_factors
        return sources - sinks

    def _fast_step(self, t, t_idx: int, y: np.ndarray, resync: bool) -> np.ndarray:
        """One fast step (dt2) of the two-timestep scheme, on y in place.

        Only the dt2 components are updated, evaluated and integrated: the slow pools do
        not move between two slow steps, so recomputing them (as a full step would) only
        to multiply them by zero is wasted. Returns the derivatives of the fast pools.

        resync: True on the first fast step after a slow step. The slow pools have just
        moved, and the fast components read some of them (Flocs reads TEPC.C through
        coupled_glue): every component is handed its new pools once.
        """
        self._update_components(self.state_layout if resync else self.fast_layout, t, t_idx, y)
        for component in self.precoupled_fast_components:
            component.get_coupled_processes_indepent_sinks_sources(t, t_idx=t_idx)

        sources, sinks = self._fast_sources, self._fast_sinks
        for component, _, pools in self.fast_layout:
            sources[pools] = component.get_sources(t, t_idx=t_idx)
        for component, _, pools in self.fast_layout:
            sinks[pools] = component.get_sinks(t, t_idx=t_idx)
        derivatives = (sources - sinks) * self.fast_dt_factors

        y[self._fast_slots] += self.used_dt * derivatives
        return derivatives

    def _compute_diagnostics(self, components) -> np.ndarray:
        """Current values of the diagnostics of `components`, concatenated."""
        diag_arrays = [comp.get_diagnostic_variables() for comp in components]
        return np.hstack(diag_arrays) if diag_arrays else np.array([])

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
        """Run model using Euler integration with pre-allocated arrays.

        Two-timestep scheme: every dt_ratio-th step is a full (slow) step, all components
        being evaluated; the steps in between are fast steps (_fast_step), where only the
        dt2 components move. The outputs are recorded at the slow steps, or at every step
        with output_fast_steps.
        """
        n_steps = len(self.dates)
        n_vars = len(self.initial_state)

        # Output rows: one every `stride` steps
        stride = self.dt_ratio if (self.two_dt and not self.output_fast_steps) else 1
        out_dates = self.dates[::stride]
        # Fast-step rows, where only the fast components' diagnostics are computed: the
        # others are left NaN and back-filled in _add_diagnostics
        self._diag_backfill = self.two_dt and stride == 1

        # Pre-allocate result arrays (NaN: the rows after an interrupted run stay NaN)
        states = np.full((len(out_dates), n_vars), np.nan, dtype=self.dtype)
        states[0] = self.initial_state

        if self.do_diagnostics:
            diagnostics = np.full((len(out_dates), len(self.diag_pool_names)), np.nan, dtype=self.dtype)
            diagnostics[0] = self._compute_diagnostics(self.diag_components)
            if self._diag_backfill:
                fast_diag_columns = np.array([col for comp in self.fast_diag_components
                                              for col in self.diag_indices[comp.name]], dtype=int)

        y = self.initial_state.astype(COMPUTE_DTYPE, copy=True)
        used_dt = self.used_dt

        # Optimization: Check for NaN every 5 days instead of every timestep
        nan_check_interval = int(5.0 / self.used_dt)
        neg_check_interval = int(self.debug_mode_check_interval / self.used_dt) if self.debug_mode else None

        previous_was_slow = True
        for t_idx, t in enumerate(self.dates[1:], start=1):
            is_slow_step = self.is_slow_dt[t_idx] if self.two_dt else True
            if is_slow_step:
                derivatives = self._compute_derivatives(t, y, t_idx=t_idx)
            else:
                derivatives = self._fast_step(t, t_idx, y, resync=previous_was_slow)
            previous_was_slow = is_slow_step

            # Check for NaN periodically and at final timestep (critical for optimization workflow)
            should_check_nan = (t_idx % nan_check_interval == 0) or (t_idx == n_steps - 1)
            if should_check_nan and np.isnan(derivatives).any():
                if self.verbose:
                    print(f'STOP MODEL: NaN values in derivatives at t_idx={t_idx}')
                self.error = True
                self.name += '-ERROR'
                break

            if is_slow_step:
                y = y + used_dt * derivatives

            # Debug mode: check for negative state variables
            if self.debug_mode:
                if (y < 0).any():
                    neg_vars = [self.pool_names[i] for i in np.where(y < 0)[0]]
                    print(f'Negative state variables at t_idx={t_idx}: {neg_vars}')

                    if (t_idx % neg_check_interval == 0 or t_idx == n_steps - 1):
                        print(f'STOP MODEL: Negative state variables at t_idx={t_idx}: {neg_vars}')
                        self.error = True
                        self.name += '-ERROR'
                        break

            if t_idx % stride == 0:
                row = t_idx // stride
                states[row] = y
                if self.do_diagnostics:
                    if is_slow_step:
                        diagnostics[row] = self._compute_diagnostics(self.diag_components)
                    elif len(fast_diag_columns):
                        diagnostics[row, fast_diag_columns] = self._compute_diagnostics(self.fast_diag_components)

            if self.verbose and t_idx % int(self.verbose_print_tstep / self.used_dt) == 0:
                print(f'Eulerian integration for t = {t}')

        self.t = out_dates
        self.y = states.T

        if self.do_diagnostics:
            self.diagnostics = diagnostics.T

    @_track_time
    def _run_ode_integration(self) -> None:
        """Run the model with an adaptive ODE solver instead of the Euler scheme.

        The components read their forcings from arrays indexed by timestep, not by time,
        so the solver's continuous t is mapped back to the nearest index below it. The
        forcings are therefore piecewise constant over a setup step, which is what the
        Euler scheme does too.

        Two limitations, both deliberate:
        - the two-timestep scheme is not supported (its dt_factors weighting only makes
          sense for a fixed step);
        - no diagnostics are produced, the solver not evaluating on the output grid.
        """
        if self.two_dt:
            raise NotImplementedError(
                'euler=False does not support the two-timestep scheme (dt2 is set). '
                'Use a single timestep (dt2=None), or keep euler=True.')

        n_steps = len(self.dates)
        used_dt = self.used_dt

        def derivatives(t, y):
            t_idx = min(int(t / used_dt), n_steps - 1)
            return self._compute_derivatives(t, y, t_idx=t_idx)

        try:
            results = integrate.solve_ivp(
                derivatives,
                self.t_span,
                self.initial_state,
                method='DOP853',
                t_eval=self.setup.t_eval
            )
            if not results.success:
                raise RuntimeError(results.message)
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
        """Process model results into a pandas DataFrame with proper units and aggregated variables.

        Every step only ADDS columns to self.df, never rebuilds it: on a 3-year run at
        dt2 = 1e-3 one full copy of the frame costs ~0.4 GB, and during an optimisation
        that per-worker memory is what bounds how many workers run in parallel. The
        price is a fragmented frame, which pandas warns about; consolidating it would be
        one more full copy, so the warning is silenced here on purpose.
        """
        # Create base DataFrame with padding if needed
        yvals = self._pad_results(self.y.T)
        self.df = pd.DataFrame(yvals, index=self.t, columns=self.pool_names)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', pd.errors.PerformanceWarning)
            self._compute_aggregate_variables(yvals)
            if self.do_diagnostics:
                self._add_diagnostics()
            self._compute_derived_variables()
            self._transform_units()
            self._cleanup_dataframe(yvals)

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
    def _compute_aggregate_variables(self, yvals: np.ndarray) -> None:
        """Compute aggregate variables, summing pool columns of the raw state array
        (not self.df.values, which would copy the whole frame)."""
        pool_index = {name: i for i, name in enumerate(self.pool_names)}
        for agg_var, components in self.aggregate_vars.items():
            self.df[agg_var] = yvals[:, [pool_index[comp] for comp in components]].sum(axis=1)

    def _evaluate_derived_expression(self, expression: str, names=None) -> Any:
        """Evaluate one `oprt` expression against the result DataFrame.

        Shared by _compute_derived_variables (at the end of a run) and by the public
        compute_derived_variables (after the fact). The two differ only in how they
        report failures, so only the evaluation itself lives here. `names` resolves the
        column names used in the expression (default: self.df).
        """
        names = self.df if names is None else names
        values = fns.eval_expr(expression, subdf=names, fulldf=names,
                               setup=self.setup, model=self)
        # Guard against division-by-zero in ratio expressions (e.g. a denominator
        # that is 0 at night) producing +/-inf.
        return pd.Series(values, index=self.df.index).replace([np.inf, -np.inf], np.nan)

    @_track_time
    def _compute_derived_variables(self) -> None:
        """Compute variables defined by expressions in output configuration."""
        names = _ModelUnitNames(self.df, self.pool_names)
        for var in set(self.output_vars) - set(self.df.columns):
            if expression := self.output_config.get(var, {}).get('oprt'):
                try:
                    self.df[var] = self._evaluate_derived_expression(expression, names)
                except KeyError:
                    if self.verbose:
                        print(f'KeyError in output preparation for {var}')
                    self.df[var] = np.nan

    def compute_derived_variables(self, var_list: Optional[List[str]] = None,
                                  output_config: Optional[Dict] = None,
                                  verbose: bool = True) -> None:
        """
        Compute derived variables from expressions, even after model instantiation.

        This method allows recalculating derived variables that were added to varinfos.py
        after the model was run, or variables that were not initially in output_vars.

        Args:
            var_list: List of specific variables to compute. If None, attempts to compute
                     all variables defined in output_config that are missing from df.
            output_config: Alternative output configuration dict. If None, uses self.output_config.
            verbose: Whether to print status messages.

        Example:
            >>> model = sim_manager.load_simulation('my_simulation')
            >>> model.compute_derived_variables(['TEPtoC', 'TEPtoChl'])
            >>> # Now you can plot these variables
        """
        config = output_config or self.output_config

        # Determine which variables to compute
        if var_list is None:
            # Compute all variables in config that have 'oprt' and are missing from df
            vars_to_compute = [
                var for var in config.keys()
                if var not in self.df.columns and config.get(var, {}).get('oprt')
            ]
        else:
            # Compute only requested variables that are missing from df
            vars_to_compute = [var for var in var_list if var not in self.df.columns]

        if not vars_to_compute:
            if verbose:
                print("All requested derived variables are already present in the model dataframe.")
            return

        if verbose:
            print(f"Computing {len(vars_to_compute)} derived variable(s): {vars_to_compute}")

        # Compute each variable
        computed = []
        failed = []
        for var in vars_to_compute:
            if expression := config.get(var, {}).get('oprt'):
                try:
                    self.df[var] = self._evaluate_derived_expression(expression)
                    computed.append(var)
                    if verbose:
                        print(f"  ✓ Successfully computed '{var}'")
                except KeyError as e:
                    failed.append((var, str(e)))
                    if verbose:
                        print(f"  ✗ Failed to compute '{var}': missing dependency {e}")
                    self.df[var] = np.nan
                except Exception as e:
                    failed.append((var, str(e)))
                    if verbose:
                        print(f"  ✗ Failed to compute '{var}': {e}")
                    self.df[var] = np.nan
            else:
                if verbose:
                    print(f"  ⚠ Variable '{var}' has no 'oprt' expression in output_config")

        if verbose and (computed or failed):
            print(f"\nSummary: {len(computed)} computed, {len(failed)} failed")

    @_track_time
    def _add_diagnostics(self) -> None:
        """Add diagnostic variables to the results DataFrame efficiently."""
        if not hasattr(self, 'diagnostics'):
            return
        ydiags = self._pad_results(self.diagnostics.T)
        diag_df = pd.DataFrame(ydiags, index=self.t, columns=self.diag_pool_names)
        if getattr(self, '_diag_backfill', False):
            # Fast-step rows: the slow components' diagnostics take their next value
            with pd.option_context('future.no_silent_downcasting', True):
                diag_df = diag_df.bfill()
        # Skip diagnostic columns that already exist in main df to avoid duplicates.
        # Added column by column: a concat would copy the whole frame.
        for col in diag_df.columns:
            if col not in self.df.columns:
                self.df[col] = diag_df[col].values

    def _transform_units(self) -> None:
        """Convert columns to output units (varinfos 'trsfrm'). Only the few columns whose
        factor is not 1 are touched, instead of multiplying (and copying) the whole frame."""
        factors = np.array([self.output_config.get(col, {}).get('trsfrm', 1)
                            for col in self.df.columns], dtype=self.dtype)
        for col, factor in zip(list(self.df.columns), factors):
            if factor != 1:
                self.df[col] = self.df[col] * factor

    def _cleanup_dataframe(self, yvals: np.ndarray) -> None:
        """
        Optionally add the model-unit columns (m-prefixed), and remove raw solver results
        to aggressively reduce memory footprint.
        """
        memory_saved = 0

        # Step 1: model-unit columns, only on request (derived expressions resolve their
        # m-names without them, cf. _ModelUnitNames)
        if self.keep_model_units:
            for i, name in enumerate(self.pool_names):
                self.df[f'm{name}'] = yvals[:, i].copy()

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
            'diagnostic_count': len(self.diag_pool_names),
            'performance': dict(self.perf_stats),
            'error_status': self.error
        }

    def _report_performance(self):
        """Report performance statistics"""
        print("\nPerformance Statistics:")
        print("-" * 80)
        for func_name, (total_time, calls) in self.perf_stats.items():
            avg_time = total_time / calls if calls else 0.
            print(f"{func_name:30s}: {total_time:8.3f}s total, {avg_time * 1000:8.3f}ms/call ({calls} calls)")


    def get_likelihood(self, obs, verbose=False, calibrated_vars=None, **kwargs):
        if self.error:
            return None
        lnl = evaluation.calculate_likelihood(self, obs, calibrated_vars=calibrated_vars, verbose=verbose, **kwargs)
        return lnl

