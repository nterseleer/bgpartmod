"""
Simulation management utilities for BGC model.
Handles saving, loading, and tracking of Simulations.
"""
import json
import os
import pickle
import re
from datetime import datetime
from enum import Enum

import dill
import numpy as np
import pandas as pd
from deepdiff import DeepDiff
from typing import Optional, List, Dict, Any, Union

from src.utils import functions as fns
from src.config_system import path_config as path_cfg

# Constants

# Optimization log columns (human-readable names)
OPTIMIZATION_LOG_COLUMNS = [
    'ID',
    '#pars',
    '#vars',
    'UserNote_PRE',
    'UserNote_POST',
    'lnl',
    '#Gen',
    'HittingLow',
    'HittingHigh',
    'Duration_h',
    'LastImprovementGen',
    'First_lnl',
    'lnl_improvement',
    'FirstValidGen',
    'BoundHit%',
    'Ref_Opt',
    'PlannedGen',
    'Date',
    'Env'
]

class SimulationTypes(Enum):
    MODEL_RUN = 1
    REFERENCES_SIMULATION = 2
    UNDEFINED = 3


# Create necessary directories
for directory in [path_cfg.SIMULATION_DIR, path_cfg.MODEL_RUNS_DIR, path_cfg.REFERENCES_SIMULATION_DIR]:
    if not os.path.exists(directory):
        os.makedirs(directory)
# Setup parameters to track in logs
SETUP_TRACKED_FIELDS = [
    # Core time settings
    'tmin', 'tmax', 'dt', 'dt2',
    # Physical parameters
    'T', 'z', 'k_att', 'kb', 'pCO2',
    # Light settings
    'light_prop', 'PARfromfile', 'lightfirst',
    # Shear settings
    'g_shear_rate', 'vary_g_shear', 'gshearfact', 'gshearper'
]

# Log columns in desired order
LOG_COLUMNS = [
                  'date',  # YYYY-MM-DD
                  'short_name',  # Name without timestamp
                  'user_notes',  # User provided notes
                  'parent_simulation',  # Parent simulation name if any
                  'configuration_diff',  # Differences from parent config
                  'setup_diff',  # Differences from parent setup
                  'status',  # Simulation status
                  'setup_start_date',  # Simulation start date
                  'setup_end_date',  # Simulation end date
              ] + [
                  # Setup fields - add all tracked fields
                  f'setup_{field}' for field in SETUP_TRACKED_FIELDS
              ] + [
                  # Configuration and runtime fields
                  'config_formulation',
                  'config_components',  # Detailed component list as JSON
                  'variables',  # Model variables
                  'runtime_info',  # Additional runtime information
                  'full_name',  # Complete name with timestamp
                  'additional_notes'  # additional notes about results, simulations or other
              ]


def create_log_file_if_not_exist() -> bool:
    if not os.path.exists(path_cfg.LOG_FILE):
        df = pd.DataFrame(columns=LOG_COLUMNS)
        df.to_csv(path_cfg.LOG_FILE, index=False)


def _load_log() -> pd.DataFrame:
    """Load or create simulation log with organized columns"""
    create_log_file_if_not_exist()
    df = pd.read_csv(path_cfg.LOG_FILE, )
    # Force column order and drop any extra columns
    return df.reindex(columns=LOG_COLUMNS)


def interactive_log_additional_note_modification(simulation_name):

    log_df = _load_log()
    target_simulation_index = _get_taget_additional_note_simulation_index(log_df, simulation_name)

    # avoid : FutureWarning: Setting an item of incompatible dtype is deprecated and will raise in a future error of pandas.
    if log_df["additional_notes"].count() == 0:
        log_df["additional_notes"] = None

    additional_notes = input(
        f'Please complete the “additional notes”  section of the simulation log for the simulation "{simulation_name}" :\n')

    if pd.isna(log_df.at[target_simulation_index[0], "additional_notes"]) or not _get_interactive_yes_no_answer("Do you want to append to the previous note (yes/no), no = overwrite : "):
        log_df.at[target_simulation_index[0], "additional_notes"] = additional_notes
        print("\033[92mAdditional notes for '{}' successfully set !\033[0m".format(simulation_name))

    else :
        previous_notes = log_df.at[target_simulation_index[0], "additional_notes"]
        log_df.at[target_simulation_index[0], "additional_notes"] = f'{previous_notes}\n{additional_notes}'
        print("\033[92mAdditional notes for '{}' successfully appended !\033[0m".format(simulation_name))

    log_df.to_csv(path_cfg.LOG_FILE, index=False)



def _get_taget_additional_note_simulation_index(log_df : pd.DataFrame, simulation_name):

    target_simulation_index = log_df.query(f"full_name == '{simulation_name}'").index
    if len(target_simulation_index) > 1:
        raise ValueError(
            f"Multiple simulations found with full_name '{simulation_name}'. Expected a unique simulation.")
    if len(target_simulation_index) == 0:
        raise ValueError(
            f"No simulations found with full_name '{simulation_name}' found in the log file.")

    return target_simulation_index



# yes = True, No = False
def _get_interactive_yes_no_answer(question : str)-> bool:
    while True:
        user_input = input(question)
        if user_input.lower() in ["yes", "y"]:
            return True
        elif user_input.lower() in ["no", "n"]:
            return False
        else:
            print("Invalid input. Please enter yes/no.")


# Optimization Log Management Functions

def get_optimization_log() -> pd.DataFrame:
    """Load or create optimization log."""
    if not os.path.exists(path_cfg.PRIVATE_OPT_LOG_FILE):
        # Create with header only
        df = pd.DataFrame(columns=OPTIMIZATION_LOG_COLUMNS)
        df.to_csv(path_cfg.PRIVATE_OPT_LOG_FILE, index=False)
        return df
    return pd.read_csv(path_cfg.PRIVATE_OPT_LOG_FILE)


def save_optimization_log(df: pd.DataFrame):
    """Save optimization log."""
    df.to_csv(path_cfg.PRIVATE_OPT_LOG_FILE, index=False)


def get_next_optimization_id() -> str:
    """Get next available optimization ID from the shared log."""
    log_df = get_optimization_log()
    if len(log_df) == 0:
        return 'OPT000'

    # Extract numbers from existing IDs
    ids = log_df['ID'].tolist()
    numbers = [int(id_str[3:]) for id_str in ids if id_str.startswith('OPT')]
    next_num = max(numbers + [-1]) + 1
    return f'OPT{next_num:03d}'


def add_optimization_to_log(opt_id: str, param_count: int, user_note: str, reference_opt: str = None,
                           boundary_hit_threshold_percent: float = 5.0, calibrated_vars_count: int = None):
    """Add new optimization entry to log."""
    log_df = get_optimization_log()

    # Detect environment with SLURM differentiation
    if "ecmwf" in os.getcwd().lower() or "copernicus" in os.getcwd().lower():
        # On ECMWF - check if batch (SLURM) or interactive (JupyterHub)
        if os.environ.get('SLURM_JOB_ID') or os.environ.get('SLURM_JOBID'):
            environment = "ecmwf_sbatch"
        else:
            environment = "ecmwf_jupyt"
    else:
        environment = "local"

    new_entry = pd.DataFrame([{
        'ID': opt_id,
        '#pars': param_count,
        '#vars': calibrated_vars_count if calibrated_vars_count is not None else '',
        'UserNote_PRE': user_note,
        'UserNote_POST': '',
        'lnl': '',
        '#Gen': '',
        'HittingLow': '',
        'HittingHigh': '',
        'Duration_h': '',
        'LastImprovementGen': '',
        'First_lnl': '',
        'lnl_improvement': '',
        'FirstValidGen': '',
        'BoundHit%': boundary_hit_threshold_percent,
        'Ref_Opt': reference_opt,
        'PlannedGen': '',
        'Date': datetime.now().strftime('%Y-%m-%d'),
        'Env': environment
    }])

    log_df = pd.concat([log_df, new_entry], ignore_index=True)
    save_optimization_log(log_df)


def update_optimization_log_convergence(
    opt_id: str,
    likelihood_score: float,
    runtime_minutes: Optional[float],
    generations_completed: int,
    generations_planned: int,
    first_valid_score: float,
    first_valid_gen: Optional[int],
    score_improvement: float,
    last_improvement_gen: int,
    params_at_lower_bound: str,
    params_at_upper_bound: str
):
    """
    Update optimization log with complete results and convergence information.

    Args:
        opt_id: Optimization ID
        likelihood_score: Final best likelihood score
        runtime_minutes: Total runtime in minutes (will be converted to hours in log), None if not available
        generations_completed: Number of generations completed
        generations_planned: Number of generations planned
        first_valid_score: First non-badlnl score
        first_valid_gen: Generation of first valid score (or None)
        score_improvement: Total improvement from first valid to final
        last_improvement_gen: Generation where best score was found
        params_at_lower_bound: Comma-separated list of parameters at lower bound
        params_at_upper_bound: Comma-separated list of parameters at upper bound
    """
    log_df = get_optimization_log()
    mask = log_df['ID'] == opt_id

    if mask.any():
        log_df.loc[mask, 'lnl'] = likelihood_score
        log_df.loc[mask, 'Duration_h'] = (runtime_minutes / 60.0) if runtime_minutes is not None else ''
        log_df.loc[mask, '#Gen'] = generations_completed
        log_df.loc[mask, 'PlannedGen'] = generations_planned
        log_df.loc[mask, 'First_lnl'] = first_valid_score
        log_df.loc[mask, 'FirstValidGen'] = first_valid_gen if first_valid_gen is not None else ''
        log_df.loc[mask, 'lnl_improvement'] = score_improvement
        log_df.loc[mask, 'LastImprovementGen'] = last_improvement_gen
        log_df.loc[mask, 'HittingLow'] = params_at_lower_bound
        log_df.loc[mask, 'HittingHigh'] = params_at_upper_bound
        save_optimization_log(log_df)


def add_post_note(item_id: str, item_type: str, post_note: str = None):
    """Unified post-note function for simulations and optimizations."""
    if item_type.lower() == 'optimization':
        _add_optimization_post_note(item_id, post_note)
    elif item_type.lower() == 'simulation':
        interactive_log_additional_note_modification(item_id)
    else:
        raise ValueError(f"Unknown item_type: {item_type}. Use 'optimization' or 'simulation'")


def _add_optimization_post_note(opt_id: str, post_note: str = None):
    """Add or update post-optimization note."""
    if post_note is None:
        print(f"\n{'='*60}")
        print(f"OPTIMIZATION {opt_id} - Post Analysis Note")
        print("="*60)
        post_note = input("Please add conclusions/observations after analyzing this optimization:\n> ")
        post_note = post_note.strip()

    log_df = get_optimization_log()
    mask = log_df['ID'] == opt_id
    if mask.any():
        log_df.loc[mask, 'UserNote_POST'] = post_note
        save_optimization_log(log_df)
        print(f"✓ Post-analysis note added for {opt_id}")
    else:
        print(f"✗ Optimization {opt_id} not found in log")


# Legacy function for backward compatibility
add_optimization_post_note = _add_optimization_post_note


def _generate_name(base_name: str = None) -> str:
    """Generate unique simulation name with timestamp"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if base_name:
        return f"{timestamp}_{base_name}"
    return timestamp


def _extract_setup_info(setup) -> Dict[str, Any]:
    """
    Extract setup information for logging.

    Args:
        setup: Setup instance

    Returns:
        Dictionary of setup information formatted for logging.
    """
    # Initialize with mandatory fields as empty strings
    setup_info = {f'setup_{field}': '' for field in SETUP_TRACKED_FIELDS}

    # Add date fields explicitly
    setup_info['setup_start_date'] = setup.start_date.strftime('%Y-%m-%d %H:%M:%S')
    setup_info['setup_end_date'] = setup.end_date.strftime('%Y-%m-%d %H:%M:%S')

    # Update with values from physical state
    physical_state = setup.get_physical_state()
    for category, attrs in physical_state.items():
        for name, value in attrs.items():
            if name in SETUP_TRACKED_FIELDS:  # Only include tracked fields
                if isinstance(value, (int, float, str, bool)):
                    setup_info[f'setup_{name}'] = value
                elif isinstance(value, (np.ndarray, pd.Series, pd.DataFrame)):
                    stats = {
                        'mean': float(np.mean(value)),
                        'min': float(np.min(value)),
                        'max': float(np.max(value))
                    }
                    setup_info[f'setup_{name}_stats'] = json.dumps(stats)

    return setup_info


def _get_setup_diff(setup1, setup2) -> str:
    """
    Get human-readable setup differences based on physical state.

    Args:
        setup1: First Setup instance
        setup2: Second Setup instance

    Returns:
        String describing differences between setups or "Idem" if identical.
    """
    if setup1 is None or setup2 is None:
        return ""

    diff_summary = []

    # Only compare attributes defined in PHYSICAL_ATTRIBUTES
    for category, attrs in setup1.PHYSICAL_ATTRIBUTES.items():
        for field in attrs:
            if not (hasattr(setup1, field) and hasattr(setup2, field)):
                continue

            val1 = getattr(setup1, field)
            val2 = getattr(setup2, field)

            if isinstance(val1, (int, float, str, bool)):
                if val1 != val2:
                    diff_summary.append(f"{field}: {val1} → {val2}")
            elif isinstance(val1, (np.ndarray, pd.Series, pd.DataFrame)):
                if not np.array_equal(np.asarray(val1), np.asarray(val2)):
                    mean1 = float(np.mean(val1))
                    mean2 = float(np.mean(val2))
                    if not np.isclose(mean1, mean2):
                        diff_summary.append(f"{field}: mean {mean1:.3g} → {mean2:.3g}")

    return "; ".join(diff_summary) if diff_summary else "Idem"


def _summarize_config(config: Dict) -> Dict:
    """Create detailed config summary for logging"""
    components = {}
    for comp_name, cfg in config.items():
        if comp_name == 'formulation':
            continue

        if 'class' in cfg:
            class_name = cfg['class'].__name__
            if class_name not in components:
                components[class_name] = []
            components[class_name].append(comp_name)

    return {
        'config_formulation': config.get('formulation', 'unknown'),
        'config_components': json.dumps(components)
    }


def _get_config_diff(config1: Dict, config2: Dict) -> str:
    """
    Get human-readable configuration differences.

    Args:
        config1: Original configuration
        config2: New configuration

    Returns:
        String with formatted parameter changes
    """
    return fns.compare_dicts(config1, config2,
                            format_output=True,
                            config_mode=True,
                            print_result=False)




def save_simulation(
        model,
        base_name: Optional[str] = None,
        sim_type: SimulationTypes = SimulationTypes.MODEL_RUN,
        parent_simulation: Optional[str] = None,
        user_notes: Optional[str] = None,
        save_full: bool = False) -> str:
    """
    Save simulation with enhanced logging and parent comparisons.
    """
    # Generate name and directory structure
    timestamp = datetime.now()
    date = timestamp.strftime('%Y-%m-%d')
    short_name = base_name or timestamp.strftime('%H%M%S')
    full_name = f"{timestamp.strftime('%Y%m%d_%H%M%S')}_{short_name}" if base_name else timestamp.strftime(
        '%Y%m%d_%H%M%S')

    sim_dir = os.path.join(
        path_cfg.REFERENCES_SIMULATION_DIR if sim_type == SimulationTypes.REFERENCES_SIMULATION else path_cfg.MODEL_RUNS_DIR,
        full_name
    )
    os.makedirs(sim_dir, exist_ok=True)

    # Load parent information if available
    parent_config = None
    parent_setup = None
    if parent_simulation:
        try:
            parent_data = load_simulation(parent_simulation)
            parent_config = getattr(parent_data, 'config', None)
            parent_setup = getattr(parent_data, 'setup', None)
        except Exception as e:
            print(f"Warning: Could not load parent simulation: {e}")

    # Prepare log entry
    setup_info = _extract_setup_info(model.setup)
    config_info = _summarize_config(model.config)

    log_entry = {
        'date': date,
        'short_name': short_name,
        'user_notes': user_notes or "No notes provided",
        'parent_simulation': parent_simulation or "",
        'configuration_diff': _get_config_diff(parent_config, model.config) if parent_simulation else "",
        'setup_diff': _get_setup_diff(parent_setup, model.setup) if parent_simulation else "",
        'status': 'completed',
        'config_formulation': config_info['config_formulation'],
        'config_components': config_info['config_components'],
        'variables': ','.join(sorted(model.df.columns)),
        'runtime_info': f"Type: {sim_type.name}",
        'full_name': full_name,
        **setup_info,  # Add setup info after known fields
    }

    # Always save setup
    with open(os.path.join(sim_dir, 'setup.pkl'), 'wb') as f:
        dill.dump(model.setup, f)

    # Save remaining files
    fns.save_human_readable_config(
        model.config,
        os.path.join(sim_dir, 'config.json')
    )
    with open(os.path.join(sim_dir, 'config.pkl'), 'wb') as f:
        pickle.dump(model.config, f)

    if save_full or sim_type == SimulationTypes.REFERENCES_SIMULATION:
        key_results = {
            'timestamps': model.df.index.values,
            'variables': {col: model.df[col].values for col in model.df.columns}
        }
        np.save(os.path.join(sim_dir, 'results.npy'), key_results)
        with open(os.path.join(sim_dir, 'model.pkl'), 'wb') as f:
            dill.dump(model, f)

    # Update log with strict column order
    create_log_file_if_not_exist()
    new_row = pd.DataFrame([{col: log_entry.get(col, '') for col in LOG_COLUMNS}])

    # log_df = pd.concat([log_df, new_row], ignore_index=True)

    new_row.to_csv(path_cfg.LOG_FILE, na_rep='NA', mode='a', index=False, header=False)

    return full_name


def run_or_load_simulation(config_dict: Dict,
                           setup: Any,
                           name: Optional[str] = None,
                           parent_simulation: Optional[str] = None,
                           user_notes: Optional[str] = None,
                           save_type: SimulationTypes = SimulationTypes.MODEL_RUN,
                           **model_kwargs) -> 'Model':
    """
    Get a simulation by either loading it or creating it.

    Args:
        config_dict: Model configuration dictionary
        setup: Model setup instance
        name: Optional name for the simulation
        parent_simulation: Optional parent simulation name
        user_notes: Optional notes about the simulation
        save_type: UNDEFINED simulation is not saved, REFERENCES_SIMULATION save as a reference simulation (full save),
        MODEL_RUN save as a model_run simulation. If the save_type from the current simulation is different
        **model_kwargs: Additional arguments for Model constructor
    """
    from src.core import model

    loaded_simulation_type = get_saved_simulation_type(name)
    loaded_simulation = None if loaded_simulation_type == SimulationTypes.UNDEFINED else \
        load_simulation(name, loaded_simulation_type)

    match loaded_simulation_type:
        case SimulationTypes.UNDEFINED:
            simulation = model.Model(config_dict, setup, name=name, **model_kwargs)
        case SimulationTypes.MODEL_RUN:
            simulation = model.Model(loaded_simulation.config, loaded_simulation.setup,
                                     name=name)  # TODO: introduce model_kwarg
        case SimulationTypes.REFERENCES_SIMULATION:
            simulation = loaded_simulation
        case _:
            raise ValueError(
                "Loading save simulation type : Simulation type should be MODEL_RUN, REFERENCES_SIMULATION or UNDEFINED")

    #
    if save_type == SimulationTypes.MODEL_RUN or (
            save_type == SimulationTypes.REFERENCES_SIMULATION and loaded_simulation_type != SimulationTypes.REFERENCES_SIMULATION):
        save_simulation(
            model=simulation,
            base_name=name if loaded_simulation_type == SimulationTypes.UNDEFINED else f're_run_of_{name}',
            sim_type=save_type,
            parent_simulation=parent_simulation,
            user_notes=user_notes
        )

    return simulation


def load_simulation(name: str, simulation_type: SimulationTypes = SimulationTypes.REFERENCES_SIMULATION) -> Any:
    """
    Load a saved simulation

    Args:
        name: Name of the simulation
        simulation_type : Type of simulation loaded. If REFERENCES_SIMULATION, load the full simulation, if MODEL_RUN, load name, setud and config
    """

    if simulation_type == SimulationTypes.REFERENCES_SIMULATION:
        ref_dir = os.path.join(path_cfg.REFERENCES_SIMULATION_DIR, name)
        if os.path.exists(ref_dir):
            model_path = os.path.join(ref_dir, 'model.pkl')
            if not os.path.exists(model_path):
                raise FileNotFoundError(f"Reference model file not found: {model_path}")
            with open(model_path, 'rb') as f:
                return dill.load(f)

    sim_dir = os.path.join(path_cfg.MODEL_RUNS_DIR, name)
    if not os.path.exists(sim_dir):
        raise FileNotFoundError(f"Simulation '{name}' not found")

    with open(os.path.join(sim_dir, 'setup.pkl'), 'rb') as f:
        setup = pickle.load(f)

    with open(os.path.join(sim_dir, 'config.pkl'), 'rb') as f:
        config = pickle.load(f)

    sim_data = type('SimulationData', (), {
        'config': config,
        'setup': setup,
        'name': name
    })
    return sim_data


def get_saved_simulation_type(name: str) -> SimulationTypes:
    """
    Determine the type of a saved simulation by checking which directory it exists in.

    Args:
        name: The name of the simulation to check

    Returns:
        SimulationTypes: The type of the simulation (REFERENCES_SIMULATION, MODEL_RUN)
                        or UNDEFINED if not found
    """
    simulation_dirs = {
        SimulationTypes.REFERENCES_SIMULATION: path_cfg.REFERENCES_SIMULATION_DIR,
        SimulationTypes.MODEL_RUN: path_cfg.MODEL_RUNS_DIR
    }

    for sim_type, sim_directory in simulation_dirs.items():
        if os.path.exists(os.path.join(sim_directory, name)):
            return sim_type

    return SimulationTypes.UNDEFINED


def list_simulations(sim_type: Optional[str] = None) -> pd.DataFrame:
    """List all saved Simulations with their metadata"""
    log_df = _load_log()
    if sim_type:
        return log_df[log_df['runtime_info'].str.contains(sim_type)]
    return log_df


def get_simulation_lineage(name: str) -> List[str]:
    """Get the full lineage of a simulation (parent chain)"""
    log_df = _load_log()
    lineage = [name]
    current = name

    while True:
        parent = log_df[log_df['full_name'] == current]['parent_simulation'].iloc[0]
        if not parent:
            break
        lineage.append(parent)
        current = parent

    return lineage


def cleanup_old_simulations(days_threshold: int = 30):
    """Clean up old regular Simulations (not references) beyond threshold"""
    current_time = datetime.now()

    for sim in os.listdir(path_cfg.MODEL_RUNS_DIR):
        sim_path = os.path.join(path_cfg.MODEL_RUNS_DIR, sim)
        if not os.path.isdir(sim_path):
            continue

        try:
            timestamp = datetime.strptime(sim.split('_')[0], "%Y%m%d")
            if (current_time - timestamp).days > days_threshold:
                for file in os.listdir(sim_path):
                    if file.endswith('.pkl'):
                        os.remove(os.path.join(sim_path, file))
        except ValueError:
            continue


def parse_parameter_changes(changes: Dict[str, Union[float, List[float], np.ndarray, Any]]) -> Dict[str, Dict]:
    """
    Parse parameter changes from 'Component+parameter' format into nested dict format.

    Args:
        changes: Dict with keys in 'Component+parameter' format and values as:
                - single values (float, int)
                - iterables (list, np.ndarray, range, etc.)

    Returns:
        Nested dict suitable for deep_update:
        {'Component': {'parameters': {'parameter': value}}}
    """
    updates = {}
    for param_key, value in changes.items():
        if '+' not in param_key:
            raise ValueError(f"Parameter key must be in format 'Component+parameter': {param_key}")

        component, param = param_key.split('+')
        if component not in updates:
            updates[component] = {'parameters': {}}

        # Convert iterables to list for consistent handling
        if not isinstance(value, (float, int)) and hasattr(value, '__iter__'):
            if isinstance(value, np.ndarray):
                value = value.tolist()
            else:
                value = list(value)

        updates[component]['parameters'][param] = value

    return updates


def run_sensitivity(base_simulation: 'Model',
                    parameter_changes: Dict[str, Union[float, List[float], np.ndarray, Any]],
                    name: Optional[str] = None,
                    user_notes: Optional[str] = None,
                    save: bool = True,
                    setup: Optional[Any] = None,
                    **kwargs) -> Union['Model', List['Model']]:
    """
    Run sensitivity analysis based on an existing simulation.
    Handles both single parameter changes and batch sensitivity analysis.

    Args:
        base_simulation: Existing Model instance to use as baseline
        parameter_changes: Dict of parameter changes where values can be either:
                         - single value: {'Component+parameter': value}
                         - list of values: {'Component+parameter': [value1, value2, ...]}
        name: Base name for the simulation(s)
        user_notes: Notes about the sensitivity run(s)
        save: Whether to save the simulation(s)
        setup: Optional Setup instance. If None, uses base_simulation.setup
        **kwargs: Additional arguments for Model constructor

    Returns:
        Single Model instance or list of Model instances depending on input
    """
    # Check if any parameter has multiple values
    has_multiple = any(isinstance(v, (list, np.ndarray)) for v in parameter_changes.values())

    # Use provided setup or fall back to base simulation setup
    modified_setup = setup if setup is not None else base_simulation.setup

    if not has_multiple:
        # Single sensitivity run
        updates = parse_parameter_changes(parameter_changes)
        new_config = fns.deep_update(base_simulation.config.copy(), updates)

        # Generate a descriptive name if none provided
        run_name = name
        if run_name is None:
            # Create a descriptive name based on the parameters and their values
            param_descriptions = []
            for param_key, value in parameter_changes.items():
                param_descriptions.append(f"{param_key.replace('+', '_')}={value:.3g}")
            run_name = f"sensitivity_{'_'.join(param_descriptions)}"

        model_instance = run_or_load_simulation(
            config_dict=new_config,
            setup=modified_setup,
            name=run_name,
            parent_simulation=base_simulation.name if hasattr(base_simulation, 'name') else None,
            user_notes=user_notes,
            save_type=SimulationTypes.UNDEFINED if not save else SimulationTypes.MODEL_RUN,
            **kwargs
        )

        model_instance.name = run_name

        return model_instance

    # Batch sensitivity run
    # Find the parameter with multiple values
    multi_param = None
    multi_values = None
    for param, values in parameter_changes.items():
        if isinstance(values, (list, np.ndarray)):
            if multi_param is not None:
                raise ValueError("Multiple parameters with lists of values not supported")
            multi_param = param
            multi_values = values

    # Generate base name if not provided
    if name is None:
        name = f"sensitivity_{multi_param.replace('+', '_')}"

    # Run Simulations for each value
    results = []
    for i, value in enumerate(multi_values):
        # Create parameter set for this run
        current_changes = parameter_changes.copy()
        current_changes[multi_param] = value

        # Create names for this run
        sim_name = f"{name}_{value:.3g}"
        sim_notes = f"{user_notes} - {multi_param}={value:.3g}" if user_notes else f"{multi_param}={value:.3g}"

        # Run simulation
        results.append(run_sensitivity(
            base_simulation,
            current_changes,
            name=sim_name,
            user_notes=sim_notes,
            setup=setup,
            **kwargs
        ))

    return results


def compare_simulations(models: Union[List, Dict[str, Any]],
                        variables: Optional[List[str]] = None,
                        plot: bool = True) -> Dict[str, Dict[str, float]]:
    """Compare multiple Simulations"""

    if not isinstance(models, dict):
        # Flatten list and create dictionary
        flattened = fns.flatten_simulation_list(models)
        models = {f"Model_{i + 1}": model for i, model in enumerate(flattened)}

    if variables is None:
        variables = sorted(set.intersection(*[set(m.df.columns) for m in models.values()]))

    differences = {}
    if plot:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(variables), 1,
                                 figsize=(10, 4 * len(variables)))
        if len(variables) == 1:
            axes = [axes]

    for i, var in enumerate(variables):
        series = {name: model.df[var] for name, model in models.items()}

        var_diffs = {}
        for name1 in models:
            for name2 in models:
                if name1 < name2:
                    diff = np.abs(series[name1] - series[name2])
                    var_diffs[f"{name1}_vs_{name2}"] = {
                        'max_diff': diff.max(),
                        'mean_diff': diff.mean(),
                        'std_diff': diff.std()
                    }

        differences[var] = var_diffs

        if plot:
            ax = axes[i]
            for name, s in series.items():
                ax.plot(s.index, s.values, label=name)
            ax.set_title(f'{var}')
            ax.legend()

    if plot:
        plt.tight_layout()

    return differences
