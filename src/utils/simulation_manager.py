"""
Simulation management utilities for BGC model.
Handles saving, loading, and tracking of simulations.
"""
import os
import dill
import pickle
import json
import pandas as pd
import numpy as np
import re
from datetime import datetime
from typing import Optional, List, Dict, Any, Union
from deepdiff import DeepDiff

from src.utils import functions as fns
# Constants
BASE_DIR = '../simulations'
SIMULATIONS_DIR = os.path.join(BASE_DIR, 'Simulations')
REFERENCES_DIR = os.path.join(BASE_DIR, 'References')
LOG_FILE = os.path.join(BASE_DIR, 'Simulations_log.xlsx')

# Create necessary directories
for directory in [BASE_DIR, SIMULATIONS_DIR, REFERENCES_DIR]:
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
    'date',              # YYYY-MM-DD
    'short_name',        # Name without timestamp
    'user_notes',        # User provided notes
    'parent_simulation', # Parent simulation name if any
    'configuration_diff',# Differences from parent config
    'setup_diff',        # Differences from parent setup
    'status',           # Simulation status
    'setup_start_date', # Simulation start date
    'setup_end_date',   # Simulation end date
] + [
    # Setup fields - add all tracked fields
    f'setup_{field}' for field in SETUP_TRACKED_FIELDS
] + [
    # Configuration and runtime fields
    'config_formulation',
    'config_components',  # Detailed component list as JSON
    'variables',         # Model variables
    'runtime_info',      # Additional runtime information
    'full_name'         # Complete name with timestamp
]

def _load_log() -> pd.DataFrame:
    """Load or create simulation log with organized columns"""
    if os.path.exists(LOG_FILE):
        df = pd.read_excel(LOG_FILE)
        # Force column order and drop any extra columns
        return df.reindex(columns=LOG_COLUMNS)

    return pd.DataFrame(columns=LOG_COLUMNS)


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
    if config1 is None or config2 is None:
        return ""

    # Don't compare 'class' keys as they cause issues with DeepDiff
    sanitized_config1 = _sanitize_config(config1)
    sanitized_config2 = _sanitize_config(config2)

    # Use DeepDiff with appropriate settings
    diff = DeepDiff(sanitized_config1, sanitized_config2,
                    ignore_order=True,
                    report_repetition=False,
                    verbose_level=2)

    if not diff:
        return "Idem"

    # Process the differences by component
    diff_by_component = {}

    # Process parameter changes
    for change_type, changes in diff.items():
        if change_type == 'values_changed':
            for path, change in changes.items():
                component, param = _extract_component_param(path)
                if component and param:
                    old_val = _format_value(change['old_value'])
                    new_val = _format_value(change['new_value'])
                    _add_change(diff_by_component, component, f"{param}: {old_val} → {new_val}")

        elif change_type == 'dictionary_item_added':
            for path in changes:
                component, param = _extract_component_param(path)
                if component and param:
                    # Try to get the value
                    try:
                        value = _get_value_from_path(config2, path)
                        value_str = _format_value(value)
                        _add_change(diff_by_component, component, f"Added {param}={value_str}")
                    except:
                        _add_change(diff_by_component, component, f"Added {param}")

        elif change_type == 'dictionary_item_removed':
            for path in changes:
                component, param = _extract_component_param(path)
                if component and param:
                    # Try to get the value
                    try:
                        value = _get_value_from_path(config1, path)
                        value_str = _format_value(value)
                        _add_change(diff_by_component, component, f"Removed {param}={value_str}")
                    except:
                        _add_change(diff_by_component, component, f"Removed {param}")

        elif change_type == 'type_changes':
            for path, change in changes.items():
                component, param = _extract_component_param(path)
                if component and param:
                    old_val = _format_value(change['old_value'])
                    new_val = _format_value(change['new_value'])
                    _add_change(diff_by_component, component, f"{param} type: {old_val} → {new_val}")

    # Construct the output string
    result = []
    for component, changes in diff_by_component.items():
        result.append(f"Changes in {component}:\n  " + "\n  ".join(changes))

    # Manual parameter comparison if DeepDiff misses something
    if not result:
        result = _manual_param_comparison(config1, config2)

    return "\n".join(result) if result else "Idem"


def _sanitize_config(config):
    """Remove problematic keys from the config"""
    result = {}
    for k, v in config.items():
        if k == 'formulation':
            result[k] = v
            continue

        result[k] = {}
        for subk, subv in v.items():
            if subk == 'class':
                # Skip class objects
                continue
            elif subk == 'parameters' and isinstance(subv, dict):
                # Deep copy parameters
                result[k][subk] = {}
                for paramk, paramv in subv.items():
                    # Convert numpy types to Python native types
                    if hasattr(paramv, 'item') and callable(getattr(paramv, 'item')):
                        result[k][subk][paramk] = paramv.item()
                    else:
                        result[k][subk][paramk] = paramv
            else:
                result[k][subk] = subv
    return result


def _extract_component_param(path):
    """Extract component and parameter from a DeepDiff path"""
    # Handle bracket notation: root['component']['parameters']['param']
    pattern = r"root\['([^']+)'\](?:\['parameters'\])?\['([^']+)'\]"
    match = re.search(pattern, path)
    if match:
        return match.group(1), match.group(2)
    return None, None


def _format_value(value):
    """Format a value for display"""
    if isinstance(value, float):
        if abs(value) < 0.001 or abs(value) > 1000:
            return f"{value:.3e}"
        return f"{value:.4g}"
    return str(value)


def _add_change(diff_dict, component, change):
    """Add a change to the component dictionary"""
    if component not in diff_dict:
        diff_dict[component] = []
    diff_dict[component].append(change)


def _get_value_from_path(config, path):
    """Get a value from a config using a DeepDiff path"""
    # Convert path like "root['component']['parameters']['param']" to value
    parts = re.findall(r"\['([^']+)'\]", path)
    if not parts or parts[0] != 'root':
        return None

    value = config
    for part in parts[1:]:
        value = value[part]
    return value


def _manual_param_comparison(config1, config2):
    """Manually compare parameters if DeepDiff didn't find anything"""
    result = []

    # Compare parameters for each component
    for component in set(config1.keys()) | set(config2.keys()):
        component_changes = []

        # Skip if component doesn't exist in both configs
        if component not in config1 or component not in config2:
            continue

        # Get parameters dictionaries
        params1 = config1[component].get('parameters', {})
        params2 = config2[component].get('parameters', {})

        # Compare parameters
        for param in set(params1.keys()) | set(params2.keys()):
            if param not in params1:
                value_str = _format_value(params2[param])
                component_changes.append(f"Added {param}={value_str}")
            elif param not in params2:
                value_str = _format_value(params1[param])
                component_changes.append(f"Removed {param}={value_str}")
            elif params1[param] != params2[param]:
                # Convert numpy values for comparison
                val1 = params1[param].item() if hasattr(params1[param], 'item') else params1[param]
                val2 = params2[param].item() if hasattr(params2[param], 'item') else params2[param]

                if val1 != val2:
                    old_val = _format_value(val1)
                    new_val = _format_value(val2)
                    component_changes.append(f"{param}: {old_val} → {new_val}")

        if component_changes:
            result.append(f"Changes in {component}:\n  " + "\n  ".join(component_changes))

    return result


def save_simulation(
        model,
        base_name: Optional[str] = None,
        sim_type: str = 'simulation',
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
        REFERENCES_DIR if sim_type == 'reference' else SIMULATIONS_DIR,
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
        'runtime_info': f"Type: {sim_type}",
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

    if save_full or sim_type == 'reference':
        key_results = {
            'timestamps': model.df.index.values,
            'variables': {col: model.df[col].values for col in model.df.columns}
        }
        np.save(os.path.join(sim_dir, 'results.npy'), key_results)
        with open(os.path.join(sim_dir, 'model.pkl'), 'wb') as f:
            dill.dump(model, f)

    # Update log with strict column order
    log_df = _load_log()
    new_row = pd.DataFrame([{col: log_entry.get(col, '') for col in LOG_COLUMNS}])
    log_df = pd.concat([log_df, new_row], ignore_index=True)
    log_df.to_excel(LOG_FILE, index=False)

    return full_name

def run_or_load_simulation(config_dict: Dict,
                   setup: Any,
                   name: Optional[str] = None,
                   parent_simulation: Optional[str] = None,
                   user_notes: Optional[str] = None,
                   save: bool = True,
                   **model_kwargs) -> 'Model':
    """
    Get a simulation by either loading it or creating it.

    Args:
        config_dict: Model configuration dictionary
        setup: Model setup instance
        name: Optional name for the simulation
        parent_simulation: Optional parent simulation name
        user_notes: Optional notes about the simulation
        save: Whether to save the simulation
        **model_kwargs: Additional arguments for Model constructor
    """
    from src.core import model

    if name and simulation_exists(name):
        return load_simulation(name)

    simulation = model.Model(config_dict, setup, name=name, **model_kwargs)

    if save:
        save_simulation(
            model=simulation,
            base_name=name,
            parent_simulation=parent_simulation,
            user_notes=user_notes
        )

    return simulation


def load_simulation(name: str, load_full: bool = False) -> Any:
    """
    Load a saved simulation

    Args:
        name: Name of the simulation
        load_full: Whether to load full pickle files
    """
    # First check if it's a reference simulation
    ref_dir = os.path.join(REFERENCES_DIR, name)
    if os.path.exists(ref_dir):
        model_path = os.path.join(ref_dir, 'model.pkl')
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Reference model file not found: {model_path}")
        with open(model_path, 'rb') as f:
            return dill.load(f)

    # Regular simulation loading logic
    sim_dir = os.path.join(SIMULATIONS_DIR, name)
    if not os.path.exists(sim_dir):
        raise FileNotFoundError(f"Simulation '{name}' not found")

    if load_full and os.path.exists(os.path.join(sim_dir, 'model.pkl')):
        with open(os.path.join(sim_dir, 'model.pkl'), 'rb') as f:
            return dill.load(f)

    # Load minimum required components
    setup = None
    setup_path = os.path.join(sim_dir, 'setup.pkl')
    if os.path.exists(setup_path):
        try:
            with open(setup_path, 'rb') as f:
                setup = dill.load(f)
        except Exception as e:
            print(f"Warning: Could not load setup from {setup_path}: {e}")

    results = np.load(os.path.join(sim_dir, 'results.npy'), allow_pickle=True).item()
    with open(os.path.join(sim_dir, 'config.pkl'), 'rb') as f:
        config = pickle.load(f)

    # Create DataFrame from results
    df = pd.DataFrame(
        {col: results['variables'][col] for col in results['variables']},
        index=results['timestamps']
    )

    # Return object with minimum required attributes
    sim_data = type('SimulationData', (), {
        'config': config,
        'setup': setup,
        'df': df,
        'name': name
    })
    return sim_data


def simulation_exists(name: str) -> bool:
    """Check if a simulation exists"""
    for sim_type in ['simulation', 'reference']:
        if os.path.exists(os.path.join(SIMULATIONS_DIR if sim_type == 'simulation' else REFERENCES_DIR, name)):
            return True
    return False


def list_simulations(sim_type: Optional[str] = None) -> pd.DataFrame:
    """List all saved simulations with their metadata"""
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
    """Clean up old regular simulations (not references) beyond threshold"""
    current_time = datetime.now()

    for sim in os.listdir(SIMULATIONS_DIR):
        sim_path = os.path.join(SIMULATIONS_DIR, sim)
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
                    setup_changes: Optional[Dict] = None,
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
        setup_changes: Optional dict of setup parameters to modify
        **kwargs: Additional arguments for Model constructor

    Returns:
        Single Model instance or list of Model instances depending on input
    """
    # Check if any parameter has multiple values
    has_multiple = any(isinstance(v, (list, np.ndarray)) for v in parameter_changes.values())

    # Create modified setup if needed
    modified_setup = base_simulation.setup

    if setup_changes:
        # Get only the initialization parameters from the original setup
        setup_params = {
            'tmin': base_simulation.setup.tmin,
            'tmax': base_simulation.setup.tmax,
            'dt': base_simulation.setup.dt,
            'dt2': base_simulation.setup.dt2,
            'dt2_s_to_d_ratio': base_simulation.setup.dt2_s_to_d_ratio,
            'start_date': base_simulation.setup.start_date,
            'PARfromfile': base_simulation.setup.PARfromfile,
            'I': base_simulation.setup.I,
            'light_prop': base_simulation.setup.light_prop,
            'lightfirst': base_simulation.setup.lightfirst,
            'T': base_simulation.setup.T - 273.15,  # Convert back to Celsius for initialization
            'k_att': base_simulation.setup.k_att,
            'z': base_simulation.setup.z,
            'pCO2': base_simulation.setup.pCO2,
            'g_shear_rate': base_simulation.setup.g_shear_rate.iloc[0,0] / (1 + base_simulation.setup.gshearfact) if
                hasattr(base_simulation.setup, 'gshearfact') else
                base_simulation.setup.g_shear_rate.iloc[0,0],
            'vary_g_shear': base_simulation.setup.vary_g_shear,
            'gshearfact': base_simulation.setup.gshearfact,
            'gshearper': base_simulation.setup.gshearper,
            'kb': base_simulation.setup.kb,
            'verbose': base_simulation.setup.verbose
        }
        # Update with requested changes
        setup_params.update(setup_changes)
        modified_setup = type(base_simulation.setup)(**setup_params)

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

        model_instance =  run_or_load_simulation(
            config_dict=new_config,
            setup=modified_setup,
            name=run_name,
            parent_simulation=base_simulation.name if hasattr(base_simulation, 'name') else None,
            user_notes=user_notes,
            save=save,
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

    # Run simulations for each value
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
            save=save,
            **kwargs
        ))

    return results

def compare_simulations(models: Union[List, Dict[str, Any]],
                        variables: Optional[List[str]] = None,
                        plot: bool = True) -> Dict[str, Dict[str, float]]:
    """Compare multiple simulations"""

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