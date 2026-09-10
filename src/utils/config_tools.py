"""Tooling for the model CONFIGURATION DICTIONARIES.

Everything that reads, writes, compares or transforms a configuration dict lives here:
merging levels (`deep_update`), applying optimised parameters (`update_config`), dropping
a currency or a set of components, and the JSON serialisation used by the simulation log.

The model's own helpers -- limitation functions, expression evaluation, coupled-term
summation -- are in `functions.py`.
"""
import json
import os
import pickle
import re
from typing import Any, Dict, List, TypeVar, Union

import numpy as np
import pandas as pd
from deepdiff import DeepDiff

from src.config_model import varinfos

KeyType = TypeVar('KeyType')

def serialize_for_json(obj, preserve_full_arrays: bool = False, circular_ref_check: bool = True) -> Any:
    """
    Serialize Python objects into JSON-compatible format, with special handling for scientific objects.

    Handles:
    - NumPy arrays and scalars
    - Pandas DataFrames, Series, timestamps
    - Python classes and types
    - Nested dictionaries and lists
    - Custom objects with __dict__

    Args:
        obj: Object to serialize
        preserve_full_arrays: If False, truncates long arrays for readability
        circular_ref_check: Enable circular reference protection for complex objects

    Returns:
        JSON-serializable version of the object
    """

    def _truncate_list(lst, max_items=4):
        """Truncate long lists/arrays for readability"""
        if not preserve_full_arrays and len(lst) > max_items:
            return lst[:2] + ["..."] + lst[-2:]
        return lst

    seen = {} if circular_ref_check else None

    def _should_track_reference(obj):
        """Determine if object should be tracked for circular references"""
        return not isinstance(obj, (int, float, bool, str, type(None)))

    def _serialize_recursive(obj, _seen=None):
        # Circular reference check only for complex objects
        if _seen is not None and _should_track_reference(obj):
            obj_id = id(obj)
            if obj_id in _seen:
                if hasattr(obj, '__dict__'):
                    return f"<circular-reference to {obj.__class__.__name__}>"
                return "<circular-reference>"
            _seen[obj_id] = True

        # Handle basic types
        if isinstance(obj, (int, float, bool, str, type(None))):
            return obj

        # Handle numpy types
        if isinstance(obj, (np.ndarray, np.generic)):
            if isinstance(obj, (np.float32, np.float64, np.int32, np.int64)):
                return obj.item()  # Convert numpy scalars to Python scalars
            return _truncate_list(obj.tolist())

        # Handle pandas types
        if isinstance(obj, pd.DataFrame):
            data_dict = obj.to_dict(orient='list')
            return {k: _truncate_list(v) for k, v in data_dict.items()}
        if isinstance(obj, pd.Series):
            return _truncate_list(obj.tolist())
        if isinstance(obj, pd.Timestamp):
            return obj.isoformat()
        if isinstance(obj, pd.DatetimeIndex):
            return _truncate_list(obj.astype(str).tolist())

        # Handle Python classes and types
        if isinstance(obj, type):
            return {
                "__type__": "class",
                "module": obj.__module__,
                "name": obj.__name__
            }

        # Handle dictionaries
        if isinstance(obj, dict):
            return {
                str(k): _serialize_recursive(v, _seen)
                for k, v in obj.items()
            }

        # Handle lists and tuples
        if isinstance(obj, (list, tuple)):
            return [_serialize_recursive(item, _seen) for item in obj]

        # Handle custom objects with __dict__
        if hasattr(obj, '__dict__'):
            class_info = {
                "__type__": "object",
                "class": obj.__class__.__name__,
                "module": obj.__class__.__module__,
            }
            try:
                # Try to get a string representation if available
                if hasattr(obj, '__str__'):
                    class_info["string_repr"] = str(obj)
                elif hasattr(obj, '__repr__'):
                    class_info["string_repr"] = repr(obj)

                # Add serialized attributes
                class_info["attributes"] = _serialize_recursive(obj.__dict__, _seen)
            except Exception as e:
                class_info["serialization_error"] = str(e)
            return class_info

        # Fall back to string representation
        try:
            return str(obj)
        except Exception:
            return f"<non-serializable: {type(obj).__name__}>"

    return _serialize_recursive(obj, seen)

def save_human_readable_config(config: Dict, filepath: str) -> None:
    """
    Save a human-readable version of the configuration dictionary.

    Args:
        config: Configuration dictionary
        filepath: Output JSON file path
    """
    serialized_config = serialize_for_json(config)

    with open(filepath, 'w') as f:
        json.dump(serialized_config, f, indent=2)

def print_dict(dict):
    print(json.dumps(serialize_for_json(dict), indent=4))

def write_dict_to_file(dict, fname,
                       fdir='optim_result/'):
    serialized_dict = serialize_for_json(dict)
    fullname = os.path.join(fdir, fname)
    with open(fullname+'.json', 'w') as file:
        json.dump(serialized_dict, file, indent=4)
    print('Dict written to ', fullname + '.json')

    with open(fullname+'.pkl', 'wb') as file:
        pickle.dump(dict, file)
    print('Dict written to ', fullname + '.pkl')

def compare_dicts(dict1: Dict, dict2: Dict,
                  format_output: bool = True,
                  config_mode: bool = False,
                  print_result: bool = False) -> Union[str, Dict]:
    """
    Compare two dictionaries and return formatted differences.

    Args:
        dict1: First dictionary
        dict2: Second dictionary
        format_output: If True returns formatted string, if False returns DeepDiff object
        config_mode: If True, uses config-specific formatting (for model configurations)
        print_result: If True, prints the result to console

    Returns:
        Formatted differences string, DeepDiff object, or prints result
    """
    if dict1 is None or dict2 is None:
        result = ""
        if print_result:
            print("One or both dictionaries are None")
        return result

    # Sanitize configs if in config mode
    if config_mode:
        dict1 = _sanitize_config_for_comparison(dict1)
        dict2 = _sanitize_config_for_comparison(dict2)

    # Use DeepDiff with appropriate settings
    diff = DeepDiff(dict1, dict2,
                    ignore_order=True,
                    report_repetition=False,
                    verbose_level=2)

    if not diff:
        result = "No differences found. The dictionaries are identical." if not config_mode else "Idem"
        if print_result:
            print(result)
        return result

    if not format_output:
        if print_result:
            print("Differences found:")
            print(diff)
        return diff

    if config_mode:
        result = _format_config_diff(diff, dict1, dict2)
    else:
        result = _format_general_diff(diff)

    if print_result:
        print("Differences found:")
        print(result)

    return result

def _sanitize_config_for_comparison(config):
    """Remove problematic keys from config for comparison"""
    result = {}
    for k, v in config.items():
        if k == 'formulation':
            result[k] = v
            continue

        result[k] = {}
        for subk, subv in v.items():
            if subk == 'class':
                continue
            elif subk == 'parameters' and isinstance(subv, dict):
                result[k][subk] = {}
                for paramk, paramv in subv.items():
                    if hasattr(paramv, 'item') and callable(getattr(paramv, 'item')):
                        result[k][subk][paramk] = paramv.item()
                    else:
                        result[k][subk][paramk] = paramv
            else:
                result[k][subk] = subv
    return result

def _format_config_diff(diff, config1, config2):
    """Format differences for model configurations"""
    diff_by_component = {}

    for change_type, changes in diff.items():
        if change_type == 'values_changed':
            for path, change in changes.items():
                component, param = _extract_component_param_from_path(path)
                if component and param:
                    old_val = _format_value_for_display(change['old_value'])
                    new_val = _format_value_for_display(change['new_value'])
                    _add_change_to_component(diff_by_component, component, f"{param}: {old_val} → {new_val}")

        elif change_type == 'dictionary_item_added':
            for path in changes:
                component, param = _extract_component_param_from_path(path)
                if component and param:
                    try:
                        value = _get_value_from_diff_path(config2, path)
                        value_str = _format_value_for_display(value)
                        _add_change_to_component(diff_by_component, component, f"Added {param}={value_str}")
                    except:
                        _add_change_to_component(diff_by_component, component, f"Added {param}")

        elif change_type == 'dictionary_item_removed':
            for path in changes:
                component, param = _extract_component_param_from_path(path)
                if component and param:
                    try:
                        value = _get_value_from_diff_path(config1, path)
                        value_str = _format_value_for_display(value)
                        _add_change_to_component(diff_by_component, component, f"Removed {param}={value_str}")
                    except:
                        _add_change_to_component(diff_by_component, component, f"Removed {param}")

        elif change_type == 'type_changes':
            for path, change in changes.items():
                component, param = _extract_component_param_from_path(path)
                if component and param:
                    old_val = _format_value_for_display(change['old_value'])
                    new_val = _format_value_for_display(change['new_value'])
                    _add_change_to_component(diff_by_component, component, f"{param} type: {old_val} → {new_val}")

    # Construct the output string
    result = []
    for component, changes in diff_by_component.items():
        result.append(f"Changes in {component}:\n  " + "\n  ".join(changes))

    return "\n".join(result) if result else "Idem"

def _format_general_diff(diff):
    """Format differences for general dictionaries"""
    result = []
    for change_type, changes in diff.items():
        result.append(f"{change_type}:")
        if isinstance(changes, dict):
            for path, change in changes.items():
                if isinstance(change, dict) and 'old_value' in change:
                    result.append(f"  {path}: {change['old_value']} → {change['new_value']}")
                else:
                    result.append(f"  {path}: {change}")
        else:
            result.append(f"  {changes}")
    return "\n".join(result)

def _extract_component_param_from_path(path):
    """Extract component and parameter from a DeepDiff path"""
    pattern = r"root\['([^']+)'\](?:\['parameters'\])?\['([^']+)'\]"
    match = re.search(pattern, path)
    if match:
        return match.group(1), match.group(2)
    return None, None

def _format_value_for_display(value):
    """Format a value for display"""
    if isinstance(value, float):
        if abs(value) < 0.001 or abs(value) > 1000:
            return f"{value:.3e}"
        return f"{value:.4g}"
    return str(value)

def _add_change_to_component(diff_dict, component, change):
    """Add a change to the component dictionary"""
    if component not in diff_dict:
        diff_dict[component] = []
    diff_dict[component].append(change)

def _get_value_from_diff_path(config, path):
    """Get a value from a config using a DeepDiff path"""
    parts = re.findall(r"\['([^']+)'\]", path)
    if not parts:
        return None

    value = config
    for part in parts[1:]:  # Skip 'root'
        value = value[part]
    return value

def print_config_diff(config1, config2, config_mode: bool = True) -> None:
    """Print what differs between two configuration dicts, grouped by component.

    The diagnostic to reach for when two runs diverge and the reason is not obvious.
    Both arguments are DICTS, not simulations: pass `sim.config` for a run, a
    configuration dict built by hand, or `setup.to_dict()` for a Setup -- in the latter
    case use config_mode=False, the component/parameter grouping does not apply.
    """
    compare_dicts(config1, config2, format_output=True,
                  config_mode=config_mode, print_result=True)

def update_config(dconf: Dict, param_dict: Dict[str, float]) -> Dict:
    """
    Update configuration with parameter values from a dictionary.

    Args:
        dconf: Base configuration dictionary
        param_dict: Dictionary mapping 'Component+parameter' to value
                   Example: {'Phy+mu_max': 1.37, 'Macroflocs+resuspension_rate': 75941.93}
                   Special prefix 'BGC+' applies parameters to all BGC components

    Returns:
        Updated configuration dictionary with parameters applied

    Raises:
        ValueError: If parameter name format is invalid (missing '+')

    Example:
        >>> new_config = update_config(base_config, {
        ...     'Phy+mu_max': 1.37,
        ...     'Phy+mortrate': 0.01,
        ...     'Macroflocs+resuspension_rate': 75941.93,
        ...     'BGC+resusp_ewma_alpha': 0.02  # Applied to all BGC components
        ... })
    """
    # BGC components that share common parameters in bgc_only mode
    BGC_COMPONENTS = ['DOCS', 'DOCL', 'DetL', 'DetS', 'BacA']
    # Mapping for BGC parameters that need renaming when applied to components
    BGC_PARAM_MAPPING = {
        'resusp_ewma_alpha': 'prescribed_resusp_ewma_alpha',
        'vertical_coupling_alpha': 'prescribed_resusp_ewma_alpha',  # backward compat (to remove when obsolete)
    }

    # Build nested update dictionary
    updates = {}
    for param_key, value in param_dict.items():
        if '+' not in param_key:
            raise ValueError(
                f"Invalid parameter name format: '{param_key}' "
                f"(expected 'Component+parameter', e.g., 'Phy+mu_max')"
            )

        component, param = param_key.split('+', 1)  # maxsplit=1 to handle '+' in param names

        # Special handling for BGC+ prefix: broadcast to all BGC components
        components = BGC_COMPONENTS if component == 'BGC' else [component]
        target_param = BGC_PARAM_MAPPING.get(param, param) if component == 'BGC' else param

        for comp in components:
            updates.setdefault(comp, {'parameters': {}})['parameters'][target_param] = value

    return deep_update(dconf, updates)

def remove_currency(dconf: Dict, currency: str = 'P') -> Dict:
    """
    Remove a specific currency (P, Si, N) from a biogeochemical model configuration.

    Removes all initialization values, couplings, aggregates, parameters, and diagnostics
    related to the specified currency. The model will treat the currency as absent (None).

    Args:
        dconf: Model configuration dictionary
        currency: Currency to remove - 'P', 'Si', or 'N' (case-insensitive)

    Returns:
        Modified configuration dictionary with currency removed

    Example:
        >>> # Remove phosphorus from the model
        >>> dconf_no_P = fns.remove_currency(dconf, 'P')
        >>>
        >>> # Remove silicon
        >>> dconf_no_Si = fns.remove_currency(dconf, 'Si')
    """
    import copy
    currency = currency.upper()

    if currency not in varinfos.CURRENCY_MAP:
        raise ValueError(f"Currency '{currency}' not supported. Use 'P', 'Si', or 'N'.")

    config = varinfos.CURRENCY_MAP[currency]
    result = copy.deepcopy(dconf)

    # Helper to ensure list
    def as_list(x): return x if isinstance(x, list) else [x]

    for comp_name, comp_config in result.items():
        if comp_name == 'formulation':
            continue

        # Remove from initialization
        if 'initialization' in comp_config and config['init_key'] in comp_config['initialization']:
            del comp_config['initialization'][config['init_key']]

        # Remove from coupling
        if 'coupling' in comp_config:
            for coup_key in as_list(config['coupling_key']):
                comp_config['coupling'].pop(coup_key, None)

            # Remove from coupled_remin_products (for detritus/DOM components)
            if 'coupled_remin_products' in comp_config['coupling']:
                remin_list = comp_config['coupling']['coupled_remin_products']
                for remin_key in as_list(config['remin_key']):
                    if remin_key in remin_list:
                        comp_config['coupling']['coupled_remin_products'] = [
                            x for x in remin_list if x != remin_key
                        ]

        # Remove from aggregates
        if 'aggregate' in comp_config:
            for agg_key in config['aggregate_keys']:
                comp_config['aggregate'].pop(agg_key, None)

        # Remove from parameters
        if 'parameters' in comp_config:
            for param_key in config['param_keys']:
                comp_config['parameters'].pop(param_key, None)

        # Remove from diagnostics
        if 'diagnostics' in comp_config:
            comp_config['diagnostics'] = [
                diag for diag in comp_config['diagnostics']
                if not any(pattern in diag for pattern in config['diag_patterns'])
            ]

    # Remove nutrient component(s) entirely
    for nut_comp in as_list(config['nutrient_component']):
        result.pop(nut_comp, None)

    return result

def remove_components(dconf: Dict, components_map: Dict) -> Dict:
    """
    Remove specific components from a biogeochemical model configuration.

    Removes all components specified in the map and cleans up external couplings
    where other components reference the removed components.

    Args:
        dconf: Model configuration dictionary
        components_map: Dictionary mapping component names to their coupling info
                       Format: {'ComponentName': {'coupled_in': {...}}}

    Returns:
        Modified configuration dictionary with components removed

    Example:
        >>> # Remove all Flocs components
        >>> dconf_no_flocs = fns.remove_components(dconf, varinfos.FLOCS_COMPONENTS_MAP)
        >>>
        >>> # Remove only Macroflocs
        >>> dconf_no_macro = fns.remove_components(dconf, {'Macroflocs': varinfos.FLOCS_COMPONENTS_MAP['Macroflocs']})
    """
    import copy
    result = copy.deepcopy(dconf)

    for comp_name, comp_info in components_map.items():
        # Remove the component itself
        result.pop(comp_name, None)

        # Clean up external couplings (where other components reference this one)
        if 'coupled_in' in comp_info:
            for external_comp, coupling_info in comp_info['coupled_in'].items():
                if external_comp in result and 'coupling' in result[external_comp]:
                    # Support both single string and list of coupling keys
                    coupling_keys = coupling_info['coupling']
                    if isinstance(coupling_keys, str):
                        coupling_keys = [coupling_keys]

                    for coupling_key in coupling_keys:
                        if coupling_key in result[external_comp]['coupling']:
                            coupled_value = result[external_comp]['coupling'][coupling_key]

                            # Handle list case (e.g., coupled_SPM: ['Microflocs', 'Micro_in_Macro'])
                            if isinstance(coupled_value, list):
                                result[external_comp]['coupling'][coupling_key] = [
                                    x for x in coupled_value if x != comp_name
                                ]
                                # Remove key if list becomes empty
                                if not result[external_comp]['coupling'][coupling_key]:
                                    del result[external_comp]['coupling'][coupling_key]

                            # Handle scalar case (e.g., coupled_aggregate: 'Macroflocs')
                            elif coupled_value == comp_name:
                                del result[external_comp]['coupling'][coupling_key]

    return result

def deep_update(mapping: Dict[KeyType, Any], *updating_mappings: Dict[KeyType, Any],
                overwrite_keys: list = None, merge_lists: list = None) -> Dict[KeyType, Any]:
    """
    Deep update dictionaries with selective overwrite behavior.

    Args:
        mapping: Base dictionary to update
        *updating_mappings: One or more dictionaries to merge in
        overwrite_keys: List of keys where dict values should be overwritten instead of merged
        merge_lists: List of keys where list values should be merged instead of overwritten.
                    Defaults to ['diagnostics'] for common use case of merging diagnostic lists.
    """
    updated_mapping = mapping.copy()
    overwrite_keys = overwrite_keys or []
    merge_lists = merge_lists or ['diagnostics']

    for updating_mapping in updating_mappings:
        for k, v in updating_mapping.items():
            if (k in updated_mapping and
                    isinstance(updated_mapping[k], dict) and
                    isinstance(v, dict) and
                    k not in overwrite_keys):
                # Recursive merge (current behavior)
                updated_mapping[k] = deep_update(updated_mapping[k], v, 
                                                overwrite_keys=overwrite_keys, 
                                                merge_lists=merge_lists)
            elif (k in updated_mapping and
                  k in merge_lists and
                  isinstance(updated_mapping[k], list) and
                  isinstance(v, list)):
                # Merge lists by extending and removing duplicates while preserving order
                merged_list = updated_mapping[k].copy()
                for item in v:
                    if item not in merged_list:
                        merged_list.append(item)
                updated_mapping[k] = merged_list
            else:
                # Direct assignment (overwrite or non-dict)
                updated_mapping[k] = v
    return updated_mapping
