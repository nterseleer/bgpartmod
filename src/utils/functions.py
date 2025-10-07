import ast
import operator as op
import os
import numpy as np
import pandas as pd
import json
import pickle
from deepdiff import DeepDiff
from typing import Any, Dict, List, Union

from core import phys
from config_model import phys_setup
from src.config_model import varinfos


def flatten_simulation_list(sims: Union[List, Any]) -> Union[List, Dict[str, Any]]:
    """
    Flatten potentially nested lists of Simulations while preserving dictionary structure if needed.

    Args:
        sims: Single simulation, list of Simulations, nested lists, or dictionary

    Returns:
        Flattened list of Simulations or dictionary mapping names to Simulations
    """
    # Handle dictionary input
    if isinstance(sims, dict):
        return sims

    # Handle non-list input
    if not isinstance(sims, list):
        return [sims]

    # Flatten nested lists
    flattened = []
    for item in sims:
        if isinstance(item, list):
            flattened.extend(flatten_simulation_list(item))
        else:
            flattened.append(item)

    return flattened

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


def load_serialized_config(filepath: str) -> Dict:
    """
    Load a serialized configuration file. Handles both JSON and pickle formats.

    Args:
        filepath: Path to configuration file (.json or .pkl)

    Returns:
        Configuration dictionary
    """
    if filepath.endswith('.json'):
        with open(filepath, 'r') as f:
            return json.load(f)
    elif filepath.endswith('.pkl'):
        with open(filepath, 'rb') as f:
            return pickle.load(f)
    else:
        raise ValueError(f"Unsupported file format: {filepath}")


def print_dict(dict):
    print(json.dumps(serialize_for_json(dict), indent=4))

def write_dict_to_file(dict, fname,
                       fdir='optim_result/',
                       serialize_objects = True,
                       truncate = True):
    # serialized_dict = nserialize(dict, serialize_objects=serialize_objects, truncate = truncate)
    serialized_dict = serialize_for_json(dict)
    fullname = os.path.join(fdir, fname)
    with open(fullname+'.json', 'w') as file:
        json.dump(serialized_dict, file, indent=4)
    print('Dict written to ', fullname + '.json')

    with open(fullname+'.pkl', 'wb') as file:
        pickle.dump(dict, file)
    print('Dict written to ', fullname + '.pkl')

def get_dict_from_file(fname,
                       fdir = 'optim_result/'):
    with open(fdir+fname+'.pkl', 'rb') as file:
        loaded_dict = pickle.load(file)
    return loaded_dict

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
    import re

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
    import re
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
    import re
    parts = re.findall(r"\['([^']+)'\]", path)
    if not parts:
        return None

    value = config
    for part in parts[1:]:  # Skip 'root'
        value = value[part]
    return value


def check_diff_two_dicts(d1, d2):
    """
    Legacy function for backward compatibility.
    Compare two dictionaries and print differences.

    Args:
        d1: First dictionary
        d2: Second dictionary
    """
    compare_dicts(d1, d2, format_output=True, print_result=True)

def update_config(dconf: Dict, param_dict: Dict[str, float]) -> Dict:
    """
    Update configuration with parameter values from a dictionary.

    Args:
        dconf: Base configuration dictionary
        param_dict: Dictionary mapping 'Component+parameter' to value
                   Example: {'Phy+mu_max': 1.37, 'Macroflocs+resuspension_rate': 75941.93}

    Returns:
        Updated configuration dictionary with parameters applied

    Raises:
        ValueError: If parameter name format is invalid (missing '+')

    Example:
        >>> new_config = update_config(base_config, {
        ...     'Phy+mu_max': 1.37,
        ...     'Phy+mortrate': 0.01,
        ...     'Macroflocs+resuspension_rate': 75941.93
        ... })
    """
    # Build nested update dictionary
    updates = {}
    for param_key, value in param_dict.items():
        if '+' not in param_key:
            raise ValueError(
                f"Invalid parameter name format: '{param_key}' "
                f"(expected 'Component+parameter', e.g., 'Phy+mu_max')"
            )

        component, param = param_key.split('+', 1)  # maxsplit=1 to handle '+' in param names
        if component not in updates:
            updates[component] = {'parameters': {}}
        updates[component]['parameters'][param] = value

    return deep_update(dconf, updates)


# From https://github.com/samuelcolvin/pydantic/blob/fd2991fe6a73819b48c906e3c3274e8e47d0f761/pydantic/utils.py#L200
# (initially from https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth)
from typing import (Any, Dict, TypeVar)
KeyType = TypeVar('KeyType')
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

def get_nested_attr(obj, attr):
    """
    Recursively gets the nested attribute of an object.

    Args:
        obj: The object to get the attribute from.
        attr: A string representing the attribute, potentially nested, e.g., 'prop2.C'.

    Returns:
        The value of the nested attribute.
    """
    attributes = attr.split('.')
    for attribute in attributes:
        obj = getattr(obj, attribute)
    return obj

def getlimT(T, A_E=4500, T_ref=288.15, boltz=False, bound_temp_to_1=True, T_max=None):
    """
    Get the temperature limitation factor from the temperature dependence

    Args:
        T: Temperature (K)
        A_E: Slope of the Arrhenius relation (=4500 K)
        T_ref: Reference temperature for A_E relation (Schartau et al 2007)
        boltz: Use Boltzmann constant if True
        bound_temp_to_1: If True, bound limitation between 0 and 1 (default: True)
        T_max: Maximum temperature for normalization (from setup.T_max)

    Returns:
        Temperature limitation factor [0,1] if bounded, unlimited if not
    """
    if boltz:
        factor = A_E / varinfos.boltz
    else:
        factor = A_E

    limT_raw = np.exp(-factor * (1/T - 1/T_ref))

    if not bound_temp_to_1:
        return limT_raw

    # Bounding requires T_max
    if T_max is None:
        raise ValueError("T_max required for temperature bounding (bound_temp_to_1=True)")

    limT_max = np.exp(-factor * (1/T_max - 1/T_ref))

    return np.minimum(limT_raw / limT_max, 1.0)

def print_coupled_attrs(obj):
    """Print all attributes of an object that start with 'coupled_'"""
    coupled_attrs = [attr for attr in dir(obj) if attr.startswith('coupled_')]
    print(f"\nCoupled attributes for {obj.__class__.__name__}:")
    for attr in coupled_attrs:
        print(f"{attr}: {getattr(obj, attr)}")


def cleantext(utxt):
    if '^' in utxt or '_' in utxt:
        return r"${}$".format(utxt)
    else:
        utxt = utxt.replace('-1', '$^{-1}$').replace('-2', '$^{-2}$').replace('-3', '$^{-3}$')
        utxt = utxt.replace('µ', '$\it{µ}$').replace('Chla', r'Chl$\mathit{a}$')
        return utxt


# import math
def get_all_contributors(contributors, term, pool=None, dkey=None):
    # dkey: if the term is a dictionary (e.g. same term for different targets)
    try:
        iterator = iter(contributors)
    except TypeError:
        if dkey is not None:
            if pool is not None:
                allcontributions = getattr(getattr(contributors, term), pool)[dkey]
            else:
                allcontributions = getattr(contributors, term)[dkey]
        else:
            if pool is not None:
                allcontributions = getattr(getattr(contributors, term), pool)
            else:
                allcontributions = getattr(contributors, term)
    else:
        if dkey is not None:
            if pool is not None:
                allcontributions = np.sum(getattr(getattr(contributor, term), pool)[dkey] for contributor in contributors)
            else:
                allcontributions = np.sum(getattr(contributor, term)[dkey] for contributor in contributors)
        else:
            if pool is not None:
                allcontributions = np.sum(getattr(getattr(contributor, term), pool) for contributor in contributors)
            else:
                allcontributions = np.sum(getattr(contributor, term) for contributor in contributors)
    return(allcontributions)

# supported operators
operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
             ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
             ast.USub: op.neg}

def eval_expr(expr, subdf, fulldf, setup=phys.Setup, varinfos=varinfos):
    # adapted from https://stackoverflow.com/questions/2371436/evaluating-a-mathematical-expression-in-a-string
    if expr.startswith('SUMALL'):
        df = fulldf[expr.replace('SUMALL(', '').replace(')', '')]
        if isinstance(df, pd.DataFrame):
            return df.sum(axis=1)
        else:
            return df
    else:
        return neweval(ast.parse(expr, mode='eval').body, subdf, fulldf, setup, varinfos)

def neweval(node, subdf, fulldf, setup, varinfos):
    if isinstance(node, ast.Num): # <number>
        return node.n
    elif isinstance(node, ast.BinOp): # <left> <operator> <right>
        # print(node.left, node.right)
        return operators[type(node.op)](neweval(node.left, subdf, fulldf, setup, varinfos),
                                        neweval(node.right, subdf, fulldf, setup, varinfos))
    elif isinstance(node, ast.UnaryOp): # <operator> <operand> e.g., -1
        return operators[type(node.op)](neweval(node.operand, subdf, fulldf, setup, varinfos))
    elif isinstance(node, ast.Name):
        if 'SUMALL' in node.id:
            df = fulldf[node.id.replace('SUMALL', '')]
            if isinstance(df, pd.DataFrame):
                return df.sum(axis=1)
            else:
                return df
        else:
            return subdf[node.id]
    elif isinstance(node, ast.Call):
        if isinstance(node.func, ast.Name): # for built-in functions like 'max'
            return eval(node.func.id)([neweval(a, subdf, fulldf, setup, varinfos) for a in node.args])
        elif isinstance(node.func, ast.Attribute): # for a call to a module etc e.g. np.exp
            # for 'np.exp', node.func.value.id = 'np' and node.func.attr = 'exp'
            return eval(node.func.value.id+'.'+node.func.attr)([neweval(a, subdf, fulldf, setup, varinfos) for a in node.args])[0]
            # [0] needed because somehow it returns a list...
    elif isinstance(node, ast.Attribute):
        # When e.g. setup.I is called...
        if node.value.id == 'setup':
            return setup.__getattribute__(node.attr)
        elif node.value.id == 'varinfos':
            return varinfos.__getattribute__(node.attr)
    else:
        print('Error with node ', node)
        print(node.value.id, node.attr)
        raise TypeError(node)

def getmodkwargs(newtmax=80, full_dia=False):
    # setup = phys.Setup(**phys.DEFAULT_SETUPS['onur22'], PARfromfile=True, tmax=newtmax, dt=1e-2, dt2=1e-3)
    setup = phys.Setup(**phys_setup.MOW1)
    modkwargs = {'setup': setup, 'verbose': True, 'do_diagnostics': True, 'full_diagnostics': full_dia}
    return modkwargs


if __name__ == "__main__":
    setup = phys.Setup()
    print(setup.__getattribute__('I'))

