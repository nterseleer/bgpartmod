import ast
import operator as op
import re
import importlib
import os
import numpy as np
import pandas as pd
import json
import pickle
from datetime import datetime
from deepdiff import DeepDiff
from typing import Any, Dict, List, Optional, Union

from . import phys
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

def check_diff_two_dicts(d1, d2):
    # Compare the original and imported dictionaries
    diff = DeepDiff(d1, d2, ignore_order=True)
    # Print the differences (if any)
    if diff:
        print("Differences found:")
        print(diff)
    else:
        print("No differences found. The dictionaries are identical.")

def update_config(dconf, parnames, parameters):
    newdict = {}
    for parname, parval in zip(parnames, parameters):
        keymod, par = parname.split('+')
        if keymod not in newdict:
            newdict[keymod] = {'parameters': {}}
        newdict[keymod]['parameters'][par] = parval
    return deep_update(dconf, newdict)


def update_config_from_optimization(base_config: Dict, best_parameters: Dict) -> Dict:
    """
    Update configuration with optimized parameters from an optimization result.

    Args:
        base_config: Base configuration dictionary
        best_parameters: Dictionary of optimized parameters (from optimization.summary['best_parameters'])

    Returns:
        Updated configuration dictionary
    """
    # Transform flat parameter dictionary into nested structure for deep_update
    updates = {}
    for param_key, value in best_parameters.items():
        if '+' not in param_key:
            continue  # Skip non-standard parameter keys

        component, param = param_key.split('+')
        if component not in updates:
            updates[component] = {'parameters': {}}
        updates[component]['parameters'][param] = value

    # Apply the updates to the base configuration
    return deep_update(base_config.copy(), updates)


# From https://github.com/samuelcolvin/pydantic/blob/fd2991fe6a73819b48c906e3c3274e8e47d0f761/pydantic/utils.py#L200
# (initially from https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth)
from typing import (Any, Dict, TypeVar)
KeyType = TypeVar('KeyType')
def deep_update(mapping: Dict[KeyType, Any], *updating_mappings: Dict[KeyType, Any]) -> Dict[KeyType, Any]:
    updated_mapping = mapping.copy()
    for updating_mapping in updating_mappings:
        for k, v in updating_mapping.items():
            if k in updated_mapping and isinstance(updated_mapping[k], dict) and isinstance(v, dict):
                updated_mapping[k] = deep_update(updated_mapping[k], v)
            else:
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

def getlimT(T, A_E=4500, T_ref=288.15, boltz=False):
    """
    Get the temperature limitation factor from the temperature dependence
    :param T: Temperature (K)
    :param A_E: Slope of the Arrhenius relation (=4500 K)
    :param T_ref: Reference temperature for A_E relation (Schartau et al 2007)
    :return:
    """
    if boltz:
        return(np.exp(-A_E/varinfos.boltz * (1 / T - 1 / T_ref)))
    else:
        return np.exp(-A_E * (1 / T - 1 / T_ref))

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

def eval_expr(expr, subdf, fulldf, setup=phys.Setup):
    # adapted from https://stackoverflow.com/questions/2371436/evaluating-a-mathematical-expression-in-a-string
    if expr.startswith('SUMALL'):
        df = fulldf[expr.replace('SUMALL(', '').replace(')', '')]
        if isinstance(df, pd.DataFrame):
            return df.sum(axis=1)
        else:
            return df
    else:
        return neweval(ast.parse(expr, mode='eval').body, subdf, fulldf, setup)

def neweval(node, subdf, fulldf, setup):
    if isinstance(node, ast.Num): # <number>
        return node.n
    elif isinstance(node, ast.BinOp): # <left> <operator> <right>
        # print(node.left, node.right)
        return operators[type(node.op)](neweval(node.left, subdf, fulldf, setup),
                                        neweval(node.right, subdf, fulldf, setup))
    elif isinstance(node, ast.UnaryOp): # <operator> <operand> e.g., -1
        return operators[type(node.op)](neweval(node.operand, subdf, fulldf, setup))
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
            return eval(node.func.id)([neweval(a, subdf, fulldf, setup) for a in node.args])
        elif isinstance(node.func, ast.Attribute): # for a call to a module etc e.g. np.exp
            # for 'np.exp', node.func.value.id = 'np' and node.func.attr = 'exp'
            return eval(node.func.value.id+'.'+node.func.attr)([neweval(a, subdf, fulldf, setup) for a in node.args])[0]
            # [0] needed because somehow it returns a list...
    elif isinstance(node, ast.Attribute):
        # When e.g. setup.I is called...
        if node.value.id == 'setup':
            return setup.__getattribute__(node.attr)
    else:
        print('Error with node ', node)
        print(node.value.id, node.attr)
        raise TypeError(node)


if __name__ == "__main__":
    setup = phys.Setup()
    print(setup.__getattribute__('I'))

