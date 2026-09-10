"""Helpers used by the model itself.

Limitation functions, nested-attribute access, summation of coupled source/sink terms,
and the small expression evaluator behind the derived output variables (`varinfos.doutput`).

Tooling for the configuration dictionaries is in `config_tools.py`.
"""
import ast
import operator as op
import numpy as np
import pandas as pd
from typing import Any, Dict, List, Union

from src.core import phys
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

def cleantext(utxt):
    if '^' in utxt or '_' in utxt:
        return r"${}$".format(utxt)
    else:
        utxt = utxt.replace('-1', '$^{-1}$').replace('-2', '$^{-2}$').replace('-3', '$^{-3}$')
        utxt = utxt.replace('µ', r'$\it{µ}$').replace('Chla', r'Chl$\mathit{a}$')
        return utxt

def get_all_contributors(contributors, term, pool=None, dkey=None):
    """
    Sum contributions from single or multiple contributors.
    Optimized: uses native sum() and single-pass iteration.

    Args:
        contributors: Single contributor or iterable of contributors
        term: Attribute name to retrieve
        pool: Optional nested attribute (e.g., 'C', 'N', 'P')
        dkey: Optional dictionary key
    """
    try:
        iter(contributors)
    except TypeError:
        # Single contributor
        result = getattr(contributors, term)
        if pool is not None:
            result = getattr(result, pool)
        return result[dkey] if dkey is not None else result

    # Multiple contributors - single-pass iteration (optimized)
    values = []
    for c in contributors:
        val = getattr(c, term)
        if pool is not None:
            val = getattr(val, pool)
        if dkey is not None:
            val = val[dkey]
        values.append(val)

    return sum(values)  # Native sum is ~90% faster than np.sum for small lists

operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
             ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
             ast.USub: op.neg}

def eval_expr(expr, subdf, fulldf, setup=phys.Setup, varinfos=varinfos, model=None):
    # adapted from https://stackoverflow.com/questions/2371436/evaluating-a-mathematical-expression-in-a-string
    if expr.startswith('SUMALL'):
        df = fulldf[expr.replace('SUMALL(', '').replace(')', '')]
        if isinstance(df, pd.DataFrame):
            return df.sum(axis=1)
        else:
            return df
    else:
        return neweval(ast.parse(expr, mode='eval').body, subdf, fulldf, setup, varinfos, model)

def neweval(node, subdf, fulldf, setup, varinfos, model=None):
    if isinstance(node, ast.Constant):  # <number> (ast.Num before Python 3.8)
        return node.value
    elif isinstance(node, ast.BinOp): # <left> <operator> <right>
            return operators[type(node.op)](neweval(node.left, subdf, fulldf, setup, varinfos, model),
                                        neweval(node.right, subdf, fulldf, setup, varinfos, model))
    elif isinstance(node, ast.UnaryOp): # <operator> <operand> e.g., -1
        return operators[type(node.op)](neweval(node.operand, subdf, fulldf, setup, varinfos, model))
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
            return eval(node.func.id)([neweval(a, subdf, fulldf, setup, varinfos, model) for a in node.args])
        elif isinstance(node.func, ast.Attribute): # for a call to a module etc e.g. np.exp
            # for 'np.exp', node.func.value.id = 'np' and node.func.attr = 'exp'
            return eval(node.func.value.id+'.'+node.func.attr)([neweval(a, subdf, fulldf, setup, varinfos, model) for a in node.args])[0]
            # [0] needed because somehow it returns a list...
    elif isinstance(node, ast.Attribute):
        # Handle attribute access for setup, varinfos, and model
        if isinstance(node.value, ast.Name):
            # Simple attribute: setup.X, varinfos.X, model.X
            if node.value.id == 'setup':
                return setup.__getattribute__(node.attr)
            elif node.value.id == 'varinfos':
                return varinfos.__getattribute__(node.attr)
            elif node.value.id == 'model' and model is not None:
                return getattr(model, node.attr)
            else:
                # Treat as a DataFrame column with dotted name (e.g., Phy_source_PP.C)
                col_name = f"{node.value.id}.{node.attr}"
                return subdf[col_name]
        elif isinstance(node.value, ast.Subscript):
            # Handle model.components['ComponentName'].attribute pattern
            if (isinstance(node.value.value, ast.Attribute) and
                isinstance(node.value.value.value, ast.Name) and
                node.value.value.value.id == 'model' and
                node.value.value.attr == 'components' and
                model is not None):
                # Extract component name from subscript
                if isinstance(node.value.slice, ast.Constant):
                    component_name = node.value.slice.value
                else:
                    raise ValueError(f"Unsupported subscript type for model.components")

                # Safely access component and attribute
                if component_name in model.components:
                    component = model.components[component_name]
                    return getattr(component, node.attr)
                else:
                    raise KeyError(f"Component '{component_name}' not found in model.components")
        # If we reach here, the ast.Attribute is not a pattern we explicitly handle
        # Let it fall through to the error handler below
        raise TypeError(f"Unsupported ast.Attribute pattern: {ast.dump(node)}")
    else:
        print('Error with node ', node)
        if hasattr(node, 'value') and hasattr(node, 'attr'):
            print(node.value.id, node.attr)
        raise TypeError(node)
