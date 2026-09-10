"""Extension point for replaying archived configurations.

A saved configuration may name parameters that have since been renamed or removed. They are
translated here, once, before the components are instantiated -- component signatures stay
strict, so a typo in a configuration still raises a TypeError.

This repository ships no translation table: the configuration it comes with
(`base_config.Kerimoglu2022`) only names current parameters. To replay your own archived
configurations, drop a `_legacy_tables.py` module next to this file, exposing:

    NON_COMPONENT_KEYS                      # top-level config keys that are not components
    translate_parameters(component, params) # -> params with legacy names translated

Without it, both are inert and configurations pass through unchanged.
"""
try:
    from src.core import _legacy_tables as _tables
except ImportError:
    _tables = None

# Top-level keys of a configuration dictionary that do not designate a component.
NON_COMPONENT_KEYS = getattr(_tables, 'NON_COMPONENT_KEYS', frozenset())


def translate_parameters(component, parameters):
    """Return `parameters` with legacy names translated, or unchanged if no table is set.

    component : component name, used in the messages.
    parameters : the 'parameters' sub-dictionary of a configuration (not modified).
    """
    if _tables is None:
        return dict(parameters)
    return _tables.translate_parameters(component, parameters)
