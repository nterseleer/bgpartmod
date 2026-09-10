"""Public fallbacks for the configuration modules a user is expected to provide.

`src.config_model.vars_to_plot`, `src.config_model.config_diagnostics` and
`src.config_model.plot_config` describe *your* choices: which variables you plot, which
diagnostics you save, how your figures are sized. They are therefore not shipped with the
framework -- you write your own and drop them in `src/config_model/`.

The modules here are the minimal versions the library falls back to when you have not.
They keep `src/utils/` importable and usable out of the box; they are not meant to be
edited in place.
"""
