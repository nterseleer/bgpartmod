"""Fallback diagnostics levels. See src/config_model/_defaults/__init__.py.

A level is nothing but ``{'<Component>': {'diagnostics': [...]}}``, merged into a model
configuration with ``fns.deep_update``. Diagnostics are stored at every timestep, so the
level is the main lever on a simulation's memory footprint.

Only what `src/utils/optimization.py` needs to work without a user-supplied module.
"""

# Attached by Optimization.get_best_model to the best model, so that a saved best model
# carries enough to be plotted. Components absent from the configuration are filtered out
# by the caller, so listing the mineral pools here is harmless in a BGC-only run.
plotting = {
    'Phy': {'diagnostics': ['kd', 'PAR_t', 'limI', 'limT', 'limNUT',
                            'PC_max', 'PC', 'source_PP.C']},
    'Microflocs': {'diagnostics': ['massconcentration']},
    'Micro_in_Macro': {'diagnostics': ['massconcentration']},
    'Macroflocs': {'diagnostics': ['diam', 'settling_vel']},
}

# EWMA state of the organic-to-floc coupling ratios. Hidden state: updated at every slow
# step in get_sink_vertical_loss, never part of the ODE pool vector, and therefore NOT
# recoverable from a saved DataFrame unless requested here. Required on any run that will
# be cropped or restarted from, otherwise the restart re-seeds the ratio from the
# instantaneous C/Nf and drifts.
smoothed_ratios = {
    'TEPC': {'diagnostics': ['smoothed_C_to_Nf_ratio']},
    'DetL': {'diagnostics': ['smoothed_C_to_Nf_ratio', 'smoothed_N_to_Nf_ratio',
                             'smoothed_P_to_Nf_ratio', 'smoothed_Si_to_Nf_ratio']},
    'DetS': {'diagnostics': ['smoothed_C_to_Nf_ratio', 'smoothed_N_to_Nf_ratio',
                             'smoothed_P_to_Nf_ratio', 'smoothed_Si_to_Nf_ratio']},
    'BacA': {'diagnostics': ['smoothed_C_to_Nf_ratio', 'smoothed_N_to_Nf_ratio',
                             'smoothed_P_to_Nf_ratio']},
}
