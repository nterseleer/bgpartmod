"""
Evaluation metrics for BGC model results, including likelihood calculations
for optimization and general validation metrics.
"""

import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Any, Tuple
from src.utils import plotting


def prepare_likelihood_data(
        model_results: pd.DataFrame,
        observations: Any,
        climatology: bool = True,
        _cached_obs: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    Streamlined data preparation for likelihood calculation.

    Args:
        model_results: Model results DataFrame
        observations: Observation data object
        climatology: Whether to use climatological means
        _cached_obs: Optional cached observation data

    Returns:
        DataFrame with merged model and observation data
    """
    # Prepare model data - direct operation on input
    if climatology and model_results.index.name != 'julian_day':
        model_data = model_results.copy()
        model_data['julian_day'] = model_data.index.dayofyear
        model_data = model_data.groupby('julian_day').mean()
    else:
        model_data = model_results

    # Use cached Observations if available
    if _cached_obs is not None:
        obs_data = _cached_obs
    else:
        obs_data = observations.df
        if climatology and obs_data.index.name != 'julian_day':
            obs_data = obs_data.copy()
            obs_data['julian_day'] = obs_data.index.dayofyear
            obs_data = obs_data.groupby('julian_day').mean()

    # Merge and return
    return pd.merge(
        model_data,
        obs_data,
        left_index=True,
        right_index=True,
        how='left',
        suffixes=('_MOD', '_OBS')
    )


def calculate_likelihood(
        model_results: pd.DataFrame,
        observations: Any,
        calibrated_vars: Optional[List[str]] = None,
        climatology: bool = True,
        plot: bool = False,
        save_plots: bool = False,
        plot_size: Tuple[int, int] = (10, 6),
        verbose: bool = True,
        veryverbose: bool = False,
        _cached_obs: Optional[pd.DataFrame] = None
) -> Optional[float]:
    """[docstring unchanged]"""
    if calibrated_vars is None:
        calibrated_vars = [
            'Phy_Chl', 'NH4_concentration', 'NO3_concentration',
            'DIP_concentration', 'DSi_concentration', 'Phy_C', 'TEPC_C'
        ]

    # Use streamlined data preparation
    merged_data = prepare_likelihood_data(
        model_results,
        observations,
        climatology,
        _cached_obs
    )

    # Create comparison plots if requested (only during analysis, not optimization)
    if plot:
        plotting.plot_results(
            model_results,
            variables=calibrated_vars,
            observations=observations,
            climatology=climatology,
            figsize=plot_size,
            save=save_plots,
            filename='likelihood_comparison'
        )

    # Calculate likelihood
    total_likelihood = 0.0
    for var in calibrated_vars:
        obs = merged_data[f'{var}_OBS']
        obs_count = obs.count()

        if obs_count > 0:
            squared_diff = np.nansum(
                (merged_data[f'{var}_MOD'] - obs) ** 2
            )
            var_likelihood = -0.5 * obs_count * np.log(squared_diff)

            if veryverbose:
                print(f'{var}: {var_likelihood}')

            if not np.isnan(var_likelihood):
                total_likelihood += var_likelihood

    if verbose:
        print(f'Total log-likelihood: {total_likelihood}')

    return total_likelihood


def calculate_rmse(
        model_results: pd.DataFrame,
        observations: Any,
        variables: Optional[List[str]] = None,
        climatology: bool = True
) -> Dict[str, float]:
    """
    Calculate Root Mean Square Error between model and Observations.

    Args:
        model_results: DataFrame containing model results
        observations: Observations object with .df attribute
        variables: List of variables to calculate RMSE for
        climatology: Whether to compare climatological means

    Returns:
        Dictionary of RMSE values by variable
    """
    if variables is None:
        variables = observations.df.columns

    model_data = model_results.copy()
    if climatology:
        model_data['julian_day'] = model_data.index.dayofyear
        model_data = model_data.groupby('julian_day').mean()

    combined_data = pd.merge(
        model_data,
        observations.df,
        left_index=True,
        right_index=True,
        how='left'
    )

    rmse_values = {}
    for var in variables:
        if var in combined_data.columns:
            squared_diff = (combined_data[f'{var}_MOD'] - combined_data[var]) ** 2
            rmse = np.sqrt(np.nanmean(squared_diff))
            rmse_values[var] = rmse

    return rmse_values

