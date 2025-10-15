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
        daily_mean: bool = True,
        _cached_obs: Optional[pd.DataFrame] = None,
        method: str = 'aggregate',  # 'aggregate' or 'interpolate'
) -> pd.DataFrame:
    """
    Streamlined data preparation for likelihood calculation.

    Args:
        model_results: Model results DataFrame
        observations: Observation data object
        daily_mean: Whether to use daily means
        method: 'aggregate' (default) aggregates model to obs periods,
                'interpolate' interpolates model to obs times
        _cached_obs: Optional cached observation data

    Returns:
        DataFrame with merged model and observation data
    """
    # Prepare model data - direct operation on input
    if daily_mean and model_results.index.name != 'julian_day':
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
        if daily_mean and obs_data.index.name != 'julian_day':
            obs_data = obs_data.copy()
            obs_data['julian_day'] = obs_data.index.dayofyear
            obs_data = obs_data.groupby('julian_day').mean()

    # Handle different temporal resolutions
    if method == 'interpolate':
        # Option 1: Interpolate model to observation times
        model_reworked = model_data.reindex(obs_data.index, method='nearest')

    else:  # method == 'aggregate' (default)
        # Option 2: Aggregate model to match observation periods

        # Detect period_days from observation index spacing
        if len(obs_data) > 1:
            obs_spacing = np.diff(obs_data.index).mean()
            period_days = int(round(obs_spacing))
        else:
            period_days = 1  # Fallback for single observation

        # Aggregate model data to observation periods
        aggregated_model = []
        for obs_time in obs_data.index:
            # Determine the period this observation represents
            period_start = obs_time - period_days / 2 + 0.5
            period_end = obs_time + period_days / 2 + 0.5

            # Find model days in this period
            period_mask = (model_data.index >= period_start) & (model_data.index <= period_end)
            period_model = model_data[period_mask]

            if len(period_model) > 0:
                # Compute mean over the period
                period_mean = period_model.mean()
                period_mean.name = obs_time
                aggregated_model.append(period_mean)

        if aggregated_model:
            model_reworked = pd.DataFrame(aggregated_model)
        else:
            # No overlap - return empty DataFrame with proper structure
            return pd.DataFrame()

    merged_data = pd.merge(
        model_reworked, obs_data,
        left_index=True, right_index=True,
        how='inner', suffixes=('_MOD', '_OBS')
    )

    return merged_data


def calculate_likelihood(
        model_results: pd.DataFrame,
        observations: Any,
        calibrated_vars: Optional[List[str]] = None,
        daily_mean: bool = True,
        plot: bool = False,
        save_plots: bool = False,
        plot_size: Tuple[int, int] = (10, 6),
        verbose: bool = True,
        veryverbose: bool = False,
        _cached_obs: Optional[pd.DataFrame] = None,
        file_name: str = 'likelihood_comparison'
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
        daily_mean,
        _cached_obs
    )

    # Create comparison plots if requested (only during analysis, not optimization)
    if plot:
        plotting.plot_results(
            model_results,
            variables=calibrated_vars,
            observations=observations,
            daily_mean=daily_mean,
            figsize=plot_size,
            save=save_plots,
            filename=file_name
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
        daily_mean: bool = True
) -> Dict[str, float]:
    """
    Calculate Root Mean Square Error between model and Observations.

    Args:
        model_results: DataFrame containing model results
        observations: Observations object with .df attribute
        variables: List of variables to calculate RMSE for
        daily_mean: Whether to compare daily means

    Returns:
        Dictionary of RMSE values by variable
    """
    if variables is None:
        variables = observations.df.columns

    model_data = model_results.copy()
    if daily_mean:
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

