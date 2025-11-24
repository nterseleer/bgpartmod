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
    Assumes both model_results and observations have DatetimeIndex.

    Args:
        model_results: Model results DataFrame with DatetimeIndex
        observations: Observation data object with DatetimeIndex
        daily_mean: Whether to use daily means
        method: 'aggregate' (default) aggregates model to obs periods,
                'interpolate' interpolates model to obs times
        _cached_obs: Optional cached observation data

    Returns:
        DataFrame with merged model and observation data
    """
    # Prepare model data with daily mean if requested
    model_data = model_results.resample('D').mean() if daily_mean else model_results

    # Use cached observations if available
    obs_data = _cached_obs if _cached_obs is not None else observations.df

    # Handle different temporal resolutions
    if method == 'interpolate':
        # Option 1: Interpolate model to observation times
        model_reworked = model_data.reindex(obs_data.index, method='nearest')

    else:  # method == 'aggregate' (default)
        # Option 2: Aggregate model to match observation periods

        # Detect spacing from observation index (in timedelta)
        if len(obs_data) > 1:
            obs_spacing = pd.Series(np.diff(obs_data.index)).mean()
        else:
            obs_spacing = pd.Timedelta(days=1)  # Fallback for single observation

        # Aggregate model data to observation periods
        aggregated_model = []
        for obs_time in obs_data.index:
            # Find model days in the period centered on this observation
            # Use semi-open interval (period_start, period_end] to match pd.cut() behavior
            # and avoid double-counting values at boundaries
            period_start = obs_time - obs_spacing / 2
            period_end = obs_time + obs_spacing / 2
            period_mask = (model_data.index > period_start) & (model_data.index <= period_end)
            period_model = model_data[period_mask]

            if len(period_model) > 0:
                # Compute mean over the period
                period_mean = period_model.mean()
                period_mean.name = obs_time
                aggregated_model.append(period_mean)

        if aggregated_model:
            model_reworked = pd.DataFrame(aggregated_model)
        else:
            # No overlap - return empty DataFrame
            return pd.DataFrame()

    merged_data = pd.merge(
        model_reworked, obs_data,
        left_index=True, right_index=True,
        how='inner', suffixes=('_MOD', '_OBS')
    )

    return merged_data


def calculate_likelihood(
        model_results: Any, #pd.DataFrame,
        observations: Any,
        calibrated_vars: Optional[List[str]] = None,
        daily_mean: bool = True,
        plot: bool = False,
        name: Optional[str] = None,
        save_plots: bool = False,
        plot_size: Tuple[int, int] = (10, 6),
        verbose: bool = True,
        veryverbose: bool = False,
        _cached_obs: Optional[pd.DataFrame] = None
) -> Optional[float]:
    """[docstring]"""
    if calibrated_vars is None:
        calibrated_vars = [
            'Phy_Chl', 'NH4_concentration', 'NO3_concentration',
            'DIP_concentration', 'DSi_concentration', 'Phy_C', 'TEPC_C'
        ]

    # Use streamlined data preparation
    merged_data = prepare_likelihood_data(
        model_results.df,
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
        daily_mean: bool = True
) -> Dict[str, float]:
    """
    Calculate Root Mean Square Error between model and Observations.
    Assumes both have DatetimeIndex.

    Args:
        model_results: DataFrame containing model results with DatetimeIndex
        observations: Observations object with .df attribute and DatetimeIndex
        variables: List of variables to calculate RMSE for
        daily_mean: Whether to compare daily means

    Returns:
        Dictionary of RMSE values by variable
    """
    if variables is None:
        variables = observations.df.columns

    model_data = model_results.resample('D').mean() if daily_mean else model_results.copy()

    combined_data = pd.merge(
        model_data,
        observations.df,
        left_index=True,
        right_index=True,
        how='left',
        suffixes=('_MOD', '_OBS')
    )

    rmse_values = {}
    for var in variables:
        # Handle both suffixed and non-suffixed columns
        mod_col = f'{var}_MOD' if f'{var}_MOD' in combined_data.columns else var
        obs_col = f'{var}_OBS' if f'{var}_OBS' in combined_data.columns else var

        if mod_col in combined_data.columns and obs_col in combined_data.columns:
            squared_diff = (combined_data[mod_col] - combined_data[obs_col]) ** 2
            rmse = np.sqrt(np.nanmean(squared_diff))
            rmse_values[var] = rmse

    return rmse_values

