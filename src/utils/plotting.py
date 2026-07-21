"""
Plotting utilities for BGC model results.
Handles both single and multiple simulation visualizations.
"""
import os
from functools import lru_cache

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import numpy as np
import pandas as pd
from typing import Union, List, Dict, Optional, Tuple, Any

from . import functions as fns
from src.config_model import varinfos
from src.config_model import vars_to_plot
from src.config_system import path_config as path_cfg
from src.utils import observations
from src.utils.plotted_variables_sets import PlottedVariablesSet

FIGURE_PATH = path_cfg.FIGURE_PATH

# Plotting constants
# Extended color palette with 14 distinct, colorblind-friendly colors
# Optimized for scientific plots with good contrast and distinguishability
DEFAULT_COLORS = [
    '#1f77b4',  # Blue
    '#ff7f0e',  # Orange
    '#2ca02c',  # Green
    '#d62728',  # Red
    '#9467bd',  # Purple
    '#8c564b',  # Brown
    '#e377c2',  # Pink
    '#17becf',  # Cyan
    '#bcbd22',  # Yellow-green
    '#7f7f7f',  # Gray
    '#c49c94',  # Tan
    '#f7b6d2',  # Light pink
    '#c5b0d5',  # Lavender
    '#9edae5',  # Light cyan
]

# DEFAULT_COLORS = [
#     'darkolivegreen',
#     'purple',
#     'darkgoldenrod',
#
# ]

# Line styles: 3 distinct styles to combine with 14 colors
# Avoiding dot-only ':' style as it becomes too discrete
# Using solid, dashed, and dash-dot for good visibility
DEFAULT_LINESTYLES = ['-', '--', '-.', ':', '-', ]
# DEFAULT_LINESTYLES = ['-', '-', '--',]

# Budget plot color palettes (cool colors for sources, warm colors for sinks)
BUDGET_SOURCE_COLORS = ['#2E86AB', '#06A77D', '#81B214', '#A7C957', '#4ECDC4', '#45B7D1']
# BUDGET_SINK_COLORS = ['#D62828', '#F77F00', '#FCBF49', '#EE6C4D', '#E63946', '#F4A261']
BUDGET_SINK_COLORS = ['#FFD300', '#E1A95F', '#E08D3C', '#E25822', '#E34234', '#8B0000']

DEFAULT_LEGEND_FONTSIZE = 8
DEFAULT_FIGURE_DPI = 500

# General plotting setup
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
# Cycle through colors first (all 14 colors with solid line),
# then repeat with different line styles
# This ensures first 14 curves are maximally distinct
# We use product of cyclers: 14 colors × 3 linestyles = 42 unique combinations
from cycler import cycler
plt.rcParams['axes.prop_cycle'] = (
    cycler('color', DEFAULT_COLORS) *
    cycler('linestyle', DEFAULT_LINESTYLES)
)
plt.rcParams.update({'font.size': 8})

# Define observation styling
OBS_STYLES = {
    'calibrated': {
        'marker': '*',
        'color': 'gray',
        'markersize': 6,
        'alpha': 0.6
    },
    'non_calibrated': {
        'marker': 'o',
        'color': 'gray',
        'markersize': 3,
        'alpha': 0.4
    }
}

# Styles for multi-window plotting (when mean_window_days is a tuple/list)
# Order: background (smoothed) -> foreground (detailed)
# Use None for color to preserve each model's original color
MULTI_WINDOW_STYLES = [
    {'linewidth': 0.5, 'alpha': 0.7, 'color': None},  # Layer 0: background
    {'linewidth': 1.5, 'alpha': 1.0, 'color': None},  # Layer 1: foreground
    {'linewidth': 2.0, 'alpha': 1.0, 'color': None},  # Layer 2: optional
]


@lru_cache(maxsize=32)
def _get_cached_model_styles(n_models: int) -> Tuple[Dict, ...]:
    """Generate cached model styles for performance."""
    return tuple({
        'color': DEFAULT_COLORS[i % len(DEFAULT_COLORS)],
        'linestyle': DEFAULT_LINESTYLES[i % len(DEFAULT_LINESTYLES)],
        'linewidth': 2.
    } for i in range(n_models))


def get_model_styles(n_models: int, custom_styles: Optional[List[Dict]] = None) -> List[Dict]:
    """
    Generate consistent model styles for plotting (cached for performance).

    Args:
        n_models: Number of models to generate styles for
        custom_styles: Optional custom styles to use instead of defaults

    Returns:
        List of style dictionaries for each model
    """
    if custom_styles and len(custom_styles) >= n_models:
        return custom_styles[:n_models]

    return list(_get_cached_model_styles(n_models))

def _filter_by_time(df: pd.DataFrame, time_filter) -> pd.DataFrame:
    """
    Filter DataFrame by time index.

    Args:
        df: DataFrame with DatetimeIndex
        time_filter: Time filtering criterion
            - None: no filtering
            - str (year): e.g., "2023" filters to that year
            - tuple: (start, end) date strings
            - slice: pandas slice for advanced filtering

    Returns:
        Filtered DataFrame
    """
    if time_filter is None:
        return df

    if isinstance(time_filter, str):
        # Assume it's a year
        return df[df.index.year == int(time_filter)]

    elif isinstance(time_filter, tuple) and len(time_filter) == 2:
        # Date range
        start, end = time_filter
        return df.loc[start:end]

    elif isinstance(time_filter, slice):
        # Direct slice
        return df.loc[time_filter]

    else:
        raise ValueError(f"Unsupported time_filter type: {type(time_filter)}")


def _expand_time_filter(time_filter, buffer_days: int):
    """
    Expand a time_filter by adding a buffer on both sides.
    This avoids edge effects when computing temporal averages.
    """
    if time_filter is None or buffer_days <= 0:
        return time_filter

    buffer = pd.Timedelta(days=buffer_days)

    if isinstance(time_filter, str):
        # Year string -> expand to tuple
        year = int(time_filter)
        start = pd.Timestamp(f"{year}-01-01") - buffer
        end = pd.Timestamp(f"{year}-12-31") + buffer
        return (start, end)

    elif isinstance(time_filter, tuple) and len(time_filter) == 2:
        start = pd.Timestamp(time_filter[0]) - buffer
        end = pd.Timestamp(time_filter[1]) + buffer
        return (start, end)

    # slice or other: return unchanged (can't easily expand)
    return time_filter


def prepare_model_obs_data(
        models: Union[Any, List[Any], pd.DataFrame, List[pd.DataFrame]],
        observations: Optional[Any] = None,
        mean_window_days: Optional[int] = 1,
        daily_mean: Optional[bool] = None,
        variables_to_plot: Optional[List[str]] = None,
        time_filter = None,
        center_trend: bool = True
) -> Tuple[List[pd.DataFrame], Optional[pd.DataFrame], Optional[pd.DataFrame], List[str]]:
    """
    Prepare model and observation data for plotting.
    Assumes both models and observations have DatetimeIndex.

    Performance optimization: Only copies necessary columns from potentially huge DataFrames.

    Args:
        models: Single model/DataFrame or list of models/DataFrames (with DatetimeIndex)
        observations: Observation data object (with DatetimeIndex)
        mean_window_days: Window size in days for temporal averaging (1 = daily, 7 = weekly, 14 = bi-weekly, etc.)
                         None or 0 disables averaging. Default: 1 (daily mean)
        daily_mean: Deprecated. Use mean_window_days=1 instead. Kept for backward compatibility.
        variables_to_plot: List of variables that will be plotted (for column standardization)
        time_filter: Time filtering criterion (None, year string, date range tuple, or slice)
        center_trend: Averaging mode when mean_window_days > 0. This is the preferred,
            more correct default across the plotting stack.
            True (default): mean_window_days-wide bin means labelled at the bin CENTRE
                (the instant the mean actually represents -- so features are not shifted
                ~half a window early as with left-edge labelling), then linearly
                interpolated onto a daily grid so the trend spans the exact data edges
                while staying smooth (one value per window -> no residual spring-neap
                ripple). Bins remain epoch-anchored, so phase alignment across simulations
                is preserved. The buffer is widened to a full window so the edge values are
                genuine two-sided (centred) means, matching the centred averaging used for
                scoring; it uses data on both sides of time_filter when available (spin-up
                before, a short run extension after -- see Setup.extend_duration).
            False: legacy epoch-anchored resampling labelled at the bin LEFT edge (sparse,
                truncates ~half a window short at each end, features shifted early). Kept as
                an explicit opt-out.

    Returns:
        Tuple of (model_data_list, merged_data, full_obs_data, model_names)
    """
    # Legacy support: daily_mean overrides mean_window_days if explicitly provided
    if daily_mean is not None:
        mean_window_days = 1 if daily_mean else None
    # Flatten nested lists of models
    models = fns.flatten_simulation_list(models)

    # Determine which columns we actually need to avoid copying massive DataFrames
    required_columns = set()
    if variables_to_plot:
        required_columns.update(variables_to_plot)

    # Extract DataFrames and names
    model_data_list = []
    model_names = []

    for i, model in enumerate(models):
        if isinstance(model, pd.DataFrame):
            full_df = model
            name = f'Model {i + 1}'
        else:
            full_df = model.df
            name = getattr(model, 'name', f'Model {i + 1}')

        # Compute on-the-fly any requested derived variable that has an 'oprt' in
        # the live output config (varinfos.doutput) but is missing from the model
        # df -- e.g. a variable added to varinfos after the simulation was built or
        # saved to a pickle. Needs the model object since eval_expr may reference
        # model./setup. Raw DataFrame inputs cannot be completed this way and fall
        # through to the NaN-fill standardization below.
        if variables_to_plot and not isinstance(model, pd.DataFrame):
            missing = [v for v in variables_to_plot if v not in full_df.columns]
            if missing and hasattr(model, 'compute_derived_variables'):
                model.compute_derived_variables(missing, output_config=varinfos.doutput,
                                                verbose=False)
                full_df = model.df  # refreshed with the newly computed column(s)

        # Performance optimization: only copy necessary columns
        if required_columns:
            # Find available columns from the required set
            available_required = [col for col in required_columns if col in full_df.columns]
            if available_required:
                model_data = full_df[available_required].copy()
            else:
                # Fallback: if no required columns found, take all (shouldn't happen normally)
                model_data = full_df.copy()
        else:
            # If no specific columns requested, copy all (backward compatibility)
            model_data = full_df.copy()

        # Sanitize +/-inf -> NaN on the plotting copy (model.df untouched). Handles
        # ratio-type variables that are inf where a denominator hits 0 (e.g.
        # Phy_limI/Phy_limI_theoretical at night, where limI is floored to 1e-6 but
        # the denominator is 0), including columns loaded from older pickles built
        # before the compute-time guard existed. Without this, a single inf in a
        # resampling window turns the whole averaged value into inf and the curve vanishes.
        model_data = model_data.replace([np.inf, -np.inf], np.nan)

        # Handle duplicate names by adding a counter
        original_name = name
        counter = 1
        while name in model_names:
            counter += 1
            name = f"{original_name} ({counter})"

        # Apply time filtering with buffer to avoid edge effects in temporal averaging.
        # center_trend needs a FULL window of neighbours on each side (so the edge bins
        # are complete and the interpolated edge values are two-sided); the plain bin
        # resample only needs half a window.
        if mean_window_days is not None and mean_window_days > 0:
            buffer_days = mean_window_days if center_trend else (mean_window_days + 1) // 2
            expanded_filter = _expand_time_filter(time_filter, buffer_days)
            if expanded_filter is not None:
                model_data = _filter_by_time(model_data, expanded_filter)
            # Apply temporal averaging.
            if center_trend:
                # Smooth biweekly-mean trend that reaches the exact data edges WITHOUT
                # exposing the residual spring-neap ripple (a mean_window_days boxcar does
                # not fully cancel the ~14.77 d beat, and a daily-dense rolling mean would
                # resolve -- i.e. show -- that residual):
                #   (1) mean_window_days-wide bin means -- one honest value per window, so
                #       there is no sub-window structure to ripple -- labelled at the bin
                #       CENTRE (the value's true time), matching the centred averaging used
                #       for scoring (evaluation.prepare_and_merge);
                #   (2) linearly interpolated onto a daily grid (pure densification: adds no
                #       new structure, identical look to connecting the bin means by lines)
                #       so the drawn line spans the whole window;
                #   (3) cropped to time_filter below.
                # With buffer_days == mean_window_days the bins bracketing each edge are
                # complete, so the interpolated first/last values are genuine two-sided means.
                binned = model_data.resample(f'{mean_window_days}D', origin='epoch').mean()
                binned.index = binned.index + pd.Timedelta(days=mean_window_days / 2.0)
                binned = binned.dropna(how='all')
                if len(binned) >= 2:
                    daily_idx = pd.date_range(binned.index[0], binned.index[-1], freq='1D')
                    model_data = binned.reindex(binned.index.union(daily_idx)) \
                                       .interpolate(method='time').reindex(daily_idx)
                else:
                    model_data = binned
            else:
                # origin='epoch' anchors the resampling bins to a fixed global grid (1970-01-01)
                # instead of the default per-series start_day. This is essential when overlaying
                # simulations that start on different dates (e.g. a full run vs a cropped restart):
                # otherwise each series gets its own bin grid and any sub-window periodicity
                # (notably the ~14.7 d spring-neap residual in bed_shear_stress / water_depth)
                # appears phase-shifted between curves even though the underlying signal is identical.
                model_data = model_data.resample(f'{mean_window_days}D', origin='epoch').mean()
            # Crop back to original time_filter
            if time_filter is not None:
                model_data = _filter_by_time(model_data, time_filter)
        elif time_filter is not None:
            model_data = _filter_by_time(model_data, time_filter)

        model_data_list.append(model_data)
        model_names.append(name)

    # Standardize columns across all DataFrames if variables_to_plot is provided
    if variables_to_plot:
        all_columns = set()
        for df in model_data_list:
            all_columns.update(df.columns)

        # Add missing columns with NaN for variables that will be plotted
        for var in variables_to_plot:
            if var not in all_columns:
                continue

            for i, df in enumerate(model_data_list):
                if var not in df.columns:
                    model_data_list[i] = df.assign(**{var: np.nan})

    if observations is None:
        return model_data_list, None, None, model_names

    # Prepare observation data (assumes DatetimeIndex, no automatic cycling)
    obs_data = observations.df.copy()

    # Apply same time filtering to observations
    if time_filter is not None:
        obs_data = _filter_by_time(obs_data, time_filter)

    # Merge with first model - the DatetimeIndex merge will handle temporal alignment automatically
    merged_data = pd.merge(
        model_data_list[0],
        obs_data,
        left_index=True,
        right_index=True,
        how='outer',
        suffixes=('_MOD', '_OBS')
    )

    return model_data_list, merged_data, obs_data, model_names


def plot_variable(
        ax: plt.Axes,
        model_data_list: List[pd.DataFrame],
        model_names: List[str],
        var_name: str,
        merged_data: Optional[pd.DataFrame] = None,
        model_styles: Optional[List[Dict]] = None,
        fill_between_data: Optional[Tuple[pd.DataFrame, pd.DataFrame]] = None,
        fill_alpha: float = 0.3,
        fill_color: str = 'magenta',
        plot_obs: bool = True,
        obs_kwargs: Optional[Dict] = None,
        obs_kwargs_calibrated: Optional[Dict] = None,
        obs_kwargs_non_calibrated: Optional[Dict] = None,
        add_labels: bool = True,
        calibrated_vars: Optional[List[str]] = None,
        show_subplot_titles: bool = True,
        apply_plt_ylim: bool = True,
        ylabel_pad: Optional[float] = None,
        ylabel_fontsize: Optional[float] = None,
        tick_label_pad: Optional[float] = None,
        tick_labelsize: Optional[float] = None,
        tick_length: Optional[float] = None,
        tick_width: Optional[float] = None,
        spine_width: Optional[float] = None
) -> None:
    """Plot a single variable with model results and optional observations."""
    # Default styles
    if model_styles is None:
        model_styles = get_model_styles(len(model_data_list))

    # Default observation kwargs
    base_obs_kwargs = obs_kwargs or {}

    # Plot fill_between variability range first (behind everything)
    if fill_between_data is not None:
        low_data, high_data = fill_between_data
        if var_name in low_data.columns and var_name in high_data.columns:
            # Use common index (intersection) for fill_between
            common_index = low_data.index.intersection(high_data.index)
            ax.fill_between(
                common_index,
                low_data.loc[common_index, var_name],
                high_data.loc[common_index, var_name],
                alpha=fill_alpha,
                color=fill_color,
                zorder=1,
                label='Uncertainty range' if add_labels else None
            )

    # Plot each model
    for model_data, name, style in zip(model_data_list, model_names, model_styles):
        if var_name in model_data.columns:
            ax.plot(model_data.index, model_data[var_name],
                    label=name if add_labels else "_" + name, **style)


    if plot_obs and merged_data is not None and f'{var_name}_OBS' in merged_data.columns:
        # Determine observation style based on calibration status
        is_calibrated = calibrated_vars is not None and var_name in calibrated_vars
        obs_style_key = 'calibrated' if is_calibrated else 'non_calibrated'
        selected_obs_style = OBS_STYLES[obs_style_key]

        # Build final style with proper priority:
        # 1. Default style from OBS_STYLES (base)
        # 2. Common obs_kwargs (overrides defaults)
        # 3. Specific obs_kwargs_calibrated/non_calibrated (highest priority)
        final_obs_style = {**selected_obs_style, **base_obs_kwargs}

        # Apply calibration-specific overrides if provided
        if is_calibrated and obs_kwargs_calibrated is not None:
            final_obs_style.update(obs_kwargs_calibrated)
        elif not is_calibrated and obs_kwargs_non_calibrated is not None:
            final_obs_style.update(obs_kwargs_non_calibrated)

        # Check if std data is available for error bars
        std_col = f"{var_name}_std"
        if std_col in merged_data.columns:
            # Use errorbar for observations with std
            # Extract errorbar-specific parameters with defaults
            elinewidth = final_obs_style.get('elinewidth', 1.5)
            capsize = final_obs_style.get('capsize', 3.)
            capthick = final_obs_style.get('capthick', 1.5)
            solid_capstyle = final_obs_style.get('solid_capstyle', 'round')

            # Filter out errorbar-specific params to avoid duplication
            errorbar_params = ['marker', 'linestyle', 's', 'elinewidth', 'capsize', 'capthick', 'solid_capstyle']
            filtered_style = {k: v for k, v in final_obs_style.items() if k not in errorbar_params}

            ax.errorbar(
                merged_data.index,
                merged_data[f'{var_name}_OBS'],
                yerr=merged_data[std_col],
                label='Observations' if add_labels else "_Used observations",
                fmt=final_obs_style['marker'],
                elinewidth=elinewidth,
                capsize=capsize,
                capthick=capthick,
                solid_capstyle=solid_capstyle,
                **filtered_style
            )
        else:
            ax.scatter(
                merged_data.index,
                merged_data[f'{var_name}_OBS'],
                label='Observations' if add_labels else "_Used observations",
                **{k: v for k, v in final_obs_style.items() if k != 'linestyle'}
            )

    # Apply ylim if specified and data exceeds bounds
    var_info = varinfos.doutput.get(var_name.lstrip('m'), {})
    plt_ylim = var_info.get('plt_ylim')
    if plt_ylim is not None and apply_plt_ylim:
        # Get actual data range from plotted elements
        y_data = []
        for line in ax.get_lines():
            y_data.extend(line.get_ydata())
        for collection in ax.collections:
            if len(collection.get_offsets()) > 0:
                y_data.extend(collection.get_offsets()[:, 1])

        if y_data:
            data_min, data_max = np.nanmin(y_data), np.nanmax(y_data)
            ylim_min, ylim_max = plt_ylim

            # Apply only if data exceeds specified bounds
            if data_min < ylim_min or data_max > ylim_max:
                ax.set_ylim(ylim_min, ylim_max)

    # Format labels and title
    clean_name = var_info.get('cleanname', None)

    # If cleanname not defined in varinfos, use var_name but escape underscores
    # to prevent matplotlib mathtext errors (double subscript issues)
    if clean_name is None:
        clean_name = var_name.replace('_', r'\_')

    units = var_info.get('munits' if var_name.startswith('m') else 'units', '')
    long_name = var_info.get('longname', None)

    # If longname not defined in varinfos, use var_name but escape underscores
    # to prevent matplotlib mathtext errors (double subscript issues)
    if long_name is None:
        long_name = var_name.replace('_', r'\_')

    # Set ylabel with optional custom padding and fontsize
    ylabel_kwargs = {}
    if ylabel_pad is not None:
        ylabel_kwargs['labelpad'] = ylabel_pad
    if ylabel_fontsize is not None:
        ylabel_kwargs['fontsize'] = ylabel_fontsize
    ax.set_ylabel(f'{fns.cleantext(clean_name)} [{fns.cleantext(units)}]', **ylabel_kwargs)

    # Adjust tick label padding, fontsize, length and width if specified
    tick_params_kwargs = {}
    if tick_label_pad is not None:
        tick_params_kwargs['pad'] = tick_label_pad
    if tick_labelsize is not None:
        tick_params_kwargs['labelsize'] = tick_labelsize
    if tick_length is not None:
        tick_params_kwargs['length'] = tick_length
    if tick_width is not None:
        tick_params_kwargs['width'] = tick_width
    if tick_params_kwargs:
        ax.tick_params(axis='both', **tick_params_kwargs)

    # Adjust spine (frame) width if specified
    if spine_width is not None:
        for spine in ax.spines.values():
            spine.set_linewidth(spine_width)

    if show_subplot_titles:
        ax.set_title(fns.cleantext(long_name))


def plot_results(
        models: Union[Any, List[Any]],
        variables: Union[List[str], PlottedVariablesSet],
        observations: Optional[Any] = observations.Obs(station='MOW1_biweekly_202509_noPhaeo'),
        calibrated_vars: Optional[List[str]] = None,
        mean_window_days: Optional[int] = 1,
        daily_mean: Optional[bool] = None,
        center_trend: bool = True,
        time_filter = None,
        ncols: Optional[int] = None,
        plot_obs = True,
        figsize: Optional[Tuple[int, int]] = None,
        model_styles: Optional[List[Dict]] = None,
        multi_window_styles: Optional[List[Dict]] = None,
        fill_between_models: Optional[Tuple[Any, Any]] = None,
        fill_alpha: float = 0.3,
        fill_color: str = 'magenta',
        subplot_labels: bool = True,
        subplot_label_fontsize: int = 10,
        label_position: str = 'top_left',
        show_subplot_titles: bool = True,
        apply_plt_ylim: bool = True,
        date_format: str = '%d/%m',
        xlabel_rotation: float = 45,
        ylabel_pad: Optional[float] = None,
        ylabel_fontsize: Optional[float] = None,
        tick_label_pad: Optional[float] = None,
        tick_labelsize: Optional[float] = None,
        tick_length: Optional[float] = None,
        tick_width: Optional[float] = None,
        spine_width: Optional[float] = None,
        tight_layout: bool = True,
        hspace: Optional[float] = None,
        wspace: Optional[float] = None,
        save: bool = False,
        filename: Optional[str] = None,
        figdir: Optional[str] = None,
        fnametimestamp: bool = True,
        **plot_kwargs
) -> Tuple[plt.Figure, np.ndarray]:
    """
    Create comprehensive plots of model results with optional observations.

    Args:
        models: Single model or list of models
        variables: List of variables to plot OR PlottedVariablesSet with embedded preferences
        observations: Optional observation data
        calibrated_vars: Calibrated variables (have different style than non-calibrated vars)
        mean_window_days: Window size in days for temporal averaging (1 = daily, 7 = weekly, 14 = bi-weekly, etc.)
        daily_mean: Deprecated. Use mean_window_days=1 instead. Kept for backward compatibility.
        center_trend: Forwarded to prepare_model_obs_data (applied to every mean_window
            layer). If True (default), each averaging window is a centre-labelled bin mean
            linearly interpolated to daily: the trend is placed at the CENTRE of its window
            (correct time position, aligned with the daily layer) instead of the bin left
            edge, and reaches the exact data edges instead of truncating ~half a window
            short at each end. Set False for the legacy left-edge sparse binning. Uses data
            on both sides of time_filter when available (spin-up before, run extension after).
        time_filter: Time filtering criterion (None, year string, date range tuple, or slice)
        ncols: Number of columns in subplot grid (auto-detected from PlottedVariablesSet if None)
        figsize: Figure size (auto-detected from PlottedVariablesSet if None, otherwise auto-calculated)
        model_styles: List of style dictionaries for each model
        fill_between_models: Optional tuple of (low_model, high_model) to create variability range with fill_between
        fill_alpha: Transparency of fill_between area (default: 0.3)
        fill_color: Color of fill_between area (default: 'gray')
        subplot_labels: Whether to add subplot labels (a, b, c, etc.)
        subplot_label_fontsize: Font size for subplot labels
        label_position: Position for subplot labels
        show_subplot_titles: Whether to show titles above each subplot
        date_format: Format string for x-axis dates (e.g., '%d/%m', '%b %d', '%Y-%m-%d')
        xlabel_rotation: Rotation angle for x-axis labels in degrees
        ylabel_pad: Padding between ylabel and tick labels (None = matplotlib default, typically 4.0)
        ylabel_fontsize: Font size for y-axis labels (None = use rcParams default)
        tick_label_pad: Padding between tick marks and tick labels (None = matplotlib default, typically 3.0)
        tick_labelsize: Font size for tick labels on both axes (None = use rcParams default)
        tick_length: Length of tick marks in points (None = matplotlib default, typically 3.5)
        tick_width: Width of tick marks in points (None = matplotlib default, typically 0.8)
        spine_width: Width of subplot frame (spines) in points (None = matplotlib default, typically 0.8)
        tight_layout: Whether to apply tight_layout for better spacing
        hspace: Height space between subplots (None = matplotlib default)
        wspace: Width space between subplots (None = matplotlib default)
        save: Whether to save the figure
        filename: Filename for saved figure (auto-generated from PlottedVariablesSet.name if None)
        figdir: directory for saved figures
        fnametimestamp: Whether to add timestamp to saved filename
        **plot_kwargs: Additional plotting parameters

    Returns:
        Figure and axes array

    Raises:
        ValueError: If input parameters are invalid
    """
    # Extract PlottedVariablesSet attributes if provided
    if isinstance(variables, PlottedVariablesSet):
        variable_set = variables
        variables_list = list(variable_set)

        # Use PlottedVariablesSet preferences if not explicitly overridden
        if ncols is None:
            ncols = variable_set.ncols
        if figsize is None:
            figsize = variable_set.figsize

        # Override legend position in plot_kwargs if specified in PlottedVariablesSet
        if 'legend_position' not in plot_kwargs and hasattr(variable_set, 'legend_position'):
            plot_kwargs['legend_position'] = variable_set.legend_position
    else:
        variables_list = variables
        # Defaults for regular list input
        if ncols is None:
            ncols = 2

    # Input validation
    if not variables_list:
        raise ValueError("Variables list cannot be empty")
    if ncols < 1:
        raise ValueError("Number of columns must be positive")
    if label_position not in ['top_left', 'top_right', 'bottom_left', 'bottom_right']:
        raise ValueError(f"Invalid label_position: {label_position}")
    # Extract calibrated_vars from models if not provided
    if calibrated_vars is None:
        # Flatten models list to check for calibrated_vars
        flat_models = fns.flatten_simulation_list(models)
        for model in flat_models:
            if hasattr(model, 'calibrated_vars') and model.calibrated_vars:
                calibrated_vars = model.calibrated_vars
                break

    # Normalize mean_window_days to list for multi-window support
    if isinstance(mean_window_days, (tuple, list)):
        mean_window_values = list(mean_window_days)
    else:
        mean_window_values = [mean_window_days]

    # Calculate grid dimensions
    nrows = (len(variables_list) + ncols - 1) // ncols
    if figsize is None:
        figsize = (5 * ncols, 4 * nrows)

    # Create figure and axes
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True)
    if nrows == 1 and ncols == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # Set intelligent figure title from PlottedVariablesSet
    if isinstance(variables, PlottedVariablesSet):
        fig.canvas.manager.set_window_title(variables.name.replace('_', ' '))

    # Extract legend-related and obs-specific kwargs before passing to plot_variable
    legend_kwargs = ['add_legend', 'legend_fontsize', 'legend_position']
    obs_specific_kwargs = ['obs_kwargs_calibrated', 'obs_kwargs_non_calibrated']
    plot_var_kwargs = {k: v for k, v in plot_kwargs.items()
                       if k not in legend_kwargs}

    # Get base model styles
    base_styles = model_styles if model_styles else get_model_styles(
        len(fns.flatten_simulation_list(models)))

    # Per-window layer styles (linewidth/alpha for each mean_window_days value).
    # Falls back to the module default; can be overridden via plot_config.
    window_styles = multi_window_styles if multi_window_styles else MULTI_WINDOW_STYLES

    # Track fill_between for legend
    fill_between_data = None

    # Plot each layer (background to foreground)
    for layer_idx, mwd in enumerate(mean_window_values):
        is_first_layer = (layer_idx == 0)
        is_last_layer = (layer_idx == len(mean_window_values) - 1)

        # Prepare data for this mean_window_days value
        layer_data_list, layer_merged, _, model_names = prepare_model_obs_data(
            models, observations, mwd, daily_mean, variables_list, time_filter,
            center_trend=center_trend
        )

        # Prepare fill_between data (only for first layer)
        layer_fill_between = None
        if is_first_layer and fill_between_models is not None:
            fill_data_list, _, _, _ = prepare_model_obs_data(
                list(fill_between_models), None, mwd, daily_mean, variables_list, time_filter,
                center_trend=center_trend
            )
            if len(fill_data_list) == 2:
                layer_fill_between = (fill_data_list[0], fill_data_list[1])
                fill_between_data = layer_fill_between  # Track for legend

        # Apply layer-specific style overrides
        # Offset index so that the last layer always uses the foreground style
        n_layers = len(mean_window_values)
        style_offset = max(0, len(window_styles) - n_layers)
        layer_style_config = window_styles[min(layer_idx + style_offset, len(window_styles) - 1)]
        layer_overrides = {k: v for k, v in layer_style_config.items() if v is not None}
        layer_styles = [{**s, **layer_overrides} for s in base_styles]

        # Plot each variable
        for i, var in enumerate(variables_list):
            if i < len(axes):
                plot_variable(
                    axes[i],
                    layer_data_list,
                    model_names,
                    var,
                    layer_merged,
                    layer_styles,
                    fill_between_data=layer_fill_between,
                    fill_alpha=fill_alpha,
                    fill_color=fill_color,
                    add_labels=(i == 0 and is_last_layer),
                    plot_obs=(plot_obs and is_first_layer),
                    calibrated_vars=calibrated_vars,
                    show_subplot_titles=(show_subplot_titles and is_last_layer),
                    apply_plt_ylim=apply_plt_ylim,
                    ylabel_pad=ylabel_pad,
                    ylabel_fontsize=ylabel_fontsize,
                    tick_label_pad=tick_label_pad,
                    tick_labelsize=tick_labelsize,
                    tick_length=tick_length,
                    tick_width=tick_width,
                    spine_width=spine_width,
                    **plot_var_kwargs
                )

    # Format axes and add subplot labels (once, after all layers)
    for i, var in enumerate(variables_list):
        if i < len(axes):
            axes[i].xaxis.set_major_formatter(mdates.DateFormatter(date_format))
            for label in axes[i].get_xticklabels():
                label.set_rotation(xlabel_rotation)
            if subplot_labels:
                add_subplot_label(axes[i], i, position=label_position, fontsize=subplot_label_fontsize)
        else:
            axes[i].set_visible(False)

    # Add global legend if needed
    if plot_kwargs.get('add_legend', True):
        # Create legend based on all models provided, ensuring consistent styling
        handles, labels = [], []

        # Get consistent styles for all models
        legend_styles = model_styles if model_styles else get_model_styles(len(model_names))

        # Create legend entries for all models
        for i, name in enumerate(model_names):
            style = legend_styles[i] if i < len(legend_styles) else legend_styles[-1]

            # Create a dummy line for the legend with the correct style
            line = plt.Line2D([0], [0], **style)
            handles.append(line)
            labels.append(name)

        # Add variability range if present
        if fill_between_data is not None:
            patch = plt.Rectangle((0, 0), 1, 1, fc=fill_color, alpha=fill_alpha)
            handles.append(patch)
            labels.append('Variability')

        # Add observations if present (check all subplots)
        if layer_merged is not None and len(axes) > 0:
            # Look for observations in any subplot
            for ax in axes:
                if hasattr(ax, 'get_legend_handles_labels'):
                    h, l = ax.get_legend_handles_labels()
                    for handle, label in zip(h, l):
                        if 'Observations' in label and 'Observations' not in labels:
                            handles.append(handle)
                            labels.append('Observations')
                            break
                    if 'Observations' in labels:
                        break

        # Create the legend
        legend_pos = plot_kwargs.get('legend_position', 'center right')

        # Set bbox_to_anchor based on legend position
        bbox_anchor_map = {
            'center right': (0.98, 0.5),
            'upper right': (0.98, 0.98),
            'lower right': (0.98, 0.02),
            'upper left': (0.02, 0.98),
            'lower left': (0.02, 0.02),
            'upper center': (0.5, 0.98),
            'lower center': (0.5, 0.02)
        }
        bbox_anchor = bbox_anchor_map.get(legend_pos, (0.98, 0.5))

        fig.legend(
            handles,
            labels,
            loc=legend_pos,
            bbox_to_anchor=bbox_anchor,
            fontsize=plot_kwargs.get('legend_fontsize', DEFAULT_LEGEND_FONTSIZE)
        )

    # Apply layout adjustments
    if tight_layout:
        plt.tight_layout()

    # Apply custom spacing if provided
    if hspace is not None or wspace is not None:
        plt.subplots_adjust(
            hspace=hspace if hspace is not None else plt.rcParams['figure.subplot.hspace'],
            wspace=wspace if wspace is not None else plt.rcParams['figure.subplot.wspace']
        )

    # Save figure using centralized function
    if save:
        if figdir is None:
            figdir = FIGURE_PATH
        # Extract PlottedVariablesSet for intelligent naming
        var_set = variables if isinstance(variables, PlottedVariablesSet) else None
        save_figure(fig, filename=filename, variable_set=var_set, figdir=figdir, add_timestamp=fnametimestamp)

    return fig, axes


# Specialized plotting functions
def plot_nutrients(model, observations=None, **kwargs):
    """Plot nutrient concentrations."""
    return plot_results(model, vars_to_plot.nutvars, observations, **kwargs)


def plot_stoichiometry(model, observations=None, **kwargs):
    """Plot stoichiometric ratios."""
    return plot_results(model, vars_to_plot.stoichioPhy, observations, **kwargs)


def plot_biomass(model, observations=None, **kwargs):
    """Plot biomass variables."""
    return plot_results(model, vars_to_plot.phyvars, observations, **kwargs)


def plot_sinks_sources(
        models: Union[Any, List[Any]],
        sources: List[str],
        sinks: List[str],
        time_filter = None,
        observations: Optional[Any] = None,
        mean_window_days: Optional[int] = 1,
        daily_mean: Optional[bool] = None,
        increase_resolution_factor: int = 2,
        figsize: Optional[Tuple[float, float]] = None,
        default_subplot_size: Tuple[float, float] = (4.6, 5.7),
        auto_adjust_figsize: bool = True,
        maxrows: int = 3,
        legend_loc: int = 0,
        legend_fontsize: int = 6,
        subplot_labels: bool = True,
        label_position: str = 'top_left',
        save: bool = False,
        filename: str = 'Sources_vs_Sinks',
        fnametimestamp: bool = False,
        **kwargs
) -> Tuple[plt.Figure, List[pd.DataFrame]]:
    """
    Plot sources and sinks as filled areas using fill_between, with a net balance line.

    Args:
        models: Single model or list of models
        sources: List of column names to be plotted as sources (positive values)
        sinks: List of column names to be plotted as sinks (negative values)
        observations: Optional observation data
        mean_window_days: Window size in days for temporal averaging (1 = daily, 7 = weekly, 14 = bi-weekly, etc.)
        daily_mean: Deprecated. Use mean_window_days=1 instead. Kept for backward compatibility.
        increase_resolution_factor: Factor to increase resolution of the datetime index
        figsize: Custom figure size (if None, calculated based on default_subplot_size)
        default_subplot_size: Default size for a single subplot when auto-calculating figsize
        auto_adjust_figsize: Whether to automatically adjust figsize based on grid dimensions
        maxrows: Maximum number of rows in the subplot grid
        legend_loc: Legend location
        legend_fontsize: Legend font size
        subplot_labels: Whether to add subplot labels (a, b, c, etc.)
        label_position: Position for subplot labels
        save: Whether to save the figure
        filename: Name for saved figure
        **kwargs: Additional plotting parameters

    Returns:
        Tuple of (figure, list of high-resolution DataFrames)
    """
    # Prepare model data
    model_data_list, _, _, model_names = prepare_model_obs_data(
        models, observations, mean_window_days, daily_mean, time_filter=time_filter,
    )

    # Create subplots
    # Calculate optimal grid dimensions
    num_models = len(model_data_list)
    if num_models <= maxrows:
        # If we have fewer models than maxrows, use a single column
        nrows, ncols = num_models, 1
    else:
        # Otherwise calculate optimal columns and rows
        ncols = (num_models + maxrows - 1) // maxrows  # Ceiling division
        # Recalculate rows to use exactly what we need
        nrows = (num_models + ncols - 1) // ncols
        nrows = min(nrows, maxrows)  # Ensure we don't exceed maxrows
        # If we have too few models for the last row, recalculate columns
        if nrows * ncols - num_models >= ncols:
            ncols = (num_models + nrows - 1) // nrows

        # Determine figure size
    if figsize is None or auto_adjust_figsize:
        # Calculate figsize based on grid dimensions
        used_figsize = (default_subplot_size[0] * ncols, default_subplot_size[1] * nrows)
    else:
        # Use user-provided figsize
        used_figsize = figsize

    # Create subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=used_figsize, sharex=True)

    # Flatten axes for easy iteration
    flat_axes = axes.flatten() if num_models > 1 else [axes]

    # Process each DataFrame
    df_high_res_list = []
    for i, (df, ax, model_name) in enumerate(zip(model_data_list, flat_axes, model_names)):
        # Validate input columns
        available_cols = set(df.columns)
        missing_sources = [col for col in sources if col not in available_cols]
        missing_sinks = [col for col in sinks if col not in available_cols]

        if missing_sources or missing_sinks:
            print(f"Warning for {model_name}: Missing columns:")
            if missing_sources:
                print(f"  Sources: {missing_sources}")
            if missing_sinks:
                print(f"  Sinks: {missing_sinks}")

        # Use only available columns
        valid_sources = [col for col in sources if col in available_cols]
        valid_sinks = [col for col in sinks if col in available_cols]

        if not valid_sources or not valid_sinks:
            print(f"Error: No valid sources or sinks for {model_name}. Skipping plot.")
            continue

        # Create high-resolution index for smoother visualization
        new_index = pd.date_range(
            start=df.index.min(),
            end=df.index.max(),
            periods=len(df.index) * increase_resolution_factor
        )
        df_high_res = pd.DataFrame(index=new_index).join(
            df[valid_sources + valid_sinks].astype(float),
            how='outer'
        ).interpolate(method='time')

        df_high_res_list.append(df_high_res)

        # Calculate cumulative sums
        cumulative_sources = df_high_res[valid_sources].cumsum(axis=1)
        cumulative_sinks = df_high_res[valid_sinks].cumsum(axis=1)

        # Calculate net balance
        netbalance = cumulative_sources.iloc[:, -1] - cumulative_sinks.iloc[:, -1]
        df_high_res['netbalance'] = netbalance

        # Plot sources (stacked positive values)
        for j, source in enumerate(valid_sources):
            if j == 0:
                ax.fill_between(
                    df_high_res.index,
                    0,
                    cumulative_sources[source],
                    color=BUDGET_SOURCE_COLORS[j % len(BUDGET_SOURCE_COLORS)],
                    label=source,
                    alpha=0.7
                )
            else:
                ax.fill_between(
                    df_high_res.index,
                    cumulative_sources[valid_sources[j - 1]],
                    cumulative_sources[source],
                    color=BUDGET_SOURCE_COLORS[j % len(BUDGET_SOURCE_COLORS)],
                    label=source,
                    alpha=0.7
                )

        # Plot sinks (stacked negative values)
        for j, sink in enumerate(valid_sinks):
            if j == 0:
                ax.fill_between(
                    df_high_res.index,
                    0,
                    -cumulative_sinks[sink],
                    color=BUDGET_SINK_COLORS[j % len(BUDGET_SINK_COLORS)],
                    label=sink,
                    alpha=0.7
                )
            else:
                ax.fill_between(
                    df_high_res.index,
                    -cumulative_sinks[valid_sinks[j - 1]],
                    -cumulative_sinks[sink],
                    color=BUDGET_SINK_COLORS[j % len(BUDGET_SINK_COLORS)],
                    label=sink,
                    alpha=0.7
                )

        # Plot net balance line
        ax.plot(
            df_high_res.index,
            netbalance.where(netbalance > 0),
            color='black',
            linewidth=1.5,
            label='Positive net growth'
        )
        ax.plot(
            df_high_res.index,
            netbalance.where(netbalance < 0),
            color='black',
            linewidth=1.5,
            linestyle='--',
            label='Negative net growth'
        )

        # Format axes
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m'))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

        # Add labels
        if model_name:
            ax.set_title(model_name)
        ax.set_ylabel(fns.cleantext('mmolC') + ' ' + fns.cleantext('m^{-3}') + ' ' + fns.cleantext('d^{-1}'))

        # Add legend only to first subplot
        if i == 0:
            ax.legend(loc=legend_loc, fontsize=legend_fontsize)

        # Add subplot label if requested
        if subplot_labels:
            add_subplot_label(ax, i, position=label_position)
    plt.tight_layout()

    # Save figure if requested
    if save:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S_')
        fname = f"{filename}_{timestamp}.png" if fnametimestamp else filename
        save_figure(fig, filename=fname)

    return fig, df_high_res_list


# Wrapper function for simplified budget plotting
def plot_budget(
        models: Union[Any, List[Any]],
        budget: Dict[str, List[str]],
        time_filter = None,
        **kwargs
) -> Tuple[plt.Figure, List[pd.DataFrame]]:
    """
    Simplified wrapper for plot_sinks_sources that takes a single budget dictionary.

    Args:
        models: Single model or list of models
        budget: Dictionary with 'sources' and 'sinks' keys
        **kwargs: Additional arguments for plot_sinks_sources

    Returns:
        Result from plot_sinks_sources
    """
    return plot_sinks_sources(
        models=models,
        sources=budget['sources'],
        sinks=budget['sinks'],
        time_filter=time_filter,
        **kwargs
    )


def plot_signed_series(
        models: Union[Any, List[Any]],
        variables: List[str],
        mean_window_days: Optional[int] = 1,
        daily_mean: Optional[bool] = None,
        time_filter=None,
        ncols: int = 2,
        figsize: Optional[Tuple[float, float]] = None,
        default_subplot_size: Tuple[float, float] = (3.0, 2.0),
        positive_color: str = 'tab:blue',
        negative_color: str = 'tab:red',
        fill_alpha: float = 0.5,
        line_color: str = '0.25',
        line_width: float = 0.7,
        show_subplot_titles: bool = True,
        subplot_labels: bool = True,
        label_position: str = 'top_left',
        subplot_label_fontsize: int = 8,
        xlabel_rotation: float = 45,
        date_format: str = '%b',
        ylabel_fontsize: Optional[float] = None,
        tick_labelsize: Optional[float] = None,
        hspace: Optional[float] = None,
        wspace: Optional[float] = None,
        legend_fontsize: float = 6,
        save: bool = False,
        filename: str = 'signed_series',
        figdir: Optional[str] = None,
        fnametimestamp: bool = False,
        **kwargs,
) -> Tuple[plt.Figure, np.ndarray]:
    """Plot each variable as a sign-coloured time series over a zero baseline.

    The area between the curve and zero is filled `positive_color` where the value
    is positive and `negative_color` where it is negative (thin outline on top), so
    net/signed diagnostics read at a glance. Designed e.g. for the floc vertical
    fluxes, where Macroflocs_settling_loss = sink_sedimentation - source_resuspension
    flips sign (>0 net deposition, <0 net resuspension); strictly one-signed inputs
    (sink_sedimentation, source_resuspension) simply come out uniformly coloured.

    Data prep (smoothing via `mean_window_days`, `time_filter`, multi-model handling)
    reuses prepare_model_obs_data, so it stays consistent with plot_results/plot_budget.
    `mean_window_days` may be the (fast, slow) tuple used by plot_results; only its
    last element is used. Extra keyword arguments (e.g. the styling keys bundled in a
    shared cmd_kwargs) are accepted and ignored.
    """
    # plot_results passes a (fast, slow) overlay tuple; signed fill needs one window.
    if isinstance(mean_window_days, (tuple, list)):
        mean_window_days = mean_window_days[-1]

    model_data_list, _, _, model_names = prepare_model_obs_data(
        models, observations=None, mean_window_days=mean_window_days,
        daily_mean=daily_mean, variables_to_plot=variables, time_filter=time_filter,
    )

    nvars = len(variables)
    nrows = -(-nvars // ncols)  # ceil
    used_figsize = figsize or (default_subplot_size[0] * ncols,
                               default_subplot_size[1] * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=used_figsize, squeeze=False)
    flat_axes = axes.flatten()

    for k, var in enumerate(variables):
        ax = flat_axes[k]
        for mi, (df, name) in enumerate(zip(model_data_list, model_names)):
            if var not in df.columns:
                ax.text(0.5, 0.5, f'{var}\nnot in output', ha='center', va='center',
                        fontsize=6, transform=ax.transAxes)
                continue
            s = pd.to_numeric(df[var], errors='coerce')
            # later models get lighter fills so overlaid series stay distinguishable
            a = fill_alpha if mi == 0 else fill_alpha * 0.5
            ax.fill_between(s.index, 0, s, where=(s >= 0), interpolate=True,
                            color=positive_color, alpha=a, linewidth=0,
                            label='positive' if (mi == 0 and k == 0) else None)
            ax.fill_between(s.index, 0, s, where=(s < 0), interpolate=True,
                            color=negative_color, alpha=a, linewidth=0,
                            label='negative' if (mi == 0 and k == 0) else None)
            ax.plot(s.index, s, color=line_color, linewidth=line_width, zorder=3)
        ax.axhline(0.0, color='0.6', linewidth=0.5, zorder=2)

        # labels/units via varinfos (same convention as plot_variable)
        var_info = varinfos.doutput.get(var.lstrip('m'), {})
        clean_name = var_info.get('cleanname') or var.replace('_', r'\_')
        units = var_info.get('units', '')
        long_name = var_info.get('longname') or var.replace('_', r'\_')
        ylab_kwargs = {'fontsize': ylabel_fontsize} if ylabel_fontsize is not None else {}
        ax.set_ylabel(f'{fns.cleantext(clean_name)} [{fns.cleantext(units)}]', **ylab_kwargs)
        if show_subplot_titles:
            ax.set_title(fns.cleantext(long_name), fontsize=ylabel_fontsize)
        ax.xaxis.set_major_formatter(mdates.DateFormatter(date_format))
        if tick_labelsize is not None:
            ax.tick_params(axis='both', labelsize=tick_labelsize)
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(xlabel_rotation)
        if subplot_labels:
            add_subplot_label(ax, k, position=label_position,
                              fontsize=subplot_label_fontsize)

    for ax in flat_axes[nvars:]:  # hide unused axes
        ax.set_visible(False)

    handles, labels = flat_axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc='upper right', fontsize=legend_fontsize,
                   framealpha=0.9)

    if hspace is not None or wspace is not None:
        fig.subplots_adjust(hspace=hspace if hspace is not None else 0.2,
                            wspace=wspace if wspace is not None else 0.2)
    else:
        fig.tight_layout()

    if save:
        fname = (f"{filename}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                 if fnametimestamp else filename)
        path = os.path.join(figdir, f'{fname}.png') if figdir else f'{fname}.png'
        fig.savefig(path, dpi=300, bbox_inches='tight')

    return fig, axes


def plot_element_distribution_stacked(model_output, element_vars, element_name=None,
                                      time_var='time', group_by_compartment=True, relative=False,
                                      stacked=True, time_filter=None, mean_window_days=None,
                                      center_trend=True, save=False,
                                      filename='element_distribution_stacked',
                                      figdir=None, fnametimestamp=True, **kwargs):
    """
    Create a stacked area plot showing element distribution across model compartments

    Parameters:
    -----------
    model_output : DataFrame
        Model output containing the variables
    element_vars : list
        List of element variables to plot (e.g., vars_to_plot.all_nitrogen_vars)
    element_name : str, optional
        Name of the element for labeling (auto-detected if None)
    time_var : str, default 'time'
        Name of time variable
    group_by_compartment : bool, default True
        If True, group variables by compartment (Phytoplankton, Heterotrophs, etc.)
        If False, plot each variable individually
    relative : bool, default False
        If True, plot relative contributions (%) instead of absolute values
    stacked : bool, default True
        If True, use stackplot (stacked area). If False, overlay individual curves.
    time_filter : optional
        Time filtering criterion (None, year string e.g. "2023", (start, end) date
        tuple, or pandas slice), consistent with the rest of the plotting module.
    mean_window_days : int, optional
        Window size in days for temporal averaging (e.g. 14 = bi-weekly mean).
        None (default) disables averaging (raw data, only time-filtered).
        Routed through prepare_model_obs_data, so it reuses the same edge-buffer
        and origin='epoch' resampling logic as plot_results.
    center_trend : bool, default True
        Forwarded to prepare_model_obs_data: when averaging, place the trend at
        its window centre and interpolate to daily (smooth, edge-reaching). Set
        False for the legacy epoch-anchored left-labelled binning.
    save : bool, default False
        Whether to save the figure
    filename : str, default 'element_distribution_stacked'
        Base filename for saved figure
    figdir : str, optional
        Directory to save figure (uses FIGURE_PATH if None)
    fnametimestamp : bool, default True
        Whether to add timestamp to filename
    **kwargs : dict
        Additional plotting arguments

    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    import matplotlib.pyplot as plt

    # Auto-detect element name if not provided
    if element_name is None:
        if any('_N' in var for var in element_vars):
            element_name = 'Nitrogen'
            unit = 'mmol N m⁻³'
        elif any('_P' in var for var in element_vars):
            element_name = 'Phosphorus'
            unit = 'mmol P m⁻³'
        elif any('_Si' in var for var in element_vars):
            element_name = 'Silicon'
            unit = 'mmol Si m⁻³'
        elif any('_C' in var for var in element_vars):
            element_name = 'Carbon'
            unit = 'mmol C m⁻³'
        else:
            element_name = 'Element'
            unit = 'mmol m⁻³'
    else:
        unit = f'mmol {element_name} m⁻³'

    # Filter available variables
    available_vars = [var for var in element_vars if var in model_output.columns]

    # Route through prepare_model_obs_data to reuse the shared time-filtering AND
    # optional temporal averaging (edge-buffer + origin='epoch' resampling), exactly
    # as plot_results / plot_kd_contributions_stacked do. mean_window_days=None keeps
    # the raw, only-time-filtered behaviour. Restrict to available_vars so only the
    # needed columns are copied out of a potentially huge df.
    if available_vars:
        model_output = prepare_model_obs_data(
            model_output, observations=None,
            mean_window_days=mean_window_days,
            variables_to_plot=available_vars,
            time_filter=time_filter,
            center_trend=center_trend,
        )[0][0]

    if group_by_compartment:
        # Define compartments based on variable naming patterns
        compartments = {
            'Phytoplankton': [var for var in available_vars if var.startswith('Phy_')],
            'Heterotrophs': [var for var in available_vars if
                             any(var.startswith(x) for x in ['BacF_', 'BacA_', 'HF_', 'Cil_'])],
            'Dissolved Organic': [var for var in available_vars if any(x in var for x in ['DOCS_', 'DOCL_', 'TEPC_'])],
            'Detritus': [var for var in available_vars if var.startswith('Det')],
            'NH4': [var for var in available_vars if 'NH4' in var and 'concentration' in var],
            'NO3': [var for var in available_vars if 'NO3' in var and 'concentration' in var],
            'DIP': [var for var in available_vars if 'DIP' in var and 'concentration' in var],
            'DSi': [var for var in available_vars if 'DSi' in var and 'concentration' in var]
        }

        # Remove empty compartments and calculate totals
        compartment_totals = {}
        for comp_name, variables in compartments.items():
            if variables:  # Only include non-empty compartments
                compartment_totals[comp_name] = sum(model_output[var] for var in variables)
    else:
        # Plot each variable individually without grouping
        compartment_totals = {}
        for var in available_vars:
            # Get clean label from varinfos if available
            var_info = varinfos.doutput.get(var.lstrip('m'), {})
            clean_label = var_info.get('cleanname', var.replace('_', ' '))
            compartment_totals[clean_label] = model_output[var]

    # Create figure
    fig, ax = plt.subplots(figsize=(6, 4))

    # Create stacked area plot with distinct colors
    if compartment_totals:
        # Convert to relative contributions if requested
        if relative:
            # Convert dict values to array for vectorized operations
            data_array = np.array([series.values for series in compartment_totals.values()])
            # Calculate total at each time point
            total = np.sum(data_array, axis=0)
            # Avoid division by zero
            total = np.where(total == 0, 1, total)
            # Convert to percentages
            data_to_plot = [(contrib / total) * 100 for contrib in data_array]
        else:
            data_to_plot = list(compartment_totals.values())

        # Use DEFAULT_COLORS for better color distinction
        colors = DEFAULT_COLORS[:len(compartment_totals)]
        if stacked:
            ax.stackplot(model_output.index,
                         *data_to_plot,
                         labels=compartment_totals.keys(),
                         colors=colors,
                         alpha=0.8)
        else:
            for data, label, color in zip(data_to_plot, compartment_totals.keys(), colors):
                ax.plot(model_output.index, data, label=label, color=color, linewidth=2, alpha=0.8)

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

    # Formatting
    ax.set_xlabel('Time [d]')

    if relative:
        ax.set_ylabel(f'{element_name} relative contribution [%]')
        ax.set_title(f'{element_name} Relative Distribution Across Model Compartments')
        ax.set_ylim(0, 100)
    else:
        ax.set_ylabel(f'{element_name} concentration [{unit}]')
        ax.set_title(f'{element_name} Distribution Across Model Compartments')
        # ax.set_ylim(0, 800)

    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save if requested
    if save:
        if figdir is None:
            figdir = FIGURE_PATH
        save_figure(fig, filename=filename, figdir=figdir, add_timestamp=fnametimestamp)

    return fig, ax


def plot_kd_contributions_stacked(model_output, kd_contrib_vars=None, relative=False,
                                  stacked=True, time_filter=None, mean_window_days=None,
                                  grouped_heterotrophs=True, grouped_detritus=True, grouped_flocs=True,
                                  center_trend=True, ax=None, add_title=True, legend_kwargs=None,
                                  total_kd_style=None, **kwargs):
    """
    Create a plot showing contributions to light attenuation coefficient (kd).

    Parameters:
    -----------
    model_output : DataFrame
        Model output containing the kd contribution variables
    kd_contrib_vars : list, optional
        List of kd contribution variables to plot. If None, uses default list.
    relative : bool, optional
        If True, plot relative contributions (%) instead of absolute values.
        Default is False (absolute contributions in m^-1).
    stacked : bool, optional
        If True, use stackplot (stacked area). If False, overlay individual curves.
        Default is True (stacked visualization).
    time_filter : optional
        Time filtering criterion (None, year string, date range tuple, or slice).
        - None: no filtering
        - str (year): e.g., "2023" filters to that year
        - tuple: (start, end) date strings
        - slice: pandas slice for advanced filtering
    mean_window_days : int, optional
        Window size in days for temporal averaging (e.g. 14 = bi-weekly mean).
        None (default) disables averaging (raw data, only time-filtered).
        Routed through prepare_model_obs_data, so it reuses the same edge-buffer
        and origin='epoch' resampling logic as plot_results.
    grouped_heterotrophs, grouped_detritus, grouped_flocs : bool, optional
        Only used when kd_contrib_vars is None (default variable set). Each controls
        one family of contributions independently; all default to True:
          - grouped_heterotrophs: merge BacA, BacF, HF, Cil -> 'kd_contrib_Heterotrophs'
          - grouped_detritus:     merge DetL, DetS          -> 'kd_contrib_Detritus'
          - grouped_flocs:        merge Microflocs,
                                  Micro_in_Macro             -> 'kd_contrib_Flocs'
        The list is assembled by vars_to_plot.build_kd_contributions_list; a False
        flag expands that family back to its per-component members. All ignored if
        kd_contrib_vars is provided explicitly.
    center_trend : bool, optional
        Forwarded to prepare_model_obs_data. If True (default), the mean_window_days
        averaging is a centre-labelled bin mean linearly interpolated to daily, so the
        trend is placed at its window centre and reaches the exact data edges while staying
        smooth (no residual spring-neap ripple, edge values two-sided); set False for the
        legacy epoch-anchored left-labelled binning.
    ax : matplotlib.axes.Axes, optional
        If provided, draw into this axis instead of creating a new figure (for
        embedding the stack as a panel in a composed/paper figure). Figure-level
        cosmetics (autofmt_xdate, tight_layout) are then left to the caller.
    add_title : bool, optional
        If True (default), set the auto-generated panel title. Pass False to
        suppress it (e.g. when the caller adds its own subplot letter).
    legend_kwargs : dict, optional
        Overrides merged into the ax.legend() call (loc, bbox_to_anchor, fontsize,
        ncol, ...). None (default) keeps the outside-right placement.
    total_kd_style : dict, optional
        Style overrides for the 'Total k_d' overlay line (absolute mode only).
        None (default) keeps the dashed black line ('k--', linewidth 2). Pass e.g.
        a solid-trend style to match a multi_window line panel.
    **kwargs : dict
        Additional plotting arguments (figsize, save, etc.)

    Returns:
    --------
    fig, ax : matplotlib figure and axis objects

    Examples:
    ---------
    >>> # Absolute contributions, stacked (default)
    >>> fig, ax = plot_kd_contributions_stacked(model.df)
    >>>
    >>> # Relative contributions (%), stacked
    >>> fig, ax = plot_kd_contributions_stacked(model.df, relative=True)
    >>>
    >>> # Relative contributions (%), overlaid (not stacked)
    >>> fig, ax = plot_kd_contributions_stacked(model.df, relative=True, stacked=False)
    >>>
    >>> # Filter to year 2023, relative %
    >>> fig, ax = plot_kd_contributions_stacked(model.df, relative=True, time_filter="2023")
    """
    import matplotlib.pyplot as plt
    from src.config_model import vars_to_plot

    # Use default list if not provided (each family grouped/expanded independently)
    if kd_contrib_vars is None:
        kd_contrib_vars = vars_to_plot.build_kd_contributions_list(
            grouped_heterotrophs=grouped_heterotrophs,
            grouped_detritus=grouped_detritus,
            grouped_flocs=grouped_flocs)

    # Each grouped series is a derived sum. When the caller passes a raw DataFrame
    # (e.g. sim.df), the on-the-fly derived-variable machinery in prepare_model_obs_data
    # does not run (it only fires for model objects), so the grouped column is absent
    # and would be silently dropped by the available_vars filter. Synthesize any
    # requested-but-missing grouped series from its constituents (single source of
    # truth: vars_to_plot.kd_contrib_groups). Summing before the downstream resampling
    # is equivalent to summing after (both linear).
    for grouped_var, components in vars_to_plot.kd_contrib_groups.items():
        if (grouped_var in kd_contrib_vars and grouped_var not in model_output.columns
                and all(c in model_output.columns for c in components)):
            model_output = model_output.copy()
            model_output[grouped_var] = model_output[components].sum(axis=1)

    # Filter to available variables
    available_vars = [var for var in kd_contrib_vars if var in model_output.columns]

    if not available_vars:
        print("Warning: No kd contribution variables found in model output")
        return None, None

    # Route through prepare_model_obs_data to reuse the time-filtering AND optional
    # temporal averaging (edge-buffer + origin='epoch' resampling) shared with
    # plot_results. mean_window_days=None keeps the raw, only-time-filtered behaviour.
    # In absolute mode also carry Phy_kd through (same smoothing) so the 'Total k_d'
    # overlay line below can still be drawn -- prepare_model_obs_data restricts the
    # output to variables_to_plot, so it would otherwise be dropped.
    prep_vars = list(available_vars)
    if not relative and 'Phy_kd' in model_output.columns and 'Phy_kd' not in prep_vars:
        prep_vars.append('Phy_kd')
    model_output = prepare_model_obs_data(
        model_output, observations=None,
        mean_window_days=mean_window_days,
        variables_to_plot=prep_vars,
        time_filter=time_filter,
        center_trend=center_trend,
    )[0][0]

    # Extract figsize from kwargs or use default
    figsize = kwargs.get('figsize', (14, 8))

    # Create figure, or draw into a caller-provided axis (composed/paper figures)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        own_figure = True
    else:
        fig = ax.figure
        own_figure = False

    # Prepare data (each variable as a separate series)
    # Ensure data is numeric and handle NaN values
    data_to_plot = []
    for var in available_vars:
        # Convert to numeric, replacing any non-numeric with NaN
        series = pd.to_numeric(model_output[var], errors='coerce')
        # Fill NaN with 0
        series = series.fillna(0)
        data_to_plot.append(series.values)

    # Convert to relative contributions if requested
    if relative:
        # Calculate total at each time point
        total = np.sum(data_to_plot, axis=0)
        # Avoid division by zero
        total = np.where(total == 0, 1, total)
        # Convert to percentages
        data_to_plot = [(contrib / total) * 100 for contrib in data_to_plot]

    # Create clean labels (remove 'kd_contrib_' prefix and format nicely)
    labels = []
    for var in available_vars:
        label = var.replace('kd_contrib_', '')
        # Special formatting for certain components
        if label == 'Micro_in_Macro':
            label = 'Microflocs in Macroflocs'
        elif label == 'Flocs':
            label = 'Mineral flocs'
        labels.append(label)

    # Ensure index is numeric/datetime (not object)
    # Convert index to numeric if it's a DatetimeIndex
    if isinstance(model_output.index, pd.DatetimeIndex):
        x_values = mdates.date2num(model_output.index)
    else:
        # Assume it's already numeric (days)
        x_values = model_output.index.values

    # Create color palette
    colors = plt.cm.tab10.colors[:len(available_vars)]
    if len(available_vars) > 10:
        colors = plt.cm.tab20.colors[:len(available_vars)]

    # Plot according to stacked option
    if stacked:
        # Stacked area plot
        ax.stackplot(x_values,
                     *data_to_plot,
                     labels=labels,
                     colors=colors,
                     alpha=0.8)
    else:
        # Overlaid line plots
        for i, (data, label, color) in enumerate(zip(data_to_plot, labels, colors)):
            ax.plot(model_output.index, data,
                   label=label,
                   color=color,
                   linewidth=2,
                   alpha=0.8)

    # Add total kd line if available and not in relative mode. Default is the dashed
    # black overlay; total_kd_style lets a caller restyle it (e.g. to the solid
    # multi_window trend used in the paper light-regime figure).
    if not relative and 'Phy_kd' in model_output.columns:
        tstyle = {'color': 'black', 'linestyle': '--', 'linewidth': 2, 'zorder': 10}
        if total_kd_style:
            tstyle.update(total_kd_style)
        ax.plot(model_output.index, model_output['Phy_kd'],
                label='Total k$_d$', **tstyle)

    # Format x-axis as dates if DatetimeIndex was used
    if isinstance(model_output.index, pd.DatetimeIndex):
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        if own_figure:
            fig.autofmt_xdate()

    # Formatting
    ax.set_xlabel('Time [d]', fontsize=12)

    if relative:
        ax.set_ylabel('Relative contribution [%]', fontsize=12)
        title_suffix = ' (Relative Contributions)'
        if stacked:
            ax.set_ylim(0, 100)
    else:
        ax.set_ylabel('Light attenuation coefficient [m$^{-1}$]', fontsize=12)
        title_suffix = ' (Absolute Contributions)'

    if add_title:
        plot_type = 'Stacked' if stacked else 'Overlaid'
        ax.set_title(f'{plot_type} Decomposition of k$_d${title_suffix}', fontsize=14)
    legend_opts = dict(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=10)
    if legend_kwargs is not None:
        legend_opts.update(legend_kwargs)
    ax.legend(**legend_opts)
    # Defer to the global rcParams['axes.grid'] switch (single lever for the paper)
    # instead of forcing the grid on.
    if plt.rcParams['axes.grid']:
        ax.grid(True, alpha=0.3)

    if own_figure:
        plt.tight_layout()

    # Save if requested
    if kwargs.get('save', False):
        savepath = kwargs.get('savepath', 'kd_contributions_stacked.png')
        fig.savefig(savepath, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {savepath}")

    return fig, ax


def plot_par_vs_depth(
        model: Any,
        depth_range: Tuple[float, float] = (0, 1.), #(0, 12.5),
        par_units: str = 'percent',
        photic_threshold: float = 0.01,
        show_minmax: bool = True,
        xscale: Optional[str] = 'linear',
        time_filter = None,
        seasons_definition: Optional[Dict[str, List[int]]] = None,
        figsize: Tuple[float, float] = (8, 6),
        save: bool = False,
        filename: str = 'PAR_vs_depth',
        figdir: Optional[str] = None,
        fnametimestamp: bool = True,
        **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot PAR intensity vs depth showing light extinction in the water column.

    Visualizes the vertical light regime by plotting PAR profiles for different
    light attenuation coefficients (kd) across seasons. Shows photic depth
    (depth where specified % of surface PAR remains).

    Args:
        model: Model simulation object with df attribute containing 'Phy_kd' and PAR data
        depth_range: Tuple of (min_depth, max_depth) in meters for the plot
        par_units: 'percent' for % of surface PAR, 'absolute' for actual PAR values
        photic_threshold: Fraction of surface PAR defining photic depth (default: 0.01 = 1%)
        show_minmax: Whether to show min and max kd curves in addition to seasonal means
        xscale: scale of the x axis (e.g. 'linear' or 'log')
        time_filter: Time filtering criterion (None, year string, date range tuple, or slice).
            - None: no filtering
            - str (year): e.g., "2023" filters to that year
            - tuple: (start, end) date strings
            - slice: pandas slice for advanced filtering
        seasons_definition: Dict mapping season names to month lists. If None, uses Silori et al. 2025:
            {'Winter': [12, 1, 2], 'Spring': [3, 4, 5], 'Summer': [6, 7, 8], 'Autumn': [9, 10, 11]}
        figsize: Figure size as (width, height)
        save: Whether to save the figure
        filename: Base filename for saved figure
        figdir: Directory to save figure (uses FIGURE_PATH if None)
        fnametimestamp: Whether to add timestamp to filename
        **kwargs: Additional plotting parameters

    Returns:
        Tuple of (figure, axes)

    Examples:
        >>> # Plot with PAR as % of surface (default)
        >>> fig, ax = plot_par_vs_depth(sim0, par_units='percent')
        >>>
        >>> # Plot with absolute PAR values
        >>> fig, ax = plot_par_vs_depth(sim0, par_units='absolute')
        >>>
        >>> # Show only seasonal means without min/max
        >>> fig, ax = plot_par_vs_depth(sim0, show_minmax=False)
        >>>
        >>> # Filter to year 2023
        >>> fig, ax = plot_par_vs_depth(sim0, time_filter="2023")
    """
    # Default season definition from Silori et al. 2025
    if seasons_definition is None:
        seasons_definition = {
            'Winter': [12, 1, 2],
            'Spring': [3, 4, 5],
            'Summer': [6, 7, 8],
            'Autumn': [9, 10, 11]
        }

    # Extract data from model
    if isinstance(model, pd.DataFrame):
        df = model
    else:
        df = model.df

    # Apply time filtering
    if time_filter is not None:
        df = _filter_by_time(df, time_filter)

    # Check required columns
    if 'Phy_kd' not in df.columns:
        raise ValueError("Model output must contain 'Phy_kd' column")

    # Get PAR data - try to find PAR column (might be in setup or derived)
    # PAR surface is typically stored in setup, but we need average values
    if hasattr(model, 'setup') and hasattr(model.setup, 'PAR_array'):
        # Get average PAR per season from setup
        par_surface_values = model.setup.PAR_array
        time_index = df.index
    else:
        # Fallback: assume constant PAR or extract from available data
        par_surface_values = np.ones(len(df)) * 100  # Default placeholder
        time_index = df.index

    # Create depth array
    depths = np.linspace(depth_range[0], depth_range[1], 1000)

    # Calculate seasonal kd statistics
    seasonal_data = {}
    for season_name, months in seasons_definition.items():
        # Filter by season
        mask = df.index.month.isin(months)
        season_df = df[mask]

        if len(season_df) > 0:
            kd_mean = season_df['Phy_kd'].mean()

            # Calculate average PAR for this season
            if hasattr(model, 'setup') and hasattr(model.setup, 'PAR_array'):
                season_par_indices = np.where(mask)[0]
                par_mean = np.mean(par_surface_values[season_par_indices])
            else:
                par_mean = 100  # Default

            # Calculate photic depth
            z_photic = -np.log(photic_threshold) / kd_mean if kd_mean > 0 else np.inf

            seasonal_data[season_name] = {
                'kd_mean': kd_mean,
                'par_mean': par_mean,
                'z_photic': z_photic
            }

    # Global min and max kd
    if show_minmax:
        kd_min = df['Phy_kd'].min()
        kd_max = df['Phy_kd'].max()
        z_photic_max = -np.log(photic_threshold) / kd_min if kd_min > 0 else np.inf
        z_photic_min = -np.log(photic_threshold) / kd_max if kd_max > 0 else np.inf

        # Use overall mean PAR for min/max curves
        if hasattr(model, 'setup') and hasattr(model.setup, 'PAR_array'):
            par_overall = np.mean(par_surface_values)
        else:
            par_overall = 100

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Define colors for seasons
    season_colors = {
        'Winter': '#3498db',
        'Spring': '#2ecc71',
        'Summer': '#f39c12',
        'Autumn': '#e74c3c'
    }

    # Plot seasonal profiles
    for season_name, data in seasonal_data.items():
        kd = data['kd_mean']
        par_0 = data['par_mean']
        z_photic = data['z_photic']

        # Calculate PAR profile: PAR(z) = PAR_0 * exp(-kd * z)
        if par_units == 'percent':
            par_profile = 100 * np.exp(-kd * depths)
            par_photic = 100 * photic_threshold
        else:  # absolute
            par_profile = par_0 * np.exp(-kd * depths)
            par_photic = par_0 * photic_threshold

        # Plot profile
        color = season_colors.get(season_name, 'gray')
        label = f'{season_name} (kd={kd:.3f} m$^{{-1}}$)'
        ax.plot(par_profile, depths, color=color, linewidth=2, label=label)

        # Mark photic depth
        if z_photic <= depth_range[1]:
            ax.plot(par_photic, z_photic, 'o', color=color, markersize=8,
                   markeredgecolor='black', markeredgewidth=1, zorder=10)
            ax.axhline(z_photic, color=color, linestyle='--', alpha=0.3, linewidth=1)

    # Plot min/max profiles
    if show_minmax:
        # Min kd (maximum light penetration)
        if par_units == 'percent':
            par_profile_max = 100 * np.exp(-kd_min * depths)
            par_photic_val = 100 * photic_threshold
        else:
            par_profile_max = par_overall * np.exp(-kd_min * depths)
            par_photic_val = par_overall * photic_threshold

        ax.plot(par_profile_max, depths, color='lightgreen', linewidth=1.5,
               linestyle=':', label=f'Min kd={kd_min:.3f} m$^{{-1}}$', alpha=0.7)
        if z_photic_max <= depth_range[1]:
            ax.plot(par_photic_val, z_photic_max, 's', color='lightgreen',
                   markersize=6, markeredgecolor='black', markeredgewidth=0.5, alpha=0.7)

        # Max kd (minimum light penetration)
        if par_units == 'percent':
            par_profile_min = 100 * np.exp(-kd_max * depths)
        else:
            par_profile_min = par_overall * np.exp(-kd_max * depths)

        ax.plot(par_profile_min, depths, color='darkred', linewidth=1.5,
               linestyle=':', label=f'Max kd={kd_max:.3f} m$^{{-1}}$', alpha=0.7)
        if z_photic_min <= depth_range[1]:
            ax.plot(par_photic_val, z_photic_min, 's', color='darkred',
                   markersize=6, markeredgecolor='black', markeredgewidth=0.5, alpha=0.7)

    # Format axes
    ax.invert_yaxis()  # Depth increases downward
    ax.set_ylabel('Depth [m]', fontsize=12)

    if par_units == 'percent':
        ax.set_xlabel('PAR [% of surface]', fontsize=12)
        # ax.set_xlim(0, 100)
    else:
        ax.set_xlabel('PAR [µmol photons m$^{-2}$ s$^{-1}$]', fontsize=12)
        # ax.set_xlim(left=0)

    ax.set_xscale(xscale)
    ax.set_ylim(depth_range[1], depth_range[0])
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='best', fontsize=9, framealpha=0.9)

    # Add photic depth reference line
    ax.axvline(100 * photic_threshold if par_units == 'percent' else par_overall * photic_threshold,
              color='gray', linestyle='-.', alpha=0.5, linewidth=1,
              label=f'Photic threshold ({photic_threshold*100:.1f}%)')

    title = f'PAR Extinction in Water Column (photic depth at {photic_threshold*100:.0f}% PAR)'
    ax.set_title(title, fontsize=13, fontweight='bold')

    plt.tight_layout()

    # Save if requested
    if save:
        if figdir is None:
            figdir = FIGURE_PATH
        save_figure(fig, filename=filename, figdir=figdir, add_timestamp=fnametimestamp)

    return fig, ax


def add_subplot_label(ax: plt.Axes,
                      index: int,
                      position: str = 'top_left',
                      fontsize: int = 10) -> None:
    """Add letter labels (a, b, c, etc.) to subplots."""
    positions = {
        'top_left': {'x': 0.02, 'y': 0.98, 'ha': 'left', 'va': 'top'},
        'top_right': {'x': 0.98, 'y': 0.98, 'ha': 'right', 'va': 'top'},
        'bottom_left': {'x': 0.02, 'y': 0.02, 'ha': 'left', 'va': 'bottom'},
        'bottom_right': {'x': 0.98, 'y': 0.02, 'ha': 'right', 'va': 'bottom'}
    }
    pos = positions.get(position, positions['top_left'])
    letter = chr(97 + index)
    ax.text(pos['x'], pos['y'], f'{letter}', transform=ax.transAxes,
            ha=pos['ha'], va=pos['va'], fontsize=fontsize)


def plot_optimization_evolution(df: pd.DataFrame,
                                costname: str = 'cost',
                                generationname: str = 'generation',
                                alpha: float = 0.2,
                                savefig: bool = False,
                                rawcost: bool = False,
                                name: Optional[str] = None,
                                ax: Optional[plt.Axes] = None) -> plt.Axes:
    """
    Plot optimization cost evolution.

    Args:
        df: Optimization results DataFrame
        costname: Name of cost column
        generationname: Name of generation column
        alpha: Scatter plot transparency
        savefig: Whether to save figure
        rawcost: Whether to plot raw cost (including bad scores)
        name: Optimization name for saving
        ax: Optional axes to plot on
    """
    if ax is None:
        _, ax = plt.subplots()

    cost_col = costname + rawcost * '_raw'

    # Plot all points
    ax.scatter(df[generationname], df[cost_col], alpha=alpha, edgecolors='none')

    # Highlight best point
    max_cost_idx = df[costname].idxmax()
    ax.plot(df[generationname][max_cost_idx],
            df[cost_col][max_cost_idx],
            'ro')
    ax.text(df[generationname][max_cost_idx],
            df[cost_col][max_cost_idx],
            f'{df[generationname][max_cost_idx]}',
            color='red')

    ax.set_xlabel('Optimization generation')
    ax.set_ylabel('Cost function (score)')

    if savefig and name:
        filename = f'{name}_opt_evol{rawcost * "_raw"}.png'
        save_figure(plt.gcf(), filename=filename)

    return ax


def plot_parameter_vs_cost(
        df: pd.DataFrame,
        param: str,
        ax: plt.Axes,
        costname: str = 'cost',
        alpha: float = 0.3,
        color: Optional[str] = None,
        label: Optional[str] = None,
        highlight_best: bool = True,
        show_reference: bool = True,
        show_labels: bool = True,
        filter_badlnl: Optional[float] = None,
        ylim_percent: Optional[float] = None
) -> plt.Axes:
    """
    Plot parameter values against optimization cost.

    Modern, flexible implementation for single or multi-optimization plots.

    Args:
        df: Optimization results DataFrame
        param: Parameter name to plot
        ax: Matplotlib axes to plot on
        costname: Name of cost column (default: 'cost')
        alpha: Scatter plot transparency
        color: Point color (None = use default blue)
        label: Label for legend (None = no label)
        highlight_best: Whether to highlight best point with star
        show_reference: Whether to show reference value line
        show_labels: Whether to set title and axis labels
        filter_badlnl: If provided, filter out rows where cost_raw == this value
        ylim_percent: If provided, zoom y-axis to top N% of cost range (default: None)

    Returns:
        Modified axes object
    """
    # Filter invalid scores if requested
    if filter_badlnl is not None:
        valid_mask = df[f'{costname}_raw'] != filter_badlnl
        df = df.loc[valid_mask]

    # Use provided color or default
    scatter_color = color if color else '#1f77b4'

    # Main scatter plot
    ax.scatter(df[param], df[costname], c=scatter_color, alpha=alpha,
               label=label, s=20, edgecolors='none')

    # Highlight best point
    if highlight_best:
        best_idx = df[costname].idxmax()
        best_value = df.loc[best_idx, param]
        best_cost = df.loc[best_idx, costname]
        ax.scatter(best_value, best_cost, marker='*', s=100,
                  c=scatter_color, edgecolors='black',
                  linewidths=0.5, zorder=4)

    # Get parameter info for reference and labels
    param_info = varinfos.ref_values.get(param, {})
    ref_value = param_info.get('reference_value')
    symbol = param_info.get('symbol', param)

    # Add reference line
    if show_reference and ref_value is not None:
        ax.axvline(ref_value, color='blue', linestyle='--',
                  alpha=0.7, linewidth=1.5)

    # Set labels if requested
    if show_labels:
        ax.set_title(fns.cleantext(symbol) if symbol else param, fontsize=10)
        ax.set_xlabel(param, fontsize=8)
        ax.set_ylabel(f'{costname} (log-likelihood)', fontsize=8)

    # Apply ylim zoom if requested
    if ylim_percent is not None:
        min_cost = df[costname].min()
        max_cost = df[costname].max()
        cost_range = max_cost - min_cost
        ylim_min = max_cost - (ylim_percent / 100) * cost_range
        ylim_max = max_cost + 0.1 * (ylim_percent / 100) * cost_range
        ax.set_ylim(ylim_min, ylim_max)

    return ax


def plot_optimization_summary(df: pd.DataFrame,
                              parameters: List[str],
                              costname: str = 'cost',
                              generationname: str = 'generation',
                              ncols: int = 3,
                              figsize: tuple = (10, 10),
                              alpha: float = 0.2,
                              savefig: bool = False,
                              rawcost: bool = False,
                              name: Optional[str] = None,
                              dateinname: bool = False,
                              ylim_percent: Optional[float] = 10.0) -> plt.Figure:
    """
    Create comprehensive optimization summary plots.

    Args:
        df: Optimization results DataFrame
        parameters: List of parameter names to plot
        costname: Name of cost column (default: 'cost')
        generationname: Name of generation column
        ncols: Number of columns in subplot grid
        figsize: Figure size
        alpha: Scatter plot transparency
        savefig: Whether to save figure
        rawcost: Whether to plot raw cost
        name: Optimization name for saving
        dateinname: Whether to include date in filename
        ylim_percent: Zoom y-axis to top N% of cost range (default: 20.0)

    Returns:
        Matplotlib figure
    """
    nrows = int(np.ceil(len(parameters) / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=figsize)
    if nrows == 1 and ncols == 1:
        axs = [axs]
    else:
        axs = axs.flatten()

    # Plot each parameter
    for i, param in enumerate(parameters):
        if i < len(parameters):
            plot_parameter_vs_cost(
                df, param,
                ax=axs[i],
                costname=costname,
                alpha=alpha,
                show_labels=True,
                highlight_best=True,
                show_reference=True,
                ylim_percent=ylim_percent
            )
        else:
            axs[i].set_visible(False)

    plt.tight_layout()

    if savefig:
        datestr = datetime.now().strftime('%Y%m%d_') if dateinname else ''
        filename = f'{datestr}{name}_pars_optim{rawcost * "_raw"}.png'
        save_figure(fig, filename=filename)

    return fig


def compare_optimizations(
        optimizations: List[Union[Any, str]],
        parameters: Optional[List[str]] = None,
        figsize: Optional[tuple] = None,
        alpha: float = 0.1,
        savefig: bool = False,
        name: str = 'optim_comparison',
        use_ref_order: bool = True
) -> plt.Figure:
    """
    Compare multiple optimizations with distribution and cost analysis.

    Creates 2-row grid with one column per parameter:
    - Row 1: Horizontal violin plots showing parameter distributions across optimizations
    - Row 2: Parameter vs cost scatter plots (all optimizations overlaid)

    Args:
        optimizations: List of Optimization objects or optimization names to load
        parameters: Parameters to compare (None = union of all optimized parameters)
        figsize: Figure size (None = auto-calculated based on number of parameters)
        alpha: Scatter plot transparency
        savefig: Whether to save figure
        name: Base filename for saving
        use_ref_order: If True, sort parameters by ref_values order; if False, alphabetically

    Returns:
        Matplotlib figure
    """
    from src.utils import optimization as optim

    # Load optimizations if given as strings or tuples (name, dir_suffix)
    opts = []
    for opt in optimizations:
        if isinstance(opt, str):
            opts.append(optim.Optimization.load_existing(opt))
        elif isinstance(opt, tuple):
            opts.append(optim.Optimization.load_existing(opt[0], dir_suffix=opt[1]))
        else:
            opts.append(opt)

    # Ensure all have processed results
    for opt in opts:
        if not hasattr(opt, 'df') or opt.df is None:
            opt.process_results()

    # Get union of all optimized parameters
    if parameters is None:
        all_params = set().union(*[set(opt.config['optimized_parameters']) for opt in opts])
        if use_ref_order:
            parameters = varinfos.sort_params_by_ref_order(all_params)
        else:
            parameters = sorted(all_params)

    n_params = len(parameters)
    n_opts = len(opts)

    # Dynamic ncols = n_params (one column per parameter)
    ncols = n_params
    nrows = 1 #2

    # Auto-calculate figsize if not provided
    if figsize is None:
        figsize = (max(16, 3.5 * ncols), 12)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)

    colors = DEFAULT_COLORS[:n_opts]

    # Track which subplot has all optimizations (for legend placement)
    legend_col = None

    for col_idx, param in enumerate(parameters):
        ax_violin = axes[0, col_idx]
        # ax_scatter = axes[1, col_idx]
        # ax_violin = axes[col_idx]

        # === ROW 1: HORIZONTAL Violin plot ===
        violin_data = []
        positions = []
        opt_names = []

        for i, opt in enumerate(opts):
            if param in opt.df.columns:
                # Filter out invalid scores (badlnl)
                badlnl = opt.config.get('badlnl', -100000.)
                valid_mask = opt.df['cost_raw'] != badlnl
                values = opt.df.loc[valid_mask, param].values

                if len(values) > 0:
                    violin_data.append(values)
                    positions.append(i)
                    opt_names.append(opt.name)

        if violin_data:
            # Create HORIZONTAL violin plot (vert=False)
            parts = ax_violin.violinplot(
                violin_data,
                positions=positions,
                vert=False,
                widths=0.7,
                showmeans=True,
                showextrema=True
            )

            # Color violin bodies
            for i, pc in enumerate(parts['bodies']):
                pc.set_facecolor(colors[positions[i]])
                pc.set_alpha(0.6)

            # Add jittered strip plot overlay
            for i, values in enumerate(violin_data):
                y_jitter = positions[i] + np.random.normal(0, 0.04, len(values))
                ax_violin.scatter(values, y_jitter, alpha=alpha, s=10,
                                c=colors[positions[i]], zorder=3, edgecolors='none')

            # Add best value stars with annotations
            for i, opt_idx in enumerate(positions):
                opt = opts[opt_idx]
                if param in opt.df.columns:
                    best_idx = opt.df['cost'].idxmax()
                    best_value = opt.df.loc[best_idx, param]
                    ax_violin.scatter(best_value, opt_idx, marker='*', s=150,
                                    c=colors[opt_idx], edgecolors='black',
                                    linewidths=0.8, zorder=5)
                    # Add value annotation above the star (smart formatting)
                    # Use scientific notation if value is very small or very large
                    if abs(best_value) < 0.01 or abs(best_value) > 1000:
                        value_text = f'{best_value:.2e}'
                    else:
                        value_text = f'{best_value:.2f}'
                    ax_violin.text(best_value, opt_idx - 0.2, value_text,
                                 fontsize=9, ha='center', va='top',
                                 color=colors[opt_idx], weight='bold',
                                 bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                          edgecolor='none', alpha=0.7))

        # Get parameter info for labels and reference
        param_info = varinfos.ref_values.get(param, {})
        ref_value = param_info.get('reference_value')
        symbol = param_info.get('symbol', param)

        # Add reference line (vertical for horizontal violins)
        if ref_value is not None:
            ax_violin.axvline(ref_value, color='blue', linestyle='--',
                            alpha=0.7, linewidth=1.5, label='Reference')

        ax_violin.set_title(fns.cleantext(symbol) if symbol else param, fontsize=10)
        ax_violin.set_yticks(range(n_opts))
        ax_violin.set_yticklabels([opt.name for opt in opts], fontsize=7)
        ax_violin.tick_params(labelsize=7)
        ax_violin.grid(axis='x', alpha=0.3)
        ax_violin.invert_yaxis()  # Invert to have first optimization at top

        # Hide x tick labels for all violins (redundant with scatter plots below)
        ax_violin.tick_params(labelbottom=False)

        # Hide y tick labels except for first column
        if col_idx > 0:
            ax_violin.tick_params(labelleft=False)

        # === ROW 2: Scatter plot (param vs cost) ===
        # for i, opt in enumerate(opts):
        #     if param in opt.df.columns:
        #         badlnl = opt.config.get('badlnl', -100000.)
        #         plot_parameter_vs_cost(
        #             opt.df, param,
        #             ax=ax_scatter,
        #             costname='cost',
        #             alpha=alpha,
        #             color=colors[i],
        #             label=opt.name,
        #             highlight_best=True,
        #             show_reference=(i == 0),  # Only first opt shows reference
        #             show_labels=False,  # We set labels manually below
        #             filter_badlnl=badlnl
        #         )
        #
        # # Set labels manually (shared across all optimizations)
        # param_info = varinfos.ref_values.get(param, {})
        # symbol = param_info.get('symbol', param)
        # ax_scatter.set_xlabel(fns.cleantext(symbol) if symbol else param, fontsize=8)
        # ax_scatter.set_ylabel('Cost (log-likelihood)', fontsize=8)
        # ax_scatter.tick_params(labelsize=7)
        # ax_scatter.grid(alpha=0.3)
        #
        # # Hide y tick labels except for first column
        # if col_idx > 0:
        #     ax_scatter.set_ylabel('')
        #     ax_scatter.tick_params(labelleft=False)
        #
        # # Track if this column has all optimizations (for legend)
        # if legend_col is None and len(positions) == n_opts:
        #     legend_col = col_idx
        #
        # # Add legend only to first column that has ALL optimizations
        # if col_idx == legend_col:
        #     ax_scatter.legend(fontsize=6, loc='best')

    plt.tight_layout()

    if savefig:
        save_figure(fig, filename=f'{name}_comparison.png')

    return fig


def compare_optim_values(
        optimizations: List[Union[Any, str, tuple]],
        parameters: Optional[List[str]] = None,
        ncols: int = 5,
        percentiles: Optional[Tuple[float, float]] = (25, 75),
        show_values: bool = True,
        marker_size: int = 120,
        figsize: Optional[tuple] = None,
        savefig: bool = False,
        name: str = 'optim_values',
        use_ref_order: bool = True
) -> plt.Figure:
    """
    Compare best parameter values across optimizations with optional uncertainty bars.

    Simple scatter plot: one subplot per parameter, one point per optimization.
    Y-axis shows optimizations (first in list at top), X-axis shows parameter value.

    Args:
        optimizations: List of Optimization objects, names, or (name, dir_suffix) tuples
        parameters: Parameters to compare (None = union of all optimized parameters)
        ncols: Maximum columns in subplot grid
        percentiles: Tuple (low, high) for uncertainty bars, None to disable
        show_values: Whether to annotate points with numeric values
        marker_size: Size of scatter markers
        figsize: Figure size (None = auto-calculated)
        savefig: Whether to save figure
        name: Base filename for saving
        use_ref_order: If True, sort parameters by ref_values order; if False, alphabetically

    Returns:
        Matplotlib figure
    """
    from src.utils import optimization as optim

    # Load optimizations
    opts = []
    for o in optimizations:
        if isinstance(o, str):
            opts.append(optim.Optimization.load_existing(o))
        elif isinstance(o, tuple):
            opts.append(optim.Optimization.load_existing(o[0], dir_suffix=o[1]))
        else:
            opts.append(o)

    for o in opts:
        if not hasattr(o, 'df') or o.df is None:
            o.process_results()

    # Get parameters
    if parameters is None:
        all_params = set().union(*[set(o.config['optimized_parameters']) for o in opts])
        if use_ref_order:
            parameters = varinfos.sort_params_by_ref_order(all_params)
        else:
            parameters = sorted(all_params)

    n_params, n_opts = len(parameters), len(opts)
    ncols = min(ncols, n_params)
    nrows = int(np.ceil(n_params / ncols))

    if figsize is None:
        figsize = (3.5 * ncols, max(4, 0.6 * n_opts) * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    colors = DEFAULT_COLORS[:n_opts]
    y_positions = np.arange(n_opts)

    for idx, param in enumerate(parameters):
        row, col = divmod(idx, ncols)
        ax = axes[row, col]

        for i, o in enumerate(opts):
            if param not in o.df.columns:
                continue

            badlnl = o.config.get('badlnl', -100000.)
            valid = o.df[o.df['cost_raw'] != badlnl]
            if valid.empty:
                continue

            # Best value
            best_idx = valid['cost'].idxmax()
            best_val = valid.loc[best_idx, param]

            # Plot point
            ax.scatter(best_val, i, s=marker_size, c=colors[i], edgecolors='black',
                       linewidths=0.8, zorder=5, label=o.name if idx == 0 else None)

            # Uncertainty bar
            if percentiles and len(valid) > 1:
                p_low, p_high = np.percentile(valid[param], percentiles)
                ax.hlines(i, p_low, p_high, colors=colors[i], linewidths=2, alpha=0.5, zorder=3)

            # Value annotation
            if show_values:
                txt = f'{best_val:.2e}' if abs(best_val) < 0.01 or abs(best_val) > 1000 else f'{best_val:.2g}'
                ax.annotate(txt, (best_val, i), xytext=(8, 0), textcoords='offset points',
                            fontsize=7, va='center', color=colors[i], zorder=10,
                            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))

        # Styling
        param_info = varinfos.ref_values.get(param, {})
        symbol = param_info.get('symbol', param)
        ax.set_title(fns.cleantext(symbol) if symbol else param, fontsize=10)
        ax.set_yticks(y_positions)
        ax.set_yticklabels([o.name for o in opts] if col == 0 else [], fontsize=8)
        ax.set_ylim(-0.5, n_opts - 0.5)
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)
        ax.tick_params(labelsize=8)

    # Hide unused subplots
    for idx in range(n_params, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row, col].set_visible(False)

    # Single legend at bottom
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc='lower center', ncol=min(n_opts, 6),
                   fontsize=8, bbox_to_anchor=(0.5, 0), frameon=True)

    plt.tight_layout(rect=[0, 0.06, 1, 0.98], h_pad=2.5)

    if savefig:
        save_figure(fig, filename=f'{name}.png')

    return fig


def create_parameter_table(df: pd.DataFrame,
                           parameters: List[tuple],
                           costname: str = 'cost',
                           fontsize: int = 12,
                           figsize: tuple = (12, 7)) -> plt.Figure:
    """
    Create parameter summary table figure.

    Args:
        df: Optimization results DataFrame
        parameters: List of (param_name, min_val, max_val) tuples
        costname: Name of cost column
        fontsize: Table font size
        figsize: Figure size
    """
    # Prepare pars_info
    pars_info = []
    for key, min_val, max_val in parameters:
        param_info = varinfos.ref_values[key]
        ref_val = param_info.get('reference_value')
        pars_info.append({
            "Name": param_info['complete_name'],
            "Symbol": fns.cleantext(param_info.get('symbol', '')),
            "Ref. value": f"{ref_val:.3g}" if ref_val is not None else "N/A",
            "Units": fns.cleantext(param_info.get('units', '')),
            "Min. value": f"{min_val:.3g}",
            "Max. value": f"{max_val:.3g}",
            "Calibrated value": f"{df[key][df[costname].idxmax()]:.3g}"
        })

    # Create table figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
    ax.axis('tight')
    ax.axis('off')

    tab_df = pd.DataFrame(pars_info)
    table = ax.table(
        cellText=tab_df.values,
        colLabels=tab_df.columns,
        loc='center',
        cellLoc='right',
        rowLoc='right',
        colLoc='center',
        bbox=(0.0, 0.0, 1., 1.)
    )

    # Style table
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)
    table.auto_set_column_width(col=list(range(len(tab_df.columns))))

    # Color header
    for i in range(len(tab_df.columns)):
        table[(0, i)].set_facecolor('#d6d6d6')

    return fig


def plot_monthly_budget_comparison(
        models: List[Any],
        var_name: str = 'Phy_source_PP.C',
        time_filter = None,
        figsize: Tuple[float, float] = (12, 6),
        colors: Optional[List[str]] = None,
        ratio_color: str = 'black',
        save: bool = False,
        filename: str = 'monthly_budget_comparison',
        figdir: Optional[str] = None,
        fnametimestamp: bool = True,
        **kwargs
) -> Tuple[plt.Figure, plt.Axes, plt.Axes]:
    """
    Plot monthly integrated budget comparison between models with dual y-axes.

    Creates a visualization with:
    - Left axis: Grouped bars showing monthly integrated values for each model
    - Right axis: Line plot showing ratio between models (model[1]/model[0])

    Args:
        models: List of model objects (with .df and .name attributes)
        var_name: Variable name to integrate monthly (default: 'Phy_source_PP.C')
        time_filter: Time filtering criterion (None, year string, date range tuple, or slice)
        figsize: Figure size (width, height)
        colors: List of colors for models (None = use DEFAULT_COLORS)
        ratio_color: Color for ratio line (default: 'black')
        save: Whether to save the figure
        filename: Base filename for saving
        figdir: Directory to save figure (uses FIGURE_PATH if None)
        fnametimestamp: Whether to add timestamp to filename
        **kwargs: Additional plotting parameters

    Returns:
        Tuple of (figure, left_axis, right_axis)

    Example:
        >>> fig, ax1, ax2 = plot_monthly_budget_comparison(
        ...     [sim_ref, sim_test],
        ...     var_name='Phy_source_PP.C',
        ...     time_filter='2023'
        ... )
    """
    # Prepare data
    model_data_list, _, _, model_names = prepare_model_obs_data(
        models, observations=None, mean_window_days=None,
        variables_to_plot=[var_name], time_filter=time_filter
    )

    n_models = len(model_data_list)
    if n_models < 1:
        raise ValueError("At least one model is required")

    # Get colors
    if colors is None:
        colors = DEFAULT_COLORS[:n_models]

    # Calculate monthly budgets for each model
    monthly_data = {name: [] for name in model_names}
    month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    for df, name in zip(model_data_list, model_names):
        # Calculate timestep
        time_diff = df.index.to_series().diff()
        dt_days = time_diff.dt.total_seconds().median() / 86400.0

        # Integrate by month
        for month in range(1, 13):
            df_month = df[df.index.month == month]
            if len(df_month) > 0 and var_name in df_month.columns:
                monthly_budget = np.sum(df_month[var_name]) * dt_days
            else:
                monthly_budget = 0.0
            monthly_data[name].append(monthly_budget)

    # Calculate ratio if multiple models
    ratio = None
    if n_models >= 2:
        ref_vals = np.array(monthly_data[model_names[0]])
        comp_vals = np.array(monthly_data[model_names[1]])
        ratio = np.divide(comp_vals, ref_vals, where=ref_vals!=0, out=np.full_like(comp_vals, np.nan))

    # Get variable info from varinfos
    var_info = varinfos.doutput.get(var_name.lstrip('m'), {})
    clean_name = var_info.get('cleanname', None)
    if clean_name is None:
        clean_name = var_name.replace('_', r'\_')
    units = var_info.get('munits' if var_name.startswith('m') else 'units', '')

    # Create figure with dual axes
    fig, ax1 = plt.subplots(figsize=figsize)
    ax2 = ax1.twinx()

    # Plot grouped bars on left axis
    x = np.arange(len(month_names))
    width = 0.8 / n_models

    for i, name in enumerate(model_names):
        offset = (i - n_models/2 + 0.5) * width
        ax1.bar(x + offset, monthly_data[name], width,
               label=name, color=colors[i], alpha=0.8)

    # Plot ratio on right axis
    if ratio is not None:
        ax2.plot(x, ratio, color=ratio_color, marker='o', linewidth=2,
                markersize=6, label='Ratio', zorder=10)
        ax2.axhline(1.0, color='gray', linestyle='--', alpha=0.5, linewidth=1)

    # Format axes
    ax1.set_xlabel('Month', fontsize=11)
    ax1.set_ylabel(f'Monthly {fns.cleantext(clean_name)}\n[{fns.cleantext(units)}]', fontsize=11, color='black')
    ax1.set_xticks(x)
    ax1.set_xticklabels(month_names)
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.grid(axis='y', alpha=0.3, zorder=0)

    if ratio is not None:
        ax2.set_ylabel(f'Ratio ({model_names[1]} / {model_names[0]})', fontsize=11, color=ratio_color)
        ax2.tick_params(axis='y', labelcolor=ratio_color)
        ax2.set_ylim(0, 1.2)

    # Legends
    ax1.legend(loc='upper left', fontsize=9)
    if ratio is not None:
        ax2.legend(loc='upper right', fontsize=9)

    plt.tight_layout()

    # Save if requested
    if save:
        if figdir is None:
            figdir = FIGURE_PATH
        save_figure(fig, filename=filename, figdir=figdir, add_timestamp=fnametimestamp)

    return fig, ax1, ax2


def save_figure(fig: plt.Figure,
                filename: Optional[str] = None,
                variable_set: Optional[PlottedVariablesSet] = None,
                add_timestamp: bool = True,
                figdir: str = FIGURE_PATH,
                dpi: int = DEFAULT_FIGURE_DPI,
                **savefig_kwargs) -> str:
    """
    Centralized figure saving with intelligent naming from PlottedVariablesSet.

    Args:
        fig: Figure to save
        filename: Explicit filename (takes priority over variable_set.name)
        variable_set: PlottedVariablesSet to extract name from
        add_timestamp: Whether to add timestamp to filename
        figdir: Directory to save figures in
        dpi: Figure DPI
        **savefig_kwargs: Additional arguments for plt.savefig

    Returns:
        Full path of saved file
    """
    # Determine filename
    if filename is not None:
        base_name = filename
    elif variable_set and hasattr(variable_set, 'name'):
        base_name = variable_set.name
    else:
        base_name = "figure"

    # Add timestamp if requested
    if add_timestamp:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S_')
        full_name = f"{timestamp}{base_name}.png"
    else:
        full_name = f"{base_name}.png" if not base_name.endswith('.png') else base_name

    # Create full path
    full_path = os.path.join(figdir, full_name)

    # Set default savefig parameters
    default_kwargs = {
        'dpi': dpi,
        'bbox_inches': 'tight',
        'facecolor': 'white',
        'edgecolor': 'none'
    }

    # Merge with user-provided kwargs
    save_kwargs = {**default_kwargs, **savefig_kwargs}

    # Save figure
    fig.savefig(full_path, **save_kwargs)

    return full_path


# ============================================================================
# ================= CARBON-FLUX NETWORK (SANKEY) =============================
# ============================================================================

# Default layout for the organic-C network produced by
# simulation_manager.carbon_flux_links. Positions are on a 0-1 canvas (x flows
# left->right along the trophic chain; external sinks are on the right column).
CARBON_SANKEY_POS = {
    'DIC_in': (0.03, 0.50), 'Phy': (0.16, 0.74),
    'DOCS': (0.35, 0.90), 'DOCL': (0.35, 0.62),
    'TEPC': (0.50, 0.78), 'DetS': (0.42, 0.32), 'DetL': (0.58, 0.15),
    'BacF': (0.71, 0.90), 'BacA': (0.71, 0.42), 'HF': (0.83, 0.66), 'Cil': (0.90, 0.56),
    'CO2': (0.95, 0.86), 'Leak': (0.95, 0.36), 'Export': (0.95, 0.13),
}
CARBON_SANKEY_CAT = {
    'Phy': 'phy', 'DOCS': 'dom', 'DOCL': 'dom', 'TEPC': 'tep', 'DetS': 'det',
    'DetL': 'det', 'BacF': 'bac', 'BacA': 'bac', 'HF': 'zoo', 'Cil': 'zoo',
    'DIC_in': 'ext', 'CO2': 'ext', 'Export': 'ext', 'Leak': 'ext',
}
CARBON_SANKEY_COL = {
    'phy': '#2e8b57', 'dom': '#4aa3df', 'tep': '#e08a1e', 'det': '#8c6d3f',
    'bac': '#d1495b', 'zoo': '#8e44ad', 'ext': '#8a8f98',
}
_CARBON_SANKEY_LAB = {'DIC_in': 'DIC', 'CO2': 'CO$_2$', 'Export': 'Export', 'Leak': 'Leak'}
_CARBON_SANKEY_LEGEND = {'phy': 'Phytoplankton', 'dom': 'DOC', 'tep': 'TEP',
                         'det': 'Detritus', 'bac': 'Bacteria', 'zoo': 'Protozoa',
                         'ext': 'External'}


_SANKEY_NODE_W = 0.026


def _sankey_ribbon(p0, p1, w, ctrl, n=28):
    """Constant-width-`w` ribbon polygon along a cubic Bezier p0->p1 whose control points
    sit a horizontal distance `ctrl` inside each endpoint, so the ribbon leaves and enters
    horizontally (perpendicular to the vertical node bars) regardless of the vertical drop."""
    x0, y0 = p0
    x1, y1 = p1
    c0 = (x0 + ctrl, y0)
    c1 = (x1 - ctrl, y1)
    t = np.linspace(0, 1, n)[:, None]
    B = ((1 - t) ** 3 * np.array(p0) + 3 * (1 - t) ** 2 * t * np.array(c0)
         + 3 * (1 - t) * t ** 2 * np.array(c1) + t ** 3 * np.array(p1))
    d = np.gradient(B, axis=0)
    nrm = np.stack([-d[:, 1], d[:, 0]], 1)
    nrm /= (np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12)
    return np.vstack([B + nrm * w / 2, (B - nrm * w / 2)[::-1]])


def _draw_carbon_sankey(ax, links, scale, pos, cat, col, lab, min_flux, title, subtitle,
                        ref_flux, mode='stacked', min_node_h=0.03, ctrl_min=0.11,
                        alpha=0.55):
    from matplotlib.path import Path as _MPath
    from matplotlib.patches import PathPatch, FancyBboxPatch, Rectangle
    ax.set_xlim(-0.04, 1.18)
    ax.set_ylim(0, 1)
    ax.axis('off')
    if title:
        ax.text(0.5, 1.05, title, ha='center', va='bottom', fontsize=12,
                fontweight='bold', transform=ax.transAxes)
    if subtitle:
        ax.text(0.5, 1.005, subtitle, ha='center', va='bottom', fontsize=7.5,
                color='0.25', transform=ax.transAxes)

    links = {k: v for k, v in links.items()
             if v >= min_flux and k[0] in pos and k[1] in pos}
    nodes = set(a for a, b in links) | set(b for a, b in links)
    out_tot = {n: sum(v for (a, b), v in links.items() if a == n) for n in nodes}
    in_tot = {n: sum(v for (a, b), v in links.items() if b == n) for n in nodes}
    hgt = {n: max(out_tot[n] * scale, in_tot[n] * scale, min_node_h) for n in nodes}

    # Anchor y of each link on the source (right edge) and target (left edge). In 'stacked'
    # mode the links are juxtaposed along each node bar (ordered by the other end's height,
    # to limit crossings) so each flux occupies a share of the bar proportional to its
    # value -- proportions read at a glance. In 'centered' mode every link attaches at the
    # node centre (ribbons overlap with transparency).
    src_anchor, tgt_anchor = {}, {}
    if mode == 'stacked':
        for n in nodes:
            cur = pos[n][1] + out_tot[n] * scale / 2
            for b, v in sorted([(b, v) for (a, b), v in links.items() if a == n],
                               key=lambda tv: -pos[tv[0]][1]):
                src_anchor[(n, b)] = cur - v * scale / 2
                cur -= v * scale
            cur = pos[n][1] + in_tot[n] * scale / 2
            for a, v in sorted([(a, v) for (a, b), v in links.items() if b == n],
                               key=lambda tv: -pos[tv[0]][1]):
                tgt_anchor[(a, n)] = cur - v * scale / 2
                cur -= v * scale
    else:
        for (a, b) in links:
            src_anchor[(a, b)] = pos[a][1]
            tgt_anchor[(a, b)] = pos[b][1]

    for (a, b), v in sorted(links.items(), key=lambda kv: -kv[1]):
        p0 = (pos[a][0] + _SANKEY_NODE_W / 2, src_anchor[(a, b)])
        p1 = (pos[b][0] - _SANKEY_NODE_W / 2, tgt_anchor[(a, b)])
        ctrl = max(0.4 * abs(p1[0] - p0[0]), ctrl_min)
        ax.add_patch(PathPatch(_MPath(_sankey_ribbon(p0, p1, v * scale, ctrl)),
                     facecolor=col[cat[a]], edgecolor='none', alpha=alpha, zorder=1))
    for n in nodes:
        x, y = pos[n]
        h = hgt[n]
        ax.add_patch(FancyBboxPatch((x - _SANKEY_NODE_W / 2, y - h / 2), _SANKEY_NODE_W, h,
                     boxstyle='round,pad=0.002,rounding_size=0.008', facecolor=col[cat[n]],
                     edgecolor='black', linewidth=0.6, zorder=3))
        if x > 0.9:  # right-hand external sinks: label outside the box (avoids clipping)
            ax.text(x + _SANKEY_NODE_W / 2 + 0.012, y, lab.get(n, n), ha='left',
                    va='center', fontsize=7.5, fontweight='bold', color=col[cat[n]],
                    zorder=4)
        else:
            ax.text(x, y, lab.get(n, n), ha='center', va='center', fontsize=6.2,
                    fontweight='bold', color='white', zorder=4,
                    rotation=90 if h < 0.05 else 0)
    if ref_flux:
        ax.add_patch(Rectangle((0.0, 0.02), 0.09, ref_flux * scale, facecolor='0.4',
                     edgecolor='none', alpha=0.5))
        ax.text(0.10, 0.02 + ref_flux * scale / 2, f'= {ref_flux:.0f}', va='center',
                fontsize=6.5, color='0.3')


def plot_carbon_sankey(links_list, names, mode='stacked', scale=None, min_flux=400,
                       ref_flux=10000.0, subtitles=None, pos=None, cat=None, col=None,
                       figsize_per_panel=None, suptitle=None, orient='vertical'):
    """Draw one carbon-flux network per simulation as side-by-side Sankey diagrams.

    Ribbon width is proportional to the annual C flux and the scale is shared across
    panels, so REF vs NO-TEP are directly comparable. Feed it the link dicts from
    simulation_manager.carbon_flux_links.

    Args:
        links_list: list of {(src, tgt): flux} dicts (one per simulation).
        names: panel titles (e.g. ['REF', 'NO-TEP']).
        mode: 'stacked' (default) juxtaposes the links along each node bar, whose height is
            proportional to the node throughput, so proportions read directly; 'centered'
            attaches every link at the node centre and overlaps them with transparency.
        scale: canvas-units per flux unit; default sizes the busiest node to ~0.17
            (stacked) or the largest link to ~0.11 (centered).
        min_flux: links below this (same units as the fluxes) are not drawn (declutter).
        ref_flux: value of the width-reference bar shown in each panel (None to hide).
        subtitles: optional {name: str} of per-panel annotation (e.g. PP / export totals).
        pos, cat, col: layout/category/colour overrides (default CARBON_SANKEY_*).
        figsize_per_panel, suptitle: figure sizing and overall title.
        orient: 'vertical' (default) stacks the panels top-to-bottom, each spanning the
            full figure width (publication page width, REF above / NO-TEP below);
            'horizontal' places them side by side.

    Returns:
        (fig, axes).
    """
    import matplotlib.patches as mpatches
    pos = pos or CARBON_SANKEY_POS
    cat = cat or CARBON_SANKEY_CAT
    col = col or CARBON_SANKEY_COL
    subtitles = subtitles or {}
    if scale is None:
        if mode == 'stacked':
            def _thru(L):
                nn = set(a for a, b in L if L[(a, b)] >= min_flux) | \
                     set(b for a, b in L if L[(a, b)] >= min_flux)
                return max((max(sum(v for (a, b), v in L.items() if a == n and v >= min_flux),
                                sum(v for (a, b), v in L.items() if b == n and v >= min_flux))
                            for n in nn), default=1.0)
            scale = 0.17 / max(_thru(L) for L in links_list)
        else:
            scale = 0.11 / max(max(L.values()) for L in links_list)
    n = len(links_list)
    if orient == 'vertical':
        fpp = figsize_per_panel or (11.0, 4.8)
        fig, axes = plt.subplots(n, 1, figsize=(fpp[0], fpp[1] * n))
    else:
        fpp = figsize_per_panel or (7.8, 7.4)
        fig, axes = plt.subplots(1, n, figsize=(fpp[0] * n, fpp[1]))
    axes = np.atleast_1d(axes)
    for ax, L, nm in zip(axes, links_list, names):
        _draw_carbon_sankey(ax, L, scale, pos, cat, col, _CARBON_SANKEY_LAB, min_flux,
                            nm, subtitles.get(nm), ref_flux, mode=mode)
    handles = [mpatches.Patch(color=col[k], label=v)
               for k, v in _CARBON_SANKEY_LEGEND.items()]
    fig.legend(handles=handles, loc='lower center', ncol=len(handles), fontsize=8,
               frameon=False)
    if suptitle is None:
        suptitle = (r'Annual carbon flux network  [mmol C m$^{-2}$ yr$^{-1}$]  '
                    r'(width $\propto$ flux)')
    fig.suptitle(suptitle, fontsize=12, y=0.99)
    fig.tight_layout(rect=[0, 0.045, 1, 0.95])
    return fig, axes


# ============================================================================
# ============ resuspension_rate SENSITIVITY (REF vs NO-TEP) =================
# ============================================================================

def plot_resusp_sensitivity(ref_df, notep_df, calib=1.323e6, figsize=(15.5, 8.4),
                            suptitle=None):
    """Six-panel resuspension_rate sensitivity of the C cycle, REF (TEP on) vs NO-TEP.

    The §5 knock-out robustness / mechanism figure. Feed two tidy DataFrames (one per
    regime), each with columns: resusp_rate, PP, e_ratio, Phy_limI, Phy_limNUT, N_tot,
    N_export (rows = one simulation per swept resuspension_rate). Produced by the sweep
    in _private/scripts (run_sensitivity over Macroflocs+resuspension_rate).

    The story it draws:
      (a) PP and (b) e-ratio vs resuspension_rate -- NO-TEP responds smoothly/monotonically
          while REF is strongly non-linear (PP peaks at the calibrated rate, export flips
          sign): the TEP-mineral coupling *amplifies* the C-cycle's sensitivity to the
          resuspension environment rather than buffering it.
      (c) REF light (limI) vs nutrient (limNUT~Si) limitation -- two competing controls:
          nutrient-limited below the calibrated rate (organic Si/N exported and not
          returned), light-limited above it (turbid), optimum at calibration.
      (d) limNUT REF vs NO-TEP -- the nutrient-limitation collapse at low resuspension
          happens *only* with TEP.
      (e) total nutrient inventory (N_tot) and (f) net nutrient export vs resuspension_rate
          -- direct evidence: with TEP the nutrient inventory is tightly slaved to the
          resuspension environment (huge swing, export at low rate / accumulation at high
          rate), whereas NO-TEP stays buffered.

    Args:
        ref_df, notep_df: per-regime DataFrames (sorted on resusp_rate internally).
        calib: calibrated resuspension_rate, marked with a dotted line.
        figsize, suptitle: figure sizing / overriding title.

    Returns:
        (fig, axs) with axs a (2, 3) array.
    """
    cref, cnot = '#2e8b57', '#d1495b'
    R = ref_df.sort_values('resusp_rate')
    N = notep_df.sort_values('resusp_rate')
    x = 'resusp_rate'
    fig, axs = plt.subplots(2, 3, figsize=figsize)

    ax = axs[0, 0]
    ax.plot(R[x], R.PP, 'o-', color=cref, label='REF (TEP on)')
    ax.plot(N[x], N.PP, 's-', color=cnot, label='NO-TEP')
    ax.set_ylabel(r'PP  [mmol C m$^{-2}$ yr$^{-1}$]')
    ax.set_title('(a) Primary production'); ax.legend(fontsize=8)

    ax = axs[0, 1]
    ax.plot(R[x], R.e_ratio, 'o-', color=cref, label='REF (TEP on)')
    ax.plot(N[x], N.e_ratio, 's-', color=cnot, label='NO-TEP')
    ax.axhline(0, color='0.7', lw=0.8)
    ax.set_ylabel('e-ratio  Export/NPP  [%]')
    ax.set_title('(b) Export efficiency'); ax.legend(fontsize=8)

    ax = axs[1, 0]
    ax.plot(R[x], R.Phy_limI, 'o-', color='#e0a11e', label='light (limI)')
    ax.plot(R[x], R.Phy_limNUT, 'o-', color='#3b6fb0',
            label=r'nutrients (limNUT $\approx$ Si)')
    ax.set_ylabel('limitation factor  [0-1]')
    ax.set_title('(c) REF: two competing controls'); ax.legend(fontsize=8)
    ax.annotate('nutrient-limited\n(Si exported)', (R[x].iloc[0], R.Phy_limNUT.iloc[0]),
                textcoords='offset points', xytext=(6, 22), fontsize=7.5, color='#3b6fb0')
    ax.annotate('light-limited\n(turbid)', (R[x].iloc[-1], R.Phy_limI.iloc[-1]),
                textcoords='offset points', xytext=(-30, 18), fontsize=7.5, color='#e0a11e')

    ax = axs[1, 1]
    ax.plot(R[x], R.Phy_limNUT, 'o-', color=cref, label='REF (TEP on)')
    ax.plot(N[x], N.Phy_limNUT, 's-', color=cnot, label='NO-TEP')
    ax.set_ylim(0, 0.8)
    ax.set_ylabel('nutrient limitation limNUT  [0-1]')
    ax.set_title('(d) Nutrient limitation collapses only with TEP'); ax.legend(fontsize=8)

    ax = axs[0, 2]
    ax.plot(R[x], R.N_tot, 'o-', color=cref, label='REF (TEP on)')
    ax.plot(N[x], N.N_tot, 's-', color=cnot, label='NO-TEP')
    ax.set_ylabel(r'total N inventory  [mmol N m$^{-3}$]')
    ax.set_title('(e) Nutrient inventory slaved to resuspension'); ax.legend(fontsize=8)

    ax = axs[1, 2]
    ax.plot(R[x], R.N_export, 'o-', color=cref, label='REF (TEP on)')
    ax.plot(N[x], N.N_export, 's-', color=cnot, label='NO-TEP')
    ax.axhline(0, color='0.7', lw=0.8)
    ax.set_ylabel(r'net N export  [mmol N m$^{-2}$ yr$^{-1}$]')
    ax.set_title('(f) Net nutrient export (>0 = lost to bed)'); ax.legend(fontsize=8)
    ax.annotate('exported', (R[x].iloc[0], R.N_export.iloc[0]),
                textcoords='offset points', xytext=(6, -14), fontsize=7.5, color=cref)
    ax.annotate('returned', (R[x].iloc[-1], R.N_export.iloc[-1]),
                textcoords='offset points', xytext=(-40, 10), fontsize=7.5, color=cref)

    for a in axs.flat:
        a.axvline(calib, ls=':', color='0.5')
        a.set_xlabel('resuspension_rate')
        a.text(calib, a.get_ylim()[0], ' calib', fontsize=7, color='0.5', va='bottom')
    if suptitle is None:
        suptitle = ('TEP-mineral coupling amplifies the C-cycle sensitivity to '
                    'resuspension\nNO-TEP: light-controlled (smooth)   |   '
                    'REF: nutrient (low rate) <-> light (high rate)')
    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    return fig, axs


def plot_tep_double_advantage(ref, notep, year=2023, days=10, figsize=(12, 10.5)):
    """Seasonal test of TEP's *double advantage* (REF vs NO-TEP at the calibrated rate).

    Phytoplankton exude TEP, which (1) glues mineral flocs -> lowers SPM/turbidity ->
    improves the light climate, and (2) fuels the microbial loop ('farming' bacteria) ->
    nutrient recycling. This figure tests how each advantage plays out over the season:

      (a) PP -- REF outproduces NO-TEP, mostly in summer.
      (b) light limitation limI -- advantage 1: clearly higher (less limiting) with TEP.
      (c) nutrient limitation lim_Si / lim_P (Droop quota) -- essentially unchanged: the
          quota formulation buffers it, so nutrient limitation is *not* the proximate PP
          control here.
      (d) dissolved DSi / DIP -- collapse post-bloom without TEP (weaker recycling).
      (e) bacterial biomass -- advantage 2: microbial loop far stronger with TEP.
      (f) nutrient regeneration (DSi+DIP remineralisation) -- ~2x higher with TEP.

    Takeaway: the two advantages act on different compartments -- light drives PP, while
    the microbial-loop advantage sustains secondary production and the nutrient inventory;
    they meet in that recycling replenishes the nutrients drawn down by the light-boosted
    production, keeping the quota limitation similar between regimes.

    Args:
        ref, notep: the two Model simulations (REF / NO-TEP) at the calibrated rate.
        year: analysis year. days: smoothing window (rolling mean, days). figsize.

    Returns:
        (fig, axs) with axs a (3, 2) array.
    """
    cref, cnot = '#2e8b57', '#d1495b'

    def tr(s, col):
        d = s.df[s.df.index.year == year]
        if col not in d.columns:
            return pd.Series(np.nan, index=d.index)
        return pd.to_numeric(d[col], errors='coerce').resample(f'{days}D',
                                                               origin='epoch').mean()

    def line(ax, cols, sumcols=False):
        for s, c, lab in [(ref, cref, 'REF'), (notep, cnot, 'NO-TEP')]:
            y = sum(tr(s, cc) for cc in cols) if sumcols else tr(s, cols)
            ax.plot(y.index, y.values, color=c, lw=1.6,
                    ls='-' if lab == 'REF' else '--', label=lab)
        ax.legend(fontsize=7.5)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b')); ax.tick_params(labelsize=7.5)

    fig, axs = plt.subplots(3, 2, figsize=figsize, sharex=True)
    line(axs[0, 0], 'Phy_source_PP.C')
    axs[0, 0].set_title('(a) Primary production', fontsize=10)
    axs[0, 0].set_ylabel(r'PP [mmol C m$^{-3}$ d$^{-1}$]', fontsize=8.5)
    line(axs[0, 1], 'Phy_limI')
    axs[0, 1].set_title('(b) Light limitation limI -- advantage 1', fontsize=10)
    axs[0, 1].set_ylabel('limI [0-1]', fontsize=8.5)

    for s, c, lab in [(ref, cref, 'REF'), (notep, cnot, 'NO-TEP')]:
        axs[1, 0].plot(tr(s, 'Phy_lim_Si').index, tr(s, 'Phy_lim_Si').values, color=c,
                       lw=1.6, ls='-' if lab == 'REF' else '--', label=f'{lab} Si')
        axs[1, 0].plot(tr(s, 'Phy_lim_P').index, tr(s, 'Phy_lim_P').values, color=c,
                       lw=1.1, ls=':', label=f'{lab} P')
    axs[1, 0].set_title('(c) Nutrient limitation (quota) -- buffered', fontsize=10)
    axs[1, 0].set_ylabel('lim_Si (solid), lim_P (dotted)', fontsize=8.5)
    axs[1, 0].legend(fontsize=7, ncol=2)
    axs[1, 0].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

    ax2 = None
    for s, c, lab in [(ref, cref, 'REF'), (notep, cnot, 'NO-TEP')]:
        axs[1, 1].plot(tr(s, 'DSi_concentration').index, tr(s, 'DSi_concentration').values,
                       color=c, lw=1.6, ls='-' if lab == 'REF' else '--', label=f'{lab} DSi')
        ax2 = axs[1, 1].twinx() if ax2 is None else ax2
        ax2.plot(tr(s, 'DIP_concentration').index, tr(s, 'DIP_concentration').values,
                 color=c, lw=1.1, ls=':', label=f'{lab} DIP')
    axs[1, 1].set_title('(d) Dissolved nutrients -- collapse without TEP', fontsize=10)
    axs[1, 1].set_ylabel(r'DSi [mmol m$^{-3}$] (solid)', fontsize=8.5)
    ax2.set_ylabel(r'DIP [mmol m$^{-3}$] (dotted)', fontsize=8.5)
    axs[1, 1].legend(fontsize=7, loc='upper right')
    axs[1, 1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

    line(axs[2, 0], ['BacF_C', 'BacA_C'], sumcols=True)
    axs[2, 0].set_title('(e) Bacterial biomass -- advantage 2 (farming)', fontsize=10)
    axs[2, 0].set_ylabel(r'Bac C [mmol m$^{-3}$]', fontsize=8.5)
    line(axs[2, 1], ['DSi_source_remineralization', 'DIP_source_remineralization'],
         sumcols=True)
    axs[2, 1].set_title('(f) Nutrient regeneration (DSi+DIP remin.)', fontsize=10)
    axs[2, 1].set_ylabel(r'regen [mmol m$^{-3}$ d$^{-1}$]', fontsize=8.5)

    fig.suptitle("TEP double advantage (REF vs NO-TEP, calibrated): light (b) drives PP; "
                 "microbial loop + recycling + dissolved nutrients (d,e,f) strongly weakened"
                 "\nwithout TEP, but nutrient limitation (c) is quota-buffered -> the 2nd "
                 "advantage does not reach PP in this model", fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return fig, axs
