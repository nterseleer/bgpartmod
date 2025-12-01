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

# Line styles: 3 distinct styles to combine with 14 colors
# Avoiding dot-only ':' style as it becomes too discrete
# Using solid, dashed, and dash-dot for good visibility
DEFAULT_LINESTYLES = ['-', '--', '-.', ':']

# Budget plot color palettes (cool colors for sources, warm colors for sinks)
BUDGET_SOURCE_COLORS = ['#2E86AB', '#06A77D', '#81B214', '#A7C957', '#4ECDC4', '#45B7D1']
BUDGET_SINK_COLORS = ['#D62828', '#F77F00', '#FCBF49', '#EE6C4D', '#E63946', '#F4A261']

DEFAULT_LEGEND_FONTSIZE = 8
DEFAULT_FIGURE_DPI = 300

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
plt.rcParams.update({'font.size': 6})

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


def prepare_model_obs_data(
        models: Union[Any, List[Any], pd.DataFrame, List[pd.DataFrame]],
        observations: Optional[Any] = None,
        daily_mean: bool = True,
        variables_to_plot: Optional[List[str]] = None,
        time_filter = None
) -> Tuple[List[pd.DataFrame], Optional[pd.DataFrame], Optional[pd.DataFrame], List[str]]:
    """
    Prepare model and observation data for plotting.
    Assumes both models and observations have DatetimeIndex.

    Performance optimization: Only copies necessary columns from potentially huge DataFrames.

    Args:
        models: Single model/DataFrame or list of models/DataFrames (with DatetimeIndex)
        observations: Observation data object (with DatetimeIndex)
        daily_mean: Whether to use daily means (uses resample)
        variables_to_plot: List of variables that will be plotted (for column standardization)
        time_filter: Time filtering criterion (None, year string, date range tuple, or slice)

    Returns:
        Tuple of (model_data_list, merged_data, full_obs_data, model_names)
    """
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

        # Handle duplicate names by adding a counter
        original_name = name
        counter = 1
        while name in model_names:
            counter += 1
            name = f"{original_name} ({counter})"

        # Apply time filtering BEFORE daily_mean to preserve temporal resolution
        if time_filter is not None:
            model_data = _filter_by_time(model_data, time_filter)

        # Apply daily_mean using resample (assumes DatetimeIndex)
        if daily_mean:
            model_data = model_data.resample('D').mean()

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
        obs_kwargs: Optional[Dict] = None,
        obs_kwargs_calibrated: Optional[Dict] = None,
        obs_kwargs_non_calibrated: Optional[Dict] = None,
        add_labels: bool = True,
        calibrated_vars: Optional[List[str]] = None,
        show_subplot_titles: bool = True,
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

    # Plot each model
    for model_data, name, style in zip(model_data_list, model_names, model_styles):
        if var_name in model_data.columns:
            ax.plot(model_data.index, model_data[var_name],
                    label=name if add_labels else "_" + name, **style)


    if merged_data is not None and f'{var_name}_OBS' in merged_data.columns:
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

    # Format labels and title
    var_info = varinfos.doutput.get(var_name.lstrip('m'), {})
    clean_name = var_info.get('cleanname', None)

    # If cleanname not defined in varinfos, use var_name but escape underscores
    # to prevent matplotlib mathtext errors (double subscript issues)
    if clean_name is None:
        clean_name = var_name.replace('_', r'\_')

    units = var_info.get('munits' if var_name.startswith('m') else 'units', '')
    long_name = var_info.get('longname', var_name)

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
        daily_mean: bool = True,
        time_filter = None,
        ncols: Optional[int] = None,
        figsize: Optional[Tuple[int, int]] = None,
        model_styles: Optional[List[Dict]] = None,
        subplot_labels: bool = True,
        subplot_label_fontsize: int = 10,
        label_position: str = 'top_left',
        show_subplot_titles: bool = True,
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
        daily_mean: Whether to use daily means
        time_filter: Time filtering criterion (None, year string, date range tuple, or slice)
        ncols: Number of columns in subplot grid (auto-detected from PlottedVariablesSet if None)
        figsize: Figure size (auto-detected from PlottedVariablesSet if None, otherwise auto-calculated)
        model_styles: List of style dictionaries for each model
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

    # Prepare data
    model_data_list, merged_data, _, model_names = prepare_model_obs_data(
        models, observations, daily_mean, variables_list, time_filter
    )

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
        fig.canvas.manager.set_window_title(variables.name.replace('_', ' ').title())

    # Extract legend-related and obs-specific kwargs before passing to plot_variable
    legend_kwargs = ['add_legend', 'legend_fontsize', 'legend_position']
    obs_specific_kwargs = ['obs_kwargs_calibrated', 'obs_kwargs_non_calibrated']
    plot_var_kwargs = {k: v for k, v in plot_kwargs.items()
                       if k not in legend_kwargs}
    
    # Plot each variable
    for i, var in enumerate(variables_list):
        if i < len(axes):
            plot_variable(
                axes[i],
                model_data_list,
                model_names,
                var,
                merged_data,
                model_styles,
                add_labels=(i == 0),  # Only add labels on the first subplot
                calibrated_vars=calibrated_vars,
                show_subplot_titles=show_subplot_titles,
                ylabel_pad=ylabel_pad,
                ylabel_fontsize=ylabel_fontsize,
                tick_label_pad=tick_label_pad,
                tick_labelsize=tick_labelsize,
                tick_length=tick_length,
                tick_width=tick_width,
                spine_width=spine_width,
                **plot_var_kwargs
            )

            # Format x-axis
            axes[i].xaxis.set_major_formatter(mdates.DateFormatter(date_format))
            for label in axes[i].get_xticklabels():
                label.set_rotation(xlabel_rotation)

            # Add subplot label if requested
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

        # Add observations if present (check all subplots)
        if merged_data is not None and len(axes) > 0:
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
        daily_mean: bool = True,
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
        daily_mean: Whether to use daily means
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
        models, observations, daily_mean, time_filter=time_filter,
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


def plot_element_distribution_stacked(model_output, element_vars, element_name=None,
                                      time_var='time', **kwargs):
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
        else:
            element_name = 'Element'
            unit = 'mmol m⁻³'
    else:
        unit = f'mmol {element_name} m⁻³'

    # Filter available variables
    available_vars = [var for var in element_vars if var in model_output.columns]

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

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Create stacked area plot
    if compartment_totals:
        ax.stackplot(model_output.index,
                     *compartment_totals.values(),
                     labels=compartment_totals.keys(),
                     alpha=0.8)

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

    # Formatting
    ax.set_xlabel('Time [d]')
    ax.set_ylabel(f'{element_name} concentration [{unit}]')
    ax.set_title(f'{element_name} Distribution Across Model Compartments')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig, ax


def plot_kd_contributions_stacked(model_output, kd_contrib_vars=None, **kwargs):
    """
    Create a stacked area plot showing contributions to light attenuation coefficient (kd).

    Parameters:
    -----------
    model_output : DataFrame
        Model output containing the kd contribution variables
    kd_contrib_vars : list, optional
        List of kd contribution variables to plot. If None, uses default list.
    **kwargs : dict
        Additional plotting arguments (figsize, save, etc.)

    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    import matplotlib.pyplot as plt
    from src.config_model import vars_to_plot

    # Use default list if not provided
    if kd_contrib_vars is None:
        kd_contrib_vars = vars_to_plot.kd_contributions_list

    # Filter to available variables
    available_vars = [var for var in kd_contrib_vars if var in model_output.columns]

    if not available_vars:
        print("Warning: No kd contribution variables found in model output")
        return None, None

    # Extract figsize from kwargs or use default
    figsize = kwargs.get('figsize', (14, 8))

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data for stacking (each variable as a separate series)
    # Ensure data is numeric and handle NaN values
    data_to_stack = []
    for var in available_vars:
        # Convert to numeric, replacing any non-numeric with NaN
        series = pd.to_numeric(model_output[var], errors='coerce')
        # Fill NaN with 0 for stacking (or you could use interpolation)
        series = series.fillna(0)
        data_to_stack.append(series.values)

    # Create clean labels (remove 'kd_contrib_' prefix and format nicely)
    labels = []
    for var in available_vars:
        label = var.replace('kd_contrib_', '')
        # Special formatting for certain components
        if label == 'Micro_in_Macro':
            label = 'Microflocs in Macroflocs'
        labels.append(label)

    # Ensure index is numeric/datetime (not object)
    # Convert index to numeric if it's a DatetimeIndex
    if isinstance(model_output.index, pd.DatetimeIndex):
        x_values = mdates.date2num(model_output.index)
    else:
        # Assume it's already numeric (days)
        x_values = model_output.index.values

    # Create stacked area plot with distinct colors
    colors = plt.cm.tab10.colors[:len(available_vars)]
    if len(available_vars) > 10:
        colors = plt.cm.tab20.colors[:len(available_vars)]

    ax.stackplot(x_values,
                 *data_to_stack,
                 labels=labels,
                 colors=colors,
                 alpha=0.8)

    # Add total kd line if available
    if 'Phy_kd' in model_output.columns:
        ax.plot(model_output.index, model_output['Phy_kd'],
                'k--', linewidth=2, label='Total k$_d$', zorder=10)

    # Formatting
    ax.set_xlabel('Time [d]', fontsize=12)
    ax.set_ylabel('Light attenuation coefficient [m$^{-1}$]', fontsize=12)
    ax.set_title('Decomposition of Light Attenuation Coefficient (k$_d$)', fontsize=14)
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=10)
    ax.grid(True, alpha=0.3)

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
        filter_badlnl: Optional[float] = None
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
                              dateinname: bool = False) -> plt.Figure:
    """
    Create comprehensive optimization summary plots.
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
                show_reference=True
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
        name: str = 'optim_comparison'
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

    Returns:
        Matplotlib figure
    """
    from src.utils import optimization as optim

    # Load optimizations if given as strings
    opts = [
        optim.Optimization.load_existing(opt) if isinstance(opt, str) else opt
        for opt in optimizations
    ]

    # Ensure all have processed results
    for opt in opts:
        if not hasattr(opt, 'df') or opt.df is None:
            opt.process_results()

    # Get union of all optimized parameters
    if parameters is None:
        parameters = sorted(set().union(*[
            set(opt.config['optimized_parameters']) for opt in opts
        ]))

    n_params = len(parameters)
    n_opts = len(opts)

    # Dynamic ncols = n_params (one column per parameter)
    ncols = n_params
    nrows = 2

    # Auto-calculate figsize if not provided
    if figsize is None:
        figsize = (max(12, 3.5 * ncols), 16)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)

    colors = DEFAULT_COLORS[:n_opts]

    # Track which subplot has all optimizations (for legend placement)
    legend_col = None

    for col_idx, param in enumerate(parameters):
        ax_violin = axes[0, col_idx]
        ax_scatter = axes[1, col_idx]

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
        for i, opt in enumerate(opts):
            if param in opt.df.columns:
                badlnl = opt.config.get('badlnl', -100000.)
                plot_parameter_vs_cost(
                    opt.df, param,
                    ax=ax_scatter,
                    costname='cost',
                    alpha=alpha,
                    color=colors[i],
                    label=opt.name,
                    highlight_best=True,
                    show_reference=(i == 0),  # Only first opt shows reference
                    show_labels=False,  # We set labels manually below
                    filter_badlnl=badlnl
                )

        # Set labels manually (shared across all optimizations)
        param_info = varinfos.ref_values.get(param, {})
        symbol = param_info.get('symbol', param)
        ax_scatter.set_xlabel(fns.cleantext(symbol) if symbol else param, fontsize=8)
        ax_scatter.set_ylabel('Cost (log-likelihood)', fontsize=8)
        ax_scatter.tick_params(labelsize=7)
        ax_scatter.grid(alpha=0.3)

        # Hide y tick labels except for first column
        if col_idx > 0:
            ax_scatter.set_ylabel('')
            ax_scatter.tick_params(labelleft=False)

        # Track if this column has all optimizations (for legend)
        if legend_col is None and len(positions) == n_opts:
            legend_col = col_idx

        # Add legend only to first column that has ALL optimizations
        if col_idx == legend_col:
            ax_scatter.legend(fontsize=6, loc='best')

    plt.tight_layout()

    if savefig:
        save_figure(fig, filename=f'{name}_comparison.png')

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
        pars_info.append({
            "Name": param_info['complete_name'],
            "Symbol": fns.cleantext(param_info.get('symbol', '')),
            "Ref. value": f"{param_info['reference_value']:.3g}",
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
