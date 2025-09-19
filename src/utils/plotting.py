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
DEFAULT_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                  '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
DEFAULT_LINESTYLES = ['-', '--', ':', '-.', (0, (3, 1, 1, 1, 1, 1))]
DEFAULT_LEGEND_FONTSIZE = 8
DEFAULT_FIGURE_DPI = 300

# General plotting setup
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams['axes.prop_cycle'] = (
    f"cycler('linestyle', {DEFAULT_LINESTYLES * 2})+"
    f"cycler('color', {DEFAULT_COLORS})"
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

def prepare_model_obs_data(
        models: Union[Any, List[Any], pd.DataFrame, List[pd.DataFrame]],
        observations: Optional[Any] = None,
        daily_mean: bool = True,
        variables_to_plot: Optional[List[str]] = None
) -> Tuple[List[pd.DataFrame], Optional[pd.DataFrame], Optional[pd.DataFrame], List[str]]:
    """
    Prepare model and observation data for plotting.

    Performance optimization: Only copies necessary columns from potentially huge DataFrames.

    Args:
        models: Single model/DataFrame or list of models/DataFrames
        observations: Observation data object
        daily_mean: Whether to use daily means
        variables_to_plot: List of variables that will be plotted (for column standardization)

    Returns:
        Tuple of (model_data_list, merged_data, full_obs_data, model_names)
    """
    # Flatten nested lists of models
    models = fns.flatten_simulation_list(models)

    # Determine which columns we actually need to avoid copying massive DataFrames
    required_columns = set()
    if variables_to_plot:
        required_columns.update(variables_to_plot)

    # Always include time-related columns for daily_mean functionality
    if daily_mean:
        required_columns.add('julian_day')

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

        if daily_mean and model_data.index.name != 'julian_day':
            model_data['julian_day'] = model_data.index.dayofyear
            model_data = model_data.groupby('julian_day').mean()

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

    # Prepare observation data
    obs_data = observations.df.copy()

    if daily_mean and obs_data.index.name != 'julian_day':
        obs_data['julian_day'] = obs_data.index.dayofyear
        obs_data = obs_data.groupby('julian_day').mean()
    elif not daily_mean and obs_data.index.name == 'julian_day':
        # Convert julian day index to datetime to match model data
        # Assumes year alignment with model data (first model's year is used)
        model_year = model_data_list[0].index[0].year
        obs_data.index = pd.to_datetime(obs_data.index.astype(int), format='%j').map(
            lambda x: x.replace(year=model_year)
        )

    # Merge with first model for the model period
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
        add_labels: bool = True,
        calibrated_vars: Optional[List[str]] = None
) -> None:
    """Plot a single variable with model results and optional observations."""
    # Default styles
    if model_styles is None:
        model_styles = get_model_styles(len(model_data_list))

    # Default observation kwargs (will be overridden by calibration-specific styles)
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

        # Merge with any user-provided obs_kwargs, giving priority to calibration-specific styles
        final_obs_style = {**base_obs_kwargs, **selected_obs_style}

        # Check if std data is available for error bars
        std_col = f"{var_name}_std"
        if std_col in merged_data.columns:
            # Use errorbar for observations with std
            ax.errorbar(
                merged_data.index,
                merged_data[f'{var_name}_OBS'],
                yerr=merged_data[std_col],
                label='Observations' if add_labels else "_Used observations",
                fmt=final_obs_style['marker'],
                elinewidth=1.5,
                capsize=3.,
                capthick=1.5,
                solid_capstyle='round',
                **{k: v for k, v in final_obs_style.items() if k not in ['marker', 'linestyle', 's']}
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
    clean_name = var_info.get('cleanname', var_name)
    units = var_info.get('munits' if var_name.startswith('m') else 'units', '')
    long_name = var_info.get('longname', var_name)

    ax.set_ylabel(f'{fns.cleantext(clean_name)} [{fns.cleantext(units)}]')
    ax.set_title(fns.cleantext(long_name))


def plot_results(
        models: Union[Any, List[Any]],
        variables: Union[List[str], PlottedVariablesSet],
        observations: Optional[Any] = observations.Obs(station='MOW1_biweekly_202509_noPhaeo'),
        calibrated_vars: Optional[List[str]] = None,
        daily_mean: bool = True,
        ncols: Optional[int] = None,
        figsize: Optional[Tuple[int, int]] = None,
        model_styles: Optional[List[Dict]] = None,
        subplot_labels: bool = True,
        label_position: str = 'top_left',
        save: bool = False,
        filename: Optional[str] = None,
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
        ncols: Number of columns in subplot grid (auto-detected from PlottedVariablesSet if None)
        figsize: Figure size (auto-detected from PlottedVariablesSet if None, otherwise auto-calculated)
        model_styles: List of style dictionaries for each model
        subplot_labels: Whether to add subplot labels (a, b, c, etc.)
        label_position: Position for subplot labels
        save: Whether to save the figure
        filename: Filename for saved figure (auto-generated from PlottedVariablesSet.name if None)
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
        models, observations, daily_mean, variables_list
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

    # Extract legend-related kwargs before passing to plot_variable
    plot_var_kwargs = {k: v for k, v in plot_kwargs.items()
                       if k not in ['add_legend', 'legend_fontsize', 'legend_position']}
    
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
                **plot_var_kwargs
            )

            # Format x-axis
            axes[i].xaxis.set_major_formatter(mdates.DateFormatter('%d/%m'))
            for label in axes[i].get_xticklabels():
                label.set_rotation(45)

            # Add subplot label if requested
            if subplot_labels:
                add_subplot_label(axes[i], i, position=label_position)
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

    # plt.tight_layout()

    # Save figure using centralized function
    if save:
        # Extract PlottedVariablesSet for intelligent naming
        var_set = variables if isinstance(variables, PlottedVariablesSet) else None
        save_figure(fig, filename=filename, variable_set=var_set, add_timestamp=fnametimestamp)

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
        models, observations, daily_mean
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
        if daily_mean:
            # For julian day indices (daily_mean=True)
            min_day = df.index.min()
            max_day = df.index.max()
            new_index = np.linspace(min_day, max_day,
                                    num=len(df.index) * increase_resolution_factor)
            df_high_res = pd.DataFrame(index=new_index).join(
                df[valid_sources + valid_sinks].astype(float),
                how='outer').interpolate(method='linear')
        else:
            # For datetime indices (daily_mean=False)
            new_index = pd.date_range(
                start=df.index.min(),
                end=df.index.max(),
                periods=len(df.index) * increase_resolution_factor
            )
            df_high_res = pd.DataFrame(index=new_index).join(
                df[valid_sources + valid_sinks].astype(float),
                how='outer').interpolate(method='time')

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
                    label=source,
                    alpha=0.7
                )
            else:
                ax.fill_between(
                    df_high_res.index,
                    cumulative_sources[valid_sources[j - 1]],
                    cumulative_sources[source],
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
                    label=sink,
                    alpha=0.7
                )
            else:
                ax.fill_between(
                    df_high_res.index,
                    -cumulative_sinks[valid_sinks[j - 1]],
                    -cumulative_sinks[sink],
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
    ax.scatter(df[generationname], df[cost_col], alpha=alpha)

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


def plot_parameter_vs_cost(df: pd.DataFrame,
                           param: str,
                           costname: str = 'cost',
                           alpha: float = 0.2,
                           rawcost: bool = False,
                           ax: Optional[plt.Axes] = None) -> plt.Axes:
    """
    Plot parameter values against optimization cost.

    Args:
        df: Optimization results DataFrame
        param: Parameter name to plot
        costname: Name of cost column
        alpha: Scatter plot transparency
        rawcost: Whether to plot raw cost
        ax: Optional axes to plot on
    """
    if ax is None:
        fig, ax = plt.subplots()
    elif not isinstance(ax, plt.Axes):
        raise TypeError("ax must be a matplotlib Axes object")

    cost_col = costname + rawcost * '_raw'

    # Main scatter plot
    ax.scatter(df[param], df[cost_col], c='#1f77b4', alpha=alpha)

    # Get parameter info from reference values
    param_info = varinfos.ref_values[param]
    complete_name = param_info['complete_name']
    symbol = param_info.get('symbol', '')
    reference_value = param_info['reference_value']

    # Set title
    ax.set_title(fns.cleantext(symbol) if symbol else complete_name)

    # Labels
    ax.set_xlabel(param)
    ax.set_ylabel(cost_col)

    # Highlight best point
    max_cost_idx = df[costname].idxmax()
    best_value = df[param][max_cost_idx]
    ax.plot(best_value, df[cost_col][max_cost_idx], 'ro')

    # Add value labels
    rounded_value = f"{best_value:.3g}"
    ax.text(best_value, df[cost_col][max_cost_idx],
            rounded_value, color='red')

    # Add reference line and value
    ax.axvline(x=reference_value, color='blue', linestyle='--')
    ax.text(reference_value, ax.get_ylim()[0], f'{reference_value}',
            color='blue', ha='center', va='bottom',
            bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

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
                costname=costname,
                alpha=alpha,
                rawcost=rawcost,
                ax=axs[i]
            )
        else:
            axs[i].set_visible(False)

    plt.tight_layout()

    if savefig:
        datestr = datetime.now().strftime('%Y%m%d_') if dateinname else ''
        filename = f'{datestr}{name}_pars_optim{rawcost * "_raw"}.png'
        save_figure(fig, filename=filename)

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
