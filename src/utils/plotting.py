"""
Plotting utilities for BGC model results.
Handles both single and multiple simulation visualizations.
"""
import os

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

FIGURE_PATH = path_cfg.FIGURE_PATH

# Pre-defined variable groups for common plotting scenarios
phy_nuts = ['Phy_C', 'Phy_Chl', 'NO3_concentration', 'NH4_concentration', 'DIN_concentration',
            'DIP_concentration', 'DSi_concentration']
phy_nuts_TEP_flocs = ['Phy_C', 'Phy_Chl', 'TEPC_C', "Microflocs_numconc",
                      'NO3_concentration', 'NH4_concentration', 'DIN_concentration', 'Macroflocs_diam', "Macroflocs_numconc",
                      'DIP_concentration', 'DSi_concentration', "Macroflocs_settling_vel", "Micro_in_Macro_numconc",
                      'SPMC']

phy_TEP_lim_sink = ['Phy_C', 'Phy_Chl', 'TEPC_C', 'Phy_mmDSi',
                    'Phy_mmNH4', 'Phy_mmNO3', 'Phy_mmDIP', "Phy_limNUT",
                    "Phy_lim_I", 'Phy_sink_lysis.C', 'Phy_sink_mortality.C', 'Phy_sink_exudation.C',
                    'Phy_sink_respiration.C','Phy_sink_aggregation.C', "TEP_to_PhyC_ratio", "Phy_source_PP.C"]


phy_PPsource_decomp = ['Phy_limNUT', 'Phy_limT', 'Phy_limI', ]
phy_C_SMS = ['Phy_source_PP.C',
             'Phy_sink_respiration.C', 'Phy_sink_exudation.C', 'Phy_sink_aggregation.C',
             'Phy_sink_ingestion.C', 'Phy_sink_lysis.C', 'Phy_sink_mortality.C' ]
phyDOM = ['Phy_C', 'DOCS_C', 'TEPC_C', 'DOCL_C']
nutvars = ['NO3_concentration', 'NH4_concentration',
           'DIP_concentration', 'DSi_concentration']
domvars = ['DOCS_C', 'DOCS_N', 'DOCS_P', 'DOCL_C',
           'DOCL_N', 'DOCL_P', 'TEPC_C', 'TEPC_N', 'TEPC_P']
detvars = ['DetS_C', 'DetS_N', 'DetS_P',
           'DetL_C', 'DetL_N', 'DetL_P']
hetvars = ['BacF_C', 'BacF_N', 'BacF_P',
           'BacA_C', 'BacA_N', 'BacA_P',
           'HF_C', 'HF_N', 'HF_P',
           'Cil_C', 'Cil_N', 'Cil_P']
QNvars = ['QN', 'BacA_QN', 'BacF_QN', 'HF_QN', 'Cil_QN']
stoichioPhy = ['Phy_QN', 'Phy_QP', 'Phy_QSi', 'CN', 'CP', 'CSi',
               'thetaC', 'CChl', 'thetaN']
phyvars = ['Phy_C', 'Phy_Chl', 'Phy_N', 'Phy_P', 'Phy_Si']
flocsvar = ["Microflocs_numconc", "Macroflocs_numconc",
            "Micro_in_Macro_numconc", 'Floc_diam']
flocsvar2 = ["Microflocs_numconc", "Macroflocs_numconc", "Micro_in_Macro_numconc",
             "Macroflocs_diam", 'Floc_diam',
             "Macroflocs_settling_vel", "Microflocs_TOT",
             "Microflocs_massconcentration", "Micro_in_Macro_massconcentration",
             'SPMC']
flocsvar3 = ["Microflocs_numconc", "Macroflocs_numconc", "Micro_in_Macro_numconc",
             # 'Floc_diam',
             "Macroflocs_diam", "Micro_in_Macro_concentration",
             "Macroflocs_concentration", "Microflocs_concentration",
             "Macroflocs_settling_vel",
             'Macroflocs_ppsource', 'Macroflocs_ffloss',
             'Macroflocs_breaksource', 'Macroflocs_settling_loss',
             "Macroflocs_massconcentration", "Microflocs_massconcentration",
             'SPMC']

flocs_tep_comprehensive = [
    "TEPC_C", "Macroflocs_diam", "Macroflocs_settling_vel", "Macroflocs_fyflocstrength",
    "Macroflocs_alpha_FF", "Microflocs_numconc", "Microflocs_massconcentration", "SPMC",
    "Macroflocs_alpha_PF", "Macroflocs_numconc", "Macroflocs_massconcentration", "SPMC", #, "Microflocs_TOT",
    "Macroflocs_alpha_PP", "Micro_in_Macro_numconc", "Micro_in_Macro_massconcentration", "Macroflocs_Ncnum"

    # "Macroflocs_volconcentration",


]

flocsrelatedvars = ["Phy_C", "TEPC_C",
                    "Microflocs_numconc", "Macroflocs_numconc",
                    "Micro_in_Macro_numconc", 'Floc_diam']
phyTEPflocs = ['Phy_C', 'TEPC_C', 'Macroflocs_numconc', 'Floc_diam']

phy_C_budget = {
    'sources': ['Phy_source_PP.C'],
    'sinks': ['Phy_sink_respiration.C', 'Phy_sink_exudation.C', 'Phy_sink_aggregation.C',
              'Phy_sink_ingestion.C', 'Phy_sink_lysis.C', 'Phy_sink_mortality.C']
}

bac_C_budget = {
    'sources': ['BacF_source_ing_C_assimilated'],
    'sinks': ['BacF_sink_respiration.C', 'BacF_sink_lysis.C',
              'BacF_sink_mortality.C', 'BacF_sink_ingestion.C']
}

dom_C_budget = {
    'sources': ['DOCS_source_exudation.C', 'DOCS_source_breakdown.C',
                'DOCS_source_aggregation.C', 'DOCS_source_sloppy_feeding.C'],
    'sinks': ['DOCS_sink_ingestion.C', 'DOCS_sink_remineralization.C',
              'DOCS_sink_breakdown.C', 'DOCS_sink_aggregation.C']
}



# General plotting setup
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams['axes.prop_cycle'] = (
    "cycler('linestyle', ['-', '--', ':', '-.', (0, (3, 1, 1, 1, 1, 1))]*2)+"
    "cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',"
    " '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])"
)
plt.rcParams.update({'font.size': 8})

# Define observation styling
OBS_STYLES = {
    'calibrated': {
        'marker': 'o',
        'color': 'grey',
        'markersize': 6,
        'alpha': 0.7
    },
    'non_calibrated': {
        'marker': 's',
        'color': 'gray',
        'markersize': 4,
        'alpha': 0.4
    }
}


def create_comparison_plots(
        mod_obs_data: pd.DataFrame,
        variables: List[str],
        observation_data: pd.DataFrame,
        calibrated_vars: Optional[List[str]] = None,
        save: bool = False,
        figsize: Tuple[int, int] = (10, 6)
) -> None:
    """
    Create comparison plots between model results and observations.

    Args:
        mod_obs_data: DataFrame with model results and corresponding observations
        variables: List of variables to plot
        observation_data: DataFrame with observation data
        calibrated_vars: Calibrated variables (have different style than non-calibrated vars)
        save: Whether to save the plots
        figsize: Size of each plot
    """
    for var in variables:
        plt.figure(figsize=figsize)

        # Plot model results
        plt.plot(
            mod_obs_data.index,
            mod_obs_data[f'{var}_MOD'],
            label=f'{var}_model'
        )

        is_calibrated = calibrated_vars is not None and var in calibrated_vars

        # Plot all observations
        plt.scatter(
            observation_data.index,
            observation_data[var],
            color='red',
            label='Observations',
            s=6
        )

        # Plot used observations
        plt.scatter(
            mod_obs_data.index,
            mod_obs_data[f'{var}_OBS'],
            color='orange',
            label='Observations',
            s=3
        )

        plt.xlabel('Date')
        plt.legend()

        if save:
            plt.savefig(f'Figs/{var}_comparison.png')
        plt.close()


def prepare_model_obs_data(
        models: Union[Any, List[Any], pd.DataFrame, List[pd.DataFrame]],
        observations: Optional[Any] = None,
        daily_mean: bool = True
) -> Tuple[List[pd.DataFrame], Optional[pd.DataFrame], Optional[pd.DataFrame], List[str]]:
    """
    Prepare model and observation data for plotting.

    Args:
        models: Single model/DataFrame or list of models/DataFrames
        observations: Observation data object
        daily_mean: Whether to use daily means

    Returns:
        Tuple of (model_data_list, merged_data, full_obs_data, model_names)
    """
    # Flatten nested lists of models
    models = fns.flatten_simulation_list(models)

    # Extract DataFrames and names
    model_data_list = []
    model_names = []
    name_counts = {}

    for i, model in enumerate(models):
        if isinstance(model, pd.DataFrame):
            model_data = model.copy()
            name = f'Model {i + 1}'
        else:
            model_data = model.df.copy()
            name = getattr(model, 'name', f'Model {i + 1}')

            # Handle duplicate names by adding a counter
            if name in name_counts:
                name_counts[name] += 1
                name = f"{name} ({name_counts[name]})"
            else:
                name_counts[name] = 1

        if daily_mean and model_data.index.name != 'julian_day':
            model_data['julian_day'] = model_data.index.dayofyear
            model_data = model_data.groupby('julian_day').mean()

        model_data_list.append(model_data)
        model_names.append(f"{name} {len(model_names) + 1}" if name == 'Model' else name)

    if observations is None:
        return model_data_list, None, None, model_names

    # Prepare observation data
    obs_data = observations.df.copy()

    if daily_mean and obs_data.index.name != 'julian_day':
        obs_data['julian_day'] = obs_data.index.dayofyear
        obs_data = obs_data.groupby('julian_day').mean()

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
        full_obs_data: Optional[pd.DataFrame] = None,
        model_styles: Optional[List[Dict]] = None,
        obs_kwargs: Optional[Dict] = None,
        show_full_obs: bool = False,
        add_labels: bool = True
) -> None:
    """Plot a single variable with model results and optional observations."""
    # Default styles
    if model_styles is None:
        # model_styles = [{'color': f'C{i}'} for i in range(len(model_data_list))]
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        linestyles = ['-', '--', ':', '-.', (0, (3, 1, 1, 1, 1, 1))]
        model_styles = []
        for i in range(len(model_data_list)):
            model_styles.append({
                'color': colors[i % len(colors)],
                'linestyle': linestyles[i % len(linestyles)]
            })

    default_obs_kwargs = {
        'marker': 'o',
        'linestyle': 'None',
        's': 4,
        'alpha': 0.6
    }
    obs_kwargs = {**default_obs_kwargs, **(obs_kwargs or {})}

    # Plot each model
    for model_data, name, style in zip(model_data_list, model_names, model_styles):
        if var_name in model_data.columns:
            ax.plot(model_data.index, model_data[var_name],
                    label=name if add_labels else "_" + name, **style)

    # if show_full_obs and full_obs_data is not None and var_name in full_obs_data.columns:
    #     # Check if std data is available for error bars
    #     std_col = f"{var_name}_std"
    #     if std_col in full_obs_data.columns:
    #         # Use errorbar for observations with std
    #         ax.errorbar(
    #             full_obs_data.index,
    #             full_obs_data[var_name],
    #             yerr=full_obs_data[std_col],
    #             color='orange',
    #             label='All observations' if add_labels else "_All observations",
    #             fmt='o',  # marker style
    #             capsize=3,
    #             **{k: v for k, v in obs_kwargs.items() if k not in ['marker', 'linestyle']}
    #         )
    #     else:
    #         ax.scatter(
    #         full_obs_data.index,
    #         full_obs_data[var_name],
    #         color='orange',
    #         label='All observations' if add_labels else "_All observations",
    #         **obs_kwargs)

    if merged_data is not None and f'{var_name}_OBS' in merged_data.columns:
        # Check if std data is available for error bars
        std_col = f"{var_name}_std"
        if std_col in merged_data.columns:
            # Use errorbar for observations with std
            ax.errorbar(
                merged_data.index,
                merged_data[f'{var_name}_OBS'],
                yerr=merged_data[std_col],
                color='grey',
                label='Observations' if add_labels else "_Used observations",
                fmt='o',  # marker style
                capsize=3,
                **{k: v for k, v in obs_kwargs.items() if k not in ['marker', 'linestyle', 's']}
            )
        else:
            ax.scatter(
            merged_data.index,
            merged_data[f'{var_name}_OBS'],
            color='red',
            label='Observations' if add_labels else "_Used observations",
            **{**obs_kwargs, 's': 6}
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
        variables: List[str],
        observations: Optional[Any] = observations.Obs(station='MOW1_202503'),
        calibrated_vars: Optional[List[str]] = None,
        daily_mean: bool = True,
        show_full_obs: bool = False,
        ncols: int = 2,
        figsize: Optional[Tuple[int, int]] = None,
        model_styles: Optional[List[Dict]] = None,
        subplot_labels: bool = True,
        label_position: str = 'top_left',
        save: bool = False,
        filename: str = None,
        fnametimestamp: bool = True,
        **plot_kwargs
) -> Tuple[plt.Figure, np.ndarray]:
    """
    Create comprehensive plots of model results with optional observations.

    Args:
        models: Single model or list of models
        variables: List of variables to plot
        observations: Optional observation data
        calibrated_vars: Calibrated variables (have different style than non-calibrated vars)
        daily_mean: Whether to use daily means
        show_full_obs: Whether to show full observation range
        ncols: Number of columns in subplot grid
        figsize: Figure size (default: auto-calculated)
        model_styles: List of style dictionaries for each model
        subplot_labels: Whether to add subplot labels (a, b, c, etc.)
        label_position: Position for subplot labels
        save: Whether to save the figure
        filename: Filename for saved figure
        **plot_kwargs: Additional plotting parameters

    Returns:
        Figure and axes array
    """
    # Prepare data
    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #     print(models[0].df['Phy_source_PP.C'])
    model_data_list, merged_data, full_obs_data, model_names = prepare_model_obs_data(
        models, observations, daily_mean
    )

    # Calculate grid dimensions
    nrows = (len(variables) + ncols - 1) // ncols
    if figsize is None:
        figsize = (5 * ncols, 4 * nrows)

    # Create figure and axes
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True)
    if nrows == 1 and ncols == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # Plot each variable
    for i, var in enumerate(variables):
        if i < len(axes):
            plot_variable(
                axes[i],
                model_data_list,
                model_names,
                var,
                merged_data,
                full_obs_data,
                model_styles,
                show_full_obs=show_full_obs,
                add_labels=(i == 0),  # Only add labels on the first subplot
                **plot_kwargs
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
        handles, labels = [], []
        # Collect unique handles and labels from all subplots
        for ax in axes:
            h, l = ax.get_legend_handles_labels()
            for handle, label in zip(h, l):
                if label not in labels:  # Only add unique labels
                    handles.append(handle)
                    labels.append(label)

        # Create a single legend for the figure
        fig.legend(
            handles,
            labels,
            loc='center right',
            bbox_to_anchor=(0.98, 0.5),
            fontsize=plot_kwargs.get('legend_fontsize', 8)
        )
    # if plot_kwargs.get('add_legend', True):
    #     handles, labels = axes[0].get_legend_handles_labels()
    #     fig.legend(
    #         handles,
    #         labels,
    #         loc='center right',
    #         bbox_to_anchor=(0.98, 0.5),
    #         fontsize=plot_kwargs.get('legend_fontsize', 8)
    #     )

    # plt.tight_layout()

    # print(save, filename, os.getcwd())
    if save:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S_')
        if filename:
            fname = f"{filename}_{timestamp}.png" if fnametimestamp else filename
        else:
            fname = f"{timestamp}_model_results.png"
        plt.savefig(os.path.join(FIGURE_PATH, fname), dpi=300, bbox_inches='tight')

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
        default_subplot_size: Tuple[float, float] = (3.6, 5.7),
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
            # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            #     print(cumulative_sources[source])

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

            continue
    plt.tight_layout()

    # Save figure if requested
    if save:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S_')
        fname = f"{filename}_{timestamp}.png" if fnametimestamp else filename
        save_figure(fig, figname=fname)

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
        fig_base_path = os.path.join(FIGURE_PATH, name)
        plt.savefig(f'{fig_base_path}_opt_evol{rawcost * "_raw"}.png')

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
        fig_base_path = os.path.join(FIGURE_PATH, f'{datestr}{name}')
        plt.savefig(f'{fig_base_path}_pars_optim{rawcost * "_raw"}.png')

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
                savetime: bool = False,
                figsdir: str = os.path.join(FIGURE_PATH, ''),
                figname: str = '',
                figsdpi: int = 200) -> None:
    """Save figure with optional timestamp."""
    timestamp = datetime.now().strftime('%Y%m%d-%H%M%S_' if savetime else '%Y%m%d_')
    fig.savefig(f"{figsdir}{timestamp}{figname}", dpi=figsdpi)
