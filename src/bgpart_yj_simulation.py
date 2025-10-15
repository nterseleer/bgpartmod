#!/usr/bin/env python3
"""
Visualization tools for bgpart_yj_simulation results comparison.

This module provides functions to compare common optimization vs independent optimization
results across multiple phytoplankton strains.
"""

import json
import sys
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Add the project root to Python path for imports
sys.path.append(os.getcwd())

from src.utils import observations, functions as fns
from src.utils.plotting import save_figure
from src.config_model import varinfos
from src.utils import phys

# Define strain names, light regimes, and their display colors
STRAINS = ['Chaetoceros', 'Lauderia', 'Odontella', 'Rhizosolenia', 'Skeletonema', 'Thalassiosira']
LIGHT_REGIMES = ['ll', 'ml', 'hl']  # Low Light, Medium Light, High Light

STRAIN_COLORS = {
    'Chaetoceros': '#1f77b4',
    'Lauderia': '#ff7f0e', 
    'Odontella': '#2ca02c',
    'Rhizosolenia': '#d62728',
    'Skeletonema': '#9467bd',
    'Thalassiosira': '#8c564b'
}

# Line styles for different light regimes
LIGHT_STYLES = {
    'll': '-',      # Solid line for Low Light
    'ml': '--',     # Dashed line for Medium Light
    'hl': '-.',     # Dash-dot line for High Light
}

LIGHT_LABELS = {
    'll': 'Low Light',
    'ml': 'Medium Light', 
    'hl': 'High Light'
}

# Observation markers for different light regimes
LIGHT_MARKERS = {
    'll': 'v',      # Triangle down for Low Light
    'ml': 'o',      # Circle for Medium Light  
    'hl': 's',      # Square for High Light
}

LIGHT_OBS_NAMES = {
    'll': 'low',
    'ml': 'medium',
    'hl': 'high'
}

CO_COLOR = '#808080'  # Gray color for common optimization


def load_simulation_data(base_dir='../bgpart_yj_simulation'):
    """Load all simulation data from bgpart_yj_simulation directory with multi-light support."""
    data = {}
    
    # Detect available light regimes by scanning directory
    available_lights = set()
    base_path = Path(base_dir)
    
    for strain in STRAINS:
        data[strain] = {}
        
        # Initialize light regime containers
        data[strain]['io_results'] = {}
        data[strain]['co_results'] = {}
        data[strain]['io_params'] = {}
        data[strain]['co_params'] = {}
        
        # Load data for each light regime
        for light in LIGHT_REGIMES:
            # Try to load independent optimization data
            io_dir = base_path / f'{strain}_io_{light}'
            if io_dir.exists():
                available_lights.add(light)
                try:
                    data[strain]['io_results'][light] = np.load(io_dir / 'results.npy', allow_pickle=True)
                    
                    # Load optimization summary
                    opt_files = list(io_dir.glob('OPT*_SUMMARY.json'))
                    if opt_files:
                        with open(opt_files[0], 'r') as f:
                            data[strain]['io_params'][light] = json.load(f)['best_parameters']
                    else:
                        data[strain]['io_params'][light] = {}
                except Exception as e:
                    print(f"Warning: Could not load {strain}_io_{light}: {e}")
            
            # Try to load common optimization data
            co_dir = base_path / f'{strain}_co_{light}'
            if co_dir.exists():
                try:
                    data[strain]['co_results'][light] = np.load(co_dir / 'results.npy', allow_pickle=True)
                except Exception as e:
                    print(f"Warning: Could not load {strain}_co_{light}: {e}")
    
    # Load common optimization parameters (if available)
    co_summary_file = base_path / 'co_SUMMARY.json'
    if co_summary_file.exists():
        try:
            with open(co_summary_file, 'r') as f:
                co_summary = json.load(f)
                for strain in STRAINS:
                    for light in available_lights:
                        data[strain]['co_params'][light] = co_summary.get('best_parameters', {})
        except Exception as e:
            print(f"Warning: Could not load common optimization parameters: {e}")
    
    # Store available light regimes in data for reference
    data['_meta'] = {'available_lights': sorted(list(available_lights))}
    
    return data


def load_observations():
    """Load observational data for each strain and light regime using the observations module."""
    obs_data = {}
    
    for strain in STRAINS:
        obs_data[strain] = {}
        
        # Load observations for each light regime
        for light_regime in LIGHT_REGIMES:
            light_name = LIGHT_OBS_NAMES[light_regime]  # ll->low, ml->medium, hl->high
            try:
                obs = observations.Obs(name=f"auria_data_{strain}_{light_name}", 
                                     station=f"{strain}_{light_name}")
                obs.format_auria_obs()
                obs_data[strain][light_regime] = obs.df
            except Exception as e:
                print(f"Warning: Could not load {strain}_{light_name} observations: {e}")
                obs_data[strain][light_regime] = None
    
    return obs_data


def _calculate_diagnostic_variables(df, variable, setup=None):
    """
    Calculate diagnostic variables using operations from varinfos.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing simulation variables
    variable : str
        Variable name to calculate (if it's a diagnostic variable)
    setup : object, optional
        Physical setup object for environmental conditions
        
    Returns:
    --------
    pd.Series or None
        Calculated variable values, or None if not a diagnostic variable
    """
    # Check if variable is in varinfos and has an 'oprt' (operation) definition
    var_info = varinfos.doutput.get(variable.lstrip('m'), {})
    oprt = var_info.get('oprt')
    
    if oprt is None:
        return None
    
    try:
        # Use functions.eval_expr to calculate the diagnostic variable
        # Default setup if not provided
        if setup is None:
            setup = phys.Setup()
            
        calculated_values = fns.eval_expr(
            oprt,
            subdf=df,
            fulldf=df,
            setup=setup,
            varinfos=varinfos
        )
        
        # Convert to pandas Series if it isn't already
        if not isinstance(calculated_values, pd.Series):
            calculated_values = pd.Series(calculated_values, index=df.index)
            
        return calculated_values
        
    except Exception as e:
        print(f"Warning: Could not calculate diagnostic variable '{variable}' with operation '{oprt}': {e}")
        return None


def _get_variable_label(var_name):
    """
    Get formatted variable label using varinfos similar to plotting.py logic.
    
    Parameters:
    -----------
    var_name : str
        Raw variable name from simulation results
    
    Returns:
    --------
    str
        Formatted label for publication
    """
    var_info = varinfos.doutput.get(var_name.lstrip('m'), {})
    clean_name = var_info.get('cleanname', var_name)
    units = var_info.get('munits' if var_name.startswith('m') else 'units', '')
    long_name = var_info.get('longname', var_name)
    
    # Create formatted label
    if clean_name != var_name and units:
        return f"{clean_name}\n[{units}]"
    elif long_name != var_name and units:
        return f"{long_name}\n[{units}]"
    elif units:
        return f"{var_name}\n[{units}]"
    else:
        return var_name


def _get_parameter_label(param_name):
    """
    Get formatted parameter label using varinfos for parameter names.
    
    Parameters:
    -----------
    param_name : str
        Raw parameter name
    
    Returns:
    --------
    str
        Formatted label for publication
    """
    param_info = varinfos.ref_values.get(param_name, {})
    symbol = param_info.get('symbol', param_name).replace("\\", "")
    units = param_info.get('units', '')
    complete_name = param_info.get('complete_name', param_name)
    
    # Create formatted label - use symbol if available, otherwise use name
    if symbol != param_name and units:
        return f"${symbol}$\n[{units}]"
    elif complete_name != param_name and units:
        # Truncate long names and add line breaks
        name_parts = complete_name.split()
        if len(name_parts) > 2:
            formatted_name = " ".join(name_parts[:2]) + "\n" + " ".join(name_parts[2:])
        else:
            formatted_name = complete_name
        return f"{formatted_name}\n[{units}]"
    elif units:
        return f"{param_name.replace('+', ' ')}\n[{units}]"
    else:
        return param_name.replace('+', ' ')


def _prepare_time_series_data(results_dict, variable, setup=None):
    """
    Convert simulation results to daily mean time series following plotting.py logic.
    Supports both direct variables and diagnostic variables that need calculation.
    
    Parameters:
    -----------
    results_dict : dict
        Results dictionary with 'timestamps' and 'variables' keys
    variable : str
        Variable name to extract or calculate
    setup : object, optional
        Physical setup object for diagnostic variable calculations
    
    Returns:
    --------
    pd.DataFrame or None
        DataFrame with julian_day index and variable column, or None if variable not found
    """
    if 'timestamps' not in results_dict or 'variables' not in results_dict:
        return None
    
    timestamps = pd.to_datetime(results_dict['timestamps'])
    
    # Check if variable exists directly in results
    if variable in results_dict['variables']:
        df = pd.DataFrame({variable: results_dict['variables'][variable]}, index=timestamps)
    else:
        # Try to calculate as diagnostic variable
        # First create a full DataFrame with all available variables
        full_df = pd.DataFrame(results_dict['variables'], index=timestamps)
        
        # Calculate diagnostic variable
        diagnostic_values = _calculate_diagnostic_variables(full_df, variable, setup)
        
        if diagnostic_values is None:
            return None
            
        df = pd.DataFrame({variable: diagnostic_values}, index=timestamps)
    
    # Convert to daily means
    df['julian_day'] = df.index.dayofyear
    return df.groupby('julian_day').mean()


def _plot_simulation_lines(ax, data, strain, variable, i, j):
    """Plot simulation lines for a given strain and variable across all light regimes."""
    available_lights = data.get('_meta', {}).get('available_lights', [])
    
    for light in available_lights:
        # Plot common optimization with gray color and solid line
        if light in data[strain]['co_results']:
            co_results = data[strain]['co_results'][light].item()
            co_daily = _prepare_time_series_data(co_results, variable)
            if co_daily is not None:
                ax.plot(co_daily.index, co_daily[variable], 
                       linestyle='-', color=CO_COLOR, 
                       alpha=0.8, linewidth=2,
                       label=f'Common opt. ({LIGHT_LABELS[light]})' if i == 0 and j == 0 else "")
        
        # Plot independent optimization with strain color and light-specific line style
        if light in data[strain]['io_results']:
            io_results = data[strain]['io_results'][light].item()
            io_daily = _prepare_time_series_data(io_results, variable)
            if io_daily is not None:
                ax.plot(io_daily.index, io_daily[variable], 
                       linestyle=LIGHT_STYLES[light], color=STRAIN_COLORS[strain], 
                       alpha=0.8, linewidth=2,
                       label=f'{strain} {LIGHT_LABELS[light]} (IO)' if i == 0 and j == 0 else "")


def _plot_observations(ax, obs_data, strain, variable, j):
    """Plot observation points for a given strain and variable across all light regimes."""
    if not obs_data or strain not in obs_data:
        return
    
    # Plot observations for each light regime with specific markers
    for light_regime in LIGHT_REGIMES:
        if light_regime not in obs_data[strain] or obs_data[strain][light_regime] is None:
            continue
            
        obs_df = obs_data[strain][light_regime]
        
        # Check if variable exists directly in observations
        if variable in obs_df.columns:
            values = obs_df[variable]
        else:
            # Try to calculate diagnostic variable from observations
            diagnostic_values = _calculate_diagnostic_variables(obs_df, variable)
            if diagnostic_values is None:
                continue
            values = diagnostic_values
        
        # Plot observations that are not NaN with light-specific markers
        mask = ~pd.isna(values)
        if mask.any():  # Only plot if there are non-NaN values
            ax.scatter(obs_df.index[mask], values.loc[mask], 
                      color=STRAIN_COLORS[strain], s=40, alpha=0.8,
                      marker=LIGHT_MARKERS[light_regime], edgecolor='white', linewidth=0.8,
                      label=f'{strain} {LIGHT_LABELS[light_regime]} obs.' if j == 0 else "")


def _format_subplot(ax, variable, strain, i, j, n_strains):
    """Apply consistent formatting to subplot."""
    if i == 0:
        # Use formatted variable name for title
        formatted_title = _get_variable_label(variable)
        ax.set_title(formatted_title, fontweight='bold', fontsize=10)
    if j == 0:
        ax.set_ylabel(strain, fontweight='bold', fontsize=11)
    if i == n_strains - 1:
        ax.set_xlabel('Julian day', fontsize=10)
    
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def _create_timeseries_legend(axes, obs_data, data):
    """Create comprehensive legend for multi-light regime time series plot."""
    available_lights = data.get('_meta', {}).get('available_lights', [])
    
    legend_handles = []
    legend_labels = []
    
    # Add common optimization legend (gray solid lines)
    if available_lights:
        legend_handles.append(plt.Line2D([0], [0], color=CO_COLOR, linestyle='-', 
                                       linewidth=3, alpha=0.8))
        legend_labels.append('Common Optimization')
    
    # Add light regime line style legends for IO
    for light in available_lights:
        legend_handles.append(plt.Line2D([0], [0], color='black', 
                                       linestyle=LIGHT_STYLES[light],
                                       linewidth=2, alpha=0.8))
        legend_labels.append(f'{LIGHT_LABELS[light]} (IO)')
    
    # Add strain color legends  
    for strain in STRAINS:
        legend_handles.append(plt.Line2D([0], [0], color=STRAIN_COLORS[strain], 
                                       linestyle='-', linewidth=3, alpha=0.8))
        legend_labels.append(f'{strain}')
    
    # Add observation marker legends for each light regime
    if obs_data:
        for light in available_lights:
            legend_handles.append(plt.Line2D([0], [0], marker=LIGHT_MARKERS[light], 
                                           color='w', markerfacecolor='black', 
                                           markersize=8, alpha=0.8, linestyle='None'))
            legend_labels.append(f'{LIGHT_LABELS[light]} obs.')
    
    return legend_handles, legend_labels


def create_timeseries_comparison_plot(variables, data, obs_data=None, figsize=None, save_fig=False):
    """
    Create a grid of subplots comparing time series between co and io optimizations.
    
    Parameters:
    -----------
    variables : list
        List of variable names to plot (columns)
    data : dict
        Simulation data loaded from load_simulation_data()
    obs_data : dict, optional
        Observational data loaded from load_observations()
    figsize : tuple, optional
        Figure size (width, height)
    save_fig : bool, optional
        Whether to save the figure
    
    Returns:
    --------
    tuple
        (figure, axes) matplotlib objects
    """
    n_strains = len(STRAINS)
    n_vars = len(variables)

    if figsize is None:
        figsize = (4 * n_vars, 3 * n_strains)

    fig, axes = plt.subplots(n_strains, n_vars, figsize=figsize, 
                            sharex=True, sharey=False)
    
    # Handle single row/column cases
    if n_strains == 1:
        axes = axes.reshape(1, -1)
    if n_vars == 1:
        axes = axes.reshape(-1, 1)
    
    # Plot data for each strain and variable
    for i, strain in enumerate(STRAINS):
        for j, var in enumerate(variables):
            ax = axes[i, j]
            
            _plot_simulation_lines(ax, data, strain, var, i, j)
            _plot_observations(ax, obs_data, strain, var, j)
            _format_subplot(ax, var, strain, i, j, n_strains)
    
    # Create and add legend with better positioning
    all_handles, all_labels = _create_timeseries_legend(axes, obs_data, data)
    
    # Calculate optimal number of columns for legend
    n_entries = len(all_labels)
    if n_entries <= 4:
        ncol = 2
        top_adjust = 0.91
    elif n_entries <= 8:
        ncol = 4
        top_adjust = 0.89
    else:
        ncol = min(6, n_entries)
        top_adjust = 0.87
    
    fig.legend(all_handles, all_labels, loc='upper center', bbox_to_anchor=(0.5, 0.98),
              ncol=ncol, frameon=True, fontsize=9, fancybox=True, shadow=True,
              edgecolor='gray', facecolor='white', framealpha=0.95)
    
    plt.tight_layout()
    plt.subplots_adjust(top=top_adjust)
    
    # Save figure if requested
    if save_fig:
        var_names = "_".join(variables[:3])  # Limit to first 3 variables for filename
        save_figure(fig, savetime=True, figname=f"timeseries_comparison_{var_names}.png")
    
    return fig, axes


def _collect_parameter_values(parameters, data):
    """Collect parameter values for boxplot visualization with multi-light support."""
    param_data = {}
    available_lights = data.get('_meta', {}).get('available_lights', [])
    
    for i, param in enumerate(parameters):
        param_values = []
        param_types = []
        param_strains = []
        param_lights = []
        
        # Collect independent optimization values for each strain and light regime
        for strain in STRAINS:
            for light in available_lights:
                if (light in data[strain]['io_params'] and 
                    param in data[strain]['io_params'][light]):
                    param_values.append(data[strain]['io_params'][light][param])
                    param_types.append('io')
                    param_strains.append(strain)
                    param_lights.append(light)
        
        # Collect common optimization values (same for all strains but may vary by light)
        for light in available_lights:
            if (light in data[STRAINS[0]]['co_params'] and 
                param in data[STRAINS[0]]['co_params'][light]):
                param_values.append(data[STRAINS[0]]['co_params'][light][param])
                param_types.append('co')
                param_strains.append('Common')
                param_lights.append(light)
        
        param_data[param] = {
            'values': param_values,
            'types': param_types,
            'strains': param_strains,
            'lights': param_lights,
            'position': i + 1
        }
    
    return param_data


def _create_parameter_legend():
    """Create legend for parameter comparison plot."""
    legend_elements = []

    for strain in STRAINS:
        legend_elements.append(plt.Line2D([0], [0], marker="o", color='w',
                                        markerfacecolor=STRAIN_COLORS[strain],
                                        markersize=8, alpha=0.8, linestyle='None',
                                        label=f'{strain} (IO)'))
    
    legend_elements.append(plt.Line2D([0], [0], color='w', marker="p",
                                    markerfacecolor=CO_COLOR, markersize=8,
                                    alpha=0.8, linestyle='None',
                                    label='Common (CO)'))
    for light in LIGHT_LABELS:
        legend_elements.append(plt.Line2D([0], [0], color='w', marker=LIGHT_MARKERS[light],
                                          markerfacecolor=CO_COLOR, markersize=8,
                                          alpha=0.8, linestyle='None',
                                          label=LIGHT_LABELS[light]))

    return legend_elements


def create_parameter_boxplot(parameters, data, figsize=None, save_fig=False):
    """
    Create boxplots comparing parameter values between co and io optimizations.
    
    Parameters:
    -----------
    parameters : list
        List of parameter names to plot
    data : dict
        Simulation data loaded from load_simulation_data()
    figsize : tuple, optional
        Figure size (width, height)
    save_fig : bool, optional
        Whether to save the figure
    
    Returns:
    --------
    tuple
        (figure, axes) matplotlib objects
    """
    n_params = len(parameters)
    
    if figsize is None:
        figsize = (4 * n_params, 6)
    
    fig, axes = plt.subplots(1, n_params, figsize=figsize, sharey=False)
    
    # Handle single parameter case
    if n_params == 1:
        axes = [axes]
    
    # Collect parameter values
    param_data = _collect_parameter_values(parameters, data)
    
    # Create each subplot with its own boxplot
    for i, param in enumerate(parameters):
        ax = axes[i]
        
        if param not in param_data:
            continue
            
        param_info = param_data[param]
        values = param_info['values']
        types = param_info['types']
        strains = param_info['strains']
        lights = param_info['lights']
        # Create boxplot for this parameter
        bp = ax.boxplot([values], positions=[1], widths=0.6,
                       patch_artist=True, showfliers=False)
        # Style boxplot
        for patch, median in zip(bp['boxes'], bp['medians']):
            plt.yticks(fontsize=12)
            patch.set_facecolor('lightgray')
            patch.set_alpha(0.5)
            median.set_color('black')


        # Add individual points
        for val, typ, strain, lights in zip(values, types, strains, lights):
            if typ == 'io':
                color = STRAIN_COLORS[strain]
                marker = LIGHT_MARKERS[lights]
                size = 60
            else:  # co
                color = CO_COLOR
                marker = "p"
                size = 80
            
            ax.scatter(1, val, color=color, s=size, alpha=0.8, 
                      marker=marker, edgecolor='white', linewidth=1, zorder=5)
        
        # Formatting for each subplot
        formatted_label = _get_parameter_label(param)
        ax.set_title(formatted_label, fontsize=14)
        ax.set_xticks([])  # Remove x-axis ticks
        ax.set_ylabel('Value', fontweight='bold', fontsize=10)
        
        ax.grid(True, alpha=0.3, axis='y')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    
    # Add shared legend
    legend_elements = _create_parameter_legend()
    fig.legend(legend_elements, [elem.get_label() for elem in legend_elements],
              loc='upper center', bbox_to_anchor=(0.5, 0.98),
              ncol=min(6, len(legend_elements)), frameon=True, fontsize=9,
              fancybox=True, shadow=True, edgecolor='gray', facecolor='white', framealpha=0.95)
    
    # # Overall title
    # fig.suptitle('Parameter Comparison: Common vs Independent (for each strain and light regime) Optimization',
    #             fontweight='bold', fontsize=14, y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.75)

    
    # Save figure if requested
    if save_fig:
        param_names = "_".join([p.replace('+', '_') for p in parameters[:3]])
        save_figure(fig, savetime=True, figname=f"parameter_comparison_{param_names}.png")
    
    return fig, axes


def main():
    """Main function to demonstrate the plotting functions."""
    
    # Load data
    print("Loading simulation data...")
    data = load_simulation_data()
    
    print("Loading observational data...")
    obs_data = load_observations()
    
    # Example variables for time series comparison (including diagnostic variables)
    example_variables = ['Phy_Chl', 'TEPC_C', "BacA_C"]
    # example_variables = ['TEPC_C', 'Phy_Chl', 'TEP_Chl_ratio', "TEP_PhyC_ratio", "Chl_PhyC_ratio", "BacA_C"]

    f_factor = 2

    print(f"Creating time series comparison plot for variables: {example_variables}")
    print("Note: 'QN' and 'CN' are diagnostic variables that will be calculated from existing data")
    fig1, axes1 = create_timeseries_comparison_plot(example_variables, data, obs_data, save_fig=True,  figsize=(7* f_factor, 5 * f_factor))
    fig1.suptitle('Time Series Comparison: Common vs Independent (for each strain and light regime) Optimization',
                  fontsize=16, fontweight='bold', y=0.99)
    plt.show()


    # # Example parameters for boxplot comparison
    # example_parameters = ['Phy+mu_max', 'Phy+alpha', 'Phy+mortrate']
    #
    # print(f"Creating parameter boxplot for parameters: {example_parameters}")
    # fig2, ax2 = create_parameter_boxplot(example_parameters, data, save_fig=True,figsize=(4.3* f_factor, 1.6 * f_factor))
    #
    # plt.show()


    print("Visualization completed!")


if __name__ == "__main__":
    main()