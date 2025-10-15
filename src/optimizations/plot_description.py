"""
Plotting utilities for multi-description optimization results.

This module provides visualization functions for optimizations with multiple descriptions.
Uses existing plotting utilities from src.utils.plotting for consistency.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Optional, Tuple

from src.utils import plotting
from src.config_model import varinfos
from src.utils import functions as fns

# Style configuration
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")


def plot_timeseries_comparison(opt,
                               variables: List[str],
                               ncols: int = 3,
                               figsize: Optional[Tuple[float, float]] = None,
                               savefig: bool = False,
                               figname: Optional[str] = None):
    """
    Create time series comparison plots across descriptions using daily means.

    Args:
        opt: Optimization object (must have processed results)
        variables: List of variable names to plot
        ncols: Number of columns for subplot grid
        figsize: Figure size (width, height)
        savefig: Whether to save the figure
        figname: Custom figure name for saving

    Returns:
        Tuple of (figure, axes)
    """
    if not opt.is_optimization_processed:
        raise ValueError("Optimization must be processed before plotting. Call opt.process_results() first.")

    n_descriptions = len(opt.descriptions)
    n_vars = len(variables)

    # Calculate layout
    nrows = n_descriptions
    actual_ncols = min(ncols, n_vars)

    if figsize is None:
        figsize = (5 * actual_ncols, 3 * nrows)

    # Create subplots
    fig, axes = plt.subplots(nrows, actual_ncols, figsize=figsize, squeeze=False, sharex=True)

    # Color palette for descriptions
    colors = sns.color_palette("husl", n_descriptions)

    # Plot each description
    for i, description in enumerate(opt.descriptions):
        desc_color = colors[i]
        desc_name = description.name

        # Skip if no best model
        if description.best_model is None:
            continue

        # Prepare model and observation data using plotting.py utilities
        model_data_list, merged_data, _, _ = plotting.prepare_model_obs_data(
            models=description.best_model,
            observations=description.observation,
            daily_mean=True,
            variables_to_plot=variables
        )

        # Get the single model DataFrame (list has only one element)
        model_data = model_data_list[0]

        # Plot each variable
        for j, var in enumerate(variables):
            if j >= actual_ncols:
                break

            ax = axes[i, j]

            # Plot model line
            if var in model_data.columns:
                ax.plot(model_data.index, model_data[var],
                       color=desc_color, linewidth=2, alpha=0.8,
                       label=desc_name)
            else:
                # Variable not available
                ax.text(0.5, 0.5, f'{var}\nNot Available',
                       ha='center', va='center', transform=ax.transAxes)
                continue

            # Plot observations if available
            if merged_data is not None and f'{var}_OBS' in merged_data.columns:
                ax.scatter(merged_data.index, merged_data[f'{var}_OBS'],
                          color=desc_color, s=50, alpha=0.8,
                          edgecolor='white', linewidth=1.5, zorder=5)

            # Formatting using varinfos
            var_info = varinfos.doutput.get(var.lstrip('m'), {})
            clean_name = var_info.get('cleanname', var)
            units = var_info.get('munits' if var.startswith('m') else 'units', '')

            if i == 0:
                title = f"{fns.cleantext(clean_name)}"
                if units:
                    title += f"\n[{fns.cleantext(units)}]"
                ax.set_title(title, fontweight='bold', fontsize=10)

            if j == 0:
                ax.set_ylabel(desc_name, fontweight='bold', fontsize=9)

            if i == nrows - 1:
                ax.set_xlabel('Julian Day', fontsize=9)

            ax.grid(True, alpha=0.3, linestyle=':')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

    # Overall title
    fig.suptitle(f'Time Series Comparison - {opt.name}\n'
                f'{n_descriptions} descriptions • {n_vars} variables',
                fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout()
    plt.subplots_adjust(top=0.93)

    if savefig:
        vars_str = "_".join(variables[:3])
        name = figname or f"{opt.name}_timeseries_{vars_str}.png"
        plotting.save_figure(fig, filename=name, add_timestamp=True)

    return fig, axes