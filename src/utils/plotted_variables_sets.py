"""
Variable set definitions for plotting.

This module contains the PlottedVariablesSet class and related utilities.
Separated from plotting.py to avoid circular imports.
"""

from typing import List, Optional, Tuple


class PlottedVariablesSet:
    """
    A container for a set of variables to be plotted together with their display preferences.

    This class encapsulates not just the list of variables, but also how they should be displayed:
    layout (ncols), sizing (figsize), naming (for automatic file saving), and positioning (legend).

    Provides backward compatibility by implementing list-like behavior through __iter__, __len__,
    and __getitem__, so existing code using list(variables) continues to work unchanged.
    """

    def __init__(
        self,
        name: str,
        variables: List[str],
        ncols: int = 2,
        figsize: Optional[Tuple[int, int]] = None,
        legend_position: str = 'center right',
        description: str = ""
    ):
        """
        Initialize a set of variables for plotting.

        Args:
            name: Identifier for this variable set (used for automatic file naming)
            variables: List of variable names to plot
            ncols: Number of columns in the subplot grid
            figsize: Figure size tuple (width, height). If None, auto-calculated
            legend_position: Position of the legend ('center right', 'upper right', etc.)
            description: Human-readable description of this variable set
        """
        self.name = name
        self.variables = variables
        self.ncols = ncols
        self.figsize = figsize
        self.legend_position = legend_position
        self.description = description

    def __iter__(self):
        """Enable iteration: for var in variable_set"""
        return iter(self.variables)

    def __len__(self):
        """Enable len(): len(variable_set)"""
        return len(self.variables)

    def __getitem__(self, key):
        """Enable indexing: variable_set[0], variable_set[1:3], etc."""
        return self.variables[key]

    def __repr__(self):
        return (f"PlottedVariablesSet(name='{self.name}', "
                f"variables={len(self.variables)} vars, ncols={self.ncols})")

    def __str__(self):
        return f"{self.name}: {len(self.variables)} variables in {self.ncols} columns"

    @property
    def nrows(self):
        """Calculate number of rows needed for the subplot grid."""
        return (len(self.variables) + self.ncols - 1) // self.ncols