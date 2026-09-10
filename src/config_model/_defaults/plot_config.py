"""Fallback figure geometry. See src/config_model/_defaults/__init__.py.

A fixed-margin grid: the axes box has the same size in every figure, and the figure size
follows from the composition. Everything is in inches, so the result is physical and a
figure saved at its exact size needs no `bbox_inches='tight'`.

Only what `src/utils/plotting.py` (`paper_grid`, `default_observations`) needs to work
without a user-supplied module.
"""

# Observation dataset used by plotting.plot_results when a call does not name one
# (station = a file name in Observations/). None = no default observations.
DEFAULT_OBS_STATION = None

PANEL_SIZE = (1.85, 1.35)   # (width, height) of one axes box, in inches

FIG_MARGINS = {
    'left': 0.57,    # ylabel + tick labels + pads
    'right': 0.15,   # overhang of the last x tick label
    'top': 0.15,
    'bottom': 0.25,  # unrotated date tick labels
    'wgap': 0.46,    # between columns: must fit the next column's ylabel
    'hgap': 0.10,    # between rows (shared x axis, so no intermediate labels)
}


def paper_figsize(nrows, ncols=1, panel=PANEL_SIZE, **margins):
    """Figure size (inches) derived from the composition -- see plotting.paper_grid."""
    m = {**FIG_MARGINS, **margins}
    pw, ph = panel
    return (m['left'] + ncols * pw + (ncols - 1) * m['wgap'] + m['right'],
            m['top'] + nrows * ph + (nrows - 1) * m['hgap'] + m['bottom'])
