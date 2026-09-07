#!/usr/bin/env python3

"""
Shared plotting functions for the uncertainty and power analyses.

This module provides the common Matplotlib configuration used to generate
all uncertainty and statistical power figures from:
"When ten cases are not enough: Uncertainty and statistical power in AI
screening evaluations"

The shared configuration ensures consistent figure dimensions, typography,
axis formatting, grid lines, and reference lines across all figures.

Author: Gwenolé Quellec
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# ============================================================================
# Configuration
# ============================================================================

FIGURE_SIZE = (5, 5)

FONT_FAMILY = "sans-serif"
FONT_SANS_SERIF = ("Arial", "DejaVu Sans")
FONT_SIZE = 10
AXES_TITLE_SIZE = 11
AXES_LABEL_SIZE = 10
LEGEND_FONT_SIZE = 9

Y_MIN = 0.0
Y_MAX = 1.0
Y_MAJOR_TICK = 0.1
Y_MINOR_TICK = 0.05

REFERENCE_SAMPLE_SIZE = 10

GRID_ALPHA = 0.3
REFERENCE_LINE_ALPHA = 0.6

LINE_STYLE = "--"
REFERENCE_LINE_COLOR = "black"


# ============================================================================
# Plotting functions
# ============================================================================

def setup_figure(
    xmax: float,
    xticks: tuple,
):
    """
    Configure the common Matplotlib figure style.

    Parameters
    ----------
    xmax : float
        Maximum value of the horizontal axis.
    xticks : tuple
        Major and minor tick intervals for the horizontal axis.
    """

    mpl.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.sans-serif": FONT_SANS_SERIF,
        "font.size": FONT_SIZE,
        "axes.titlesize": AXES_TITLE_SIZE,
        "axes.labelsize": AXES_LABEL_SIZE,
        "legend.fontsize": LEGEND_FONT_SIZE,
    })

    plt.figure(figsize=FIGURE_SIZE)

    ax = plt.gca()

    ax.set_xlim(0, xmax)
    ax.set_ylim(Y_MIN, Y_MAX)

    ax.axvline(
        REFERENCE_SAMPLE_SIZE,
        linestyle=LINE_STYLE,
        color=REFERENCE_LINE_COLOR,
        alpha=REFERENCE_LINE_ALPHA,
    )

    ax.xaxis.set_major_locator(
        ticker.MultipleLocator(xticks[0])
    )
    ax.xaxis.set_minor_locator(
        ticker.MultipleLocator(xticks[1])
    )
    ax.yaxis.set_major_locator(
        ticker.MultipleLocator(Y_MAJOR_TICK)
    )
    ax.yaxis.set_minor_locator(
        ticker.MultipleLocator(Y_MINOR_TICK)
    )

    ax.grid(
        which="both",
        linestyle=LINE_STYLE,
        alpha=GRID_ALPHA,
    )


def plot_curves(
    curves: dict,
    xlabel: str,
    ylabel: str,
    filename: str,
    xmax: float,
    xticks: tuple,
    hline: float | None = None,
):
    """
    Plot a collection of curves using the common figure style.

    Parameters
    ----------
    curves : dict
        Mapping from curve labels to (x, y) data pairs.
    xlabel : str
        Label of the horizontal axis.
    ylabel : str
        Label of the vertical axis.
    filename : str
        Output figure filename.
    xmax : float
        Maximum value of the horizontal axis.
    xticks : tuple
        Major and minor tick intervals for the horizontal axis.
    hline : float or None, optional
        Horizontal reference line. No line is drawn if None.
    """

    setup_figure(
        xmax=xmax,
        xticks=xticks,
    )

    for label, (x, y) in curves.items():
        plt.plot(
            x,
            y,
            label=label,
        )

    if hline is not None:
        plt.axhline(
            hline,
            linestyle=LINE_STYLE,
            color=REFERENCE_LINE_COLOR,
            alpha=REFERENCE_LINE_ALPHA,
        )

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()

    plt.tight_layout()
    plt.savefig(filename)
