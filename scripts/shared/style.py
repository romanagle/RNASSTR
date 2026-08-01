# shared/style.py

"""
Centralized plotting style for the RNASSTR revision pipeline.

Implements:
- Feng Zhang-inspired palette
- Arial font
- SVG-friendly settings
- KDE styling
- consistent axis formatting
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

# ============================================================
# Feng-style color palette
# ============================================================

FENG_PINK = "#FD4078"
FENG_GREY = "#6E6E6E"
FENG_BLUE = "#51ADE0"
FENG_GOLD = "#E6A53A"

SPLIT_COLORS = {
    "train": FENG_GREY,
    "val": FENG_BLUE,
    "test": FENG_PINK,
}

# ============================================================
# Matplotlib global settings
# ============================================================

mpl.rcParams["font.family"] = "Arial"

mpl.rcParams["svg.fonttype"] = "none"

mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False

mpl.rcParams["axes.linewidth"] = 1.2

mpl.rcParams["xtick.major.width"] = 1.2
mpl.rcParams["ytick.major.width"] = 1.2

mpl.rcParams["xtick.major.size"] = 5
mpl.rcParams["ytick.major.size"] = 5

mpl.rcParams["legend.frameon"] = False

mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["savefig.dpi"] = 300

mpl.rcParams["savefig.bbox"] = "tight"

# ============================================================
# Common plotting helpers
# ============================================================

KDE_ALPHA = 0.25
LINEWIDTH = 2.0

GRID_ALPHA = 0.2

DEFAULT_FIGSIZE = (5, 5)

# ============================================================
# Helper functions
# ============================================================

def setup_axis(ax):
    """
    Apply consistent axis formatting.
    """
    ax.grid(True, alpha=GRID_ALPHA)
    ax.tick_params(axis='both', labelsize=11)

    return ax


def save_figure(fig, out_prefix):
    """
    Save figure as both SVG and PNG.
    """
    fig.savefig(f"{out_prefix}.svg")
    fig.savefig(f"{out_prefix}.png")
