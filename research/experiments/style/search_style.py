"""The house style of the search case studies' figures (research/experiments/*-search-*).

Warm paper, plum ink, serif type, hairlines, and a family of softened colours for the algorithms, so that every
figure of the gallery's search section reads as one family. A figure is 1600 x 1000 pixels, saved as JPEG into
assets/gallery/ as search-NN-subject.jpg.

Use:
    import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
    from search_style import *
    fig, ax = figure()
    ...
    save(fig, "search-01-subject")
"""

from __future__ import annotations

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[3]
GALLERY = ROOT / "assets" / "gallery"

PAPER = "#FBF6F1"   # warm paper
PANEL = "#F6EEE8"   # a panel on the paper
INK = "#3B2F3A"     # plum graphite
MUTED = "#8A7D86"   # secondary ink
HAIRLINE = "#D9CCD3"

# One colour per algorithm family, the same in every figure.
ALGORITHM = {
    "random": "#9A9097",        # graphite: the baseline
    "(1+1)-ES": "#C98B6B",      # sand
    "CMA-ES": "#8E3B78",        # deep orchid
    "IPOP-CMA-ES": "#B5629E",   # orchid
    "BIPOP-CMA-ES": "#CA8DBE",  # dusty orchid
    "DE": "#4E8F8B",            # teal
    "SHADE": "#2F6E6A",         # deep teal
    "PSO": "#D49A3A",           # amber
    "CE": "#6C8EBF",            # dusty blue
    "BO": "#D9738B",            # rose
    "GA": "#7FA37A",            # sage
    "(1+1) EA": "#A67C52",      # umber
    "ACO": "#C9A227",           # ochre
    "SA": "#5C7FA3",            # slate blue
    "ILS": "#6B3E75",           # plum
    "NSGA-II": "#8E3B78",
    "NSGA-III": "#4E8F8B",
    "MOEA/D": "#D49A3A",
    "SMS-EMOA": "#6C8EBF",
    "MAP-Elites": "#8E3B78",
    "CVT-MAP-Elites": "#4E8F8B",
    "novelty": "#D49A3A",
    "GP": "#8E3B78",
    "GP untyped": "#9A9097",
    "Grover": "#8E3B78",
    "QIEA": "#CA8DBE",
    "PBIL": "#4E8F8B",
}
CYCLE = ["#8E3B78", "#4E8F8B", "#D49A3A", "#6C8EBF", "#D9738B", "#7FA37A", "#C98B6B", "#6B3E75", "#9A9097"]

# Sequential maps from the paper: orchid (values), and a cool-to-warm diverging map for signed quantities.
ORCHID = LinearSegmentedColormap.from_list("orchid", ["#FBF6F1", "#EBD3E3", "#CA8DBE", "#8E3B78", "#3B1F3A"])
DUSK = LinearSegmentedColormap.from_list("dusk", ["#1F2A44", "#4E6E9A", "#A9B9CF", "#FBF6F1", "#E7B9A1", "#C2665A", "#6B2333"])
TERRAIN = LinearSegmentedColormap.from_list("terrain_warm", ["#3B2F3A", "#6B3E75", "#B5629E", "#E7B9A1", "#F6EEE8"])


def apply() -> None:
    matplotlib.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Source Serif 4", "Source Serif Pro", "Georgia", "Cambria", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "font.size": 13,
        "axes.titlesize": 15,
        "axes.titleweight": "regular",
        "axes.labelsize": 13,
        "axes.edgecolor": HAIRLINE,
        "axes.linewidth": 0.8,
        "axes.labelcolor": INK,
        "axes.facecolor": PAPER,
        "axes.prop_cycle": matplotlib.cycler(color=CYCLE),
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": HAIRLINE,
        "grid.linewidth": 0.5,
        "grid.alpha": 0.7,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "text.color": INK,
        "figure.facecolor": PAPER,
        "savefig.facecolor": PAPER,
        "legend.frameon": False,
        "legend.fontsize": 11,
        "lines.linewidth": 1.8,
    })


apply()


def figure(rows: int = 1, cols: int = 1, **kw):
    """A 1600 x 1000 figure (16 x 10 inches at 100 dpi)."""
    fig, ax = plt.subplots(rows, cols, figsize=(16, 10), dpi=100, **kw)
    return fig, ax


def title(fig, text: str, subtitle: str | None = None) -> None:
    """A title at the top left, in the gallery's manner, with an optional line beneath."""
    fig.text(0.04, 0.955, text, fontsize=22, color=INK, ha="left", va="top")
    if subtitle:
        fig.text(0.04, 0.915, subtitle, fontsize=13, color=MUTED, ha="left", va="top")


def source(fig, text: str) -> None:
    """The figure's provenance at the bottom left: what was run, on what, how."""
    fig.text(0.04, 0.02, text, fontsize=10, color=MUTED, ha="left", va="bottom")


def save(fig, name: str) -> pathlib.Path:
    """Saves assets/gallery/<name>.jpg at 1600 x 1000."""
    GALLERY.mkdir(parents=True, exist_ok=True)
    path = GALLERY / f"{name}.jpg"
    fig.savefig(path, dpi=100, pil_kwargs={"quality": 90, "optimize": True})
    plt.close(fig)
    return path
