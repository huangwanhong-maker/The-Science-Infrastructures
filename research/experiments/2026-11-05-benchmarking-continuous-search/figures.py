"""Figures of the case study 'benchmarking continuous search', drawn from results/ only.

    .venv/Scripts/python.exe research/experiments/2026-11-05-benchmarking-continuous-search/figures.py

writes assets/gallery/search-06-ecdf-by-dimension.jpg (runtime distributions per dimension, all functions pooled),
search-07-ecdf-by-function.jpg (per function at n = 10), search-08-ert.jpg (expected running time against dimension)
and search-09-solved.jpg (which algorithm solves which function in which dimension).
"""

from __future__ import annotations

import csv
import math
import pathlib
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F403
import matplotlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# The house serif, with fallbacks for the glyphs it lacks (the superscript minus).
_have = {f.name for f in font_manager.fontManager.ttflist}
matplotlib.rcParams["font.family"] = [f for f in ["Source Serif 4", "Source Serif Pro", "Georgia", "Cambria", "DejaVu Serif"] if f in _have]

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
RECORD = "research/experiments/2026-11-05-benchmarking-continuous-search"
VERSION = "6b3cc88"

ALGORITHMS = ["random", "(1+1)-ES", "CMA-ES", "IPOP-CMA-ES", "BIPOP-CMA-ES", "DE", "SHADE", "PSO", "SPSO 2011", "CE"]
FUNCTIONS = ["sphere", "ellipsoid", "rotated ellipsoid", "rosenbrock", "rastrigin", "ackley", "griewank", "schwefel"]
TITLES = {"sphere": "Sphere", "ellipsoid": "Ellipsoid (separable)", "rotated ellipsoid": "Rotated ellipsoid", "rosenbrock": "Rosenbrock",
          "rastrigin": "Rastrigin", "ackley": "Ackley", "griewank": "Griewank", "schwefel": "Schwefel (separable)"}
SHORT = {"sphere": "sphere", "ellipsoid": "ellipsoid", "rotated ellipsoid": "rot. ellipsoid", "rosenbrock": "Rosenbrock",
         "rastrigin": "Rastrigin", "ackley": "Ackley", "griewank": "Griewank", "schwefel": "Schwefel"}
DIMENSIONS = [2, 5, 10, 20]

# The swarm of SPSO 2011 shares PSO's colour, dashed; every other colour is the house's.
COLOUR = dict(ALGORITHM)
COLOUR["SPSO 2011"] = ALGORITHM["PSO"]
STYLE = {a: "-" for a in ALGORITHMS}
STYLE["SPSO 2011"] = (0, (4, 2))
WIDTH = {a: 1.8 for a in ALGORITHMS}
WIDTH["random"] = 1.4
SUITE = "8 functions × 15 instances (shifted, rotated) × 11 targets f − f* ∈ {10², …, 10⁻⁸}"


def read_csv(name):
    with open(RESULTS / name, newline="") as fh:
        return list(csv.DictReader(fh))


def num(s):
    return float(s) if s not in ("", None) else math.nan


def sup(e):
    return str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))


def ecdfs():
    data = defaultdict(list)
    for r in read_csv("ecdf.csv"):
        data[(r["function"], int(r["n"]), r["algorithm"])].append((num(r["budget_per_n"]), num(r["fraction"])))
    return {k: np.array(v) for k, v in data.items()}


def spread_labels(ax, items, x, gap=0.04, lo=0.0, hi=1.0, size=11):
    """Labels at the ends of curves on a linear axis, at least `gap` apart."""
    items = sorted(items, key=lambda t: t[1])
    ys = [max(lo, y) for _, y, _ in items]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + gap)
    if ys and ys[-1] > hi:
        ys[-1] = hi
        for i in range(len(ys) - 2, -1, -1):
            ys[i] = min(ys[i], ys[i + 1] - gap)
    for (text, _, colour), y in zip(items, ys):
        ax.text(x, y, text, color=colour, fontsize=size, va="center", ha="left")


def ecdf_axis(ax, right=10 ** 5.6):
    ax.set_xscale("log")
    ax.set_xlim(1, right)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xticks([1, 10, 100, 1000, 10000])
    ax.set_xticklabels(["1", "10", "10²", "10³", "10⁴"])
    ax.axvline(1e4, color=HAIRLINE, lw=0.8)
    ax.grid(False, axis="x")
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])


def legend_row(fig, y, x0=0.04, size=12.5):
    """The algorithms' names in their colours, in one line, in place of a boxed legend."""
    x = x0
    renderer = fig.canvas.get_renderer()
    for a in ALGORITHMS:
        t = fig.text(x, y, a, color=COLOUR[a], fontsize=size, ha="left", va="center")
        if a == "SPSO 2011":
            t.set_text("SPSO 2011 (dashed)")
        w = t.get_window_extent(renderer=renderer).width / fig.bbox.width
        x += w + 0.018


def figure_ecdf_by_dimension(E):
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), dpi=100)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.83, bottom=0.1, wspace=0.16, hspace=0.3)
    for ax, n in zip(axes.flat, DIMENSIONS):
        ends = []
        for a in ALGORITHMS:
            d = E[("all", n, a)]
            ax.step(d[:, 0], d[:, 1], where="post", color=COLOUR[a], ls=STYLE[a], lw=WIDTH[a])
            ends.append((a, d[-1, 1], COLOUR[a]))
        ecdf_axis(ax)
        spread_labels(ax, ends, 1e4 * 1.25, gap=0.047, hi=1.05, size=10.5)
        ax.set_title(f"n = {n}", loc="left", fontsize=16, color=INK)
        ax.set_ylabel("fraction of (run, target) pairs")
    for ax in axes[1]:
        ax.set_xlabel("evaluations / n")
    title(fig, "Runtime distributions, by dimension",
          "For each budget, the fraction of (run, target) pairs solved: 8 functions, 15 instances each, 11 targets f − f* from 10² to 10⁻⁸, 1320 pairs per algorithm "
          "and dimension.\nThe budget is 10⁴·n evaluations (the hairline). A curve that ends low has failed on some functions or targets for good.")
    source(fig, f"heuresis {VERSION} · {RECORD} · empirical runtime distributions by heuresis::analysis::runtime_distribution, as COCO's · from results/ecdf.csv")
    return save(fig, "search-06-ecdf-by-dimension")


def figure_ecdf_by_function(E, n=10):
    fig, axes = plt.subplots(2, 4, figsize=(16, 10), dpi=100, sharey=True)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.79, bottom=0.1, wspace=0.08, hspace=0.32)
    for ax, f in zip(axes.flat, FUNCTIONS):
        for a in ALGORITHMS:
            d = E[(f, n, a)]
            ax.step(d[:, 0], d[:, 1], where="post", color=COLOUR[a], ls=STYLE[a], lw=WIDTH[a] * 0.9)
        ecdf_axis(ax, right=1.6e4)
        ax.set_title(TITLES[f], loc="left", fontsize=14, color=INK)
    for ax in axes[:, 0]:
        ax.set_ylabel("fraction of (run, target) pairs")
    for ax in axes[1]:
        ax.set_xlabel("evaluations / n")
    title(fig, f"Runtime distributions on each function, n = {n}",
          "15 instances × 11 targets f − f* from 10² to 10⁻⁸ for each algorithm and function; the budget is 10⁵ evaluations.")
    legend_row(fig, 0.855)
    source(fig, f"heuresis {VERSION} · {RECORD} · heuresis::analysis::runtime_distribution · from results/ecdf.csv")
    return save(fig, "search-07-ecdf-by-function")


def figure_ert(target="1e-08"):
    rows = read_csv("ert.csv")
    E = {}
    for r in rows:
        if r["target"] == target:
            E[(r["function"], int(r["n"]), r["algorithm"])] = (num(r["ert"]), int(r["successes"]))
    fig, axes = plt.subplots(2, 4, figsize=(16, 10), dpi=100, sharey=True)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.79, bottom=0.1, wspace=0.08, hspace=0.32)
    for ax, f in zip(axes.flat, FUNCTIONS):
        for a in ALGORITHMS:
            ys = [E[(f, n, a)][0] / n for n in DIMENSIONS]  # +∞ where nothing was solved: the line breaks there
            if not any(math.isfinite(y) for y in ys):
                continue
            ax.plot(DIMENSIONS, [y if math.isfinite(y) else np.nan for y in ys], color=COLOUR[a], ls=STYLE[a], lw=WIDTH[a] * 0.9)
            for n, y in zip(DIMENSIONS, ys):
                if math.isfinite(y):
                    full = E[(f, n, a)][1] == 15
                    ax.scatter([n], [y], s=26, color=COLOUR[a] if full else PAPER, edgecolor=COLOUR[a], lw=1.2, zorder=4)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(DIMENSIONS)
        ax.set_xticklabels([str(n) for n in DIMENSIONS])
        ax.minorticks_off()
        ax.set_xlim(1.7, 24)
        ax.set_ylim(10, 3e6)
        ax.axhline(1e4, color=MUTED, lw=0.7, ls=(0, (2, 3)))
        ax.grid(False, axis="x")
        ax.set_title(TITLES[f], loc="left", fontsize=14, color=INK)
    for ax in axes[:, 0]:
        ax.set_ylabel("ERT / n to f − f* ≤ 10⁻⁸")
    for ax in axes[1]:
        ax.set_xlabel("dimension n")
    axes[0, 1].text(2.0, 1.25e4, "the budget, 10⁴·n", color=MUTED, fontsize=10, va="bottom")
    title(fig, "Expected running time against dimension",
          "Evaluations per dimension that independent restarts would need to reach f − f* ≤ 10⁻⁸ (ERT: all evaluations spent over all runs, divided by the successes),"
          "\n15 instances each. A filled point: every instance solved; an open one: some failed, so ERT counts their full budgets. No point: no instance solved.")
    legend_row(fig, 0.855)
    source(fig, f"heuresis {VERSION} · {RECORD} · heuresis::analysis::expected_running_time · from results/ert.csv")
    return save(fig, "search-08-ert")


def figure_solved(target="1e-08"):
    ert = {}
    for r in read_csv("ert.csv"):
        if r["target"] == target:
            ert[(r["function"], int(r["n"]), r["algorithm"])] = int(r["successes"])
    gaps = defaultdict(list)
    for r in read_csv("runs.csv"):
        gaps[(r["function"], int(r["n"]), r["algorithm"])].append(num(r["gap"]))
    fig = plt.figure(figsize=(16, 10), dpi=100)
    ax = fig.add_axes([0.13, 0.12, 0.85, 0.66])
    cmap = LinearSegmentedColormap.from_list("solved", [PANEL, "#EBD3E3", "#CA8DBE", "#8E3B78"])
    cols = [(f, n) for f in FUNCTIONS for n in DIMENSIONS]
    M = np.array([[ert[(f, n, a)] / 15 for (f, n) in cols] for a in ALGORITHMS])
    ax.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    for i, a in enumerate(ALGORITHMS):
        for j, (f, n) in enumerate(cols):
            k = ert[(f, n, a)]
            if k > 0:
                ax.text(j, i, str(k), ha="center", va="center", fontsize=11, color=PAPER if k >= 10 else INK)
            else:
                g = np.nanmedian(gaps[(f, n, a)])
                e = math.floor(math.log10(g)) if g > 0 else None
                ax.text(j, i, f"10{sup(e)}" if e is not None else "0", ha="center", va="center", fontsize=9, color=MUTED)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([str(n) for f, n in cols], fontsize=10, color=MUTED)
    ax.set_yticks(range(len(ALGORITHMS)))
    ax.set_yticklabels(ALGORITHMS, fontsize=12.5)
    for t, a in zip(ax.get_yticklabels(), ALGORITHMS):
        t.set_color(COLOUR[a])
    ax.tick_params(length=0)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    for j in range(len(FUNCTIONS) + 1):
        ax.axvline(j * 4 - 0.5, color=PAPER, lw=4)
    for i in range(len(ALGORITHMS) + 1):
        ax.axhline(i - 0.5, color=PAPER, lw=1.5)
    for j, f in enumerate(FUNCTIONS):
        ax.text(j * 4 + 1.5, -0.85, SHORT[f], ha="center", va="bottom", fontsize=13, color=INK)
    ax.text(-0.75, len(ALGORITHMS) - 0.5 + 0.2, "n", ha="right", va="top", fontsize=11, color=MUTED, style="italic")
    title(fig, "Which algorithm solves which function",
          "In each cell, how many of 15 instances reached f − f* ≤ 10⁻⁸ within 10⁴·n evaluations, shaded from none to all. Where none did, the median gap f − f*"
          "\nthat the runs reached instead, as a power of ten. Columns: the functions, each in 2, 5, 10 and 20 dimensions.")
    source(fig, f"heuresis {VERSION} · {RECORD} · successes and gaps from results/ert.csv and results/runs.csv")
    return save(fig, "search-09-solved")


if __name__ == "__main__":
    E = ecdfs()
    print(figure_ecdf_by_dimension(E))
    print(figure_ecdf_by_function(E))
    print(figure_ert())
    print(figure_solved())
