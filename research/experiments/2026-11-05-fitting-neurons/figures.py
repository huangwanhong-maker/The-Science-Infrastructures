"""Figures of fitting neurons and the repertoire of their firing (README.md §7), drawn from results/ only.

    .venv/Scripts/python.exe research/experiments/2026-11-05-fitting-neurons/figures.py
"""

import csv
import json
import math
import pathlib
import sys
from collections import defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F403

import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.colors import ListedColormap, BoundaryNorm  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
VERSION = "heuresis 6b3cc88 · neuroscience 585f4d7 · root 4153cb9"
ALGOS = ["CMA-ES", "DE", "PSO", "random"]
STYLE = {"CMA-ES": "-", "DE": "-", "PSO": (0, (6, 2)), "random": (0, (2, 2))}
TRUTH = {"a_nS": 2.0, "tau_w_ms": 300.0, "b_pA": 60.0, "V_r_mV": -58.0}
FLOOR = 1e-7  # where an exact fit (error 0) is drawn on a logarithmic axis


def rows(name):
    with open(RESULTS / name, newline="") as f:
        return list(csv.DictReader(f))


def panel_title(ax, text, sub=None):
    ax.set_title(text, loc="left", fontsize=15, color=INK, pad=24 if sub else 10)
    if sub:
        ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=11, color=MUTED, ha="left", va="bottom")


def steps(progress, budget):
    """The best-so-far as a step function on 1 … budget, from the improvements (evaluation, value)."""
    x = np.arange(1, budget + 1)
    y = np.full(budget, np.nan)
    for e, v in progress:
        y[int(e) - 1:] = max(v, FLOOR)
    return x, y


# ---- search-40: convergence on the two objectives ------------------------------------------------------------

def convergence_figure():
    conv = rows("convergence.csv")
    fits = rows("fits.csv")
    by = defaultdict(list)
    for r in conv:
        by[(r["objective"], r["algorithm"], int(r["seed"]))].append((float(r["evaluations"]), float(r["best"])))
    budget = 1000
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Fitting an adaptive neuron to its own behaviour, by four searches",
          "AdEx with four hidden parameters (a, τw, b, Vr); five runs of 1000 evaluations each, the best error found so far; thin lines single runs, thick the median")
    gs = fig.add_gridspec(1, 2, left=0.065, right=0.97, top=0.80, bottom=0.12, wspace=0.18)
    for col, (objective, label, sub) in enumerate((("features", "The rates: f–I curves and their adaptation",
                                                    "mean of log²((f + 1)/(f* + 1)) over initial and steady rates at six currents"),
                                                   ("spikes", "The spike times", "mean over six currents of van Rossum's squared distance (τ = 20 ms), relative to the target's norm"))):
        ax = fig.add_subplot(gs[0, col])
        for alg in ALGOS:
            colour = ALGORITHM[alg]
            curves = []
            for seed in range(1, 6):
                x, y = steps(sorted(by[(objective, alg, seed)]), budget)
                curves.append(y)
                ax.plot(x, y, color=colour, linewidth=0.8, alpha=0.35, linestyle=STYLE[alg])
            med = np.nanmedian(np.vstack(curves), axis=0)
            ax.plot(np.arange(1, budget + 1), med, color=colour, linewidth=2.4, linestyle=STYLE[alg], label=alg)
        if objective == "spikes":
            # The two-stage fit: CMA-ES on the spike times from each CMA-ES fit of the rates, after its evaluations.
            stage1 = {int(r["seed"]): float(r["evaluations"]) for r in fits if r["objective"] == "features" and r["algorithm"] == "CMA-ES"}
            rc = defaultdict(list)
            for r in rows("refine_convergence.csv"):
                rc[int(r["seed"])].append((float(r["evaluations"]), float(r["best"])))
            for seed, pr in rc.items():
                pr.sort()
                x = [stage1[seed] + e for e, _ in pr] + [stage1[seed] + 600]
                y = [v for _, v in pr] + [pr[-1][1]]
                ax.step(x, y, where="post", color=ALGORITHM["IPOP-CMA-ES"], linewidth=1.6, alpha=0.9, label="two stages" if seed == 1 else None)
            ax.text(40, 1.2e-2, "two stages: CMA-ES on the rates (its own run, left),\nthen CMA-ES on the spike times from that fit,\nσ 0.05 of the box, 600 evaluations more", color=ALGORITHM["IPOP-CMA-ES"],
                    fontsize=11, ha="left", va="center")
            ax.set_xlim(0, 1820)
            ax.set_ylim(2e-3, 2.5)
        else:
            ax.axhline(FLOOR, color=HAIRLINE, linewidth=0.8)
            ax.text(995, FLOOR / 1.4, "exact: all twelve rates equal on the 0.1 ms grid", color=MUTED, fontsize=10.5, ha="right", va="top")
            ax.set_xlim(0, 1000)
            ax.set_ylim(FLOOR / 3, 3)
        ax.set_yscale("log")
        ax.set_xlabel("evaluations (each six simulated seconds)")
        ax.set_ylabel("best error so far")
        ax.legend(loc="upper right", fontsize=11, handlelength=3)
        panel_title(ax, label, sub)
    source(fig, f"neuroscience::snn::models::Adex on the domain's engine (Dormand–Prince, 0.1 ms grid); heuresis strategies::Cma, DifferentialEvolution, swarms::ParticleSwarm, RandomSearch.  {VERSION}")
    save(fig, "search-40-fitting-neurons")


# ---- search-41: identifiability ---------------------------------------------------------------------------

def identifiability_figure():
    land = rows("landscape.csv")
    fits = rows("fits.csv")
    refine = rows("refine.csv")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "What the rates leave open and the spike times pin down",
          "Error on two planes through the hidden parameters (the other two held at the truth), and the parameters each fit recovered")
    gs = fig.add_gridspec(2, 3, left=0.06, right=0.97, top=0.82, bottom=0.11, wspace=0.32, hspace=0.42, width_ratios=[1, 1, 1.15])
    cmap = TERRAIN
    exact = [r for r in fits if r["objective"] == "features" and float(r["error_features"]) == 0.0]
    two = [r for r in refine if r["stage"] == "spike times"]
    for row, (plane, xkey, xlabel) in enumerate((("a-b", "a_nS", "a (nS)"), ("tau_w-b", "tau_w_ms", "τw (ms)"))):
        r = [x for x in land if x["plane"] == plane]
        xs = sorted(set(float(x[xkey]) for x in r))
        ys = sorted(set(float(x["b_pA"]) for x in r))
        for col, (key, name) in enumerate((("error_features", "rates"), ("error_spikes", "spike times"))):
            ax = fig.add_subplot(gs[row, col])
            Z = np.full((len(ys), len(xs)), np.nan)
            ix = {v: i for i, v in enumerate(xs)}
            iy = {v: i for i, v in enumerate(ys)}
            for x in r:
                Z[iy[float(x["b_pA"])], ix[float(x[xkey])]] = max(float(x[key]), 1e-5)
            L = np.log10(Z)
            im = ax.pcolormesh(xs, ys, L, cmap=cmap, vmin=-4 if key == "error_features" else -2, vmax=0.5 if key == "error_features" else 0.7, shading="nearest", rasterized=True)
            if key == "error_features":
                ax.contour(xs, ys, L, levels=[-2], colors=[PAPER], linewidths=0.8)
            ax.plot(TRUTH[xkey], TRUTH["b_pA"], marker="*", markersize=15, color="#FFFFFF", markeredgecolor=INK, markeredgewidth=0.8, linestyle="none", zorder=5)
            if key == "error_features":
                ax.plot([float(x[xkey]) for x in exact], [float(x["b_pA"]) for x in exact], "o", markersize=6, markerfacecolor="none", markeredgecolor="#FFFFFF",
                        markeredgewidth=1.4, zorder=4)
            else:
                ax.plot([float(x[xkey]) for x in two], [float(x["b_pA"]) for x in two], "o", markersize=6, markerfacecolor="none", markeredgecolor="#FFFFFF",
                        markeredgewidth=1.4, zorder=4)
            if plane == "tau_w-b":
                from matplotlib.ticker import NullFormatter, NullLocator
                ax.set_xscale("log")
                ax.xaxis.set_minor_locator(NullLocator())
                ax.xaxis.set_minor_formatter(NullFormatter())
                ax.set_xticks([20, 50, 100, 200, 500])
                ax.set_xticklabels(["20", "50", "100", "200", "500"])
            ax.set_xlabel(xlabel)
            ax.set_ylabel("b (pA)")
            ax.grid(False)
            cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
            cb.set_label("log₁₀ error", color=MUTED, fontsize=10)
            cb.outline.set_visible(False)
            cb.ax.tick_params(labelsize=9, colors=MUTED)
            other = "τw and Vr at the truth" if plane == "a-b" else "a and Vr at the truth"
            note = "; minimum finer than the grid" if key == "error_spikes" else "; dashed: 0.01"
            panel_title(ax, f"{name}: {('a' if plane == 'a-b' else 'τw')} against b", f"41 × 41, {other}{note}")
    # The recovered parameters: each fit's relative error in each parameter.
    ax = fig.add_subplot(gs[:, 2])
    params = [("a_nS", "a"), ("tau_w_ms", "τw"), ("b_pA", "b"), ("V_r_mV", "Vr")]
    groups = [("rates, exact fits", exact, ALGORITHM["DE"], "o", -0.18),
              ("two stages", two, ALGORITHM["IPOP-CMA-ES"], "D", 0.0),
              ("spike times alone (best: DE, seed 3)", [r for r in fits if r["objective"] == "spikes" and r["algorithm"] == "DE" and r["seed"] == "3"], ALGORITHM["CMA-ES"], "s", 0.18)]
    for gname, members, colour, marker, dx in groups:
        for j, (key, _) in enumerate(params):
            vals = [abs(float(m[key]) - TRUTH[key]) / abs(TRUTH[key]) for m in members]
            jitter = np.linspace(-0.05, 0.05, len(vals)) if len(vals) > 1 else [0]
            ax.plot(j + dx + np.array(jitter), np.maximum(vals, 1e-5), marker, color=colour, markersize=7, markerfacecolor=PAPER if marker != "s" else colour,
                    markeredgewidth=1.4, linestyle="none", label=gname if j == 0 else None)
    ax.set_yscale("log")
    ax.set_ylim(5e-6, 2)
    ax.set_xticks(range(4))
    ax.set_xticklabels([p[1] for p in params], fontsize=14)
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylabel("relative error of the recovered value")
    ax.legend(loc="lower left", fontsize=10.5)
    failed = [m for m in two if float(m["error_spikes"]) > 0.1]
    for m in failed:
        ax.annotate(f"seed {m['seed']}: drifted to another basin", xy=(0.03, abs(float(m["a_nS"]) - 2) / 2), xytext=(0.3, 0.75), fontsize=10.5, color=MUTED,
                    arrowprops=dict(arrowstyle="-", color=HAIRLINE, linewidth=0.8))
    n_exact = len(exact)
    panel_title(ax, "The parameters recovered", f"{n_exact} fits of the rates at error 0; 5 two-stage fits")
    source(fig, f"Stars: the truth (a 2 nS, τw 300 ms, b 60 pA, Vr −58 mV). Rings: the fits' values on the plane's two axes (their other parameters are not those of the plane).  {VERSION}")
    save(fig, "search-41-identifiability")


# ---- search-42: the repertoire ---------------------------------------------------------------------------------

CLASSES = ["tonic", "fast", "adapting", "initial burst", "bursting"]
CLASS_COLOUR = {"tonic": "#4E8F8B", "fast": "#CA8DBE", "adapting": "#8E3B78", "initial burst": "#3B2F3A", "bursting": "#D49A3A"}


def repertoire_figure():
    el = rows("elites.csv")
    cov = rows("coverage.csv")
    rep = json.load(open(RESULTS / "repertoire.json"))
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The repertoire of one model neuron: MAP-Elites over six of its parameters",
          "Each cell of (rate, irregularity) holds the cell found nearest Naud et al.'s tonic neuron that fires so under 500 pA, named by its firing class")
    gs = fig.add_gridspec(2, 2, left=0.06, right=0.97, top=0.82, bottom=0.11, wspace=0.2, hspace=0.45, width_ratios=[1.55, 1])
    ax = fig.add_subplot(gs[:, 0])
    lo, hi, res = (0.3, 0.0), (2.5, 1.5), (44, 30)
    wx, wy = (hi[0] - lo[0]) / res[0], (hi[1] - lo[1]) / res[1]
    for r in el:
        if r["algorithm"] != "MAP-Elites":
            continue
        c = int(r["cell"])
        i, j = c % res[0], c // res[0]
        ax.add_patch(plt.Rectangle((lo[0] + i * wx, lo[1] + j * wy), wx * 0.92, wy * 0.92, color=CLASS_COLOUR[r["class"]], linewidth=0))
    # random search's cells that MAP-Elites did not fill, as hairline outlines
    me_cells = {int(r["cell"]) for r in el if r["algorithm"] == "MAP-Elites"}
    for r in el:
        if r["algorithm"] != "random" or int(r["cell"]) in me_cells:
            continue
        c = int(r["cell"])
        i, j = c % res[0], c // res[0]
        ax.add_patch(plt.Rectangle((lo[0] + i * wx, lo[1] + j * wy), wx * 0.92, wy * 0.92, fill=False, edgecolor=INK, linewidth=0.6))
    ex = {r["class"]: r for r in rows("exemplars.csv")}
    for k, r in ex.items():
        ax.plot(float(r["log_rate"]), float(r["cv"]), marker="o", markersize=11, markerfacecolor="none", markeredgecolor=INK, markeredgewidth=1.6)
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], 1.9)
    ax.axhline(hi[1], color=HAIRLINE, linewidth=0.8)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
    ticks = [2, 5, 10, 20, 50, 100, 200]
    ax.set_xticks([math.log10(t) for t in ticks])
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_xlabel("mean firing rate over 1 s (Hz, logarithmic)")
    ax.set_ylabel("coefficient of variation of the interspike intervals")
    ax.grid(False)
    counts = rep["classes_map_elites"]
    handles = [plt.Rectangle((0, 0), 1, 1, color=CLASS_COLOUR[k], label=f"{k} ({counts.get(k, 0)} cells)") for k in CLASSES]
    handles.append(plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor=INK, linewidth=0.6, label="found by random search only"))
    handles.append(Line2D([], [], marker="o", markersize=10, markerfacecolor="none", markeredgecolor=INK, linestyle="none", label="the exemplars drawn in search-43"))
    ax.legend(handles=handles, loc="upper left", fontsize=10.5, ncol=3, columnspacing=1.2, handlelength=1.4)
    panel_title(ax, "The map, MAP-Elites after 20,000 evaluations",
                f"{rep['map_elites']['filled']} of {rep['cells']} cells filled; below 10 Hz a rate is a whole number of spikes in the second, so columns stand apart")

    ax = fig.add_subplot(gs[0, 1])
    for alg, colour, ls in (("MAP-Elites", ALGORITHM["MAP-Elites"], "-"), ("random", ALGORITHM["random"], (0, (2, 2)))):
        r = [x for x in cov if x["algorithm"] == alg]
        ax.plot([int(x["evaluations"]) for x in r], [100 * float(x["coverage"]) for x in r], color=colour, linewidth=2.0, linestyle=ls)
        last = r[-1]
        ax.text(20500, 100 * float(last["coverage"]), f"{alg} {100 * float(last['coverage']):.1f}%", color=colour, fontsize=11, va="center", ha="left")
    ax.set_xlim(0, 27000)
    ax.set_ylim(0, 85)
    ax.set_xticks([0, 5000, 10000, 15000, 20000])
    ax.set_xlabel("evaluations (one simulated second each)")
    ax.set_ylabel("cells filled (%)")
    panel_title(ax, "Coverage", "the same grid; random search's samples offered to an archive")

    ax = fig.add_subplot(gs[1, 1])
    cm, cr = rep["classes_map_elites"], rep["classes_random"]
    y = np.arange(len(CLASSES))
    ax.barh(y + 0.19, [cm.get(k, 0) for k in CLASSES], height=0.36, color=[CLASS_COLOUR[k] for k in CLASSES])
    ax.barh(y - 0.19, [cr.get(k, 0) for k in CLASSES], height=0.36, color=[CLASS_COLOUR[k] for k in CLASSES], alpha=0.35)
    for i, k in enumerate(CLASSES):
        ax.text(cm.get(k, 0) + 6, i + 0.19, f"{cm.get(k, 0)}", va="center", fontsize=10.5, color=INK)
        ax.text(cr.get(k, 0) + 6, i - 0.19, f"{cr.get(k, 0)}", va="center", fontsize=10.5, color=MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels(CLASSES)
    ax.invert_yaxis()
    ax.set_xlim(0, 600)
    ax.set_xlabel("cells of the class")
    ax.grid(axis="y", visible=False)
    panel_title(ax, "Cells by class", "dark: MAP-Elites; pale: random search")
    source(fig, f"heuresis diversity::MapElites on a 44 × 30 grid (Gaussian mutation σ 0.05, iso+line), quality the distance in the unit box from the reference; AdEx on the domain's engine.  {VERSION}")
    save(fig, "search-42-firing-repertoire")


# ---- search-43: the exemplars' potentials ------------------------------------------------------------------------

def traces_figure():
    tr = rows("traces.csv")
    sp = rows("trace_spikes.csv")
    ex = {r["class"]: r for r in rows("exemplars.csv")}
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Five ways one model neuron fires under the same current",
          "The membrane potential over 600 ms under 500 pA of the MAP-Elites elite nearest the median behaviour of each class; spikes drawn to 0 mV")
    gs = fig.add_gridspec(len(CLASSES), 1, left=0.06, right=0.70, top=0.84, bottom=0.10, hspace=0.25)
    for i, k in enumerate(CLASSES):
        ax = fig.add_subplot(gs[i, 0])
        r = [x for x in tr if x["class"] == k]
        t = np.array([float(x["t_ms"]) for x in r])
        v = np.array([float(x["V_mV"]) for x in r])
        ax.plot(t, v, color=CLASS_COLOUR[k], linewidth=1.0)
        s = [float(x["t_ms"]) for x in sp if x["class"] == k]
        ax.vlines(s, -50, 0, color=CLASS_COLOUR[k], linewidth=0.9)
        ax.set_xlim(0, 600)
        ax.set_ylim(-75, 5)
        ax.set_yticks([-70, -35, 0])
        ax.grid(False)
        ax.spines["left"].set_color(HAIRLINE)
        if i < len(CLASSES) - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("time (ms)")
        e = ex[k]
        rate = 10 ** float(e["log_rate"])
        mid = ax.get_position().y0 + ax.get_position().height * 0.5
        fig.text(0.72, mid + 0.052, k, fontsize=14, color=CLASS_COLOUR[k], va="center", ha="left")
        fig.text(0.72, mid - 0.012,
                 f"{rate:.0f} Hz, CV {float(e['cv']):.2f}\n"
                 f"C {float(e['C_pF']):.0f} pF · gL {float(e['g_L_nS']):.1f} nS · a {float(e['a_nS']):.2f} nS\n"
                 f"τw {float(e['tau_w_ms']):.0f} ms · b {float(e['b_pA']):.0f} pA · Vr {float(e['V_r_mV']):.1f} mV\n"
                 f"{int(e['elites'])} cells of the class",
                 fontsize=10.5, color=INK, va="center", ha="left", linespacing=1.45)
    fig.text(0.022, 0.47, "membrane potential (mV)", rotation=90, fontsize=13, color=INK, va="center", ha="center")
    source(fig, f"AdEx (E_L −70 mV, V_T −50 mV, Δ_T 2 mV, cut at 0 mV) on the domain's engine, recorded every 0.3 ms; parameters as found by heuresis diversity::MapElites.  {VERSION}")
    save(fig, "search-43-firing-classes")


if __name__ == "__main__":
    convergence_figure()
    identifiability_figure()
    repertoire_figure()
    traces_figure()
