"""Figures of the study 'Rugged landscapes', drawn from results/ only.

    .venv/Scripts/python.exe research/experiments/2026-11-05-rugged-landscapes/figures.py

Writes assets/gallery/search-27 ... search-29.
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys
from collections import defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F403

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
N = 20
ALGORITHMS = ["random", "(1+1) EA", "GA", "SA", "tabu", "cGA", "UMDA"]
COLOUR = {
    "random": ALGORITHM["random"],
    "(1+1) EA": ALGORITHM["(1+1) EA"],
    "GA": ALGORITHM["GA"],
    "SA": ALGORITHM["SA"],
    "tabu": "#8E3B78",
    "cGA": ALGORITHM["PBIL"],
    "UMDA": ALGORITHM["PSO"],
}
LABEL = {"random": "random search", "tabu": "tabu search", "SA": "simulated annealing", "cGA": "compact GA"}


def rows(name: str) -> list[dict]:
    with open(RESULTS / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def provenance() -> str:
    v = json.loads((RESULTS / "versions.json").read_text(encoding="utf-8"))
    return (f"heuresis {v['heuresis']} (foundations/heuresis), root {v['root']}; NK landscapes of heuresis::benchmarks, n = 20, "
            "five instances per k, enumerated whole; research/experiments/2026-11-05-rugged-landscapes")


def landscapes():
    out = defaultdict(list)
    for r in rows("landscapes.csv"):
        out[int(r["k"])].append({k: float(v) for k, v in r.items()})
    return dict(sorted(out.items()))


# ------------------------------------------------------------------------------------------------------------------
# search-27: ruggedness against k.

def finite(k):
    """Weinberger's ρ(1) with each table's 2ᵐ entries a finite sample (m = k + 1): the mean squared change of a
    table entry under a flip of one of its inputs is 2σ², while the tables' empirical variance has expectation
    σ²(1 − 2⁻ᵐ), so the deficit (k + 1)/n grows by 2ᵐ/(2ᵐ − 1). At k = 0 it gives 1 − 2/n, the elementary value."""
    k = np.asarray(k, dtype=float)
    m = 2.0 ** (k + 1)
    return 1 - (k + 1) / N * m / (m - 1)

def figure_ruggedness():
    L = landscapes()
    ks = np.array(list(L))
    fig, axes = figure(2, 3)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.83, bottom=0.10, wspace=0.28, hspace=0.38)
    title(fig, "From smooth to rugged: NK landscapes of twenty bits, counted whole",
          "Each dot is one instance, all 2²⁰ points counted, and the line their mean; k is the number of other bits on which each bit's contribution depends")

    def dots(ax, key, transform=lambda v: v):
        for k, ls_ in L.items():
            v = [transform(l[key]) for l in ls_]
            ax.scatter([k] * len(v), v, s=14, color=MUTED, alpha=0.6, linewidths=0, zorder=2)
        mean = [np.mean([transform(l[key]) for l in L[k]]) for k in ks]
        ax.plot(ks, mean, color=INK, lw=1.4, zorder=3)
        return mean

    ax = axes[0, 0]
    dots(ax, "local_optima")
    ax.set_yscale("log")
    ax.axhline(2 ** N / (N + 1), color="#C2665A", lw=1, ls="--")
    ax.text(0.4, 2 ** N / (N + 1) * 0.6, "uncorrelated: $2^n/(n+1)$ = 49 932", color="#C2665A", fontsize=10, ha="left", va="top")
    ax.set_title("Local optima under one-bit flips", loc="left")
    ax.set_xlabel("k")

    ax = axes[0, 1]
    for k, ls_ in L.items():
        ax.scatter([k] * len(ls_), [l["rho1_exact"] for l in ls_], s=14, color=MUTED, alpha=0.6, linewidths=0, zorder=2)
        ax.scatter([k + 0.25] * len(ls_), [l["rho1_walk"] for l in ls_], s=14, marker="x", color=ALGORITHM["CMA-ES"], alpha=0.7, linewidths=0.9, zorder=2)
    kk = np.linspace(0, 19, 100)
    ax.plot(kk, 1 - (kk + 1) / N, color="#C2665A", lw=1, ls="--", label="Weinberger: 1 − (k+1)/n")
    ax.plot(ks, finite(ks), color="#C2665A", lw=1, ls=":", label="with finite tables: deficit × $2^m/(2^m-1)$, $m = k+1$")
    ax.plot(ks, [np.mean([l["rho1_exact"] for l in L[k]]) for k in ks], color=INK, lw=1.4, label="exact, from the whole landscape")
    ax.scatter([], [], marker="x", color=ALGORITHM["CMA-ES"], label="a random walk of 2 × 10⁵ steps")
    ax.set_title("One-step correlation ρ(1)", loc="left")
    ax.set_xlabel("k")
    ax.legend(loc="lower left", fontsize=9.5)

    ax = axes[0, 2]
    dots(ax, "length_exact")
    kl = kk[kk < 18.5]
    ax.plot(kl, -1 / np.log(1 - (kl + 1) / N), color="#C2665A", lw=1, ls="--")
    kf = ks[ks < 19]
    ax.plot(kf, -1 / np.log(finite(kf)), color="#C2665A", lw=1, ls=":")
    ax.text(4.2, 12, "−1/ln ρ(1): dashed Weinberger's,\ndotted with finite tables", color="#C2665A", fontsize=10)
    ax.set_title("Correlation length (steps)", loc="left")
    ax.set_xlabel("k")
    ax.set_ylim(0, 21)

    ax = axes[1, 0]
    dots(ax, "fdc")
    ax.axhline(0, color=HAIRLINE, lw=0.8)
    ax.set_title("Fitness–distance correlation", loc="left")
    ax.set_xlabel("k")

    ax = axes[1, 1]
    dots(ax, "global_basin")
    ax.set_yscale("log")
    ax.set_title("Share of the space in the global optimum's basin", loc="left")
    ax.set_xlabel("k")

    ax = axes[1, 2]
    A = defaultdict(list)
    for r in rows("autocorrelation.csv"):
        A[int(r["k"])].append((int(r["lag"]), float(r["rho_mean"])))
    shown = [0, 2, 5, 10, 19]
    shades = [ORCHID(0.35 + 0.6 * i / (len(shown) - 1)) for i in range(len(shown))]
    for k, c in zip(shown, shades):
        lag, rho = np.array(A[k]).T
        ax.plot(lag, rho, color=c, lw=1.6, label=f"k = {k}")
        ax.plot(lag, finite(k) ** lag, color=c, lw=0.8, ls=":")
    ax.set_title("Autocorrelation along the walk", loc="left")
    ax.set_xlabel("lag s (steps)")
    ax.set_ylabel(r"ρ(s); dotted: $\rho(1)^s$ with finite tables", fontsize=11)
    ax.legend(fontsize=9.5, loc="upper right")
    ax.set_ylim(-0.1, 1.02)
    source(fig, provenance())
    return save(fig, "search-27-ruggedness")


# ------------------------------------------------------------------------------------------------------------------
# search-28: performance against k.

def position(k):
    """The k axis of the performance figure, broken between 10 and 19."""
    return np.where(np.asarray(k) == 19, 12.0, np.asarray(k, dtype=float))


def broken_axis(ax):
    ax.set_xticks([0, 2, 4, 6, 8, 10, 12])
    ax.set_xticklabels(["0", "2", "4", "6", "8", "10", "19"])
    ax.set_xlim(-0.6, 12.6)
    ax.axvspan(10.7, 11.3, color=PAPER, zorder=6)
    ax.text(11.0, 0.0, "//", transform=ax.get_xaxis_transform(), ha="center", va="center", fontsize=11, color=MUTED, zorder=7)
    ax.set_xlabel("k")


def figure_performance():
    P = defaultdict(dict)
    for r in rows("performance.csv"):
        P[(r["algorithm"], float(r["budget"]))][int(r["k"])] = r
    fig, axes = figure(2, 2)
    fig.subplots_adjust(left=0.07, right=0.80, top=0.84, bottom=0.115, wspace=0.2, hspace=0.38)
    title(fig, "Seven searches against ruggedness, at fixed budgets",
          "Over five instances and ten seeds per k: the gap F* − F of the best point found, its rank among the 2²⁰ points, and the share of runs that found the optimum")
    offsets = np.linspace(-0.3, 0.3, len(ALGORITHMS))
    for j, alg in enumerate(ALGORITHMS):
        d = P[(alg, 5000.0)]
        ks = np.array(sorted(d))
        x = position(ks)
        q1 = np.array([float(d[k]["gap_q1"]) for k in ks])
        md = np.array([float(d[k]["gap_median"]) for k in ks])
        q3 = np.array([float(d[k]["gap_q3"]) for k in ks])
        ax = axes[0, 0]
        ax.vlines(x + offsets[j], q1, q3, color=COLOUR[alg], lw=1.1, alpha=0.8)
        ax.scatter(x + offsets[j], md, s=16, color=COLOUR[alg], zorder=3, linewidths=0, label=LABEL.get(alg, alg))
        smooth = ks <= 10
        ax = axes[0, 1]
        rank = np.array([float(d[k]["rank_median"]) + 1 for k in ks])
        ax.plot(x[smooth], rank[smooth], color=COLOUR[alg], lw=1.6, marker="o", ms=3.5)
        ax.plot(x[~smooth], rank[~smooth], color=COLOUR[alg], lw=0, marker="o", ms=4.5)
        for ax, b in [(axes[1, 0], 5000.0), (axes[1, 1], 50000.0)]:
            e = P[(alg, b)]
            y = np.array([float(e[k]["success"]) for k in ks])
            ax.plot(x[smooth], y[smooth], color=COLOUR[alg], lw=1.6, marker="o", ms=3.5)
            ax.plot(x[~smooth], y[~smooth], color=COLOUR[alg], lw=0, marker="o", ms=4.5)
    axes[0, 0].set_title("Gap to the optimum, 5 000 evaluations (median, quartiles)", loc="left", fontsize=14)
    axes[0, 0].set_ylabel("F* − F(best)")
    axes[0, 0].set_ylim(-0.003, 0.11)
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_title("Median rank of the best point found (1 = the optimum)", loc="left", fontsize=14)
    med = 2 ** N * (1 - 0.5 ** (1 / 5000)) + 1
    axes[0, 1].axhline(med, color=MUTED, lw=0.8, ls=":")
    axes[0, 1].text(0.0, med * 1.35, r"random search's median, $2^{20}(1 - 2^{-1/5000}) + 1 \approx 146$", fontsize=9.5, color=MUTED)
    axes[1, 0].set_title("Runs that found the optimum, 5 000 evaluations", loc="left", fontsize=14)
    axes[1, 1].set_title("Runs that found the optimum, 50 000 evaluations", loc="left", fontsize=14)
    for ax in axes[1]:
        ax.set_ylim(-0.03, 1.03)
    for ax in axes.flat:
        broken_axis(ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.815, 0.84), fontsize=11, markerscale=1.8, title="searches", title_fontsize=11)
    source(fig, provenance().split("; research")[0] + "; budgets of 5 000 (0.48% of the space) and 50 000 evaluations")
    return save(fig, "search-28-performance-against-k")


# ------------------------------------------------------------------------------------------------------------------
# search-29: the local optima of three landscapes and their basins.

def figure_optima():
    O = defaultdict(list)
    for r in rows("optima.csv"):
        O[int(r["k"])].append((float(r["fitness"]), float(r["basin"]), float(r["distance"])))
    L = landscapes()
    fig, axes = figure(1, 3)
    fig.subplots_adjust(left=0.06, right=0.92, top=0.80, bottom=0.12, wspace=0.22)
    title(fig, "The local optima of three landscapes, and their basins",
          "Instance 0 at k = 2, 6 and 10: every local optimum under one-bit flips, by its fitness and the size of its basin under steepest descent")
    cmap = LinearSegmentedColormap.from_list("dist", ["#3B1F3A", "#8E3B78", "#CA8DBE", "#E7B9A1", "#D9C3D0"])
    norm = Normalize(0, 12)
    for ax, k in zip(axes, sorted(O)):
        a = np.array(O[k])
        order = np.argsort(-a[:, 2])
        sc = ax.scatter(a[order, 0], a[order, 1], c=a[order, 2], cmap=cmap, norm=norm, s=16 if len(a) < 500 else 7, linewidths=0, alpha=0.9)
        g = a[np.argmin(a[:, 2])]
        ax.scatter([g[0]], [g[1]], marker="*", s=220, facecolor="none", edgecolor=INK, linewidths=1.1, zorder=5)
        ax.set_yscale("log")
        ax.set_xlabel("fitness F of the local optimum")
        inst = [l for l in L[k] if l["instance"] == 0][0]
        ax.set_title(f"k = {k}: {len(a)} local optima", loc="left")
        ax.text(0.03, 0.97, f"global optimum (star): basin {int(g[1])} points\ncorrelation of F with log basin size {inst['fitness_basin_correlation']:.2f}",
                transform=ax.transAxes, fontsize=10, color=MUTED, va="top")
    axes[0].set_ylabel("points in the basin")
    cax = fig.add_axes([0.935, 0.12, 0.012, 0.68])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("Hamming distance to the global optimum")
    cb.outline.set_edgecolor(HAIRLINE)
    source(fig, provenance().split("; research")[0] + "; basins by steepest descent, ties to the lower bit (heuresis::landscapes::basins)")
    return save(fig, "search-29-local-optima-and-basins")


if __name__ == "__main__":
    for f in (figure_ruggedness, figure_performance, figure_optima):
        print(f())
