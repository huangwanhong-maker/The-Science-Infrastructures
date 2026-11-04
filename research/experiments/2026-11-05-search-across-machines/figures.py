"""Figures of the search-across-machines case study, drawn from results/ in the house style of the search studies.

    .venv/Scripts/python.exe research/experiments/2026-11-05-search-across-machines/figures.py

Every number is read from results/, written by machines.cpp. Two theoretical curves are computed here, from closed
forms stated in the README: the expected takeover time of binary tournaments on a ring and on a complete graph (exact
Markov chains), and the expected makespan of a generational search on lognormal workers (the expected maximum of 32
lognormal durations, by quadrature).
"""

from __future__ import annotations

import csv
import json
import math
import pathlib
import sys
from collections import defaultdict

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.ticker import MaxNLocator
from scipy import integrate, stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F401,F403
from search_style import ALGORITHM, CYCLE, HAIRLINE, INK, MUTED, ORCHID, PANEL, PAPER, figure, save, source, title  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE / "results"

KIND = {"ring": "#8E3B78", "torus": "#4E8F8B", "torus (Moore)": "#6C8EBF", "hypercube": "#D49A3A", "complete": "#7FA37A"}
SCHEME = {"generational": "#8E3B78", "steady-state": "#4E8F8B", "islands": "#D49A3A"}
SCHEME_NAME = {"generational": "generational GA (μ = λ = 32)", "steady-state": "steady-state GA (μ = 32)", "islands": "32 asynchronous islands (μ = 8, ring)"}


def rows(name):
    with open(RESULTS / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ---- theory -------------------------------------------------------------------------------------------------

def ring_takeover(k: int, p: float = 5 / 9) -> float:
    """The expected takeover time of binary tournaments (with replacement) in the closed neighbourhood of a ring of k.

    The best spreads as two fronts. A cell beside a front draws two of its three closed neighbours and takes the best
    unless both draws miss the best neighbour: probability 1 − (2/3)² = 5/9. While two or more cells remain, each front
    advances independently with probability p; the last cell, between both fronts, falls with 1 − (1/3)² = 8/9.
    E[g] over the gap g of cells not yet taken solves E[g] = 1 + p² E[g−2] + 2p(1−p) E[g−1] + (1−p)² E[g].
    """
    E = [0.0, 9 / 8]
    for g in range(2, k):
        E.append((1 + p * p * E[g - 2] + 2 * p * (1 - p) * E[g - 1]) / (1 - (1 - p) ** 2))
    return E[k - 1]


def complete_takeover(N: int) -> float:
    """The expected takeover time of binary tournaments on a complete graph of N (a Markov chain on the copies k).

    A cell without the best draws two of the N cells and takes the best with probability 1 − (1 − k/N)², independently
    of the others, so k' = k + Bin(N − k, 1 − (1 − k/N)²); E[T_N] = 0 and E[T_k] = (1 + Σ_{k'>k} P(k, k') E[T_k'])/(1 − P(k, k)).
    """
    E = np.zeros(N + 1)
    for k in range(N - 1, 0, -1):
        q = 1 - (1 - k / N) ** 2
        j = np.arange(0, N - k + 1)
        pmf = stats.binom.pmf(j, N - k, q)
        E[k] = (1 + np.dot(pmf[1:], E[k + 1:N + 1])) / (1 - pmf[0])
    return float(E[1])


def expected_max_lognormal(n: int, sigma: float) -> float:
    """E[max of n independent lognormal durations of mean 1], ∫₀^∞ (1 − F(x)ⁿ) dx."""
    if sigma == 0:
        return 1.0
    mu = -sigma * sigma / 2
    F = lambda x: stats.lognorm.cdf(x, s=sigma, scale=math.exp(mu))  # noqa: E731
    val, _ = integrate.quad(lambda x: 1 - F(x) ** n, 0, np.inf, limit=400)
    return val


# ---- figures ------------------------------------------------------------------------------------------------

def wave_image(values):
    """Cells by generation: the best (0) in deep plum; the others from rose (good) to paper (worst)."""
    v = np.array(values, dtype=float).T  # cells × generations
    worst = 255.0
    lo, hi = np.array(to_rgb("#F6EEE8")), np.array(to_rgb("#DDB3CB"))
    frac = ((worst - v) / (worst - 1)) ** 2
    img = lo + (hi - lo) * frac[..., None]
    img[v == 0] = to_rgb("#4A2346")
    return img


def figure_waves():
    d = json.loads((RESULTS / "takeover_waves.json").read_text(encoding="utf-8"))
    kinds = ["ring", "torus", "hypercube", "complete"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "How the best individual takes over a structured population",
          "Selection alone on graphs of 256 cells: each cell, ordered by its distance from the first best (top), generation by generation")
    for row, sel in enumerate(["best", "tournament"]):
        for col, k in enumerate(kinds):
            ax = fig.add_axes([0.08 + col * 0.23, 0.535 - row * 0.385, 0.20, 0.29])
            w = d[k][sel]
            img = wave_image(w["values"])
            T = len(w["values"]) - 1
            ax.imshow(img, aspect="auto", extent=[-0.5, T + 0.5, 255.5, -0.5], interpolation="nearest")
            dist = np.array(d[k]["distance"])
            ax.step(dist, np.arange(256), where="mid", color="#D49A3A", lw=1.2)
            ax.grid(False)
            ax.set_xlim(-0.5, T + 0.5)
            if row == 0:
                ax.set_title(f"{d[k]['name']}\neccentricity {d[k]['eccentricity']}", loc="left", fontsize=13, color=INK)
            ax.text(0.98, 0.04, f"takeover {w['takeover']}", transform=ax.transAxes, ha="right", fontsize=11, color=INK,
                    bbox=dict(boxstyle="square,pad=0.25", fc=PAPER, ec="none", alpha=0.85))
            if col == 0:
                ax.set_ylabel("cell, by distance")
                fig.text(0.018, 0.535 - row * 0.385 + 0.145, "best of the closed\nneighbourhood" if sel == "best" else "binary tournaments\nin the neighbourhood",
                         rotation=90, ha="center", va="center", fontsize=13, color=INK)
            else:
                ax.set_yticklabels([])
            ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
            ax.set_xlabel("generation")
    fig.text(0.04, 0.052, "Deep plum: cells holding the best; rose to paper: the other values, good to worst. Amber: each cell's distance from the start, where the wave "
             "arrives under the first rule; under tournaments it lags.", fontsize=10.5, color=MUTED)
    source(fig, "heuresis::populations::CellularEA without variation (synchronous updates); distances by breadth-first search (morphe::combinatorics). One run of each.")
    return save(fig, "search-33-takeover-waves")


def figure_takeover_times():
    r = rows("takeover_times.csv")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Takeover time is the eccentricity, and tournaments are slower by a factor of their own",
          "Generations until the best fills the graph, on 25 graphs: deterministic best-of-neighbourhood selection, and binary tournaments over 20 seeds")
    a = fig.add_axes([0.06, 0.12, 0.40, 0.70])
    b = fig.add_axes([0.56, 0.12, 0.40, 0.70])
    graphs = defaultdict(list)
    meta = {}
    for x in r:
        graphs[(x["graph"], x["selection"])].append(int(x["takeover"]))
        meta[x["graph"]] = (x["kind"], int(x["cells"]), int(x["eccentricity"]))
    a.plot([0.8, 60], [0.8, 60], color=INK, lw=0.8)
    a.text(40, 33, "takeover = eccentricity", rotation=38, fontsize=10.5, color=INK, ha="center")
    a.plot([4, 140], [4 * 1.8, 140 * 1.8], color=KIND["ring"], lw=0.8, ls="--")
    a.text(70, 70 * 1.8 * 1.12, "9/5 × eccentricity", rotation=38, fontsize=10.5, color=KIND["ring"], ha="center")
    seen = set()
    exact = 0
    for g, (kind, cells, ecc) in meta.items():
        best = graphs[(g, "best")][0]
        exact += best == ecc
        t = np.array(graphs[(g, "tournament")])
        col = KIND[kind]
        a.plot([ecc], [best], "D", ms=7, mfc=PAPER, mec=col, mew=1.3, zorder=4)
        a.errorbar([ecc], [t.mean()], yerr=[[t.mean() - t.min()], [t.max() - t.mean()]], fmt="o", ms=6, color=col, elinewidth=0.9, capsize=2, zorder=3,
                   label=kind if kind not in seen else None)
        seen.add(kind)
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlim(0.8, 200)
    a.set_ylim(0.8, 400)
    a.set_xlabel("eccentricity of the starting cell")
    a.set_ylabel("takeover time (generations)")
    a.set_title(f"diamonds: best of the neighbourhood, on the line in {exact} of {len(meta)} graphs; dots: tournaments, mean and range", loc="left", fontsize=11.5, color=MUTED)
    a.legend(loc="upper left", title="graph", title_fontsize=11)

    # Against the number of cells, with the exact expectations where the chain is simple.
    for kind in ["ring", "torus", "torus (Moore)", "hypercube", "complete"]:
        gs = sorted([g for g, m in meta.items() if m[0] == kind], key=lambda g: meta[g][1])
        if kind == "torus":
            gs = [g for g in gs if meta[g][1] in (16, 64, 256, 1024) and g.split()[-1] in ("4x4", "8x8", "16x16", "32x32")]
        n = [meta[g][1] for g in gs]
        m = [np.mean(graphs[(g, "tournament")]) for g in gs]
        b.plot(n, m, "o", ms=6, color=KIND[kind], label=kind)
        b.plot(n, [graphs[(g, "best")][0] for g in gs], "D", ms=5, mfc=PAPER, mec=KIND[kind], mew=1.0)
    ks = np.array([16, 32, 64, 128, 256])
    b.plot(ks, [ring_takeover(k) for k in ks], color=KIND["ring"], lw=1.2, label="ring: exact expectation")
    Ns = np.array([16, 32, 64, 128, 256, 512, 1024])
    b.plot(Ns, [complete_takeover(int(N)) for N in Ns], color=KIND["complete"], lw=1.2, label="complete: exact expectation")
    b.set_xscale("log", base=2)
    b.set_yscale("log")
    b.set_xlabel("cells")
    b.set_ylabel("takeover time (generations)")
    b.set_title("against the size (dots: tournaments; diamonds: best); lines: exact expectations", loc="left", fontsize=11.5, color=MUTED)
    b.legend(loc="upper left", fontsize=10.5)
    source(fig, "heuresis::populations::CellularEA, selection alone; eccentricities by morphe::combinatorics; expectations of the ring's two fronts and the complete graph's copies computed in figures.py.")
    return save(fig, "search-34-takeover-times")


def figure_synchronisation():
    r = rows("synchronisation.csv")
    ident = rows("identity.csv")
    summ = rows("identity_summary.csv")
    sigmas = sorted({float(x["sigma"]) for x in r if x["straggler"] == "0"})
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "What waiting costs: one search on 32 workers of uneven speed",
          "A GA on Rastrigin's function in 20 dimensions, 32 000 evaluations, lognormal durations of mean 1 and spread σ; one straggler ten times slower")
    axes = [fig.add_axes([0.06, 0.50, 0.27, 0.29]), fig.add_axes([0.40, 0.50, 0.27, 0.29]), fig.add_axes([0.06, 0.09, 0.27, 0.29]), fig.add_axes([0.40, 0.09, 0.27, 0.29])]
    xs = sigmas + [sigmas[-1] + 0.45]
    labels = [f"{s:g}" for s in sigmas] + ["straggler"]

    def series(scheme, key):
        med, lo, hi = [], [], []
        for s in sigmas:
            v = np.array([float(x[key]) for x in r if x["scheme"] == scheme and x["straggler"] == "0" and abs(float(x["sigma"]) - s) < 1e-9])
            med.append(np.median(v))
            lo.append(np.percentile(v, 25))
            hi.append(np.percentile(v, 75))
        v = np.array([float(x[key]) for x in r if x["scheme"] == scheme and x["straggler"] == "1"])
        return np.array(med), np.array(lo), np.array(hi), (np.median(v), np.percentile(v, 25), np.percentile(v, 75))

    def panel(ax, key, ylabel, heading, log=True):
        for sch in ["generational", "steady-state", "islands"]:
            med, lo, hi, st = series(sch, key)
            ok = np.isfinite(med)
            band = ok & np.isfinite(hi)
            ax.fill_between(np.array(sigmas)[band], lo[band], hi[band], color=SCHEME[sch], alpha=0.15, lw=0)
            ax.plot(np.array(sigmas)[ok], med[ok], "o-", ms=4, color=SCHEME[sch], lw=1.6, label=SCHEME_NAME[sch])
            if np.isfinite(st[0]):
                up = st[2] - st[0] if np.isfinite(st[2]) else 0
                ax.errorbar([xs[-1]], [st[0]], yerr=[[st[0] - st[1]], [up]], fmt="s", ms=6, color=SCHEME[sch], capsize=2, elinewidth=0.9)
        if log:
            ax.set_yscale("log")
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, fontsize=10)
        ax.axvline((sigmas[-1] + xs[-1]) / 2, color=HAIRLINE, lw=0.8)
        ax.set_xlabel("spread σ of the durations")
        ax.set_ylabel(ylabel)
        ax.set_title(heading, loc="left", fontsize=12.5, color=INK)

    ideal = 32000 / 32
    theory = [ideal * expected_max_lognormal(32, s) for s in sigmas]
    axes[0].plot(sigmas, theory, color="#E7C9DC", lw=7, solid_capstyle="round", zorder=0)
    panel(axes[0], "makespan", "simulated wall time", "Time to spend the budget")
    axes[0].axhline(ideal, color=INK, lw=0.6, ls=":")
    axes[0].text(1.05, theory[4] * 1.5, "under it, broad and pale:\n1000 × E[max of 32 durations]", fontsize=9.5, color=SCHEME["generational"], ha="right", va="bottom")
    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_, loc="upper left", bbox_to_anchor=(0.04, 0.885), ncol=3, fontsize=11, handlelength=2.2)
    panel(axes[1], "best", "best value found (median, quartiles)", "Quality at the budget (the generational runs alike at every σ)")
    panel(axes[2], "time_to_10", "simulated wall time", "Time until f ≤ 10")
    panel(axes[3], "time_to_1", "simulated wall time", "Time until f ≤ 1 (generational: 6 of 10 runs; islands: none)")

    # The identity.
    ax = fig.add_axes([0.74, 0.09, 0.23, 0.70])
    gen = defaultdict(list)
    for x in ident:
        gen[(x["scheme"], x["placement"])].append((int(x["generation"]), float(x["best"])))
    places = list(dict.fromkeys(x["placement"] for x in ident))
    styles = ["-", "--", ":", "-.", (0, (5, 2, 1, 2)), (0, (1, 3))]
    for i, p in enumerate(places):
        g = sorted(gen[("generational", p)])
        ax.plot([a + 1 for a, _ in g], [b for _, b in g], color=SCHEME["generational"], lw=2.6 - 0.35 * i, ls=styles[i], alpha=0.9, label=p)
        s = sorted(gen[("steady-state", p)])
        ax.plot([a + 1 for a, _ in s], [b for _, b in s], color=SCHEME["steady-state"], lw=0.9, alpha=0.75)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("generation (32 evaluations)")
    ax.set_ylabel("best of the population")
    ax.set_title("The protocol's identity", loc="left", fontsize=12.5, color=INK)
    gsum = [x for x in summ if x["scheme"] == "generational"]
    ssum = [x for x in summ if x["scheme"] == "steady-state"]
    same = sum(int(x["identical"]) for x in gsum)
    worst_s = max(float(x["largest_difference_from_serial"]) for x in ssum)
    ax.legend(loc="lower left", fontsize=9, title="generational GA, placed:", title_fontsize=9.5)
    ax.text(0.03, 0.50, f"Plum: the generational GA on {len(gsum)} placements,\nidentical bit for bit in {same} of {len(gsum)} over all\n{gsum[0]['generations']} generations. Teal: the steady state\non the same placements, as far apart as {worst_s:.0f}.",
            transform=ax.transAxes, fontsize=10, color=INK, va="top")
    source(fig, "heuresis::search::run_asynchronous and populations::run_islands_asynchronous (discrete-event workers); 10 seeds; the pale curve is computed in figures.py.")
    return save(fig, "search-35-synchronisation")


def figure_islands():
    s = rows("islands_summary.csv")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Islands against a deceptive problem: how often, and how connected",
          "16 islands of a GA (μ = 20) on 8 concatenated traps of order 5, 128 000 evaluations: the share of runs that find the optimum")
    tops = ["isolated", "ring", "torus 4 × 4", "hypercube", "complete"]
    intervals = sorted({int(x["interval"]) for x in s if x["topology"] not in ("panmictic", "isolated")})
    grid = np.full((len(tops), len(intervals)), np.nan)
    ann = {}
    for x in s:
        if x["topology"] == "panmictic":
            pan = x
            continue
        i = tops.index(x["topology"])
        if x["topology"] == "isolated":
            for j in range(len(intervals)):
                grid[i, j] = float(x["rate"])
            ann[(i, None)] = x
        else:
            j = intervals.index(int(x["interval"]))
            grid[i, j] = float(x["rate"])
            ann[(i, j)] = x
    ax = fig.add_axes([0.10, 0.14, 0.46, 0.66])
    cmap = LinearSegmentedColormap.from_list("rate", ["#FBF6F1", "#EBD3E3", "#CA8DBE", "#8E3B78", "#3B1F3A"])
    im = ax.imshow(grid, cmap=cmap, vmin=0, vmax=max(0.6, np.nanmax(grid)), aspect="auto")
    for (i, j), x in ann.items():
        txt = f"{x['successes']}/{x['runs']}"
        v = float(x["rate"])
        colour = PAPER if v > 0.35 else INK
        if j is None:
            ax.text((len(intervals) - 1) / 2, i, f"{txt}   (no migration: the interval does not apply)", ha="center", va="center", fontsize=12, color=colour)
        else:
            ax.text(j, i, txt, ha="center", va="center", fontsize=12.5, color=colour)
    ax.set_xticks(range(len(intervals)))
    ax.set_xticklabels([str(v) for v in intervals])
    ax.set_yticks(range(len(tops)))
    ax.set_yticklabels(tops)
    ax.set_xlabel("generations between migrations (one emigrant along each edge)")
    ax.grid(False)
    ax.set_title(f"one panmictic GA of 320 instead: {pan['successes']}/{pan['runs']}", loc="left", fontsize=12.5, color=MUTED)
    cax = fig.add_axes([0.57, 0.14, 0.008, 0.66])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("success rate")
    cb.outline.set_edgecolor(HAIRLINE)

    b = fig.add_axes([0.68, 0.14, 0.29, 0.66])
    cols = {"ring": "#8E3B78", "torus 4 × 4": "#4E8F8B", "hypercube": "#D49A3A", "complete": "#6C8EBF"}
    for k, t in enumerate(["ring", "torus 4 × 4", "hypercube", "complete"]):
        xx = [x for x in s if x["topology"] == t]
        xx.sort(key=lambda x: int(x["interval"]))
        iv = np.array([int(x["interval"]) for x in xx]) * (1 + 0.04 * (k - 1.5))
        est = np.array([float(x["rate"]) for x in xx])
        lo = np.array([float(x["lower"]) for x in xx])
        hi = np.array([float(x["upper"]) for x in xx])
        b.errorbar(iv, est, yerr=[est - lo, hi - est], fmt="o-", ms=5, color=cols[t], capsize=2, elinewidth=0.8, lw=1.4, label=t)
    iso = [x for x in s if x["topology"] == "isolated"][0]
    b.axhline(float(iso["rate"]), color=MUTED, lw=0.9, ls=":", label="isolated")
    b.axhline(float(pan["rate"]), color=INK, lw=0.9, ls="--", label="panmictic GA of 320")
    b.set_xscale("log")
    b.set_xticks(intervals)
    b.set_xticklabels([str(v) for v in intervals])
    b.set_ylim(-0.02, 1.0)
    b.set_xlabel("generations between migrations")
    b.set_ylabel("success rate, with Clopper–Pearson 95 % intervals")
    b.legend(loc="upper right", fontsize=10)
    source(fig, "heuresis::populations::run_islands (synchronous migration, best out, worst in), genetic::GeneticAlgorithm with two-point crossover and bit flips at 1/n; benchmarks::trap.")
    return save(fig, "search-36-island-topologies")


def main():
    print(figure_waves())
    print(figure_takeover_times())
    if (RESULTS / "identity_summary.csv").exists():
        print(figure_synchronisation())
    if (RESULTS / "islands_summary.csv").exists():
        print(figure_islands())


if __name__ == "__main__":
    main()
