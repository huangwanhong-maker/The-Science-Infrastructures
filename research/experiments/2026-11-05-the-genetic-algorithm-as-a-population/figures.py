"""Figures of the genetic algorithm as a population (README.md §7), drawn from results/ only.

    .venv/Scripts/python.exe research/experiments/2026-11-05-the-genetic-algorithm-as-a-population/figures.py
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

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
VERSION = "heuresis 6b3cc88 · root 4153cb9"


def rows(name):
    with open(RESULTS / name, newline="") as f:
        return list(csv.DictReader(f))


def panel_title(ax, text, sub=None):
    ax.set_title(text, loc="left", fontsize=15, color=INK, pad=24 if sub else 10)
    if sub:
        ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=11, color=MUTED, ha="left", va="bottom")


# ---- search-44: one generation's law, and the fixation of one copy ---------------------------------------------

def fixation_figure():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "A genetic algorithm with roulette selection is a Wright–Fisher population",
          "One generation's law of the carriers' count, and the probability that one copy of a variant fixes, against the exact chain and Kimura's diffusion")
    gs = fig.add_gridspec(1, 2, left=0.06, right=0.97, top=0.80, bottom=0.12, wspace=0.22, width_ratios=[1, 1.35])

    # (a) one generation from k = 10 carriers of fitness 1.5 among M = 40.
    ax = fig.add_subplot(gs[0, 0])
    law = rows("one_generation.csv")
    styles = {"roulette": (ALGORITHM["CMA-ES"], "o"), "tournament 2": (ALGORITHM["DE"], "s"), "universal sampling": (ALGORITHM["PSO"], "D")}
    width = 0.28
    for i, (scheme, (colour, marker)) in enumerate(styles.items()):
        r = [x for x in law if x["scheme"] == scheme]
        j = np.array([int(x["count"]) for x in r])
        h = np.array([float(x["frequency"]) for x in r])
        t = np.array([float(x["theory"]) for x in r])
        keep = (h > 0) | (t > 1e-4)
        off = (i - 1) * width
        ax.bar(j[keep] + off, np.minimum(h[keep], 0.195), width=width * 0.9, color=colour, alpha=0.55, linewidth=0)
        inside = t[keep] < 0.195
        ax.plot((j[keep] + off)[inside], t[keep][inside], marker, color=colour, markersize=6, markerfacecolor=PAPER, markeredgewidth=1.4, linestyle="none")
        if scheme == "universal sampling":
            for jj, hh, tt in zip(j[keep], h[keep], t[keep]):
                left = jj == min(j[keep])
                ax.text(jj + off + (-0.35 if left else 0.35), 0.172, f"{hh:.3f}\n({tt:.3f})", color="#9A6A1F", fontsize=10, ha="right" if left else "left", va="center")
    summ = {x["scheme"]: x for x in rows("one_generation_summary.csv")}
    ax.text(2.8, 0.205, "roulette: Binomial(40, 1/3)\ntotal variation " + f"{float(summ['roulette']['total_variation']):.4f}", color=ALGORITHM["CMA-ES"], fontsize=11, ha="left", va="top")
    ax.text(29.5, 0.205, "tournament of two: Binomial(40, 7/16)\ntotal variation " + f"{float(summ['tournament 2']['total_variation']):.4f}", color=ALGORITHM["DE"], fontsize=11, ha="right", va="top")
    ax.text(29.5, 0.152, "universal sampling: 13 or 14 copies only\n(bars cut at 0.195: measured, exact in brackets)\ntotal variation " + f"{float(summ['universal sampling']['total_variation']):.4f}",
            color="#9A6A1F", fontsize=11, ha="right", va="top")
    ax.set_xlim(2, 30)
    ax.set_ylim(0, 0.235)
    ax.set_xlabel("carriers among the 40 children")
    ax.set_ylabel("probability")
    panel_title(ax, "One generation", "k = 10 carriers of fitness 1.5 among M = 40; 200,000 generations; bars measured, rings exact")

    # (b) the fixation probability of one copy, scaled: M·u against S = 2Ms.
    ax = fig.add_subplot(gs[0, 1])
    fx = rows("fixation.csv")
    colours = {20: ALGORITHM["CMA-ES"], 50: ALGORITHM["DE"], 100: ALGORITHM["PSO"], 200: ALGORITHM["CE"]}
    markers = {20: "o", 50: "s", 100: "D", 200: "^"}
    S_dense = np.linspace(-2.2, 16.5, 300)
    limit = np.where(np.abs(S_dense) < 1e-9, 1.0, S_dense / -np.expm1(-S_dense))
    ax.plot(S_dense, limit, color=INK, linewidth=1.2, linestyle=(0, (5, 4)), zorder=1)
    ax.text(12.6, 13.6, "diffusion limit as $s \\to 0$:\n$Mu = S\\,/\\,(1 - e^{-S})$", color=INK, fontsize=11, ha="right", va="bottom")
    for M in (20, 50, 100, 200):
        r = sorted([x for x in fx if int(x["M"]) == M], key=lambda x: float(x["S"]))
        S = np.array([float(x["S"]) for x in r])
        ex = np.array([float(x["u_exact"]) for x in r]) * M
        ga = np.array([float(x["u_ga"]) for x in r]) * M
        se = np.array([float(x["se_ga"]) for x in r]) * M
        dx = {20: -0.24, 50: -0.08, 100: 0.08, 200: 0.24}[M]
        ax.plot(S, ex, color=colours[M], linewidth=1.4, alpha=0.9, zorder=2)
        ax.errorbar(S + dx, ga, yerr=2 * se, fmt=markers[M], color=colours[M], markersize=7, markerfacecolor=PAPER, markeredgewidth=1.5,
                    elinewidth=1.0, capsize=0, zorder=3, label=f"M = {M}")
    ax.annotate("M = 20, s = 0.4:\nGA 0.506, exact 0.507,\nKimura 0.551", xy=(16, 0.506 * 20), xytext=(13.2, 4.2), fontsize=11, color=MUTED, ha="left",
                arrowprops=dict(arrowstyle="-", color=HAIRLINE, linewidth=0.8))
    ax.set_xlabel("S = 2Ms, the scaled selection coefficient")
    ax.set_ylabel("M · u, fixation relative to neutrality")
    ax.set_xlim(-2.8, 17)
    ax.set_ylim(0, 17)
    leg = ax.legend(loc="upper left", title="GA, ±2 s.e.; line: exact chain", fontsize=11, title_fontsize=11)
    leg._legend_box.align = "left"
    panel_title(ax, "Fixation of one copy", "roulette GA, 20,000 runs a point (40,000 for S < 0), no mutation or crossover")
    source(fig, f"heuresis selection::roulette, universal_sampling, tournament; the exact chain by calx's dense solve; Kimura's u from evolution::population.  {VERSION}")
    save(fig, "search-44-ga-as-wright-fisher")


# ---- search-45: drift and takeover ------------------------------------------------------------------------------

def drift_takeover_figure():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "What each selection does to a population: drift without fitness, takeover with it",
          "Heterozygosity of a neutral label under five schemes against 1 − c per generation, and the time for one best individual to fill 500")
    gs = fig.add_gridspec(1, 2, left=0.06, right=0.97, top=0.80, bottom=0.12, wspace=0.2)

    ax = fig.add_subplot(gs[0, 0])
    dr = rows("drift.csv")
    order = [("universal sampling", ALGORITHM["PSO"], "D", "universal sampling, c = 0"),
             ("roulette", ALGORITHM["CMA-ES"], "o", "roulette, c = 1/M"),
             ("tournament 2 fresh ties", ALGORITHM["CE"], "v", "tournament 2, fresh ties, c = 1/M"),
             ("tournament 2", ALGORITHM["DE"], "s", "tournament 2, ties to index, c = (4M²−1)/3M³"),
             ("truncation 1/2", ALGORITHM["ILS"], "^", "truncation to M/2, c = 2/M")]
    for scheme, colour, marker, label in order:
        r = [x for x in dr if x["scheme"] == scheme]
        g = np.array([int(x["generation"]) for x in r])
        h = np.array([float(x["mean_H"]) for x in r])
        t = np.array([float(x["theory_H"]) for x in r])
        ax.plot(g, t, color=colour, linewidth=1.2, alpha=0.8)
        sel = g % 15 == 0
        ax.plot(g[sel], h[sel], marker, color=colour, markersize=6, markerfacecolor=PAPER, markeredgewidth=1.4, linestyle="none", label=label)
    dom = [x for x in dr if x["scheme"] == "domain Wright-Fisher"]
    g = np.array([int(x["generation"]) for x in dom])
    h = np.array([float(x["mean_H"]) for x in dom])
    sel = (g % 15 == 7)
    ax.plot(g[sel], h[sel], "x", color=MUTED, markersize=6, label="evolution domain's Wright–Fisher, 2N = 50")
    ax.set_yscale("log")
    ax.set_ylim(2e-4, 0.8)
    ax.set_xlim(0, 300)
    ax.set_xlabel("generation")
    ax.set_ylabel("mean heterozygosity 2p(1 − p)")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.93), fontsize=10.5)
    panel_title(ax, "Drift of a neutral label, M = 50", "4,000 runs a scheme; markers measured, lines (1 − c)^t / 2 with c = Σ q_j²")

    ax = fig.add_subplot(gs[0, 1])
    tk = rows("takeover.csv")
    times = {x["case"]: x for x in rows("takeover_times.csv")}
    cases = [("truncation 1/4", ALGORITHM["ILS"], "^"), ("tournament 4", ALGORITHM["CE"], "v"), ("roulette s = 4", ALGORITHM["PSO"], "D"),
             ("truncation 1/2", "#B5629E", "<"), ("tournament 2", ALGORITHM["DE"], "s"), ("roulette s = 1", ALGORITHM["CMA-ES"], "o")]
    for case, colour, marker in cases:
        r = [x for x in tk if x["case"] == case]
        g = np.array([int(x["generation"]) for x in r])
        m = np.array([float(x["mean_fraction_fixed"]) for x in r])
        d = np.array([float(x["deterministic"]) for x in r])
        keep = g <= 30
        ax.plot(g[keep], d[keep], color=colour, linewidth=1.2, alpha=0.8, linestyle=(0, (4, 3)))
        ax.plot(g[keep], m[keep], marker=marker, color=colour, markersize=6, markerfacecolor=PAPER, markeredgewidth=1.4, linewidth=1.6)
    # labels at the right, in the order the curves reach the top
    label_y = {"truncation 1/4": 0.40, "tournament 4": 0.33, "roulette s = 4": 0.26, "truncation 1/2": 0.19, "tournament 2": 0.12, "roulette s = 1": 0.05}
    for case, colour, marker in cases:
        t = times[case]
        ax.text(30.5, label_y[case], f"{case}   {float(t['mean_time']):.1f} ± {float(t['sd_time']):.1f}   ({int(t['fixed'])}/1000)",
                color=colour if colour not in (ALGORITHM["PSO"],) else "#9A6A1F", fontsize=11, ha="right", va="center")
    ax.text(30.5, 0.48, "takeover, generations (fixed of 1000 runs)", color=MUTED, fontsize=11, ha="right", va="center")
    ax.set_xlim(0, 31)
    ax.set_ylim(0, 1.03)
    ax.set_xlabel("generation")
    ax.set_ylabel("share of the best, among runs in which it fixed")
    inv = json.load(open(RESULTS / "takeover_invariance.json"))
    same = all(inv.values())
    panel_title(ax, "Takeover from one copy, M = 500",
                "markers measured, dashed deterministic; tournaments and truncation " + ("identical at s = 1 and 4" if same else "differ between s = 1 and 4"))
    source(fig, f"heuresis selection::tournament, truncation, roulette, universal_sampling; evolution::population::wright_fisher_replicates.  {VERSION}")
    save(fig, "search-45-drift-and-takeover")


# ---- search-46: mutation–selection balance on OneMax ------------------------------------------------------------

def balance_figure():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Mutation against selection on OneMax: where a genetic algorithm settles",
          "The equilibrium share of ones in 100 bits under standard bit mutation at rate χ/n, against the infinite population's fixed point")
    gs = fig.add_gridspec(1, 2, left=0.06, right=0.97, top=0.80, bottom=0.12, wspace=0.2, width_ratios=[1.25, 1])

    ax = fig.add_subplot(gs[0, 0])
    bal = rows("balance.csv")
    inf = rows("balance_infinite.csv")
    groups = defaultdict(list)
    for x in bal:
        groups[(x["scheme"], int(x["M"]), float(x["chi"]))].append(float(x["equilibrium"]))
    for scheme, colour in (("roulette", ALGORITHM["CMA-ES"]), ("tournament 2", ALGORITHM["DE"])):
        r = [x for x in inf if x["scheme"] == scheme]
        chi = np.array([float(x["chi"]) for x in r])
        f = np.array([float(x["infinite"]) for x in r])
        ax.plot(chi, f, color=colour, linewidth=2.0)
        if scheme == "roulette":
            le = np.array([float(x["linkage_equilibrium"]) for x in r])
            ax.plot(chi, le, color=colour, linewidth=1.0, linestyle=(0, (2, 2)))
        for M, marker, dx in ((20, "o", 0.93), (100, "s", 1.0), (500, "D", 1.07)):
            ks = sorted(k for k in groups if k[0] == scheme and k[1] == M)
            c = np.array([k[2] for k in ks]) * dx
            m = np.array([np.mean(groups[k]) for k in ks])
            lo = np.array([min(groups[k]) for k in ks])
            hi = np.array([max(groups[k]) for k in ks])
            ax.errorbar(c, m, yerr=[m - lo, hi - m], fmt=marker, color=colour, markersize=6 if M != 500 else 6.5, markerfacecolor=PAPER if M != 500 else colour,
                        markeredgewidth=1.3, elinewidth=0.9, capsize=0, alpha=0.95)
    ax.set_xscale("log")
    ax.set_xticks([0.25, 0.5, 1, 2, 4, 8])
    ax.set_xticklabels(["¼", "½", "1", "2", "4", "8"])
    ax.set_xlim(0.08, 10)
    ax.set_ylim(0.48, 1.01)
    ax.set_xlabel("χ, the expected number of bits flipped in a child")
    ax.set_ylabel("share of ones at equilibrium")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=MUTED, linewidth=2.0, label="infinite population (exact on the classes)"),
               Line2D([], [], color=MUTED, linewidth=1.0, linestyle=(0, (2, 2)), label="roulette, independent loci: (1 + χ)/(1 + 2χ)"),
               Line2D([], [], color=MUTED, marker="o", markerfacecolor=PAPER, linestyle="none", label="M = 20"),
               Line2D([], [], color=MUTED, marker="s", markerfacecolor=PAPER, linestyle="none", label="M = 100"),
               Line2D([], [], color=MUTED, marker="D", linestyle="none", label="M = 500")]
    ax.legend(handles=handles, loc="upper right", fontsize=10.5)
    ax.text(0.11, 0.985, "tournament of two", color=ALGORITHM["DE"], fontsize=12, va="top")
    ax.text(0.095, 0.865, "roulette", color=ALGORITHM["CMA-ES"], fontsize=12, va="top")
    panel_title(ax, "Equilibria", "mean over generations 1000–2000; 4 seeds (range shown); filled: M = 500")

    ax = fig.add_subplot(gs[0, 1])
    tr = rows("balance_trajectories.csv")
    for scheme, colour, ls in (("roulette", ALGORITHM["CMA-ES"], "-"), ("tournament 2", ALGORITHM["DE"], "-")):
        for chi, alpha in ((0.25, 1.0), (1.0, 0.75), (4.0, 0.5)):
            r = [x for x in tr if x["scheme"] == scheme and abs(float(x["chi"]) - chi) < 1e-9]
            g = np.array([int(x["generation"]) for x in r])
            m = np.array([float(x["mean_fraction"]) for x in r])
            ax.plot(g, m, color=colour, linewidth=1.3, alpha=alpha, linestyle=ls)
            ax.text(2010, m[-1], f"χ = {chi:g}", color=colour, fontsize=10.5, va="center", ha="left", alpha=max(alpha, 0.8))
    ax.text(600, 0.975, "tournament of two (the library's GA)", color=ALGORITHM["DE"], fontsize=11)
    ax.text(600, 0.715, "roulette", color=ALGORITHM["CMA-ES"], fontsize=11)
    ax.set_xlim(0, 2300)
    ax.set_ylim(0.48, 1.01)
    ax.set_xlabel("generation")
    ax.set_ylabel("population's mean share of ones")
    panel_title(ax, "Trajectories, M = 100, seed 1", "from random strings; no crossover, no elites")
    source(fig, f"heuresis variation::standard_bit_mutation, selection::roulette, genetic::GeneticAlgorithm; the classes' quasispecies iterated to its fixed point.  {VERSION}")
    save(fig, "search-46-mutation-selection-balance")


if __name__ == "__main__":
    fixation_figure()
    drift_takeover_figure()
    balance_figure()
