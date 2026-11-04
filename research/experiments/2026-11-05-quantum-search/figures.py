"""Figures of the quantum-search case study, drawn from results/ in the house style of the search case studies.

    .venv/Scripts/python.exe research/experiments/2026-11-05-quantum-search/figures.py

Every number is read from results/, written by quantum.cpp. The device is never simulated: its measurements are
drawn from the exact law of amplitude amplification (heuresis::quantum::EmulatedGrover).
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys
from collections import defaultdict

import numpy as np
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F401,F403
from search_style import ALGORITHM, CYCLE, HAIRLINE, INK, MUTED, ORCHID, PANEL, PAPER, figure, save, source, title  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE / "results"
EMULATED = "heuresis::quantum, the device's law emulated exactly (no quantum computer simulated)"


def sp(x):
    """A whole number with spaces between its thousands."""
    return f"{x:,.0f}".replace(",", " ")


def rows(name):
    with open(RESULTS / name, newline="") as f:
        return list(csv.DictReader(f))


def figure_law():
    law = rows("grover_law.csv")
    em = rows("grover_emulated.csv")
    N = int(law[0]["N"])
    ts = sorted({int(r["t"]) for r in law}, reverse=True)
    shades = ["#E7B9A1", "#D9738B", "#B5629E", "#8E3B78", "#5C3A75", "#3B2F3A"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Grover's law, and an emulator that keeps it",
          f"The probability that a measurement after r iterations returns one of t marked points among N = {N}: sin²((2r + 1)θ), sin²θ = t/N")
    a = fig.add_axes([0.06, 0.12, 0.50, 0.70])
    b = fig.add_axes([0.64, 0.12, 0.32, 0.70])
    phi = np.linspace(0, np.pi, 400)
    b.plot(phi, np.sin(phi) ** 2, color=HAIRLINE, lw=6, solid_capstyle="round", zorder=1)
    b.plot(phi, np.sin(phi) ** 2, color=INK, lw=1.0, zorder=2)
    for t, col in zip(ts, shades):
        theta = np.arcsin(np.sqrt(t / N))
        e = [r for r in em if int(r["t"]) == t]
        x = np.array([int(r["r"]) for r in e])
        f = np.array([float(r["frequency"]) for r in e])
        if t <= 64:  # the slower rotations; the fast ones are read on the angle's panel
            rr = np.linspace(0, 60, 1200)
            a.plot(rr, np.sin((2 * rr + 1) * theta) ** 2, color=col, lw=1.0, alpha=0.5)
            pts = [r for r in law if int(r["t"]) == t]
            a.plot([int(r["r"]) for r in pts], [float(r["probability"]) for r in pts], "o", ms=3.2, color=col, alpha=0.9)
            a.scatter(x, f, s=46, facecolors="none", edgecolors=col, linewidths=1.1, zorder=4)
        ang = ((2 * x + 1) * theta) % np.pi
        b.scatter(ang, f, s=30, facecolors="none", edgecolors=col, linewidths=1.0, zorder=3, label=f"t/N = 1/{N // t}")
        rstar = int(np.floor(np.pi / (4 * theta)))
        if t <= 64:
            a.annotate(f"r* = {rstar}", xy=(rstar, np.sin((2 * rstar + 1) * theta) ** 2), xytext=(rstar, 1.06), ha="center", fontsize=10, color=col,
                       arrowprops=dict(arrowstyle="-", color=col, lw=0.6))
    a.set_xlim(-1, 61)
    a.set_ylim(-0.02, 1.12)
    a.set_xlabel("Grover iterations r (oracle calls)")
    a.set_ylabel("probability of a marked point")
    a.set_title("against the iterations, t/N from 1/64 to 1/4096: the law (dots) and the emulator's frequencies (rings)", loc="left", fontsize=12.5, color=MUTED)
    b.set_xlim(0, np.pi)
    b.set_xticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
    b.set_xticklabels(["0", "π/4", "π/2", "3π/4", "π"])
    b.set_ylim(-0.02, 1.05)
    b.set_xlabel("rotation angle (2r + 1)θ, modulo π")
    b.set_title("against the angle: one curve, sin²", loc="left", fontsize=12.5, color=MUTED)
    b.legend(loc="lower center", fontsize=10, ncol=2, handletextpad=0.3, columnspacing=1.0)
    # The emulator's agreement with the law, as a statistic.
    z = np.array([float(r["z"]) for r in em if 1e-9 < float(r["probability"]) < 1 - 1e-9])
    chi2 = float(np.sum(z ** 2))
    p = float(stats.chi2.sf(chi2, len(z)))
    K = int(em[0]["measurements"])
    b.text(0.5, 0.50, f"{len(z)} cells of {sp(K)} measurements\nΣz² = {chi2:.1f} on {len(z)} d.f., p = {p:.2f}\nlargest |z| = {np.abs(z).max():.2f}",
           transform=b.transAxes, fontsize=10.5, color=INK, va="center", ha="center", linespacing=1.5)
    source(fig, f"{EMULATED}: EmulatedGrover on N = {N} points with f(i) = i, the first t marked; frequencies from keyed streams.")
    return save(fig, "search-37-grover-law")


def figure_adaptive():
    sc = rows("gas_scaling.csv")
    bound = rows("gas_bound.csv")
    desc = rows("gas_descent.csv")
    qubo = rows("gas_qubo.csv")
    q = json.loads((RESULTS / "qubo.json").read_text())
    by = defaultdict(list)
    meas = defaultdict(list)
    rnd = defaultdict(list)
    for r in sc:
        if r["method"] == "grover":
            by[int(r["N"])].append(float(r["oracle_calls"]))
            meas[int(r["N"])].append(float(r["measurements"]))
        else:
            rnd[int(r["N"])].append(float(r["measurements"]))
    Ns = np.array(sorted(by))
    mean = np.array([np.mean(by[n]) for n in Ns])
    lo = np.array([np.percentile(by[n], 10) for n in Ns])
    hi = np.array([np.percentile(by[n], 90) for n in Ns])
    mm = np.array([np.mean(meas[n]) for n in Ns])
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Grover adaptive search, counted in oracle calls",
          "Oracle calls until the threshold holds the unique minimum, against the size of the space, beside classical random search")
    a = fig.add_axes([0.06, 0.12, 0.47, 0.70])
    col = ALGORITHM["Grover"]
    a.fill_between(Ns, lo, hi, color=col, alpha=0.15, lw=0)
    a.plot(Ns, mean, "o-", color=col, ms=5, lw=1.8, label="Grover adaptive search: mean oracle calls (band: 10–90 %)")
    a.plot(Ns, mm, "--", color=col, lw=1.2, alpha=0.8, label="its measurements, each checked classically")
    bN = np.array([int(r["N"]) for r in bound])
    a.plot(bN, [float(r["dh_expected_bound"]) for r in bound], color=INK, lw=1.2, label="Dürr and Høyer's bound, (45/4)√N + (7/10) lg² N")
    a.plot(Ns, Ns, color=ALGORITHM["random"], lw=1.6, label="classical random search: expectation N")
    rN = np.array(sorted(rnd))
    a.plot(rN, [np.mean(rnd[n]) for n in rN], "s", ms=6, mfc="none", mec=ALGORITHM["random"], mew=1.2, label="random search, measured (heuresis's RandomSearch)")
    big = Ns >= 4096
    slope, icpt = np.polyfit(np.log(Ns[big]), np.log(mean[big]), 1)
    a.set_xscale("log", base=2)
    a.set_yscale("log")
    a.set_xlabel("points in the space N")
    a.set_ylabel("oracle calls (Grover iterations) or evaluations")
    a.legend(loc="upper left", fontsize=10.5)
    a.text(Ns[-1], mean[-1] * 0.55, f"slope {slope:.2f}\nfor N ≥ 2¹²", fontsize=11, color=col, ha="right", va="top")
    a.set_xlim(Ns[0] / 1.4, Ns[-1] * 1.4)

    b = fig.add_axes([0.62, 0.12, 0.34, 0.70])
    runs = sorted({int(r["run"]) for r in desc})
    for i, run in enumerate(runs):
        d = [r for r in desc if int(r["run"]) == run]
        x = [1] + [max(1, int(r["oracle_calls"])) for r in d]
        y = [2 ** q["bits"]] + [int(r["rank"]) for r in d]
        b.step(x, y, where="post", color=CYCLE[i % 8], lw=1.1, alpha=0.85)
    calls = np.array([float(r["oracle_calls"]) for r in qubo])
    b.axhline(1, color=INK, lw=0.8)
    b.set_xscale("log")
    b.set_yscale("log")
    b.set_xlabel("oracle calls so far")
    b.set_ylabel("rank of the threshold among the 2²⁰ points (1: the minimum)")
    b.set_title(f"the threshold's descent on a QUBO of {q['bits']} bits, {len(runs)} runs", loc="left", fontsize=12.5, color=MUTED)
    N = 2 ** q["bits"]
    b.text(0.97, 0.97, f"To the minimum, over {len(calls)} runs: median {sp(np.median(calls))} oracle calls\n(the bound {sp(11.25 * np.sqrt(N) + 0.7 * np.log2(N) ** 2)}; random search expects {sp(N)} evaluations)",
           transform=b.transAxes, fontsize=10.5, color=INK, ha="right", va="top", bbox=dict(boxstyle="square,pad=0.4", fc=PAPER, ec=HAIRLINE, lw=0.6))
    source(fig, f"{EMULATED}: 200 runs at each N (random search 100), each landscape with a unique minimum. The emulator evaluates all N points once, a cost apart.")
    return save(fig, "search-38-grover-adaptive-search"), slope


def figure_inspired():
    d = json.loads((RESULTS / "inspired_agreement.json").read_text())
    conv = rows("inspired_convergence.csv")
    hit = rows("inspired_hitting.csv")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "A quantum-inspired search is a product law",
          f"Han and Kim's QIEA beside PBIL on a QUBO of {d['bits']} bits: the probability of reading each bit as the optimum's, generation by generation")
    G = d["generations"]
    for k, (name, key, x0) in enumerate([(r"QIEA: $\sin^2\varphi_i$ of the first individual's Q-bits", "qiea", 0.06), (r"PBIL: the product law's $p_i$", "pbil", 0.06)]):
        ax = fig.add_axes([x0, 0.55 - k * 0.30, 0.52, 0.25])
        m = np.array(d[key]).T
        im = ax.imshow(m, aspect="auto", cmap=ORCHID, vmin=0, vmax=1, extent=[0.5, G + 0.5, d["bits"] - 0.5, -0.5], interpolation="nearest")
        ax.set_ylabel("bit")
        ax.set_title(name, loc="left", fontsize=12.5, color=INK)
        ax.grid(False)
        if k == 1:
            ax.set_xlabel("generation (10 evaluations each)")
        else:
            ax.set_xticklabels([])
    cax = fig.add_axes([0.59, 0.25, 0.008, 0.55])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("P(bit i reads as the optimum's)")
    cb.outline.set_edgecolor(HAIRLINE)

    c = fig.add_axes([0.745, 0.47, 0.225, 0.33])
    for alg, col in (("QIEA", ALGORITHM["QIEA"]), ("PBIL", ALGORITHM["PBIL"])):
        rr = [r for r in conv if r["algorithm"] == alg]
        ev = sorted({int(r["evaluations"]) for r in rr})
        table = defaultdict(list)
        for r in rr:
            table[int(r["evaluations"])].append(max(float(r["excess"]), 1e-3))
        med = [np.median(table[e]) for e in ev]
        c.fill_between(ev, [np.percentile(table[e], 25) for e in ev], [np.percentile(table[e], 75) for e in ev], color=col, alpha=0.18, lw=0)
        c.plot(ev, med, color=col, lw=1.8, label=alg)
    c.set_xscale("log")
    c.set_yscale("log")
    c.set_ylim(8e-4, 30)
    c.set_yticks([1e-3, 1e-2, 1e-1, 1, 10])
    c.set_yticklabels(["0", "0.01", "0.1", "1", "10"])
    c.set_xlabel("evaluations")
    c.set_ylabel("excess over the minimum (0 at the floor)")
    c.set_title("the best so far (median, quartiles)", loc="left", fontsize=12.5, color=MUTED)
    c.legend(loc="lower left")
    hq = [float(r["evaluations_to_minimum"]) for r in hit if r["algorithm"] == "QIEA"]
    hp = [float(r["evaluations_to_minimum"]) for r in hit if r["algorithm"] == "PBIL"]
    notes = [
        f"QIEA's {sp(d['observed_bits'])} observed bits, against the",
        "product-Bernoulli family's reference draws at",
        r"$p_i = \sin^2\varphi_i$ on the same keyed streams:",
        f"{d['differing_bits']} differ.",
        "",
        "Unentangled Q-bits measured in the computational",
        "basis are a product Bernoulli law; the rotation",
        "gate only moves its parameters.",
        "",
        f"Minimum found: QIEA {d['qiea_found']} of {d['seeds']} runs (median",
        f"{np.median(hq):.0f} evaluations), PBIL {d['pbil_found']} of {d['seeds']} ({np.median(hp):.0f}).",
    ]
    fig.text(0.70, 0.37, "\n".join(notes), fontsize=11, color=INK, va="top", linespacing=1.45)
    source(fig, f"heuresis::quantum::QuantumInspired (10 individuals, global migration every 20 generations) and distributions' PBIL (λ = 10, rate 0.1); {d['seeds']} seeds; minimum by enumerating 2²⁰ points.")
    return save(fig, "search-39-quantum-inspired")


def main():
    print(figure_law())
    p, slope = figure_adaptive()
    print(p, f"slope {slope:.3f}")
    print(figure_inspired())


if __name__ == "__main__":
    main()
