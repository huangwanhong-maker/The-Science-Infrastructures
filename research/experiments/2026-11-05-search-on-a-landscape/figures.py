"""Figures of the case study 'search on a landscape', drawn from results/ only.

    .venv/Scripts/python.exe research/experiments/2026-11-05-search-on-a-landscape/figures.py

writes assets/gallery/search-01-rosenbrock.jpg, search-02-rastrigin.jpg, search-03-branin.jpg (small multiples of
six algorithms on each landscape), search-04-anytime.jpg (best value against evaluations over 31 seeds) and
search-05-spread.jpg (the distance to the optimum and the spread of each search).
"""

from __future__ import annotations

import csv
import json
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
from matplotlib.patches import Circle, Ellipse  # noqa: E402

# The house serif, with fallbacks for the glyphs it lacks (the superscript minus, the subscript i).
_have = {f.name for f in font_manager.fontManager.ttflist}
matplotlib.rcParams["font.family"] = [f for f in ["Source Serif 4", "Source Serif Pro", "Georgia", "Cambria", "DejaVu Serif"] if f in _have]

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
TRAJ = json.loads((RESULTS / "trajectories.json").read_text())
LAND = json.loads((RESULTS / "landscapes.json").read_text())
VERSION = TRAJ["meta"]["heuresis"]
RECORD = "research/experiments/2026-11-05-search-on-a-landscape"

COLUMNS = ["(1+1)-ES", "CMA-ES", "DE", "PSO", "CE", "BO"]
STAGES = [50, 300, 1500]   # evaluations, for every algorithm but BO
BO_STAGES = [6, 30, 100]   # Bayesian optimisation's: its initial design, and two later surrogates
FILL = LinearSegmentedColormap.from_list("fill", ["#D3ADC1", "#E2C7D3", "#EEDDE1", "#F6ECEA", "#FBF6F1"])

# Per landscape: how its gap f − f* is shaded, where its panels' labels go (away from the start and the minima), ticks.
SPEC = {
    "rosenbrock": dict(title="Six searches in Rosenbrock's valley",
                       formula="f(x) = 100(x₂ − x₁²)² + (1 − x₁)² on [−2, 2] × [−1, 3], minimum 0 at (1, 1)",
                       transform=lambda g: np.log10(np.maximum(g, 1e-4)), fill=np.arange(-2.0, 3.76, 0.5), lines=np.arange(-2.0, 3.76, 0.5),
                       scale="log₁₀ f", corner="right", xticks=[-1, 0, 1], yticks=[0, 1, 2]),
    "rastrigin": dict(title="Six searches on Rastrigin's egg crate",
                      formula="f(x) = 20 + Σ (xᵢ² − 10 cos 2πxᵢ) on [−5.12, 5.12]², minimum 0 at the origin among 120 local minima",
                      transform=lambda g: g, fill=np.arange(0, 81, 5), lines=np.arange(0, 81, 10),
                      scale="f", corner="left", xticks=[-4, -2, 0, 2, 4], yticks=[-4, -2, 0, 2, 4]),
    "branin": dict(title="Six searches on Branin's three minima",
                   formula="Branin's function on [−5, 10] × [0, 15], minimum 0.397887 at (−π, 12.275), (π, 2.275) and (3π, 2.475)",
                   transform=lambda g: np.log10(np.maximum(g, 1e-3)), fill=np.arange(-1.5, 2.6, 0.5), lines=np.arange(-1.5, 2.6, 0.5),
                   scale="log₁₀ (f − f*)", corner="left", xticks=[0, 5], yticks=[5, 10]),
}


def values(seq):
    return np.array([np.nan if v is None else v for v in seq], dtype=float)


def grid_axes(name, n):
    L = LAND["functions"][name]
    return np.linspace(L["lower"][0], L["upper"][0], n), np.linspace(L["lower"][1], L["upper"][1], n)


def sup(e):
    return str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))


def fmt(v):
    if v is None or not math.isfinite(v):
        return "–"
    if v <= 0:
        return "0"
    e = math.floor(math.log10(v))
    if -2 <= e <= 2:
        return f"{v:.3g}"
    return f"{v / 10 ** e:.1f}×10{sup(e)}"


def stage_index(gens, budget):
    k = 0
    for i, g in enumerate(gens):
        if g["e"] <= budget:
            k = i
    return k


def ellipse(ax, centre, cov, k, **kw):
    a = np.array([[cov[0], cov[1]], [cov[1], cov[2]]])
    w, V = np.linalg.eigh(a)
    w = np.maximum(w, 0)
    angle = math.degrees(math.atan2(V[1, 1], V[0, 1]))
    ax.add_patch(Ellipse(centre, 2 * k * math.sqrt(w[1]), 2 * k * math.sqrt(w[0]), angle=angle, fill=False, **kw))


def frame(ax, name, show_x, show_y):
    L, spec = LAND["functions"][name], SPEC[name]
    ax.set_xlim(L["lower"][0], L["upper"][0])
    ax.set_ylim(L["lower"][1], L["upper"][1])
    ax.set_aspect("equal")
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color(HAIRLINE)
        s.set_linewidth(0.8)
    ax.set_xticks(spec["xticks"])
    ax.set_yticks(spec["yticks"])
    ax.tick_params(labelsize=9.5, length=2, colors=MUTED)
    if not show_x:
        ax.set_xticklabels([])
    if not show_y:
        ax.set_yticklabels([])


def backdrop(ax, name, gap=None):
    spec = SPEC[name]
    if gap is None:
        gap = values(LAND["functions"][name]["values"])
    n = int(round(math.sqrt(gap.size)))
    x, y = grid_axes(name, n)
    z = spec["transform"](gap.reshape(n, n))
    ax.contourf(x, y, z, levels=spec["fill"], cmap=FILL, extend="both", zorder=0)
    ax.contour(x, y, z, levels=spec["lines"], colors=MUTED, linewidths=0.35, alpha=0.5, zorder=1)
    for m in LAND["functions"][name]["minimisers"]:
        ax.plot(*m, marker="+", ms=8, mew=1.2, color=INK, zorder=6)


def label(ax, name, top, bottom):
    kw = dict(transform=ax.transAxes, fontsize=9.5, color=INK, zorder=10,
              bbox=dict(boxstyle="square,pad=0.18", fc=PAPER, ec="none", alpha=0.85))
    right = SPEC[name]["corner"] == "right"
    ax.text(0.03, 0.97, top, ha="left", va="top", **kw)
    ax.text(0.97 if right else 0.03, 0.03, bottom, ha="right" if right else "left", va="bottom", **kw)


def draw(ax, name, algorithm, stage):
    run = TRAJ["functions"][name]["algorithms"][algorithm]
    gens = run["generations"]
    colour = ALGORITHM[algorithm]
    start = TRAJ["functions"][name]["start"]
    budget = (BO_STAGES if algorithm == "BO" else STAGES)[stage]
    k = stage_index(gens, budget)
    g = gens[k]
    if algorithm == "BO":
        snaps = {s["count"]: s for s in run["gp"]["snapshots"]}
        backdrop(ax, name, values(snaps[min(snaps, key=lambda c: abs(c - budget))]["mean"]))
    else:
        backdrop(ax, name)
    if algorithm in ("(1+1)-ES", "CMA-ES", "CE"):
        ax.plot(*start, marker="x", ms=6, mew=1.0, color=INK, alpha=0.7, zorder=4)

    if algorithm == "(1+1)-ES":
        path = []
        for h in gens[: k + 1]:
            if not path or h["m"] != path[-1]:
                path.append(h["m"])
        p = np.array(path)
        ax.plot(p[:, 0], p[:, 1], color=colour, lw=1.2, zorder=3)
        recent = np.array([h["pts"][0] for h in gens[max(0, k - 30): k + 1]])
        ax.scatter(recent[:, 0], recent[:, 1], s=10, facecolor="none", edgecolor=colour, lw=0.8, zorder=4)
        ax.add_patch(Circle(g["m"], g["s"], fill=False, ec=colour, lw=1.4, zorder=5))
        ax.scatter(*g["m"], s=28, color=colour, edgecolor=INK, lw=0.6, zorder=6)
    elif algorithm == "CMA-ES":
        means = np.array([h["m"] for h in gens[: k + 1]])
        ax.plot(means[:, 0], means[:, 1], color=colour, lw=1.0, alpha=0.8, zorder=3)
        cov = [g["s"] ** 2 * c for c in g["C"]]
        ellipse(ax, g["m"], cov, 1, ec=colour, lw=1.5, zorder=5)
        ellipse(ax, g["m"], cov, 2, ec=colour, lw=0.7, ls=(0, (3, 2)), zorder=5)
        pts = np.array(g["pts"])
        ax.scatter(pts[:, 0], pts[:, 1], s=15, color=colour, alpha=0.8, lw=0, zorder=5)
        ax.scatter(*g["m"], s=28, color=colour, edgecolor=INK, lw=0.6, zorder=6)
    elif algorithm == "DE":
        pop = np.array(g["pop"])
        pf = values(g["pf"])
        ax.scatter(pop[:, 0], pop[:, 1], s=18, color=colour, alpha=0.9, lw=0, zorder=5)
        ax.scatter(*pop[int(np.nanargmin(pf))], s=64, facecolor="none", edgecolor=INK, lw=0.8, zorder=6)
    elif algorithm == "PSO":
        pts = np.array(g["pts"])
        v = np.array(g["v"])
        tail = pts if k == 0 else pts - v  # the first generation's velocities are those it will move with
        ax.quiver(tail[:, 0], tail[:, 1], v[:, 0], v[:, 1], angles="xy", scale_units="xy", scale=1,
                  color=colour, width=0.006, headwidth=3.5, headlength=4, alpha=0.8, zorder=4)
        ax.scatter(pts[:, 0], pts[:, 1], s=13, color=colour, edgecolor=INK, lw=0.3, zorder=5)
    elif algorithm == "CE":
        pts = np.array(g["pts"])
        elite = np.argsort(values(g["f"]), kind="stable")[:10]
        ax.scatter(pts[:, 0], pts[:, 1], s=9, color=colour, alpha=0.35, lw=0, zorder=4)
        ax.scatter(pts[elite, 0], pts[elite, 1], s=15, color=colour, lw=0, zorder=5)
        ellipse(ax, g["m"], g["C"], 1, ec=colour, lw=1.5, zorder=5)
        ellipse(ax, g["m"], g["C"], 2, ec=colour, lw=0.7, ls=(0, (3, 2)), zorder=5)
        ax.scatter(*g["m"], s=28, color=colour, edgecolor=INK, lw=0.6, zorder=6)
    elif algorithm == "BO":
        seen = np.array([p for h in gens[: k + 1] for p in h["pts"]])
        ax.scatter(seen[:, 0], seen[:, 1], s=13, color=colour, edgecolor=INK, lw=0.3, zorder=5)
        ax.scatter(*seen[-1], s=64, facecolor="none", edgecolor=INK, lw=0.8, zorder=6)

    ended = k == len(gens) - 1 and run["ending"] != "budget" and g["e"] < budget
    top = f"{g['e']:.0f} evaluations · its surrogate" if algorithm == "BO" else f"gen {k + 1} · {g['e']:.0f} evaluations" + (" · ended" if ended else "")
    label(ax, name, top, f"best f − f* = {fmt(g['gap'])}")


def smallmultiples(name, number):
    spec = SPEC[name]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    left, right, top, bottom = 0.085, 0.985, 0.835, 0.075
    wgap, hgap = 0.006, 0.012
    cols, rows = len(COLUMNS), len(STAGES)
    w = (right - left - (cols - 1) * wgap) / cols
    h = (top - bottom - (rows - 1) * hgap) / rows
    for r in range(rows):
        for c, algorithm in enumerate(COLUMNS):
            ax = fig.add_axes([left + c * (w + wgap), top - (r + 1) * h - r * hgap, w, h])
            draw(ax, name, algorithm, r)
            frame(ax, name, show_x=r == rows - 1, show_y=c == 0)
            if r == 0:
                ax.set_title(algorithm, fontsize=14, color=ALGORITHM[algorithm], pad=6, loc="center")
        fig.text(left - 0.045, top - (r + 0.5) * h - r * hgap, f"≈ {STAGES[r]}\nevaluations\n(BO: {BO_STAGES[r]})",
                 ha="center", va="center", fontsize=11, color=MUTED, linespacing=1.3)
    title(fig, spec["title"],
          "Seed 0 of each search; the (1+1)-ES, CMA-ES and CE start at ×, DE, PSO and BO at uniform points of the box. ES: path, recent offspring, 1-σ circle. "
          "CMA-ES, CE: mean,\nsamples, 1- and 2-σ ellipses (CE's elite solid). DE: population, best ringed. PSO (the Standard PSO 2011): particles, velocities. "
          "BO: its samples over its surrogate's mean, the latest ringed.")
    source(fig, f"heuresis {VERSION} · {RECORD} · {spec['formula']} · shaded by {spec['scale']} · from results/trajectories.json")
    return save(fig, f"search-0{number}-{name}")


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def num(s):
    return float(s) if s not in ("", None) else math.nan


def end_labels(ax, items, x, min_gap=0.42, lo=None, hi=None):
    """Labels at the ends of curves on a log y axis, spread apart by at least min_gap decades."""
    items = sorted(items, key=lambda t: t[1])
    ys = [math.log10(max(y, 1e-30)) for _, y, _ in items]
    if lo is not None:
        ys = [max(y, lo) for y in ys]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + min_gap)
    if hi is not None and ys and ys[-1] > hi:
        ys[-1] = hi
        for i in range(len(ys) - 2, -1, -1):
            ys[i] = min(ys[i], ys[i + 1] - min_gap)
    for (text, _, colour), y in zip(items, ys):
        ax.text(x, 10 ** y, text, color=colour, fontsize=11, va="center", ha="left")


def log_axis(ax, right):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, right)
    ax.set_xticks([1, 10, 100, 1000])
    ax.set_xticklabels(["1", "10", "100", "1000"])
    ax.axvline(2000, color=HAIRLINE, lw=0.8)
    ax.text(2000, 1.0, "2000", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=9.5, color=MUTED)
    ax.grid(False, axis="x")


NAMES = {"rosenbrock": "Rosenbrock", "rastrigin": "Rastrigin", "branin": "Branin"}


def anytime():
    data = defaultdict(lambda: defaultdict(list))
    for r in read_csv(RESULTS / "anytime.csv"):
        data[r["function"]][r["algorithm"]].append((num(r["budget"]), num(r["q25"]), num(r["q50"]), num(r["q75"])))
    order = ["random", "(1+1)-ES", "CMA-ES", "DE", "PSO", "CE", "BO"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 10), dpi=100)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.80, bottom=0.11, wspace=0.2)
    floor = 1e-13
    for ax, f in zip(axes, ["rosenbrock", "rastrigin", "branin"]):
        ends = []
        for a in order:
            d = np.array(data[f][a])
            b, q25, q50, q75 = d[:, 0], *(np.clip(d[:, i], floor, None) for i in (1, 2, 3))
            c = ALGORITHM[a]
            ax.fill_between(b, q25, q75, color=c, alpha=0.13, lw=0, step="post")
            ax.step(b, q50, where="post", color=c, lw=1.9 if a != "random" else 1.4)
            ends.append((a, q50[-1], c))
        log_axis(ax, 2000 * 5.5)
        ax.set_ylim(floor / 3, 3e3)
        ax.axhline(1e-12, color=MUTED, lw=0.7, ls=(0, (2, 3)))
        ax.set_title(NAMES[f], loc="left", fontsize=16, color=INK, pad=14)
        ax.set_xlabel("evaluations")
        end_labels(ax, [e for e in ends if e[0] != "BO"], 2000 * 1.15, lo=-12.6, hi=3.2)
        for a, y, c in ends:
            if a == "BO":
                ax.text(100 * 1.15, max(y, floor), "BO", color=c, fontsize=11, va="center", ha="left")
    axes[0].set_ylabel("best f − f*: median, and quartiles as a band")
    title(fig, "The best value found, against evaluations",
          "Median and quartiles over 31 seeds of each search on each two-dimensional landscape; runs end at f − f* ≤ 10⁻¹² (dotted).\n"
          "Bayesian optimisation is given 100 evaluations, the others 2000; random search samples the box uniformly; PSO is the Standard PSO 2011.")
    source(fig, f"heuresis {VERSION} · {RECORD} · quartiles by heuresis::analysis::anytime_quantile (type 7) · from results/anytime.csv")
    return save(fig, "search-04-anytime")


def spread():
    data = defaultdict(lambda: defaultdict(list))
    for r in read_csv(RESULTS / "spread.csv"):
        data[r["function"]][r["algorithm"]].append([num(r[k]) for k in ("budget", "distance_q25", "distance_q50", "distance_q75", "spread_q25", "spread_q50", "spread_q75")])
    order = ["(1+1)-ES", "CMA-ES", "DE", "PSO", "CE", "BO"]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), dpi=100)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.82, bottom=0.115, wspace=0.2, hspace=0.16)
    for j, f in enumerate(["rosenbrock", "rastrigin", "branin"]):
        for i, (lo, mid, hi, what) in enumerate([(1, 2, 3, "distance"), (4, 5, 6, "spread")]):
            ax = axes[i, j]
            ends = []
            for a in order:
                if what == "spread" and a == "BO":
                    continue
                d = np.array(data[f][a])
                ok = np.isfinite(d[:, mid])
                c = ALGORITHM[a]
                ax.fill_between(d[ok, 0], d[ok, lo], d[ok, hi], color=c, alpha=0.12, lw=0, step="post")
                ax.step(d[ok, 0], d[ok, mid], where="post", color=c, lw=1.8)
                ends.append((a, d[ok, mid][-1], c, d[ok, 0][-1]))
            log_axis(ax, 2000 * 6)
            floor = 1e-8 if what == "distance" else 1e-10
            ax.set_ylim(floor, 30)
            end_labels(ax, [(a, y, c) for a, y, c, b in ends if a != "BO"], 2000 * 1.15, lo=math.log10(floor) + 0.2, hi=1.2)
            for a, y, c, b in ends:
                if a == "BO":
                    ax.text(b * 1.15, y, "BO", color=c, fontsize=11, va="center", ha="left")
            if i == 0:
                ax.set_title(NAMES[f], loc="left", fontsize=16, color=INK, pad=14)
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("evaluations")
    axes[0, 0].set_ylabel("distance of the best point\nto the nearest minimiser")
    axes[1, 0].set_ylabel("spread of the search:\nr.m.s. distance of a candidate")
    title(fig, "How far each search is from the optimum, and how widely it searches",
          "Medians, with quartiles as bands, over 31 seeds. Spread: σ√n for the (1+1)-ES, σ√tr C for CMA-ES, √tr C for CE, the r.m.s. distance from the centroid for\n"
          "DE's population and PSO's particles. A search whose spread falls while its distance does not has converged where the optimum is not.")
    source(fig, f"heuresis {VERSION} · {RECORD} · the state after the last generation completed within each budget · from results/spread.csv")
    return save(fig, "search-05-spread")


if __name__ == "__main__":
    for i, name in enumerate(["rosenbrock", "rastrigin", "branin"], start=1):
        print(smallmultiples(name, i))
    print(anytime())
    print(spread())
