"""Figures of the study 'Fronts in three objectives', drawn from results/ only.

    .venv/Scripts/python.exe research/experiments/2026-11-05-fronts-in-three-objectives/figures.py

Writes assets/gallery/search-20 ... search-23.
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
from matplotlib.ticker import FuncFormatter  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"

ALGORITHMS = ["NSGA-II", "NSGA-III", "MOEA/D-Tch", "MOEA/D-PBI", "SMS-EMOA"]
COLOUR = {
    "NSGA-II": ALGORITHM["NSGA-II"],
    "NSGA-III": ALGORITHM["NSGA-III"],
    "MOEA/D-Tch": ALGORITHM["MOEA/D"],
    "MOEA/D-PBI": "#D9738B",
    "SMS-EMOA": ALGORITHM["SMS-EMOA"],
}
LABEL = {"MOEA/D-Tch": "MOEA/D, Tchebycheff", "MOEA/D-PBI": "MOEA/D, PBI"}
PROVENANCE = "heuresis {h} (foundations/heuresis), root {r}; 15 seeds per algorithm and problem; research/experiments/2026-11-05-fronts-in-three-objectives"


def rows(name: str) -> list[dict]:
    with open(RESULTS / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def versions() -> dict:
    return json.loads((RESULTS / "versions.json").read_text(encoding="utf-8"))


def provenance() -> str:
    v = versions()
    return PROVENANCE.format(h=v["heuresis"], r=v["root"])


def fronts() -> dict:
    out = defaultdict(list)
    for r in rows("fronts.csv"):
        out[r["problem"]].append((int(r["component"]), float(r["f1"]), float(r["f2"]), float(r["f3"])))
    return {k: np.array(v) for k, v in out.items()}


def populations() -> dict:
    out = defaultdict(list)
    for r in rows("populations.csv"):
        out[(r["problem"], r["algorithm"])].append((float(r["f1"]), float(r["f2"]), float(r["f3"])))
    return {k: np.array(v) for k, v in out.items()}


def finals() -> dict:
    out = defaultdict(list)
    for r in rows("final.csv"):
        out[(r["problem"], r["algorithm"])].append(r)
    return out


def style_3d(ax, lim, zlim=None, labels=("$f_1$", "$f_2$", "$f_3$")):
    """Three-dimensional axes in the house style: hairline box and grid, no tick labels, the objectives named."""
    ax.set_facecolor(PAPER)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((0.984, 0.965, 0.945, 0.0))
        axis.line.set_color(HAIRLINE)
        axis._axinfo["grid"]["color"] = HAIRLINE
        axis._axinfo["grid"]["linewidth"] = 0.4
        axis._axinfo["tick"]["color"] = (0, 0, 0, 0)
        axis.set_ticklabels([])
    ax.set_xlim(*lim[0])
    ax.set_ylim(*lim[1])
    ax.set_zlim(*(zlim or lim[2]))
    ax.set_xlabel(labels[0], fontsize=11, labelpad=-12, color=MUTED)
    ax.set_ylabel(labels[1], fontsize=11, labelpad=-12, color=MUTED)
    ax.set_zlabel(labels[2], fontsize=11, labelpad=-12, color=MUTED)


def heading(fig, x, y, alg, note):
    """An algorithm's name in its colour, and a line of medians beneath it."""
    fig.text(x, y, LABEL.get(alg, alg), fontsize=13, color=COLOUR[alg], ha="left", va="bottom")
    fig.text(x, y - 0.006, note, fontsize=9.5, color=MUTED, ha="left", va="top")


def median(xs):
    return float(np.median(np.array(xs, dtype=float)))


# ------------------------------------------------------------------------------------------------------------------
# search-20: final populations on DTLZ1 and DTLZ2 in three dimensions, against the true fronts.

def figure_fronts_dtlz():
    F = fronts()
    P = populations()
    fin = finals()
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Five searches on the fronts of DTLZ1 and DTLZ2",
          "The final population of the median run (by hypervolume) after 40 000 evaluations, against the true front drawn as a surface")
    problems = [("DTLZ1", "DTLZ1: the simplex $f_1+f_2+f_3=1/2$, behind $11^5-1$ local fronts", [(0, 0.55), (0, 0.55), (0, 0.55)]),
                ("DTLZ2", "DTLZ2: the positive octant of the unit sphere", [(0, 1.05), (0, 1.05), (0, 1.05)])]
    for row, (prob, text, lim) in enumerate(problems):
        tri = F[prob]
        top = 0.86 - row * 0.425
        fig.text(0.04, top, text, fontsize=14, color=INK, ha="left", va="bottom")
        for col, alg in enumerate(ALGORITHMS):
            runs = fin[(prob, alg)]
            hv = median([float(r["hv_ratio"]) for r in runs])
            on = median([int(r["on_front"]) for r in runs])
            size = median([int(r["front_size"]) for r in runs])
            heading(fig, 0.045 + col * 0.19, top - 0.055, alg, f"HV/HV* {hv:.3f} · {on:.0f} of {size:.0f} on the front")
            ax = fig.add_axes([0.02 + col * 0.19, top - 0.405, 0.2, 0.345], projection="3d")
            style_3d(ax, lim)
            ax.plot_trisurf(tri[:, 1], tri[:, 2], tri[:, 3], color="#E6D3DE", alpha=0.55, linewidth=0, antialiased=True, shade=False)
            pts = P[(prob, alg)]
            inside = np.all(pts <= np.array([l[1] for l in lim]), axis=1)
            q = pts[inside]
            ax.scatter(q[:, 0], q[:, 1], q[:, 2], s=11, color=COLOUR[alg], depthshade=False, linewidths=0, alpha=0.95)
            ax.view_init(elev=24, azim=40)
            if (~inside).sum():
                ax.text2D(0.5, 0.02, f"{(~inside).sum()} point beyond the box" if (~inside).sum() == 1 else f"{(~inside).sum()} points beyond the box",
                          transform=ax.transAxes, fontsize=9, color=MUTED, ha="center")
    source(fig, provenance() + "; the notes are medians over the 15 seeds")
    return save(fig, "search-20-fronts-dtlz1-dtlz2")


# ------------------------------------------------------------------------------------------------------------------
# search-21: the disconnected fronts, DTLZ7 (four patches) and ZDT3 (five segments).

def figure_disconnected():
    F = fronts()
    P = populations()
    fin = finals()
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Disconnected fronts: who covers every piece",
          "Final populations of the median runs on DTLZ7 (four patches, 40 000 evaluations) and ZDT3 (five segments, 25 000 evaluations)")
    d7 = F["DTLZ7"]
    fig.text(0.04, 0.855, r"DTLZ7 in three objectives: $f_3 = 6 - \varphi(f_1) - \varphi(f_2)$, $\varphi(t) = t\,(1 + \sin 3\pi t)$, on four patches", fontsize=14, color=INK, ha="left", va="bottom")
    for col, alg in enumerate(ALGORITHMS):
        runs = fin[("DTLZ7", alg)]
        comps = np.array([[int(x) for x in r["components"].split("|")] for r in runs])
        whole = int(((comps > 0).sum(axis=1) == 4).sum())
        med = np.median(comps, axis=0)
        hv = median([float(r["hv_ratio"]) for r in runs])
        per_patch = " · ".join(f"{m:.0f}" for m in med)
        heading(fig, 0.045 + col * 0.19, 0.80, alg, f"HV/HV* {hv:.3f} · all four patches in {whole} of 15 runs\npoints per patch {per_patch}")
        ax = fig.add_axes([0.02 + col * 0.19, 0.42, 0.2, 0.33], projection="3d")
        style_3d(ax, [(0, 1), (0, 1), (2.5, 6.1)])
        for c in range(4):
            q = d7[d7[:, 0] == c]
            ax.plot_trisurf(q[:, 1], q[:, 2], q[:, 3], color="#E2CCD9", alpha=0.6, linewidth=0, shade=False)
        pts = P[("DTLZ7", alg)]
        inside = pts[:, 2] <= 6.1
        q = pts[inside]
        ax.scatter(q[:, 0], q[:, 1], q[:, 2], s=11, color=COLOUR[alg], depthshade=False, linewidths=0)
        ax.view_init(elev=32, azim=-60)
    z3 = F["ZDT3"]
    fig.text(0.04, 0.37, r"ZDT3 in two objectives: $f_2 = 1 - \sqrt{f_1} - f_1 \sin 10\pi f_1$, on five segments", fontsize=14, color=INK, ha="left", va="bottom")
    for col, alg in enumerate(ALGORITHMS):
        runs = fin[("ZDT3", alg)]
        comps = np.array([[int(x) for x in r["components"].split("|")] for r in runs])
        whole = int(((comps > 0).sum(axis=1) == 5).sum())
        hv = median([float(r["hv_ratio"]) for r in runs])
        heading(fig, 0.045 + col * 0.19, 0.325, alg, f"HV/HV* {hv:.4f} · all five in {whole} of 15 runs")
        ax = fig.add_axes([0.045 + col * 0.19, 0.105, 0.165, 0.175])
        for c in np.unique(z3[:, 0]):
            q = z3[z3[:, 0] == c]
            ax.plot(q[:, 1], q[:, 2], color="#DCC6D3", lw=5, solid_capstyle="round", zorder=1)
        pts = P[("ZDT3", alg)]
        ax.scatter(pts[:, 0], pts[:, 1], s=10, color=COLOUR[alg], zorder=2, linewidths=0)
        ax.set_xlim(-0.03, 0.9)
        ax.set_ylim(-0.85, 1.05)
        ax.tick_params(labelsize=9)
        ax.set_xlabel("$f_1$", fontsize=11, color=MUTED, labelpad=1)
        if col == 0:
            ax.set_ylabel("$f_2$", fontsize=11, color=MUTED)
    source(fig, provenance() + "; a point counts on a piece when it lies within 0.01 of it; medians over 15 seeds")
    return save(fig, "search-21-disconnected-fronts")


# ------------------------------------------------------------------------------------------------------------------
# search-22: hypervolume and IGD against evaluations.

def figure_convergence():
    C = defaultdict(list)
    for r in rows("curves.csv"):
        C[(r["problem"], r["algorithm"])].append(r)
    ref = json.loads((RESULTS / "reference.json").read_text(encoding="utf-8"))
    hv_true = {p["name"]: p["hv_true"] for p in ref["problems"]}
    problems = ["DTLZ1", "DTLZ2", "DTLZ7", "ZDT1", "ZDT3"]
    fig, axes = figure(2, 5)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.83, bottom=0.12, wspace=0.28, hspace=0.32)
    title(fig, "How fast each search approaches the front",
          "Median over 15 seeds (line) and the interquartile range (band); hypervolume as a share of the true front's, IGD to a dense sample of the front")
    for col, prob in enumerate(problems):
        for alg in ALGORITHMS:
            rs = C[(prob, alg)]
            e = np.array([float(r["evaluations"]) for r in rs]) / 1000
            for row, (lo, mid, hi) in enumerate([("hv_ratio_q1", "hv_ratio_median", "hv_ratio_q3"), ("igd_q1", "igd_median", "igd_q3")]):
                ax = axes[row, col]
                y = np.array([float(r[mid]) for r in rs])
                ax.fill_between(e, [float(r[lo]) for r in rs], [float(r[hi]) for r in rs], color=COLOUR[alg], alpha=0.15, linewidth=0)
                ax.plot(e, y, color=COLOUR[alg], lw=1.6, label=LABEL.get(alg, alg))
        axes[0, col].set_title(f"{prob}", loc="left", fontsize=15)
        axes[0, col].set_ylim(0.0 if prob == "DTLZ1" else 0.5, 1.0)
        axes[1, col].set_yscale("log")
        axes[1, col].yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        axes[1, col].yaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
        axes[1, col].set_xlabel("evaluations (thousands)")
        for row in range(2):
            axes[row, col].set_xlim(0, e.max())
        # The same measures of a natural set of about μ points lying on the front (reference.json), where there is one.
        p = next(q for q in ref["problems"] if q["name"] == prob)
        if "hv_set_ratio" in p:
            axes[0, col].axhline(p["hv_set_ratio"], color=INK, lw=0.8, ls=(0, (3, 3)), zorder=0)
            axes[1, col].axhline(p["igd_set"], color=INK, lw=0.8, ls=(0, (3, 3)), zorder=0)
    axes[0, 0].set_ylabel("hypervolume / HV of the front")
    axes[1, 0].set_ylabel("IGD")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles.append(plt.Line2D([], [], color=INK, lw=0.8, ls=(0, (3, 3))))
    labels.append("μ points on the front, evenly set")
    axes[0, 1].legend(handles, labels, loc="lower right", fontsize=9.5)
    source(fig, provenance() + "; reference points (1, 1, 1) for DTLZ1 and DTLZ2, (1, 1, 7) for DTLZ7, (1, 1) for ZDT")
    return save(fig, "search-22-convergence")


# ------------------------------------------------------------------------------------------------------------------
# search-23: how the population spreads over the directions, on the simplex.

def barycentric(p):
    s = p.sum(axis=1, keepdims=True)
    w = p / s
    x = w[:, 1] + 0.5 * w[:, 2]
    y = np.sqrt(3) / 2 * w[:, 2]
    return x, y


def figure_spread():
    P = populations()
    fin = finals()
    D = np.array([[float(r["w1"]), float(r["w2"]), float(r["w3"])] for r in rows("directions.csv")])
    dx, dy = barycentric(D)
    fig, axes = figure(2, 5)
    fig.subplots_adjust(left=0.04, right=0.99, top=0.79, bottom=0.10, wspace=0.05, hspace=0.12)
    title(fig, "Where each search puts its points among the 91 directions",
          "Final populations of the median runs, each point projected along its ray from the origin onto the simplex; rings mark the Das–Dennis directions (p = 12)")
    for row, prob in enumerate(["DTLZ1", "DTLZ2"]):
        for col, alg in enumerate(ALGORITHMS):
            ax = axes[row, col]
            ax.set_aspect("equal")
            ax.axis("off")
            ax.plot([0, 1, 0.5, 0], [0, 0, np.sqrt(3) / 2, 0], color=HAIRLINE, lw=0.8)
            ax.scatter(dx, dy, s=34, facecolors="none", edgecolors="#C9B8C3", linewidths=0.8, zorder=1)
            pts = P[(prob, alg)]
            near = pts[np.all(np.isfinite(pts), axis=1) & (pts.sum(axis=1) > 0)]
            x, y = barycentric(near)
            ax.scatter(x, y, s=10, color=COLOUR[alg], zorder=2, linewidths=0)
            runs = fin[(prob, alg)]
            occ = median([int(r["occupied"]) for r in runs])
            if row == 0:
                ax.set_title(LABEL.get(alg, alg), fontsize=13, color=COLOUR[alg], loc="left")
            ax.text(1.0, 0.86, f"{occ:.0f} of 91\noccupied", ha="right", va="top", fontsize=10, color=MUTED, transform=ax.transData)
        axes[row, 0].text(-0.06, 0.45, prob, rotation=90, fontsize=15, color=INK, ha="center", va="center", transform=axes[row, 0].transData)
    source(fig, provenance().split('; research')[0] + "; a direction is occupied when a point within 0.01 of the true front is nearest its ray (heuresis's NSGA-III association); median over 15 seeds")
    return save(fig, "search-23-spread-on-the-simplex")


if __name__ == "__main__":
    for f in (figure_fronts_dtlz, figure_disconnected, figure_convergence, figure_spread):
        print(f())
