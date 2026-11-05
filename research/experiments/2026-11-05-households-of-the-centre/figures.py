"""The figures of "The households of the centre": urban-science's housing on the City of Melbourne, drawn from results/
in the house style of the case studies, the place in the Urban Science Workbench's sea blue and sand.

    .venv/Scripts/python.exe research/experiments/2026-11-05-households-of-the-centre/figures.py [--only 10,11]

Writes assets/gallery/agents-10-... to agents-19-... (1600 x 1000).
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F401,F403

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"

SEA, SAND, DEEP = "#5B8FB9", "#C98B6B", "#1F3B57"
TYPE_COLOUR = ["#C2703F", "#3F7FB5", "#9A4A8A"]          # house, apartment, student (validated; markers differ too)
TYPE_NAME = ["houses and townhouses", "apartments", "student apartments"]
TYPE_MARK = ["s", "o", "^"]
CLASS_COLOUR = ["#C98420", "#3E9B6A", "#3366AA"]         # low, middle, high income
CLASS_NAME = ["low income (< A$650 a week)", "middle (A$650-1,999)", "high (A$2,000 or more)"]
BLOCK_FILL, BLOCK_EDGE = "#EFE7DF", "#D9CCD3"
LOW_MAP = LinearSegmentedColormap.from_list("low", ["#FBF6F1", "#F3DDB8", "#DFAA5E", "#B0741E", "#5E3A10"])
CHANGE_MAP = LinearSegmentedColormap.from_list("change", ["#2E5E86", "#8FB3D2", "#F4EEE8", "#E1A877", "#9C4F1F"])
SCEN_COLOUR = {"baseline": "#7A6F78", "towers under construction": "#3F7FB5", "every approved tower": "#C2703F", "student apartments": "#9A4A8A"}
CREDIT = "Data: CLUE and the Development Activity Monitor, City of Melbourne; Rental Report, Homes Victoria; Source: ABS, 2021 Census and ASGS (all CC BY 4.0)"


def load(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def km(v):
    return np.asarray(v, float) / 1000.0


def map_axes(ax, blocks, fill=None, values=None, cmap=None, norm=None, edge=BLOCK_EDGE, lw=0.4):
    polys = [km(b["ring"]) for b in blocks]
    if values is None:
        pc = PolyCollection(polys, facecolors=fill or BLOCK_FILL, edgecolors=edge, linewidths=lw)
    else:
        pc = PolyCollection(polys, array=np.asarray(values, float), cmap=cmap, norm=norm, edgecolors=edge, linewidths=lw)
    ax.add_collection(pc)
    xs = np.concatenate([p[:, 0] for p in polys])
    ys = np.concatenate([p[:, 1] for p in polys])
    ax.set_xlim(xs.min() - 0.1, xs.max() + 0.1)
    ax.set_ylim(ys.min() - 0.1, ys.max() + 0.1)
    ax.set_aspect("equal")
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    return pc


def scale_bar(ax, x0, y0, length=1.0):
    ax.plot([x0, x0 + length], [y0, y0], color=INK, lw=1.2)
    ax.text(x0 + length / 2, y0 - 0.12, f"{length:g} km", ha="center", va="top", fontsize=10, color=MUTED)


def quarter_years(labels):
    return np.array([int(q[:4]) + (int(q[-1]) - 0.5) / 4 for q in labels])


def nanarr(v):
    return np.array([np.nan if x is None else x for x in v], float)


# ---- 10 ----------------------------------------------------------------------------------------------------------

def fig10():
    st = load("stock.json")
    b = np.array(st["buildings"], float)
    blocks = st["blocks"]
    areas = st["areas"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.25, 1], left=0.03, right=0.97, top=0.86, bottom=0.08, wspace=0.22, hspace=0.38)
    ax = fig.add_subplot(gs[:, 0])
    map_axes(ax, blocks)
    for t in range(3):
        n = b[:, 7 + t]
        keep = n > 0
        ax.scatter(km(b[keep, 0]), km(b[keep, 1]), s=2.2 * np.sqrt(n[keep]) + 1, marker=TYPE_MARK[t], color=TYPE_COLOUR[t], alpha=0.75,
                   linewidths=0, label=f"{TYPE_NAME[t]}: {int(n.sum()):,}")
    an = km(st["anchor"])
    ax.plot(an[0], an[1], marker="+", color=INK, ms=12, mew=1.5)
    ax.text(an[0] + 0.12, an[1] + 0.1, "Town Hall", fontsize=10, color=INK)
    ax.legend(loc="lower left", fontsize=11, title="dwellings in 2024, by building", title_fontsize=11, markerscale=1.6)
    lim = ax.get_xlim(), ax.get_ylim()
    scale_bar(ax, lim[0][0] + 0.3, lim[1][1] - 0.3)
    # Right: the stock through time.
    tot = np.array(st["totals_by_year_type_area"], float)  # year × (type · area)
    A = len(areas)
    years = np.arange(2002, 2025)
    ax2 = fig.add_subplot(gs[0, 1])
    by_type = np.stack([tot[:, t * A:(t + 1) * A].sum(axis=1) for t in range(3)])
    ax2.stackplot(years, by_type[[1, 2, 0]], colors=[TYPE_COLOUR[1], TYPE_COLOUR[2], TYPE_COLOUR[0]], alpha=0.9, edgecolor=PAPER, linewidth=1.2)
    ax2.set_xlim(2002, 2024)
    ax2.set_xticks([2002, 2006, 2010, 2014, 2018, 2022, 2024])
    ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v)}"))
    ax2.set_ylabel("dwellings")
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v / 1000:.0f}k"))
    ax2.text(2003, by_type[1, 2] * 0.45, f"apartments: {int(by_type[1, 0]):,} to {int(by_type[1, -1]):,}", fontsize=11, color="white")
    top = by_type.sum(axis=0)
    ax2.annotate(f"student apartments: {int(by_type[2, 0]):,} to {int(by_type[2, -1]):,}", (2019, by_type[1, 17] + by_type[2, 17] / 2), xytext=(2003, 78000),
                 fontsize=10, color=INK, arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.8})
    ax2.annotate(f"houses and townhouses: {int(by_type[0, 0]):,} to {int(by_type[0, -1]):,}", (2017, top[15] - by_type[0, 15] / 2), xytext=(2003, 98000),
                 fontsize=10, color=INK, arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.8})
    ax2.set_title("the stock by type, CLUE's census years", loc="left", fontsize=13, color=INK)
    ax3 = fig.add_subplot(gs[1, 1])
    by_area0 = tot[0].reshape(3, A).sum(axis=0)
    by_area1 = tot[-1].reshape(3, A).sum(axis=0)
    order = [i for i in np.argsort(by_area1)[::-1] if by_area1[i] >= 50]
    y = np.arange(len(order))
    ax3.barh(y, by_area1[order], color=SEA, height=0.62, label="2024")
    ax3.barh(y, by_area0[order], color=DEEP, height=0.62, label="2002")
    short = {"West Melbourne (Residential)": "West Melbourne (res.)", "Melbourne (Remainder)": "Melbourne (rest)"}
    ax3.set_yticks(y)
    ax3.set_yticklabels([short.get(areas[i], areas[i]) for i in order], fontsize=10)
    ax3.invert_yaxis()
    ax3.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v / 1000:.0f}k"))
    ax3.grid(axis="y", visible=False)
    for i, k in enumerate(order):
        ax3.text(by_area1[k] + 400, i, f"×{by_area1[k] / max(by_area0[k], 1):.1f}", va="center", fontsize=9, color=MUTED)
    ax3.legend(loc="lower right", fontsize=10)
    ax3.set_title("dwellings by CLUE small area, 2002 and 2024", loc="left", fontsize=13, color=INK)
    title(fig, "The centre's dwellings", f"The City of Melbourne's {int(by_type[:, -1].sum()):,} dwellings of 2024 in {len(b):,} buildings on CLUE's {len(blocks)} blocks: "
          f"{int(by_type[:, 0].sum()):,} in 2002, apartments nearly all of the growth")
    source(fig, "CLUE residential dwellings by building and type, 2002-2024. " + CREDIT)
    save(fig, "agents-10-the-centres-dwellings")


# ---- 11 ----------------------------------------------------------------------------------------------------------

def two_way_demean(y, g, q, iters=200):
    r = y.astype(float).copy()
    for _ in range(iters):
        for key in (g, q):
            m = {k: r[key == k].mean() for k in np.unique(key)}
            r = r - np.array([m[k] for k in key])
    return r


def fig11():
    c = load("calibration.json")
    fit = c["inverse_demand"]
    rows = np.array(c["panel"]["rows"], float)
    seg, q, lr, ls = rows[:, 0].astype(int), rows[:, 1].astype(int), rows[:, 2], rows[:, 3]
    ry, rx = two_way_demean(lr, seg, q), two_way_demean(ls, seg, q)
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.6, 1], left=0.06, right=0.97, top=0.84, bottom=0.11, wspace=0.2, hspace=0.32)
    ax = fig.add_subplot(gs[:, 0])
    names = fit["segments"]
    palette = ["#3F7FB5", "#9A4A8A", "#3E9B6A", "#C2703F", "#E3A86B", "#C98420", "#8C6A4F", "#7A6F78"]
    for k, name in enumerate(names):
        m = seg == k
        ax.scatter(rx[m], ry[m], s=14, color=palette[k % len(palette)], alpha=0.7, linewidths=0, label=name)
    xx = np.linspace(rx.min(), rx.max(), 10)
    ax.plot(xx, fit["slope"] * xx, color=INK, lw=2)
    ax.text(-0.95, 0.19, f"slope {fit['slope']:.3f}, clustered s.e. {fit['se_cluster']:.3f}\nα = −1/slope = {fit['alpha']:.2f} ± {fit['alpha_se']:.2f}", fontsize=12,
            color=INK, va="top")
    ax.set_xlabel("log stock, less its segment's and its quarter's means")
    ax.set_ylabel("log rent of new lettings, less the same")
    ax.legend(fontsize=9, loc="lower left", ncol=2, markerscale=1.5)
    ax.set_title("the inverse demand of eight segments, 2003Q1-2025Q3", loc="left", fontsize=13)
    qs = c["quarters"]
    t = quarter_years(qs)
    d = nanarr(c["demand_index"])
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(t, d, color=SEA, lw=2.2)
    ax.set_yscale("log")
    ax.set_yticks([0.25, 0.5, 1, 2])
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_ylabel("2003Q1 = 1")
    ax.axvspan(2020.2, 2022.0, color=HAIRLINE, alpha=0.5, lw=0)
    ax.text(2020.25, 2.25, "borders closed to\ninternational students", fontsize=10, color=MUTED, va="top")
    i19, i21 = qs.index("2019Q3"), qs.index("2021Q3")
    ax.annotate(f"{d[i21]:.2f}", (t[i21], d[i21]), xytext=(t[i21] - 3.2, d[i21] * 0.95), fontsize=11, color=INK, arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.8})
    ax.annotate(f"{d[i19]:.2f}", (t[i19], d[i19]), xytext=(t[i19] - 3.2, d[i19] * 1.1), fontsize=11, color=INK, arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.8})
    ax.set_xlim(2003, 2025.8)
    ax.set_ylim(0.2, 2.6)
    ax.set_title("the demand the rents reveal: exp(α (quarter effect − log metropolitan rent))", loc="left", fontsize=12)
    ax = fig.add_subplot(gs[1, 1])
    metro = nanarr(c["metro_index"])
    ax.plot(t, metro, color=SAND, lw=2)
    ax.set_xlim(2003, 2025.8)
    ax.set_ylabel("A$ a week")
    ax.set_title("the outside option's rent: the metropolis's other regions", loc="left", fontsize=12)
    title(fig, "What rents reveal", f"A segment whose stock grows against the others' lets for less: {abs(fit['slope']):.2f} % for each 1 % of stock, a semi-elasticity "
          f"of choice in rent α = {fit['alpha']:.1f}")
    source(fig, f"Two-way fixed effects on {int(fit['observations'])} segment-quarters (tyche regression, HC1 and segment-clustered errors). " + CREDIT)
    save(fig, "agents-11-what-rents-reveal")


# ---- 12 ----------------------------------------------------------------------------------------------------------

def fig12():
    h = load("history.json")
    t = quarter_years(h["ticks"])
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), dpi=100, sharex=True, gridspec_kw={"left": 0.06, "right": 0.98, "top": 0.84, "bottom": 0.1, "wspace": 0.14, "hspace": 0.28})
    for ax, g in zip(axes.flat, h["groups"]):
        ax.plot(t, nanarr(g["observed"]), color=INK, lw=2.2, label="the report (new lettings)")
        for k, m in enumerate(g["model"]):
            ax.plot(t, nanarr(m), color=SEA, lw=1.0, alpha=0.75, label="the model, four seeds" if k == 0 else None)
        ax.set_title(g["group"].replace("-", "–"), loc="left", fontsize=13)
        ax.set_xlim(2002.5, 2025.8)
    axes[0, 0].legend(fontsize=10, loc="upper left")
    for ax in axes[:, 0]:
        ax.set_ylabel("median rent of new lettings, A$ a week")
    rm = np.array(h["rmse"])
    title(fig, "Rents of new lettings, reported and simulated", f"Each suburb group's flats, moving annual medians: the agents calibrated at 2002's rents, run on CLUE's stock and "
          f"the revealed demand (log error: centre {rm[:, 0].mean():.3f}, groups {rm[:, 1].mean():.3f})")
    source(fig, f"κ = {h['kappa']}, π = {h['power']}; newcomers {h['arrivals0']:.0f} a quarter in 2002. " + CREDIT)
    save(fig, "agents-12-rents-reported-and-simulated")


# ---- 13 ----------------------------------------------------------------------------------------------------------

def fig13():
    c = load("census-2021.json")
    rows = c["rows"]
    rings = [km(r["ring"]) for r in rows]
    cen = np.array([r["census"] for r in rows], float)
    mod = np.array([[m for m in r["model"]] for r in rows], float)  # sa1 × seed × col
    mod = np.nan_to_num(mod)
    mm = mod.mean(axis=1)
    def share(a, i, cols=(0, 1, 2)):
        tot = a[:, list(cols)].sum(axis=1)
        return np.where(tot > 0, a[:, i] / np.maximum(tot, 1e-9), np.nan), tot
    low_c, nc = share(cen, 0)
    low_m, nm = share(mm, 0)
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1, 1, 0.75], left=0.02, right=0.98, top=0.86, bottom=0.14, wspace=0.12, hspace=0.38)
    norm = Normalize(0, 0.6)
    for k, (vals, n, label) in enumerate([(low_c, nc, "the 2021 Census"), (low_m, nm, "the model, 2021Q3 (four seeds)")]):
        ax = fig.add_subplot(gs[:, k])
        v = np.where(n >= 20, vals, np.nan)
        pc = PolyCollection(rings, array=np.ma.masked_invalid(v), cmap=LOW_MAP, norm=norm, edgecolors="white", linewidths=0.3)
        pc.cmap.set_bad("#EEE7E1")
        ax.add_collection(pc)
        xs = np.concatenate([p[:, 0] for p in rings]); ys = np.concatenate([p[:, 1] for p in rings])
        ax.set_xlim(xs.min(), xs.max()); ax.set_ylim(ys.min(), ys.max()); ax.set_aspect("equal"); ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(label, loc="left", fontsize=13)
        if k == 0:
            scale_bar(ax, xs.min() + 0.2, ys.min() + 0.3)
    cax = fig.add_axes([0.12, 0.1, 0.3, 0.015])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=LOW_MAP), cax=cax, orientation="horizontal")
    cb.set_label("share of households with a low income (SA1s of 20 households or more)", fontsize=10, color=INK)
    cb.outline.set_visible(False)
    # Right: by small area.
    areas = sorted({r["area"] for r in rows})
    ax = fig.add_subplot(gs[0, 2])
    ax2 = fig.add_subplot(gs[1, 2])
    for a in areas:
        m = np.array([r["area"] == a for r in rows])
        cs, ms = cen[m].sum(axis=0), mm[m].sum(axis=0)
        if cs[:3].sum() < 200:
            continue
        for i, col in enumerate(CLASS_COLOUR):
            ax.scatter(cs[i] / cs[:3].sum(), ms[i] / ms[:3].sum(), color=col, s=26, linewidths=0, alpha=0.85)
        ax2.scatter(cs[3] / (cs[3] + cs[4]), ms[3] / max(ms[3] + ms[4], 1), color=SEA, s=30, linewidths=0)
        pass
    for a_ in (ax, ax2):
        a_.plot([0, 1], [0, 1], color=HAIRLINE, lw=1)
        a_.set_xlim(0, 1); a_.set_ylim(0, 1)
    ax.set_xlabel("the Census", fontsize=11, labelpad=2); ax.set_ylabel("the model", fontsize=11)
    ax.set_title("each small area's income classes", loc="left", fontsize=12)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=col, label=n.split(" (")[0]) for col, n in zip(CLASS_COLOUR, CLASS_NAME)], fontsize=9, loc="upper left")
    ax2.set_xlabel("the Census", fontsize=11); ax2.set_ylabel("the model", fontsize=11)
    ax2.set_title("each small area's share of renters", loc="left", fontsize=12)
    ax2.set_xlim(0.3, 1); ax2.set_ylim(0.3, 1)
    dm = np.array(c["d_model"])
    title(fig, "The 2021 Census and the model", f"Where the low incomes live, by SA1: dissimilarity of low from high incomes {c['d_census']:.3f} in the Census, "
          f"{dm.mean():.3f} ± {dm.std():.3f} in the model, its income power π fitted to it")
    source(fig, "The Census of 10 August 2021, taken in a lockdown, against the model's 2021Q3. " + CREDIT)
    save(fig, "agents-13-the-census-and-the-model")


# ---- 14, 15 ------------------------------------------------------------------------------------------------------

def fig14():
    s = load("shocks.json")
    t = quarter_years(s["ticks"])
    fig, axes = plt.subplots(1, 3, figsize=(16, 10), dpi=100, gridspec_kw={"left": 0.06, "right": 0.98, "top": 0.76, "bottom": 0.09, "wspace": 0.25})
    base = s["scenarios"][0]
    for sc in s["scenarios"]:
        col = SCEN_COLOUR[sc["name"]]
        ls = "--" if sc["name"] == "baseline" else "-"
        cen = np.array([nanarr(x["centre"]) for x in sc["seeds"]])
        stu = np.array([nanarr(x["student"]) for x in sc["seeds"]])
        res = np.array([np.array(x["residents"]).sum(axis=1) for x in sc["seeds"]])
        low = np.array([np.array(x["residents"])[:, 0] for x in sc["seeds"]])
        lab = sc["name"] + ("" if sc["name"] == "baseline" else f" (+{sc['added_apartments']:,.0f} apartments, +{sc['added_students']:,.0f} student)")
        axes[0].plot(t, np.nanmean(cen, axis=0), color=col, lw=2, ls=ls, label=lab)
        axes[0].fill_between(t, np.nanmin(cen, axis=0), np.nanmax(cen, axis=0), color=col, alpha=0.12, lw=0)
        axes[1].plot(t, np.nanmean(stu, axis=0), color=col, lw=2, ls=ls)
        axes[2].plot(t, np.nanmean(low, axis=0), color=col, lw=2, ls=ls)
    axes[0].set_ylabel("the centre's flats: median rent of new lettings, A$ a week")
    axes[1].set_ylabel("student apartments' rent, A$ a week")
    axes[2].set_ylabel("low-income households resident")
    axes[0].set_title("rents of new lettings", loc="left", fontsize=13)
    axes[1].set_title("student apartments", loc="left", fontsize=13)
    axes[2].set_title("who stays: low incomes", loc="left", fontsize=13)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(0.055, 0.895), ncol=2, fontsize=11)
    for ax in axes:
        ax.set_xlim(t[0], t[-1])
    title(fig, "A supply shock, in time", "The Development Activity Monitor's pipeline completed in the model, 2026-2030, against no new supply beyond 2026's completions: "
          "mean of four seeds, the band their range")
    source(fig, f"Demand and the metropolis's rents carried at their 2010-2019 growth ({s['g_demand']:+.3f} and {s['g_metro']:+.3f} a year, assumptions). CLUE, DAM, Homes Victoria, ABS (CC BY 4.0).")
    save(fig, "agents-14-a-supply-shock-in-time")


def fig15():
    s = load("shocks.json")
    blocks = s["blocks"]
    base = np.array([np.array(x["blocks"]) for x in s["scenarios"][0]["seeds"]]).mean(axis=0)
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 1, 0.55], left=0.02, right=0.98, top=0.84, bottom=0.1, wspace=0.08)
    for k, name in enumerate(["every approved tower", "student apartments"]):
        sc = next(x for x in s["scenarios"] if x["name"] == name)
        b = np.array([np.array(x["blocks"]) for x in sc["seeds"]]).mean(axis=0)
        d = b[:, 0] - base[:, 0]
        ax = fig.add_subplot(gs[0, k])
        lim = max(5, np.percentile(np.abs(d), 99))
        map_axes(ax, blocks, values=d, cmap=CHANGE_MAP, norm=TwoSlopeNorm(0, -lim, lim), edge="white", lw=0.3)
        devs = [x for x in s["developments"] if x["status"] in ("UNDER CONSTRUCTION", "APPROVED")]
        key = "dwellings" if k == 0 else "student"
        xs = km([x["x"] for x in devs]); ys = km([x["y"] for x in devs]); n = np.array([x[key] for x in devs])
        m = n > 0
        ax.scatter(xs[m], ys[m], s=np.sqrt(n[m]) * 3, facecolors="none", edgecolors=INK, linewidths=0.8)
        ax.set_title(f"{name}: low-income households, against the baseline, 2030Q4", loc="left", fontsize=12)
        cax = ax.inset_axes([0.05, 0.03, 0.4, 0.02])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=TwoSlopeNorm(0, -lim, lim), cmap=CHANGE_MAP), cax=cax, orientation="horizontal")
        cb.set_label("households a block (circles: the pipeline's buildings)", fontsize=9)
        cb.outline.set_visible(False)
    ax = fig.add_subplot(gs[0, 2])
    names = [x["name"] for x in s["scenarios"][1:]]
    y = np.arange(len(names))
    left = np.zeros(len(names))
    for i in range(3):
        v = np.array([np.mean([x["in_new"][i] for x in sc["seeds"]]) for sc in s["scenarios"][1:]])
        ax.barh(y, v, left=left, color=CLASS_COLOUR[i], height=0.55, label=CLASS_NAME[i], edgecolor=PAPER, linewidth=1.5)
        left += v
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=11); ax.invert_yaxis(); ax.grid(axis="y", visible=False)
    ax.set_xlabel("households in the pipeline's buildings, 2030Q4")
    ax.legend(fontsize=9, loc="lower right")
    ax.set_title("who moves into the new buildings", loc="left", fontsize=12)
    title(fig, "A supply shock, in space", "Where the low incomes are held or lost when the towers or the student apartments are built (mean of four seeds)")
    source(fig, CREDIT)
    save(fig, "agents-15-a-supply-shock-in-space")


# ---- 16, 17 ------------------------------------------------------------------------------------------------------

def fig16():
    s = load("schelling.json")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = GridSpec(2, 4, figure=fig, height_ratios=[1, 1.1], left=0.05, right=0.98, top=0.84, bottom=0.06, wspace=0.12, hspace=0.3)
    ax = fig.add_subplot(gs[0, :])
    cols = [SEA, "#9A4A8A"]
    for L, col in zip(s["lattices"], cols):
        f = L["threshold"]
        n = np.array(L["segregation_nomos"]); h = np.array(L["segregation_housing"])
        ax.plot(np.arange(len(n)), n, color=col, lw=2.2, label=f"nomos's Schelling, f = {f:.2f}")
        ax.plot(np.arange(len(h)), h, ls="", marker="o", ms=6, mfc="none", mec=INK, mew=0.9, label=f"the households in Schelling's limit, f = {f:.2f}")
        ax.text(len(n) - 0.5, n[-1] + 0.01, f"{n[-1]:.3f}; {int(L['ticks_differing'])} of 40 ticks differ", ha="right", fontsize=11, color=INK)
    ax.set_xlabel("tick"); ax.set_ylabel("mean share of like neighbours")
    ax.legend(fontsize=10, loc="lower right", ncol=2)
    ax.set_title("one trajectory, two programs: draw for draw", loc="left", fontsize=13)
    cmap = LinearSegmentedColormap.from_list("two", ["#F4EEE8", "#C98420", "#3366AA"])
    for k, (L, which) in enumerate([(s["lattices"][0], "start"), (s["lattices"][0], "end"), (s["lattices"][1], "start"), (s["lattices"][1], "end")]):
        a = fig.add_subplot(gs[1, k])
        w = int(L["width"])
        grid = np.array(L[which], float).reshape(w, w) + 1
        a.imshow(grid, cmap=cmap, vmin=0, vmax=2, interpolation="nearest")
        a.set_xticks([]); a.set_yticks([]); a.grid(False)
        for sp in a.spines.values():
            sp.set_visible(False)
        a.set_title(f"f = {L['threshold']:.2f}, {'placed' if which == 'start' else 'after 40 ticks'}", loc="left", fontsize=12)
    title(fig, "Schelling's limit, exactly", "With no rent, no type, no distance, a step utility at the threshold f and the best dwelling taken, the households' model on a 40 × 40 torus "
          "is nomos's Schelling, tick by tick")
    source(fig, "1,200 households of two classes on 1,600 cells, Moore's neighbourhood, random sequential order, seed 2026 (urban::housing, nomos::agents).")
    save(fig, "agents-16-schellings-limit-exactly")


def fig17():
    s = load("schelling.json")
    st = load("stock.json")
    blocks = st["blocks"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.3, 1], left=0.06, right=0.98, top=0.86, bottom=0.1, wspace=0.08, hspace=0.42)
    norm = Normalize(0, 0.8)
    panels = [(np.array(s["model_blocks"]), "the calibrated model, 2021Q3")] + [(np.array(L["blocks"]), f"Schelling's limit by block, f = {L['threshold']:.2f}") for L in s["city_limits"]]
    for k, (b, label) in enumerate(panels):
        ax = fig.add_subplot(gs[0, k])
        tot = b.sum(axis=1)
        v = np.where(tot >= 10, b[:, 0] / np.maximum(tot, 1), np.nan)
        pc = map_axes(ax, blocks, values=np.ma.masked_invalid(v), cmap=LOW_MAP, norm=norm, edge="white", lw=0.25)
        pc.cmap.set_bad("#EEE7E1")
        ax.set_title(label, loc="left", fontsize=12)
    cax = fig.add_axes([0.66, 0.47, 0.3, 0.012])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=LOW_MAP), cax=cax, orientation="horizontal")
    cb.set_label("share of a block's households with a low income", fontsize=9); cb.outline.set_visible(False)
    ax = fig.add_subplot(gs[1, :2])
    et = s["etas"]
    ax.plot([e["eta"] for e in et], [e["d"] for e in et], color=SEA, lw=2, marker="o", label="the calibrated model with a taste η for like neighbours")
    ax.axhline(s["d_census"], color=INK, lw=1.2, ls="--")
    ax.text(0.1, s["d_census"] - 0.06, f"the 2021 Census {s['d_census']:.3f}", fontsize=11, color=INK)
    for L, col in zip(s["city_limits"], [SAND, "#9A4A8A"]):
        ax.axhline(L["d"][-1], color=col, lw=1.2)
        ax.text(0.1, L["d"][-1] + 0.015, f"Schelling's limit, f = {L['threshold']:.2f}: {L['d'][-1]:.3f}", fontsize=11, color=col)
    ax.set_xlabel("η, utility of an all-alike block"); ax.set_ylabel("dissimilarity of low from high incomes, SA1s")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=10, loc="lower right")
    ax.set_title("how far the city is from the limit", loc="left", fontsize=12)
    lim = s["city_limits"]
    ax = fig.add_subplot(gs[1, 2])
    for L, col in zip(s["city_limits"], [SAND, "#9A4A8A"]):
        ax.plot(L["d"], color=col, lw=2, label=f"D, f = {L['threshold']:.2f}")
        ax.plot(L["content"], color=col, lw=1.2, ls=":", label=f"content, f = {L['threshold']:.2f}")
    ax.set_xlabel("tick of the limit"); ax.set_ylim(0, 1.02); ax.set_xticks(range(len(lim[0]["d"]) if False else 4))
    ax.legend(fontsize=9, loc="lower right")
    ax.set_title("the limit's own path", loc="left", fontsize=12)
    lim = s["city_limits"]
    frac = (s["d_census"] - s["d_model"]) if False else None
    title(fig, "How far the city is from Schelling's limit", f"Sorting by price alone gives the Census's dissimilarity {s['d_census']:.3f} (model {s['d_model']:.3f}); "
          f"sorting by neighbours alone, from 2021's households, gives {lim[0]['d'][-1]:.3f} at f = 1/3 and {lim[1]['d'][-1]:.3f} at f = 1/2")
    source(fig, "Blocks as neighbourhoods; three income classes; random sequential order. " + CREDIT)
    save(fig, "agents-17-how-far-from-schellings-limit")


# ---- 18, 19 ------------------------------------------------------------------------------------------------------

def fig18():
    p = load("planner.json")
    X, Y = int(p["X"]), int(p["Y"])
    tab = np.array(p["table"], float)
    rent = tab[:, 2].reshape(X, Y)
    plan = np.array(p["plan"], float)
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.05, 1], height_ratios=[1, 1], left=0.07, right=0.97, top=0.84, bottom=0.12, wspace=0.28, hspace=0.42)
    ax = fig.add_subplot(gs[:, 0])
    cmap = LinearSegmentedColormap.from_list("r", ["#E8EEF3", "#A9C3DA", "#5B8FB9", "#2E5E86", "#1F3B57"])
    im = ax.imshow(rent, origin="lower", cmap=cmap, aspect="auto")
    for x in range(X):
        for y in range(Y):
            a = int(plan[x, y])
            ax.text(y, x, ["·", "1", "2"][a], ha="center", va="center", fontsize=14 if a else 16, color="white" if rent[x, y] > np.median(rent) else INK)
    u = p["unit"]
    ax.set_xticks(range(Y)); ax.set_xticklabels([f"{(y - int(p['y_now'])) * 10:+d}%" for y in range(Y)], fontsize=9)
    ax.set_yticks(range(X)); ax.set_yticklabels([f"+{x * u / 1000:.1f}k" for x in range(X)], fontsize=9)
    ax.set_xlabel("demand: the stock it would fill at 2025Q3's rent, against 2024's apartments")
    ax.set_ylabel("apartments beyond 2024's stock")
    ax.grid(False)
    ax.set_title("HSVI's plan where demand is known: units permitted (cells: the clearing rent, A$ a week)", loc="left", fontsize=12)
    cax = ax.inset_axes([1.02, 0.0, 0.025, 0.6])
    cb = fig.colorbar(im, cax=cax, orientation="vertical")
    cb.outline.set_visible(False)
    pol = p["policies"]
    names = [q["name"] for q in pol]
    vals = np.array([q["value"] for q in pol]); err = np.array([q["error"] for q in pol])
    colours = ["#7A6F78", "#C2703F", "#C98420", "#3E9B6A", "#3F7FB5", "#9A4A8A"]
    yy = np.arange(len(names))
    ax = fig.add_subplot(gs[0, 1])
    ax.barh(yy, vals / 1000, color=colours, height=0.55, edgecolor=PAPER, linewidth=1.5)
    for i, q in enumerate(pol):
        a = q["actions"]
        ax.text(vals[i] / 1000 + 0.6, i, f"{vals[i] / 1000:.2f}   (permits 0/1/2 a year: {a[0]:.2f} / {a[1]:.2f} / {a[2]:.2f})", va="center", fontsize=9, color=INK)
    ax.set_yticks(yy); ax.set_yticklabels(names, fontsize=11); ax.invert_yaxis(); ax.grid(axis="y", visible=False)
    ax.set_xlim(0, vals.max() / 1000 * 1.75)
    ax.set_xlabel("discounted welfare, A$ billion")
    ax.set_title("six policies, simulated (3,000 episodes of 150 years)", loc="left", fontsize=13)
    ax = fig.add_subplot(gs[1, 1])
    keep = [0, 3, 4, 5]
    ref = vals[0]
    for k, i in enumerate(keep):
        ax.errorbar((vals[i] - ref), k, xerr=2 * err[i], fmt="o", color=colours[i], ms=8, capsize=4, lw=1.4)
        ax.text((vals[i] - ref) + 2 * err[i] + 12, k, f"{vals[i] - ref:+.0f}", va="center", fontsize=10, color=INK)
    lo, hi = p["hsvi"][0] - ref, p["hsvi"][1] - ref
    ax.axvspan(lo, hi, color="#9A4A8A", alpha=0.25, lw=0)
    ax.text(lo - 8, -0.55, f"HSVI's interval for the optimum: [{lo:+.1f}, {hi:+.1f}]", fontsize=10, color="#9A4A8A", ha="right")
    ax.set_yticks(range(len(keep))); ax.set_yticklabels([names[i] for i in keep], fontsize=11); ax.invert_yaxis(); ax.grid(axis="y", visible=False)
    ax.set_ylim(len(keep) - 0.5, -0.9)
    ax.set_xlabel("against never permitting, A$ million (two standard errors)")
    ax.set_title("the best of them, closer", loc="left", fontsize=13)
    title(fig, "The planner's permits under uncertain demand", f"A MOMDP: the stock seen, demand hidden and read through the band of the rent; a unit is "
          f"{u:,.0f} apartments at A$ {p['cost'] / 1000:.2f} billion, completed the next year")
    source(fig, f"Welfare: the households' surplus and the rents, from urban::housing::market; demand's yearly steps from the revealed demand; discount {p['discount']}. "
           "Costs, bands and the reading's error are assumptions.")
    save(fig, "agents-18-the-planners-permits")


def fig19():
    r = load("rent-spectra.json")
    p = load("planner.json")
    fig, axes = plt.subplots(1, 2, figsize=(16, 10), dpi=100, gridspec_kw={"left": 0.06, "right": 0.98, "top": 0.83, "bottom": 0.11, "wspace": 0.22})
    ax = axes[0]
    for k, s in enumerate(r["spectra"]):
        f = np.array(s["spectrum"]["frequency"]); d = np.array(s["spectrum"]["density"])
        if k == 0:
            ax.plot(f[1:], d[1:], color=INK, lw=2.2, label=s["name"])
        else:
            ax.plot(f[1:], d[1:], color=SEA, lw=1, alpha=0.8, label="the model, four seeds" if k == 1 else None)
    ax.set_yscale("log")
    ax.set_xlabel("cycles a year"); ax.set_ylabel("spectral density of the quarterly growth of log rent")
    ax.axvline(1, color=HAIRLINE, lw=1)
    ax.text(1.02, ax.get_ylim()[1] * 0.5 if ax.get_ylim()[1] > 0 else 1, "the year", fontsize=10, color=MUTED)
    ax.legend(fontsize=11, loc="lower left")
    ax.set_title("the centre's rents, 2003-2025: reported and simulated", loc="left", fontsize=13)
    ax = axes[1]
    colours = ["#7A6F78", "#C2703F", "#C98420", "#3E9B6A", "#3F7FB5", "#9A4A8A"]
    for q, col in zip(p["policies"], colours):
        f = np.array(q["frequency"]); d = np.array(q["density"])
        ax.plot(1 / f[1:], d[1:], color=col, lw=1.8, label=f"{q['name']}: variance {d[1:].sum() * (f[1] - f[0]):.5f}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("period, years"); ax.set_ylabel("spectral density of log rent")
    ax.legend(fontsize=10, loc="lower right")
    ax.set_title("the cycles of rent each planner's policy makes", loc="left", fontsize=13)
    title(fig, "Rents and their cycles", "Left: multitaper spectra of the quarterly growth of the centre's rents. Right: the periodogram of the clearing rent, year by year, under each "
          "policy of the planner (400 simulated paths of 100 years)")
    source(fig, "Thomson's multitaper (NW 3, five tapers) and the periodogram, calx signal. " + CREDIT)
    save(fig, "agents-19-rents-and-their-cycles")


FIGURES = {10: fig10, 11: fig11, 12: fig12, 13: fig13, 14: fig14, 15: fig15, 16: fig16, 17: fig17, 18: fig18, 19: fig19}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    which = [int(x) for x in args.only.split(",") if x] or list(FIGURES)
    for k in which:
        if k in FIGURES:
            try:
                FIGURES[k]()
                print(f"agents-{k}: drawn")
            except FileNotFoundError as e:
                print(f"agents-{k}: skipped ({e.filename} missing)")
