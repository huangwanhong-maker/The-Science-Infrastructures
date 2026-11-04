"""The figures of the case study "A district's economy without a model", drawn from results/ (and the shared problem's
results), in the house style of the search case studies (research/experiments/style/search_style.py), the place in
the Urban Science Workbench's sea blue and sand.

    .venv/Scripts/python.exe research/experiments/2026-11-05-a-district-economy-without-a-model/figures.py [--only 62,63]

Writes assets/gallery/search-62-... to search-69-... (1600 x 1000). The maps need the shared problem's uncommitted
results/district-geometry.json (made by its calibrate from the data library's cache).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import sys
from collections import defaultdict

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
PROBLEM = HERE.parent / "2026-11-05-a-district-economy-the-problem"
sys.path.insert(0, str(HERE.parent / "style"))
from search_style import *  # noqa: E402,F401,F403

import importlib.util  # noqa: E402

# The shared problem's drawing of the district (its figures.py), loaded under its own name.
_spec = importlib.util.spec_from_file_location("district_economy_problem_figures", PROBLEM / "figures.py")
_problem = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_problem)
Frame, KIND_COLOUR, KINDS, SEA, SAND = _problem.Frame, _problem.KIND_COLOUR, _problem.KINDS, _problem.SEA, _problem.SAND
district_axes, draw_base, north_arrow, street_labels = _problem.district_axes, _problem.draw_base, _problem.north_arrow, _problem.street_labels

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize, TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

RESULTS = HERE / "results"
PRESULTS = PROBLEM / "results"
PROVENANCE = "Data: City of Melbourne, Pedestrian Counting System and CLUE (CC BY 4.0); OpenStreetMap contributors (ODbL 1.0)"
FLOW = LinearSegmentedColormap.from_list("flow", ["#EEF2F6", "#A9C3DA", "#5B8FB9", "#2F5A86", "#1F2A44"])
ALG = {"GA": ALGORITHM["GA"], "GA from as built": "#4A6B46", "CMA-ES": ALGORITHM["CMA-ES"], "MAP-Elites": "#4E8F8B", "NSGA-II": ALGORITHM["MOEA/D"], "random": ALGORITHM["random"], "BO (REMBO)": ALGORITHM["BO"]}
ORDER = ["GA", "CMA-ES", "MAP-Elites", "NSGA-II", "random", "BO (REMBO)", "GA from as built"]
LARGE = ["GA", "CMA-ES", "MAP-Elites", "NSGA-II", "random", "GA from as built"]
DAYS = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]


def load(name, where=RESULTS):
    return json.loads((where / name).read_text(encoding="utf-8"))


def read_csv(name):
    with (RESULTS / name).open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def week_axis(ax, label=True):
    ax.set_xlim(0, 168)
    ax.set_xticks(np.arange(0, 169, 24), [""] * 8)
    if label:
        for d, nm in enumerate(DAYS):
            ax.text(24 * d + 12, -0.05, nm, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=10.5, color=MUTED)


def short_name(a):
    """An algorithm's name for a tick: the warm-started GA on two lines."""
    return "GA from\nas built" if a == "GA from as built" else a


def plans_by_key():
    return {p["key"]: p for p in load("plans.json")["plans"]}


# ---- 62: the people without a model -------------------------------------------------------------------------------

def figure_62():
    G, S, people = load("district-geometry.json", PRESULTS), load("sensors.json", PRESULTS), load("people.json")
    T = Frame(G["quad"])
    fig = plt.figure(figsize=(16, 10), dpi=100)
    k = people["kriging"]
    title(fig, "The people without a model: footfall interpolated from the counters",
          f"Kriging ln(1 + count) in the plane by tyche's Gaussian conditioning: an exponential covariance of length {k['length_m']:.0f} m, half of it a nugget, chosen by the likelihood")
    seg = {s["id"]: s for s in G["segments"]}
    hours = people["segments_hours"]
    norm = LogNorm(vmin=30, vmax=5000)
    sens_h = [8, 13, 19, 5 * 24 + 14]
    boxes = [[0.02, 0.52, 0.33, 0.36], [0.36, 0.52, 0.33, 0.36], [0.02, 0.135, 0.33, 0.36], [0.36, 0.135, 0.33, 0.36]]
    for q in range(4):
        ax = fig.add_axes(boxes[q])
        draw_base(ax, T, G, ["#EFE9E3"] * len(G["frontages"]), edge="#E2D9D3")
        lines, vals = [], []
        for sid, v in people["segments"].items():
            if sid in seg:
                lines.append(T(seg[sid]["line"])), vals.append(max(1.0, v[q]))
        widths = [1.2 + 2.6 * (math.log10(max(v, 30)) - 1.5) for v in vals]
        ax.add_collection(LineCollection(lines, colors=FLOW(norm(vals)), linewidths=widths, zorder=4, capstyle="round"))
        district_axes(ax, T, G, pad=(60, 30))
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        for s in S["sensors"]:
            p = T([[s["x"], s["y"]]])[0]
            if xlim[0] <= p[0] <= xlim[1] and ylim[0] <= p[1] <= ylim[1]:
                ax.scatter([p[0]], [p[1]], marker="s", s=60, color=FLOW(norm(max(1, s["mean_week"][sens_h[q]]))), edgecolors=INK, linewidths=0.9, zorder=6)
        ax.set_title(hours[q], fontsize=13, loc="left", color=INK)
    cax = fig.add_axes([0.06, 0.085, 0.25, 0.013])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=FLOW)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=10)
    fig.text(0.33, 0.091, "people an hour: kriged at the segments' middles (lines), measured at the counters (squares)", fontsize=10.5, color=INK, va="center")
    cb.outline.set_visible(False)

    # Leave one out.
    ax = fig.add_axes([0.76, 0.52, 0.21, 0.33])
    ob = np.array([x["observed_week"] for x in k["loo"]])
    pr = np.array([x["predicted_week"] for x in k["loo"]])
    ax.loglog([2e3, 5e5], [2e3, 5e5], color=HAIRLINE, lw=1)
    ax.scatter(ob, pr, s=26, color=SEA, edgecolors=PAPER, linewidths=0.6)
    ax.set_xlim(2e3, 5e5), ax.set_ylim(2e3, 5e5)
    ax.set_xlabel("observed, people a week")
    ax.set_ylabel("predicted from the others")
    ax.set_title("each counter left out in turn", fontsize=12.5, loc="left")
    ax.text(0.04, 0.96, f"RMSE {k['loo_rmse']:.2f} in ln(1 + count)\nR² {k['loo_r2']:.2f} over all hours\nmedian ratio of weeks {k['loo_weekly_median_ratio']:.2f}",
            transform=ax.transAxes, va="top", fontsize=10.5, color=INK)
    # The likelihood over the grid.
    ax2 = fig.add_axes([0.76, 0.12, 0.21, 0.28])
    grid = np.array(k["grid"])
    for eta in sorted(set(grid[:, 1])):
        g = grid[grid[:, 1] == eta]
        ax2.semilogx(g[:, 0], g[:, 2] - grid[:, 2].max(), color=ORCHID(0.25 + eta), lw=1.4, label=f"η = {eta:.2f}" if eta in (0.0, 0.2, 0.5) else None)
    ax2.set_ylim(-400, 20)
    ax2.axvline(k["length_m"], color=INK, lw=0.8, ls=(0, (3, 2)))
    ax2.set_xlabel("length ℓ of the covariance (m)")
    ax2.set_ylabel("log-likelihood − its maximum")
    ax2.set_title("the profile likelihood, 168 hours", fontsize=12.5, loc="left")
    ax2.legend(fontsize=9.5, loc="lower right")
    source(fig, PROVENANCE + ". Hours of the mean week of the reference year, 27 October 2025 to 25 October 2026; Flinders Lane at the left, Bourke Street at the right.")
    save(fig, "search-62-people-without-a-model")


# ---- 63: the searches compared ------------------------------------------------------------------------------------

def spaced(values, gap):
    """Label positions near the values, at least `gap` apart, in the values' order."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = list(values)
    for a, b in zip(order, order[1:]):
        if out[b] - out[a] < gap:
            out[b] = out[a] + gap
    return out


def figure_63():
    rows = read_csv("anytime.csv")
    summ = load("summary.json")
    v0 = summ["as_built"]["value_aud_week"] / 1000
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Six searches of one black box, at equal budgets, and a GA started from the plan as built",
          "The best weekly value found against the evaluations spent, five seeds each: the median and the range; the plan as built dashed")
    for b, box, label in ((20000, [0.07, 0.12, 0.52, 0.72], "20,000 evaluations"), (400, [0.68, 0.12, 0.27, 0.72], "400 evaluations")):
        ax = fig.add_axes(box)
        by = defaultdict(lambda: defaultdict(list))
        for r in rows:
            if int(r["budget"]) != b:
                continue
            by[r["algorithm"]][int(r["evaluations"])].append(float(r["best_value"]) / 1000)
        ends, names = [], []
        for alg in ORDER:
            if alg not in by:
                continue
            e = np.array(sorted(by[alg]))
            v = [np.array(by[alg][x]) for x in e]
            med = np.array([np.median(x) for x in v])
            lo = np.array([x.min() for x in v])
            hi = np.array([x.max() for x in v])
            warm = alg == "GA from as built"
            ax.fill_between(e, lo, hi, color=ALG[alg], alpha=0.08 if warm else 0.15, lw=0)
            ax.plot(e, med, color=ALG[alg], lw=2.0 if warm else 2.2, ls=(0, (4, 2)) if warm else "-", label=alg)
            ends.append(med[-1]), names.append(alg)
        span = ax.get_ylim()[1] - ax.get_ylim()[0]
        for y, alg, yv in zip(spaced(ends, 0.035 * span), names, ends):
            ax.text(b * 1.015, y, f"{yv:.0f}", color=ALG[alg], fontsize=10.5, va="center")
        ax.axhline(v0, color=INK, lw=1, ls=(0, (5, 3)))
        ax.text(b * 0.02, v0 + 0.012 * span, f"as built, {v0:.1f}", fontsize=10, color=INK)
        ax.set_xlabel("evaluations")
        ax.set_ylabel("best value found, A$ thousand a week")
        ax.set_title(label, fontsize=13.5, loc="left")
        ax.set_xlim(0, b * 1.08)
        ax.legend(loc="upper left" if b == 20000 else "lower right", fontsize=10.5 if b == 20000 else 10, ncol=1)
    source(fig, "heuresis on four threads with keyed draws: GA on the mixed encoding; separable CMA-ES on the 621-dimensional relaxation; MAP-Elites over the evening and weekend shares;\n"
                "NSGA-II on value, spread and mix; random search; BO on a random 10-dimensional embedding (REMBO). Every search but the dashed GA starts from random plans.")
    save(fig, "search-63-searches-compared")


# ---- 64: the best plan on the map -----------------------------------------------------------------------------------

def map_plan(ax, T, G, plan, title_text, max_spend):
    kinds = [f["kind"] for f in plan["frontages"]]
    draw_base(ax, T, G, [KIND_COLOUR[k] for k in kinds])
    anchors = T([g["anchor"] for g in G["frontages"]])
    sp = np.array([f["spending_week"] for f in plan["frontages"]])
    ax.scatter(anchors[:, 0], anchors[:, 1], s=6 + 560 * sp / max_spend, c=[KIND_COLOUR[k] for k in kinds], edgecolors=INK, linewidths=0.7, zorder=6, alpha=0.95)
    for k, l in enumerate(G["links"]):
        if plan["plan"]["links"][k]:
            ab = T([l["a"], l["b"]])
            ax.plot(ab[:, 0], ab[:, 1], color=SAND, lw=3, ls=(0, (3, 1.6)), zorder=7)
    district_axes(ax, T, G, pad=(50, 25))
    ax.set_title(title_text, fontsize=13, loc="left")


def figure_64():
    G = load("district-geometry.json", PRESULTS)
    T = Frame(G["quad"])
    pl = plans_by_key()
    a, b = pl["as-built"], pl["best"]
    oa, ob = a["outcome"], b["outcome"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The best plan found without a model of the people",
          f"{b['from'].replace('20000', '20,000')}: A$ {ob['value_aud_week'] / 1000:.0f}k a week against {oa['value_aud_week'] / 1000:.1f}k as built; each frontage sized by the spending it captures")
    m = max(max(f["spending_week"] for f in a["frontages"]), max(f["spending_week"] for f in b["frontages"]))
    ax1 = fig.add_axes([0.01, 0.485, 0.58, 0.37])
    map_plan(ax1, T, G, a, f"as built: {sum(1 for f in a['frontages'] if f['kind'] != 'vacant')} frontages occupied", m)
    street_labels(ax1, T, G)
    ax2 = fig.add_axes([0.01, 0.075, 0.58, 0.37])
    map_plan(ax2, T, G, b, f"the best plan: {sum(1 for f in b['frontages'] if f['kind'] != 'vacant')} occupied, {ob['changed']:.0f} changed, {ob['links']:.0f} links made", m)
    street_labels(ax2, T, G)
    xl, yl = ax2.get_xlim(), ax2.get_ylim()
    north_arrow(ax2, T, (xl[1] - 25, yl[0] + 25))
    handles = [Patch(facecolor=KIND_COLOUR[k], edgecolor=PAPER, label=k) for k in KINDS]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.62, 0.87), ncol=3, fontsize=10.5, title="kind of ground floor", title_fontsize=11, alignment="left")
    # The value's parts.
    ax3 = fig.add_axes([0.66, 0.44, 0.31, 0.28])
    parts = [("margin", "margin_aud_week", 1), ("staffing", "staffing_aud_week", -1), ("outgoings", "outgoings_aud_week", -1), ("refits", "refits_aud_week", -1), ("links", "links_aud_week", -1)]
    x = np.arange(len(parts) + 1)
    for j, (o, col, off) in enumerate(((oa, MUTED, -0.18), (ob, "#8E3B78", 0.18))):
        vals = [s * o[k] / 1000 for _, k, s in parts] + [o["value_aud_week"] / 1000]
        ax3.bar(x + off, vals, width=0.34, color=col, alpha=0.85, label="as built" if j == 0 else "best")
    ax3.set_xticks(x, [p[0] for p in parts] + ["value"], fontsize=10.5)
    ax3.axhline(0, color=INK, lw=0.8)
    ax3.set_ylabel("A$ thousand a week")
    ax3.legend(fontsize=10.5, loc="upper right")
    ax3.set_title("the value's parts", fontsize=12.5, loc="left")
    # Frontages by kind.
    ax4 = fig.add_axes([0.66, 0.12, 0.31, 0.23])
    ka, kb = oa["frontages_by_kind"], ob["frontages_by_kind"]
    y = np.arange(len(KINDS))
    ax4.barh(y - 0.18, ka, height=0.34, color=[KIND_COLOUR[k] for k in KINDS], alpha=0.45, label="as built")
    ax4.barh(y + 0.18, kb, height=0.34, color=[KIND_COLOUR[k] for k in KINDS], label="best")
    for i in range(len(KINDS)):
        ax4.text(max(ka[i], kb[i]) + 0.6, i, f"{ka[i]:.0f} to {kb[i]:.0f}", va="center", fontsize=10, color=INK)
    ax4.set_yticks(y, KINDS, fontsize=10.5)
    ax4.invert_yaxis()
    ax4.set_xlabel("frontages (pale: as built)")
    ax4.grid(axis="y", visible=False)
    ax4.set_xlim(0, max(max(ka), max(kb)) + 8)
    source(fig, PROVENANCE + ". The value: the shared problem's, under study (1)'s flows. A study of methods on a real place, not a planning proposal.")
    save(fig, "search-64-best-plan-map")


# ---- 65: the week of the best plan ------------------------------------------------------------------------------------

def raster(ax, plan, P):
    n = len(plan["frontages"])
    img = np.ones((n, 168, 3)) * np.array(matplotlib.colors.to_rgb(PAPER))
    from matplotlib.colors import to_rgb
    for f, fr in enumerate(plan["frontages"]):
        k = fr["kind"]
        if k == "vacant":
            continue
        s1, l1, s2, l2 = fr["hours"]
        c = np.array(to_rgb(KIND_COLOUR[k]))
        for h in range(168):
            t = h % 24
            s, l = (s2, l2) if h >= 120 else (s1, l1)
            if l > 0 and s <= t < min(24, s + l):
                img[f, h] = c
            else:
                img[f, h] = 0.82 * np.array(to_rgb(PAPER)) + 0.18 * c
    ax.imshow(img, aspect="auto", interpolation="nearest", extent=(0, 168, n, 0))
    ax.grid(False)
    for d in range(1, 7):
        ax.axvline(24 * d, color=PAPER, lw=1.2)


import matplotlib  # noqa: E402


def figure_65():
    pl = plans_by_key()
    a, b = pl["as-built"], pl["best"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The week of the best plan: when each frontage opens, and what the district takes",
          "Each row a frontage, south to north along the area; coloured when open, pale when its kind is held but shut. Below, the district's spending and its people through the week")
    for j, (pln, box, lab) in enumerate(((a, [0.06, 0.53, 0.43, 0.33], "as built"), (b, [0.54, 0.53, 0.43, 0.33], "the best plan"))):
        ax = fig.add_axes(box)
        raster(ax, pln, None)
        week_axis(ax)
        ax.set_ylabel("frontage, south to north" if j == 0 else "")
        ax.set_title(lab, fontsize=13, loc="left")
    ax = fig.add_axes([0.06, 0.10, 0.87, 0.30])
    h = np.arange(168) + 0.5
    for pln, col, lab in ((a, MUTED, "as built"), (b, "#8E3B78", "the best plan")):
        s = np.array(pln["outcome"]["spending_by_hour"]) / 1000
        ax.fill_between(h, 0, s, step="mid", color=col, alpha=0.18 if pln is a else 0.25, lw=0)
        ax.step(h, s, where="mid", color=col, lw=1.6, label=f"spending, {lab}: A$ {sum(s):.0f}k a week")
    ax.set_ylabel("spending, A$ thousand an hour")
    week_axis(ax)
    ax2 = ax.twinx()
    ax2.grid(False)
    fa, fb = np.array(a["outcome"]["footfall_by_hour"]), np.array(b["outcome"]["footfall_by_hour"])
    ax2.plot(h, fa / 1000, color=SEA, lw=1.2, ls=(0, (4, 2)), label="people passing the frontages, as built")
    ax2.plot(h, fb / 1000, color=SEA, lw=1.6, label=f"… under the best plan ({fb.sum() / fa.sum() - 1:+.1%})")
    ax2.set_ylabel("people passing, thousands an hour", color=SEA)
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_color(HAIRLINE)
    l1, n1 = ax.get_legend_handles_labels()
    l2, n2 = ax2.get_legend_handles_labels()
    ax.legend(l1 + l2, n1 + n2, loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=10.5, ncol=4)
    ax.set_ylim(0, None)
    ax2.set_ylim(0, None)
    source(fig, PROVENANCE + ". Civil hours of the mean week, Monday 00:00 first.")
    save(fig, "search-65-week-of-the-best-plan")


# ---- 66: the week's lines of the spending -----------------------------------------------------------------------------

def figure_66():
    pl = plans_by_key()
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The spending in the frequency domain: the week's lines",
          "By Parseval the power of the spending's 168 hours is the sum of its lines: the mean (its spread), the day's harmonics, and the week's other lines (weekday against weekend)")
    ax = fig.add_axes([0.06, 0.15, 0.52, 0.69])
    keys = [("as-built", MUTED, "as built"), ("best", "#8E3B78", "the best value found"), ("spread", "#4E8F8B", "the most spread of NSGA-II's front worth as much as built")]
    for i, (k, col, lab) in enumerate(keys):
        p = np.array(pl[k]["outcome"]["power"])
        sh = p / p.sum()
        m = np.arange(len(sh))
        keep = m <= 49
        ax.vlines(m[keep] + (i - 1) * 0.27, 1e-7, np.maximum(sh[keep], 1e-7), color=col, lw=2.6, alpha=0.95, label=lab)
    for j in range(1, 8):
        ax.axvspan(7 * j - 0.5, 7 * j + 0.5, color=HAIRLINE, alpha=0.35, lw=0, zorder=0)
        ax.text(7 * j, 2.2, f"D{j}", ha="center", fontsize=10.5, color=INK)
    ax.text(0, 2.2, "mean", ha="center", fontsize=10.5, color=INK)
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1.6)
    ax.set_xlim(-1.2, 49.8)
    ax.set_xlabel("line of the week, m (cycles a week): the day's harmonics shaded, the week's other lines between them")
    ax.set_ylabel("share of the spending's power")
    ax.legend(loc="upper right", fontsize=11)
    # Shares.
    ax2 = fig.add_axes([0.66, 0.53, 0.31, 0.31])
    names = ["mean (spread)", "the day's lines", "the week's other lines"]
    for i, (k, col, lab) in enumerate(keys):
        o = pl[k]["outcome"]
        v = [o["spread"], o["daily_share"], o["weekly_share"]]
        ax2.barh(np.arange(3) + (i - 1) * 0.25, v, height=0.24, color=col, label=lab)
        for q in range(3):
            ax2.text(v[q] + 0.01, q + (i - 1) * 0.25, f"{v[q]:.2f}", va="center", fontsize=9.5, color=col)
    ax2.set_yticks(range(3), names)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 1)
    ax2.set_title("shares of the power", fontsize=12.5, loc="left")
    ax2.grid(axis="y", visible=False)
    # The mean days.
    ax3 = fig.add_axes([0.66, 0.15, 0.31, 0.27])
    for k, col, lab in keys:
        s = np.array(pl[k]["outcome"]["spending_by_hour"]).reshape(7, 24)
        wd, we = s[:5].mean(axis=0), s[5:].mean(axis=0)
        tot = s.sum()
        ax3.plot(np.arange(24) + 0.5, 1000 * wd / tot, color=col, lw=2)
        ax3.plot(np.arange(24) + 0.5, 1000 * we / tot, color=col, lw=1.4, ls=(0, (3, 2)))
    ax3.set_xlim(0, 24)
    ax3.set_xticks([0, 6, 12, 18, 24])
    ax3.set_xlabel("civil hour: weekday (full) and weekend day (dashed)")
    ax3.set_ylabel("‰ of the week's spending")
    ax3.set_title("the mean days", fontsize=12.5, loc="left")
    o0, o1, o2 = pl["as-built"]["outcome"], pl["best"]["outcome"], pl["spread"]["outcome"]
    fig.text(0.06, 0.055, f"Effective hours (168 × spread): as built {o0['effective_hours']:.0f}, best {o1['effective_hours']:.0f}, most spread {o2['effective_hours']:.0f}; "
                          f"evening share {o0['evening_share']:.2f}, {o1['evening_share']:.2f}, {o2['evening_share']:.2f}; weekend share {o0['weekend_share']:.2f}, {o1['weekend_share']:.2f}, {o2['weekend_share']:.2f}.",
             fontsize=11, color=INK)
    source(fig, PROVENANCE + ". The lines are the exact discrete Fourier transform of the mean week's 168 hours (economy.hpp, week_power).")
    save(fig, "search-66-spending-lines")


# ---- 67: a repertoire of evenings and weekends ---------------------------------------------------------------------

def figure_67():
    rep = load("repertoire.json")
    scores = read_csv("repertoire-scores.csv")
    summ = load("summary.json")
    v0 = summ["as_built"]["value_aud_week"]
    n = int(rep["grid"]["resolution"])
    hi = rep["grid"]["high"]
    M = np.full((n, n), np.nan)
    for c in rep["cells"]:
        M[int(c["weekend_bin"]), int(c["evening_bin"])] = c["value"] / 1000
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "A repertoire of evenings and weekends",
          "MAP-Elites keeps the best plan for each share of the spending taken in the evening (18–24 h) and at the weekend: how much value each balance of the week can keep")
    ax = fig.add_axes([0.07, 0.12, 0.48, 0.72])
    vmax = np.nanmax(M)
    norm = TwoSlopeNorm(vmin=min(np.nanmin(M), -1), vcenter=0, vmax=vmax)
    im = ax.imshow(M, origin="lower", extent=(0, hi, 0, hi), cmap=DUSK.reversed(), norm=norm, aspect="equal")
    ax.grid(False)
    o = summ["as_built"]
    ax.scatter([o["evening_share"]], [o["weekend_share"]], marker="*", s=260, color=PAPER, edgecolors=INK, linewidths=1.2, zorder=5)
    box = dict(boxstyle="round,pad=0.25", facecolor=PAPER, edgecolor="none", alpha=0.9)
    ax.annotate(f"as built, {v0 / 1000:.1f}k", (o["evening_share"], o["weekend_share"]), xytext=(-12, -22), textcoords="offset points", fontsize=10.5, color=INK, ha="right", bbox=box)
    b = summ["best"]
    ax.scatter([b["evening_share"]], [b["weekend_share"]], marker="D", s=70, color="#8E3B78", edgecolors=PAPER, linewidths=1.0, zorder=5)
    ax.annotate(f"the best value found, {b['value_aud_week'] / 1000:.0f}k", (b["evening_share"], b["weekend_share"]), xytext=(8, -26), textcoords="offset points", fontsize=10.5, color=INK, bbox=box)
    ax.axhline(2 / 7, color=INK, lw=0.6, ls=(0, (2, 3)))
    ax.text(hi * 0.99, 2 / 7 + 0.008, "2/7, the weekend's share of the hours", ha="right", fontsize=9.5, color=INK, bbox=box)
    ax.set_xlabel("evening share of the spending (18–24 h)")
    ax.set_ylabel("weekend share of the spending")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("best value in the cell, A$ thousand a week")
    cb.outline.set_visible(False)
    # QD of all that each algorithm evaluated, binned in the same grid.
    ax2 = fig.add_axes([0.66, 0.50, 0.31, 0.34])
    ax3 = fig.add_axes([0.66, 0.12, 0.31, 0.28])
    for i, alg in enumerate([a for a in LARGE if a in {r["algorithm"] for r in scores}]):
        q = [float(r["qd_score_kaud"]) for r in scores if r["algorithm"] == alg and r["budget"] == "20000"]
        c = [float(r["coverage"]) for r in scores if r["algorithm"] == alg and r["budget"] == "20000"]
        ax2.scatter([i] * len(q), q, color=ALG[alg], s=36, alpha=0.85)
        ax2.plot([i - 0.25, i + 0.25], [np.median(q)] * 2, color=ALG[alg], lw=2.4)
        ax3.scatter([i] * len(c), c, color=ALG[alg], s=36, alpha=0.85)
        ax3.plot([i - 0.25, i + 0.25], [np.median(c)] * 2, color=ALG[alg], lw=2.4)
    present = [a for a in LARGE if a in {r["algorithm"] for r in scores}]
    for a_ in (ax2, ax3):
        a_.set_xticks(range(len(present)), [short_name(a) for a in present], fontsize=10)
        a_.grid(axis="x", visible=False)
    ax2.set_ylabel("QD-score: Σ over cells of max(0, best value), A$k")
    ax2.set_title("every evaluation of each run, in the same grid (20,000; five seeds)", fontsize=12, loc="left")
    ax3.set_ylabel("cells reached")
    ax3.set_ylim(0, 1)
    source(fig, "MAP-Elites (heuresis diversity) on a 12 × 12 grid of [0, 0.6]², seed 1, 20,000 evaluations, half bred by crossover. " + PROVENANCE + ".")
    save(fig, "search-67-repertoire-evenings-weekends")


# ---- 68: value, spread and mix -------------------------------------------------------------------------------------------

def figure_68():
    fr = read_csv("fronts.csv")
    hvs = read_csv("hypervolume.csv")
    summ = load("summary.json")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Value against spread and mix: the fronts of all that was evaluated",
          "The non-dominated plans in value, spread over the week and mix of kinds among every evaluation of each run (20,000, five seeds); their colour the mix")
    ax = fig.add_axes([0.07, 0.12, 0.52, 0.72])
    rnd = [r for r in fr if r["algorithm"] == "random"]
    ax.scatter([float(r["spread"]) for r in rnd], [float(r["value"]) / 1000 for r in rnd], s=10, color=ALG["random"], alpha=0.5, label="random search's front")
    sc = None
    for alg, marker, size in (("NSGA-II", "o", 18), ("GA", "s", 16), ("CMA-ES", "^", 22), ("GA from as built", "D", 16)):
        ns = [r for r in fr if r["algorithm"] == alg]
        sc = ax.scatter([float(r["spread"]) for r in ns], [float(r["value"]) / 1000 for r in ns], s=size, marker=marker, c=[float(r["mix"]) for r in ns], cmap=ORCHID, vmin=0, vmax=1.6,
                        edgecolors="none", label=f"{alg}'s fronts")
    o = summ["as_built"]
    ax.scatter([o["spread"]], [o["value_aud_week"] / 1000], marker="*", s=300, color=PAPER, edgecolors=INK, linewidths=1.2, zorder=5, label=f"as built (mix {o['mix_nats']:.2f})")
    ax.set_xlabel("spread of the spending over the week (1/(1 + CV²))")
    ax.set_ylabel("value, A$ thousand a week")
    ax.legend(loc="lower left", fontsize=10.5)
    cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("mix H(K | block), nats")
    cb.outline.set_visible(False)
    ax2 = fig.add_axes([0.70, 0.40, 0.27, 0.44])
    present = [a for a in LARGE if a in {r["algorithm"] for r in hvs}]
    for i, alg in enumerate(present):
        v = [float(r["hypervolume"]) for r in hvs if r["algorithm"] == alg]
        ax2.scatter([i] * len(v), v, color=ALG[alg], s=40, alpha=0.85)
        ax2.plot([i - 0.25, i + 0.25], [np.median(v)] * 2, color=ALG[alg], lw=2.4)
        ax2.text(i, max(v) + 0.01, f"{np.median(v):.3f}", ha="center", fontsize=10, color=ALG[alg])
    ax2.set_xticks(range(len(present)), [short_name(a) for a in present], fontsize=10)
    ax2.grid(axis="x", visible=False)
    ax2.set_ylabel("hypervolume, normalised (reference 1.1)")
    ax2.set_title("hypervolume of all each run evaluated", fontsize=12.5, loc="left")
    fig.text(0.70, 0.30, "The value's axis runs from the best value found to the plan\nas built: a plan worse than as built adds nothing.\n"
                         "Spread and mix run from the union's least to its greatest.", fontsize=10.5, color=MUTED, va="top")
    source(fig, "heuresis multiobjective: exact hypervolume (WFG). " + PROVENANCE + ".")
    save(fig, "search-68-value-spread-mix")


# ---- 69: where the model-free plan is blind -----------------------------------------------------------------------------

def figure_69():
    G, P = load("district-geometry.json", PRESULTS), load("problem.json", PRESULTS)
    bl = load("blindness.json")
    T = Frame(G["quad"])
    links = bl["links"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    fp = bl["footfall"]
    title(fig, "Where a plan without a model is blind",
          "In the data alone a link changes nobody's way: the value sees its cost and the catalysis of frontages it brings nearer. On the network, people's shortest walks would take it")
    ax = fig.add_axes([0.01, 0.30, 0.58, 0.58])
    best = max(links, key=lambda l: l["share_of_walks_through_it"])
    change = np.array(best["own_segment_change"])
    norm = TwoSlopeNorm(vmin=-0.6, vcenter=0, vmax=0.6)
    cmap = DUSK.reversed()
    draw_base(ax, T, G, [cmap(norm(c)) for c in change], edge=PAPER)
    for k, l in enumerate(G["links"]):
        ab = T([l["a"], l["b"]])
        sh = links[k]["share_of_walks_through_it"]
        ax.plot(ab[:, 0], ab[:, 1], color=SAND, lw=1.2 + 40 * sh, solid_capstyle="round", zorder=7, alpha=0.9)
        mid = ab.mean(axis=0)
        ax.text(mid[0], mid[1] + 10, f"{100 * sh:.1f}%", fontsize=9.5, color=INK, ha="center", zorder=8,
                bbox=dict(boxstyle="round,pad=0.2", facecolor=PAPER, edgecolor="none", alpha=0.85))
    district_axes(ax, T, G)
    street_labels(ax, T, G)
    xl, yl = ax.get_xlim(), ax.get_ylim()
    north_arrow(ax, T, (xl[0] + 25, yl[0] + 40))
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cax = fig.add_axes([0.06, 0.25, 0.25, 0.015])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label(f"change in the walks past each frontage if {' to '.join(best['joins'])} were made", fontsize=10)
    cb.outline.set_visible(False)
    fig.text(0.06, 0.19, "Links in sand, their width and label the share of the people's shortest walks between frontages (pairs weighted by\n"
                         "their passers) that would pass through them. The value of study (1) sees none of it.", fontsize=10.5, color=MUTED, va="top")
    # Bars: what the network says, against what the value says.
    ax2 = fig.add_axes([0.71, 0.52, 0.26, 0.34])
    y = np.arange(len(links))
    names = [f"{l['joins'][0]} – {l['joins'][1]}" for l in links]
    ax2.barh(y, [100 * l["share_of_walks_through_it"] for l in links], color=SAND, height=0.6)
    ax2.set_yticks(y, names, fontsize=9.5)
    ax2.invert_yaxis()
    ax2.set_xlabel("% of the walks between frontages through it")
    ax2.set_title("what the network says", fontsize=12.5, loc="left")
    ax2.grid(axis="y", visible=False)
    ax3 = fig.add_axes([0.71, 0.11, 0.26, 0.30])
    ax3.barh(y, [l["value_change_aud_week"] / 1000 for l in links], color=MUTED, height=0.6)
    ax3.set_yticks(y, names, fontsize=9.5)
    ax3.invert_yaxis()
    ax3.set_xlabel("change in value with the link made, A$ thousand a week")
    ax3.set_title("what the value without a model says", fontsize=12.5, loc="left")
    ax3.grid(axis="y", visible=False)
    fig.text(0.06, 0.135, f"Nor is anybody conserved: under the best plan {100 * fp['change']:+.1f}% more people pass the frontages than as built,\n"
                          "drawn from nowhere by the empirical rule, and taken from no other street.", fontsize=11, color=INK, va="top")
    source(fig, "Shortest walks by Dijkstra (the designing study's walking.hpp) on the urban-science domain's district with each link made alone. " + PROVENANCE + ".")
    save(fig, "search-69-where-it-is-blind")


FIGURES = {62: figure_62, 63: figure_63, 64: figure_64, 65: figure_65, 66: figure_66, 67: figure_67, 68: figure_68, 69: figure_69}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    which = [int(x) for x in args.only.split(",") if x] or list(FIGURES)
    for k in which:
        FIGURES[k]()
        print("drawn", k)
