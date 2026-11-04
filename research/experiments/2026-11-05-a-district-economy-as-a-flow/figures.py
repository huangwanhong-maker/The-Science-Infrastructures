"""The figures of study (2) of "A district's economy": people as a flow on the district's streets, drawn from results/
(and the shared problem's uncommitted geometry) in the house style of the search case studies, the place in the Urban
Science Workbench's sea blue and sand.

    .venv/Scripts/python.exe research/experiments/2026-11-05-a-district-economy-as-a-flow/figures.py [--only 70,71] [--results DIR]

Writes assets/gallery/search-70-... to search-79-... (1600 x 1000).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F401,F403

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Patch  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
PROBLEM = HERE.parent / "2026-11-05-a-district-economy-the-problem" / "results"

SEA = "#5B8FB9"
SAND = "#C98B6B"
DEEP = "#1F3B57"
CONTEXT = "#EEE7DF"
PARCEL = "#E2D8CE"
NETWORK = "#CDBFB6"
KIND_COLOUR = {0: "#E9E2DA", 1: "#5B8FB9", 2: "#C98B6B", 3: "#7FA37A", 4: "#9A9097", 5: "#8E3B78"}
KIND_NAME = {0: "vacant", 1: "retail", 2: "food and drink", 3: "services", 4: "office lobby", 5: "culture"}
FLOW_MAP = LinearSegmentedColormap.from_list("people", ["#E8EEF3", "#A9C3DA", "#5B8FB9", "#2E5E86", "#1F3B57"])
CREST_MAP = LinearSegmentedColormap.from_list("crest", ["#2F3E5C", "#5B8FB9", "#A9C3DA", "#E9B97F", "#C98B6B", "#8E3B78", "#2F3E5C"])
CREST_SEQ = LinearSegmentedColormap.from_list("crest_seq", ["#2E5E86", "#5B8FB9", "#9DBFD6", "#E9C9A3", "#C98B6B", "#8E3B78", "#4A2140"])
CIRC_MAP = LinearSegmentedColormap.from_list("circ", ["#2E5E86", "#7FA6C7", "#F4EEE8", "#E3A985", "#A85A3A"])
ALGO = {
    "GA": ALGORITHM["GA"], "CMA-ES": ALGORITHM["CMA-ES"], "random": ALGORITHM["random"],
    "GA, marginal kernels": "#D49A3A", "BR": SEA, "BR + (1+1) EA": "#1F3B57", "GA from as built": "#4E8F8B",
}
ORDER = ["random", "CMA-ES", "GA", "GA from as built", "GA, marginal kernels", "BR", "BR + (1+1) EA"]
CREDIT = "Data: OpenStreetMap contributors (ODbL 1.0); City of Melbourne pedestrian counts and CLUE (CC BY 4.0)"


def load(path):
    return json.loads(pathlib.Path(path).read_text(encoding="utf-8"))


def read_csv(path):
    with pathlib.Path(path).open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


class Data:
    def __init__(self, results: pathlib.Path):
        self.cal = load(results / "calibration.json")
        self.plans = load(results / "plans.json")
        self.maps = load(results / "flow-maps.json")
        self.summary = load(results / "summary.json")
        self.runs = read_csv(results / "runs.csv")
        self.anytime = read_csv(results / "anytime.csv")
        if (results / "runs-extra.csv").exists():
            self.runs += read_csv(results / "runs-extra.csv")
            self.anytime += read_csv(results / "anytime-extra.csv")
        self.geo = load(PROBLEM / "district-geometry.json")
        self.problem = load(PROBLEM / "problem.json")
        self.edges = self.maps["edges"]
        self.lines = [np.array(e["line"]) for e in self.edges]
        self.nodes = np.array(self.maps["nodes"])
        self.plan = {p["plan"]: p for p in self.plans["plans"]}
        self.map = {p["plan"]: p for p in self.maps["plans"]}

    def present(self, plan_name):
        links = self.map[plan_name]["links"]
        return [e["link"] < 0 or links[int(e["link"])] == 1 for e in self.edges]


# ---- drawing helpers ------------------------------------------------------------------------------------------

def frame(ax, D, pad=12):
    xs = np.concatenate([l[:, 0] for l in D.lines])
    ys = np.concatenate([l[:, 1] for l in D.lines])
    ax.set_xlim(xs.min() - pad, xs.max() + pad)
    ax.set_ylim(ys.min() - pad, ys.max() + pad)
    ax.set_aspect("equal")
    ax.axis("off")


def buildings(ax, D, parcels=True, alpha=1.0):
    ctx = [np.array(c["footprint"]) for c in D.geo["context"]]
    ax.add_collection(PolyCollection(ctx, facecolors=CONTEXT, edgecolors="none", alpha=alpha, zorder=1))
    if parcels:
        ps = [np.array(f["footprint"]) for f in D.geo["frontages"]]
        ax.add_collection(PolyCollection(ps, facecolors=PARCEL, edgecolors="#D3C6BB", linewidths=0.4, alpha=alpha, zorder=1))


def network(ax, D, present=None, colour=NETWORK, lw=0.6, zorder=2):
    segs = [l for i, l in enumerate(D.lines) if (present is None or present[i]) and D.edges[i]["link"] < 0]
    ax.add_collection(LineCollection(segs, colors=colour, linewidths=lw, zorder=zorder, capstyle="round"))


def links(ax, D, made=None, all_candidates=True, lw=2.2):
    for k, l in enumerate(D.geo["links"]):
        a, b = l["a"], l["b"]
        on = made is not None and made[k] == 1
        if not (on or all_candidates):
            continue
        ax.plot([a[0], b[0]], [a[1], b[1]], color=SAND if on else "#6B3E75", lw=lw if on else 1.6, ls=(0, (3, 1.5)) if l["kind"] == "lane carried through" else (0, (1, 1.2)),
                zorder=6, solid_capstyle="round", dash_capstyle="round")


def flow_lines(ax, D, values, present, vmax, width_max=9.0, cmap=FLOW_MAP, zorder=3, power=0.7):
    segs, ws, cs = [], [], []
    for i, l in enumerate(D.lines):
        if not present[i]:
            continue
        v = max(0.0, values[i])
        segs.append(l)
        ws.append(0.35 + width_max * (min(v, vmax) / vmax) ** power)
        cs.append(cmap(min(1.0, (v / vmax) ** power)))
    ax.add_collection(LineCollection(segs, colors=cs, linewidths=ws, zorder=zorder, capstyle="round", joinstyle="round"))


def midpoint_direction(line):
    d = np.r_[0, np.cumsum(np.hypot(*np.diff(line, axis=0).T))]
    L = d[-1]
    s = L / 2
    k = max(0, min(len(line) - 2, np.searchsorted(d, s) - 1))
    t = (s - d[k]) / max(1e-9, d[k + 1] - d[k])
    p = line[k] + t * (line[k + 1] - line[k])
    u = line[k + 1] - line[k]
    u = u / max(1e-9, np.hypot(*u))
    return p, u, L


def net_arrows(ax, D, net, present, scale, colour=DEEP, min_value=None, zorder=7, max_len=26):
    """An arrow at each segment's middle along its net flow (u → v positive), its length and width with |net|."""
    vmax = scale
    for i, l in enumerate(D.lines):
        if not present[i]:
            continue
        v = net[i]
        if min_value is not None and abs(v) < min_value:
            continue
        p, u, L = midpoint_direction(l)
        if L < 12:
            continue
        a = min(1.0, abs(v) / vmax)
        length = min(L * 0.8, 6 + max_len * a ** 0.6)
        d = u * (1 if v > 0 else -1)
        start, end = p - d * length / 2, p + d * length / 2
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=4 + 9 * a ** 0.5, lw=0.5 + 1.6 * a ** 0.6, color=colour, zorder=zorder, alpha=0.95))


def face_rings(D, faces):
    rings = []
    for f in faces:
        ring = []
        for e, s in f["edges"]:
            l = D.lines[int(e)]
            l = l if s > 0 else l[::-1]
            ring.extend(l[:-1].tolist())
        rings.append(np.array(ring))
    return rings


def legend_lines(ax, entries, loc="lower left", title=None, fontsize=10, **kw):
    handles = [Line2D([0], [0], color=c, lw=w, ls=ls) for (c, w, ls, _) in entries]
    ax.legend(handles, [t for (*_, t) in entries], loc=loc, frameon=True, framealpha=0.92, facecolor=PAPER, edgecolor=HAIRLINE, fontsize=fontsize, title=title, title_fontsize=fontsize, **kw)


def colourbar(fig, rect, cmap, vmin, vmax, label, ticks=None, ticklabels=None):
    cax = fig.add_axes(rect)
    sm = plt.cm.ScalarMappable(norm=Normalize(vmin, vmax), cmap=cmap)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.outline.set_edgecolor(HAIRLINE)
    cb.ax.tick_params(labelsize=10, colors=MUTED)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabels is not None:
        cb.set_ticklabels(ticklabels)
    cb.set_label(label, fontsize=11, color=INK)
    return cb


def fmt_aud(v):
    return f"A${v / 1000:,.0f}k" if abs(v) >= 1000 else f"A${v:,.0f}"


def day_axis(ax, ymax=None):
    for d in range(1, 7):
        ax.axvline(24 * d, color=HAIRLINE, lw=0.8, zorder=0)
    ax.set_xticks([12 + 24 * d for d in range(7)])
    ax.set_xticklabels(["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"])
    ax.set_xlim(0, 168)
    ax.grid(axis="x", visible=False)


# ---- 70: the district as a network ---------------------------------------------------------------------------

def fig70(D):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The district as a network the people flow on",
          "Segments drawn by walkable width; floor on the junctions; CLUE's establishments; the City's counters; the entrances and their counters")
    ax = fig.add_axes([0.02, 0.06, 0.6, 0.83])
    buildings(ax, D)
    widths = [e["width_m"] for e in D.edges]
    segs = [l for i, l in enumerate(D.lines) if D.edges[i]["link"] < 0]
    ws = [0.4 + 0.55 * w for i, w in enumerate(widths) if D.edges[i]["link"] < 0]
    cs = [SEA if D.edges[i]["kind"] != "laneway" else SAND for i in range(len(D.edges)) if D.edges[i]["link"] < 0]
    ax.add_collection(LineCollection(segs, colors=cs, linewidths=ws, zorder=3, capstyle="round", alpha=0.85))
    links(ax, D, made=None)
    nodes = D.nodes
    floor = nodes[:, 2] + nodes[:, 3]
    ok = floor > 0
    ax.scatter(nodes[ok, 0], nodes[ok, 1], s=3 + floor[ok] / 900, color=DEEP, alpha=0.4, lw=0, zorder=4)
    att = D.maps["attractors"]
    ax.scatter([nodes[a["node"], 0] for a in att], [nodes[a["node"], 1] for a in att], s=[3 + 2.2 * sum(a["count"]) for a in att], facecolors="none",
               edgecolors="#8E3B78", lw=0.7, alpha=0.8, zorder=5)
    cal = {c["id"]: c for c in D.cal["counters"]}
    for c in D.maps["counters"]:
        x, y = c["at"]
        ax.scatter([x], [y], marker="^", s=70, color=SAND, edgecolor=INK, lw=0.6, zorder=8)
    for e in D.maps["entrances"]:
        x, y = nodes[e["node"], 0], nodes[e["node"], 1]
        week = sum(e["gross"])
        ax.scatter([x], [y], marker="o", s=18 + week / 6000, color="white", edgecolor=SEA if e["source"] != "default" else MUTED, lw=1.4, zorder=8)
    frame(ax, D)
    handles = [Line2D([0], [0], color=SEA, lw=6.6, label="street, 9 m walkable (Swanston 12 m)"), Line2D([0], [0], color=SEA, lw=2.6, label="little street, 4 m"),
               Line2D([0], [0], color=SAND, lw=2.3, label="laneway, 3.5 m"), Line2D([0], [0], color="#6B3E75", lw=1.6, ls=(0, (3, 1.5)), label="candidate link: lane carried through"),
               Line2D([0], [0], color="#6B3E75", lw=1.6, ls=(0, (1, 1.2)), label="candidate link: passage"),
               Line2D([0], [0], marker="o", ls="none", color=DEEP, mec="none", alpha=0.45, ms=8, label="floor on a junction"),
               Line2D([0], [0], marker="o", ls="none", mfc="none", mec="#8E3B78", ms=7, label="CLUE establishments (trading)"),
               Line2D([0], [0], marker="^", ls="none", color=SAND, mec=INK, ms=8, label="the City's counter"),
               Line2D([0], [0], marker="o", ls="none", mfc="white", mec=SEA, mew=1.4, ms=8, label="entrance (size: its weekly flow)")]
    ax.legend(handles=handles, loc="lower left", fontsize=9.2, frameon=True, framealpha=0.94, facecolor=PAPER, edgecolor=HAIRLINE, numpoints=1)
    # The right side: what the network holds, and the route choice's calibration.
    t = fig.add_axes([0.64, 0.52, 0.33, 0.37])
    t.axis("off")
    s = D.summary
    att_tot = D.cal["attractors"]
    lines = [
        ("Network", f"{len(D.edges)} segments, {len(D.nodes)} junctions (with every link's)"),
        ("Entrances", f"{len(D.maps['entrances'])} dead ends on the box's edge, fed by counters outside it"),
        ("Floor", f"{D.cal['floor_m2']['work'] / 1e6:.2f} M m² of work, {D.cal['floor_m2']['home'] / 1e6:.2f} M m² of homes"),
        ("Attractors", f"{att_tot['retail']:.0f} retail, {att_tot['food and drink']:.0f} food, {att_tot['services']:.0f} services, {att_tot['culture']:.0f} culture"),
        ("Counters", f"{len(D.cal['counters'])} inside the box, {sum(1 for c in D.cal['counters'] if c.get('directed'))} read each way"),
        ("Route choice", f"θ = {D.cal['theta_per_m']:.3f} per metre; spectral radius {D.cal['radius']['all_links']:.2f}"),
    ]
    for k, (a, b) in enumerate(lines):
        t.text(0.0, 1 - k * 0.15, a, fontsize=13, color=INK, weight="bold", va="top")
        t.text(0.0, 1 - k * 0.15 - 0.065, b, fontsize=11.5, color=MUTED, va="top")
    g = fig.add_axes([0.68, 0.12, 0.27, 0.3])
    grid = [x for x in D.cal["theta_grid"] if not x.get("diverges")]
    th = [x["theta_per_m"] for x in grid]
    g.plot(th, [x["r2"] for x in grid], "-o", color=SEA, ms=5, label="gross, fitted")
    g.plot(th, [x["loo_r2"] for x in grid], "-o", color=SAND, ms=5, label="gross, each counter left out")
    g.axvline(D.cal["theta_per_m"], color=MUTED, lw=0.8, ls=":")
    g.set_xscale("log")
    g.set_xticks(th)
    g.set_xticklabels([f"{t:g}" for t in th], fontsize=9.5)
    g.xaxis.set_minor_formatter(plt.NullFormatter())
    g.set_xlabel("θ, the route choice's sensitivity to length (per m)")
    g.set_ylabel("R² of the counters' mean weeks")
    g.legend(loc="lower right", fontsize=9.5)
    g.set_title("Calibrating the walkers' choice of way", fontsize=12, loc="left")
    source(fig, CREDIT + ". Walkable widths by Hoddle's classes of street (assumptions).")
    return fig


# ---- 71: calibration ----------------------------------------------------------------------------------------------

def fig71(D):
    cs = sorted(D.cal["counters"], key=lambda c: -sum(c["observed"]))
    n = len(cs)
    cols = 4
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(16, 10), dpi=100, sharex=True)
    fig.subplots_adjust(left=0.05, right=0.985, top=0.86, bottom=0.08, hspace=0.42, wspace=0.18)
    title(fig, "The counters' mean weeks, and the flow fitted to them",
          f"Hour by hour by non-negative least squares on each way; gross R² {D.cal['r2']:.2f} fitted, {D.cal['loo_r2']:.2f} with each counter left out (θ = {D.cal['theta_per_m']:.3f}/m)")
    h = np.arange(168) + 0.5
    for k, ax in enumerate(axes.flat):
        if k >= n:
            ax.axis("off")
            continue
        c = cs[k]
        obs, fit, loo = np.array(c["observed"]), np.array(c["fitted"]), np.array(c["loo"])
        ax.fill_between(h, 0, obs, color=SEA, alpha=0.22, lw=0, step="mid")
        ax.plot(h, obs, color=SEA, lw=1.0, drawstyle="steps-mid")
        ax.plot(h, fit, color=SAND, lw=1.4)
        ax.plot(h, loo, color=INK, lw=0.8, ls=(0, (2, 1.5)), alpha=0.8)
        day_axis(ax)
        ax.tick_params(labelsize=8.5)
        name = c["name"].replace("(", "").replace(")", "")
        ax.set_title(f"{name[:30]}", fontsize=10.5, loc="left")
        ax.text(0.99, 0.95, f"R² {c['r2']:.2f} · out {c['loo_r2']:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=MUTED)
        ax.set_ylim(0, max(obs.max(), fit.max()) * 1.18)
    handles = [Patch(color=SEA, alpha=0.4, label="observed (both ways, people an hour)"), Line2D([0], [0], color=SAND, lw=1.6, label="the model, fitted"),
               Line2D([0], [0], color=INK, lw=0.9, ls=(0, (2, 1.5)), label="the model with this counter left out")]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.985, 0.995), ncol=1, fontsize=10.5)
    source(fig, "City of Melbourne pedestrian counts (CC BY 4.0), 52 weeks from 27 October 2025, public holidays left out; each way read by the urban-science domain's current meter. Sorted by weekly count.")
    return fig


# ---- 72: people on the streets ---------------------------------------------------------------------------------

def fig72(D):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "People on the streets, as built", "The walkers' gross flow on every segment, both ways, at four hours of the mean week, and what the streams of trips are")
    p = D.map["as built"]
    present = D.present("as built")
    panels = [("gross_wed_08", "Wednesday 08:00"), ("gross_wed_12", "Wednesday 12:00"), ("gross_wed_17", "Wednesday 17:00"), ("gross_sat_14", "Saturday 14:00")]
    vmax = np.percentile([v for k, _ in panels for v, ok in zip(p[k], present) if ok], 99.5)
    for i, (key, label) in enumerate(panels):
        r, c = divmod(i, 2)
        ax = fig.add_axes([0.015 + c * 0.31, 0.49 - r * 0.41, 0.3, 0.37])
        buildings(ax, D, alpha=0.8)
        flow_lines(ax, D, p[key], present, vmax, width_max=7)
        frame(ax, D, pad=4)
        tot = sum(v for v, ok in zip(p[key], present) if ok)
        ax.set_title(label, fontsize=13.5, loc="left", color=INK)
        ax.text(0.99, 1.005, f"{tot / 1000:,.0f}k walker-segments an hour", transform=ax.transAxes, fontsize=9.5, color=MUTED, ha="right", va="bottom")
    colourbar(fig, [0.67, 0.13, 0.29, 0.016], FLOW_MAP, 0, vmax, "people an hour on the segment, both ways", ticks=np.linspace(0, vmax, 5),
              ticklabels=[f"{v:,.0f}" for v in np.linspace(0, vmax, 5)])
    # The streams: each mode's share of the week's walker-segments.
    a = fig.add_axes([0.73, 0.52, 0.24, 0.31])
    mg = p["mode_week_gross"]
    names = [n for n in mg if sum(mg[n]) > 0]
    tot = sum(sum(mg[n]) for n in names)
    order = sorted(names, key=lambda n: -sum(mg[n]))
    cols = {"through": SEA, "local": SAND, "arrive (work)": "#7FA37A", "leave (work)": "#4E8F8B", "visit": "#8E3B78", "arrive (home)": "#D49A3A", "leave (home)": "#9A9097"}
    a.barh(range(len(order)), [sum(mg[n]) / tot for n in order], color=[cols.get(n, MUTED) for n in order], height=0.6)
    for k, n in enumerate(order):
        a.text(sum(mg[n]) / tot + 0.01, k, f"{sum(mg[n]) / tot:.0%}", va="center", fontsize=10.5, color=INK)
    a.set_yticks(range(len(order)))
    a.set_yticklabels(order, fontsize=11)
    a.invert_yaxis()
    a.set_xlim(0, 1)
    a.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
    a.set_title("What the walkers are doing: each stream's share\nof the week's walking, as fitted", fontsize=12, loc="left")
    a.grid(axis="y", visible=False)
    t = fig.add_axes([0.66, 0.22, 0.32, 0.22])
    t.axis("off")
    txt = ("through: from one entrance of the box to another\n"
           "local: between the buildings, by floor and distance\n"
           "arrive and leave: between the entrances and the floor\n"
           "visit: to the open shops, cafés, services and venues\n    and back, set at a tenth of the hour's trips\n\n"
           "Each stream is conserved and chooses its way by a logit\nover every walk, wider streets preferred as a fluid's\n"
           "flux prefers a wider channel. The counters fix the\namplitudes hour by hour; the patterns are the model's.")
    t.text(0, 1, txt, fontsize=11, color=INK, va="top", linespacing=1.45)
    source(fig, CREDIT + ".")
    return fig


# ---- 73: circulation -----------------------------------------------------------------------------------------------

def fig73(D, plan="as built", label="Wednesday 08:00"):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    m = D.map[plan]["moments"][label]
    sh = m["shares"]
    key = {"Wednesday 08:00": "gross_wed_08", "Wednesday 12:00": "gross_wed_12", "Wednesday 17:00": "gross_wed_17", "Saturday 14:00": "gross_sat_14"}[label]
    present = D.present(plan)
    gross = sum(v for v, ok in zip(D.map[plan][key], present) if ok)
    net = sum(abs(v) for v, ok in zip(m["net"], present) if ok)
    title(fig, "Where the morning's people are driven, and where they eddy",
          f"{label}, {plan}: the net flow (Σ|net| is {net / gross:.0%} of Σ gross) by the economics domain's circulation; "
          f"energy {sh[0]:.0%} driven, {sh[1]:.0%} round blocks, {sh[2]:.0%} round holes")
    parts = [("driven", "driven: up the attraction, as a Darcy fluid would go"), ("local", "round the blocks: the eddies"), ("global", "round the blocks the domain leaves as holes")]
    big = max(abs(v) for v in m["net"])
    for i, (key2, lab) in enumerate(parts):
        ax = fig.add_axes([0.01 + i * 0.33, 0.3, 0.32, 0.58])
        buildings(ax, D, parcels=True, alpha=0.6)
        network(ax, D, present, lw=0.5)
        if key2 == "local":
            rings = face_rings(D, m["faces"])
            circ = np.array([f["circulation"] for f in m["faces"]])
            cm = max(1e-9, np.abs(circ).max())
            ax.add_collection(PolyCollection(rings, facecolors=[CIRC_MAP(0.5 + 0.5 * c / cm) for c in circ], edgecolors="none", alpha=0.8, zorder=2))
        vals = np.array(m[key2])
        pmax = max(1e-9, np.abs(vals).max())
        net_arrows(ax, D, vals, present, scale=pmax * 0.5, colour=DEEP if key2 != "local" else "#7A3A22", min_value=pmax * 0.04)
        frame(ax, D, pad=4)
        ax.set_title(lab, fontsize=12.5, loc="left")
        ax.text(0.99, -0.01, f"largest {pmax:,.0f} people an hour" + ("" if key2 == "driven" else f" (arrows × {big / pmax:,.0f} beside the driven)"),
                transform=ax.transAxes, ha="right", va="top", fontsize=9.5, color=MUTED)
    ax = fig.add_axes([0.06, 0.07, 0.88, 0.16])
    shh = np.array(D.plan[plan]["circulation_by_hour"])
    h = np.arange(168) + 0.5
    E = shh[:, 3]
    for k, (c, lab2) in enumerate([(SEA, "driven"), (SAND, "round the blocks"), ("#8E3B78", "round the holes")]):
        ax.plot(h, np.maximum(1e-2, shh[:, k] * E), color=c, lw=1.6, drawstyle="steps-mid", label=lab2)
    day_axis(ax)
    ax.set_yscale("log")
    ax.set_ylim(1, E.max() * 3)
    ax.set_ylabel("energy Σ F²/c", fontsize=10.5)
    w = (shh[:, :3] * E[:, None]).sum(axis=0) / E.sum()
    ax.legend(loc="lower right", ncol=3, fontsize=10, bbox_to_anchor=(1, 1.02), frameon=False)
    ax.set_title(f"The net flow's energy by part through the mean week; weighted by energy, {w[0]:.1%} driven, {w[1]:.1%} round blocks, {w[2]:.1%} round holes", fontsize=12, loc="left")
    source(fig, CREDIT + ". Conductance = width/length. Blocks: orange turn counter-clockwise, blue clockwise.")
    return fig


# ---- 74: the tide of the district ----------------------------------------------------------------------------

def fig74(D):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The tide of the district: its daily line on every street",
          "The D1 line (one cycle a day) of the gross flow, as built: its amplitude as width, the hour of its crest as colour; solved as phasors of the sources")
    p = D.map["as built"]
    present = D.present("as built")
    amp = np.array(p["d1_amplitude"])
    crest = np.array(p["d1_crest_hour"]) + 0.5  # the values stand for the middle of each hour
    ax = fig.add_axes([0.01, 0.14, 0.55, 0.75])
    buildings(ax, D, alpha=0.7)
    vmax = np.percentile(amp[np.array(present)], 98)
    pc = crest[np.array(present) & (amp > 0)]
    lo, hi = np.floor(np.percentile(pc, 3) * 2) / 2, np.ceil(np.percentile(pc, 97) * 2) / 2
    segs, ws, cs = [], [], []
    for i in np.argsort(amp):
        l = D.lines[i]
        if not present[i] or amp[i] <= 0:
            continue
        segs.append(l)
        ws.append(0.5 + 9 * min(1, amp[i] / vmax) ** 0.7)
        cs.append(CREST_SEQ(np.clip((crest[i] - lo) / (hi - lo), 0, 1)))
    ax.add_collection(LineCollection(segs, colors=cs, linewidths=ws, zorder=3, capstyle="round"))
    frame(ax, D)
    tk = np.linspace(lo, hi, 5)
    colourbar(fig, [0.12, 0.09, 0.33, 0.016], CREST_SEQ, lo, hi, "the hour of the daily crest of the walkers", ticks=tk,
              ticklabels=[f"{int(t):02d}:{int(round((t % 1) * 60)):02d}" for t in tk])
    a1 = fig.add_axes([0.63, 0.56, 0.33, 0.31])
    a2 = fig.add_axes([0.63, 0.12, 0.33, 0.31])
    for c in D.cal["counters"]:
        l7 = c["lines"]["7"]
        a1.scatter(l7["observed_amplitude"], l7["model_amplitude"], s=46, color=SEA, edgecolor=INK, lw=0.4, zorder=3)
        a2.scatter(l7["observed_crest_hour"] + 0.5, l7["model_crest_hour"] + 0.5, s=46, color=SAND, edgecolor=INK, lw=0.4, zorder=3)
    lim = max(max(c["lines"]["7"]["observed_amplitude"], c["lines"]["7"]["model_amplitude"]) for c in D.cal["counters"]) * 1.1
    a1.plot([0, lim], [0, lim], color=MUTED, lw=0.8, ls=":")
    a1.set_xlim(0, lim), a1.set_ylim(0, lim)
    a1.set_xlabel("observed amplitude of D1 (people an hour)"), a1.set_ylabel("modelled")
    a1.set_title("At the 16 counters: the daily line's amplitude", fontsize=12, loc="left")
    obs = [c["lines"]["7"]["observed_crest_hour"] + 0.5 for c in D.cal["counters"]]
    mod = [c["lines"]["7"]["model_crest_hour"] + 0.5 for c in D.cal["counters"]]
    lo2, hi2 = math.floor(min(obs + mod)) - 0.5, math.ceil(max(obs + mod)) + 0.5
    a2.plot([lo2, hi2], [lo2, hi2], color=MUTED, lw=0.8, ls=":")
    a2.set_xlim(lo2, hi2), a2.set_ylim(lo2, hi2)
    a2.set_xlabel("observed hour of the daily crest"), a2.set_ylabel("modelled")
    a2.set_title("... and the hour of its crest", fontsize=12, loc="left")
    source(fig, CREDIT + ". Linear in its sources, each line of the week is the sum of the sources' lines times their responses.")
    return fig


# ---- 75: the searches compared ------------------------------------------------------------------------------------

def fig75(D):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    large = max(int(float(r["budget"])) for r in D.runs)
    title(fig, "What the structure bought: the searches at equal budgets",
          f"The value of the best plan found against evaluations of the value (168 margins of an hour each), five seeds; as built {D.summary['as_built']['value_aud_week']:,.0f} A$/week")
    ax = fig.add_axes([0.07, 0.12, 0.52, 0.74])
    curves = {}
    for r in D.anytime:
        if int(float(r["budget"])) != large:
            continue
        curves.setdefault((r["algorithm"], r["seed"]), []).append((float(r["evaluations"]), float(r["best_value"])))
    grid = np.unique(np.round(np.geomspace(50, large, 240)))
    for algo in ORDER:
        ys = []
        for (a, seed), pts in curves.items():
            if a != algo:
                continue
            pts = sorted(pts)
            e = np.array([q[0] for q in pts])
            v = np.array([q[1] for q in pts])
            y = np.full(len(grid), np.nan)
            for k, g in enumerate(grid):
                idx = np.searchsorted(e, g, side="right") - 1
                if idx >= 0:
                    y[k] = v[idx]
            ys.append(y)
        if not ys:
            continue
        Y = np.array(ys)
        ok = ~np.all(np.isnan(Y), axis=0)
        ax.fill_between(grid[ok], np.nanmin(Y[:, ok], axis=0), np.nanmax(Y[:, ok], axis=0), color=ALGO[algo], alpha=0.16, lw=0)
        ax.plot(grid[ok], np.nanmedian(Y[:, ok], axis=0), color=ALGO[algo], lw=2.4, label=algo)
    ax.axhline(D.summary["as_built"]["value_aud_week"], color=MUTED, lw=1, ls=":")
    ax.text(55, D.summary["as_built"]["value_aud_week"], "as built", color=MUTED, va="bottom", fontsize=10)
    ax.set_xscale("log")
    ax.set_xlim(50, large)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel("evaluations of the value")
    ax.set_ylabel("best value found, median and range of five seeds (A$ a week)")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v / 1000:,.0f}k" if v else "0"))
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.62), fontsize=11, frameon=True, framealpha=0.92, facecolor=PAPER, edgecolor=HAIRLINE)
    b = fig.add_axes([0.73, 0.12, 0.25, 0.74])
    budgets = sorted({int(float(r["budget"])) for r in D.runs})
    for j, algo in enumerate(ORDER):
        for bi, budget in enumerate(budgets):
            vals = [float(r["best_value"]) for r in D.runs if r["algorithm"] == algo and int(float(r["budget"])) == budget]
            if not vals:
                continue
            y = j + (0.18 if bi == 0 else -0.18)
            b.scatter(vals, [y] * len(vals), s=32, color=ALGO[algo], alpha=0.45 if bi == 0 else 1.0, edgecolor=INK if bi else "none", lw=0.4, zorder=3)
            b.plot([np.mean(vals)] * 2, [y - 0.14, y + 0.14], color=INK, lw=1.4)
    b.set_yticks(range(len(ORDER)))
    b.set_yticklabels(ORDER, fontsize=11)
    b.axvline(D.summary["as_built"]["value_aud_week"], color=MUTED, lw=1, ls=":")
    b.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v / 1000:,.0f}k" if v else "0"))
    b.set_xlabel("best value of each seed (A$ a week)")
    b.set_title(f"Pale: {budgets[0]:,} evaluations; dark: {budgets[-1]:,}", fontsize=11.5, loc="left")
    b.grid(axis="y", visible=False)
    source(fig, "heuresis: GA, separable CMA-ES and random search as study (1) runs them; BR: exact best responses of each frontage (6 evaluations each) with a descent over the links; GA with kernels biased by each frontage's marginal value.")
    return fig


# ---- 76: the best plan ----------------------------------------------------------------------------------------------

def kinds_panel(ax, D, genotype, outcome, label, spend_flux=True):
    buildings(ax, D, parcels=False, alpha=0.9)
    network(ax, D, None, lw=0.5)
    kinds = genotype["kind"]
    polys = [np.array(f["footprint"]) for f in D.geo["frontages"]]
    ax.add_collection(PolyCollection(polys, facecolors=[KIND_COLOUR[int(k)] for k in kinds], edgecolors=PAPER, linewidths=0.6, zorder=3))
    links(ax, D, made=genotype["links"], all_candidates=False, lw=3)
    if spend_flux:
        spend = outcome["spending_by_frontage"]
        smax = max(spend) or 1
        for f, fr in enumerate(D.geo["frontages"]):
            if spend[f] <= 0:
                continue
            a = np.array(fr["anchor"])
            c = np.array(fr["centre"])
            d = c - a
            n = np.hypot(*d)
            if n < 1e-6:
                continue
            d = d / n
            L = 6 + 22 * (spend[f] / smax) ** 0.6
            ax.add_patch(FancyArrowPatch(a - d * L * 0.6, a + d * L * 0.4, arrowstyle="-|>", mutation_scale=5 + 8 * (spend[f] / smax) ** 0.5, lw=0.6 + 2 * (spend[f] / smax) ** 0.6,
                                         color=INK, zorder=6, alpha=0.9))
    q = np.array(D.geo["quad"])
    ax.set_xlim(q[:, 0].min() - 40, q[:, 0].max() + 40)
    ax.set_ylim(q[:, 1].min() - 30, q[:, 1].max() + 30)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(label, fontsize=13, loc="left")


def fig76(D):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    a, b, k = D.plan["as built"], D.plan["best"], D.plan["best black box"]
    title(fig, "The plans: what each frontage holds, and the spending it takes in",
          f"Arrows: the week's spending as flux from the street into each frontage. As built {a['outcome']['value_aud_week']:,.0f}, the best black-box plan {k['outcome']['value_aud_week']:,.0f}, "
          f"the best plan {b['outcome']['value_aud_week']:,.0f} A$/week")
    for i, (p, label) in enumerate([(a, "as built"), (k, f"best black box ({D.summary['best_black_box_algorithm']})"), (b, f"best ({D.summary['best_algorithm']})")]):
        ax = fig.add_axes([0.01 + i * 0.33, 0.12, 0.32, 0.74])
        kinds_panel(ax, D, p["genotype"], p["outcome"], label)
        o = p["outcome"]
        ax.text(0.02, 0.0, f"spending {o['spending_aud_week'] / 1000:,.0f}k · open {o['open_frontages']:.0f} · links {o['links']:.0f}\nmix {o['mix_nats']:.2f} nats · access {o['access_m']:.0f} m · spread {o['spread']:.2f}",
                transform=ax.transAxes, fontsize=10, color=MUTED, va="top")
    handles = [Patch(color=KIND_COLOUR[k], label=KIND_NAME[k]) for k in range(6)] + [Line2D([0], [0], color=SAND, lw=3, ls=(0, (3, 1.5)), label="link made")]
    fig.legend(handles=handles, loc="lower center", ncol=7, fontsize=11, bbox_to_anchor=(0.5, 0.035))
    source(fig, CREDIT + ". The value is the shared problem's (economy.hpp); the passers are this study's flow.")
    return fig


# ---- 77: the week ---------------------------------------------------------------------------------------------------

def fig77(D):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The week: people, spending and the frontages open", "Hour by hour through the mean week, as built and under the best plans; and the spending's lines of the week")
    h = np.arange(168) + 0.5
    names = [("as built", MUTED, (0, (2, 1.5))), ("best black box", ALGO["GA"], "-"), ("best", SEA, "-")]
    ax1 = fig.add_axes([0.06, 0.66, 0.9, 0.2])
    ax2 = fig.add_axes([0.06, 0.40, 0.9, 0.2])
    ax3 = fig.add_axes([0.06, 0.08, 0.42, 0.24])
    ax4 = fig.add_axes([0.55, 0.08, 0.41, 0.24])
    for n, c, ls in names:
        p = D.plan[n]
        ax1.plot(h, np.array(p["passers_by_hour"]) / 1000, color=c, lw=1.6, ls=ls, label=n)
        ax2.plot(h, np.array(p["outcome"]["spending_by_hour"]) / 1000, color=c, lw=1.6, ls=ls, label=n)
        ax3.plot(h, p["open_trading_by_hour"], color=c, lw=1.6, ls=ls, drawstyle="steps-mid", label=n)
        pw = np.array(p["outcome"]["power"])
        m = np.arange(len(pw))
        tot = pw.sum()
        off = {"as built": -0.25, "best black box": 0.0, "best": 0.25}[n]
        ax4.vlines(m[1:43] + off, 1e-6, np.maximum(1e-6, pw[1:43] / tot), color=c, lw=1.6)
        ax4.scatter(m[1:43] + off, np.maximum(1e-6, pw[1:43] / tot), s=10, color=c, zorder=3)
    for a in (ax1, ax2, ax3):
        day_axis(a)
    ax1.set_ylabel("thousand passings\nan hour, 61 frontages", fontsize=10.5)
    ax2.set_ylabel("thousand A$ an hour", fontsize=10.5)
    ax3.set_ylabel("trading frontages open", fontsize=10.5)
    ax1.set_title("The people passing the frontages", fontsize=12, loc="left")
    ax2.set_title("The spending they leave", fontsize=12, loc="left")
    ax1.legend(loc="upper left", ncol=3, fontsize=10)
    for k in range(1, 7):
        ax4.axvline(7 * k, color=HAIRLINE, lw=0.6, zorder=0)
    ax4.set_xlim(0.5, 42)
    ax4.set_xticks([1, 7, 14, 21, 28, 35, 42])
    ax4.set_xticklabels(["W1", "D1", "D2", "D3", "D4", "D5", "D6"])
    ax4.set_yscale("log")
    ax4.set_ylim(1e-5, 1)
    ax4.set_ylabel("share of the power", fontsize=10.5)
    ax4.set_title("The spending's lines (Parseval; the mean line left out)", fontsize=12, loc="left")
    source(fig, "The week's 168 hours from Monday 00:00 civil; lines m = 0…84 folded, D_j at m = 7j; the value, spread and lines are the shared problem's measures.")
    return fig


# ---- 78: what a link does -------------------------------------------------------------------------------------------

def fig78(D):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    L = D.plans["links"]
    title(fig, "What a link does: the people it carries, the frontages it passes by, and the value",
          "Each candidate link made alone with the plan as built (no frontage changed), and each flipped from the best plan with the frontages re-optimised")
    ax = fig.add_axes([0.25, 0.14, 0.38, 0.72])
    y = np.arange(len(L))
    on = np.array([l["on_link_week"] for l in L]) / 1000
    gain = np.array([l["frontage_passers_gained_week"] for l in L]) / 1000
    lost = np.array([l["frontage_passers_lost_week"] for l in L]) / 1000
    ax.barh(y + 0.2, on, height=0.36, color=SEA, label="walkers on the link")
    ax.barh(y - 0.2, gain, height=0.36, color="#7FA37A", label="passings gained at frontages")
    ax.barh(y - 0.2, lost, height=0.36, color=SAND, label="passings lost at frontages")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{l['joins'][0]} to {l['joins'][1]} ({'passage' if l['kind'] == 'passage' else 'lane'}, {l['length_m']:.0f} m)" for l in L], fontsize=10.5)
    ax.axvline(0, color=INK, lw=0.8)
    ax.set_xlabel("thousand people a week")
    ax.legend(loc="lower left", fontsize=10)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    b = fig.add_axes([0.68, 0.14, 0.29, 0.72])
    v1 = np.array([l["value_change_as_built"] for l in L])
    v2 = np.array([l["value_change_from_best_flipped_and_settled"] for l in L])
    cost = np.array([l["weekly_cost_aud"] for l in L])
    b.barh(y + 0.2, v1, height=0.36, color=MUTED, label="made with the plan as built")
    b.barh(y - 0.2, v2, height=0.36, color=[SEA if l["in_best_plan"] else SAND for l in L], label="flipped from the best plan, re-optimised")
    b.scatter(-cost, y + 0.2, marker="|", s=140, color=INK, zorder=4, label="minus the link's weekly cost")
    b.axvline(0, color=INK, lw=0.8)
    b.set_yticks(y)
    b.set_yticklabels([])
    b.invert_yaxis()
    b.set_xlabel("change in value (A$ a week)")
    b.legend(loc="lower left", fontsize=9.5)
    b.grid(axis="y", visible=False)
    b.set_title("Blue: in the best plan, so its removal; sand: its addition", fontsize=11, loc="left")
    source(fig, "A link reroutes every stream at once (Woodbury's update of the walks): it adds passers on some frontages' segments and takes them from others.")
    return fig


# ---- 79: arrangement -----------------------------------------------------------------------------------------------

def fig79(D):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    arr = {a["plan"]: a for a in D.plans["arrangement"]}
    b = D.plan["best"]
    title(fig, "The arrangement the best plan uses: catalysis, competition and party walls",
          f"Catalysis adds {arr['best']['value'] - arr['best']['value_without_catalysis']:,.0f} of the best plan's {arr['best']['value']:,.0f} A$/week; "
          f"{arr['best']['same_kind_same_segment_pairs']:.0f} pairs of one kind share a segment (as built {arr['as built']['same_kind_same_segment_pairs']:.0f})")
    ax = fig.add_axes([0.01, 0.1, 0.58, 0.78])
    kinds = b["genotype"]["kind"]
    buildings(ax, D, parcels=False, alpha=0.9)
    network(ax, D, None, lw=0.5)
    polys = [np.array(f["footprint"]) for f in D.geo["frontages"]]
    ax.add_collection(PolyCollection(polys, facecolors=[KIND_COLOUR[int(k)] for k in kinds], edgecolors=PAPER, linewidths=0.6, zorder=3, alpha=0.9))
    centres = np.array([f["centre"] for f in D.geo["frontages"]])
    rel = D.problem["relations"]
    aff = D.problem["affinity"]["a"]
    trades = lambda k: int(k) in (1, 2, 3, 5)
    for f, g, d in rel["reach_within_3_coupling_lengths"]:
        f, g = int(f), int(g)
        if not (trades(kinds[f]) and trades(kinds[g])):
            continue
        a = aff[int(kinds[g])][int(kinds[f])] + aff[int(kinds[f])][int(kinds[g])]
        if a <= 0:
            continue
        ax.plot(*zip(centres[f], centres[g]), color="#8E3B78", lw=0.4 + 1.4 * math.exp(-d / 42.5) * a, alpha=0.55, zorder=4)
    for f, g, area in rel["party_walls"]:
        f, g = int(f), int(g)
        same = kinds[f] == kinds[g] and trades(kinds[f])
        ax.plot(*zip(centres[f], centres[g]), color=INK if same else MUTED, lw=1.6 if same else 0.6, ls="-" if same else ":", zorder=5)
    q = np.array(D.geo["quad"])
    ax.set_xlim(q[:, 0].min() - 30, q[:, 0].max() + 30)
    ax.set_ylim(q[:, 1].min() - 25, q[:, 1].max() + 25)
    ax.set_aspect("equal")
    ax.axis("off")
    handles = [Patch(color=KIND_COLOUR[k], label=KIND_NAME[k]) for k in range(6)] + [
        Line2D([0], [0], color="#8E3B78", lw=1.4, label="catalysing pair within reach"), Line2D([0], [0], color=INK, lw=1.6, label="party wall, one trading kind"),
        Line2D([0], [0], color=MUTED, lw=0.8, ls=":", label="party wall, other")]
    ax.legend(handles=handles, loc="lower left", fontsize=9.5, frameon=True, framealpha=0.92, facecolor=PAPER, edgecolor=HAIRLINE)
    c = fig.add_axes([0.66, 0.55, 0.31, 0.3])
    names = ["as built", "best black box", "best"]
    v = [arr[n]["value"] for n in names]
    v0 = [arr[n]["value_without_catalysis"] for n in names]
    y = np.arange(3)
    c.barh(y + 0.18, v0, color=MUTED, height=0.34, label="without catalysis (the same flows)")
    c.barh(y - 0.18, v, color="#8E3B78", height=0.34, label="with it")
    c.axvline(0, color=INK, lw=0.8)
    c.set_yticks(y)
    c.set_yticklabels(names)
    c.invert_yaxis()
    c.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x / 1000:,.0f}k" if x else "0"))
    c.set_xlabel("value (A$ a week)")
    c.legend(loc="upper center", bbox_to_anchor=(0.5, 1.2), ncol=2, fontsize=9.5, frameon=False)
    c.grid(axis="y", visible=False)
    d = fig.add_axes([0.66, 0.12, 0.31, 0.3])
    keys = [("same_kind_same_segment_pairs", "pairs of one kind\non one segment"), ("positively_catalysing_pairs_in_reach", "catalysing pairs\nwithin reach"), ("party_walls_same_kind", "party walls between\none trading kind")]
    w = 0.26
    for i, n in enumerate(names):
        d.bar(np.arange(3) + (i - 1) * w, [arr[n][k] for k, _ in keys], width=w, color=[MUTED, ALGO["GA"], SEA][i], label=n)
    d.set_xticks(range(3))
    d.set_xticklabels([t for _, t in keys], fontsize=10)
    d.legend(fontsize=9.5)
    d.grid(axis="x", visible=False)
    source(fig, "Relations from the shared problem (§2.3): party walls from architecture's massing, reach within 3ξ = 127.6 m, affinities ln CLQ of CLUE 2024.")
    return fig


FIGURES = {70: ("search-70-the-district-as-a-network", fig70), 71: ("search-71-counters-and-flow", fig71), 72: ("search-72-people-on-the-streets", fig72),
           73: ("search-73-circulation", fig73), 74: ("search-74-the-tide-of-the-district", fig74), 75: ("search-75-what-the-structure-bought", fig75),
           76: ("search-76-the-plans-and-their-spending", fig76), 77: ("search-77-the-week-of-the-plans", fig77), 78: ("search-78-what-a-link-does", fig78),
           79: ("search-79-the-arrangement-used", fig79)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--results", default=str(HERE / "results"))
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    D = Data(pathlib.Path(a.results))
    only = {int(x) for x in a.only.split(",") if x}
    for k, (name, f) in FIGURES.items():
        if only and k not in only:
            continue
        fig = f(D)
        if a.out:
            pathlib.Path(a.out).mkdir(parents=True, exist_ok=True)
            fig.savefig(pathlib.Path(a.out) / f"{name}.jpg", dpi=100, pil_kwargs={"quality": 90})
            plt.close(fig)
            print(pathlib.Path(a.out) / f"{name}.jpg")
        else:
            print(save(fig, name))


if __name__ == "__main__":
    main()
