"""The figures of "The walkers of the Hoddle Grid": agents-01 to agents-09, from results/ (and results/maps.json, the
geometry, which walkers.cpp writes and which is not committed).

House style: research/experiments/style/search_style.py. The place in the Urban Science Workbench's colours: a gentle
sea blue (#5B8FB9) for the walkers' water and sand (#C98B6B) for the shore; deeper steps of both for series, validated
for colour vision deficiencies with the dataviz validator (walkers #2F6FA6, flow #B4643A, plum #8E3B78, green #3D8C5A).

    .venv/Scripts/python.exe research/experiments/2026-11-05-walkers-of-the-hoddle-grid/figures.py
"""

from __future__ import annotations

import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "style"))
from search_style import HAIRLINE, INK, MUTED, PANEL, PAPER, figure, save, source, title  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize, TwoSlopeNorm  # noqa: E402

R = HERE / "results"
SEA, SAND = "#5B8FB9", "#C98B6B"
WALKERS, FLOW, PLUM, GREEN = "#2F6FA6", "#B4643A", "#8E3B78", "#3D8C5A"
OBSERVED = INK
SEA_MAP = LinearSegmentedColormap.from_list("sea", ["#EEF3F7", "#C9DAE8", "#8FB3D3", "#5B8FB9", "#2F6FA6", "#1D3F63"])
SEA_LIGHT = LinearSegmentedColormap.from_list("sea_light", ["#F1ECE7", "#DCE6EE", "#B7CDE0", "#8FB3D3", "#6C9AC2"])
SHORE = LinearSegmentedColormap.from_list("shore", ["#7A3E1E", "#B4643A", "#DDB79F", "#EDE8E3", "#B9CFE2", "#5B8FB9", "#1D4F7C"])
DAYS = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
SRC = "City of Melbourne pedestrian counts (CC BY 4.0); OpenStreetMap (ODbL) via the data library's city model; urban-science movement on nomos agents and tyche decision"


def load(name):
    with open(R / name, encoding="utf-8") as f:
        return json.load(f)


MAPS = load("maps.json")
SEG = MAPS["segments"]
BOX = MAPS["box"]


def lines():
    return [np.array(s["line"]) for s in SEG]


def base_map(ax, values=None, cmap=SEA_MAP, norm=None, widths=None, under=True, lw_scale=0.55):
    """The district's segments, coloured by a value per segment (or in shore grey), widths by walkable width."""
    ls = lines()
    w = np.array([s["width"] for s in SEG], dtype=float)
    lw = (0.8 + lw_scale * w) if widths is None else widths
    if under:
        ax.add_collection(LineCollection(ls, colors="#E7DED6", linewidths=lw + 2.2, capstyle="round", zorder=1))
    if values is None:
        ax.add_collection(LineCollection(ls, colors="#D8CEC6", linewidths=lw, capstyle="round", zorder=2))
    else:
        lc = LineCollection(ls, array=np.asarray(values, dtype=float), cmap=cmap, norm=norm, linewidths=lw, capstyle="round", zorder=2)
        ax.add_collection(lc)
        return lc
    return None


def frame(ax, pad=10):
    ax.set_xlim(BOX[0] - pad, BOX[2] + pad)
    ax.set_ylim(BOX[1] - pad, BOX[3] + pad)
    ax.set_aspect("equal")
    ax.set_xticks([]), ax.set_yticks([])
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def street_labels(ax, names=("Swanston Street", "Elizabeth Street", "Russell Street", "Collins Street", "Bourke Street", "Little Bourke Street", "Lonsdale Street",
                             "Little Collins Street", "Flinders Lane", "Flinders Street")):
    """One label per named street, at the middle of its longest segment, along it."""
    for n in names:
        best = None
        for s in SEG:
            if s["name"] != n:
                continue
            l = np.array(s["line"])
            L = np.sum(np.hypot(*np.diff(l, axis=0).T))
            if best is None or L > best[0]:
                best = (L, l)
        if best is None:
            continue
        l = best[1]
        m = l[len(l) // 2] if len(l) > 2 else l.mean(axis=0)
        d = l[-1] - l[0]
        ang = math.degrees(math.atan2(d[1], d[0]))
        if ang > 90:
            ang -= 180
        if ang < -90:
            ang += 180
        ax.text(m[0], m[1], n, fontsize=8.5, color=MUTED, rotation=ang, rotation_mode="anchor", ha="center", va="center", zorder=6,
                bbox=dict(boxstyle="square,pad=0.08", fc=PAPER, ec="none", alpha=0.75))


def week_axis(ax, label=True):
    ax.set_xlim(0, 168)
    ax.set_xticks(range(0, 169, 24))
    ax.set_xticklabels([""] * 8)
    for d in range(7):
        ax.text(24 * d + 12, -0.06, DAYS[d], transform=ax.get_xaxis_transform(), ha="center", va="top", color=MUTED, fontsize=11)
    if label:
        ax.set_xlabel("hour of the mean week (civil time)", labelpad=18)


# ======================================================================================================================
def fig01():
    W = load("week.json")
    import csv
    rows = list(csv.DictReader(open(R / "week-segments.csv", encoding="utf-8")))
    passes = np.array([float(r["walkers_week_by_entry"]) for r in rows])
    fig, axes = figure(1, 2)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.14, wspace=0.03)
    norm = LogNorm(vmin=max(1, np.percentile(passes[passes > 0], 3)), vmax=passes.max())
    snaps = MAPS["snapshots"]
    for ax, (when, pts) in zip(axes, snaps.items()):
        lc = base_map(ax, np.maximum(passes, 1), cmap=SEA_LIGHT, norm=norm, lw_scale=0.7)
        p = np.array([[q[0], q[1]] for q in pts])
        st = np.array([q[3] for q in pts])
        walking = st == 1
        ax.scatter(p[walking, 0], p[walking, 1], s=3.2, color="#1D3F63", alpha=0.75, lw=0, zorder=4)
        ax.scatter(p[~walking, 0], p[~walking, 1], s=6.0, color=PLUM, alpha=0.9, lw=0, zorder=5)
        frame(ax)
        street_labels(ax)
        ax.set_title(f"{when}: {len(pts):,} walkers on the streets", loc="left", color=INK)
        for c in MAPS["counters"]:
            ax.plot(c["x"], c["y"], marker="s", ms=5.5, mfc=PAPER, mec=FLOW, mew=1.4, zorder=7)
    cax = fig.add_axes([0.3, 0.105, 0.4, 0.018])
    cb = fig.colorbar(lc, cax=cax, orientation="horizontal")
    cb.set_label("walkers entering each segment over the week (both ways)", color=INK)
    cb.outline.set_visible(False)
    t = W["totals"]
    title(fig, "Walkers of a mean week, each an agent",
          f"{t['arrived']:,.0f} trips in the week, each with its origin, its target and a way chosen at every junction by the logit over walks.\n"
          f"Dots: walkers walking (deep sea) and waiting at a signal or for room (plum); squares: the 16 counters the trips were calibrated to.")
    source(fig, SRC + f"; one week at 2 s steps, {W['run_s']:.0f} s on {W['threads']:.0f} threads")
    save(fig, "agents-01-the-walkers-of-a-week")


def fig02():
    W = load("week.json")
    import csv
    rows = list(csv.DictReader(open(R / "week-segments.csv", encoding="utf-8")))
    flow = np.array([float(r["flow_week"]) for r in rows])
    dep = np.array([float(r["walkers_week_by_departure"]) for r in rows])
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(2, 3, left=0.06, right=0.97, top=0.84, bottom=0.12, wspace=0.28, hspace=0.42)
    ax = fig.add_subplot(gs[:, 0])
    ok = flow > 0
    ax.loglog(flow[ok], dep[ok], "o", ms=4, color=WALKERS, alpha=0.6, mec="none")
    lim = [max(1, flow[ok].min() * 0.7), flow.max() * 1.4]
    ax.plot(lim, lim, color=FLOW, lw=1.4)
    ax.set_xlim(lim), ax.set_ylim(lim)
    ax.set_xlabel("the logit flow over the week (walker-segments)")
    ax.set_ylabel("the walkers' entries over the week, by hour of departure")
    ax.set_title("each segment: the walkers on the flow", loc="left")
    mm = W["micro_macro"]
    ax.text(0.04, 0.96, f"{mm['walker_entries']:,.0f} entries\nagainst {mm['flow_entries']:,.0f}", transform=ax.transAxes, va="top", color=INK, fontsize=12)
    # The standardised deviations of every state-hour, against N(0, 1).
    ax2 = fig.add_subplot(gs[0, 1:])
    z = None
    ax2.set_title(f"every segment-way and hour with μ ≥ 5 ({mm['state_hours']:,.0f}): (n − μ)/√μ", loc="left")
    xs = np.linspace(-5, 5, 400)
    m, v = mm["mean_z"], mm["var_z"]
    hist = np.array(mm.get("z_histogram", []), dtype=float)
    if hist.size:
        edges = np.linspace(-5, 5, hist.size + 1)
        dens = hist / (hist.sum() * 0.25)
        ax2.bar(edges[:-1] + 0.125, dens, width=0.23, color=SEA, alpha=0.85, label=f"the walkers: mean {m:+.3f}, variance {v:.3f}")
    ax2.plot(xs, np.exp(-xs ** 2 / 2) / math.sqrt(2 * math.pi), color=FLOW, lw=2, label="standard normal (Corollary 3.5)")
    ax2.set_xlabel("standardised deviation")
    ax2.set_yticks([])
    ax2.legend(loc="upper right")
    ax2.text(0.02, 0.9, f"beyond ±3: {100 * mm['beyond_3']:.2f}% (normal: 0.27%)\nχ²/dof {mm['chi2_per_dof']:.3f}", transform=ax2.transAxes, va="top", fontsize=12, color=INK)
    ax3 = fig.add_subplot(gs[1, 1:])
    e = W["entries_by_hour"]
    f = np.array(e["flow"]); a = np.array(e["walkers_by_entry"]); d = np.array(e["walkers_by_departure"])
    hrs = np.arange(168) + 0.5
    ax3.plot(hrs, (a - f) / 1000, color=WALKERS, lw=1.6, label="by the hour they enter a segment (the walking takes time)")
    ax3.plot(hrs, (d - f) / 1000, color=PLUM, lw=1.2, label="by the hour they set out (exact in expectation)")
    ax3.axhline(0, color=FLOW, lw=1)
    ax3.set_ylabel("walkers − flow (thousand entries an hour)")
    ax3.set_title("hour by hour: the lag of walking against the quasi-static flow", loc="left")
    week_axis(ax3)
    ax3.legend(loc="lower right", fontsize=10.5)
    title(fig, "The walkers' mean is the logit flow", "Their route choice is the flow study's law as a Markov chain on directed segments; counted by departure the two agree exactly in expectation")
    source(fig, "City of Melbourne pedestrian counts (CC BY 4.0); OpenStreetMap (ODbL); the flow: research/experiments/2026-11-05-a-district-economy-as-a-flow, its calibration reproduced")
    save(fig, "agents-02-the-walkers-and-the-flow")


def fig03():
    W = load("week.json")
    cs = {c["id"]: c for c in W["counters"]}
    pick = [4, 185, 19, 21, 53, 72]
    fig, axes = figure(3, 2, sharex=True)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.82, bottom=0.08, hspace=0.62, wspace=0.16)
    for ax, cid in zip(axes.ravel(), pick):
        c = cs[cid]
        h = np.arange(168) + 0.5
        ax.fill_between(h, 0, c["observed"], color="#E9E1DA", lw=0, step=None)
        ax.plot(h, c["observed"], color=OBSERVED, lw=1.2, label="counted (mean week, 2025–26)")
        ax.plot(h, c["flow"], color=FLOW, lw=1.6, label="the logit flow")
        ax.plot(h, c["walkers"], color=WALKERS, lw=1.1, alpha=0.9, label="the walkers (one simulated week)")
        ax.set_title(f"{c['name']} — R² of the walkers {c['r2_walkers_observed']:.2f}, of the flow {c['r2_flow_observed']:.2f}", loc="left", fontsize=12.5)
        ax.set_ylim(0, max(max(c["observed"]), max(c["walkers"])) * 1.12)
        week_axis(ax, label=False)
        ax.set_ylabel("people an hour")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", ncol=3, fontsize=11, bbox_to_anchor=(0.035, 0.895))
    r = W["r2"]
    title(fig, "The counters' week, counted and walked",
          f"All 16 counters × 168 hours: R² of the walkers against the counts {r['walkers_observed']:.3f}, of the flow {r['flow_observed']:.3f}; of the walkers against the flow {r['walkers_flow']:.3f}")
    source(fig, SRC)
    save(fig, "agents-03-the-counters-week")


def fig04():
    W = load("week.json")
    x = np.array(W["present_by_minute"])
    n = len(x)
    X = np.fft.rfft(x - x.mean())
    P = np.abs(X) ** 2
    f = np.fft.rfftfreq(n, d=1 / (24 * 60))   # cycles a day
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(2, 2, left=0.06, right=0.97, top=0.84, bottom=0.11, wspace=0.22, hspace=0.42, height_ratios=[1, 1.15])
    ax = fig.add_subplot(gs[0, :])
    Pn = P[1:] / P[1:].max()
    fk = f[1:]
    daily = np.isclose(fk % 1, 0, atol=1e-6) | np.isclose(fk % 1, 1, atol=1e-6)
    ax.vlines(fk[~daily], 1e-8, Pn[~daily], color=SEA, lw=1.6)
    ax.vlines(fk[daily], 1e-8, Pn[daily], color=WALKERS, lw=3)
    ax.plot(fk[daily], Pn[daily], "o", color=WALKERS, ms=5)
    ax.set_yscale("log")
    ax.set_ylim(1e-7, 3)
    ax.set_xlim(0, 12.3)
    ax.set_xlabel("frequency (cycles a day)")
    ax.set_ylabel("power (of its largest)")
    ax.set_title("the spectrum of the walkers on the district's streets, minute by minute over the week", loc="left")
    for k, name, dx in [(1 / 7, "W1", 0), (1, "D1", 0), (2, "D2", 0), (3, "D3", 0), (4, "D4", 0), (5, "D5", 0), (6, "D6", 0), (6 / 7, "D1−W1", -16), (8 / 7, "D1+W1", 16)]:
        i = np.argmin(np.abs(f - k))
        ax.annotate(name, (f[i], P[i] / P[1:].max()), xytext=(dx, 7), textcoords="offset points", ha="center", fontsize=10, color=INK)
    ax.text(0.99, 0.95, "dark: the day and its harmonics; light: the week's lines and their sidebands", transform=ax.transAxes, ha="right", va="top", fontsize=10.5, color=MUTED)
    # The daily line at each counter: amplitude and crest, counted against walked.
    cs = W["counters"]
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    for c in cs:
        o, w, fl = c["lines"]["observed"]["D1"], c["lines"]["walkers"]["D1"], c["lines"]["flow"]["D1"]
        ax2.plot([o["amplitude"]], [w["amplitude"]], "o", color=WALKERS, ms=7)
        ax2.plot([o["amplitude"]], [fl["amplitude"]], "s", color=FLOW, ms=5, mfc="none", mew=1.3)
        ax3.plot([o["crest_h"]], [w["crest_h"]], "o", color=WALKERS, ms=7)
        ax3.plot([o["crest_h"]], [fl["crest_h"]], "s", color=FLOW, ms=5, mfc="none", mew=1.3)
    m = max(max(c["lines"]["observed"]["D1"]["amplitude"] for c in cs), max(c["lines"]["walkers"]["D1"]["amplitude"] for c in cs)) * 1.1
    ax2.plot([0, m], [0, m], color=HAIRLINE, lw=1)
    ax2.set_xlim(0, m), ax2.set_ylim(0, m)
    ax2.set_xlabel("counted D1 amplitude (people an hour)")
    ax2.set_ylabel("walked D1 amplitude")
    ax2.set_title("the daily line's amplitude at the 16 counters", loc="left")
    ax2.plot([], [], "o", color=WALKERS, label="the walkers"), ax2.plot([], [], "s", color=FLOW, mfc="none", label="the logit flow")
    ax2.legend(loc="upper left")
    ax3.plot([12, 19], [12, 19], color=HAIRLINE, lw=1)
    ax3.set_xlim(12.5, 18.5), ax3.set_ylim(12.5, 18.5)
    ax3.set_xlabel("counted crest of D1 (hour of the day)")
    ax3.set_ylabel("walked crest")
    ax3.set_title("its crest: one tide for the whole district", loc="left")
    title(fig, "The week's lines", "The walkers keep the counted week's frequencies — the day, its harmonics and the week's modulation of it — but, like the flow, one hour of high water for every street")
    source(fig, SRC)
    save(fig, "agents-04-the-lines-of-the-week")


# ======================================================================================================================
def fig05():
    C = load("closure.json")
    import csv
    rows = list(csv.DictReader(open(R / "closure-segments.csv", encoding="utf-8")))
    wo = np.array([float(r["walkers_open"]) for r in rows]); wc = np.array([float(r["walkers_closed"]) for r in rows])
    fo = np.array([float(r["flow_open"]) for r in rows]); fc = np.array([float(r["flow_closed"]) for r in rows])
    d = wc - wo
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(2, 2, left=0.03, right=0.97, top=0.82, bottom=0.08, wspace=0.12, hspace=0.35, width_ratios=[1.35, 1])
    ax = fig.add_subplot(gs[:, 0])
    lim = np.percentile(np.abs(d), 99.5)
    lc = base_map(ax, np.clip(d, -lim, lim), cmap=SHORE, norm=TwoSlopeNorm(0, -lim, lim))
    blk = set(int(e) for e in C["edges"])
    for e in blk:
        l = np.array(SEG[e]["line"])
        ax.plot(l[:, 0], l[:, 1], color=PLUM, lw=2, ls=(0, (2, 2)), zorder=5)
    frame(ax)
    street_labels(ax)
    cb = fig.colorbar(lc, ax=ax, orientation="horizontal", fraction=0.04, pad=0.02, aspect=40)
    cb.set_label("walkers entering each segment on the Sunday, closed − open", color=INK)
    cb.outline.set_visible(False)
    ax.set_title("where the walkers go instead (the closed block dashed)", loc="left")
    ax2 = fig.add_subplot(gs[0, 1])
    ok = np.abs(fc - fo) > 1
    ax2.plot(fc[ok] - fo[ok], d[ok], "o", color=WALKERS, ms=4.5, alpha=0.7, mec="none")
    m = np.abs(fc - fo).max() * 1.05
    ax2.plot([-m, m], [-m, m], color=FLOW, lw=1.3)
    ax2.set_xlabel("change in the logit flow (exact)")
    ax2.set_ylabel("change in the walkers' entries")
    ax2.set_title("the walkers' change against the flow's, segment by segment", loc="left", fontsize=13)
    ax3 = fig.add_subplot(gs[1, 1])
    po, pc = np.array(C["person_hours_by_hour"]["open"]), np.array(C["person_hours_by_hour"]["closed"])
    h = np.arange(24) + 0.5
    ax3.bar(h, pc - po, width=0.8, color=FLOW)
    ax3.axhline(0, color=INK, lw=0.8)
    ax3.set_xlim(0, 24), ax3.set_xticks(range(0, 25, 3))
    ax3.set_xlabel("hour of the Sunday")
    ax3.set_ylabel("person-hours on the streets, closed − open")
    ax3.set_title("what the closure costs the walkers, hour by hour", loc="left", fontsize=13)
    o, c = C["open"], C["closed"]
    title(fig, "Little Bourke Street closed to walkers passing through",
          f"As the City's permit closed it for the Lunar New Year (27–29 January 2017), on the mean Sunday. Of the {C['flow_on_block']:,.0f} passes the street carries, "
          f"{C['flow_moved_elsewhere']:,.0f} move to other streets;\nits own walkers keep it. Mean trip {o['mean_travel_s']:.0f} s open, {c['mean_travel_s']:.0f} s closed: "
          f"{c['person_hours'] - o['person_hours']:+,.0f} person-hours on the day.")
    source(fig, "City of Melbourne counts and event permits 2014–2018 (CC BY 4.0); OpenStreetMap (ODbL); common random numbers: every walker draws alike in both runs")
    save(fig, "agents-05-a-street-closed")


def fig06():
    F = load("festival.json")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(2, 3, left=0.06, right=0.97, top=0.82, bottom=0.1, wspace=0.28, hspace=0.45)
    fitted = [c for c in F["counters"] if c["fitted"]]
    order = sorted(fitted, key=lambda c: -np.nansum(np.array(c["observed_excess"], dtype=float)))[:3]
    h = np.arange(24) + 0.5
    for i, c in enumerate(order):
        ax = fig.add_subplot(gs[0, i])
        ob = np.array([np.nan if v is None else v for v in c["observed_2017"]], dtype=float)
        od = np.array([np.nan if v is None else v for v in c["ordinary_sundays_2017"]], dtype=float)
        ax.fill_between(h, 0, od, color="#E9E1DA", lw=0, label="ordinary Sundays about it (mean of six)")
        ax.plot(h, ob, color=OBSERVED, lw=1.6, label="Sunday 29 January 2017")
        ex = np.array(c["walkers_excess"], dtype=float)
        ax.plot(h, od + ex, color=WALKERS, lw=1.8, label="ordinary + the walkers' excess")
        ax.set_title(c["name"], loc="left", fontsize=12.5)
        ax.set_xlim(6, 24), ax.set_xticks(range(6, 25, 3))
        ax.set_xlabel("hour")
        if i == 0:
            ax.set_ylabel("people an hour")
            handles, labels = ax.get_legend_handles_labels()
            fig.legend(handles, labels, loc="upper left", ncol=3, fontsize=11, bbox_to_anchor=(0.035, 0.895))
    ax = fig.add_subplot(gs[1, :2])
    a = np.array(F["arrivals_by_hour"])
    ae = np.array(F.get("arrivals_from_entrances", a)); ab = np.array(F.get("arrivals_from_buildings", np.zeros(24)))
    ax.bar(h, ae, width=0.8, color=SEA, label="from the district's entrances")
    ax.bar(h, ab, width=0.8, bottom=ae, color=SAND, label="from its buildings")
    ax.legend(loc="upper right")
    ax.set_xlim(6, 24), ax.set_xticks(range(6, 25, 1))
    ax.set_xlabel("hour of arrival")
    ax.set_ylabel("visitors arriving an hour")
    ax.set_title(f"the festival's visitors, fitted: {a.sum():,.0f} arrive, each dwelling {F['dwell_h'] * 60:.0f} minutes on average", loc="left")
    ax2 = fig.add_subplot(gs[1, 2])
    g = F["dwell_grid"]
    ax2.plot([q["dwell_h"] * 60 for q in g], [q["r2"] for q in g], "o-", color=PLUM)
    ax2.set_xlabel("mean dwell (minutes)")
    ax2.set_ylabel("R² of the counters' excess")
    ax2.set_title("the dwell, chosen by the fit", loc="left")
    best = max(F["dwell_grid"], key=lambda q: q["r2"])["r2"]
    title(fig, "The Lunar New Year's crowd, fitted to the counters",
          f"Visitors from the district's entrances and buildings to Little Bourke Street and back, their arrivals and dwell fitted by NNLS to the excess at {F['fit_counters']:.0f} counters: R² {best:.2f}")
    source(fig, "City of Melbourne pedestrian counts of 2017 (CC BY 4.0); ordinary Sundays 8, 15, 22 January and 5, 12, 19 February 2017; OpenStreetMap (ODbL)")
    save(fig, "agents-06-the-festival-crowd")


def fig07():
    F = load("festival.json")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(2, 2, left=0.05, right=0.97, top=0.82, bottom=0.09, wspace=0.3, hspace=0.42, width_ratios=[1.25, 1])
    ax = fig.add_subplot(gs[:, 0])
    along = F["along"]
    L = np.array([s["length"] for s in along])
    edges = np.concatenate([[0], np.cumsum(L)])
    rho = np.array([np.array(s["walkers"]) / (s["length"] * F["block_width_m"]) for s in along])  # walkers per m², minute by minute
    t = np.arange(rho.shape[1] + 1) / 60
    top = float(np.nanmax(rho))
    pc = ax.pcolormesh(edges, t, rho.T, cmap=SEA_MAP, vmin=0, vmax=max(0.7, 1.05 * top), shading="flat")
    ax.set_ylim(8, 24)
    ax.set_xlabel("along Little Bourke Street from Swanston Street (m)", labelpad=2)
    ax.set_ylabel("hour of Sunday 29 January")
    ax.set_title("the crowd's density along the street, minute by minute", loc="left")
    ax.grid(False)
    cb = fig.colorbar(pc, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("walkers per m², over a segment and a minute (line: Fruin's level C from 0.43)")
    cb.outline.set_visible(False)
    for lev, name in [(0.43, "Fruin C"), (0.72, "Fruin D")]:
        if lev <= max(0.7, 1.05 * top):
            cb.ax.axhline(lev, color=FLOW, lw=1.6)
    ax2 = fig.add_subplot(gs[0, 1])
    bw = np.array(F["block_walkers_by_minute"])
    ax2.plot(np.arange(len(bw)) / 60, bw, color=WALKERS, lw=1.4, label="walkers on the street")
    le = np.sum([np.array(s["at_level_e"]) for s in along], axis=0)
    ax2.plot(np.arange(len(le)) / 60, le, color=FLOW, lw=1.4, label="of them at Fruin's level E or worse")
    ax2.set_xlim(6, 24), ax2.set_xticks(range(6, 25, 3))
    ax2.set_ylabel("walkers")
    ax2.legend(loc="upper right", fontsize=10)
    ax2.set_title("how many, and in what crush", loc="left", fontsize=13)
    ax3 = fig.add_subplot(gs[1, 1])
    street = set(int(s["edge"]) for s in along)
    signals = sorted([q for q in F["queues"] if int(q["state"]) // 2 not in street], key=lambda q: -max(q["longest"]))[:2]
    own = sorted([q for q in F["queues"] if int(q["state"]) // 2 in street], key=lambda q: -max(q["longest"]))[:2]
    for q, col, lab in zip(signals + own, [PLUM, GREEN, FLOW, WALKERS],
                           ["at a signal, into Bourke Street", "at a signal, the next busiest", "into Little Bourke Street (its busiest mouth)", "into Little Bourke Street (the next)"]):
        ax3.plot(np.arange(24) + 0.5, q["longest"], color=col, lw=1.6, marker="o", ms=3.5, label=lab if "Bourke Street" not in lab or q["name"] == "Bourke Street" or "Little" in lab else f"at a signal, into {q['name']}")
    ax3.set_xlim(6, 24), ax3.set_xticks(range(6, 25, 3))
    ax3.set_xlabel("hour")
    ax3.set_ylabel("longest queue in the hour (walkers)")
    ax3.legend(loc="upper left", fontsize=9.5)
    ax3.set_title("queues: at the signals of the main streets, not on the festival's street", loc="left", fontsize=12.5)
    tt = F["totals"]
    title(fig, "Densities and queues of the festival's crowd",
          f"An ordinary Sunday's walkers and the fitted visitors together, the street given to walkers at {F['block_width_m']:.0f} m (an assumption): at most {int(bw.max())} on it, "
          f"{top:.2f} per m² at most over a segment and a minute;\nin the district {tt['los_e_hours']:,.0f} person-hours at Fruin's level E or worse, {tt['los_f_hours']:,.1f} at F")
    source(fig, SRC + "; Fruin 1971's levels of service; Weidmann 1993's speeds and jam density")
    save(fig, "agents-07-densities-and-queues")


def fig08():
    G = load("regime.json")
    cases = G["cases"]
    pols = ["open", "one-way", "closed", "POMCP"]
    colours = {"open": WALKERS, "one-way": PLUM, "closed": FLOW, "POMCP": GREEN}
    scales = G["scales"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(2, len(scales), left=0.08, right=0.97, top=0.8, bottom=0.13, wspace=0.25, hspace=0.42, height_ratios=[1.2, 1])
    for j, sc in enumerate(scales):
        ax = fig.add_subplot(gs[0, j])
        base = {}
        for p in pols:
            v = [c["value"] for c in cases if c["scale"] == sc and c["policy"] == p]
            base[p] = v
        o = np.array(base["open"])
        for i, p in enumerate(pols):
            d = np.array(base[p]) - o
            ax.bar(i, d.mean(), color=colours[p], width=0.62)
            ax.plot([i] * len(d), d, "o", color=INK, ms=4, zorder=4)
        ax.axhline(0, color=INK, lw=0.8)
        ax.grid(False), ax.grid(axis="y")
        ax.set_xticks(range(4)), ax.set_xticklabels(pols)
        ax.set_title(f"the crowd at {sc:g}× its fitted size", loc="left", fontsize=13)
        if j == 0:
            ax.set_ylabel("value against holding the street open\n(person-hours saved, crush weighted ×10)")
        ax2 = fig.add_subplot(gs[1, j])
        days = [c for c in cases if c["scale"] == sc and c["policy"] == "POMCP"]
        for k, c in enumerate(days):
            for q, a in enumerate(c["actions"]):
                t0 = G["start_h"] + q * G["epoch_s"] / 3600
                ax2.add_patch(plt.Rectangle((t0, k - 0.36), G["epoch_s"] / 3600, 0.72, color=colours[pols[int(a)]], lw=0))
        ax2.set_xlim(G["start_h"], G["start_h"] + G["epochs"] * G["epoch_s"] / 3600)
        ax2.set_ylim(-0.6, len(days) - 0.4)
        ax2.set_yticks(range(len(days)))
        ax2.set_yticklabels([f"day {c['day']:.0f}" + (f"\n{c['restarts']} restart" + ("s" if c["restarts"] > 1 else "") if c.get("restarts") else "") for c in days])
        ax2.set_xlabel("hour")
        ax2.grid(False)
        ax2.set_title("the regimes POMCP chose", loc="left", fontsize=12)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colours[p]) for p in pols[:3]]
    fig.legend(handles + [plt.Rectangle((0, 0), 1, 1, color=GREEN)], ["open", "one-way (from Swanston Street)", "closed to walkers passing through", "POMCP, deciding every half hour"],
               loc="upper left", ncol=4, frameon=False, bbox_to_anchor=(0.035, 0.895))
    title(fig, "The street's regime decided online",
          "Every half hour, POMCP plans on the walkers themselves, its belief about the crowd's size carried by particles and the counts at two counters; against holding one regime all day")
    source(fig, f"POMCP: {G['pomcp']['simulations']} simulations a decision, depth {G['pomcp']['depth']}, {G['pomcp']['particles']} particles, its walkers at {G['pomcp']['model_step_s']:g} s steps. "
                f"Reward: −(person-hours) − 10 × (person-hours at Fruin's level F). Common random numbers across the policies of a day.")
    save(fig, "agents-08-the-regime-decided")


def fig09():
    F = load("festival.json")
    rec = F["record"]
    ev = rec["list"]
    mjd0 = min(e["from_mjd"] for e in ev)
    day0 = math.floor(mjd0 - 13 / 24) + 13 / 24   # Sunday 00:00 AEDT
    rows, names = {}, {f"street:{s['edge']}": s["name"] or s["kind"] for s in SEG}
    for e in ev:
        for c in e["couplings"]:
            key, _, name = c.partition("|")
            rows.setdefault(key, []).append(e)
            if name:
                names[key] = name
    dur = {k: sum(e["to_mjd"] - e["from_mjd"] for e in rows[k]) for k in rows}
    crowd_keys = sorted([k for k in rows if any(e["kind"] == "crowd" for e in rows[k])], key=lambda k: -dur[k])[:10]
    queue_keys = sorted([k for k in rows if k not in crowd_keys], key=lambda k: -dur[k])[:14]
    keys = crowd_keys + queue_keys
    fig, ax = figure()
    fig.subplots_adjust(left=0.24, right=0.97, top=0.84, bottom=0.16)
    col = {"crowd": FLOW, "queue": PLUM, "open": WALKERS, "closed": FLOW, "one-way": PLUM}
    for i, k in enumerate(keys):
        for e in rows[k]:
            a, b = (e["from_mjd"] - day0) * 24, (e["to_mjd"] - day0) * 24
            if e["kind"] == "crowd":
                ax.add_patch(plt.Rectangle((a, i - 0.32), max(b - a, 0.05), 0.64, color=FLOW, lw=0))
            elif e["kind"] == "queue":
                ax.add_patch(plt.Rectangle((a, i - 0.18), max(b - a, 0.05), 0.36, color=PLUM, lw=0, alpha=0.85))
            else:
                ax.add_patch(plt.Rectangle((a, i - 0.42), max(b - a, 0.05), 0.84, fill=False, ec=WALKERS, lw=1.2))
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([("the corner of " if k.startswith("junction:") else "") + names.get(k, k) for k in keys], fontsize=9.5)
    ax.set_xlim(6, 24), ax.set_xticks(range(6, 25, 2))
    ax.set_ylim(-0.8, len(keys) - 0.2)
    ax.invert_yaxis()
    ax.set_xlabel("hour of Sunday 29 January 2017 (simulated)", labelpad=2)
    ax.grid(axis="y", visible=False)
    ax.add_patch(plt.Rectangle((0, 0), 0, 0, color=FLOW, label="crowd: ten walkers or more at Fruin's level E, minute after minute"))
    ax.add_patch(plt.Rectangle((0, 0), 0, 0, color=PLUM, label="queue: twenty or more waiting at a street's mouth in the hour"))
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.07), ncol=2, fontsize=10.5)
    k = rec["kinds"]
    title(fig, "The festival day as a record of simulated events",
          f"Written in urban science's record (events/record.hpp): {rec['events']:.0f} events — {k.get('crowd', 0)} crowds, {k.get('queue', 0)} queues, "
          f"{k.get('count', 0)} hourly counts — every one with provenance \"simulated\"; the record is {'sound' if rec['sound'] else 'NOT sound'}")
    source(fig, "urban-science movement, written as urban::events::Record; never drawn as observed: the workbench tells a simulated record by its provenance. Rows: the street spaces with crowds, then the junctions with the longest queues")
    save(fig, "agents-09-the-simulated-record")


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
    for w in which:
        name = f"fig0{w}"
        try:
            globals()[name]()
            print("drew", name)
        except FileNotFoundError as e:
            print("skipped", name, "(", e, ")")
