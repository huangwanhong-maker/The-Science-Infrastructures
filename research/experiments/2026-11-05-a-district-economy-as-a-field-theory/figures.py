"""Figures of the case study "A district's economy as a field theory" (README §9), in the house style of the search
case studies (research/experiments/style/search_style.py): search-80 ... search-89, 1600 x 1000, into assets/gallery/.

Reads results/ (calibration.json, field.json, search.json; field-geometry.json, the ODbL geometry, not committed and
made again by field_theory), and the shared problem's results (problem.json, sensors.json, district-geometry.json).

    .venv/Scripts/python.exe research/experiments/2026-11-05-a-district-economy-as-a-field-theory/figures.py
"""

from __future__ import annotations

import json
import os
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import (ALGORITHM, DUSK, HAIRLINE, INK, MUTED, ORCHID, PANEL, PAPER, figure, save, source, title)  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RES = pathlib.Path(os.environ.get("FIELD_RESULTS", str(HERE / "results")))
PROB = HERE.parent / "2026-11-05-a-district-economy-the-problem" / "results"


def load(path: pathlib.Path):
    return json.loads(path.read_text(encoding="utf-8"))


CAL = load(RES / "calibration.json")
FLD = load(RES / "field.json")
SRCH = load(RES / "search.json")
GEO = load(RES / "field-geometry.json")
PG = load(PROB / "district-geometry.json")
PROBLEM = load(PROB / "problem.json")
SENSORS = load(PROB / "sensors.json")

KINDS = ["vacant", "retail", "food and drink", "services", "office lobby", "culture"]
KIND_COLOUR = ["#DCD3D9", "#5C88B4", "#8E3B78", "#7FA37A", "#C98B6B", "#D49A3A"]
CREDIT = "Data: OpenStreetMap contributors (ODbL 1.0); City of Melbourne, CLUE and Pedestrian Counting System (CC BY 4.0). A study of methods on a real place, not a planning proposal."
DAYS = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
WEEK = 168

# ---- the frame: the grid's streets square to the page -----------------------------------------------------------------
Q = np.array(PG["quad"])
_d = Q[1] - Q[0]
THETA = -math.atan2(_d[1], _d[0])
CENTRE = Q.mean(axis=0)
_R = np.array([[math.cos(THETA), -math.sin(THETA)], [math.sin(THETA), math.cos(THETA)]])


def rot(p) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    return (p - CENTRE) @ _R.T


VERT = rot(np.array(GEO["vertices"]))
EDGES = GEO["edges"]
LINES = [rot(np.array(e["line"])) for e in EDGES]
ENDS = np.array([e["ends"] for e in EDGES])
CLS = np.array([e["class"] for e in EDGES])
FEET = rot(np.array(GEO["frontage_feet"]))
FRONT = PROBLEM["frontages"]
NF = len(FRONT)
FOOT = [rot(np.array(f["footprint"])) for f in PG["frontages"]]
CONTEXT = [rot(np.array(c["footprint"])) for c in PG["context"]]
COUNTERS = {c["id"]: c for c in GEO["counters"]}
LINKS = [(rot(l["a"]), rot(l["b"]), l["kind"]) for l in GEO["links"]]
XLIM = (VERT[:, 0].min() - 15, VERT[:, 0].max() + 15)
YLIM = (VERT[:, 1].min() - 15, VERT[:, 1].max() + 15)
_QR = rot(Q)
AREA_X = (_QR[:, 0].min() - 45, _QR[:, 0].max() + 45)
AREA_Y = (_QR[:, 1].min() - 45, _QR[:, 1].max() + 45)
THOUSANDS = FuncFormatter(lambda x, _: f"{x / 1000:,.0f}k")


def base_map(ax, context=True, frontages=True, alpha=1.0, area=False):
    ax.set_facecolor(PAPER)
    if context:
        ax.add_collection(PolyCollection(CONTEXT, facecolor="#EFE6E1", edgecolor="#E2D6D0", linewidth=0.3, zorder=0))
    if frontages:
        ax.add_collection(PolyCollection(FOOT, facecolor="#E7DBD6", edgecolor="#D8C9C3", linewidth=0.4, zorder=0.5, alpha=alpha))
    ax.set_xlim(*(AREA_X if area else XLIM)), ax.set_ylim(*(AREA_Y if area else YLIM))
    ax.set_aspect("equal")
    ax.set_xticks([]), ax.set_yticks([])
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


def network(ax, values=None, cmap=ORCHID, norm=None, widths=None, colour=HAIRLINE, z=2):
    if values is None:
        lc = LineCollection(LINES, colors=colour, linewidths=0.8 if widths is None else widths, zorder=z, capstyle="round")
    else:
        lc = LineCollection(LINES, array=np.asarray(values), cmap=cmap, norm=norm, linewidths=widths if widths is not None else 2.0, zorder=z, capstyle="round")
    ax.add_collection(lc)
    return lc


def edge_values(vertex_values):
    v = np.asarray(vertex_values)
    return 0.5 * (v[ENDS[:, 0]] + v[ENDS[:, 1]])


def scale_bar(ax, length=100, area=False):
    x0, y0 = ((AREA_X[0] + 12, AREA_Y[0] + 10) if area else (XLIM[0] + 20, YLIM[0] + 18))
    ax.plot([x0, x0 + length], [y0, y0], color=INK, lw=1.2, zorder=10)
    ax.text(x0 + length / 2, y0 + 8, f"{length} m", ha="center", va="bottom", fontsize=9, color=MUTED, zorder=10)


def hours_axis(ax, label=True):
    ax.set_xlim(0, WEEK)
    ax.set_xticks([24 * d + 12 for d in range(7)])
    ax.set_xticklabels(DAYS)
    for d in range(1, 7):
        ax.axvline(24 * d, color=HAIRLINE, lw=0.6, zorder=0)
    ax.grid(False, axis="x")


def plan_kinds(plan):
    return [int(k) for k in plan["kind"]]


def open_mask(plan, f):
    k = plan["kind"][f]
    if k == 0:
        return np.zeros(WEEK, bool)
    h = plan["hours"][4 * f: 4 * f + 4]
    out = np.zeros(WEEK, bool)
    for t in range(WEEK):
        wd = (t // 24) < 5
        s, l = (h[0], h[1]) if wd else (h[2], h[3])
        tt = t % 24
        out[t] = l > 0 and s <= tt < min(24, s + l)
    return out


PLANS = {p["name"]: p for p in SRCH["plans"]}
AS_BUILT = PLANS["as built"]
BEST_HB = PLANS["best heat bath"]
BEST_BB = next(p for p in SRCH["plans"] if p["name"].startswith("best black box"))


# ---- 80: the people as a field ------------------------------------------------------------------------------------------
def fig80():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The people of the district as a field", "Passers an hour on a footpath: the stationary path of the action, calibrated to the City's counters by Gaussian conditioning line by line")
    ax = fig.add_axes([0.02, 0.05, 0.60, 0.84])
    base_map(ax)
    vals = edge_values(FLD["vertex_weekly_mean"])
    norm = LogNorm(vmin=30, vmax=3000)
    w = 1.0 + 2.8 * np.clip(np.log10(np.maximum(vals, 30) / 30) / 2, 0, 1)
    lc = network(ax, np.maximum(vals, 30), ORCHID, norm, w)
    sens = {s["id"]: s for s in SENSORS["sensors"]}
    for cid, c in COUNTERS.items():
        p = rot(c["at"])
        mw = np.mean(sens[cid]["mean_week"])
        ax.scatter([p[0]], [p[1]], s=95, marker="s", c=[ORCHID(norm(max(mw, 30)))], edgecolors=INK, linewidths=1.0, zorder=6)
    ax.scatter(FEET[:, 0], FEET[:, 1], s=10, color=INK, zorder=7, linewidths=0)
    scale_bar(ax)
    cax = fig.add_axes([0.20, 0.072, 0.18, 0.012])
    cb = fig.colorbar(lc, cax=cax, orientation="horizontal")
    cb.set_label("passers an hour, weekly mean", color=INK, fontsize=10.5)
    cb.ax.xaxis.set_label_position("top")
    cb.outline.set_visible(False)
    ax.legend(handles=[Line2D([], [], marker="s", ls="", mfc=ORCHID(0.6), mec=INK, ms=9, label="a counter, coloured by its observed mean"),
                       Line2D([], [], marker="o", ls="", color=INK, ms=4, label="a frontage's foot on its street")], loc="upper left", fontsize=10.5)
    ax2 = fig.add_axes([0.69, 0.52, 0.28, 0.34])
    obs = np.array([np.mean(sens[c["id"]]["mean_week"]) * WEEK for c in CAL["loo"]["counters"]])
    loo = np.array([sum(c["predicted"]) for c in CAL["loo"]["counters"]])
    nn = np.array([sum(c["predicted"]) for c in CAL["nearest_counter"]["counters"]])
    ax2.loglog([1e4, 1e6], [1e4, 1e6], color=HAIRLINE, lw=1)
    ax2.scatter(obs, nn, s=34, color=ALGORITHM["random"], label="the nearest other counter", zorder=3, alpha=0.85)
    ax2.scatter(obs, np.maximum(loo, 1.6e4), s=46, color="#8E3B78", label="the field, its site left out", zorder=4)
    ax2.set_xlim(3e4, 4e5), ax2.set_ylim(1.5e4, 8e5)
    for a in (ax2.xaxis, ax2.yaxis):
        a.set_major_formatter(THOUSANDS)
        a.set_minor_formatter(FuncFormatter(lambda x, _: ""))
    ax2.set_xticks([3e4, 1e5, 3e5]), ax2.set_yticks([3e4, 1e5, 3e5])
    ax2.set_xlabel("observed weekly passers")
    ax2.set_ylabel("predicted")
    ax2.set_title("leave one site out: the 17 counters on the network", loc="left", fontsize=13)
    ax2.legend(loc="upper left", fontsize=10)
    ch = CAL["chosen"]
    ax3 = fig.add_axes([0.69, 0.06, 0.28, 0.38])
    ax3.axis("off")
    rows = [
        ("screening length ξ = √(D/μ)", f"{ch['length_m']:.0f} m, D = {ch['diffusivity_m2_per_h']/3600:.1f} m²/s"),
        ("leaving rate μ (assumed)", f"{ch['leave_per_hour']:.0f} an hour"),
        ("ways' conductances", f"little streets ×{ch['little']:g}, lanes ×{ch['lane']:g}"),
        ("nugget τ²", f"{ch['nugget']:g}, relative to the level"),
        ("source groups", f"{len(ch['groups'])}: edges, buildings, ways"),
        ("error left out (normalised)", f"{CAL['loo']['normalised_sse']:.3f}   nearest counter {CAL['nearest_counter']['normalised_sse']:.2f}"),
        ("median R² of a week, left out", f"{CAL['loo']['median_r2']:.2f}   nearest counter {CAL['nearest_counter']['median_r2']:.2f}"),
        ("weekly totals, RMS log error", f"{CAL['loo']['log_total_rms']:.2f}   nearest counter {CAL['nearest_counter']['log_total_rms']:.2f}"),
        ("in sample, median R²", f"{CAL['in_sample']['median_r2']:.3f}"),
    ]
    for i, (k, v) in enumerate(rows):
        y = 1 - i * 0.111
        ax3.text(0.0, y, k, fontsize=10.5, color=MUTED, va="top", transform=ax3.transAxes)
        ax3.text(0.0, y - 0.047, v, fontsize=11.5, color=INK, va="top", transform=ax3.transAxes)
    source(fig, CREDIT)
    save(fig, "search-80-the-people-as-a-field")

# ---- 81: the modes ------------------------------------------------------------------------------------------------------
def fig81():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The district's modes: Hodge Laplacians of every degree", "Δ₀'s first modes above its kernel on the lattice, and the spectrum of Δ₁ split into gradients (Δ₀'s), curls (Δ₂'s) and harmonic currents")
    rates = np.array(FLD["mode_rates_per_hour"])
    b0, b1, b2 = FLD["hodge"]["betti"]
    modes = np.array(FLD["modes"])
    mu = CAL["chosen"]["leave_per_hour"]
    Dm = CAL["chosen"]["diffusivity_m2_per_h"]
    for i in range(6):
        k = min(b0 + i, len(modes) - 1)
        ax = fig.add_axes([0.015 + 0.205 * (i % 3), 0.47 - 0.43 * (i // 3), 0.20, 0.40])
        base_map(ax, context=False, frontages=False)
        network(ax, None, colour="#EDE4E8", widths=2.6, z=1)
        v = edge_values(modes[k])
        m = np.max(np.abs(v)) or 1
        network(ax, v, DUSK, TwoSlopeNorm(0, -m, m), 2.2)
        lam = rates[k] / Dm
        wl = 2 * math.pi / math.sqrt(lam) if lam > 0 else float("inf")
        ax.set_title(f"mode {k}: Dλ = {rates[k]:.2f} an hour; 2π/√λ = {wl:,.0f} m", fontsize=11, loc="left")
    ax = fig.add_axes([0.665, 0.11, 0.315, 0.74])
    d0 = np.array(FLD["hodge"]["delta0"])
    d1 = np.array(FLD["hodge"]["delta1"])
    d2 = np.array(FLD["hodge"]["delta2"])
    top = max(d1.max(), d0.max())
    nz0 = np.sort(d0[d0 > 1e-9 * top])
    nz2 = np.sort(d2[d2 > 1e-9 * max(d2.max(), 1e-300)])
    kind = []
    for x in d1:
        if x <= 1e-9 * top:
            kind.append(2)
        else:
            e0 = np.min(np.abs(nz0 - x)) / x if len(nz0) else 1
            e2 = np.min(np.abs(nz2 - x)) / x if len(nz2) else 1
            kind.append(0 if e0 <= e2 else 1)
    kind = np.array(kind)
    rank = np.arange(len(d1))
    floor = 1e-3
    cols = {0: "#8E3B78", 1: "#D49A3A", 2: "#4E8F8B"}
    labs = {0: f"gradients: Δ₀'s nonzero modes ({(kind == 0).sum()})", 1: f"curls round the blocks: Δ₂'s ({(kind == 1).sum()})", 2: f"harmonic, round the holes: b₁ = {(kind == 2).sum()}"}
    for c in (0, 1, 2):
        sel = kind == c
        ax.scatter(rank[sel], np.maximum(d1[sel] * Dm, floor), s=7 if c == 0 else 22, color=cols[c], label=labs[c], zorder=3 if c else 2)
    ax.axhline(mu, color=MUTED, lw=0.8, ls="--")
    ax.text(len(d1) * 0.98, mu * 0.55, "μ: below it, leaving screens the mode", fontsize=9.5, color=MUTED, ha="right", va="top")
    ax.set_yscale("log")
    ax.set_ylim(floor * 0.7, 3e5)
    ax.set_xlabel("rank of Δ₁'s eigenvalue")
    ax.set_ylabel("D × eigenvalue (an hour); zeros drawn at 10$^{-3}$")
    ax.legend(loc="upper left", fontsize=10)
    ax.set_title(f"Δ₁ on the {len(d1):,} pieces: b₀ = {b0}, b₁ = {b1}, b₂ = {b2}", loc="left", fontsize=13)
    source(fig, CREDIT.replace(" A study of methods on a real place, not a planning proposal.", "") + " physis::discretisation::hodge_laplacian on urban::landscape's district, subdivided.")
    save(fig, "search-81-the-districts-modes")

# ---- 82: space, time and frequency ----------------------------------------------------------------------------------------
def fig82():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Space, time and frequency: the field in its modes", "The power of the calibrated field by spatial mode and line of the week, the action's filter, and the weeks it makes at four frontages")
    bins = FLD["spectrum"]["bins"]
    P = np.array([b["power"] for b in bins])
    ax = fig.add_axes([0.06, 0.52, 0.40, 0.34])
    im = ax.imshow(np.log10(np.maximum(P[:, :50] / P.max(), 1e-9)), aspect="auto", origin="lower", cmap=ORCHID, extent=[-0.5, 49.5, -0.5, len(bins) - 0.5], vmin=-7, vmax=0)
    labels = []
    for b in bins:
        lo, hi = b["modes"]
        labels.append("kernel" if lo == 0 else f"{lo}–{hi}")
    ax.set_yticks(range(len(bins)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xticks([0, 7, 14, 21, 28, 35, 42, 49])
    ax.set_xticklabels(["W0", "D1", "D2", "D3", "D4", "D5", "D6", "D7"], fontsize=9)
    ax.set_xlabel("line of the week, m (cycles a week; D_j = 7j, a day's harmonics)")
    ax.set_ylabel("Δ₀'s modes, by octave")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
    cb.set_label("log₁₀ share of the largest", color=INK)
    cb.outline.set_visible(False)
    ax.set_title("the field's power |c_k(m)|²", loc="left", fontsize=13)
    # The filter |t_k(m)| = 1/|iω + Dλ + μ|, relative to the static.
    ax2 = fig.add_axes([0.56, 0.52, 0.40, 0.34])
    mu = CAL["chosen"]["leave_per_hour"]
    rates = np.array(FLD["mode_rates_per_hour"])
    m = np.arange(85)
    om = 2 * np.pi * m / WEEK
    for b, col in zip([1, 3, 5, 7], ["#8E3B78", "#4E8F8B", "#D49A3A", "#6C8EBF"]):
        lo, hi = bins[b]["modes"]
        r = rates[(lo + hi) // 2]
        t = (r + mu) / np.sqrt(om ** 2 + (r + mu) ** 2)
        ph = np.degrees(np.arctan2(om, r + mu))
        ax2.plot(m, t, color=col, lw=1.8, label=f"modes {lo}–{hi}, Dλ + μ = {r + mu:.0f}/h")
    ax2.set_ylim(0.6, 1.02)
    ax2.set_xlabel("line of the week, m")
    ax2.set_ylabel("|t_k(m)| / t_k(0)")
    ax2.set_title("the action's filter: a lag of 1/(Dλ + μ), at most 10 minutes", loc="left", fontsize=13)
    ax2.legend(loc="lower left", fontsize=10)
    # The frontages' weeks.
    ax3 = fig.add_axes([0.06, 0.08, 0.62, 0.34])
    P0 = np.array(FLD["frontage_passers"])
    tot = P0.sum(1)
    pick = [int(np.argmax(tot)), int(np.argsort(tot)[len(tot) * 3 // 4]), int(np.argsort(tot)[len(tot) // 2]), int(np.argsort(tot)[len(tot) // 6])]
    for f, col in zip(pick, ["#8E3B78", "#4E8F8B", "#D49A3A", "#6C8EBF"]):
        nm = FRONT[f]["name"] if FRONT[f]["name"] and not FRONT[f]["name"].startswith("way/") else f"frontage {f}"
        ax3.plot(np.arange(WEEK) + 0.5, P0[f], color=col, lw=1.6, label=f"{nm}, {FRONT[f]['street']}: {P0[f].sum():,.0f} a week")
    hours_axis(ax3)
    ax3.set_ylabel("passers an hour")
    ax3.set_title("the stationary path at four frontages, as built", loc="left", fontsize=13)
    ax3.legend(loc="upper left", fontsize=9, ncol=2)
    ax4 = fig.add_axes([0.74, 0.08, 0.23, 0.34])
    bl = np.array(FLD["spectrum"]["by_line"])
    share = bl / bl.sum()
    cols = ["#8E3B78" if (i % 7 == 0 and i > 0) else ("#3B2F3A" if i == 0 else "#6C8EBF") for i in range(85)]
    ax4.bar(np.arange(85), np.maximum(share, 1e-7), color=cols, width=0.8)
    ax4.set_yscale("log")
    ax4.set_ylim(1e-6, 1)
    ax4.set_xlabel("line, m")
    ax4.set_ylabel("share of the power")
    ax4.set_title("its lines: mean, day, week", loc="left", fontsize=13)
    source(fig, CREDIT)
    save(fig, "search-82-space-time-and-frequency")


# ---- 83: people flow -----------------------------------------------------------------------------------------------------
def fig83():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "People flow: the conserved number and the currents", "Noether's charge of the action, the people on the district's streets; Fick's current on the network; the counters' directions against it")
    N = np.array(FLD["noether"]["people_on_the_streets"])
    S = np.array(FLD["noether"]["inflow_per_hour"])
    ax = fig.add_axes([0.06, 0.58, 0.40, 0.28])
    ax.plot(np.arange(WEEK) + 0.5, N, color="#8E3B78", lw=2.2, label="people on the streets, N(t)")
    ax.plot(np.arange(WEEK) + 0.5, S / CAL["chosen"]["leave_per_hour"], color="#F6EEE8", lw=0.9, ls=(0, (2, 2)), label="inflow / μ: the static limit, indistinguishable")
    hours_axis(ax)
    ax.set_ylabel("people")
    ax.legend(loc="upper left", fontsize=10)
    ax.set_title(f"dN/dt = Σ s − μN, at every line to {FLD['noether']['largest_relative_residual']:.0e}", loc="left", fontsize=13)
    ax2 = fig.add_axes([0.50, 0.06, 0.48, 0.84])
    base_map(ax2)
    vals = edge_values(FLD["vertex_hours"]["Saturday 13:00"])
    network(ax2, np.maximum(vals, 30), ORCHID, LogNorm(30, 3000), 1.6)
    J = np.array(FLD["current_saturday_13"])
    mids, dirs = [], []
    for e, l in enumerate(LINES):
        a, b = VERT[ENDS[e, 0]], VERT[ENDS[e, 1]]
        n = np.linalg.norm(b - a)
        if n < 5:
            continue
        mids.append(0.5 * (a + b))
        dirs.append((b - a) / n * J[e])
    mids, dirs = np.array(mids), np.array(dirs)
    mag = np.linalg.norm(dirs, axis=1)
    keep = mag > np.percentile(mag, 75)
    keep &= (np.arange(len(mag)) % 2 == 0)
    ax2.quiver(mids[keep, 0], mids[keep, 1], dirs[keep, 0], dirs[keep, 1], color=INK, angles="xy", scale_units="xy",
               scale=np.percentile(mag[keep], 90) / 28, width=0.0024, zorder=6, headwidth=4.5, headlength=5, alpha=0.85)
    scale_bar(ax2)
    ax2.set_title(f"Fick's current at Saturday 13:00 over the field (the quarter strongest; the largest {np.max(np.abs(J)):,.0f} an hour)", loc="left", fontsize=12)
    dl = FLD["directions"]["counters"]
    best = sorted(dl, key=lambda d: -d["rms_observed"])[:3]
    ax3 = fig.add_axes([0.06, 0.08, 0.40, 0.38])
    for d, col in zip(best, ["#8E3B78", "#4E8F8B", "#D49A3A"]):
        ax3.plot(np.arange(WEEK) + 0.5, d["observed"], color=col, lw=1.5, label=f"{d['name'][:32]} (along {d['along']})")
        ax3.plot(np.arange(WEEK) + 0.5, d["model"], color=col, lw=1.6, ls=":")
    hours_axis(ax3)
    ax3.axhline(0, color=MUTED, lw=0.6)
    ax3.set_ylabel("net walkers an hour on the footpath")
    ax3.set_ylim(-1050, 1150)
    ax3.legend(loc="upper left", fontsize=9, ncol=2, title="observed (solid), and Fick's current on that footpath (dotted)", title_fontsize=9)
    cors = [d["correlation"] for d in dl]
    ax3.set_title(f"the counters' directions: mean correlation with Fick's current {np.mean(cors):.2f}", loc="left", fontsize=13)
    source(fig, CREDIT)
    save(fig, "search-83-people-flow")

# ---- 84: calibration ---------------------------------------------------------------------------------------------------
def fig84():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Calibration to the counters, line by line", "Each counter's mean week (ink), the field's prediction there with its site left out (orchid), and the nearest other counter's week (grey)")
    sens = {s["id"]: s for s in SENSORS["sensors"]}
    cs = CAL["loo"]["counters"]
    nn = {c["id"]: c for c in CAL["nearest_counter"]["counters"]}

    def short(name):
        name = name.replace(" St ", " St. ").replace("Street", "St")
        return name if len(name) <= 30 else name[:29].rsplit(" ", 1)[0] + "…"

    for i, c in enumerate(cs):
        ax = fig.add_axes([0.045 + 0.158 * (i % 6), 0.635 - 0.265 * (i // 6), 0.13, 0.19])
        y = np.array(sens[c["id"]]["mean_week"])
        ax.plot(np.arange(WEEK), nn[c["id"]]["predicted"], color=ALGORITHM["random"], lw=0.8, alpha=0.8)
        ax.plot(np.arange(WEEK), y, color=INK, lw=1.0)
        ax.plot(np.arange(WEEK), c["predicted"], color="#8E3B78", lw=1.2)
        ax.set_xlim(0, WEEK)
        ax.set_xticks([])
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))
        ax.tick_params(labelsize=7.5, pad=1)
        ax.set_title(f"{short(c['name'])}\nR² {c['r2']:.2f}; total ×{c['total_ratio']:.2f}", fontsize=9, loc="left")
    ax = fig.add_axes([0.045 + 0.158 * 5, 0.635 - 0.265 * 2, 0.13, 0.19])
    grid = CAL["grid"]
    xs = sorted({g["length_m"] for g in grid})
    for variant, col in [((False, False), "#8E3B78"), ((False, True), "#4E8F8B"), ((True, False), "#D49A3A"), ((True, True), "#6C8EBF")]:
        ys = [min(g["loo_nsse"] for g in grid if g["length_m"] == x and (g["edge_from_counters"], g["trams"]) == variant) for x in xs]
        ax.plot(xs, ys, color=col, lw=1.4, marker="o", ms=3)
    ax.set_xscale("log")
    ax.set_xticks([30, 100, 300])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.0f}"))
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda x, _: ""))
    ax.set_xlabel("ξ (m)", fontsize=9)
    ax.tick_params(labelsize=7.5)
    ax.set_title("the grid: error left out,\nbest over the rest", fontsize=9, loc="left")
    fig.legend(handles=[Line2D([], [], color=INK, lw=1.2, label="observed mean week"), Line2D([], [], color="#8E3B78", lw=1.4, label="the field, its site left out"),
                        Line2D([], [], color=ALGORITHM["random"], lw=1, label="the nearest other counter"),
                        Line2D([], [], color="#8E3B78", marker="o", label="grid: edges by width"), Line2D([], [], color="#4E8F8B", marker="o", label="… and tram streets"),
                        Line2D([], [], color="#D49A3A", marker="o", label="edges fed by outside counters"), Line2D([], [], color="#6C8EBF", marker="o", label="… and tram streets")],
               loc="lower center", ncol=7, fontsize=9.5, bbox_to_anchor=(0.44, 0.04))
    source(fig, CREDIT)
    save(fig, "search-84-calibration")

# ---- 85: the Euclidean measure ------------------------------------------------------------------------------------------
def fig85():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    fl = SRCH["fluctuations"]
    title(fig, "The Euclidean measure: fluctuations and uncertainty", "How a real week departs from the mean week; how far the counters' local residuals are correlated, against the action; and what the counters leave unknown")
    ax = fig.add_axes([0.06, 0.57, 0.40, 0.29])
    sd = np.array(fl["common_sd_by_hour"])
    ids = set(c["id"] for c in GEO["counters"])
    tot = np.sum([s["mean_week"] for s in SENSORS["sensors"] if s["id"] in ids], axis=0)
    busy = tot > 0.15 * tot.max()
    x = np.arange(WEEK) + 0.5
    ax.plot(x, np.where(busy, 100 * sd, np.nan), color="#8E3B78", lw=1.8)
    ax.plot(x, np.where(~busy, np.minimum(100 * sd, 100), np.nan), color="#CA8DBE", lw=0.9, ls=":")
    hours_axis(ax)
    ax.set_ylim(0, 60)
    ax.set_ylabel("% of the mean week")
    ax.set_title("the district's common deviation: its s.d. over 52 weeks (dotted: the night's few walkers)", loc="left", fontsize=12)
    ax2 = fig.add_axes([0.06, 0.09, 0.40, 0.36])
    pr = np.array([p[3] for p in fl["pairs"]])
    ob = np.array([p[4] for p in fl["pairs"]])
    dist = np.array([p[2] for p in fl["pairs"]])
    ax2.scatter(np.maximum(dist, 3), ob, s=26, color="#4E8F8B", label="observed, the local residual", zorder=3)
    ax2.scatter(np.maximum(dist, 3), pr, s=22, color="#8E3B78", marker="D", label="the action's, static line", zorder=4)
    ax2.axhline(0, color=MUTED, lw=0.6)
    ax2.set_xscale("log")
    ax2.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax2.set_xlabel("distance between the two counters (m)")
    ax2.set_ylabel("correlation")
    ax2.legend(loc="upper right", fontsize=10)
    ax2.set_title(f"136 pairs of counters: the two agree at r = {fl['correlation_agreement']:.2f}", loc="left", fontsize=13)
    ax3 = fig.add_axes([0.51, 0.06, 0.47, 0.84])
    base_map(ax3)
    network(ax3, None, colour="#D9CCD3", widths=1.0)
    mean = np.array(fl["frontage_weekly_total"])
    sdp = np.array(fl["frontage_weekly_total_posterior_sd"])
    rel = sdp / np.maximum(mean, 1)
    norm = Normalize(0, 1.0)
    sc = ax3.scatter(FEET[:, 0], FEET[:, 1], s=25 + 150 * np.sqrt(mean / mean.max()), c=np.minimum(rel, 1), cmap=ORCHID, norm=norm, edgecolors=INK, linewidths=0.5, zorder=6)
    for cid, c in COUNTERS.items():
        p = rot(c["at"])
        ax3.scatter([p[0]], [p[1]], s=55, marker="s", facecolor="none", edgecolor=INK, linewidths=1.0, zorder=7)
    cax = fig.add_axes([0.75, 0.105, 0.20, 0.014])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label("posterior s.d. of weekly passers / their mean", color=INK, fontsize=10.5)
    cb.outline.set_visible(False)
    ax3.set_title(f"what the counters (squares) leave unknown at each frontage; median {np.median(rel):.2f} of the mean", loc="left", fontsize=12)
    scale_bar(ax3)
    source(fig, CREDIT + f" T_ζ = {fl['T_zeta']:.4f}.")
    save(fig, "search-85-the-euclidean-measure")

# ---- 86: the searches ----------------------------------------------------------------------------------------------------
FAMILY = {"GA": "GA", "GA from as built": "GA", "CMA-ES": "CMA-ES", "SA": "SA", "UMDA": "PBIL", "random": "random",
          "heat bath (from as built)": "ILS", "heat bath (from a random plan)": "BO"}
ORDER = ["heat bath (from as built)", "heat bath (from a random plan)", "UMDA", "SA", "GA from as built", "GA", "CMA-ES", "random"]
COLOURS = {"heat bath (from as built)": "#8E3B78", "heat bath (from a random plan)": "#D9738B", "UMDA": "#4E8F8B", "SA": "#5C7FA3", "GA from as built": "#7FA37A",
           "GA": "#A6C29F", "CMA-ES": "#D49A3A", "random": "#9A9097"}


def fig86():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    budget = SRCH["budget_evaluations"]
    title(fig, "What the action's structure buys: searches at equal budgets", f"Five seeds each, {budget:,.0f} evaluations; the heat bath's local work counted at the measured cost of one evaluation ({SRCH['evaluation_ms']:.2f} ms)")
    ax = fig.add_axes([0.07, 0.10, 0.55, 0.76])
    v0 = AS_BUILT["outcome"]["value_aud_week"]
    for name in ORDER:
        rs = [r for r in SRCH["runs"] if r["algorithm"] == name]
        for i, r in enumerate(rs):
            c = np.array(r["curve"])
            if len(c) == 0:
                continue
            x = np.concatenate([c[:, 0], [max(c[-1, 0], r["evaluations"])]])
            y = np.concatenate([c[:, 1], [c[-1, 1]]])
            ax.step(x, y / 1000, where="post", color=COLOURS[name], lw=1.6 if i == 0 else 1.0, alpha=1 if i == 0 else 0.55, label=name if i == 0 else None)
    ax.axhline(v0 / 1000, color=INK, lw=0.8, ls="--")
    ax.text(1.2, v0 / 1000 + 12, "the plan as built", fontsize=10, color=INK)
    bb = max(r["value"] for r in SRCH["runs"] if not r["algorithm"].startswith("heat"))
    ax.axhline(bb / 1000, color="#4E8F8B", lw=0.8, ls=":")
    hits = []
    for r in SRCH["runs"]:
        if r["algorithm"].startswith("heat"):
            hits.append(next((x for x, y in r["curve"] if y >= bb), np.nan))
    ax.text(1.2, bb / 1000 + 12, f"the best of the 30 black-box runs, A$ {bb:,.0f}: the heat bath passes it after {np.nanmin(hits):,.0f}–{np.nanmax(hits):,.0f} evaluations' worth", fontsize=10, color="#4E8F8B")
    ax.set_xscale("log")
    ax.set_xlim(1, budget * 1.1)
    ax.set_xlabel("evaluations (the heat bath: their equivalents)")
    ax.set_ylabel("best value found (A$ thousand a week)")
    ax.legend(loc="lower right", fontsize=10)
    ax.set_title("anytime: the best value so far, each seed", loc="left", fontsize=13)
    ax2 = fig.add_axes([0.68, 0.10, 0.30, 0.76])
    for i, name in enumerate(ORDER):
        rs = [r for r in SRCH["runs"] if r["algorithm"] == name]
        vals = np.array([r["value"] for r in rs]) / 1000
        ev = np.mean([r["evaluations"] for r in rs])
        sec = np.mean([r["seconds"] for r in rs])
        y = len(ORDER) - 1 - i
        ax2.scatter(vals, np.full(len(vals), y), color=COLOURS[name], s=40, zorder=3)
        ax2.plot([vals.min(), vals.max()], [y, y], color=COLOURS[name], lw=1.2)
        ax2.text(ax2.get_xlim()[0] if False else -95, y + 0.28, f"{name}: {vals.mean():.0f}k ± {vals.std(ddof=1):.1f}k; {ev:,.0f} ev., {sec:.1f} s", fontsize=9.5, color=INK)
    ax2.axvline(v0 / 1000, color=INK, lw=0.8, ls="--")
    ax2.set_yticks([])
    ax2.set_xlim(-100, 760)
    ax2.set_ylim(-0.6, len(ORDER))
    ax2.set_xlabel("best value (A$ thousand a week)")
    ax2.set_title("the end of each run", loc="left", fontsize=13)
    source(fig, "The value is the shared problem's (economy.hpp) on the field's passers; heuresis's GA, CMA-ES (separable), SA, UMDA and random search as study (1) runs them, and the study's heat bath.")
    save(fig, "search-86-searches-at-equal-budgets")


# ---- 87: the plan found, and the relations it uses -------------------------------------------------------------------------


def draw_plan(ax, plan, title_text):
    base_map(ax, frontages=False, area=True)
    network(ax, None, colour="#D9CCD3", widths=1.2)
    kinds = plan_kinds(plan["plan"])
    ax.add_collection(PolyCollection(FOOT, facecolor=[KIND_COLOUR[k] for k in kinds], edgecolor=PAPER, linewidth=0.8, zorder=3))
    for j, made in enumerate(plan["plan"]["links"]):
        if made:
            a, b, kind = LINKS[j]
            ax.plot([a[0], b[0]], [a[1], b[1]], color="#C2477F", lw=2.6, ls=(0, (3, 1.5)), zorder=8, solid_capstyle="round")
    ax.set_title(title_text, loc="left", fontsize=12.5)


def fig87():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The plan the heat bath found, and the relations it uses", "The 61 frontages' kinds as built and in the best plan; the catalysis and competition between them in that plan; the links made, dashed")
    v0 = AS_BUILT["outcome"]["value_aud_week"]
    vb = BEST_HB["outcome"]["value_aud_week"]
    ax1 = fig.add_axes([0.01, 0.12, 0.32, 0.74])
    draw_plan(ax1, AS_BUILT, f"as built: A$ {v0:,.0f} a week")
    ax2 = fig.add_axes([0.34, 0.12, 0.32, 0.74])
    draw_plan(ax2, BEST_HB, f"the heat bath's best: A$ {vb:,.0f} a week")
    scale_bar(ax2, 50, area=True)
    ax3 = fig.add_axes([0.67, 0.12, 0.32, 0.74])
    base_map(ax3, frontages=False, area=True)
    network(ax3, None, colour="#E2D6D0", widths=1.0)
    kinds = plan_kinds(BEST_HB["plan"])
    ax3.add_collection(PolyCollection(FOOT, facecolor=[KIND_COLOUR[k] for k in kinds], edgecolor=PAPER, linewidth=0.6, zorder=3, alpha=0.45))
    A = np.array(PROBLEM["affinity"]["a"])
    cen = np.array([np.mean(f, axis=0) for f in FOOT])
    xi = PROBLEM["parameters"]["catalysis"]["coupling_length_m"]
    for f, g, d in PROBLEM["relations"]["reach_within_3_coupling_lengths"]:
        f, g = int(f), int(g)
        kf, kg = kinds[f], kinds[g]
        if kf in (0, 4) or kg in (0, 4):
            continue
        w = math.exp(-d / xi)
        a = 0.5 * (A[kg][kf] + A[kf][kg])
        if abs(a * w) < 0.03:
            continue
        ax3.plot([cen[f, 0], cen[g, 0]], [cen[f, 1], cen[g, 1]], color="#5E8C57" if a > 0 else "#C2477F", lw=5 * abs(a) * w, alpha=0.75, zorder=4, solid_capstyle="round")
    for s_f in FRONT:
        f = s_f["index"]
        for g in s_f["same_segment"]:
            g = int(g)
            if g > f and kinds[f] == kinds[g] and kinds[f] in (1, 2, 3, 5):
                ax3.plot([cen[f, 0], cen[g, 0]], [cen[f, 1], cen[g, 1]], color="#D49A3A", lw=1.6, ls=(0, (1, 1.2)), zorder=5)
    for a, b, area in PROBLEM["relations"]["party_walls"]:
        a, b = int(a), int(b)
        ax3.plot([cen[a, 0], cen[b, 0]], [cen[a, 1], cen[b, 1]], color=INK, lw=0.6, alpha=0.55, zorder=4.5)
    ax3.set_title("its relations", loc="left", fontsize=12.5)
    fig.legend(handles=[Patch(facecolor=KIND_COLOUR[k], label=KINDS[k]) for k in range(6)] + [Line2D([], [], color="#C2477F", lw=2.2, ls="--", label="a link made"),
               Line2D([], [], color="#5E8C57", lw=3, label="catalysis +, by strength"), Line2D([], [], color="#C2477F", lw=3, label="catalysis −"),
               Line2D([], [], color="#D49A3A", lw=1.6, ls=":", label="competitors on one segment"), Line2D([], [], color=INK, lw=0.8, label="a party wall")],
               loc="lower center", ncol=6, fontsize=10.5, bbox_to_anchor=(0.5, 0.025))
    source(fig, CREDIT)
    save(fig, "search-87-the-plan-and-its-relations")

# ---- 88: the week of the plan ---------------------------------------------------------------------------------------------
def fig88():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The week of the plan: time and frequency", "Spending and frontages open by hour, the passers the plan draws, and the week's lines of its spending")
    plans = [(AS_BUILT, "as built", INK), (BEST_HB, "the heat bath's best", "#8E3B78"), (BEST_BB, BEST_BB["name"].replace("best black box", "the best black box"), "#4E8F8B")]
    ax = fig.add_axes([0.06, 0.58, 0.62, 0.28])
    for p, lab, col in plans:
        ax.plot(np.arange(WEEK) + 0.5, np.array(p["outcome"]["spending_by_hour"]) / 1000, color=col, lw=1.7, label=f"{lab}: A$ {p['outcome']['spending_aud_week']/1000:,.0f}k spent a week")
    hours_axis(ax)
    ax.set_ylabel("spending (A$ thousand an hour)")
    ax.set_ylim(0, 66)
    ax.legend(loc="upper left", fontsize=10, ncol=3, frameon=False)
    ax.set_title("spending by hour", loc="left", fontsize=13)
    ax2 = fig.add_axes([0.06, 0.32, 0.62, 0.18])
    for p, lab, col in plans:
        cnt = np.zeros(WEEK)
        for f in range(NF):
            if p["plan"]["kind"][f] in (1, 2, 3, 5):
                cnt += open_mask(p["plan"], f)
        ax2.step(np.arange(WEEK + 1), np.concatenate([cnt, [cnt[-1]]]), where="post", color=col, lw=1.4)
    hours_axis(ax2)
    ax2.set_ylabel("trading frontages open")
    ax3 = fig.add_axes([0.06, 0.07, 0.62, 0.18])
    P0 = np.array(AS_BUILT["passers"]).sum(0)
    P1 = np.array(BEST_HB["passers"]).sum(0)
    ax3.fill_between(np.arange(WEEK) + 0.5, 0, P1 - P0, color="#8E3B78", alpha=0.6, step="mid")
    hours_axis(ax3)
    ax3.axhline(0, color=MUTED, lw=0.6)
    ax3.set_ylabel("passers drawn")
    ax3.set_title(f"the plan's attraction: passers at the 61 frontages, the best plan less as built ({(P1.sum() - P0.sum()) / P0.sum() * 100:+.2f}% of the week)", loc="left", fontsize=11.5)
    ax4 = fig.add_axes([0.75, 0.32, 0.22, 0.54])
    for p, lab, col in plans:
        pw = np.array(p["outcome"]["power"])
        ax4.semilogy(np.arange(len(pw)), np.maximum(pw / pw.sum(), 1e-7), color=col, lw=1.2, marker="o", ms=2.5)
    for j in range(1, 7):
        ax4.axvline(7 * j, color=HAIRLINE, lw=0.6)
    ax4.set_xlim(-1, 50)
    ax4.set_ylim(1e-5, 1)
    ax4.set_xlabel("line of the week, m")
    ax4.set_ylabel("share of the spending's power")
    ax4.set_title("the spending's lines (Parseval)", loc="left", fontsize=13)
    ax5 = fig.add_axes([0.75, 0.07, 0.22, 0.17])
    ax5.axis("off")
    for i, (p, lab, col) in enumerate(plans):
        o = p["outcome"]
        ax5.text(0, 0.95 - 0.33 * i, f"{lab}\nspread {o['spread']:.3f} ({o['effective_hours']:.0f} h), evening {o['evening_share']:.2f}, weekend {o['weekend_share']:.2f}, mix {o['mix_nats']:.2f}, access {o['access_m']:.0f} m",
                 fontsize=9, color=col, va="top", transform=ax5.transAxes, wrap=True)
    source(fig, CREDIT)
    save(fig, "search-88-the-week-of-the-plan")


# ---- 89: links by linear response, plans by free energy -----------------------------------------------------------------
def fig89():
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Links by linear response, plans by free energy", "Each link's effect on the passers from the action's propagator alone, against the lattice with the link and the domain's rebuilt complex; and the plans ranked by F_T = E V − Var V / 2T")
    L = SRCH["links"]
    ax = fig.add_axes([0.20, 0.12, 0.34, 0.72])
    y = np.arange(len(L))
    ex = np.array([l["weekly_passers_change"]["exact"] for l in L])
    li = np.array([l["weekly_passers_change"]["linear_response"] for l in L])
    rb = np.array([l["weekly_passers_change"]["domain_rebuilt"] for l in L])
    ax.barh(y + 0.27, rb / 1000, height=0.25, color=ALGORITHM["random"], label="the domain's complex rebuilt, subdivided alike")
    ax.barh(y, ex / 1000, height=0.25, color="#8E3B78", label="the lattice with the link (Sherman–Morrison, exact)")
    ax.barh(y - 0.27, li / 1000, height=0.25, color="#CA8DBE", label="first order: the connected correlation (FDT)")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{l['joins'][0]} – {l['joins'][1]}" for l in L], fontsize=9.5)
    ax.axvline(0, color=MUTED, lw=0.6)
    ax.set_xlabel("change of the 61 frontages' weekly passers (thousand)")
    ax.legend(loc="lower left", fontsize=9, bbox_to_anchor=(0.0, -0.005))
    cors = [l["frontage_correlation"]["linear_response"] for l in L]
    corr_rb = [l["frontage_correlation"]["domain_rebuilt"] for l in L]
    ax.set_title(f"each link alone, from the plan as built (median r with the exact:\nfirst order {np.median(cors):.2f}, rebuilt {np.median(corr_rb):.2f})", loc="left", fontsize=11)
    ax2 = fig.add_axes([0.62, 0.12, 0.25, 0.72])
    Ts = np.logspace(3.5, 7, 60)
    shown = [p for p in SRCH["plans"] if p["name"] in ("as built", "best heat bath")] + [BEST_BB]
    for p in SRCH["plans"]:
        if p["name"].startswith("best ") and p["name"] not in ("best heat bath",) and "black box" not in p["name"] and all(p["plan"] != q["plan"] for q in shown):
            shown.append(p)
    for p in shown:
        col = INK if p["name"] == "as built" else ("#8E3B78" if p["name"] == "best heat bath" else ("#4E8F8B" if p is BEST_BB else "#B3A9AF"))
        ev = p["outcome"]["value_aud_week"]
        var = p["week_sd"] ** 2 + p["posterior_sd"] ** 2
        fe = ev - var / (2 * Ts)
        ax2.plot(Ts, fe / 1000, color=col, lw=2.0 if col != "#B3A9AF" else 1.2)
        pb = p["difference_from_as_built"]["probability_better"]
        lab = p["name"].replace("best black box", "black box") + ("" if pb is None else f"\nP(better than as built) {pb:.3f}")
        ax2.text(1.15e7, fe[-1] / 1000, lab, fontsize=9, color=col, va="center")
    ax2.set_xscale("log")
    ax2.set_xlim(Ts[0], Ts[-1])
    ax2.set_ylim(-200, 760)
    ax2.set_xlabel("temperature T (A$ a week): risk tolerance")
    ax2.set_ylabel("free energy F_T (A$ thousand a week)")
    ax2.set_title("free energies, under the week's fluctuations and the posterior", loc="left", fontsize=11)
    source(fig, CREDIT)
    save(fig, "search-89-links-and-free-energies")

if __name__ == "__main__":
    which = sys.argv[1:] or ["80", "81", "82", "83", "84", "85", "86", "87", "88", "89"]
    for w in which:
        globals()[f"fig{w}"]()
        print("drew", w)
