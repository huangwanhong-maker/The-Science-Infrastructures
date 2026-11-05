"""Figures of the case study "Owners and occupiers" (README.md), drawn from results/ in the house style of
research/experiments/style/search_style.py and saved as assets/gallery/agents-20-... to agents-29-... (1600 x 1000).

The maps draw the Hoddle Grid's buildings from the data library's occupancy product (their footprints from the city
model, OpenStreetMap's) and the city model's other buildings as context, turned so that the grid's streets run across
the page; north is marked. Places are named only as the data library's naming policy allows (streets by their names).

Colours. The place in the workbench's sea blue and sand; the uses in four hues validated for colour vision on the paper
(office #4A84BA, retail and hospitality #9A3F82, residential #D8952A, vacant #C2552F; all pairs pass the dataviz
validator), other occupied uses in graphite, works hatched. The models: observed in ink, the full model in sea blue, the
model without its neighbours in amber, the chain alone in graphite, dashed.

Usage (from the repository's root):
    python research/experiments/2026-11-05-owners-and-occupiers/figures.py [results directory or -] [figures, e.g. 20,24]
"""

from __future__ import annotations

import csv
import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "style"))
from search_style import HAIRLINE, INK, MUTED, PANEL, PAPER, figure, save, title  # noqa: E402
from search_style import source as _source  # noqa: E402
import textwrap  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = HERE.parents[2]
RESULTS = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != "-" else HERE / "results"
PRODUCT = ROOT / "data" / "cache" / "products" / "occupancy-hoddle-grid.json"
CITY = ROOT / "data" / "cache" / "products" / "city-model-melbourne.json"

SEA, SAND = "#5B8FB9", "#C98B6B"          # the Urban Science Workbench's place
OFFICE, RETAIL, RESID, VACANT = "#4A84BA", "#9A3F82", "#D8952A", "#C2552F"
OTHER, WORKS_FACE, CONTEXT = "#B7AEB3", "#F1D9CC", "#E7DFDA"
OBS, FULL, NONB, CHAIN = INK, "#4A84BA", "#D8952A", "#9A9097"
GROUPS = [("Office", OFFICE), ("Retail and hospitality", RETAIL), ("Residential", RESID), ("Other occupied", OTHER),
          ("Vacant: unused", VACANT), ("Vacant: under works", WORKS_FACE)]
RUST = LinearSegmentedColormap.from_list("rust", ["#FBF6F1", "#F1D2C2", "#DE9A7C", "#C2552F", "#7A2E18"])
SEAMAP = LinearSegmentedColormap.from_list("sea", ["#FBF6F1", "#D3E2EE", "#8DB3D3", "#4A84BA", "#24496B"])
SOURCE = "CLUE (City of Melbourne, CC BY 4.0), 2002-2024; buildings © OpenStreetMap contributors (ODbL); urban-science occupancy on nomos agents, tyche decision, physis nested chains."


def source(fig, text: str) -> None:
    """The provenance, wrapped to the page's width."""
    _source(fig, "\n".join(textwrap.wrap(text, 215)))


def group_of(use: str) -> int:
    if use == "Unoccupied - Unused":
        return 4
    if use.startswith("Unoccupied"):
        return 5
    if use == "Office":
        return 0
    if use.startswith("Retail") or use.startswith("Entertainment") or use in ("Commercial Accommodation", "Performances, Conferences, Ceremonies", "Public Display Area"):
        return 1
    if use in ("Residential Apartment", "House/Townhouse", "Student Accommodation", "Institutional Accommodation"):
        return 2
    return 3


def load(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def buildings():
    with (RESULTS / "buildings.csv").open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ---- the map ------------------------------------------------------------------------------------------------------------

class Grid:
    """The Hoddle Grid turned so that Collins Street runs across the page."""

    def __init__(self):
        P = json.loads(PRODUCT.read_text(encoding="utf-8"))
        self.P = P
        city = json.loads(CITY.read_text(encoding="utf-8"))
        best, angle = 0, 0.0
        for r in city["roads"]:
            if r.get("name") == "Collins Street" and len(r["line"]) > 1:
                a, b = np.array(r["line"][0]), np.array(r["line"][-1])
                d = float(np.hypot(*(b - a)))
                if d > best:
                    best, angle = d, math.atan2(b[1] - a[1], b[0] - a[0])
        if angle > math.pi / 2:
            angle -= math.pi
        if angle < -math.pi / 2:
            angle += math.pi
        self.angle = angle
        self.c, self.s = math.cos(-angle), math.sin(-angle)
        props = P["properties"]["rows"]
        pts = np.array([r[3] for r in props], float)
        self.origin = pts.mean(axis=0)
        self.rings = {b: np.array(r, float) for b, r in P["footprints"]["rows"]}
        self.building_of = {r[0]: r[4] for r in props}
        box = P["box"]["local"]
        self.context = []
        for b in city["buildings"]:
            f = b.get("footprint") or []
            if len(f) < 3 or b["id"] in self.rings:
                continue
            a = np.array(f, float)
            if a[:, 0].max() < box[0] or a[:, 0].min() > box[2] or a[:, 1].max() < box[1] or a[:, 1].min() > box[3]:
                continue
            self.context.append(self.turn(a))
        self.blocks = [(blk, self.turn(np.array(r, float))) for blk, _, r in P["blocks"]["rows"]]
        t = self.turn(pts)
        lo, hi = np.percentile(t, 0.5, axis=0), np.percentile(t, 99.5, axis=0)
        self.core = (np.percentile(t[:, 0], 3), np.percentile(t[:, 0], 97), np.percentile(t[:, 1], 6), np.percentile(t[:, 1], 94))
        pad = 40
        self.extent = (lo[0] - 150, hi[0] + pad, lo[1] - 70, hi[1] + pad)
        self.names = []
        for r in city["roads"]:
            if r.get("name") in ("Collins Street", "Bourke Street", "Flinders Street", "Lonsdale Street", "La Trobe Street", "Elizabeth Street",
                                 "Swanston Street", "Spencer Street", "Spring Street", "William Street", "Queen Street", "Russell Street", "King Street",
                                 "Exhibition Street") and len(r["line"]) > 1:
                self.names.append((r["name"], self.turn(np.array(r["line"], float))))

    def turn(self, a):
        a = np.asarray(a, float) - getattr(self, "origin", np.zeros(2))
        return np.column_stack([self.c * a[:, 0] - self.s * a[:, 1], self.s * a[:, 0] + self.c * a[:, 1]])

    def ring_of(self, prop):
        b = self.building_of.get(prop)
        return self.turn(self.rings[b]) if b in self.rings else None

    def base(self, ax, blocks=True):
        ax.set_facecolor(PAPER)
        ax.add_collection(PolyCollection(self.context, facecolors=CONTEXT, edgecolors="none", zorder=1))
        if blocks:
            ax.add_collection(PolyCollection([r for _, r in self.blocks], facecolors="none", edgecolors=HAIRLINE, linewidths=0.6, zorder=0.5))
        ax.set_xlim(self.extent[0], self.extent[1])
        ax.set_ylim(self.extent[2], self.extent[3])
        ax.set_aspect("equal")
        ax.set_xticks([]), ax.set_yticks([])
        ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)

    def north(self, ax, x=0.985, y=0.86):
        # North in the turned frame: the local +y axis turned by −angle.
        dx, dy = -self.s, self.c
        ax.annotate("", xy=(x + 0.025 * dx, y + 0.045 * dy), xytext=(x, y), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1))
        ax.text(x + 0.04 * dx, y + 0.07 * dy, "N", transform=ax.transAxes, color=MUTED, fontsize=10, ha="center", va="center")

    def street_labels(self, ax, across=("Flinders Street", "Collins Street", "Bourke Street", "La Trobe Street"),
                      along=("Spencer Street", "William Street", "Swanston Street", "Exhibition Street")):
        """Streets across the page named in the left margin, at the height where they meet the grid's west edge; streets
        up the page named beneath the grid, where they meet its south edge."""
        done = set()
        x0, x1, y0, y1 = self.core
        for name, line in self.names:
            if name in done:
                continue
            d = line[-1] - line[0]
            horizontal = abs(d[0]) > abs(d[1])
            if name in across and horizontal:
                order = np.argsort(line[:, 0])
                xs, ys = line[order, 0], line[order, 1]
                if xs[0] > x0 + 60 or xs[-1] < x0:
                    continue
                y = float(np.interp(x0, xs, ys))
                ax.text(self.extent[0] + 4, y, name, fontsize=9, color=MUTED, ha="left", va="center", zorder=5, style="italic",
                        bbox=dict(facecolor=PAPER, edgecolor="none", alpha=0.85, pad=1.2))
                done.add(name)
            elif name in along and not horizontal:
                order = np.argsort(line[:, 1])
                xs, ys = line[order, 0], line[order, 1]
                if ys[0] > y0 + 60 or ys[-1] < y0:
                    continue
                x = float(np.interp(y0, ys, xs))
                ax.text(x, self.extent[2] + 6, name, fontsize=9, color=MUTED, ha="center", va="bottom", zorder=5, style="italic",
                        bbox=dict(facecolor=PAPER, edgecolor="none", alpha=0.85, pad=1.2))
                done.add(name)

    def polygons(self, ax, props, faces, hatch=None, edge="none", lw=0.0, z=2, values=None):
        """Footprints filled; with values, drawn in their ascending order, so that where several properties share one
        building the greatest is seen."""
        polys, fc = [], []
        pairs = list(zip(props, faces))
        if values is not None:
            pairs = [pf for _, pf in sorted(zip(values, pairs), key=lambda t: t[0])]
        for p, f in pairs:
            r = self.ring_of(p)
            if r is not None:
                polys.append(r), fc.append(f)
        ax.add_collection(PolyCollection(polys, facecolors=fc, edgecolors=edge, linewidths=lw, zorder=z, hatch=hatch))


def legend_groups(fig, y=0.075):
    handles = [Patch(facecolor=c, edgecolor=(VACANT if i == 5 else "none"), hatch=("////" if i == 5 else None), label=n) for i, (n, c) in enumerate(GROUPS)]
    fig.legend(handles=handles, loc="lower center", ncol=6, bbox_to_anchor=(0.5, y), fontsize=11.5, handlelength=1.4, columnspacing=1.6)


def uses_at(P, year):
    out = {}
    for y, prop, use, *_ in P["censuses"]["rows"]:
        if y == year:
            out[prop] = use
    return out


# ---- figures ------------------------------------------------------------------------------------------------------------

def fig_uses(G, memory):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The uses of the Hoddle Grid's buildings, 2002 and 2024",
          f"Each building's predominant use at CLUE's census, by four groups of its families; {memory['hoddle_grid']['units']:.0f} buildings, "
          f"{memory['hoddle_grid']['transitions']:.0f} transitions between censuses")
    for k, year in enumerate((2002, 2024)):
        ax = fig.add_axes([0.03 + 0.485 * k, 0.17, 0.47, 0.68])
        G.base(ax)
        U = uses_at(G.P, year)
        props = list(U)
        gs = [group_of(U[p]) for p in props]
        G.polygons(ax, [p for p, g in zip(props, gs) if g != 5], [GROUPS[g][1] for g in gs if g != 5])
        G.polygons(ax, [p for p, g in zip(props, gs) if g == 5], [WORKS_FACE] * sum(g == 5 for g in gs), hatch="////", edge=VACANT, lw=0.3)
        counts = np.bincount(gs, minlength=6)
        ax.text(0.0, 1.02, f"{year}", transform=ax.transAxes, fontsize=18, color=INK, ha="left", va="bottom")
        ax.text(0.10, 1.025, f"{len(props)} buildings: {counts[0]} office, {counts[4]} unused, {counts[5]} under works", transform=ax.transAxes,
                fontsize=11, color=MUTED, ha="left", va="bottom")
        G.street_labels(ax)
        if k == 1:
            G.north(ax)
    legend_groups(fig)
    m = memory["hoddle_grid"]
    source(fig, f"The nested chain fitted to every transition: a building half keeps its family for {m['half_life'][0]:.1f} ± {m['half_life_error'][0]:.1f} years and its use for "
                f"{m['half_life'][1]:.1f} ± {m['half_life_error'][1]:.1f}. " + SOURCE)
    return save(fig, "agents-20-the-uses-of-the-grid")


def fig_record(G, memory):
    S = load("simulated-record-district.json")
    # The district's properties, ordered along the grid's long axis.
    ids = [r[0] for r in S["properties"]["rows"]]
    at = {r[0]: r[1] for r in S["properties"]["rows"]}
    obs = {}
    for y, prop, use, *_ in G.P["censuses"]["rows"]:
        if prop in at:
            obs.setdefault(prop, {})[y] = use
    sim = {}
    for y, prop, use in S["censuses"]["rows"]:
        sim.setdefault(prop, {})[y] = use
    order = sorted(ids, key=lambda p: (group_of(next(iter(obs.get(p, {2002: "Office"}).values()))), G.turn(np.array([at[p]]))[0, 0]))
    order = [p for p in order if p in obs]
    years = list(range(2002, 2025))
    palette = [OFFICE, RETAIL, RESID, OTHER, VACANT, "#E9B8A3"]

    def image(rec):
        img = np.ones((len(order), len(years), 3))
        for i, p in enumerate(order):
            for j, y in enumerate(years):
                u = rec.get(p, {}).get(y)
                if u is None:
                    img[i, j] = (0.984, 0.965, 0.945)
                else:
                    h = palette[group_of(u)].lstrip("#")
                    img[i, j] = tuple(int(h[k:k + 2], 16) / 255 for k in (0, 2, 4))
        return img

    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The record and a simulated history, building by building",
          f"Event Spacetime's district ({len(order)} buildings, rows, grouped by their first use): CLUE's censuses (left) and one history the calibrated model "
          "simulated from the same first censuses (right)")
    for k, (rec, label) in enumerate(((obs, "Observed: CLUE's censuses"), (sim, "Simulated: provenance \"simulated\", never observed"))):
        ax = fig.add_axes([0.06 + 0.47 * k, 0.14, 0.42, 0.66])
        ax.imshow(image(rec), aspect="auto", interpolation="nearest", extent=(2001.5, 2024.5, len(order), 0))
        changes = sum(1 for p in order for a, b in zip(years, years[1:]) if rec.get(p, {}).get(a) and rec.get(p, {}).get(b) and rec[p][a] != rec[p][b])
        still = sum(1 for p in order if len(set(rec.get(p, {}).values())) <= 1) / len(order)
        ax.text(0.0, 1.06, label, transform=ax.transAxes, fontsize=14, color=INK if k == 0 else FULL, ha="left", va="bottom")
        ax.text(0.0, 1.015, f"{changes} changes of use between consecutive censuses; {still * 100:.0f} % of the buildings never change", transform=ax.transAxes,
                fontsize=11, color=MUTED, ha="left", va="bottom")
        ax.set_xticks([2002, 2006, 2010, 2014, 2018, 2022, 2024])
        ax.set_yticks([])
        ax.grid(False)
        if k == 0:
            ax.set_ylabel("buildings")
        if k == 1:
            for s in ax.spines.values():
                s.set_visible(True), s.set_color(FULL), s.set_linewidth(1.6), s.set_linestyle((0, (4, 3)))
    handles = [Patch(facecolor=c, label=n) for n, c in zip([g[0] for g in GROUPS], palette)]
    fig.legend(handles=handles, loc="lower center", ncol=6, bbox_to_anchor=(0.5, 0.045), fontsize=11.5, handlelength=1.4)
    m = S["model"]
    source(fig, f"The simulated history: seed {m['seed']:.0f}, λ0 {m['lambda'][0]:.3f}, λ1 {m['lambda'][1]:.3f}, ν {m['vacate']:.3f}, η {m['let']:.3f}, ρ {m['returning']:.2f}, "
                f"βn {m['beta_neighbours']:.2f}, βd {m['beta_demand']:.2f} a year (results/simulated-record-district.json, written in the city-events record's tables). " + SOURCE)
    return save(fig, "agents-21-the-record-and-a-simulated-history")


SPECS = [("modulated", "the full model", FULL, "-"), ("without_neighbours", "without the neighbours", NONB, "-"), ("chain_alone", "the chain alone", CHAIN, "--")]


GLOSS = {
    "ln λ̂_0": ("the chain's λ0, fitted (log)", "the nested chain on every transition"),
    "ln λ̂_1": ("the chain's λ1, fitted (log)", "physis decision 0020"),
    "falling vacant": ("falling vacant", "occupied, then vacant, a year on"),
    "let": ("let", "vacant, then occupied, a year on"),
    "falling into works": ("falling into works", "occupied, then under works"),
    "returning to the family": ("returning to the family", "vacancies ending in the family before"),
    "neighbours' pull": ("the neighbours' pull", "a new use's affinity to its neighbours"),
    "demand: falling vacant": ("demand: falling vacant", "slope on the demand index"),
    "demand: let": ("demand: let", "slope on the demand index"),
    "moving between occupied families": ("moving between families", "occupied to another occupied family"),
}


def fig_calibration(cal):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Calibration by indirect inference", "The auxiliary model's ten statistics on CLUE's record against their means over panels simulated at each specification's estimate, "
                                                    "in standard deviations of a simulated panel")
    names = cal["modulated"]["psi_names"]
    ax = fig.add_axes([0.27, 0.11, 0.42, 0.70])
    yy = np.arange(len(names))[::-1]
    handles = []
    for k, (key, label, colour, ls) in enumerate(SPECS):
        c = cal[key]
        z = [(o - m) / s if s else 0 for o, m, s in zip(c["psi_observed"], c["psi_mean"], c["psi_sd"])]
        zc = np.clip(z, -12, 12)
        h = ax.scatter(zc, yy + (1 - k) * 0.22, s=70, color=colour, zorder=3, label=f"{label}: J {c['J']:.1f}", marker="o" if ls == "-" else "D",
                       edgecolors=PAPER, linewidths=1.2)
        handles.append(h)
        for zi, zr, y in zip(z, zc, yy + (1 - k) * 0.22):
            if abs(zi) > 12:
                ax.text(zr + (0.45 if zi > 0 else -0.45), y, f"{zi:.0f}", fontsize=9.5, color=MUTED, va="center", ha="left" if zi > 0 else "right")
    ax.axvspan(-2, 2, color=PANEL, zorder=0)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_yticks([])
    for y, n in zip(yy, names):
        main, gloss = GLOSS.get(n, (n, ""))
        ax.text(-0.02, y + 0.12, main, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=12.5, color=INK)
        ax.text(-0.02, y - 0.22, gloss, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=10, color=MUTED)
    ax.set_xlim(-14.5, 14.5)
    ax.set_ylim(-0.7, len(names) - 0.3)
    ax.set_xlabel("(observed − simulated) / simulated standard deviation; the band is ±2")
    ax.grid(axis="y", visible=False)
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.48, 0.875), ncol=3, fontsize=11.5)
    # The estimates.
    tx = fig.add_axes([0.72, 0.11, 0.26, 0.70])
    tx.axis("off")
    c = cal["modulated"]
    th, er = c["theta"], c["error"]
    rows = [("λ0, succeeded out of the family", math.exp(th[0]), math.exp(th[0]) * er[0], "a year"),
            ("λ1, succeeded within it", math.exp(th[1]), math.exp(th[1]) * er[1], "a year"),
            ("ν, left with no successor", math.exp(th[2]), math.exp(th[2]) * er[2], "a year"),
            ("η, a vacancy let", math.exp(th[3]), math.exp(th[3]) * er[3], "a year"),
            ("ρ, let back into its family", c["returning"], c["returning"] * (1 - c["returning"]) * er[4], ""),
            ("βn, the neighbours' pull", th[5], er[5], ""),
            ("βd, demand", th[6], er[6], "")]
    tx.text(0, 1.0, "The full model's estimates", fontsize=15, color=INK, va="top")
    for i, (n, v, e, u) in enumerate(rows):
        y = 0.92 - i * 0.085
        tx.text(0, y, n, fontsize=12, color=INK, va="top")
        tx.text(0, y - 0.036, f"{v:.3f} ± {e:.3f} {u}", fontsize=12.5, color=FULL, va="top")
    y = 0.92 - len(rows) * 0.085
    lines = [f"Hansen's statistic {c['hansen']:.1f} on {c['degrees']:.0f} degrees",
             f"of freedom (Σ from {c['replicates_sigma']:.0f} replicates): rejected,",
             "mainly by falling into works.",
             f"Without the neighbours {cal['without_neighbours']['hansen']:.0f} on {cal['without_neighbours']['degrees']:.0f};",
             f"the chain alone {cal['chain_alone']['hansen']:.0f} on {cal['chain_alone']['degrees']:.0f}. The neighbours'",
             f"distance metric {cal['distance_neighbours']:.1f} on 1 degree of freedom."]
    tx.text(0, y, chr(10).join(lines),
            fontsize=11, color=MUTED, va="top", linespacing=1.5)
    source(fig, f"Levenberg-Marquardt on the binding function, then Nelder-Mead; {cal['replicates']} replicates of common random numbers a value of J; "
                "errors by Gouriéroux, Monfort and Renault (1993) with a diagonal weight. " + SOURCE)
    return save(fig, "agents-22-calibration-by-indirect-inference")


def band(ax, x, b, colour, label, ls="-", alpha=0.18):
    ax.fill_between(x, b["q05"], b["q95"], color=colour, alpha=alpha, lw=0)
    ax.plot(x, b["q50"], color=colour, lw=2, ls=ls, label=label)


def fig_uses_in_time(ens, uses):
    years = np.arange(2002, 2025)
    names = [u["use"] for u in uses["uses"]]
    obs = np.array(ens["observed"]["counts"])  # years × uses
    fig, axes = figure(1, 2, gridspec_kw={"width_ratios": [1.25, 1]})
    fig.subplots_adjust(left=0.06, right=0.98, top=0.83, bottom=0.12, wspace=0.18)
    title(fig, "The uses in time: the record against the calibrated model",
          "Left: the share of each group of uses among the buildings standing, CLUE (points) and the full model (median, 5-95 % of 200 histories). "
          "Right: the vacant share, three models")
    ax = axes[0]
    tot = obs.sum(axis=1)
    sim = np.array(ens["modulated"]["use_counts"]).T  # years × uses (mean)
    for g, (n, c) in enumerate(GROUPS[:5]):
        idx = [k for k, u in enumerate(names) if group_of(u) == g]
        o = obs[:, idx].sum(axis=1) / tot
        s = sim[:, idx].sum(axis=1) / sim.sum(axis=1)
        ax.plot(years, s, color=c if g != 3 else "#8C8288", lw=2.2)
        ax.scatter(years, o, color=c if g != 3 else "#8C8288", s=26, zorder=3, edgecolors=PAPER, linewidths=0.8)
        ax.text(2024.4, s[-1], n, color=INK, fontsize=11, va="center")
    ax.set_xlim(2001.5, 2029)
    ax.set_ylabel("share of the buildings")
    ax.set_xticks([2002, 2006, 2010, 2014, 2018, 2022, 2024])
    ax.text(0.01, 0.98, "points: observed; lines: the full model's mean", transform=ax.transAxes, fontsize=10.5, color=MUTED, va="top")
    ax = axes[1]
    for key, label, colour, ls in SPECS:
        band(ax, years, ens[key]["vacant_share"], colour, label, ls)
    ax.plot(years, ens["observed"]["vacant_share"], color=OBS, lw=0, marker="o", ms=6, label="observed (CLUE)")
    ax.set_ylabel("share vacant (unused or under works)")
    ax.set_xticks([2002, 2006, 2010, 2014, 2018, 2022, 2024])
    ax.legend(loc="upper left", fontsize=11)
    v = ens["observed"]["vacant_share"]
    ax.annotate("from 2020 the record's vacancy" + chr(10) + "jumps; no model holds the jump", xy=(2020, v[18]), xytext=(2013.2, 0.19), fontsize=10.5, color=MUTED,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    source(fig, "The full model's histories start from each building's first census and are observed at its censuses, as the record is. " + SOURCE)
    return save(fig, "agents-23-the-uses-in-time")


def fig_where(G, rows):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Where vacancy persists", "Years each building stood vacant at CLUE's censuses 2002-2024: observed (left) and expected under the full model (right, the mean of "
                                         "200 histories from the same first censuses)")
    vmax = 12
    for k, key in enumerate(("observed_vacant", "expected_vacant_modulated")):
        ax = fig.add_axes([0.03 + 0.485 * k, 0.27, 0.47, 0.58])
        G.base(ax)
        props = [r["id"] for r in rows]
        vals = [min(float(r[key]), vmax) / vmax for r in rows]
        G.polygons(ax, props, [RUST(v) if v > 0 else "#EFE7E2" for v in vals], values=vals)
        ax.text(0.0, 1.02, "Observed" if k == 0 else "Expected, the full model", transform=ax.transAxes, fontsize=16, color=INK, ha="left", va="bottom")
        G.street_labels(ax)
        if k == 1:
            G.north(ax)
    cax = fig.add_axes([0.35, 0.215, 0.30, 0.018])
    sm = plt.cm.ScalarMappable(cmap=RUST, norm=plt.Normalize(0, vmax))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("censuses vacant, 2002-2024 (12 or more at the top)", color=INK)
    cb.outline.set_visible(False)
    # Below: observed against expected, by zone and in all.
    ax = fig.add_axes([0.04, 0.075, 0.92, 0.06])
    o = np.array([float(r["observed_vacant"]) for r in rows])
    e = np.array([float(r["expected_vacant_modulated"]) for r in rows])
    ax.axis("off")
    r = np.corrcoef(o, e)[0, 1]
    ever = (o > 0).mean()
    long = (o >= 5).sum()
    ax.text(0.5, 0.9, f"Observed and expected years vacant, building by building, correlate at {r:.2f}. {ever * 100:.0f} % of the buildings were vacant at some census; "
                      f"{long} for five censuses or more.", fontsize=12.5, color=INK, va="top", ha="center")
    source(fig, "The expected years carry the model's memory of each building's first use, its neighbours and demand; what they miss is persistence the model does not hold. " + SOURCE)
    return save(fig, "agents-24-where-vacancy-persists")


def fig_when(ens, uses, G):
    years = np.arange(2002, 2025)
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Vacancy in time, and how long it lasts", "Left: vacant floor space in the grid's blocks (CLUE's floor space by use, a lower bound) and vacant buildings, "
                                                          "with the demand index. Right: the chance a vacancy is let within a year, by its age")
    ax = fig.add_axes([0.06, 0.34, 0.40, 0.50])
    axd = fig.add_axes([0.06, 0.12, 0.40, 0.15], sharex=ax)
    axes = [ax, fig.add_axes([0.55, 0.12, 0.43, 0.72])]
    fs = G.P["floor_space"]["rows"]
    blocks_in = {r[2] for r in G.P["properties"]["rows"]}
    share = []
    for y in years:
        tot = sum(r[2] or 0 for r in fs if r[0] == y and r[1] in blocks_in)
        vac = sum(sum(v or 0 for v in r[3:8]) for r in fs if r[0] == y and r[1] in blocks_in)
        share.append(vac / tot if tot else np.nan)
    ax.plot(years, share, color=SAND, lw=2.4, marker="s", ms=5, label="vacant floor space (blocks)")
    ax.plot(years, ens["observed"]["vacant_share"], color=VACANT, lw=2.4, marker="o", ms=5, label="vacant buildings (predominant use)")
    ax.set_ylabel("share vacant")
    ax.legend(loc="upper left", fontsize=11)
    plt.setp(ax.get_xticklabels(), visible=False)
    d = uses["demand"]
    axd.bar(years, d, color=[SEA if v > 0 else "#B9CBDA" for v in d], width=0.7)
    axd.axhline(0, color=MUTED, lw=0.6)
    axd.set_ylabel("demand", color=SEA)
    axd.set_xticks([2002, 2006, 2010, 2014, 2018, 2022, 2024])
    axd.set_title("demand: the occupied establishments outside the box, ln, standardised", loc="left", fontsize=10.5, color=MUTED)
    ax = axes[1]
    obs = ens["observed"]
    k = np.arange(1, 12)
    risk, ex = np.array(obs["at_risk"])[1:12], np.array(obs["exits"])[1:12]
    p = np.where(risk > 0, ex / np.maximum(risk, 1), np.nan)
    se = np.sqrt(np.maximum(p * (1 - p), 1e-12) / np.maximum(risk, 1))
    for key, label, colour, ls in SPECS:
        b = ens[key]["exit_by_age"]
        bb = {q: np.array(b[q], float)[1:12] for q in ("q05", "q50", "q95")}
        band(ax, k, bb, colour, label, ls)
    ax.errorbar(k, p, yerr=1.96 * se, color=OBS, fmt="o", ms=6, capsize=0, lw=1.2, label="observed (95 %)")
    for x, n in zip(k, risk):
        ax.text(x, -0.03, f"{n:.0f}", fontsize=8.5, color=MUTED, ha="center", va="top")
    ax.set_xlabel("years vacant so far (spells whose start is seen; under each, how many)")
    ax.set_ylabel("chance of being let within the year")
    ax.set_ylim(-0.08, 0.6)
    ax.set_xticks(k)
    ax.legend(loc="upper right", fontsize=11)
    source(fig, "The record's chance of letting is low after one year and high after two: CLUE surveys a property at least every two years, carrying its use forward between. Vacant floor space: the five unoccupied columns of floor-space-by-use-by-block "
                "over the blocks of the box's buildings (small cells suppressed). " + SOURCE)
    return save(fig, "agents-25-vacancy-in-time-and-how-long")


def fig_why(ens):
    fig, axes = figure(1, 3)
    fig.subplots_adjust(left=0.06, right=0.98, top=0.82, bottom=0.14, wspace=0.28)
    title(fig, "Why vacancy persists: the neighbours", "Being let (left) and falling vacant (middle) within a year, by the share of the six nearest buildings vacant; "
                                                        "right, vacancy's clustering by year")
    labels = ["none", "up to 1/5", "1/5 to 2/5", "over 2/5"]
    obs = ens["observed"]
    for ax, (rk, ek, bk, ylab) in zip(axes[:2], (("nb_vac_risk", "nb_vac_exit", "exit_by_neighbours", "chance of being let"),
                                                 ("nb_occ_risk", "nb_occ_fall", "fall_by_neighbours", "chance of falling vacant"))):
        x = np.arange(4)
        for j, (key, label, colour, ls) in enumerate(SPECS):
            b = ens[key][bk]
            m, lo, hi = np.array(b["q50"], float), np.array(b["q05"], float), np.array(b["q95"], float)
            ax.errorbar(x + (j - 1) * 0.16 + 0.08, m, yerr=[m - lo, hi - m], fmt="D" if ls == "--" else "s", color=colour, ms=7, lw=1.4, capsize=0, label=label)
        r, e = np.array(obs[rk], float), np.array(obs[ek], float)
        p = e / np.maximum(r, 1)
        se = np.sqrt(np.maximum(p * (1 - p), 1e-12) / np.maximum(r, 1))
        ax.errorbar(x - 0.24, p, yerr=1.96 * se, fmt="o", color=OBS, ms=7, lw=1.4, capsize=0, label="observed (95 %)")
        for xi, n in zip(x, r):
            ax.text(xi, ax.get_ylim()[0], f"{n:.0f} pairs", fontsize=9, color=MUTED, ha="center", va="bottom")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_xlabel("neighbours vacant")
        ax.set_ylabel(ylab)
    axes[0].legend(loc="center left", bbox_to_anchor=(0.0, 0.40), fontsize=10.5)
    ax = axes[2]
    years = np.arange(2002, 2025)
    for key, label, colour, ls in SPECS:
        band(ax, years, ens[key]["clustering"], colour, label, ls)
    ax.plot(years, obs["clustering"], color=OBS, marker="o", lw=0, ms=6, label="observed")
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_ylabel("vacant share about a vacancy, less about an occupied building")
    ax.text(0.98, 0.97, "observed: points", transform=ax.transAxes, ha="right", va="top", fontsize=10.5, color=MUTED)
    ax.set_xticks([2002, 2008, 2014, 2020, 2024])
    source(fig, "Neighbours: each building's six nearest (assumption). The full model's owners choose a use by a logit with the neighbours' affinity (ln CLQ of the 2002 census, 60 m). " + SOURCE)
    return save(fig, "agents-26-why-vacancy-persists")


def fig_frequency(ens):
    fig, axes = figure(1, 2)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.76, bottom=0.12, wspace=0.22)
    title(fig, "Vacancy in frequency", "The pooled periodogram of each building's vacancy over the 23 censuses (buildings seen every year), against 200 simulated "
                                       "records of each model and the nested chain's own spectrum")
    n = 23
    f = np.arange(1, n // 2 + 1) / n
    for ax, (key, name, sp) in zip(axes, (("periodogram_vacant", "vacant (unused or under works)", "vacant"), ("periodogram_lettable", "vacant and unused", "lettable"))):
        for k, label, colour, ls in SPECS:
            band(ax, f, ens[k][key], colour, label, ls)
        ax.plot(f, ens["spectrum"][sp], color=MUTED, lw=1.4, ls=":", label="the chain's spectrum (closed form)")
        ax.plot(f, ens["observed"][key], color=OBS, marker="o", ms=6, lw=1.2, label="observed")
        ax.set_yscale("log")
        ax.set_xlabel("frequency (cycles a year)")
        ax.set_ylabel("power (variance a cycle a year)")
        ax.set_title(f"A building's indicator: {name}", loc="left", fontsize=14)
        sec = ax.secondary_xaxis("top", functions=(lambda x: 1 / np.maximum(x, 1e-3), lambda p: 1 / np.maximum(p, 1e-3)))
        sec.set_xticks([23, 10, 5, 3, 2])
        sec.set_xlabel("period (years)", color=MUTED)
    axes[0].legend(loc="upper right", fontsize=10.5)
    mu = ens["spectrum"]["mu"]
    source(fig, f"Closed form: Σ_m c_m (1 − r²)/(1 − 2r cos ω + r²), r = e^(−μ_m), with the fitted memories μ = {mu[1]:.3f}, {mu[2]:.3f} a year; the family's indicator lies in one mode, "
                f"the use's in two. {ens['observed']['complete']:.0f} buildings seen every year. " + SOURCE)
    return save(fig, "agents-27-vacancy-in-frequency")


RULES = [("the POMDP's ranked policy", FULL, "o"), ("every vacancy", RESID, "s"), ("the longest vacant", RETAIL, "^"), ("by neighbours", VACANT, "D")]


def fig_interventions(plan):
    fig, axes = figure(1, 2)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.82, bottom=0.12, wspace=0.22)
    title(fig, "Interventions on vacant space under a budget", "Left: discounted unused building-years averted in the planner's model, under the Lagrangian bound. "
                                                               "Right: building-years averted in the agent-based model, 2025-2034")
    P = plan["modulated"]["budgets"]
    B = np.array([b["budget"] for b in P])
    ax = axes[0]
    ceiling = np.array([b["bound"] - b["nothing"] for b in P])
    ax.fill_between(B, ceiling, ceiling * 0 + ceiling.max() * 1.25, color=PANEL, lw=0)
    ax.plot(B, ceiling, color=INK, lw=2, label="the Lagrangian bound: no plan averts more")
    for name, colour, marker in RULES:
        v, e = [], []
        for b in P:
            r = next(x for x in b["decoupled_averted"] if x["rule"] == name)
            v.append(r["years"]), e.append(r["error"])
        ax.errorbar(B, v, yerr=1.96 * np.array(e), color=colour, marker=marker, ms=8, lw=1.8, capsize=0, label=name)
        if name == RULES[0][0]:
            for x, y, c in zip(B, v, ceiling):
                ax.text(x * 0.96, y + 0.025 * ceiling.max(), f"{100 * y / c:.0f} %", fontsize=10.5, color=colour, va="bottom", ha="right")
    ax.set_xscale("log")
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xticks(B)
    ax.set_xticklabels([f"{b:.0f}" for b in B])
    ax.set_xlabel("budget (incentive-years a year)")
    ax.set_ylabel("discounted unused building-years averted")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left", fontsize=10.5)
    ax = axes[1]
    for name, colour, marker in RULES:
        v, e = [], []
        for b in P:
            r = next(x for x in b["abm_averted"] if x["rule"] == name)
            v.append(r["years"]), e.append(r["error"])
        ax.errorbar(B * (1 + 0.04 * RULES.index((name, colour, marker))), v, yerr=1.96 * np.array(e), color=colour, marker=marker, ms=8, lw=1.8, capsize=0, label=name)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xscale("log")
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xticks(B)
    ax.set_xticklabels([f"{b:.0f}" for b in B])
    ax.set_xlabel("budget (incentive-years a year)")
    ax.set_ylabel("unused building-years averted, 2025-2034 (paired, 95 %)")
    ax.legend(loc="upper left", fontsize=10.5)
    a = plan["assumptions"]
    h = P[0].get("horizon_decoupled", 0)
    source(fig, f"Left: the policy for {h:.0f} years then nothing, against doing nothing on the same draws ({P[0].get('episodes_decoupled', 0):.0f} episodes); the percentages are of the bound's ceiling. "
                f"Right: {P[0].get('episodes_abm', 0):.0f} paired episodes. Assumptions: a census surveys a building with chance {a['survey']:.1f} a year; an incentive costs 1, an inspection {a['cost'][1]:.1f}; it doubles a vacancy's letting "
                f"(κ = {a['incentive']:.0f}); discount {a['discount']:.2f} a year. The planner's model holds each building's neighbours at the 2024 census. " + SOURCE)
    return save(fig, "agents-28-interventions-under-a-budget")


def fig_planner(G, plan, rows, cal):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "What the neighbours add: where the planner acts, and how the models compare",
          "Left: the ten buildings each policy enrols in 2025 (budget ten), over the planner's belief that each is unused. Right: the models' fit and the policies' effect")
    ax = fig.add_axes([0.03, 0.19, 0.55, 0.66])
    G.base(ax)
    pm = plan["modulated"]
    bv = pm["belief_vacant"]
    props = [r["id"] for r in rows]
    # rows are in the order of buildings.csv, which omits buildings never seen; map by id through the product's order.
    order = [r[0] for r in G.P["properties"]["rows"]]
    idx = {p: i for i, p in enumerate(order)}
    vals = [bv[idx[p]] if bv[idx[p]] is not None else 0 for p in props]
    G.polygons(ax, props, [SEAMAP(min(1, v / 0.8)) if v > 0.02 else "#EFE7E2" for v in vals], values=vals)
    act = pm["first_actions"]["the POMDP's ranked policy"]
    ev = pm["first_actions"]["every vacancy"]
    pts_p, pts_e, pts_i = [], [], []
    point = {r[0]: r[3] for r in G.P["properties"]["rows"]}
    for p in props:
        a, b = act[idx[p]], ev[idx[p]]
        c = G.turn(np.array([point[p]], float))[0]  # the property's own point, inside its footprint
        if a == 2:
            pts_p.append(c)
        if a == 1:
            pts_i.append(c)
        if b == 2:
            pts_e.append(c)
    if pts_e:
        e = np.array(pts_e)
        ax.scatter(e[:, 0], e[:, 1], s=110, facecolors="none", edgecolors=RESID, linewidths=2, zorder=6, label=f"every vacancy enrols ({len(pts_e)})")
    if pts_p:
        q = np.array(pts_p)
        ax.scatter(q[:, 0], q[:, 1], s=30, color=VACANT, edgecolors=PAPER, linewidths=0.8, zorder=7, label=f"the POMDP enrols ({len(pts_p)})")
    if pts_i:
        q = np.array(pts_i)
        ax.scatter(q[:, 0], q[:, 1], s=40, marker="x", color=INK, linewidths=1.5, zorder=7, label=f"the POMDP inspects ({len(pts_i)})")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.01), ncol=3, fontsize=11)
    G.street_labels(ax)
    G.north(ax)
    cax = fig.add_axes([0.36, 0.115, 0.20, 0.014])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=SEAMAP, norm=plt.Normalize(0, 0.8)), cax=cax, orientation="horizontal")
    cb.set_label("belief in 2025 that the building is unused", color=INK, fontsize=11)
    cb.ax.xaxis.set_label_position("top")
    cb.outline.set_visible(False)
    # Right top: J by specification.
    ax = fig.add_axes([0.66, 0.58, 0.31, 0.24])
    names = ["the full model", "without the neighbours", "the chain alone"]
    Js = [cal[k]["J"] for k, *_ in SPECS]
    ax.barh(range(3)[::-1], Js, color=[FULL, NONB, CHAIN], height=0.55)
    ax.set_xscale("log")
    for i, (n, j) in enumerate(zip(names, Js)):
        ax.text(j * 1.1, 2 - i, f"{n}: J {j:.1f}", va="center", fontsize=11, color=INK)
    ax.set_yticks([])
    ax.set_xlim(1, max(Js) * 40)
    ax.set_xlabel("J, the distance from the record (log)")
    ax.grid(axis="y", visible=False)
    # Right bottom: averted years at budget 10 under the two models.
    ax = fig.add_axes([0.75, 0.15, 0.22, 0.32])
    for k, (key, label, colour) in enumerate((("modulated", "the full model", FULL), ("without_neighbours", "without the neighbours", NONB))):
        b = next(x for x in plan[key]["budgets"] if x["budget"] == 10)
        v = [next(x for x in b["abm_averted"] if x["rule"] == n)["years"] for n, *_ in RULES]
        e = [next(x for x in b["abm_averted"] if x["rule"] == n)["error"] for n, *_ in RULES]
        ax.barh(np.arange(4)[::-1] + (0.2 if k == 0 else -0.2), v, xerr=1.96 * np.array(e), height=0.36, color=colour, label=label, error_kw={"ecolor": MUTED, "lw": 1})
    ax.set_yticks(np.arange(4)[::-1])
    ax.set_yticklabels([n for n, *_ in RULES], fontsize=11)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xlabel("building-years averted, 2025-2034, budget 10")
    ax.legend(loc="upper center", bbox_to_anchor=(0.4, 1.16), ncol=2, fontsize=10.5)
    ax.grid(axis="y", visible=False)
    source(fig, "The planner's belief: its lumped POMDP (occupied, unused, under works) from the 2024 census, a year on; every candidate was seen unused in 2024, so the policy ranks them "
                "by their models, and its advantage grows as surveys skip buildings. " + SOURCE)
    return save(fig, "agents-29-what-the-neighbours-add")


def main():
    G = Grid()
    memory, cal, ens, uses, plan = load("memory.json"), load("calibration.json"), load("ensembles.json"), load("uses.json"), load("planning.json")
    rows = buildings()
    figs = {20: lambda: fig_uses(G, memory), 21: lambda: fig_record(G, memory), 22: lambda: fig_calibration(cal), 23: lambda: fig_uses_in_time(ens, uses),
            24: lambda: fig_where(G, rows), 25: lambda: fig_when(ens, uses, G), 26: lambda: fig_why(ens), 27: lambda: fig_frequency(ens),
            28: lambda: fig_interventions(plan), 29: lambda: fig_planner(G, plan, rows, cal)}
    only = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else sorted(figs)
    for k in only:
        print(figs[k]())


if __name__ == "__main__":
    main()
