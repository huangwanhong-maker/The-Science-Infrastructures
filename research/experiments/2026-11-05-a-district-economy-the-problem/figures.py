"""The figures of the shared problem "A district's economy", in the house style of the search case studies
(research/experiments/style/search_style.py), the place in the Urban Science Workbench's sea blue and sand.

    .venv/Scripts/python.exe research/experiments/2026-11-05-a-district-economy-the-problem/figures.py

Draws from results/problem.json, results/sensors.json and the uncommitted results/district-geometry.json (made by
calibrate from the cache), and writes assets/gallery/search-60-... and search-61-... (1600 x 1000).
"""

from __future__ import annotations

import json
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F401,F403

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"

# The kinds' colours, shared by the three studies: the place's sea blue and sand among the house's family.
KIND_COLOUR = {"vacant": "#DCD3D7", "retail": "#5B8FB9", "food and drink": "#8E3B78", "services": "#7FA37A", "office lobby": "#C98B6B", "culture": "#D49A3A"}
KINDS = list(KIND_COLOUR)
SEA, SAND = "#5B8FB9", "#C98B6B"
CONTEXT = "#ECE6DF"
STREET = "#CFC4CB"
PROVENANCE = "Data: OpenStreetMap contributors (ODbL 1.0); City of Melbourne, CLUE and Pedestrian Counting System (CC BY 4.0)"


def load(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


class Frame:
    """The map turned so that the area's long axis (Flinders Lane to Bourke Street, along Swanston Street) runs left to
    right: a rotation of the local frame about the area's centre."""

    def __init__(self, quad):
        q = np.array(quad, float)
        self.c = q.mean(axis=0)
        s0, s1 = (q[0] + q[1]) / 2, (q[2] + q[3]) / 2
        self.theta = math.atan2(s1[1] - s0[1], s1[0] - s0[0])
        ct, st = math.cos(-self.theta), math.sin(-self.theta)
        self.R = np.array([[ct, -st], [st, ct]])

    def __call__(self, pts):
        a = np.asarray(pts, float)
        return (a - self.c) @ self.R.T

    def north(self):
        """The direction of north on the turned map."""
        return self.R @ np.array([0.0, 1.0])


def district_axes(ax, T, G, pad=(70, 40)):
    q = T(G["quad"])
    ax.set_xlim(q[:, 0].min() - pad[0], q[:, 0].max() + pad[0])
    ax.set_ylim(q[:, 1].min() - pad[1], q[:, 1].max() + pad[1])
    ax.set_aspect("equal")
    ax.axis("off")


def north_arrow(ax, T, at):
    n = T.north() * 22
    ax.annotate("", xy=(at[0] + n[0], at[1] + n[1]), xytext=at, arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.0))
    ax.text(at[0] + 1.35 * n[0], at[1] + 1.35 * n[1], "N", color=MUTED, fontsize=10, ha="center", va="center")


def draw_base(ax, T, G, frontage_colour=None, alpha=1.0, edge=PAPER):
    ctx = [T(c["footprint"]) for c in G["context"]]
    ax.add_collection(PolyCollection(ctx, facecolor=CONTEXT, edgecolor="#E2D9D3", linewidth=0.4, zorder=1))
    for sp in G["open_spaces"]:
        ax.add_patch(plt.Polygon(T(sp["ring"]), facecolor="#DCE8D3", edgecolor="none", zorder=1))
    segs = [T(sg["line"]) for sg in G["segments"]]
    widths = [1.8 if sg["type"] == "street" else (1.5 if sg["type"] == "mall" else 0.8) for sg in G["segments"]]
    ax.add_collection(LineCollection(segs, colors=STREET, linewidths=widths, zorder=2))
    if frontage_colour is not None:
        polys = [T(f["footprint"]) for f in G["frontages"]]
        ax.add_collection(PolyCollection(polys, facecolor=frontage_colour, edgecolor=edge, linewidth=0.6, alpha=alpha, zorder=3))


STREET_LABELS = {"Swanston Street": (-205, 0, 90), "Russell Street": (205, 0, 90), "Flinders Lane": (0, -130, 0), "Bourke Street": (0, 128, 0)}


def street_labels(ax, T, G):
    """The four streets that bound the area, beside its edges (in the turned frame: Flinders Lane at the left)."""
    q = T(G["quad"])
    cx, cy = q[:, 0].mean(), q[:, 1].mean()
    w, h = (q[:, 0].max() - q[:, 0].min()) / 2, (q[:, 1].max() - q[:, 1].min()) / 2
    ax.text(cx - w - 22, cy, "Flinders Lane", rotation=90, ha="center", va="center", fontsize=10, color=MUTED, style="italic", zorder=9)
    ax.text(cx + w + 22, cy, "Bourke Street", rotation=90, ha="center", va="center", fontsize=10, color=MUTED, style="italic", zorder=9)
    ax.text(cx, cy + h + 16, "Swanston Street", ha="center", va="center", fontsize=10, color=MUTED, style="italic", zorder=9)
    ax.text(cx, cy - h - 16, "Russell Street", ha="center", va="center", fontsize=10, color=MUTED, style="italic", zorder=9)


def figure_60():
    P, G, S = load("problem.json"), load("district-geometry.json"), load("sensors.json")
    fr = P["frontages"]
    T = Frame(G["quad"])
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "A district's economy: the problem the three studies share",
          "61 ground-floor frontages between Swanston and Russell Streets, Flinders Lane and Bourke Street: their kinds as built, how they touch, what they can be")
    ax = fig.add_axes([0.02, 0.30, 0.58, 0.58])
    colours = [KIND_COLOUR[f["kind_as_built"]] for f in fr]
    draw_base(ax, T, G, colours)
    anchors = T([g["anchor"] for g in G["frontages"]])
    centres = T([g["centre"] for g in G["frontages"]])
    # Reach: frontages within 60 m of each other along the network (faint); party walls (solid, centre to centre).
    reach = [(anchors[a], anchors[b]) for a, b, d in P["relations"]["reach_within_3_coupling_lengths"] if d <= 60]
    ax.add_collection(LineCollection(reach, colors=MUTED, linewidths=0.4, alpha=0.5, zorder=4))
    party = [(centres[a], centres[b]) for a, b, _ in P["relations"]["party_walls"]]
    ax.add_collection(LineCollection(party, colors=INK, linewidths=1.0, alpha=0.8, zorder=5))
    ax.scatter(anchors[:, 0], anchors[:, 1], s=[10 + 0.9 * f["length_m"] for f in fr], c=colours, edgecolors=INK, linewidths=0.6, zorder=6)
    for l in G["links"]:
        ls = (0, (3, 1.6)) if l["kind"] == "lane carried through" else (0, (1, 1.2))
        ab = T([l["a"], l["b"]])
        ax.plot(ab[:, 0], ab[:, 1], color=SAND, lw=2.6, ls=ls, zorder=7)
    district_axes(ax, T, G)
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    for sn in S["sensors"]:
        q = T([[sn["x"], sn["y"]]])[0]
        if xlim[0] <= q[0] <= xlim[1] and ylim[0] <= q[1] <= ylim[1]:
            ax.scatter([q[0]], [q[1]], marker="s", s=44, color=SEA, edgecolors=PAPER, linewidths=0.9, zorder=8)
    quad = T(G["quad"] + [G["quad"][0]])
    ax.plot(quad[:, 0], quad[:, 1], color=INK, lw=0.7, ls=(0, (6, 3)), alpha=0.45, zorder=2)
    street_labels(ax, T, G)
    north_arrow(ax, T, (xlim[1] - 30, ylim[0] + 30))
    handles = [Patch(facecolor=KIND_COLOUR[k], edgecolor=PAPER, label=f"{k} ({sum(1 for f in fr if f['kind_as_built'] == k)})") for k in KINDS if k != "culture"]
    lh = [Line2D([], [], color=INK, lw=1.0, label="a party wall (54 pairs)"), Line2D([], [], color=MUTED, lw=0.6, alpha=0.7, label="within 60 m along the network"),
          Line2D([], [], color=SAND, lw=2.6, ls=(0, (3, 1.6)), label="a lane carried through (8)"), Line2D([], [], color=SAND, lw=2.6, ls=(0, (1, 1.2)), label="a passage (3)"),
          Line2D([], [], marker="s", color=SEA, lw=0, markersize=7, label="a pedestrian counter"),
          Line2D([], [], marker="o", color=PAPER, markeredgecolor=INK, lw=0, markersize=7, label="a frontage (size: its length)")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.05, 0.27), fontsize=11, title="ground floors as built (CLUE 2024)", title_fontsize=11.5, alignment="left")
    fig.legend(handles=lh, loc="upper left", bbox_to_anchor=(0.27, 0.27), fontsize=11, title="relations and interventions", title_fontsize=11.5, alignment="left")

    # The affinities.
    a = np.array(P["affinity"]["clq"], float)
    ax2 = fig.add_axes([0.66, 0.54, 0.28, 0.31])
    la = np.log(a)
    im = ax2.imshow(la, cmap=DUSK.reversed(), vmin=-1, vmax=1)
    short = ["vacant", "retail", "food &\ndrink", "services", "office\nlobby", "culture"]
    ax2.set_xticks(range(6), short, fontsize=10)
    ax2.set_yticks(range(6), short, fontsize=10)
    ax2.grid(False)
    for i in range(6):
        for j in range(6):
            ax2.text(j, i, f"{a[i, j]:.2f}", ha="center", va="center", fontsize=9.5, color=PAPER if abs(la[i, j]) > 0.6 else INK)
    ax2.set_title("co-location quotient, CLUE 2024, within 50 m", fontsize=12.5, loc="left", pad=8)
    ax2.set_xlabel("kind j, found among the neighbours of kind i (rows)", fontsize=11)
    cb = fig.colorbar(im, ax=ax2, fraction=0.04, pad=0.02)
    cb.set_label("ln CLQ", fontsize=10)
    cb.outline.set_visible(False)

    # The hours as built and the value's ingredients.
    ax3 = fig.add_axes([0.71, 0.09, 0.25, 0.22])
    kinds = ["retail", "food and drink", "services", "office lobby", "culture"]
    for i, k in enumerate(kinds):
        h = P["kinds"][KINDS.index(k)]["hours_as_built"]
        for j, (key, off, alpha) in enumerate((("weekday", -0.17, 1.0), ("weekend", 0.17, 0.55))):
            w = h[key]
            if w and w[1] > w[0]:
                ax3.barh(i + off, w[1] - w[0], left=w[0], height=0.3, color=KIND_COLOUR[k], alpha=alpha, edgecolor="none")
    ax3.set_yticks(range(len(kinds)), kinds, fontsize=10.5)
    ax3.invert_yaxis()
    ax3.set_xlim(0, 24)
    ax3.set_xticks([0, 6, 9, 12, 15, 18, 21, 24])
    ax3.set_xlabel("civil hour: weekday (full) and weekend (pale)", fontsize=11)
    ax3.set_title("hours as built: OpenStreetMap's medians", fontsize=12.5, loc="left")
    ax3.grid(axis="y", visible=False)
    fig.text(0.635, 0.425, "the value, A$ a week, summed over the frontages' open hours:", fontsize=11.5, color=MUTED)
    fig.text(0.635, 0.395, "V = Σ passers × chance to stop × competition × spend × margin", fontsize=12, color=INK)
    fig.text(0.635, 0.367, "      − staffing − outgoings − refits − links", fontsize=12, color=INK)
    source(fig, PROVENANCE + ". A study of methods on a real place, not a planning proposal.")
    save(fig, "search-60-district-economy-the-problem")


def figure_61():
    P, S = load("problem.json"), load("sensors.json")
    G = load("district-geometry.json")
    xs = [p[0] for p in G["quad"]]
    ys = [p[1] for p in G["quad"]]
    box = P  # noqa: F841
    sensors = S["sensors"]
    near = [s for s in sensors if min(xs) - 120 <= s["x"] <= max(xs) + 120 and min(ys) - 120 <= s["y"] <= max(ys) + 120]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The people of the district: a tide with one crest a day",
          "The City's counters over the 52 weeks from 27 October 2025, read as tides by the urban-science domain: the mean week, its lines, its spectrum, its response")
    hours = np.arange(168) + 0.5
    # The mean weeks.
    ax = fig.add_axes([0.06, 0.53, 0.52, 0.33])
    for s in near:
        w = np.array(s["mean_week"])
        ax.plot(hours, w, color=SEA, lw=0.8, alpha=0.45)
    hi = next(s for s in sensors if s["name"] == "Town Hall (West)")
    ax.plot(hours, hi["mean_week"], color="#8E3B78", lw=2.0, label="Town Hall (West), Swanston Street")
    ru = next(s for s in sensors if s["name"].startswith("155-161 Russell"))
    ax.plot(hours, ru["mean_week"], color=SAND, lw=2.0, label="155–161 Russell Street")
    ax.plot([], [], color=SEA, lw=0.8, alpha=0.6, label=f"the {len(near)} counters about the area")
    ax.set_xlim(0, 168)
    ax.set_xticks(np.arange(0, 169, 24), [""] * 8)
    for d, nm in enumerate(["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]):
        ax.text(24 * d + 12, -0.06, nm, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=11, color=MUTED)
    ax.set_ylabel("people an hour (mean week)")
    ax.set_ylim(0, 4900)
    ax.legend(loc="upper left", fontsize=10.5, ncol=3, bbox_to_anchor=(0.0, 1.0))
    ax.set_title("the mean week: the week's constituents alone, public holidays left out", fontsize=12.5, loc="left")

    # The week's lines of the Town Hall's mean week.
    ax2 = fig.add_axes([0.65, 0.53, 0.31, 0.33])
    w = np.array(hi["mean_week"])
    X = np.fft.rfft(w) / 168
    pw = np.abs(X) ** 2
    pw[1:-1] *= 2
    m = np.arange(len(pw))
    share = pw / pw.sum()
    daily = (m % 7 == 0) & (m > 0)
    ax2.vlines(m[~daily & (m > 0)], 0, share[~daily & (m > 0)], color=SEA, lw=2)
    ax2.vlines(m[daily], 0, share[daily], color="#8E3B78", lw=3)
    ax2.set_yscale("log")
    ax2.set_ylim(1e-6, 1)
    ax2.set_xlim(-1, 50)
    for j in range(1, 7):
        ax2.text(7 * j, share[7 * j] * 1.8, f"D{j}", ha="center", fontsize=10, color="#8E3B78")
    ax2.text(0.8, share[0] * 0.55, f"the mean, {share[0]:.0%}", fontsize=10, color=INK)
    ax2.vlines([0], 0, share[0], color=INK, lw=3)
    ax2.set_xlabel("line of the week, m (cycles a week)")
    ax2.set_ylabel("share of the power")
    ax2.set_title("Town Hall (West): the week's lines (Parseval)", fontsize=12.5, loc="left")

    # Spectra.
    ax3 = fig.add_axes([0.06, 0.12, 0.52, 0.30])
    spec = {s["name"]: s for s in S["spectra"]}
    for name, colour in (("Town Hall (West)", "#8E3B78"), ("155-161 Russell Street", SAND), ("Flinders Ln -Degraves St (South)", "#4E8F8B")):
        if name not in spec:
            continue
        s = spec[name]
        f = np.array(s["cycles_per_day"])
        ax3.semilogy(f, s["density"], color=colour, lw=1.0, label=name)
        ax3.semilogy(f, s["level"], color=colour, lw=0.8, ls=(0, (3, 2)), alpha=0.8)
    for j in range(1, 5):
        ax3.axvline(j, color=HAIRLINE, lw=0.8, zorder=0)
        ax3.text(j, 1.02, f"D{j}", transform=ax3.get_xaxis_transform(), ha="center", fontsize=10, color=MUTED)
    ax3.set_xlim(0, 4.2)
    ax3.set_xlabel("cycles a day")
    ax3.set_ylabel("density of ln(1 + count)")
    ax3.legend(loc="upper right", fontsize=10)
    ax3.set_title("multitaper spectra (NW = 4, 7 tapers) and the red-noise level (dashed)", fontsize=12.5, loc="left", pad=16)

    # The response.
    ax4 = fig.add_axes([0.65, 0.12, 0.31, 0.30])
    ins = {s["id"] for s in near}
    for s in sensors:
        r = s["response"]
        if "rain" not in r:
            continue
        c = "#8E3B78" if s["id"] in ins else MUTED
        ax4.scatter(r["rain"]["percent"], r["holiday"]["percent"], s=36 if s["id"] in ins else 20, color=c, alpha=0.9 if s["id"] in ins else 0.55, edgecolors="none")
    for s in sensors:
        label = {"Collins Place (North)": ("Collins Place (North)", (8, -2)), "Hamer Hall Entrance (COM Pole 1544)": ("Hamer Hall", (8, -4)),
                 "Town Hall (West)": ("Town Hall (West)", (-6, 8)), "Spencer St-Collins St (North)": ("Spencer St–Collins St", (-10, 8))}.get(s["name"])
        if s["name"].startswith("Hamer Hall"):
            label = ("Hamer Hall", (8, -4))
        if label:
            r = s["response"]
            ax4.annotate(label[0], (r["rain"]["percent"], r["holiday"]["percent"]), xytext=label[1], textcoords="offset points", fontsize=9.5, color=INK,
                         ha="right" if label[1][0] < 0 else "left")
    ax4.axhline(0, color=HAIRLINE, lw=1)
    ax4.set_xlabel("rain: % change in the counts a mm an hour")
    ax4.set_ylabel("a public holiday: % change")
    ax4.set_title("the response (errors clustered by day)", fontsize=12.5, loc="left")
    ax4.scatter([], [], color="#8E3B78", s=36, label="about the area")
    ax4.scatter([], [], color=MUTED, s=20, label="elsewhere in the Hoddle Grid")
    ax4.legend(loc="upper right", fontsize=10)
    source(fig, "Data: City of Melbourne, Pedestrian Counting System (the City's open data, CC BY 4.0); Open-Meteo from ERA5 (CC BY 4.0); Nager.Date. "
                "urban::rhythm (tide, spectrum) and urban::response on calx and tyche.")
    save(fig, "search-61-the-people-of-the-district")


if __name__ == "__main__":
    figure_60()
    figure_61()
    print("drawn")
