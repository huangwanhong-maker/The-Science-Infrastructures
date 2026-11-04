"""The figures of the case study "Designing a district", drawn from results/ only, in the house style of the search case
studies (research/experiments/style/search_style.py), the place in the Urban Science Workbench's sea blue and sand.

    .venv/Scripts/python.exe research/experiments/2026-11-05-designing-a-district/figures.py [--only 50,51]

Writes assets/gallery/search-50-... to search-59-... (1600 x 1000).

The districts are drawn in an oblique parallel projection by the painter's algorithm: the ground and the shadows first,
then each building (its walls facing the viewer, then its roof), buildings from the farthest to the nearest. Walls and
roofs are lit by the Sun of the drawing's instant (a Lambertian term over an ambient one). Shadows are drawn on the
ground only: the prisms' shadows, swept from each footprint along the Sun's direction, as the study's sunlight model
casts them; shadows falling on other buildings are not drawn.
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
from matplotlib.collections import PolyCollection, LineCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, to_rgb  # noqa: E402
from matplotlib.patches import Polygon as MplPolygon, Rectangle  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"

# The place: the workbench's sea blue and sand, with the house's orchid and a pale sage for the public ground.
USE_COLOUR = {0: "#C98B6B", 1: "#5B8FB9", 2: "#8E3B78", 3: "#B3D19C"}   # dwellings, offices, shops, open space
USE_NAME = {0: "dwellings", 1: "offices and other workplaces", 2: "shops, food and entertainment", 3: "open space"}
CONTEXT = "#ECE6DF"
GROUND = PAPER
SHADOW = "#CBD0DA"
SHADOW_GREEN = "#91B28C"   # open ground in shadow
SEA = "#5B8FB9"
SAND = "#C98B6B"
LINK_STYLE = {"lane carried through": (0, (3, 1.6)), "passage": (0, (1, 1.2))}
SUN_MAP = LinearSegmentedColormap.from_list("winter_sun", ["#2F3E5C", "#5B8FB9", "#A9C3DA", "#F1E3CF", "#E9B97F", "#C98B6B"])
PROVENANCE = "Data: OpenStreetMap contributors (ODbL 1.0); City of Melbourne, CLUE (CC BY 4.0)"


def load(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def read_csv(name):
    with (RESULTS / name).open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ---- the oblique view -------------------------------------------------------------------------------------------

class View:
    """An orthographic view from azimuth `toward` (the direction the viewer looks, degrees from north through east)
    and elevation `elevation` above the horizon."""

    def __init__(self, toward_deg: float, elevation_deg: float):
        a, e = math.radians(toward_deg), math.radians(elevation_deg)
        self.f = np.array([math.sin(a), math.cos(a)])      # forward, horizontal
        self.r = np.array([math.cos(a), -math.sin(a)])     # right
        self.se, self.ce = math.sin(e), math.cos(e)

    def xy(self, p, z=0.0):
        p = np.asarray(p, float)
        u = p[..., 0] * self.r[0] + p[..., 1] * self.r[1]
        v = p[..., 0] * self.f[0] + p[..., 1] * self.f[1]
        return np.stack([u, v * self.se + np.asarray(z) * self.ce], axis=-1)

    def depth(self, p, z=0.0):
        p = np.asarray(p, float)
        v = p[..., 0] * self.f[0] + p[..., 1] * self.f[1]
        return v * self.ce - np.asarray(z) * self.se


def lit(colour, normal, sun, alpha=1.0, ambient=0.70):
    """A face's colour under the Sun: an ambient term and the rest of the direct light on it."""
    k = max(0.0, float(np.dot(normal, sun)))
    c = np.array(to_rgb(colour))
    return tuple(np.clip(c * (ambient + (1 - ambient) * k) + 0.05 * k, 0, 1)) + (alpha,)


def sun_vector(alt, az):
    return np.array([math.cos(alt) * math.sin(az), math.cos(alt) * math.cos(az), math.sin(alt)])


def shadow_polygons(ring, h, alt, az):
    """The ground shadow of a prism: the footprint swept along the shadow's direction (quads of each edge, and the
    footprint carried to the shadow's end)."""
    if h <= 0 or alt <= 0:
        return []
    L = h / math.tan(alt)
    s = np.array([-math.sin(az), -math.cos(az)]) * L
    r = np.asarray(ring, float)
    out = [r + s]
    for k in range(len(r)):
        a, b = r[k], r[(k + 1) % len(r)]
        out.append(np.array([a, b, b + s, a + s]))
    return out


def draw_district(ax, inputs, heights, uses, view, sun=None, context=True, shadows=True, passages=None, edge=0.25,
                  extent=None, ground_spaces=True):
    """Draws the district: the context as built in pale stone, cut away in front of the area (as an architect's model
    is cut, so that the area shows), the parcels by their use at the design's heights (a parcel given to open space
    drawn as open ground), the existing open spaces, the links made."""
    sunv = sun_vector(*sun) if sun else np.array([0.3, 0.6, 0.74])
    polys_ground, colours_ground = [], []
    # Open spaces of the city model inside the area, and the parcels given to open space.
    if ground_spaces:
        for s in inputs["open_spaces"]:
            polys_ground.append(view.xy(np.array(s["ring"])))
            colours_ground.append(USE_COLOUR[3])
    for i, p in enumerate(inputs["parcels_ordered"]):
        if uses[i] == 3:
            polys_ground.append(view.xy(np.array(p["footprint"])))
            colours_ground.append(USE_COLOUR[3])
    buildings = []
    if context:
        near = float(np.min(view.depth(np.array(inputs["redesign_quad"])))) - 12
        for b in inputs["context_drawn"]:
            ring = np.array(b["footprint"])
            if float(np.min(view.depth(ring))) >= near:
                buildings.append((ring, b["height_m"], CONTEXT, 0.99))
    for i, p in enumerate(inputs["parcels_ordered"]):
        if uses[i] != 3:
            buildings.append((np.array(p["footprint"]), heights[i], USE_COLOUR[uses[i]], 1.0))
    # Shadows on the ground.
    shadow_polys = []
    if shadows and sun:
        for ring, h, _, _ in buildings:
            for q in shadow_polygons(ring, h, *sun):
                shadow_polys.append(view.xy(q))
    if shadow_polys:
        ax.add_collection(PolyCollection(shadow_polys, facecolors=SHADOW, edgecolors=SHADOW, linewidths=0.4, zorder=1))
    if polys_ground:
        ax.add_collection(PolyCollection(polys_ground, facecolors=colours_ground, edgecolors="none", zorder=2))
        if shadow_polys:
            # The shadows that fall on open ground, in a shaded green, clipped to the open ground.
            from matplotlib.path import Path
            from matplotlib.patches import PathPatch
            clip = PathPatch(Path.make_compound_path(*[Path(np.vstack([q, q[:1]]), closed=True) for q in polys_ground]), fc="none", ec="none", transform=ax.transData)
            ax.add_patch(clip)
            shaded = PolyCollection(shadow_polys, facecolors=SHADOW_GREEN, edgecolors=SHADOW_GREEN, linewidths=0.4, zorder=2.5)
            ax.add_collection(shaded)
            shaded.set_clip_path(clip)
    # Buildings, farthest first.
    order = sorted(range(len(buildings)), key=lambda k: -float(np.min(view.depth(buildings[k][0]))))
    faces, fc, ec, lw = [], [], [], []
    for k in order:
        ring, h, colour, alpha = buildings[k]
        if h <= 0:
            continue
        n = len(ring)
        walls = []
        for j in range(n):
            a, b = ring[j], ring[(j + 1) % n]
            d = b - a
            L = math.hypot(*d)
            if L < 1e-6:
                continue
            normal = np.array([d[1] / L, -d[0] / L])        # outward for a counter-clockwise ring
            if float(np.dot(normal, view.f)) >= 0:
                continue                                   # faces away from the viewer
            depth = float(np.mean(view.depth(np.array([a, b]))))
            quad = np.array([view.xy(a, 0), view.xy(b, 0), view.xy(b, h), view.xy(a, h)])
            walls.append((depth, quad, lit(colour, np.array([normal[0], normal[1], 0.0]), sunv, alpha, 0.70 if alpha == 1.0 else 0.88)))
        walls.sort(key=lambda w: -w[0])
        for _, quad, c in walls:
            faces.append(quad), fc.append(c), ec.append((0.23, 0.18, 0.23, 0.2 if alpha == 1.0 else 0.14)), lw.append(edge)
        faces.append(view.xy(ring, h)), fc.append(lit(colour, np.array([0, 0, 1.0]), sunv, alpha, 0.70 if alpha == 1.0 else 0.92)), ec.append((0.23, 0.18, 0.23, 0.28 if alpha == 1.0 else 0.16)), lw.append(edge)
    ax.add_collection(PolyCollection(faces, facecolors=fc, edgecolors=ec, linewidths=lw, zorder=4))
    # The links made, at ground level, drawn over the buildings they pass beneath (they are arcades).
    for a, b in passages or []:
        P = view.xy(np.array([a, b]))
        ax.plot(P[:, 0], P[:, 1], color=PAPER, lw=3.2, zorder=5, solid_capstyle="round", alpha=0.9)
        ax.plot(P[:, 0], P[:, 1], color=INK, lw=1.5, ls=(0, (2, 1.2)), zorder=6, solid_capstyle="butt")
    if extent is not None:
        ax.set_xlim(extent[0], extent[1]), ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal")
    ax.axis("off")


def prepare_inputs():
    """The inputs with the parcels in the genotype's order, and the context to draw (within a ring about the area)."""
    inputs = load("district-inputs.json")
    cells = load("district-cells.json")
    by_id = {p["id"]: p for p in inputs["parcels"]}
    inputs["parcels_ordered"] = [by_id[i] for i in cells["parcel_order"]]
    quad = np.array(inputs["redesign_quad"])
    c = quad.mean(axis=0)
    R = 1.0 * max(np.hypot(*(quad - c).T))
    inputs["context_drawn"] = [b for b in inputs["context"] if math.hypot(*(np.mean(b["footprint"], axis=0) - c)) < R + 30]
    inputs["cells"] = cells
    return inputs


def design_of(d):
    return d["storeys"], d["uses"], d["links"]


def heights_of(inputs, storeys, uses):
    return [0.0 if u == 3 else s * inputs["storey_m"] for s, u in zip(storeys, uses)]


def as_built(inputs):
    st = [p["floors"] for p in inputs["parcels_ordered"]]
    us = [p["use"] for p in inputs["parcels_ordered"]]
    return st, us


def extent_of(inputs, view, pad=8, ztop=60):
    quad = np.array(inputs["redesign_quad"])
    c = quad.mean(axis=0)
    ring = c + 1.18 * (quad - c)
    P = np.vstack([view.xy(ring, 0), view.xy(ring, ztop)])
    return [P[:, 0].min() - pad, P[:, 0].max() + pad, P[:, 1].min() - pad, P[:, 1].max() + pad]


# ---- the figures ------------------------------------------------------------------------------------------------

def fig_as_built(inputs, summary):
    view = View(40, 38)
    sun = inputs["sun"]["draw"]
    sun = (math.radians(sun["altitude_deg"]), math.radians(sun["azimuth_deg"]))
    fig = plt.figure(figsize=(16, 10), dpi=100)
    ax = fig.add_axes([0.0, 0.05, 0.72, 0.83])
    st, us = as_built(inputs)
    draw_district(ax, inputs, heights_of(inputs, st, us), us, view, sun=sun, extent=extent_of(inputs, view, ztop=70))
    title(fig, "The district as built", "Sixty-one buildings between Swanston and Russell Streets, Flinders Lane and Bourke Street, by their use; shadows at noon on 21 June")
    ab = summary["as_built"]
    lx = 0.74
    fig.text(lx, 0.80, "Uses (CLUE's latest census)", fontsize=13, color=INK)
    for k, u in enumerate([0, 1, 2, 3]):
        y = 0.765 - 0.035 * k
        fig.patches.append(Rectangle((lx, y - 0.006), 0.014, 0.022, transform=fig.transFigure, color=USE_COLOUR[u]))
        n = summary["uses_as_built"][u]
        label = USE_NAME[u] + (f"  ({n:.0f} buildings)" if u < 3 else "  (City Square and small parks)")
        fig.text(lx + 0.022, y, label, fontsize=11, color=INK, va="bottom")
    rows = [("floor area", f"{ab['floor_area_m2']:,.0f} m²"), ("plot ratio", f"{ab['plot_ratio']:.1f}"),
            ("winter sun on public ground", f"{ab['winter_sun_m2h']:,.0f} m²·h"),
            ("heating, 21 June", f"{1e3 * ab['heating_kwh_per_m2_day']:.0f} Wh per m² of floor"),
            ("walk, home to shop + open space", f"{ab['walk_m']:.0f} m"), ("mix of uses H(U | X)", f"{ab['mix_nats']:.3f} nats"),
            ("homes / shops of the floor", f"{100 * ab['dwelling_share']:.1f}% / {100 * ab['shop_share']:.1f}%"),
            ("open space of the land", f"{100 * ab['open_share']:.1f}%")]
    fig.text(lx, 0.57, "Its measures", fontsize=13, color=INK)
    for k, (a, b) in enumerate(rows):
        y = 0.535 - 0.032 * k
        fig.text(lx, y, a, fontsize=11, color=MUTED)
        fig.text(0.975, y, b, fontsize=11, color=INK, ha="right")
    fig.text(lx, 0.24, "Heights: CLUE's floors above ground × 3.2 m where the census\njoins the building, else the city model's height. The context\n"
             "(pale) stands as built; only the sixty-one parcels are redesigned.", fontsize=10, color=MUTED, va="top")
    source(fig, "Drawn from results/district-inputs.json and summary.json. " + PROVENANCE + "; the Sun from the astronomy domain. The context cut away in front of the area.")
    return fig


MEASURES = [  # (fronts.csv column, label, unit, larger is better)
    ("floor_area_m2", "floor area", "thousand m²", True),
    ("winter_sun_m2h", "winter sun", "thousand m²·h", True),
    ("heating_kwh_per_m2_day", "heating", "Wh/m², 21 June", False),
    ("walk_m", "walk to shop + open space", "m", False),
    ("mix_nats", "mix of uses H(U | X)", "nats", True),
]
SCALE = {"floor_area_m2": 1e-3, "winter_sun_m2h": 1e-3, "heating_kwh_per_m2_day": 1e3, "walk_m": 1.0, "mix_nats": 1.0}
SUMMARY_KEY = {"floor_area_m2": "floor_area_m2", "winter_sun_m2h": "winter_sun_m2h", "heating_kwh_per_m2_day": "heating_kwh_per_m2_day", "walk_m": "walk_m", "mix_nats": "mix_nats"}
# The chosen designs, in colours kept apart from the algorithms' (NSGA-III teal, MOEA/D amber, NSGA-II orchid).
CHOSEN_COLOUR = {"most-floor": "#5B8FB9", "most-sun": "#C2665A", "least-heating": "#6B3E75", "shortest-walks": "#C98B6B", "most-mixed": "#B5629E", "balanced": "#3B2F3A"}


def front_rows(alg="NSGA-III", seed="1"):
    return [r for r in read_csv("fronts.csv") if r["algorithm"] == alg and r["seed"] == seed]


def value(row, key):
    return float(row[key]) * SCALE[key]


def centroid(ring):
    r = np.asarray(ring, float)
    x, y = r[:, 0], r[:, 1]
    w = x * np.roll(y, -1) - np.roll(x, -1) * y
    a = w.sum() / 2
    return np.array([((x + np.roll(x, -1)) * w).sum() / (6 * a), ((y + np.roll(y, -1)) * w).sum() / (6 * a)])


def contains(ring, p):
    x, y = p
    c = False
    n = len(ring)
    for i in range(n):
        x1, y1 = ring[i]
        x2, y2 = ring[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            c = not c
    return c


def plan_axes(ax, inputs, pad=24):
    quad = np.array(inputs["redesign_quad"])
    c = quad.mean(axis=0)
    ring = c + 1.12 * (quad - c)
    ax.set_xlim(ring[:, 0].min() - pad, ring[:, 0].max() + pad)
    ax.set_ylim(ring[:, 1].min() - pad, ring[:, 1].max() + pad)
    ax.set_aspect("equal")
    ax.axis("off")


def north_arrow(ax, x, y, size=18):
    ax.annotate("", xy=(x, y + size), xytext=(x, y), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0))
    ax.text(x, y + size + 3, "N", ha="center", va="bottom", fontsize=11, color=INK)


def scale_bar(ax, x, y, length=50):
    ax.plot([x, x + length], [y, y], color=INK, lw=1.2, solid_capstyle="butt")
    ax.text(x + length / 2, y - 3, f"{length} m", ha="center", va="top", fontsize=10, color=INK)


def fig_problem(inputs, summary):
    """The decisions on a plan of the area: the parcels with their storey bounds, the walking network as the domain's
    complex, the candidate passage; and the winter sun on the ground as built."""
    cells = inputs["cells"]
    maps = load("sun-maps.json")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The decisions, and the winter sun as built",
          "Left: each parcel's storeys bound and the walking network as the urban-science domain builds it. Right: hours of sun on the ground, 10:00–14:00 on 21 June")
    ax = fig.add_axes([0.03, 0.08, 0.46, 0.80])
    bx = fig.add_axes([0.51, 0.08, 0.46, 0.80])
    for a in (ax, bx):
        plan_axes(a, inputs)
    # Left: the area's blocks, 2-cells tinted and holes hatched; segments by kind; parcels by bound; the context pale.
    quad_ring = inputs["redesign_quad"]
    for h in cells["holes"]:
        if contains(quad_ring, np.mean(h, axis=0)):
            ax.add_patch(MplPolygon(h, closed=True, fc="#F6EFE8", ec="#D9C7B8", hatch="////", lw=0, zorder=0))
    for b in cells["blocks"]:
        if contains(quad_ring, np.mean(b["ring"], axis=0)):
            ax.add_patch(MplPolygon(b["ring"], closed=True, fc="#EFE2D3", ec="none", zorder=0))
    for b in inputs["context_drawn"]:
        ax.add_patch(MplPolygon(b["footprint"], closed=True, fc=CONTEXT, ec="#D6CEC6", lw=0.3, zorder=1))
    bound_cmap = LinearSegmentedColormap.from_list("bound", ["#F1F5F9", "#C4D6E6", "#8FB2D1", "#5B8FB9"])
    bounds = cells["storeys_bound"]
    vmax = max(bounds)
    for p, bd in zip(inputs["parcels_ordered"], bounds):
        ax.add_patch(MplPolygon(p["footprint"], closed=True, fc=bound_cmap(bd / vmax), ec="#FFFFFF", lw=0.5, zorder=2))
    kind_style = {"street": (SEA, 2.2), "laneway": (SAND, 1.4), "mall": ("#7FA37A", 2.2), "lane": (INK, 1.6)}
    for s in cells["segments"]:
        L = np.array(s["line"])
        c, w = kind_style.get(s["type"], (SEA, 1.6))
        ax.plot(L[:, 0], L[:, 1], color=c, lw=w, zorder=3, solid_capstyle="round", alpha=0.9)
    for k, q in enumerate(cells["links"]):
        P = np.array([q["a"], q["b"]])
        ls = LINK_STYLE[q["kind"]]
        ax.plot(P[:, 0], P[:, 1], color=INK, lw=2.0, ls=ls, zorder=4, solid_capstyle="butt")
        ax.plot(P[0, 0], P[0, 1], "o", ms=3.2, color=INK, zorder=5)
    quad = np.array(inputs["redesign_quad"] + [inputs["redesign_quad"][0]])
    ax.plot(quad[:, 0], quad[:, 1], color=INK, lw=0.6, ls=(0, (1, 2)), zorder=4)
    # A legend for the left panel.
    xs, ys = ax.get_xlim(), ax.get_ylim()
    lx, ly = xs[0] + 8, ys[0] + 8
    ax.add_patch(Rectangle((lx - 5, ly - 5), 172, 84, fc=PAPER, ec=HAIRLINE, lw=0.6, alpha=0.94, zorder=5.5))
    entries = [("streets", SEA, "line"), ("laneways", SAND, "line"), ("pedestrian malls", "#7FA37A", "line"), ("blocks of the area that are 2-cells", "#EFE2D3", "patch"), ("blocks the domain leaves as holes", None, "hatch"),
               ("a blind lane carried through (candidate)", "lane carried through", "link"), ("a passage through a block (candidate)", "passage", "link")]
    for k, (lab, c, kind) in enumerate(entries):
        y = ly + 11 * (len(entries) - 1 - k)
        if kind == "link":
            ax.plot([lx, lx + 14], [y + 3, y + 3], color=INK, lw=2.0, ls=LINK_STYLE[c], zorder=6)
        elif kind == "line":
            ax.plot([lx, lx + 14], [y + 3, y + 3], color=c, lw=2.2, zorder=6)
        elif kind == "patch":
            ax.add_patch(Rectangle((lx, y), 14, 7, fc=c, ec=HAIRLINE, lw=0.5, zorder=6))
        else:
            ax.add_patch(Rectangle((lx, y), 14, 7, fc="#F6EFE8", ec="#D9C7B8", hatch="////", lw=0.5, zorder=6))
        ax.text(lx + 18, y + 3, lab, fontsize=10, color=INK, va="center", zorder=6)
    cax = fig.add_axes([0.345, 0.115, 0.12, 0.014])
    sm = plt.cm.ScalarMappable(cmap=bound_cmap, norm=plt.Normalize(0, vmax))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=9, colors=MUTED, length=2)
    cax.set_title("bound on a parcel's storeys", fontsize=10, color=INK, pad=4)
    # Right: the sun-hours map as built.
    pts = np.array(maps["points"])
    hours = np.array([np.nan if v is None else v for v in maps["hours"]["as-built"]])
    for b in inputs["context_drawn"]:
        bx.add_patch(MplPolygon(b["footprint"], closed=True, fc=CONTEXT, ec="#D6CEC6", lw=0.3, zorder=1))
    st, us = as_built(inputs)
    for p, u in zip(inputs["parcels_ordered"], us):
        if u != 3:
            bx.add_patch(MplPolygon(p["footprint"], closed=True, fc="#E9E2DA", ec="#CFC5BC", lw=0.4, zorder=1))
    ok = np.isfinite(hours)
    d = maps["spacing_m"] / 2
    squares = [np.array([[x - d, y - d], [x + d, y - d], [x + d, y + d], [x - d, y + d]]) for x, y in pts[ok]]
    total = sum(s[2] for s in maps["suns"])
    bx.add_collection(PolyCollection(squares, facecolors=SUN_MAP(hours[ok] / total), edgecolors="none", zorder=2))
    bx.plot(quad[:, 0], quad[:, 1], color=INK, lw=0.6, ls=(0, (1, 2)), zorder=4)
    cax2 = fig.add_axes([0.81, 0.115, 0.13, 0.014])
    sm2 = plt.cm.ScalarMappable(cmap=SUN_MAP, norm=plt.Normalize(0, total))
    cb2 = fig.colorbar(sm2, cax=cax2, orientation="horizontal")
    cb2.outline.set_visible(False)
    cb2.ax.tick_params(labelsize=9, colors=MUTED, length=2)
    cax2.set_title("hours of sun, 10:00–14:00", fontsize=10, color=INK, pad=4)
    xs, ys = bx.get_xlim(), bx.get_ylim()
    north_arrow(bx, xs[0] + 14, ys[0] + 18)
    scale_bar(bx, xs[0] + 34, ys[0] + 14)
    ab = summary["as_built"]
    bx.text(xs[0] + 8, ys[1] - 6, f"public ground sunlit: {ab['winter_sun_m2h']:,.0f} m²·h", fontsize=11, color=INK, va="top")
    source(fig, f"Drawn from results/district-cells.json and sun-maps.json: {len(pts):,} samples of unbuilt ground every {maps['spacing_m']:.0f} m, the study's own shadow model, "
                "the Sun from the astronomy domain.\nStorey bound: the lesser of 48 (the 90th percentile of the Hoddle Grid's completed developments) and 8 √(footprint) / 3.2 m. "
                + PROVENANCE + ".")
    return fig


def fig_parallel(inputs, summary):
    rows = front_rows()
    chosen = load("front-designs.json")["chosen"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The five-objective front", "NSGA-III's final population (seed 1): every line a design, each measure on its own axis, better upward; the district as built dashed")
    ax = fig.add_axes([0.06, 0.19, 0.88, 0.6])
    keys = [m[0] for m in MEASURES]
    V = np.array([[value(r, k) for k in keys] for r in rows])
    ab = summary["as_built"]
    ref = np.array([ab[SUMMARY_KEY[k]] * SCALE[k] for k in keys])
    allv = np.vstack([V, ref])
    lo, hi = allv.min(axis=0), allv.max(axis=0)
    span = np.where(hi > lo, hi - lo, 1)

    def y_of(v):
        y = (v - lo) / span
        better_up = np.array([m[3] for m in MEASURES])
        return np.where(better_up, y, 1 - y)

    xs = np.arange(len(keys))
    Y = y_of(V)
    for y in Y:
        ax.plot(xs, y, color=HAIRLINE, lw=0.9, alpha=0.9, zorder=1)
    for c in chosen:
        y = Y[int(c["index"])]
        ax.plot(xs, y, color=CHOSEN_COLOUR[c["key"]], lw=2.6, zorder=3, label=c["title"])
    ax.plot(xs, y_of(ref), color=INK, lw=1.6, ls=(0, (4, 2)), zorder=4, label="the district as built")
    for k, (key, lab, unit, up) in enumerate(MEASURES):
        ax.axvline(k, color=MUTED, lw=0.8, zorder=0)
        top, bottom = (hi[k], lo[k]) if up else (lo[k], hi[k])
        fmt = (lambda v: f"{v:,.0f}") if key in ("floor_area_m2", "winter_sun_m2h", "walk_m") else (lambda v: f"{v:.1f}" if key == "heating_kwh_per_m2_day" else f"{v:.2f}")
        ax.text(k, 1.04, fmt(top), ha="center", va="bottom", fontsize=10.5, color=INK)
        ax.text(k, -0.04, fmt(bottom), ha="center", va="top", fontsize=10.5, color=INK)
        ax.text(k, -0.1, f"{lab}\n({unit})", ha="center", va="top", fontsize=12, color=INK)
    ax.set_xlim(-0.35, len(keys) - 0.65)
    ax.set_ylim(-0.02, 1.02)
    ax.axis("off")
    fig.legend(*ax.get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(0.035, 0.885), ncol=7, fontsize=11, handlelength=2.2, columnspacing=1.6)
    source(fig, f"Drawn from results/fronts.csv and front-designs.json: {len(rows)} feasible non-dominated designs after "
                f"{summary['generations']:.0f} generations of {summary['population']:.0f}. Heating in Wh per m² of floor on 21 June; walking from homes, weighted by their floor.")
    return fig


def fig_pairs(inputs, summary):
    rows = front_rows()
    rows_m = front_rows("MOEA/D")
    chosen = load("front-designs.json")["chosen"]
    keys = [m[0] for m in MEASURES]
    V = np.array([[value(r, k) for k in keys] for r in rows])
    W = np.array([[value(r, k) for k in keys] for r in rows_m]) if rows_m else None
    ab = summary["as_built"]
    ref = np.array([ab[SUMMARY_KEY[k]] * SCALE[k] for k in keys])
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Trade-offs, two measures at a time", "The five-objective fronts of NSGA-III (teal) and MOEA/D (amber), seed 1; the chosen designs ringed; the district as built a star")
    n = len(keys)
    gx0, gy0, gw, gh = 0.07, 0.13, 0.64, 0.74
    for i in range(1, n):
        for j in range(i):
            ax = fig.add_axes([gx0 + j * gw / (n - 1), gy0 + (n - 1 - i) * gh / (n - 1), gw / (n - 1) * 0.9, gh / (n - 1) * 0.86])
            if W is not None:
                ax.scatter(W[:, j], W[:, i], s=9, color=ALGORITHM["MOEA/D"], alpha=0.5, lw=0, zorder=1)
            ax.scatter(V[:, j], V[:, i], s=9, color=ALGORITHM["NSGA-III"], alpha=0.65, lw=0, zorder=2)
            for c in chosen:
                v = V[int(c["index"])]
                ax.scatter([v[j]], [v[i]], s=46, facecolors="none", edgecolors=CHOSEN_COLOUR[c["key"]], lw=1.6, zorder=3)
            ax.scatter([ref[j]], [ref[i]], s=90, marker="*", color=INK, zorder=4)
            ax.tick_params(labelsize=8.5)
            ax.grid(True, lw=0.4)
            if i == n - 1:
                ax.set_xlabel(MEASURES[j][1] + f"\n({MEASURES[j][2]})", fontsize=10.5)
            else:
                ax.set_xticklabels([])
            if j == 0:
                ax.set_ylabel(MEASURES[i][1] + f"\n({MEASURES[i][2]})", fontsize=10.5)
            else:
                ax.set_yticklabels([])
    # The chosen designs' key.
    key = fig.add_axes([0.75, 0.3, 0.24, 0.57])
    key.set_xlim(0, 1), key.set_ylim(0, 1), key.axis("off")
    key.text(0, 0.98, "The chosen designs", fontsize=13, color=INK, va="top")
    for k, c in enumerate(chosen):
        y = 0.86 - 0.105 * k
        key.scatter([0.03], [y], s=60, facecolors="none", edgecolors=CHOSEN_COLOUR[c["key"]], lw=1.8)
        r = rows[int(c["index"])]
        key.text(0.08, y, c["title"], fontsize=11.5, color=INK, va="center")
        key.text(0.08, y - 0.045, f"{value(r, 'floor_area_m2'):.0f}k m², sun {value(r, 'winter_sun_m2h'):.1f}k m²·h, walk {value(r, 'walk_m'):.0f} m",
                 fontsize=9.5, color=MUTED, va="center")
    y = 0.86 - 0.105 * len(chosen)
    key.scatter([0.03], [y], s=110, marker="*", color=INK)
    key.text(0.08, y, "the district as built", fontsize=11.5, color=INK, va="center")
    key.scatter([0.03], [y - 0.08], s=14, color=ALGORITHM["NSGA-III"]), key.text(0.08, y - 0.08, "NSGA-III's front", fontsize=10.5, color=INK, va="center")
    key.scatter([0.03], [y - 0.13], s=14, color=ALGORITHM["MOEA/D"]), key.text(0.08, y - 0.13, "MOEA/D's front", fontsize=10.5, color=INK, va="center")
    fig.text(0.75, 0.27, "Each panel plots two of the five measures for\nevery design of the two fronts. A design that\nlooks dominated in a panel is not: it is better\n"
             "in a measure the panel does not show.", fontsize=10.5, color=MUTED, va="top")
    source(fig, "Drawn from results/fronts.csv and front-designs.json. Larger is better for floor area, winter sun and mix; smaller for heating and walking.")
    return fig


def measure_line(m):
    return (f"{m['floor_area_m2'] / 1e3:.0f}k m² floor · {m['winter_sun_m2h'] / 1e3:.1f}k m²·h sun · {1e3 * m['heating_kwh_per_m2_day']:.0f} Wh/m² heat\n"
            f"walk {m['walk_m']:.0f} m · mix {m['mix_nats']:.2f} · open {100 * m['open_share']:.0f}% · homes {100 * m['dwelling_share']:.0f}%")


def fig_chosen(inputs, summary):
    designs = {d["key"]: d for d in load("designs.json")["designs"]}
    fd = load("front-designs.json")
    view = View(40, 38)
    sun = inputs["sun"]["draw"]
    sun = (math.radians(sun["altitude_deg"]), math.radians(sun["azimuth_deg"]))
    keys = [c["key"] for c in fd["chosen"]]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Six designs from the front", "Each the best of the front in one measure, and one balanced among all five; shadows at noon on 21 June")
    ext = extent_of(inputs, view, pad=4, ztop=120)
    passages = inputs["cells"]["links"]
    for k, key in enumerate(keys):
        r, c = divmod(k, 3)
        ax = fig.add_axes([0.01 + c * 0.33, 0.475 - r * 0.395, 0.32, 0.325])
        idx = next(int(x["index"]) for x in fd["chosen"] if x["key"] == key)
        st, us, ps = design_of(fd["nsga3"][idx]["design"])
        opened = [(q["a"], q["b"]) for q, o in zip(passages, ps) if o]
        draw_district(ax, inputs, heights_of(inputs, st, us), us, view, sun=sun, extent=ext, passages=opened, edge=0.18)
        d = designs[key]
        fig.text(0.02 + c * 0.33, 0.875 - r * 0.395, d["title"][0].upper() + d["title"][1:], fontsize=14, color=CHOSEN_COLOUR[key], va="top")
        fig.text(0.02 + c * 0.33, 0.848 - r * 0.395, measure_line(d["measures"]), fontsize=9.5, color=MUTED, va="top")
    lx = 0.02
    for k, u in enumerate([0, 1, 2, 3]):
        x = lx + 0.16 * k
        fig.patches.append(Rectangle((x, 0.052), 0.012, 0.018, transform=fig.transFigure, color=USE_COLOUR[u]))
        fig.text(x + 0.017, 0.061, USE_NAME[u], fontsize=10.5, color=INK, va="center")
    source(fig, "Drawn from results/designs.json and front-designs.json (NSGA-III, seed 1). The context pale and as built, cut away in front; a link made dashed. One view and scale in every panel.")
    return fig


def fig_two(inputs, summary):
    rows2 = [r for r in read_csv("fronts.csv") if r["algorithm"] == "NSGA-II"]
    rows3 = front_rows()
    ab = summary["as_built"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Two objectives: floor against winter sun", "NSGA-II asked only these two (three seeds), beside the five-objective front of NSGA-III seen in the same plane")
    ax = fig.add_axes([0.07, 0.12, 0.55, 0.72])
    x3 = [value(r, "floor_area_m2") for r in rows3]
    y3 = [value(r, "winter_sun_m2h") for r in rows3]
    ax.scatter(x3, y3, s=16, color=ALGORITHM["NSGA-III"], alpha=0.5, lw=0, label="NSGA-III, five objectives (seed 1)", zorder=1)
    shades = {"1": "#8E3B78", "2": "#B5629E", "3": "#D79BC8"}
    for seed in ("1", "2", "3"):
        rs = sorted([r for r in rows2 if r["seed"] == seed], key=lambda r: float(r["floor_area_m2"]))
        ax.plot([value(r, "floor_area_m2") for r in rs], [value(r, "winter_sun_m2h") for r in rs], color=shades[seed], lw=1.6, marker="o", ms=3.2,
                label=f"NSGA-II, two objectives (seed {seed})", zorder=3)
    ax.scatter([ab["floor_area_m2"] / 1e3], [ab["winter_sun_m2h"] / 1e3], s=200, marker="*", color=INK, zorder=5, label="the district as built")
    ax.set_xlabel("floor area (thousand m²)")
    ax.set_ylabel("winter sun on the public ground (thousand m²·h)")
    ax.legend(loc="upper right", fontsize=10.5)
    h2 = summary["hypervolume_two_objectives"]
    ax.text(0.02, 0.04, "hypervolume in this plane (normalised, reference 1.1):\n"
            f"NSGA-II {np.mean(h2['NSGA-II']):.3f} ± {np.std(h2['NSGA-II']):.3f};  NSGA-III's projection {np.mean(h2['NSGA-III projected']):.3f} ± {np.std(h2['NSGA-III projected']):.3f}",
            transform=ax.transAxes, fontsize=10.5, color=INK)
    # Thumbnails of three NSGA-II designs along its front (seed 1).
    fd = load("front-designs.json")
    n2 = sorted(fd["nsga2"], key=lambda d: d["objectives"][0])
    picks = [n2[0], n2[len(n2) // 2], n2[-1]] if len(n2) >= 3 else n2
    view = View(40, 38)
    ext = extent_of(inputs, view, pad=4, ztop=120)
    passages = inputs["cells"]["links"]
    for k, d in enumerate(picks):
        bx = fig.add_axes([0.64, 0.585 - 0.262 * k, 0.35, 0.235])
        st, us, ps = design_of(d["design"])
        opened = [(q["a"], q["b"]) for q, o in zip(passages, ps) if o]
        draw_district(bx, inputs, heights_of(inputs, st, us), us, view, sun=None, shadows=False, extent=ext, passages=opened, edge=0.12)
        floor, sunv = -d["objectives"][0] * ab["floor_area_m2"] / 1e3, -d["objectives"][1] * ab["winter_sun_m2h"] / 1e3
        fig.text(0.645, 0.585 - 0.262 * k + 0.245, f"{k + 1}.  {floor:.0f}k m² of floor, {sunv:.1f}k m²·h of winter sun", fontsize=11, color=INK, va="bottom")
        ax.scatter([floor], [sunv], s=90, facecolors="none", edgecolors=INK, lw=1.2, zorder=6)
        ax.annotate(str(k + 1), (floor, sunv), xytext=(6, 6), textcoords="offset points", fontsize=10, color=INK)
    source(fig, "Drawn from results/fronts.csv, front-designs.json and summary.json. Same constraints in both searches: open space ≥ 15% of the land; homes ≥ 20% and shops 5–20% of the floor.")
    return fig


def fig_hypervolume(inputs, summary):
    rows = read_csv("hypervolume.csv")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "How far each search got", "Hypervolume of the feasible front found, five objectives normalised by the union of the final fronts, reference point 1.1; three seeds each")
    ax = fig.add_axes([0.07, 0.12, 0.58, 0.72])
    colours = {"NSGA-III": ALGORITHM["NSGA-III"], "MOEA/D": ALGORITHM["MOEA/D"], "random": ALGORITHM["random"]}
    for alg in ("NSGA-III", "MOEA/D", "random"):
        for seed in ("1", "2", "3"):
            rs = [r for r in rows if r["algorithm"] == alg and r["seed"] == seed]
            ax.plot([float(r["evaluations"]) / 1e3 for r in rs], [float(r["hypervolume"]) for r in rs], color=colours[alg], lw=1.8, marker="o", ms=3.5,
                    label=alg if seed == "1" else None, alpha=0.9)
    ax.set_xlabel("evaluations (thousands)")
    ax.set_ylabel("hypervolume (normalised)")
    ax.legend(loc="lower right")
    bx = fig.add_axes([0.72, 0.12, 0.25, 0.72])
    final = summary["hypervolume_final"]
    for k, alg in enumerate(("NSGA-III", "MOEA/D", "random")):
        v = final.get(alg, [])
        bx.scatter([k] * len(v), v, s=60, color=colours[alg], zorder=3)
        if v:
            bx.plot([k - 0.25, k + 0.25], [np.mean(v)] * 2, color=colours[alg], lw=2)
            bx.text(k, max(v) + 0.004, f"{np.mean(v):.3f}", ha="center", va="bottom", fontsize=11, color=INK)
    bx.set_xticks(range(3))
    bx.set_xticklabels(["NSGA-III", "MOEA/D", "random"])
    bx.set_xlim(-0.6, 2.6)
    bx.set_title("at the end", fontsize=13, color=INK, loc="left")
    runs = summary["runs"]
    ev = runs[0]["evaluations"]
    source(fig, f"Drawn from results/hypervolume.csv and summary.json: {ev:,.0f} evaluations a run; random search's front is its whole feasible non-dominated archive, the others' their final populations' fronts.")
    return fig


def fig_atlas(inputs, summary):
    me = load("map-elites.json")
    g = me["grid"]
    ni, nj = int(g["resolution"][0]), int(g["resolution"][1])
    pr, op = g["plot_ratio"], g["open_share"]
    cells = {(int(c["i"]), int(c["j"])): c for c in me["cells"]}
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "A repertoire of districts", "MAP-Elites: in each cell of plot ratio and open share, the design that keeps the most winter sun on the public ground")
    gx0, gy0, gw, gh = 0.08, 0.1, 0.84, 0.76
    view = View(40, 42)
    ext = extent_of(inputs, view, pad=-14, ztop=95)
    quals = [-c["quality"] for c in me["cells"] if c["feasible"]]
    qn = plt.Normalize(min(quals), max(quals))
    tint = LinearSegmentedColormap.from_list("tint", ["#F6EEE8", "#F3E1C9", "#EBC79B"])
    ax_bg = fig.add_axes([gx0, gy0, gw, gh])
    ax_bg.set_xlim(pr[0], pr[1]), ax_bg.set_ylim(op[0], op[1])
    for spine in ("top", "right"):
        ax_bg.spines[spine].set_visible(False)
    ax_bg.grid(False)
    ax_bg.set_xlabel("plot ratio (floor area / land)")
    ax_bg.set_ylabel("open space (share of the land)")
    dx, dy = (pr[1] - pr[0]) / ni, (op[1] - op[0]) / nj
    for (i, j), c in cells.items():
        ax_bg.add_patch(Rectangle((pr[0] + i * dx, op[0] + j * dy), dx, dy, fc=tint(qn(-c["quality"])) if c["feasible"] else "#EEE8E4", ec="#FFFFFF", lw=1.0))
    for (i, j), c in cells.items():
        x0 = gx0 + gw * i / ni
        y0 = gy0 + gh * j / nj
        ax = fig.add_axes([x0 + 0.002, y0 + 0.004, gw / ni - 0.004, gh / nj - 0.008])
        st, us, ps = design_of(c["design"])
        draw_district(ax, inputs, heights_of(inputs, st, us), us, view, sun=None, shadows=False, context=False, extent=ext, edge=0.05)
    ab = summary["as_built"]
    ov = fig.add_axes([gx0, gy0, gw, gh])
    ov.set_xlim(pr[0], pr[1]), ov.set_ylim(op[0], op[1])
    ov.axis("off")
    ov.patch.set_alpha(0)
    ov.scatter([ab["plot_ratio"]], [ab["open_share"]], s=200, marker="*", color=INK, edgecolors=PAPER, linewidths=1.0, zorder=10, clip_on=False)
    ov.annotate("as built", (ab["plot_ratio"], ab["open_share"]), xytext=(9, 6), textcoords="offset points", fontsize=11, color=INK,
                bbox=dict(boxstyle="round,pad=0.2", fc=PAPER, ec="none", alpha=0.9))
    source(fig, f"Drawn from results/map-elites.json: {me['budget']:,.0f} evaluations, {100 * me['coverage']:.0f}% of the {ni * nj} cells filled; a cell's tint its elite's winter sun "
                "(darker is sunnier). Homes ≥ 20% and shops 5–20% of the floor; context not drawn.")
    return fig


def fig_repertoire_map(inputs, summary):
    me = load("map-elites.json")
    prog = read_csv("map-elites-progress.csv")
    g = me["grid"]
    ni, nj = int(g["resolution"][0]), int(g["resolution"][1])
    pr, op = g["plot_ratio"], g["open_share"]
    ab = summary["as_built"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "What the repertoire says", "Left: the most winter sun the search found at each plot ratio and open share, relative to the district as built. Right: the search's progress")
    ax = fig.add_axes([0.07, 0.12, 0.5, 0.72])
    Z = np.full((nj, ni), np.nan)
    for c in me["cells"]:
        if c["feasible"]:
            Z[int(c["j"]), int(c["i"])] = -c["quality"]
    from matplotlib.colors import TwoSlopeNorm
    finite = Z[np.isfinite(Z)]
    norm = TwoSlopeNorm(vcenter=1.0, vmin=min(0.99, float(finite.min())), vmax=max(1.01, float(finite.max())))
    im = ax.imshow(Z, origin="lower", extent=[pr[0], pr[1], op[0], op[1]], aspect="auto", cmap=SUN_MAP, norm=norm, interpolation="nearest")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.outline.set_visible(False)
    cb.set_label("winter sun, relative to the district as built", fontsize=11)
    ax.scatter([ab["plot_ratio"]], [ab["open_share"]], s=180, marker="*", color=INK, zorder=5)
    ax.annotate("as built", (ab["plot_ratio"], ab["open_share"]), xytext=(8, -4), textcoords="offset points", fontsize=11, color=INK)
    rows = front_rows()
    ax.scatter([float(r["plot_ratio"]) for r in rows], [float(r["open_share"]) for r in rows], s=10, color="#FFFFFF", edgecolors=INK, lw=0.5, zorder=4, label="NSGA-III's front")
    ax.set_xlim(pr[0], pr[1]), ax.set_ylim(op[0], op[1])
    beyond = sum(1 for r in rows if not (pr[0] <= float(r["plot_ratio"]) <= pr[1] and op[0] <= float(r["open_share"]) <= op[1]))
    ax.set_xlabel("plot ratio (floor area / land)")
    ax.set_ylabel("open space (share of the land)")
    ax.legend(loc="upper right", fontsize=10.5)
    bx = fig.add_axes([0.68, 0.52, 0.28, 0.32])
    E = [float(r["evaluations"]) / 1e3 for r in prog]
    bx.plot(E, [100 * float(r["coverage"]) for r in prog], color=ALGORITHM["MAP-Elites"], lw=2)
    bx.set_ylabel("cells filled (%)")
    bx.set_xlabel("evaluations (thousands)")
    cx = fig.add_axes([0.68, 0.12, 0.28, 0.28])
    cx.plot(E, [float(r["qd_score"]) for r in prog], color="#4E8F8B", lw=2)
    cx.set_ylabel("QD-score")
    cx.set_xlabel("evaluations (thousands)")
    source(fig, f"Drawn from results/map-elites.json, map-elites-progress.csv and fronts.csv ({beyond} of NSGA-III's {len(rows)} designs lie beyond the grid). "
                "QD-score: the sum over filled feasible cells of the elite's winter sun relative to the district as built.")
    return fig


def fig_where(inputs, summary):
    fd = load("front-designs.json")
    designs = [d["design"] for d in fd["nsga3"]]
    n = len(inputs["parcels_ordered"])
    S = np.array([d["storeys"] for d in designs], float)
    U = np.array([d["uses"] for d in designs])
    H = np.where(U == 3, np.nan, S)
    med = np.nanmedian(np.where(np.isnan(H), np.nan, H), axis=0)
    share = [(U[:, i][:, None] == np.arange(4)[None, :]).mean(axis=0) for i in range(n)]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Where the front puts height, and what it gives to the public",
          "Over NSGA-III's front (seed 1): left, each parcel's median storeys when built; right, its most frequent use and how often the front gives it to open space")
    ax = fig.add_axes([0.03, 0.15, 0.46, 0.74])
    bx = fig.add_axes([0.51, 0.15, 0.46, 0.74])
    for a in (ax, bx):
        plan_axes(a, inputs)
        for b in inputs["context_drawn"]:
            a.add_patch(MplPolygon(b["footprint"], closed=True, fc=CONTEXT, ec="#D6CEC6", lw=0.3, zorder=1))
        for s in inputs["open_spaces"]:
            a.add_patch(MplPolygon(s["ring"], closed=True, fc=USE_COLOUR[3], ec="none", zorder=1))
    hcm = LinearSegmentedColormap.from_list("height", ["#F3EDE6", "#E7C9A6", "#C98B6B", "#8A4E3B", "#3B2F3A"])
    vmax = np.nanmax(med) if np.isfinite(med).any() else 1
    for p, m in zip(inputs["parcels_ordered"], med):
        ax.add_patch(MplPolygon(p["footprint"], closed=True, fc=hcm(m / vmax) if np.isfinite(m) else USE_COLOUR[3], ec="#FFFFFF", lw=0.5, zorder=2))
    cax = fig.add_axes([0.16, 0.085, 0.2, 0.016])
    sm = plt.cm.ScalarMappable(cmap=hcm, norm=plt.Normalize(0, vmax))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=9, colors=MUTED, length=2)
    cax.set_title("median storeys when built", fontsize=10, color=INK, pad=4)
    for p, sh in zip(inputs["parcels_ordered"], share):
        u = int(np.argmax(sh))
        bx.add_patch(MplPolygon(p["footprint"], closed=True, fc=USE_COLOUR[u], ec="#FFFFFF", lw=0.5, alpha=0.35 + 0.65 * float(sh[u]), zorder=2))
        if sh[3] >= 0.35 and p["area_m2"] >= 250:
            c = centroid(p["footprint"])
            bx.text(c[0], c[1], f"{100 * sh[3]:.0f}%", ha="center", va="center", fontsize=9, color=INK, zorder=3)
    for k, u in enumerate([0, 1, 2, 3]):
        x = 0.53 + 0.095 * k
        fig.patches.append(Rectangle((x, 0.083), 0.012, 0.018, transform=fig.transFigure, color=USE_COLOUR[u], zorder=6))
        fig.text(x + 0.017, 0.092, USE_NAME[u].replace(" and other workplaces", "").replace(", food and entertainment", ""), fontsize=10.5, color=INK, va="center", zorder=6)
    fig.text(0.53, 0.115, "The most frequent use over the front, stronger where it agrees; the figures: how often it gives the parcel to open space",
             fontsize=10, color=MUTED, va="bottom")
    xs, ys = ax.get_xlim(), ax.get_ylim()
    north_arrow(ax, xs[0] + 14, ys[1] - 40)
    scale_bar(ax, xs[0] + 34, ys[1] - 36)
    source(fig, f"Drawn from results/front-designs.json: {len(designs)} designs of NSGA-III's final front (seed 1).")
    return fig


FIGURES = {
    50: ("search-50-district-as-built", fig_as_built),
    51: ("search-51-decisions-and-sun", fig_problem),
    52: ("search-52-front-parallel", fig_parallel),
    53: ("search-53-front-pairs", fig_pairs),
    54: ("search-54-six-designs", fig_chosen),
    55: ("search-55-two-objectives", fig_two),
    56: ("search-56-hypervolume", fig_hypervolume),
    57: ("search-57-repertoire-atlas", fig_atlas),
    58: ("search-58-repertoire-map", fig_repertoire_map),
    59: ("search-59-where-height-goes", fig_where),
}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--results", default=str(RESULTS), help="the results directory (results/ by default)")
    args = ap.parse_args()
    RESULTS = pathlib.Path(args.results)
    only = {int(x) for x in args.only.split(",") if x}
    inputs = prepare_inputs()
    summary = load("summary.json")
    for k, (name, f) in FIGURES.items():
        if only and k not in only:
            continue
        fig = f(inputs, summary)
        print(save(fig, name))
