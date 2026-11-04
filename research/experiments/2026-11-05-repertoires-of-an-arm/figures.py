"""Figures of the study 'Repertoires of an arm', drawn from results/ only.

    .venv/Scripts/python.exe research/experiments/2026-11-05-repertoires-of-an-arm/figures.py

Writes assets/gallery/search-24 ... search-26.
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
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm  # noqa: E402
from matplotlib.patches import PathPatch  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
RES = 50
JOINTS = 12

COLOUR = {
    "MAP-Elites, iso+line": ALGORITHM["MAP-Elites"],
    "MAP-Elites, Gaussian 0.03": "#E7B9A1",
    "MAP-Elites, Gaussian 0.1": "#D9738B",
    "MAP-Elites, Gaussian 0.3": "#A8577E",
    "CVT-MAP-Elites, iso+line": ALGORITHM["CVT-MAP-Elites"],
    "novelty search, Gaussian 0.1": ALGORITHM["novelty"],
    "random search": ALGORITHM["random"],
    "CMA-ES, IPOP": "#6C8EBF",
}
ORDER = list(COLOUR)
# The angle spread of an elite, its standard deviation √var, on a log scale shared by every map.
SD = LogNorm(vmin=0.01, vmax=1.5)
MAP = LinearSegmentedColormap.from_list("spread", ["#2E1A2D", "#6B2F5E", "#A8577E", "#CF97B5", "#E6C3D3", "#F1DDE5"])


def rows(name: str) -> list[dict]:
    with open(RESULTS / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def provenance() -> str:
    v = json.loads((RESULTS / "versions.json").read_text(encoding="utf-8"))
    return (f"heuresis {v['heuresis']} (foundations/heuresis), root {v['root']}; heuresis::benchmarks::planar_arm with 12 joints in [−π, π]; "
            "research/experiments/2026-11-05-repertoires-of-an-arm")


def disc(ax, lw=0.9, mask=False):
    if mask:
        # paper over everything outside the disc: a square with a circular hole
        t = np.linspace(0, 2 * np.pi, 400)
        square = [(-0.05, -0.05), (1.05, -0.05), (1.05, 1.05), (-0.05, 1.05), (-0.05, -0.05)]
        circle = list(zip(0.5 + 0.5 * np.cos(t[::-1]), 0.5 + 0.5 * np.sin(t[::-1])))
        verts = square + circle + [circle[0]]
        codes = [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY] + [Path.MOVETO] + [Path.LINETO] * (len(circle) - 1) + [Path.CLOSEPOLY]
        ax.add_patch(PathPatch(Path(verts, codes), facecolor=PAPER, edgecolor="none", zorder=3))
    ax.add_patch(Circle((0.5, 0.5), 0.5, fill=False, edgecolor=MUTED, lw=lw, zorder=4))
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.grid(False)


def grid_map(ax, elites):
    img = np.full((RES, RES), np.nan)
    for e in elites:
        img[int(e["j"]), int(e["i"])] = np.sqrt(max(float(e["variance"]), 1e-12))
    ax.imshow(img, origin="lower", extent=(0, 1, 0, 1), cmap=MAP, norm=SD, interpolation="nearest", zorder=2)
    disc(ax)


def final_medians():
    F = defaultdict(list)
    for r in rows("final.csv"):
        F[r["algorithm"]].append(r)
    return F


# ------------------------------------------------------------------------------------------------------------------
# search-24: maps of the elites over the descriptor space.

def figure_maps():
    E = defaultdict(list)
    for r in rows("elites.csv"):
        E[r["algorithm"]].append(r)
    F = final_medians()
    ref = json.loads((RESULTS / "reference.json").read_text(encoding="utf-8"))
    fig, axes = figure(2, 3)
    fig.subplots_adjust(left=0.03, right=0.88, top=0.80, bottom=0.05, wspace=0.05, hspace=0.24)
    title(fig, "Six searches illuminating the reach of a twelve-jointed arm",
          "Where the tip can be put (the disc), and the evenest arm each search found for each place: the colour is the spread of its joint angles")
    panels = ["MAP-Elites, iso+line", "MAP-Elites, Gaussian 0.1", "CVT-MAP-Elites, iso+line", "novelty search, Gaussian 0.1", "random search", "CMA-ES, IPOP"]
    for ax, alg in zip(axes.flat, panels):
        if alg.startswith("CVT"):
            cvt = rows("cvt.csv")
            pts = np.array([[float(c["x"]), float(c["y"])] for c in cvt])
            polys, vals = voronoi_cells(pts)
            filled = np.array([int(c["filled"]) for c in cvt]).astype(bool)
            sd = np.sqrt(np.maximum([float(c["variance"]) for c in cvt], 1e-12))
            colours = [MAP(SD(s)) if f else PANEL for s, f in zip(sd, filled)]
            ax.add_collection(PolyCollection([polys[i] for i in range(len(polys))], facecolors=colours, edgecolors=PAPER, linewidths=0.3, zorder=2))
            disc(ax, mask=True)
            own = np.median([float(r["own_coverage"]) for r in F[alg]])
            note = f"{100 * own:.1f}% of its {ref['cvt_sites']} cells"
        else:
            grid_map(ax, E[alg])
            note = ""
        cov = np.median([float(r["coverage"]) for r in F[alg]])
        qd = np.median([float(r["qd"]) for r in F[alg]])
        ax.text(0.0, 1.085, alg, transform=ax.transAxes, fontsize=14, color=COLOUR[alg], va="bottom", ha="left")
        text = f"coverage {100 * cov:.1f}% · QD-score {qd:.3f}"
        if note:
            text += " · " + note
        ax.text(0.0, 1.03, text, transform=ax.transAxes, fontsize=10, color=MUTED, va="bottom", ha="left")
    cax = fig.add_axes([0.9, 0.18, 0.014, 0.55])
    sm = plt.cm.ScalarMappable(norm=SD, cmap=MAP)
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("standard deviation of the 12 joint angles (rad)")
    cb.outline.set_edgecolor(HAIRLINE)
    source(fig, provenance().split("; research")[0] + "; 100 000 evaluations each; 50 × 50 cells; medians over 15 seeds, the run shown the median seed's")
    return save(fig, "search-24-repertoire-maps")


def voronoi_cells(pts):
    """The Voronoi cells of the sites, clipped to the unit square (by a dense nearest-site raster turned into polygons)."""
    from scipy.spatial import Voronoi

    far = np.array([[-10, -10], [-10, 11], [11, -10], [11, 11]])
    vor = Voronoi(np.vstack([pts, far]))
    polys = []
    for i in range(len(pts)):
        region = vor.regions[vor.point_region[i]]
        poly = np.array([vor.vertices[v] for v in region if v >= 0])
        polys.append(clip_square(poly))
    return polys, None


def clip_square(poly):
    def clip(p, inside, intersect):
        out = []
        for a, b in zip(p, np.roll(p, -1, axis=0)):
            if inside(b):
                if not inside(a):
                    out.append(intersect(a, b))
                out.append(b)
            elif inside(a):
                out.append(intersect(a, b))
        return np.array(out) if out else np.zeros((0, 2))

    for axis, bound, keep_less in [(0, 0.0, False), (0, 1.0, True), (1, 0.0, False), (1, 1.0, True)]:
        if len(poly) == 0:
            break
        inside = (lambda q, a=axis, b=bound, l=keep_less: q[a] <= b if l else q[a] >= b)

        def intersect(a, b, ax_=axis, bd=bound):
            t = (bd - a[ax_]) / (b[ax_] - a[ax_])
            return a + t * (b - a)
        poly = clip(poly, inside, intersect)
    return poly


# ------------------------------------------------------------------------------------------------------------------
# search-25: arms for chosen places.

def joints(theta):
    a = np.cumsum(theta)
    L = 1.0 / JOINTS
    x = np.concatenate([[0], np.cumsum(L * np.cos(a))])
    y = np.concatenate([[0], np.cumsum(L * np.sin(a))])
    return (x + 1) / 2, (y + 1) / 2  # descriptor coordinates


def figure_arms():
    A = rows("arms.csv")
    E = defaultdict(list)
    for r in rows("elites.csv"):
        E[r["algorithm"]].append(r)
    radii = [0.06, 0.16, 0.26, 0.36, 0.46]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "A repertoire of arms: one for each place the tip can reach",
          "Elites of the median MAP-Elites run (iso+line, 100 000 evaluations) for places on rings about the base, drawn joint by joint")
    ax = fig.add_axes([0.03, 0.06, 0.50, 0.80])
    grid_map(ax, E["MAP-Elites, iso+line"])
    ax.images[0].set_alpha(0.25)
    for r in A:
        if int(r["ring"]) not in (2, 4):
            continue
        th = np.array([float(r[f"theta{q + 1}"]) for q in range(JOINTS)])
        x, y = joints(th)
        ax.plot(x, y, color="#6B3E75", lw=0.75, alpha=0.9, zorder=5)
        ax.scatter(x[-1:], y[-1:], s=10, color="#6B3E75", zorder=6, linewidths=0)
    for rad in (0.26, 0.46):
        ax.add_patch(Circle((0.5, 0.5), rad, fill=False, edgecolor=MUTED, lw=0.6, ls=(0, (2, 3)), zorder=4))
    ax.scatter([0.5], [0.5], s=30, color=INK, zorder=7)
    ax.text(0.0, 1.0, "32 places on two rings (dotted), radii 0.26 and 0.46 of the unit square; the arm's reach is 0.5", transform=ax.transAxes, fontsize=10.5, color=MUTED, va="bottom")
    # Small multiples on a common scale: one arm per panel, for a place on each ring.
    picks = []
    for ring in range(5):
        for d in (1, 9, 5):
            for r in A:
                if int(r["ring"]) == ring and int(r["direction"]) == d:
                    picks.append(r)
    picks = picks[:15]
    for idx, r in enumerate(picks):
        col, row = divmod(idx, 3)
        size = 0.083
        sub = fig.add_axes([0.555 + col * 0.088, 0.62 - row * 0.255, size, size * 1.6])
        th = np.array([float(r[f"theta{q + 1}"]) for q in range(JOINTS)])
        x, y = joints(th)
        sub.add_patch(Circle((0.5, 0.5), 0.5, fill=False, edgecolor=HAIRLINE, lw=0.6))
        sub.add_patch(Circle((0.5, 0.5), radii[int(r["ring"])], fill=False, edgecolor=MUTED, lw=0.5, ls=(0, (2, 3))))
        sub.plot(x, y, color="#6B3E75", lw=1.1)
        sub.scatter(x[1:-1], y[1:-1], s=5, color="#6B3E75", zorder=3, linewidths=0)
        sub.scatter([0.5], [0.5], s=12, color=INK, zorder=4)
        sub.scatter(x[-1:], y[-1:], s=18, facecolor=PAPER, edgecolor="#6B3E75", zorder=4, linewidths=1)
        sub.set_xlim(-0.02, 1.02)
        sub.set_ylim(-0.02, 1.02)
        sub.set_aspect("equal")
        sub.axis("off")
        sd = np.sqrt(float(r["variance"]))
        sub.text(0.5, -0.06, f"s.d. {sd:.3f} rad", transform=sub.transAxes, fontsize=9, color=MUTED, ha="center", va="top")
        if row == 0:
            sub.text(0.5, 1.10, f"ring {int(r['ring']) + 1}\nradius {radii[int(r['ring'])]:.2f}", transform=sub.transAxes, fontsize=10, color=INK, ha="center", va="bottom")
    fig.text(0.555, 0.835, "One arm for each of three places on every ring: base (dot), tip (ring), joints", fontsize=10.5, color=MUTED, ha="left", va="bottom")
    source(fig, provenance().split("; research")[0] + "; link length 1/24 in these coordinates; the faint map is the run's elites (search-24)")
    return save(fig, "search-25-arms")


# ------------------------------------------------------------------------------------------------------------------
# search-26: coverage and QD-score over evaluations.

def figure_curves():
    C = defaultdict(list)
    for r in rows("curves.csv"):
        C[r["algorithm"]].append(r)
    fig, axes = figure(1, 3)
    fig.subplots_adjust(left=0.06, right=0.98, top=0.80, bottom=0.12, wspace=0.26)
    title(fig, "How fast each search fills the map, and how good its arms become",
          "On the common grid of 50 × 50 cells: median over 15 seeds (line) and the interquartile range (band); MAP-Elites samples at random for its first 1 000")
    for alg in ORDER:
        rs = C[alg]
        e = np.array([float(r["evaluations"]) for r in rs])
        for ax, (lo, mid, hi) in zip(axes[:2], [("coverage_q1", "coverage_median", "coverage_q3"), ("qd_q1", "qd_median", "qd_q3")]):
            ax.fill_between(e, [float(r[lo]) for r in rs], [float(r[hi]) for r in rs], color=COLOUR[alg], alpha=0.14, linewidth=0)
            ax.plot(e, [float(r[mid]) for r in rs], color=COLOUR[alg], lw=1.7, label=alg)
        axes[2].plot(e, np.sqrt([float(r["mean_variance_median"]) for r in rs]), color=COLOUR[alg], lw=1.7)
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("evaluations")
        ax.set_xlim(100, 1e5)
    axes[0].set_title("Coverage of the reachable cells", loc="left")
    axes[0].set_ylim(0, 1.0)
    axes[1].set_title("QD-score, normalised", loc="left")
    axes[1].set_ylim(0, 1.0)
    axes[2].set_title("Typical spread of an elite's angles", loc="left")
    axes[2].set_ylabel("root of the elites' mean variance (rad)")
    axes[2].set_yscale("log")
    axes[2].set_yticks([0.2, 0.3, 0.5, 1.0, 1.5])
    axes[2].yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
    axes[2].yaxis.set_minor_formatter(plt.NullFormatter())
    axes[0].legend(loc="upper left", fontsize=9.5)
    source(fig, provenance().split("; research")[0] + "; QD-score = Σ (π² − var) over the filled cells, divided by π² × the 2 040 reachable cells")
    return save(fig, "search-26-coverage-and-qd")


if __name__ == "__main__":
    for f in (figure_maps, figure_arms, figure_curves):
        print(f())
