"""Figures of the walk past Melbourne's sensors, drawn from results/ in the house style of the search case studies.

    .venv/Scripts/python.exe research/experiments/2026-11-05-a-walk-past-melbournes-sensors/figures.py [results/]

The numbers come from results/ alone. The streets drawn beneath the tours are not results: they are read from the
data library's city model of Melbourne (data/cache/products/city-model-melbourne.json, built by
data/products/street_records.py from OpenStreetMap, ODbL 1.0), whose local ENU frame the program's sensors.csv
shares; without that file the tours are drawn on plain paper.
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys

import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F401,F403
from search_style import ALGORITHM, HAIRLINE, INK, MUTED, PANEL, PAPER, ROOT, figure, save, source, title  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE / "results"
CITY = ROOT / "data" / "cache" / "products" / "city-model-melbourne.json"

NAMES = {
    "random": "Random tours",
    "Plackett-Luce": "Plackett–Luce EDA",
    "GA-CX": "GA, cycle crossover",
    "GA-PMX": "GA, partially mapped",
    "GA-OX": "GA, order crossover",
    "MMAS": "MAX–MIN ant system",
    "SA": "Simulated annealing",
    "ILS": "Iterated local search",
}
ORDER = ["random", "Plackett-Luce", "GA-CX", "GA-PMX", "GA-OX", "MMAS", "SA", "ILS"]
COLOUR = {
    "random": ALGORITHM["random"],
    "Plackett-Luce": ALGORITHM["DE"],
    "GA-CX": "#A9C29F",
    "GA-PMX": ALGORITHM["GA"],
    "GA-OX": "#4F7A4B",
    "MMAS": ALGORITHM["ACO"],
    "SA": ALGORITHM["SA"],
    "ILS": ALGORITHM["ILS"],
}
TOUR = "#8E3B78"
# Darker inks for thin tours on a faint map.
DARK = {"random": "#8A7D86", "GA-CX": "#7E9C74", "GA-PMX": "#5F8559", "GA-OX": "#3F6A3B", "MMAS": "#B08A1C"}
CREDIT_SENSORS = "Sensors: City of Melbourne, Pedestrian Counting System (CC BY 4.0), retrieved 30 Oct 2026"
CREDIT_OSM = "streets © OpenStreetMap contributors (ODbL 1.0)"


def read_sensors():
    with open(RESULTS / "sensors.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    return np.array([[float(r["east_m"]), float(r["north_m"])] for r in rows])


def read_city():
    if not CITY.exists():
        return None
    d = json.loads(CITY.read_text(encoding="utf-8"))
    major = {"motorway", "trunk", "primary", "secondary", "tertiary", "motorway_link", "trunk_link", "primary_link", "secondary_link", "tertiary_link"}
    minor = {"residential", "unclassified", "living_street", "service", "pedestrian", "busway"}
    roads = {"major": [], "minor": [], "path": []}
    for r in d["roads"]:
        k = r["kind"]
        group = "major" if k in major else "minor" if k in minor else "path" if k in {"footway", "path", "cycleway", "steps"} else None
        if group:
            roads[group].append(np.array(r["line"]))
    spaces = {"water": [], "park": []}
    for s in d["spaces"]:
        if s["kind"] == "water":
            spaces["water"].append(np.array(s["ring"]))
        elif s["kind"] in ("park", "green"):
            spaces["park"].append(np.array(s["ring"]))
    return roads, spaces


def draw_city(ax, city, strength=1.0, faint=False):
    if city is None:
        return
    roads, spaces = city
    ax.add_collection(PolyCollection(spaces["water"], facecolors="#E2E8EC" if faint else "#DCE4EA", edgecolors="none", zorder=0))
    ax.add_collection(PolyCollection(spaces["park"], facecolors="#EEEFE5" if faint else "#E9EBDD", edgecolors="none", zorder=0))
    if not faint:
        ax.add_collection(LineCollection(roads["path"], colors="#E6DCDF", linewidths=0.25 * strength, zorder=1))
    ax.add_collection(LineCollection(roads["minor"], colors="#E9E0E3" if faint else "#DDD0D6", linewidths=0.45 * strength, zorder=1))
    ax.add_collection(LineCollection(roads["major"], colors="#DED2D8" if faint else "#CDBDC6", linewidths=0.9 * strength, zorder=1))


def draw_tour(ax, xy, tour, colour=TOUR, width=1.8, dots=18):
    t = np.array(tour + tour[:1])
    ax.plot(xy[t, 0], xy[t, 1], color=colour, lw=width, solid_joinstyle="round", zorder=3)
    ax.scatter(xy[:, 0], xy[:, 1], s=dots, color=PAPER, edgecolors=INK, linewidths=0.7, zorder=4)


def frame_map(ax, xy, pad=250):
    lo, hi = xy.min(axis=0) - pad, xy.max(axis=0) + pad
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_facecolor(PAPER)


def scale_bar(ax, xy, length=1000):
    lo = xy.min(axis=0)
    x0, y0 = lo[0] + 50, lo[1] - 120
    ax.plot([x0, x0 + length], [y0, y0], color=INK, lw=1.2, solid_capstyle="butt", zorder=5)
    for x in (x0, x0 + length):
        ax.plot([x, x], [y0 - 25, y0 + 25], color=INK, lw=1.0, zorder=5)
    ax.text(x0 + length / 2, y0 + 45, f"{length / 1000:g} km", ha="center", va="bottom", fontsize=10, color=INK)
    # North, beside the scale (the local frame's y axis points to geodetic north at the origin).
    xn = x0 + length + 250
    ax.annotate("", xy=(xn, y0 + 160), xytext=(xn, y0 - 40), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0), zorder=5)
    ax.text(xn + 70, y0 + 60, "N", ha="left", va="center", fontsize=11, color=INK)


def km(x):
    return f"{x / 1000:,.2f} km".replace(",", " ")


def figure_tour(xy, best, city):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "A walk past Melbourne's pedestrian sensors",
          f"The shortest closed walk known through all {best['sensors']} counters, {km(best['length_m'])} along geodesics on WGS 84")
    ax = fig.add_axes([0.03, 0.07, 0.60, 0.81])
    draw_city(ax, city)
    draw_tour(ax, xy, best["tour"], width=2.0, dots=22)
    frame_map(ax, xy)
    scale_bar(ax, xy)

    # The evidence for the best known tour: its restarts.
    side = fig.add_axes([0.70, 0.50, 0.26, 0.30])
    lengths = np.array(best["restart_lengths_m"]) / 1000
    order = np.sort(lengths)
    excess = 100 * (order / order[0] - 1)
    side.scatter(np.arange(1, len(order) + 1), excess, s=26, color=TOUR, edgecolors=PAPER, linewidths=0.6, zorder=3)
    side.set_xlabel("restart, by its result")
    side.set_ylabel("excess over the best (%)")
    side.set_title(f"{best['restarts_reaching_best']} of {best['restarts']} restarts end at {best['length_m'] / 1000:.4f} km", loc="left", fontsize=13, color=INK)
    side.set_xlim(0, len(order) + 1)
    side.set_ylim(-0.05, max(0.5, excess.max() * 1.15))

    notes = [
        "Iterated local search on the whole problem: 2-opt with",
        "moves evaluated from four distances, double-bridge kicks,",
        f"{best['kicks']:,} kicks from each of {best['restarts']} random starts.".replace(",", " "),
        "The tour is the best known, not proved optimal.",
        "",
        "Distances are Karney's geodesics on the WGS 84 ellipsoid",
        "(the geodesy domain). In the city model's local frame the",
        f"same walk measures {km(best['length_in_local_frame_m'])}: chords differ from",
        f"geodesics by at most {best['geodesic_vs_chord_max_abs_m'] * 1000:.1f} mm over every pair.",
    ]
    fig.text(0.70, 0.40, "\n".join(notes), fontsize=11.5, color=INK, va="top", linespacing=1.55)
    source(fig, f"{CREDIT_SENSORS}; {CREDIT_OSM}, from the data library's city model. Search: heuresis.")
    return save(fig, "search-30-melbourne-best-tour")


def figure_multiples(xy, tours, final, city):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    budget = tours["budget"]
    words = "a million" if budget == 1e6 else f"{budget:,.0f}".replace(",", " ")
    title(fig, f"Each search's best walk after {words} evaluations",
          f"The best of each algorithm's runs at {budget:,.0f} evaluations; its length and its excess over the best known walk".replace(",", " "))
    best = tours["best_known_m"]
    medians = {}
    for a in ORDER:
        v = [r["best_m"] for r in final if r["algorithm"] == a]
        if v:
            medians[a] = (np.median(v), len(v))
    for i, a in enumerate(ORDER):
        if a not in tours["algorithms"]:
            continue
        r, c = divmod(i, 4)
        ax = fig.add_axes([0.03 + c * 0.24, 0.47 - r * 0.40, 0.22, 0.36])
        draw_city(ax, city, strength=0.7, faint=True)
        t = tours["algorithms"][a]
        draw_tour(ax, xy, t["tour"], colour=DARK.get(a, COLOUR[a]), width=1.2, dots=5)
        frame_map(ax, xy, pad=150)
        L = t["length_m"]
        ax.set_title(NAMES[a], loc="left", fontsize=14, color=INK, pad=4)
        m, n = medians[a]
        ax.text(0.0, -0.035, f"{km(L)}   +{100 * (L / best - 1):.1f} %   (median of {n}: +{100 * (m / best - 1):.1f} %)", transform=ax.transAxes, fontsize=10.5, color=MUTED, va="top")
    seeds = sorted({(r["algorithm"], r["seed"]) for r in final})
    n_seeds = max(sum(1 for a, _ in seeds if a == b) for b in ORDER)
    n_ant = sum(1 for a, _ in seeds if a == "MMAS")
    source(fig, f"Each panel: the shortest walk among the algorithm's runs ({n_seeds} seeds; the ant system {n_ant}). {CREDIT_SENSORS}; {CREDIT_OSM}.")
    return save(fig, "search-31-melbourne-tours-by-algorithm")


def figure_convergence(convergence, best_known, timing):
    fig, ax = figure()
    fig.subplots_adjust(left=0.08, right=0.80, top=0.84, bottom=0.11)
    title(fig, "How fast each search closes on the best walk",
          "Excess length over the best known walk against evaluations: the median of the runs, and the band between their quartiles")
    labels = []
    for a in ORDER:
        rows = [r for r in convergence if r["algorithm"] == a]
        if not rows:
            continue
        evals = sorted({r["evaluations"] for r in rows})
        seeds = sorted({r["seed"] for r in rows})
        table = {(r["seed"], r["evaluations"]): r["best_m"] for r in rows}
        e_ok, med, q1, q3 = [], [], [], []
        for e in evals:
            v = [table[(s, e)] for s in seeds if (s, e) in table]
            if len(v) < len(seeds):
                continue
            x = 100 * (np.array(v) / best_known - 1)
            e_ok.append(e)
            med.append(np.median(x))
            q1.append(np.percentile(x, 25))
            q3.append(np.percentile(x, 75))
        col = COLOUR[a]
        ax.fill_between(e_ok, q1, q3, color=col, alpha=0.18, lw=0)
        ax.plot(e_ok, med, color=col, lw=2.0)
        labels.append([np.log10(med[-1]), NAMES[a], DARK.get(a, col), e_ok[-1], med[-1]])
    # Labels at the curves' ends, kept at least 0.07 decades apart.
    labels.sort()
    for i in range(1, len(labels)):
        labels[i][0] = max(labels[i][0], labels[i - 1][0] + 0.085)
    for y, name, col, e, m in labels:
        ax.text(e * 1.12, 10 ** y, f"{name}  {m:.1f} %", color=col, fontsize=12, va="center", ha="left")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("evaluations (tour lengths computed)")
    ax.set_ylabel("excess over the best known walk (%)")
    ax.set_xlim(100, 1e6)
    ax.set_yticks([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000])
    ax.set_yticklabels(["1", "2", "5", "10", "20", "50", "100", "200", "500", "1000"])
    ax.set_ylim(0.8, 1500)
    seeds = timing.get("seeds", 10)
    source(fig, f"{seeds} seeds each (the ant system {timing.get('ant_seeds', seeds)}), heuresis's drivers and anytime curves; best known {best_known / 1000:.4f} km. {CREDIT_SENSORS}.")
    return save(fig, "search-32-melbourne-convergence")


def main():
    xy = read_sensors()
    best = json.loads((RESULTS / "best_known.json").read_text())
    tours = json.loads((RESULTS / "tours.json").read_text())
    timing = json.loads((RESULTS / "timing.json").read_text())
    with open(RESULTS / "final.csv", newline="") as f:
        final = [{**r, "best_m": float(r["best_m"])} for r in csv.DictReader(f)]
    with open(RESULTS / "convergence.csv", newline="") as f:
        convergence = [{"algorithm": r["algorithm"], "seed": int(r["seed"]), "evaluations": float(r["evaluations"]), "best_m": float(r["best_m"])} for r in csv.DictReader(f)]
    city = read_city()
    for p in (figure_tour(xy, best, city), figure_multiples(xy, tours, final, city), figure_convergence(convergence, best["length_m"], timing)):
        print(p)


if __name__ == "__main__":
    main()
