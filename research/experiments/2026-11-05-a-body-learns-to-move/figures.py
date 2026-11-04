"""Figures of the case study `2026-11-05-a-body-learns-to-move`, drawn from results/ only.

    search-14-rabbit-convergence  five searches tuning the rabbit's gains: the best feasible cost of transport over evaluations
    search-15-rabbit-budget       where each search spends its trials: falls, failures, infeasible and feasible gaits; the final costs
    search-16-rabbit-gains        the gains each search settled on, against the build's defaults
    search-17-rabbit-front        speed against the power drawn: the front of NSGA-II and SMS-EMOA
    search-18-rabbit-strobes      four gaits of the front, seen from the side
    search-19-rabbit-before-after the default gains and the tuned ones at 0.6 m/s
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F403

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
ALGOS = ["random", "(1+1)-ES", "CMA-ES", "DE", "PSO"]
MASS, MOON = 1.69, 1.62  # the rabbit's mass (robotics/legged/rabbit.hpp) and the Moon's gravity at its site (m/s²)
SOURCE = "The robotics domain's rabbit (robotics/legged, ADR 0058) on the Moon, semi-implicit Euler at 0.5 ms, its controller at 500 Hz; heuresis (6b3cc88)."
REASONS = ["fell", "tipped over", "stopped hopping", "refused", "diverged"]
REASON_COLOUR = {"fell": "#6B2333", "tipped over": "#C2665A", "stopped hopping": "#E7B9A1", "refused": "#6C8EBF", "diverged": "#3B2F3A"}


def load(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def grid_of(runs):
    return np.array([c[0] for c in runs[0]["checkpoints"]], dtype=float)


def series(run, column):
    return np.array([np.inf if c[column] is None else c[column] for c in run["checkpoints"]], dtype=float)


# ---- search-14: convergence ------------------------------------------------------------------------------------------

def convergence():
    d = load("tune.json")
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(16, 10), dpi=100, sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    fig.subplots_adjust(left=0.08, right=0.8, top=0.83, bottom=0.09, hspace=0.12)
    t = d["task"]
    title(fig, "Five searches teach the rabbit to hop at 0.6 m/s",
          f"Ten gains of its controller; the cost of transport over the last {t['window']:.0f} s of {t['hop']:.0f}, the speed within 3 cm/s, the apex within 2 cm of {t['height']*100:.0f} cm, the pitch within 20°; {d['seeds']} seeds each")
    ends = []
    for a in ALGOS:
        rs = [r for r in d["runs"] if r["algorithm"] == a]
        n = min(len(r["checkpoints"]) for r in rs)
        x = grid_of(rs)[:n]
        Y = np.array([series(r, 1)[:n] for r in rs])
        found = np.isfinite(Y)
        med = np.array([np.median(c) if np.isfinite(np.median(c)) else np.nan for c in Y.T])
        lo = np.array([np.quantile(c, 0.25) for c in Y.T])
        hi = np.array([np.quantile(c, 0.75) for c in Y.T])
        colour = ALGORITHM[a]
        ax.plot(x, med, color=colour, lw=2.4, label=a)
        ax.fill_between(x, np.where(np.isfinite(lo), lo, np.nan), np.where(np.isfinite(hi), hi, np.nan), color=colour, alpha=0.12, lw=0)
        bx.plot(x, found.mean(axis=0), color=colour, lw=2.2)
        ends.append((med[-1], a, colour))
    ax.set_ylabel("best feasible cost of transport (median, quartiles)")
    ax.set_ylim(1.3, 2.0)
    ends.sort()
    placed = -np.inf
    for yv, a, colour in ends:
        placed = max(yv, placed + 0.03)
        ax.text(x[-1] * 1.015, placed, f"{a}  {yv:.3f}", color=colour, fontsize=12, va="center")
    ax.text(0.99, 0.03, f"the build's default gains reach {d['default']['speed']:.2f} m/s at a cost of {d['default']['cot']:.2f}, pitching to {d['default']['pitch_deg']:.0f} degrees: infeasible",
            transform=ax.transAxes, fontsize=12, color=MUTED, va="bottom", ha="right")
    bx.set_ylabel("seeds with a\nfeasible gait")
    bx.set_ylim(-0.03, 1.05)
    bx.set_xlabel("trials (evaluations)")
    source(fig, "Cost of transport: electrical energy over m g d. CMA-ES and the (1+1)-ES start at a uniform point of the box, step 0.25 of its width; DE (rand/1/bin) and PSO (inertia, ring) with 20 members. heuresis (6b3cc88).")
    save(fig, "search-14-rabbit-convergence")


# ---- search-15: where the trials go ---------------------------------------------------------------------------------

def budget():
    d = load("tune.json")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], left=0.07, right=0.97, top=0.83, bottom=0.17, wspace=0.22)
    ax, bx = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    title(fig, "Where the trials go", "Each search's trials by what became of them, pooled over its seeds; and each seed's best feasible cost of transport at the end")
    y = np.arange(len(ALGOS))[::-1]
    for yi, a in zip(y, ALGOS):
        rs = [r for r in d["runs"] if r["algorithm"] == a]
        total = sum(r["evaluations"] for r in rs)
        left = 0.0
        parts = [(k, sum(r["reasons"].get(k, 0) for r in rs)) for k in REASONS]
        failures = sum(v for _, v in parts)
        feasible = sum(r["feasible"] for r in rs)
        parts += [("infeasible", total - failures - feasible), ("feasible", feasible)]
        for k, v in parts:
            w = v / total
            colour = REASON_COLOUR.get(k, {"infeasible": "#D9CCD3", "feasible": ALGORITHM[a]}.get(k))
            ax.barh(yi, w, left=left, color=colour, height=0.62, edgecolor=PAPER, lw=0.8)
            if w > 0.06:
                ax.text(left + w / 2, yi, f"{w*100:.0f}%", ha="center", va="center", fontsize=11, color=PAPER if k in ("fell", "diverged", "feasible") else INK)
            left += w
    ax.set_yticks(y, ALGOS)
    ax.set_xlim(0, 1)
    ax.set_xlabel("share of the trials")
    ax.grid(axis="y", visible=False)
    present = [k for k in REASONS if any(r["reasons"].get(k, 0) for r in d["runs"])]
    handles = [plt.Rectangle((0, 0), 1, 1, color=REASON_COLOUR[k]) for k in present] + [plt.Rectangle((0, 0), 1, 1, color="#D9CCD3"), plt.Rectangle((0, 0), 1, 1, color=ALGORITHM["CMA-ES"])]
    ax.legend(handles, present + ["hopped, outside the bounds asked", "feasible (in the search's colour)"], loc="upper left", bbox_to_anchor=(0.0, -0.1), ncol=4, fontsize=11)
    ax.set_title("what became of the trials", loc="left", fontsize=14, color=INK)
    for yi, a in zip(y, ALGOS):
        vals = [r["best_value"] for r in d["runs"] if r["algorithm"] == a and not r["best_failed"] and (r["best_violation"] or 0) == 0]
        jitter = np.linspace(-0.12, 0.12, len(vals)) if len(vals) > 1 else [0]
        bx.scatter(vals, yi + np.array(jitter), s=55, color=ALGORITHM[a], zorder=3)
        if vals:
            bx.plot([np.median(vals)] * 2, [yi - 0.25, yi + 0.25], color=INK, lw=1.4)
        missing = sum(1 for r in d["runs"] if r["algorithm"] == a) - len(vals)
        if missing:
            bx.text(1.555, yi - 0.3, f"{missing} seed: no feasible gait in {d['budget']:.0f} trials", fontsize=11, color=ALGORITHM[a], va="center")
    bx.set_yticks(y, ALGOS)
    bx.set_xlabel("best feasible cost of transport at the end")
    bx.set_title("each seed's result (bar: the median)", loc="left", fontsize=14, color=INK)
    bx.grid(axis="y", visible=False)
    source(fig, f"{d['budget']:.0f} trials a run. Fell: the torso touched the ground; tipped over: the pitch beyond 60 degrees. No hop was refused or stopped. heuresis (6b3cc88).")
    save(fig, "search-15-rabbit-budget")


# ---- search-16: the gains ------------------------------------------------------------------------------------------

def gains():
    d = load("tune.json")
    params = d["parameters"]
    fig, ax = plt.subplots(figsize=(16, 10), dpi=100)
    fig.subplots_adjust(left=0.17, right=0.97, top=0.83, bottom=0.15)
    title(fig, "The gains the searches settled on", "Each run's best feasible gains on the unit box of the search (0: the lower end of the range, 1: the upper); the build's defaults ringed")
    names = [p["name"].replace("_", " ") for p in params]
    y = np.arange(len(params))[::-1]
    offsets = np.linspace(-0.3, 0.3, len(ALGOS))
    for off, a in zip(offsets, ALGOS):
        rs = [r for r in d["runs"] if r["algorithm"] == a and not r["best_failed"] and (r["best_violation"] or 0) == 0]
        for r in rs:
            ax.scatter(r["best_x"], y + off, s=26, color=ALGORITHM[a], alpha=0.85, zorder=3)
        ax.scatter([], [], s=40, color=ALGORITHM[a], label=a)
    ax.scatter(d["default_x"], y, s=260, facecolor="none", edgecolor=INK, lw=1.5, zorder=4, label="the build's default")
    best = d["runs"][d["best_run"]]["best_x"]
    ax.scatter(best, y, s=140, marker="D", facecolor="none", edgecolor=ALGORITHM["BO"], lw=1.8, zorder=4, label="the best of all runs")
    labels = []
    for p in params:
        u = f" ({p['unit']})" if p["unit"] else ""
        lo, hi = p["lower"], p["upper"]
        if p["unit"] == "rad":
            lo, hi, u = np.degrees(lo), np.degrees(hi), " (°)"
        labels.append(f"{p['name'].replace('_', ' ')}{u}\n{lo:g} to {hi:g}")
    ax.set_yticks(y, labels, fontsize=11)
    for yy in y[:-1]:
        ax.axhline(yy - 0.5, color=HAIRLINE, lw=0.6)
    ax.set_xlim(-0.03, 1.03)
    ax.set_xlabel("position in the range searched")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.075), ncol=7, fontsize=11)
    source(fig, "Ten fields of robotics::legged::Gains (legged/controller.hpp): the virtual leg, the ankle in a hop, attitude, Raibert's placement and its integral, the apex's servo, the thrust's cap. heuresis (6b3cc88).")
    save(fig, "search-16-rabbit-gains")


# ---- search-17: the front ------------------------------------------------------------------------------------------

def front():
    d = load("front.json")
    explore = load("explore.json")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.7, 1], left=0.07, right=0.97, top=0.83, bottom=0.11, wspace=0.22)
    ax, bx = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    title(fig, "Faster costs power", f"Speed against the mean electrical power, both searched with the gains and the speed asked: the non-dominated gaits of {d['seeds']} runs each of NSGA-II and SMS-EMOA")
    v = np.linspace(0.02, 1.3, 50)
    for cot in (1.0, 1.5, 2, 3, 5):
        ax.plot(v, cot * MASS * MOON * v, color=HAIRLINE, lw=0.9, zorder=1)
        vv = min(1.2, 5.75 / (cot * MASS * MOON))
        ax.text(vv, cot * MASS * MOON * vv, f"CoT {cot:g} ", color=MUTED, fontsize=11, va="top", ha="right")
    for r in d["runs"]:
        f = np.array(sorted(r["front"]))
        if len(f):
            ax.plot(f[:, 0], f[:, 1], color=ALGORITHM[r["algorithm"]], lw=0.9, alpha=0.45, drawstyle="steps-post")
    u = d["union_front"]
    for a in ("NSGA-II", "SMS-EMOA"):
        pts = [p for p in u if p["from"] == a]
        ax.scatter([p["speed"] for p in pts], [p["power"] for p in pts], s=48, color=ALGORITHM[a], zorder=3, label=f"{a} (on the union's front)")
    dflt = [e for e in explore if abs(e["asked_height"] - 0.12) < 1e-9 and not e["trial"]["failed"] and e["asked_speed"] > 0]
    ax.scatter([e["trial"]["speed"] for e in dflt], [e["trial"]["power"] for e in dflt], s=60, facecolor=PAPER, edgecolor=INK, lw=1.3, zorder=3, label="the default gains, asked 0.1 to 0.6 m/s")
    ax.set_xlim(0, 1.25)
    ax.set_ylim(0, 6)
    ax.set_xlabel("speed reached (m/s)")
    ax.set_ylabel("mean electrical power (W)")
    ax.legend(loc="lower right")
    ax.set_title("the front, and lines of equal cost of transport", loc="left", fontsize=14, color=INK)
    for r in d["runs"]:
        h = np.array(r["hypervolume"])
        bx.plot(h[:, 0], h[:, 1], color=ALGORITHM[r["algorithm"]], lw=1.8, alpha=0.9)
    for a in ("NSGA-II", "SMS-EMOA"):
        bx.plot([], [], color=ALGORITHM[a], lw=2, label=a)
    bx.legend(loc="lower right")
    bx.set_xlabel("trials")
    bx.set_ylabel(f"hypervolume (reference: 0 m/s, {d['reference'][1]:g} W)")
    bx.set_title("each run's hypervolume", loc="left", fontsize=14, color=INK)
    source(fig, "Constraints: the apex within 2 cm of 12 cm, the pitch within 20 degrees; falls fail. Population 40, SBX and polynomial mutation. The rabbit on the Moon, semi-implicit Euler; heuresis (6b3cc88).")
    save(fig, "search-17-rabbit-front")


# ---- strobes ---------------------------------------------------------------------------------------------------------

def draw_pose(ax, p, x_shift, colour, alpha, lw=1.6):
    P = np.array(p)
    P[:, 0] -= x_shift
    ax.add_patch(Polygon(P[:4], closed=True, facecolor=colour, edgecolor=colour, alpha=alpha * 0.2, lw=lw * 0.6))
    segs = [[P[4], P[5]], [P[5], P[6]], [P[6], P[7]], [P[8], P[7]], [P[9], P[10]], [P[10], P[11]]]
    ax.add_collection(LineCollection(segs, colors=[colour], linewidths=lw, alpha=alpha, capstyle="round"))


def pitch_of(p):
    c = np.array(p[:4])
    d = c[1] - c[0]
    return np.degrees(np.arctan2(d[1], d[0]))


def strobe(ax, poses, colour, t_from, t_to, every=4, label=None):
    sel = [q for q in poses if t_from <= q["t"] <= t_to]
    if not sel:
        return
    x0 = sel[0]["p"][12][0]
    for k, q in enumerate(sel):
        if k % every:
            continue
        a = 0.25 + 0.75 * k / max(1, len(sel) - 1)
        draw_pose(ax, q["p"], x0, colour, a)
    com = np.array([q["p"][12] for q in sel])
    ax.plot(com[:, 0] - x0, com[:, 1], color=colour, lw=0.8, alpha=0.6, ls=(0, (2, 2)))
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_aspect("equal")
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_yticks([])
    if label:
        ax.text(0.0, 1.02, label, transform=ax.transAxes, fontsize=12, color=INK, va="bottom")


def strobes():
    ms = load("front_motions.json")
    fig, axes = plt.subplots(len(ms), 1, figsize=(16, 10), dpi=100)
    fig.subplots_adjust(left=0.05, right=0.98, top=0.86, bottom=0.1, hspace=0.6)
    title(fig, "Four gaits of the front", "Each row 2.5 s of the window, the body every 0.2 s, later poses darker; the dashed line the centre of mass; the ground at 0")
    cmap = [ALGORITHM["SA"], ALGORITHM["DE"], ALGORITHM["PSO"], ALGORITHM["BO"]]
    for ax, m, colour in zip(np.atleast_1d(axes), ms, cmap):
        t = m["trial"]
        lab = f"asked {m['asked']:.2f} m/s: reached {t['speed']:.2f} m/s at {t['power']:.2f} W, cost of transport {t['cot']:.2f}, apex {t['apex']*100:.1f} cm, pitch within {t['pitch_deg']:.0f}°"
        end = m["poses"][-1]["t"]
        strobe(ax, m["poses"], colour, end - 2.5, end, every=6, label=lab)
        ax.set_xlim(-0.35, 1.0 * 2.5 + 0.4)
        ax.set_ylim(-0.02, 0.42)
        ax.set_xticks(np.arange(0, 3.0, 0.5))
        ax.tick_params(axis="x", labelsize=10)
    np.atleast_1d(axes)[-1].set_xlabel("distance from the row's first pose (m)")
    source(fig, "Side view of the planar model: the torso, the hind leg (thigh, shank, foot) and the foreleg. The rabbit on the Moon, semi-implicit Euler; heuresis (6b3cc88).")
    save(fig, "search-18-rabbit-strobes")


def before_after():
    m = load("tune_motions.json")
    d = load("tune.json")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1.1], left=0.06, right=0.98, top=0.85, bottom=0.1, hspace=0.55)
    title(fig, "Before and after: asked to hop at 0.6 m/s", "The default gains and the best tuned ones on the same body and ground: strobes from 1.6 s to 4.1 s, as the default pitches over gathering speed")
    for i, (key, colour, name) in enumerate((("default", MUTED, "default gains"), ("best", ALGORITHM["CMA-ES"], "tuned gains"))):
        ax = fig.add_subplot(gs[i, 0])
        t = d[key]
        lab = f"{name}: {t['speed']:.2f} m/s, cost of transport {t['cot']:.2f}, {t['power']:.2f} W, apex {t['apex']*100:.1f} cm, pitch within {t['pitch_deg']:.0f}°"
        end = m[key][-1]["t"]
        strobe(ax, m[key], colour, 1.6, 4.1, every=6, label=lab)
        ax.set_xlim(-0.35, 1.9)
        ax.set_ylim(-0.02, 0.45)
    cx = fig.add_subplot(gs[2, 0])
    for key, colour, name in (("default", MUTED, "default gains"), ("best", ALGORITHM["CMA-ES"], "tuned gains")):
        ts = np.array([q["t"] for q in m[key]])
        ps = np.array([pitch_of(q["p"]) for q in m[key]])
        cx.plot(ts, ps, color=colour, lw=1.8, label=name)
    cx.axhspan(-20, 20, color=HAIRLINE, alpha=0.35, lw=0)
    cx.text(0.1, 21, "the bound asked, ±20°", color=MUTED, fontsize=11)
    cx.set_xlabel("time (s); the hop asked at 0.5 s", labelpad=1)
    cx.set_ylabel("pitch (°)")
    cx.legend(loc="lower left", ncol=2)
    source(fig, "The tuned gains: the best feasible point of all runs of search-14 (CMA-ES, seed 2005). The rabbit on the Moon, semi-implicit Euler; heuresis (6b3cc88).")
    save(fig, "search-19-rabbit-before-after")


if __name__ == "__main__":
    convergence()
    budget()
    gains()
    front()
    strobes()
    before_after()
