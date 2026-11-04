"""Figures of the case study `2026-11-05-keplers-third-law`, drawn from results/ only.

    search-10-kepler-law      the law on log-log data, and its residuals body by body
    search-11-kepler-front    the front of formulas, error against size, typed and untyped
    search-12-kepler-success  how often, and how soon, each search meets the law
    search-13-kepler-units    a law or a coincidence of the units: each formula's error before and after a change of units
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F403

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
SOURCE = "JPL DE440s (NAIF, retrieved 2026-10-21) read by the astronomy domain's SPK reader; heuresis programs (6b3cc88); 30 seeds of 30 000 evaluations each."

TYPED, UNTYPED = ALGORITHM["GP"], ALGORITHM["GP untyped"]
PLANET, MOON = "#8E3B78", "#4E8F8B"


def orbits():
    with open(RESULTS / "orbits.csv", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def runs():
    return json.loads((RESULTS / "gp_runs.json").read_text(encoding="utf-8"))


def fronts():
    return json.loads((RESULTS / "gp_fronts.json").read_text(encoding="utf-8"))


AU3_YR2 = r"au$^3$ yr$^{-2}$"

# ---- formulas: the prefix terms of results/ written for the page --------------------------------------------------

ARITY = {"+": 2, "-": 2, "*": 2, "/": 2, "sq": 1, "cube": 1, "sqrt": 1}
VARIABLES = {"a", "GM", "e"}


def parse(text: str):
    """A prefix term, as heuresis's to_string writes it, as nested tuples (head, args...)."""
    pos = 0

    def term():
        nonlocal pos
        start = pos
        while pos < len(text) and text[pos] not in "(),":
            pos += 1
        head = text[start:pos].strip()
        if pos < len(text) and text[pos] == "(":
            pos += 1
            args = []
            while True:
                args.append(term())
                while text[pos] == " ":
                    pos += 1
                if text[pos] == ",":
                    pos += 1
                    while text[pos] == " ":
                        pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
            return (head, *args)
        return (head,)

    return term()


def has_variable(t) -> bool:
    return t[0] in VARIABLES or any(has_variable(a) for a in t[1:])


def value(t) -> float:
    h = t[0]
    if len(t) == 1:
        return np.pi if h == "pi" else float(h)
    x = [value(a) for a in t[1:]]
    return {"+": lambda: x[0] + x[1], "-": lambda: x[0] - x[1], "*": lambda: x[0] * x[1], "/": lambda: x[0] / x[1],
            "sq": lambda: x[0] ** 2, "cube": lambda: x[0] ** 3, "sqrt": lambda: np.sqrt(x[0])}[h]()


def number(v: float) -> str:
    for k, name in ((2 * np.pi, r"2\pi"), (np.pi, r"\pi"), (4 * np.pi**2, r"4\pi^2")):
        if abs(v - k) < 1e-12:
            return name
    if abs(v - round(v)) < 1e-12 and abs(v) < 1e6:
        return f"{v:.0f}"
    return f"{v:.4g}"


def tex(t, top=True):
    """(text, precedence): 0 a sum, 1 a product or quotient, 2 an atom or a power. Constant subterms are folded."""
    h = t[0]
    if len(t) == 1:
        return ({"pi": r"\pi", "GM": "GM"}.get(h, h) if h in VARIABLES or h == "pi" else number(float(h)), 2)
    if not has_variable(t):
        v = value(t)
        return number(v), (0 if v < 0 else 2)
    if ARITY[h] == 1:
        x, p = tex(t[1], False)
        if h == "sqrt":
            return rf"\sqrt{{{x}}}", 2
        bare = p == 2 and "sqrt" not in x and "^" not in x
        return (x if bare else f"({x})") + ("^2" if h == "sq" else "^3"), 2
    (x, px), (y, py) = tex(t[1], False), tex(t[2], False)
    if h == "+":
        return f"{x} + {y}", 0
    if h == "-":
        return f"{x} - " + (f"({y})" if py == 0 else y), 0
    if h == "/":
        return (f"({x})" if px == 0 else x) + "/" + (f"({y})" if py < 2 else y), 1
    if not has_variable(t[2]) and has_variable(t[1]):  # a number first
        (x, px), (y, py) = (y, py), (x, px)
    return (f"({x})" if px == 0 else x) + r"\," + (f"({y})" if py == 0 else y), 1


def formula(prefix: str) -> str:
    return "$" + tex(parse(prefix))[0] + "$"


def sci(x: float, digits: int = 1) -> str:
    """x as a × 10^n in mathtext."""
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10**e
    if abs(round(m, digits) - 10) < 1e-9:
        m, e = 1.0, e + 1
    return rf"{m:.{digits}f}\times10^{{{e}}}"


# ---- search-10: the law ------------------------------------------------------------------------------------------

def law():
    rows = orbits()
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 1], left=0.075, right=0.97, top=0.82, bottom=0.1, hspace=0.42, wspace=0.18)
    ax1, ax2, ax3 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
    title(fig, "Kepler's third law in the ephemeris", "Ten orbits of JPL's DE440s over 1850–2150: the periods measured, the semi-major axes averaged, two centres on one line")
    a = np.array([float(r["a_mean_m"]) for r in rows])
    T = np.array([float(r["period_s"]) for r in rows])
    mu = np.array([float(r["mu_m3s2"]) for r in rows])
    names = [r["body"] for r in rows]
    is_moon = np.array([n == "Moon" for n in names])
    day = 86400.0
    year = 365.25 * day

    # T against a: the planets on a line of slope 3/2, the Moon far below it.
    xs = np.geomspace(a.min() / 3, a.max() * 2, 50)
    ax1.plot(xs / 1.495978707e11, 2 * np.pi * np.sqrt(xs**3 / 1.32712440041279419e20) / year, color=HAIRLINE, lw=1.4, zorder=1)
    ax1.scatter(a[~is_moon] / 1.495978707e11, T[~is_moon] / year, s=60, color=PLANET, zorder=3, label="planets about the Sun")
    ax1.scatter(a[is_moon] / 1.495978707e11, T[is_moon] / year, s=60, color=MOON, zorder=3, label="the Moon about the Earth")
    offsets = {"Mercury": (8, -4), "Venus": (-44, 4), "Earth": (8, -10), "Mars": (-34, 6), "Jupiter": (8, -6), "Saturn": (8, -6), "Uranus": (8, -6),
               "Neptune": (-58, 2), "Pluto": (8, -8), "Moon": (8, 2)}
    for n, x, y in zip(names, a / 1.495978707e11, T / year):
        ax1.annotate(n, (x, y), xytext=offsets[n], textcoords="offset points", fontsize=11, color=INK)
    ax1.set_xscale("log"), ax1.set_yscale("log")
    ax1.set_xlabel("mean semi-major axis $a$ (au)")
    ax1.set_ylabel("sidereal period $T$ (years)")
    ax1.set_title(r"$T$ against $a$: one line for one centre", loc="left", fontsize=14, color=INK)
    ax1.text(0.03, 0.93, r"the line: $T = 2\pi\sqrt{a^3/GM_\odot}$", transform=ax1.transAxes, fontsize=12, color=MUTED)
    ax1.legend(loc="lower right")

    # T against √(a³/GM): every orbit on T = 2π x.
    x = np.sqrt(a**3 / mu)
    xs = np.geomspace(x.min() / 3, x.max() * 2, 50)
    ax2.plot(xs, 2 * np.pi * xs, color=HAIRLINE, lw=1.4, zorder=1)
    ax2.scatter(x[~is_moon], T[~is_moon], s=60, color=PLANET, zorder=3)
    ax2.scatter(x[is_moon], T[is_moon], s=60, color=MOON, zorder=3)
    for n, xx, y in zip(names, x, T):
        if n in ("Moon", "Mercury", "Earth", "Jupiter", "Neptune"):
            ax2.annotate(n, (xx, y), xytext=(7, -4), textcoords="offset points", fontsize=11, color=INK)
    ax2.set_xscale("log"), ax2.set_yscale("log")
    ax2.set_xlabel(r"$\sqrt{a^3/GM}$ (s), with each orbit's $GM = G(M + m)$")
    ax2.set_ylabel("$T$ (s)")
    ax2.set_title(r"$T$ against $\sqrt{a^3/GM}$: one line for both centres", loc="left", fontsize=14, color=INK)
    ax2.text(0.03, 0.93, r"the line: $T = 2\pi\sqrt{a^3/GM}$", transform=ax2.transAxes, fontsize=12, color=MUTED)

    # Residuals, by their size on a log scale: Jacobi (filled) and heliocentric (hollow); the sign by the marker.
    res = np.array([float(r["residual"]) for r in rows])
    helio = np.array([float(r["helio_residual"]) if r["helio_residual"] else np.nan for r in rows])
    pos = np.arange(len(rows))
    colours = [MOON if m else PLANET for m in is_moon]
    for p, r, c in zip(pos, res, colours):
        ax3.vlines(p, 1e-7, abs(r), color=c, lw=2.0, alpha=0.5)
        ax3.scatter([p], [abs(r)], s=70, marker="^" if r > 0 else "v", color=c, zorder=3)
    for p, r in zip(pos, helio):
        if np.isfinite(r):
            ax3.scatter([p + 0.12], [abs(r)], s=60, marker="^" if r > 0 else "v", facecolor=PAPER, edgecolor=MUTED, lw=1.3, zorder=3)
    ax3.scatter([], [], s=60, marker="^", color=PLANET, label="Jacobi coordinates (the Moon: geocentric), the period longer than the law's")
    ax3.scatter([], [], s=60, marker="v", color=PLANET, label="shorter")
    ax3.scatter([], [], s=60, marker="v", facecolor=PAPER, edgecolor=MUTED, lw=1.3, label="heliocentric")
    ax3.set_yscale("log")
    ax3.set_ylim(3e-7, 0.15)
    ax3.set_xticks(pos, names)
    ax3.set_ylabel(r"$|\,T\,/\,2\pi\sqrt{a^3/GM} - 1\,|$")
    ax3.set_title("What the two-body law leaves: the Sun's pull on the Moon, Jupiter and Saturn's great inequality; heliocentric, the planets inside each orbit", loc="left", fontsize=13, color=INK)
    ax3.legend(loc="upper left", ncol=3, fontsize=11)
    source(fig, "Periods: the slope of the J2000 ecliptic longitude (20 revolutions or more) or its mean return time; a: the time mean of the osculating axis. GM: DE440 (Park et al. 2021). DE440s from NAIF, retrieved 2026-10-21.")
    save(fig, "search-10-kepler-law")


# ---- search-11: the front of formulas ------------------------------------------------------------------------------

def front():
    rows = [r for r in fronts() if r["problem"] == "two centres" and r["units"] == "au, yr"]
    law_err = runs()["law_error"]["two centres"]
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], left=0.07, right=0.98, top=0.84, bottom=0.07, wspace=0.08, hspace=0.22)
    axes = [fig.add_subplot(gs[0, 0])]
    axes.append(fig.add_subplot(gs[0, 1], sharey=axes[0]))
    tables = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    for t in tables:
        t.axis("off")
    title(fig, "The front of formulas", "Two centres, in Kepler's units (au, year): each formula is the most accurate of its size met in 30 runs; the non-dominated set of error and size")
    # The formulas listed in each panel's empty corner, numbered at their points.
    picks = {True: [5, 7, 8, 10, 12, 15, 18], False: [1, 3, 5, 8, 10, 12, 21, 22]}
    for ax, tab, typed in zip(axes, tables, (True, False)):
        pts = sorted([r for r in rows if r["typed"] == typed], key=lambda r: r["size"])
        xs = [r["size"] for r in pts]
        ys = [r["error"] for r in pts]
        colour = TYPED if typed else UNTYPED
        ax.step(xs, ys, where="post", color=colour, lw=1.0, alpha=0.6)
        ax.scatter(xs, ys, s=45, color=colour, zorder=3)
        for r in pts:
            if r["law"]:
                ax.scatter([r["size"]], [r["error"]], s=230, facecolor="none", edgecolor=ALGORITHM["BO"], lw=1.8, zorder=4)
        ax.axhline(law_err, color=MUTED, lw=0.8, ls=(0, (4, 3)))
        ax.text(41.5, law_err * 0.78, f"the law's own error on these data, ${sci(law_err)}$", ha="right", va="top", fontsize=11, color=MUTED)
        ax.set_yscale("log")
        ax.set_xlim(0, 42)
        ax.set_ylim(3e-4, 2.5)
        ax.set_xlabel("size of the formula (nodes)", labelpad=2)
        shown = [r for r in pts if r["size"] in picks[typed]]
        tab.text(0.0, 1.0, "size", transform=tab.transAxes, fontsize=11, color=MUTED, va="top")
        tab.text(0.08, 1.0, "formula", transform=tab.transAxes, fontsize=11, color=MUTED, va="top")
        tab.text(0.805 if not typed else 0.80, 1.0, "error", transform=tab.transAxes, fontsize=11, color=MUTED, va="top")
        if not typed:
            tab.text(0.915, 1.0, "in SI", transform=tab.transAxes, fontsize=11, color=MUTED, va="top")
        for k, r in enumerate(shown):
            y = 0.86 - k * 0.118
            ax.annotate(str(r["size"]), (r["size"], r["error"]), xytext=(0, 9), textcoords="offset points", fontsize=10.5, color=colour if typed else INK, ha="center")
            c = ALGORITHM["BO"] if r["law"] else INK
            tab.text(0.0, y, str(r["size"]), transform=tab.transAxes, fontsize=12, color=c, va="top")
            tab.text(0.08, y, formula(r["formula"]), transform=tab.transAxes, fontsize=12, color=c, va="top")
            tab.text(0.805 if not typed else 0.80, y, f"${sci(r['error'])}$", transform=tab.transAxes, fontsize=12, color=c, va="top")
            if not typed:
                e2 = r["error_other_units"]
                tab.text(0.915, y, f"${sci(e2)}$" if e2 is not None else "not finite", transform=tab.transAxes, fontsize=12, color=c, va="top")
        ax.set_title("typed by dimension" if typed else "untyped: every symbol a bare number", loc="left", fontsize=15, color=colour if typed else INK)
    plt.setp(axes[1].get_yticklabels(), visible=False)
    axes[0].set_ylabel("RMS relative error of $T$")
    axes[0].text(0.98, 0.95, r"ringed, and rose in the table: the law, $2\pi\sqrt{a^3/GM}$ in one of its forms", transform=axes[0].transAxes, fontsize=11.5, color=ALGORITHM["BO"], va="top", ha="right")
    source(fig, "Variables a (au), GM (" + AU3_YR2 + ") and the eccentricity e; constants 2, π and c in [−1, 1]; +, −, ×, ÷, square, cube, √. Error: RMS of T̂/T − 1 over the nine planets and the Moon. " + SOURCE.split(";", 1)[1].strip())
    save(fig, "search-11-kepler-front")


# ---- search-12: success over evaluations ---------------------------------------------------------------------------

def success():
    d = runs()
    budget, seeds = d["budget"], d["seeds"]
    grid = np.geomspace(100, budget, 300)
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), dpi=100, sharex=True)
    fig.subplots_adjust(left=0.08, right=0.97, top=0.83, bottom=0.1, hspace=0.25, wspace=0.2)
    title(fig, "Typing by dimension finds the law; the untyped search finds the units", f"Above: the share of {seeds} seeds that have met the law (its form, within 1 % on the data, and still within 1 % after a change of units). Below: the median best error")
    styles = [(True, "au, yr", TYPED, "-", "typed, in au and years"), (True, "SI", TYPED, (0, (1, 2)), "typed, in SI"),
              (False, "au, yr", UNTYPED, "-", "untyped, in au and years"), (False, "SI", ALGORITHM["(1+1)-ES"], "-", "untyped, in SI")]
    for col, problem in enumerate(("two centres", "one centre")):
        ax, bx = axes[0, col], axes[1, col]
        for typed, units, colour, ls, label in styles:
            rs = [r for r in d["runs"] if r["problem"] == problem and r["units"] == units and r["typed"] == typed]
            hits = np.array([r["law_at"] if r["law_at"] is not None and r["law_at"] >= 0 else np.inf for r in rs])
            share = [(hits <= g).mean() for g in grid]
            ax.plot(grid, share, color=colour, ls=ls, lw=2.4 if ls == "-" else 2.0, label=label)
            # median best error, the progress being a staircase of improvements
            curves = []
            for r in rs:
                pr = np.array(r["progress"]) if r["progress"] else np.zeros((0, 2))
                y = np.full(grid.shape, np.nan)
                for k, g in enumerate(grid):
                    m = pr[:, 0] <= g if len(pr) else []
                    if len(pr) and m.any():
                        y[k] = pr[m, 1][-1]
                curves.append(y)
            med = np.nanmedian(np.array(curves), axis=0)
            bx.plot(grid, med, color=colour, ls=ls, lw=2.4 if ls == "-" else 2.0)
        ax.set_xscale("log")
        ax.set_ylim(-0.03, 1.05)
        ax.set_ylabel("share of seeds")
        ax.set_title({"two centres": "two centres: the planets and the Moon, variables a, GM, e", "one centre": "one centre: the planets alone, variables a, e and a constant k"}[problem], loc="left", fontsize=14, color=INK)
        bx.set_xscale("log"), bx.set_yscale("log")
        bx.axhline(d["law_error"][problem], color=MUTED, lw=0.8, ls=(0, (4, 3)))
        bx.text(2.2e4, d["law_error"][problem] * 1.25, "the law's error on the data", fontsize=11, color=MUTED, ha="right", va="bottom")
        bx.set_xlabel("evaluations")
        bx.set_ylabel("median best RMS error")
    axes[0, 0].legend(loc="center right")
    source(fig, "Typed, the two systems of units meet the law at the same evaluation on 59 of 60 seeds (its terms are covariant; rounding breaks some ties). Population 500, tournaments of 7, lexicographic parsimony. " + SOURCE.split(";", 1)[1].strip())
    save(fig, "search-12-kepler-success")


# ---- search-13: units -----------------------------------------------------------------------------------------------

def units():
    rows = fronts()
    fig = plt.figure(figsize=(16, 10), dpi=100)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1], left=0.075, right=0.97, top=0.83, bottom=0.12, wspace=0.24)
    ax, bx = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    title(fig, "A law, or a coincidence of the units", "Every formula of the fronts, evaluated again after the units change (au and years to SI, or back), its constants converted by their declared dimension")
    lim = (2e-4, 3)
    for typed, colour, label in ((True, TYPED, "typed by dimension"), (False, UNTYPED, "untyped")):
        for problem, marker in (("two centres", "o"), ("one centre", "s")):
            rs = [r for r in rows if r["typed"] == typed and r["problem"] == problem]
            x = np.array([r["error"] for r in rs])
            y = np.array([r["error_other_units"] if r["error_other_units"] is not None else 1e30 for r in rs])
            ax.scatter(x, np.minimum(y, 1e19), s=50, marker=marker, color=colour, alpha=0.85, label=f"{label}, {problem}", zorder=3)
    ax.plot(lim, lim, color=HAIRLINE, lw=1.2)
    ax.text(0.03, 0.0025, "on the diagonal: unchanged, a law", fontsize=12, color=MUTED, ha="left")
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlim(*lim)
    ax.set_ylim(1e-4, 1e20)
    ax.set_xlabel("RMS relative error in the units searched")
    ax.set_ylabel("RMS relative error after the change of units")
    ax.set_title("the error before and after", loc="left", fontsize=14, color=INK)
    ax.legend(loc="upper right", fontsize=11, ncol=1)
    # The GM read from the one-centre typed front, against the planets' own μ in Jacobi coordinates.
    d = runs()
    gm = d["gm_sun_au3_yr2"]
    year, au = 365.25 * 86400.0, 149597870700.0
    mus = [(o["body"], float(o["mu_m3s2"]) * year**2 / au**3) for o in orbits() if o["body"] != "Moon"]
    rs = sorted([r for r in rows if r["typed"] and r["problem"] == "one centre" and r["units"] == "au, yr" and r["law"]], key=lambda r: r["size"])
    lo, hi = min(m for _, m in mus), max(m for _, m in mus)
    bx.axvspan(lo, hi, color=HAIRLINE, alpha=0.45, lw=0)
    for name, m in mus:
        bx.plot([m, m], [-1.05, -0.75], color=PLANET, lw=1.2)
    bx.text(lo - 0.002, -1.2, "Mercury to Mars", fontsize=10.5, color=PLANET, ha="right", va="center")
    bx.text(hi + 0.002, -1.2, "Jupiter to Pluto", fontsize=10.5, color=PLANET, ha="left", va="center")
    for i, r in enumerate(rs):
        x = max(r["implied"], 39.447)
        bx.scatter([x], [i], s=60, color=TYPED, zorder=3, marker="o" if x == r["implied"] else "<")
        if x != r["implied"]:
            bx.text(x + 0.003, i, f"{r['implied']:.2f}", va="center", fontsize=11, color=TYPED)
        bx.text(39.582, i, f"size {r['size']}", ha="right", va="center", fontsize=11, color=MUTED)
    bx.axvline(gm, color=MOON, lw=1.4)
    bx.text(gm - 0.002, len(rs) - 0.35, "the Sun's GM (DE440)\n" + f"{gm:.4f} " + AU3_YR2, color=MOON, fontsize=11.5, va="bottom", ha="right")
    bx.text(hi + 0.002, len(rs) - 0.35, "the planets' μ:\nthe Sun and the\nplanets inside", color=PLANET, fontsize=11.5, va="bottom", ha="left")
    bx.set_yticks([])
    bx.set_ylim(-1.9, len(rs) + 1.0)
    bx.set_xlim(39.44, 39.585)
    bx.set_xlabel(r"the GM a formula implies, $4\pi^2 a^3/T(a)^2$ at $a$ = 1 au (" + AU3_YR2 + ")")
    bx.set_title("one centre, typed: the constant k becomes a GM", loc="left", fontsize=14, color=INK)
    source(fig, "An untyped constant is a bare number and cannot change with the units; a typed one carries its dimension. " + SOURCE.split(";", 1)[1].strip())
    save(fig, "search-13-kepler-units")


if __name__ == "__main__":
    law()
    front()
    success()
    units()
