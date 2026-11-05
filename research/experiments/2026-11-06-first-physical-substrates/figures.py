"""The figures of "The first physical substrates": one example of each simulator of the physical-computing domain's
fluids, chemical and molecular families (ADR 0076), drawn from results/ in the house style of the research records.

    .venv/Scripts/python.exe research/experiments/2026-11-06-first-physical-substrates/figures.py [--only 1,5]

Writes assets/gallery/physcomp-01-... to physcomp-16-... (1600 x 1000).
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F401,F403

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, ListedColormap  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"

# A family of colours for the three families: water's blues and sand, chemistry's plum and amber, DNA's teal and sage.
SEA, DEEP, SAND, FOAM = "#5B8FB9", "#1F3B57", "#C98B6B", "#A9C3DA"
PLUM, ORCHID, AMBER, ROSE = "#6B3E75", "#B5629E", "#D49A3A", "#D9738B"
TEAL, SAGE, MOSS = "#4E8F8B", "#7FA37A", "#2F6E6A"
WATER = LinearSegmentedColormap.from_list("water", ["#1F3B57", "#5B8FB9", "#F4EEE8", "#E3A985", "#A85A3A"])
BZ = LinearSegmentedColormap.from_list("bz", ["#FBF6F1", "#EBD3E3", "#CA8DBE", "#8E3B78", "#3B1F3A"])
CRED = "physical-computing domain (ADR 0076), its simulators on physis; research/experiments/2026-11-06-first-physical-substrates"


def load(name):
    return json.loads((RESULTS / f"{name}.json").read_text(encoding="utf-8"))


def arr(values):
    """An array with JSON's nulls (thresholds never reached, written as null) as NaN, which matplotlib leaves out."""
    return np.array([np.nan if v is None else v for v in values], dtype=float)


def panel_title(ax, text):
    ax.set_title(text, loc="left", fontsize=14, color=INK, pad=10)


# ---- 1. convection -----------------------------------------------------------------------------------------------------
def fig_convection():
    d = load("convection")
    fig, (a, b) = figure(1, 2)
    fig.subplots_adjust(left=0.07, right=0.97, top=0.80, bottom=0.11, wspace=0.25)
    title(fig, "A layer heated from below: the onset of convection, and Lorenz's chaos in its rolls",
          "The marginal Rayleigh number by shooting, for three pairs of plates; the cell's heat flux against its heating, by Lorenz's truncation")
    for key, colour, label in [("rigid_rigid", DEEP, "rigid plates"), ("rigid_free", SEA, "rigid below, free above"), ("free_free", SAND, "free plates")]:
        a.plot(d["a"], d[key], color=colour, label=label)
        k, r = d[key + "_critical"]
        a.plot([k], [r], "o", color=colour, ms=8)
        a.annotate(f"Ra$_c$ = {r:.2f}\n at a = {k:.3f}", (k, r), textcoords="offset points", xytext=(8, -34), fontsize=11, color=colour)
    a.set_yscale("log"), a.set_xlabel("horizontal wavenumber a (times the depth)"), a.set_ylabel("marginal Rayleigh number")
    a.set_ylim(500, 20000), a.legend(loc="upper center")
    panel_title(a, "Onset: Chandrasekhar's 1707.762, 1100.65 and 27π⁴/4")
    r = np.array(d["r"])
    b.fill_between(r, d["nusselt_min"], d["nusselt_max"], color=FOAM, alpha=0.8, lw=0, label="range of the Nusselt number")
    b.plot(r, d["nusselt_max"], color=DEEP, lw=0.8)
    b.plot(r, d["nusselt_min"], color=DEEP, lw=0.8)
    rr = np.linspace(1, d["hopf_r"], 100)
    b.plot(rr, 1 + 2 * (rr - 1) / rr, color=SAND, lw=2.2, label="steady rolls, 1 + 2(r − 1)/r")
    b.axvline(d["hopf_r"], color=PLUM, lw=1, ls="--")
    b.set_ylim(0.9, 4.95)
    b.text(d["hopf_r"] - 0.4, 4.62, f"Hopf point r$_H$ = {d['hopf_r']:.2f} (σ = Pr = {d['prandtl']:.2f})", color=PLUM, fontsize=11, ha="right")
    b.set_xlabel(f"heating r = ΔT / ΔT$_c$   (ΔT$_c$ = {d['dTc']:.4f} K for 1 cm of water)"), b.set_ylabel("Nusselt number (heat carried / conducted)")
    b.legend(loc="lower center")
    panel_title(b, "The cell's states: conduction, steady rolls, then chaos")
    ins = b.inset_axes([0.04, 0.55, 0.34, 0.34])
    ins.plot(d["orbit_x"], d["orbit_z"], color=DEEP, lw=0.4)
    ins.set_xticks([]), ins.set_yticks([]), ins.set_title("rolls' circulation X against Z at r = 28", fontsize=10, color=MUTED)
    source(fig, CRED + "; physis::continuum::critical_rayleigh, lorenz_parameters; fluids::ConvectionCell")
    save(fig, "physcomp-01-a-layer-heated-from-below")


# ---- 2. rotating tank -------------------------------------------------------------------------------------------------
def fig_rotating():
    d = load("rotating")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The general circulation in a tank: baroclinic waves in a rotating annulus of water",
          f"Two-layer quasi-geostrophic flow in the unrolled annulus (Θ = {d['theta']:.2f}, Ta = {d['taylor']:.2e}); from noise, its most unstable wavenumber m = {int(d['most_unstable'])} takes over")
    ax = fig.add_axes([0.03, 0.08, 0.50, 0.78])
    nx, ny = int(d["nx"]), int(d["ny"])
    psi = np.array(d["psi1"]).reshape(ny, nx)[: ny // 2 + 1, :]  # the channel: y in [0, L]
    R, L = d["radius"], d["gap"]
    th = np.linspace(0, 2 * np.pi, nx + 1)
    rad = np.linspace(R - L / 2, R + L / 2, psi.shape[0])
    T, Rr = np.meshgrid(th, rad)
    P = np.concatenate([psi, psi[:, :1]], axis=1)
    lim = np.abs(P).max()
    ax.pcolormesh(Rr * np.cos(T), Rr * np.sin(T), P, cmap=WATER, vmin=-lim, vmax=lim, shading="gouraud")
    ax.add_patch(plt.Circle((0, 0), R - L / 2, color=PAPER, zorder=3))
    ax.add_patch(plt.Circle((0, 0), R + L / 2, fill=False, color=INK, lw=1, zorder=3))
    ax.add_patch(plt.Circle((0, 0), R - L / 2, fill=False, color=INK, lw=1, zorder=4))
    ax.set_aspect("equal"), ax.axis("off")
    ax.text(0, 0, "cold\ninner wall", ha="center", va="center", color=SEA, fontsize=12, zorder=5)
    ax.text(0, -(R + L / 2) * 1.08, "warm outer wall", ha="center", color=SAND, fontsize=12)
    ax.set_title(f"the upper layer's streamfunction: {int(d['dominant'])} waves", color=INK, fontsize=14, loc="left")
    g = fig.add_axes([0.60, 0.52, 0.37, 0.32])
    for i, (th_, colour) in enumerate(zip(d["thetas"], [DEEP, SEA, SAND, ROSE])):
        g.plot(d["m"], d["growth_inviscid"][i], "-", color=colour, label=f"Θ = {th_}")
        g.plot(d["m"], d["growth_ekman"][i], "--", color=colour, lw=1.2)
    g.set_xlabel("azimuthal wavenumber m"), g.set_ylabel("linear growth rate (1/s)"), g.legend(ncol=2, fontsize=10)
    g.set_title("growth by wavenumber: inviscid (solid), with Ekman friction (dashed)", loc="left", fontsize=12)
    e = fig.add_axes([0.60, 0.10, 0.37, 0.30])
    e.semilogy(d["time"], d["energy"], color=DEEP)
    e.set_xlabel("time (s)"), e.set_ylabel("energy (m²/s²)")
    e.set_title(f"from noise: growth at {d['rate']:.3f}/s, then saturation", loc="left", fontsize=12)
    source(fig, CRED + "; physis::continuum::TwoLayerQG (channel by images); fluids::RotatingTank (Hide and Mason 1975)")
    save(fig, "physcomp-02-the-general-circulation-in-a-tank")


# ---- 3. Faraday -------------------------------------------------------------------------------------------------------
def fig_faraday():
    d = load("faraday")
    fig, (a, b) = figure(1, 2)
    fig.subplots_adjust(left=0.07, right=0.93, top=0.80, bottom=0.11, wspace=0.3)
    title(fig, "A layer of water on a shaker: the Faraday instability, mode by mode",
          "Each surface mode a damped Mathieu oscillator; its threshold from the Floquet multipliers, its response at half the drive's frequency")
    k = np.array(d["k"]) / d["resonant_k"]
    a.plot(k, arr(d["threshold"]) / 9.80665, color=SEA, label="Floquet threshold (none below 4 g where absent)")
    tip = arr(d["small_damping"])[np.argmin(np.abs(k - 1))] / 9.80665
    a.plot([1], [tip], "D", color=SAND, ms=9, label="the tip for small damping, 4γω₀/(k tanh kh)")
    a.axvline(1, color=MUTED, lw=0.8)
    a.set_xlabel("wavenumber k / k$_r$  (k$_r$: ω₀ = Ω/2)"), a.set_ylabel("threshold acceleration (g)"), a.legend()
    panel_title(a, "The first tongue at 40 Hz: its tip at resonance")
    b.plot(d["drive_hz"], arr(d["onset_acceleration"]) / 9.80665, "o-", color=DEEP, label="onset acceleration")
    b.set_xlabel("drive frequency (Hz)"), b.set_ylabel("onset acceleration (g)")
    c = b.twinx()
    c.plot(d["drive_hz"], arr(d["onset_wavelength"]) * 1000, "s--", color=SAND, label="wavelength at onset")
    c.set_ylabel("wavelength at onset (mm)", color=SAND), c.grid(False)
    b.legend(loc="upper left"), c.legend(loc="upper right")
    panel_title(b, "Onset against the drive (bulk damping only: lower than a real dish)")
    source(fig, CRED + "; physis::dynamics::parametric_oscillator, physis::continuum::gravity_capillary_frequency; fluids::FaradayLayer")
    save(fig, "physcomp-03-a-layer-of-water-on-a-shaker")


# ---- 4. walkers -------------------------------------------------------------------------------------------------------
def fig_walkers():
    d = load("walkers")
    fig, (a, b) = figure(1, 2)
    fig.subplots_adjust(left=0.07, right=0.97, top=0.80, bottom=0.11, wspace=0.25)
    title(fig, "A droplet walking on its own waves: hydrodynamic pilot-wave dynamics",
          "Oza, Rosales and Bush's stroboscopic model: walking begins at β = 2, at u = √(2β − 1 − √(1 + 4β))/√2")
    a.plot(d["beta_theory"], d["speed_theory"], color=SAND, label="steady speed, closed form")
    a.plot(d["beta"], d["speed"], "o", color=DEEP, ms=8, label="simulated, after 80 memory times")
    a.axvline(2, color=MUTED, lw=0.8, ls="--")
    a.set_xlabel("memory β"), a.set_ylabel("walking speed"), a.legend()
    panel_title(a, "A pitchfork from rest at β = 2")
    b.plot(d["orbit_x"], d["orbit_y"], color=TEAL, lw=0.7)
    b.plot(d["orbit_x"][0], d["orbit_y"][0], "o", color=SAND)
    b.set_aspect("equal"), b.set_xlabel("x (1/k$_F$)"), b.set_ylabel("y (1/k$_F$)")
    panel_title(b, f"In a harmonic well (k = {d['well']}), β = {d['orbit_beta']}: the walker's orbits")
    source(fig, CRED + "; fluids::Walker (Oza, Rosales and Bush 2013); quantum-like statistics are claims to be checked (Andersen et al. 2015)")
    save(fig, "physcomp-04-a-droplet-walking-on-its-waves")


# ---- 5. tank ----------------------------------------------------------------------------------------------------------
def fig_tank():
    d = load("tank")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "A tank of water as a reservoir: two motors' waves, a camera, and a memory that fades",
          "The damped long-wave equation of physis's Basin; the camera's pixels are the reservoir's state, read by a linear readout")
    nx, ny = int(d["nx"]), int(d["ny"])
    s = np.array(d["surface"]).reshape(ny, nx)
    ax = fig.add_axes([0.04, 0.42, 0.44, 0.42])
    lim = np.abs(s).max()
    ax.imshow(s, cmap=WATER, vmin=-lim, vmax=lim, origin="lower", extent=[0, d["length"], 0, d["width"]])
    ax.plot([0.08, 0.32], [0.15, 0.15], "o", color=INK, ms=6)
    ax.set_title("the surface after 14 s (motors at 2.0 and 2.6 Hz, black)", loc="left", fontsize=12)
    ax.set_xlabel("x (m)"), ax.set_ylabel("y (m)")
    cx, cy = int(d["cx"]), int(d["cy"])
    c = fig.add_axes([0.54, 0.42, 0.42, 0.42])
    c.imshow(np.array(d["binary"]).reshape(cy, cx), cmap=ListedColormap([PANEL, DEEP]), origin="lower")
    c.set_xticks([]), c.set_yticks([])
    c.set_title(f"what the camera sees: {cx} × {cy} pixels, thresholded", loc="left", fontsize=12)
    e = fig.add_axes([0.07, 0.12, 0.88, 0.20])
    t, g = np.array(d["esp_time"]), np.array(d["esp_gap"])
    e.semilogy(t, g / g[0], color=DEEP, label="gap between two runs under the same inputs")
    e.semilogy(t, np.exp(-d["damping"] * t), color=SAND, ls="--", label="e$^{−γt}$")
    e.set_xlabel("time (s)"), e.set_ylabel("gap / first gap"), e.legend(loc="upper right")
    e.set_title("the echo state property: the start is forgotten at the damping rate (its volume, the mean level, never is)", loc="left", fontsize=12)
    source(fig, CRED + "; physis::continuum::Basin; fluids::WaveTank (Fernando and Sojakka 2003)")
    save(fig, "physcomp-05-a-tank-of-water-as-a-reservoir")


# ---- 6. capillary maze ------------------------------------------------------------------------------------------------
def fig_capillary():
    d = load("capillary")
    W = int(d["width"])
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "A maze solved by flow: capillaries carry liquid only along the way through",
          "A perfect maze of 20 × 20 cells as tubes of Hagen–Poiseuille conductance; one Kirchhoff solve, from the entrance to the exit")
    ax = fig.add_axes([0.25, 0.06, 0.50, 0.80])
    q = np.abs(np.array(d["flow"]))
    peak = q.max()
    segs, cols, widths = [], [], []
    for a, b, f in zip(d["a"], d["b"], q):
        a, b = int(a), int(b)
        segs.append([(a % W, a // W), (b % W, b // W)])
        on = f > 1e-9 * peak
        cols.append(SEA if on else HAIRLINE)
        widths.append(5 if on else 2)
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=widths, capstyle="round"))
    ax.plot([0], [0], "o", color=SAND, ms=14), ax.plot([W - 1], [W - 1], "o", color=PLUM, ms=14)
    ax.text(-0.6, -0.9, "entrance", color=SAND, fontsize=12), ax.text(W - 2.5, W - 0.2, "exit", color=PLUM, fontsize=12)
    ax.set_xlim(-1, W), ax.set_ylim(-1.5, W + 0.5), ax.set_aspect("equal"), ax.axis("off")
    n_on = int((q > 1e-9 * peak).sum())
    fig.text(0.78, 0.5, f"{n_on} of {len(q)} tubes carry flow:\nexactly the solution.\n\nA dead end cannot\nconserve a non-zero flow,\nso it carries none\n(|Q| < $10^{-9}$ of the peak).",
             fontsize=13, color=INK, va="center")
    source(fig, CRED + "; physis::discretisation::solve_kirchhoff; fluids::CapillaryNetwork (after Fuerstman et al. 2003)")
    save(fig, "physcomp-06-a-maze-solved-by-flow")


# ---- 7. foam ----------------------------------------------------------------------------------------------------------
def fig_foam():
    d = load("foam")
    n = int(d["size"])
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "A foam coarsening: cells of few sides vanish, cells of many grow",
          "The large-q Potts model of a two-dimensional foam, from a Voronoi tessellation; von Neumann's law, dA/dt = k (n − 6), measured")
    rng = np.random.default_rng(2)
    palette = np.array([to for to in [PANEL, "#EBD3E3", "#E2D8CE", "#DCE6EF", "#E8EEF3", "#F1E3D3", "#E3EDE0"]])
    for i, key in enumerate(["s0", "s1", "s2"]):
        ax = fig.add_axes([0.03 + 0.235 * i, 0.30, 0.21, 0.55])
        s = np.array(d[key]).reshape(n, n).astype(int)
        colours = rng.integers(0, len(palette), int(s.max()) + 1)
        img = np.zeros((n, n, 3))
        from matplotlib.colors import to_rgb
        lut = np.array([to_rgb(palette[c]) for c in colours])
        img = lut[s]
        edge = (s != np.roll(s, 1, 0)) | (s != np.roll(s, 1, 1))
        img[edge] = to_rgb(PLUM)
        ax.imshow(img, origin="lower"), ax.set_xticks([]), ax.set_yticks([])
        label = {"s0": f"start: {int(d['cells0'])} cells", "s1": f"{int(d['sweeps1'])} sweeps: {int(d['cells1'])} cells", "s2": f"{int(d['sweeps2'])} sweeps: {int(d['cells2'])} cells"}[key]
        ax.set_title(label, loc="left", fontsize=12)
    v = fig.add_axes([0.78, 0.30, 0.19, 0.55])
    sides, rate, count = np.array(d["sides"]), np.array(d["rate"]), np.array(d["count"])
    keep = count >= 5
    v.scatter(sides[keep], rate[keep], s=20 + 4 * count[keep], color=SEA, alpha=0.8)
    if keep.sum() > 2:
        w = count[keep]
        p = np.polyfit(sides[keep], rate[keep], 1, w=np.sqrt(w))
        xs = np.linspace(sides[keep].min(), sides[keep].max(), 10)
        v.plot(xs, np.polyval(p, xs), color=SAND)
        v.text(0.05, 0.92, f"zero at n = {-p[1] / p[0]:.2f}", transform=v.transAxes, color=SAND, fontsize=12)
    v.axhline(0, color=MUTED, lw=0.8), v.axvline(6, color=MUTED, lw=0.8, ls="--")
    v.set_xlabel("sides n"), v.set_ylabel("mean dA/dt (sites per sweep)")
    v.set_title("von Neumann's law", loc="left", fontsize=12)
    fig.text(0.04, 0.13, "A cell's growth depends on its number of sides alone: below six it shrinks and vanishes, above six it grows, whatever its shape (von Neumann 1952).\n"
             "Here the law is statistical, read from the cells alive at 40 sweeps (vanished cells counted with area zero), so the line crosses zero near, not at, six;\n"
             "the lattice's anisotropy is reduced by the Moore neighbourhood and a small temperature, not removed.", fontsize=12, color=INK)
    source(fig, CRED + "; physis::statistical_mechanics::PottsLattice; fluids::Foam (Glazier, Anderson and Grest 1990)")
    save(fig, "physcomp-07-a-foam-coarsening")


# ---- 8. bubbles and Phillips ------------------------------------------------------------------------------------------
def fig_logic_and_economy():
    b = load("bubbles")
    p = load("phillips")
    fig, (a, c) = figure(1, 2)
    fig.subplots_adjust(left=0.06, right=0.97, top=0.80, bottom=0.12, wspace=0.25)
    title(fig, "Water as logic and as an economy: a bubble gate, and Phillips's hydraulic machine",
          "Bubbles that raise their channel's resistance divert the next (AND and OR in one junction); a tank whose outflow is income")
    t = np.array(b["truth"]).reshape(4, 4)
    labels = [f"A={int(r[0])}, B={int(r[1])}" for r in t]
    xs = np.arange(4)
    a.bar(xs - 0.18, t[:, 2], 0.34, color=SEA, label="bubbles at the OR outlet")
    a.bar(xs + 0.18, t[:, 3], 0.34, color=SAND, label="bubbles at the AND outlet")
    a.set_xticks(xs, labels), a.set_ylabel("bubbles out"), a.set_ylim(0, 1.6), a.legend(loc="upper left")
    panel_title(a, "The truth table of one junction: a lone bubble takes OR; two together take both")
    tt = np.array(p["time"]) / p["tau"]
    c.plot(tt, np.array(p["income"]) * 1e4, color=DEEP, label="income Y = h/R")
    inj = np.array(p["injections"])
    c.plot(tt, inj / (1 - p["propensity"]) * 1e4, color=SAND, ls="--", label="the multiplier's I/(1 − c)")
    c.plot(tt, inj * 1e4, color=MUTED, lw=1, label="injections I")
    c.set_xlabel("time (time constants A R/(1 − c))"), c.set_ylabel("flow ($10^{-4}$ m³/s)"), c.legend()
    panel_title(c, f"Phillips's machine (c = {p['propensity']}): income settles at the multiplier")
    source(fig, CRED + "; fluids::BubbleCircuit (provisional; Prakash and Gershenfeld 2007); fluids::PhillipsMachine (Phillips 1950)")
    save(fig, "physcomp-08-water-as-logic-and-as-an-economy")


# ---- 9. Oregonator ----------------------------------------------------------------------------------------------------
def fig_oregonator():
    d = load("oregonator")
    fig, (a, b) = figure(1, 2)
    fig.subplots_adjust(left=0.07, right=0.97, top=0.80, bottom=0.11, wspace=0.25)
    title(fig, "The Belousov–Zhabotinsky medium, stirred: a window of oscillation between two Hopf points",
          f"The two-variable Oregonator (ε = {d['epsilon']}, q = {d['q']}); its steady state, its Hopf points by two methods, its oscillation's range")
    f, us, st = np.array(d["f"]), np.array(d["steady"]), np.array(d["stable"]).astype(bool)
    a.fill_between(f, d["u_min"], d["u_max"], where=~st, color="#EBD3E3", lw=0, label="oscillation's range of u")
    a.semilogy(f[st], us[st], ".", color=PLUM, ms=5, label="steady state, stable")
    a.semilogy(f[~st], us[~st], ".", color=ROSE, ms=3, label="steady state, unstable")
    for h in d["hopf"]:
        a.axvline(h, color=AMBER, lw=1.2, ls="--")
        a.text(h + 0.03, 0.6, f"Hopf\nf = {h:.4f}", color=AMBER, fontsize=11)
    a.set_xlabel("stoichiometric factor f"), a.set_ylabel("activator u (HBrO₂, scaled)"), a.legend(loc="lower left")
    panel_title(a, "Oscillatory between the Hopf points, excitable beyond")
    b.semilogy(d["t"], d["u"], color=PLUM, label="activator u")
    b.semilogy(d["t"], d["v"], color=AMBER, label="oxidised catalyst v")
    b.set_xlabel("time (scaled)"), b.set_ylabel("concentration (scaled)"), b.legend()
    panel_title(b, "A relaxation oscillation at f = 1")
    source(fig, CRED + "; chemical::Oregonator (Tyson and Fife 1980); Hopf points by the trace and by physis::dynamics::continue_equilibria")
    save(fig, "physcomp-09-the-bz-medium-stirred")


def fig_spiral():
    d = load("spiral")
    n = int(d["n"])
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "The BZ medium in a dish: a broken wave curls about its free end and keeps turning",
          f"Excitable Oregonator kinetics (ε = {d['epsilon']}, f = {d['f']:.3f}) in a {n} × {n} dish with no-flux walls, the catalyst immobile; stepped implicitly in the stiff kinetics")
    top = 0.9 * max(max(f) for f in d["frames"])
    for i, (frame, t) in enumerate(zip(d["frames"], d["times"])):
        ax = fig.add_axes([0.02 + 0.245 * i, 0.20, 0.225, 0.60])
        ax.imshow(np.array(frame).reshape(n, n), cmap=BZ, vmin=0, vmax=top, origin="lower")
        ax.set_xticks([]), ax.set_yticks([])
        ax.set_title(f"t = {t:.0f}: {int(d['excited'][i])} cells excited", loc="left", fontsize=13)
    fig.text(0.04, 0.07, "Activator u, dark where the medium is excited. The wave was cut at the middle with a refractory region behind it; its free end cannot advance into that\n"
             "region and turns about it. A plane wave would leave through the walls; this one keeps turning to t = 40. The dish is about one wavelength across, so it\n"
             "holds one turn of the spiral: a spiral needs a medium wider than its wavelength, some T/ε front widths, several times this dish at ε = 0.04.", fontsize=12, color=INK)
    source(fig, CRED + "; chemical::BZMedium::dish on physis::chemistry::FieldKinetics::step_split (Strang splitting, SDIRK2)")
    save(fig, "physcomp-10-a-broken-wave-curls-into-a-spiral")


def fig_array():
    d = load("array")
    fig, (a, b) = figure(2, 1, sharex=False)
    fig.subplots_adjust(left=0.07, right=0.97, top=0.82, bottom=0.08, hspace=0.45)
    title(fig, "Two cells of the BZ medium, 5% apart in their recipes: uncoupled they drift, coupled they lock",
          "Stirred Oregonators exchanging the activator through their walls: the programmable chemical computer's coupled cells")
    for ax, key, label in [(a, "uncoupled", "uncoupled: each keeps its own period, and their peaks drift apart"), (b, "coupled", "coupled through the activator: one period, a fixed lag")]:
        ax.plot(d[key + "_t"], d[key + "_a"], color=PLUM, label="cell 1")
        ax.plot(d[key + "_t"], d[key + "_b"], color=AMBER, label="cell 2 (f × 1.05)")
        ax.set_ylabel("catalyst v"), ax.legend(loc="upper right")
        panel_title(ax, label)
    b.set_xlabel("time (scaled)")
    source(fig, CRED + "; chemical::BZArray (after Parrilla-Gutierrez et al. 2020)")
    save(fig, "physcomp-11-two-cells-lock")


def fig_programs():
    d = load("programs")
    fig, (a, b) = figure(1, 2)
    fig.subplots_adjust(left=0.07, right=0.96, top=0.80, bottom=0.11, wspace=0.28)
    title(fig, "Reaction networks as programs: Lorenz's flow on dual rails, and max(x₁, x₂) with probability one",
          f"A polynomial ODE compiled to {int(d['reactions'])} reactions of mass action; a function on counts computed exactly by its Markov chain")
    t = np.array(d["t"])
    a.plot(t, d["x_plus"], color=PLUM, lw=1.2, label="rail $x^{+}$")
    a.plot(t, d["x_minus"], color=AMBER, lw=1.2, label="rail $x^{-}$")
    a.plot(t, np.array(d["x_plus"]) - np.array(d["x_minus"]), color=DEEP, lw=2.4, label="$x^{+} - x^{-}$")
    a.plot(t, d["x_lorenz"], color=SAND, lw=1.2, ls="--", label="Lorenz's x, integrated directly")
    a.set_xlabel("time"), a.set_ylabel("concentration"), a.legend(loc="upper left", fontsize=10)
    panel_title(a, f"The difference carries the flow (largest gap {max(d['error']):.1e} to t = 4)")
    tab = np.array(d["max_table"]).reshape(-1, 4)
    grid = np.zeros((5, 5))
    for x1, x2, y, prob in tab:
        grid[int(x2), int(x1)] = y
    b.imshow(grid, cmap=BZ, origin="lower")
    for x1, x2, y, prob in tab:
        b.text(x1, x2, f"{int(y)}\n(p = {prob:.0f})", ha="center", va="center", fontsize=10, color=PAPER if y > 2 else INK)
    b.set_xlabel("x₁ molecules"), b.set_ylabel("x₂ molecules"), b.grid(False)
    panel_title(b, "max(x₁, x₂): the output's whole law, exactly (rationals)")
    source(fig, CRED + "; physis::chemistry::dual_rail, absorption; chemical::max_program (Chen, Doty and Soloveichik 2014)")
    save(fig, "physcomp-12-reaction-networks-as-programs")


def fig_electrodes():
    d = load("electrodes")
    fig, ax = figure()
    fig.subplots_adjust(left=0.08, right=0.96, top=0.80, bottom=0.11)
    title(fig, "Six hundred dissolving electrodes: Kuramoto's transition to synchrony",
          "The electrodes' phase reduction, natural frequencies of Lorentzian spread γ, coupled through a common resistance")
    k = np.linspace(0.2, 3.0, 300)
    ax.plot(k, np.where(k > 1, np.sqrt(np.clip(1 - 1 / k, 0, None)), 0), color=SAND, lw=2.4, label="r = √(1 − 2γ/K), as N grows")
    ax.plot(d["K_over_Kc"], d["r"], "o", color=PLUM, ms=9, label="simulated, N = 600 (mean over the last 100 time units)")
    ax.axvline(1, color=MUTED, lw=0.8, ls="--")
    ax.set_xlabel("coupling K / K$_c$   (K$_c$ = 2γ)"), ax.set_ylabel("order parameter r"), ax.legend(loc="upper left")
    source(fig, CRED + "; physis::dynamics::Kuramoto; chemical::ElectrodeArray (provisional; Kiss, Zhai and Hudson 2002)")
    save(fig, "physcomp-13-dissolving-electrodes-synchronise")


def fig_strands():
    d = load("strands")
    fig, (a, b) = figure(1, 2)
    fig.subplots_adjust(left=0.07, right=0.97, top=0.80, bottom=0.11, wspace=0.25)
    title(fig, "Any reaction network in DNA: strand displacement approaches its formal network as its fuel grows",
          "X $\\rightarrow$ Y, X + Y $\\rightarrow$ Z compiled to gates, fuels and cascades (Soloveichik, Seelig and Winfree 2010), integrated implicitly")
    a.loglog(d["c_max"], d["error"], "o-", color=TEAL, ms=9, label="largest error to t = 2")
    c = np.array(d["c_max"])
    a.loglog(c, d["error"][0] * c[0] / c, color=MUTED, ls="--", label="slope −1")
    a.set_xlabel("fuel C$_{max}$ (q$_{max}$ = 10 C$_{max}$)"), a.set_ylabel("|DSD − formal|"), a.legend()
    panel_title(a, "Error falls as 1/C$_{max}$: fuel in excess, cascades fast")
    t = d["t"]
    for key, colour, lab in [("x", TEAL, "X"), ("y", SAGE, "Y"), ("z", AMBER, "Z")]:
        b.plot(t, d["formal_" + key], color=colour, lw=2, label=f"{lab}, formal")
        b.plot(t, d["dsd_" + key], "o", color=colour, ms=5, label=f"{lab}, DNA (C$_{{max}}$ = 30)")
    b.set_xlabel("time"), b.set_ylabel("concentration"), b.legend(ncol=2, fontsize=10)
    panel_title(b, "The formal species against their DNA realisation, at a small fuel")
    source(fig, CRED + "; molecular::compile_strand_displacement; physis::chemistry::integrate_stiff")
    save(fig, "physcomp-14-any-network-in-dna")


def fig_seesaw():
    d = load("seesaw")
    fig, (a, b) = figure(1, 2, sharey=True)
    fig.subplots_adjust(left=0.07, right=0.97, top=0.80, bottom=0.11, wspace=0.1)
    title(fig, "Seesaw gates of DNA: a threshold makes OR, a higher one makes AND",
          "An integrating gate sums the inputs; an amplifying gate with a threshold and fuel releases its output above it (Qian and Winfree 2011)")
    labels = ["no input", "x₁", "x₂", "x₁ and x₂"]
    colours = [MUTED, TEAL, SAGE, AMBER]
    for ax, key, name in [(a, "or", "threshold 0.6: OR"), (b, "and", "threshold 1.2: AND")]:
        for curve, lab, col, ls in zip(d[key], labels, colours, ["-", "-", ":", "-"]):
            ax.plot(d["t"], curve, color=col, ls=ls, lw=2.2, label=lab)
        ax.set_xlabel("time"), panel_title(ax, name)
    a.set_ylabel("output strand y (units)"), a.legend(loc="center right")
    source(fig, CRED + "; molecular::two_input_gate, by mass action (provisional)")
    save(fig, "physcomp-15-seesaw-gates")


def fig_pen():
    d = load("pen")
    fig, (a, b) = figure(1, 2)
    fig.subplots_adjust(left=0.07, right=0.97, top=0.80, bottom=0.11, wspace=0.25)
    title(fig, "A molecular ecosystem: DNA prey and predators oscillate past a Hopf point, and the prey invades as a front",
          f"The PEN toolbox's predator–prey in Rosenzweig and MacArthur's form; Hopf at a prey template of K$_H$ = {d['hopf_template']:.3f}")
    a.plot(d["below_prey"], d["below_predator"], color=SAGE, label=f"K = {d['below_K']:.2f} (below K$_H$): a spiral to rest")
    a.plot(d["above_prey"], d["above_predator"], color=TEAL, label=f"K = {d['above_K']:.2f} (above): a limit cycle")
    a.set_xlabel("prey strand N"), a.set_ylabel("predator strand P"), a.legend()
    panel_title(a, "Enriching the prey's template destabilises the ecosystem")
    x = np.arange(len(d["profiles"][0])) * d["spacing"]
    for i, prof in enumerate(d["profiles"]):
        b.plot(x, prof, color=[SAGE, TEAL, MOSS, AMBER, ROSE][i], label=f"t = {5 * (i + 1)}")
    b.set_xlabel("position"), b.set_ylabel("prey N"), b.legend()
    panel_title(b, "The prey alone invades at 2√(D r), less Bramson's correction")
    source(fig, CRED + "; molecular::PenPredatorPrey (Fujii and Rondelez 2013; Padirac et al. 2013); physis::dynamics::RosenzweigMacArthur")
    save(fig, "physcomp-16-a-molecular-ecosystem")


FIGURES = [fig_convection, fig_rotating, fig_faraday, fig_walkers, fig_tank, fig_capillary, fig_foam, fig_logic_and_economy,
           fig_oregonator, fig_spiral, fig_array, fig_programs, fig_electrodes, fig_strands, fig_seesaw, fig_pen]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    only = {int(x) for x in args.only.split(",") if x}
    for i, f in enumerate(FIGURES, 1):
        if only and i not in only:
            continue
        try:
            f()
            print(f"figure {i}: {f.__name__}")
        except FileNotFoundError as e:
            print(f"figure {i}: {f.__name__} skipped, no results ({e.filename})")
