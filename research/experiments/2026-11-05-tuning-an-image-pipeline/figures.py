"""Figures of tuning an image pipeline on phantoms (README.md §7), drawn from results/ only.

    .venv/Scripts/python.exe research/experiments/2026-11-05-tuning-an-image-pipeline/figures.py
"""

import csv
import pathlib
import sys
from collections import OrderedDict, defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "style"))
from search_style import *  # noqa: E402,F403

import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
VERSION = "heuresis 6b3cc88 · eikon 10b108a · root 4153cb9"
HAND = "#C98B6B"  # sand: the hand-set chain


def rows(name):
    with open(RESULTS / name, newline="") as f:
        return list(csv.DictReader(f))


def panel_title(ax, text, sub=None):
    ax.set_title(text, loc="left", fontsize=15, color=INK, pad=24 if sub else 10)
    if sub:
        ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=11, color=MUTED, ha="left", va="bottom")


def nondominated(points):
    out = []
    for i, (a, b) in enumerate(points):
        if not any((c >= a and d >= b and (c > a or d > b)) for c, d in points):
            out.append(i)
    return out


# ---- search-47: the front, and what it does on phantoms it has not seen --------------------------------------

def front_figure():
    ev = rows("evaluations.csv")
    held = rows("held_out.csv")
    ft = rows("front_test.csv")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Tuning a vessel chain for overlap and for centrelines, and what the tuning does not carry",
          "Eight parameters of Frangi's vesselness, hysteresis and growth, on two phantoms; the non-dominated chains then held to three phantoms not seen")
    gs = fig.add_gridspec(1, 2, left=0.065, right=0.97, top=0.80, bottom=0.12, wspace=0.22, width_ratios=[1.15, 1])

    ax = fig.add_subplot(gs[0, 0])
    for alg, colour, z in (("random", ALGORITHM["random"], 1), ("NSGA-II", ALGORITHM["NSGA-II"], 2)):
        r = [x for x in ev if x["algorithm"] == alg]
        d = np.array([float(x["dice"]) for x in r])
        c = np.array([float(x["cldice"]) for x in r])
        ax.plot(d, c, "o", color=colour, markersize=3.5, alpha=0.28, markeredgewidth=0, zorder=z)
        nd = nondominated(list(zip(d, c)))
        pts = sorted(set((d[i], c[i]) for i in nd))
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o" if alg == "NSGA-II" else "s", color=colour, markersize=9, markerfacecolor=PAPER,
                markeredgewidth=1.8, zorder=4, label=f"{alg}: its non-dominated chains ({len(nd)}; {len(pts)} distinct scores)")
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=colour, linewidth=0.9, zorder=3)
    cm = defaultdict(list)
    for x in held:
        if x["chain"].startswith("CMA-ES") and x["set"] == "training":
            cm[x["chain"]].append(x)
    hs = [x for x in held if x["chain"] == "hand-set" and x["set"] == "training"]
    hd, hc = np.mean([float(x["dice"]) for x in hs]), np.mean([float(x["cldice"]) for x in hs])
    ax.plot(hd, hc, "D", color=HAND, markersize=10, markeredgecolor=INK, markeredgewidth=0.8, zorder=5)
    ax.annotate("hand-set chain: Frangi's α = β = 0.5,\nthresholds 0.1–0.3; it grows into the balls", xy=(hd, hc), xytext=(0.65, 0.918), fontsize=11, color=INK,
                arrowprops=dict(arrowstyle="-", color=HAIRLINE, linewidth=0.8))
    for name, xs in cm.items():
        ax.plot(np.mean([float(x["dice"]) for x in xs]), np.mean([float(x["cldice"]) for x in xs]), "^", color=ALGORITHM["CMA-ES"], markersize=9,
                markerfacecolor=ALGORITHM["CMA-ES"], markeredgecolor=PAPER, zorder=6, label="CMA-ES's best of each seed (one score)" if name.endswith("1") else None)
    ax.set_xlim(0.55, 0.95)
    ax.set_ylim(0.90, 1.0)
    # the front's corner, enlarged
    box = (0.855, 0.936, 0.986, 0.9945)
    ax.add_patch(plt.Rectangle((box[0], box[2]), box[1] - box[0], box[3] - box[2], fill=False, edgecolor=MUTED, linewidth=0.8, zorder=7))
    ins = ax.inset_axes([0.1, 0.64, 0.4, 0.32])
    for alg, colour, mk in (("random", ALGORITHM["random"], "s"), ("NSGA-II", ALGORITHM["NSGA-II"], "o")):
        r = [x for x in ev if x["algorithm"] == alg]
        d = np.array([float(x["dice"]) for x in r])
        c = np.array([float(x["cldice"]) for x in r])
        ins.plot(d, c, "o", color=colour, markersize=3, alpha=0.3, markeredgewidth=0)
        nd = nondominated(list(zip(d, c)))
        pts = sorted(set((d[i], c[i]) for i in nd))
        ins.plot([q[0] for q in pts], [q[1] for q in pts], mk, color=colour, markersize=8, markerfacecolor=PAPER, markeredgewidth=1.6)
    for name, xs in cm.items():
        ins.plot(np.mean([float(x["dice"]) for x in xs]), np.mean([float(x["cldice"]) for x in xs]), "^", color=ALGORITHM["CMA-ES"], markersize=8, markeredgecolor=PAPER)
    ins.set_xlim(box[0], box[1])
    ins.set_ylim(box[2], box[3])
    ins.tick_params(labelsize=9)
    ins.set_facecolor(PAPER)
    for sp in ins.spines.values():
        sp.set_visible(True)
        sp.set_color(MUTED)
        sp.set_linewidth(0.8)
    ins.set_title("the corner, enlarged", fontsize=10.5, color=MUTED, loc="left", pad=4)
    ax.set_xlabel("Dice: overlap of the volumes (mean of two training phantoms)")
    ax.set_ylabel("clDice: overlap of the centrelines")
    ax.legend(loc="lower left", fontsize=10.5)
    hv = {x["algorithm"]: float(x["hypervolume"]) for x in rows("hypervolume.csv")}
    panel_title(ax, "On the training phantoms", f"800 evaluations each; hypervolume to (0, 0): NSGA-II {hv['NSGA-II']:.4f}, random {hv['random']:.4f}")

    ax = fig.add_subplot(gs[0, 1])
    chains = OrderedDict()
    for x in ft:
        chains.setdefault((x["algorithm"], x["index"]), []).append(x)
    for alg, colour, marker in (("NSGA-II", ALGORITHM["NSGA-II"], "o"), ("random", ALGORITHM["random"], "s")):
        tr, worst, single = [], [], []
        for (a, _), xs in chains.items():
            if a != alg:
                continue
            tr.append(float(xs[0]["train_dice"]))
            worst.append(min(float(y["dice"]) for y in xs))
            single.append(xs[0]["sigma_lo"] == xs[0]["sigma_hi"])
        tr, worst, single = np.array(tr), np.array(worst), np.array(single)
        # identical points stacked: their count sets the size
        cnt = defaultdict(int)
        for a, b in zip(tr, worst):
            cnt[(round(a, 4), round(b, 4))] += 1
        for (a, b), n in cnt.items():
            ax.plot(a, b, marker, color=colour, markersize=6 + 2.2 * np.sqrt(n), markerfacecolor=PAPER, markeredgewidth=1.6, alpha=0.95)
            if n >= 5:
                ax.text(a + 0.0012, b + 0.03, f"×{n}", fontsize=10.5, color=colour, ha="left", va="bottom")
        failed = int(np.sum(worst == 0))
        ax.plot([], [], marker, color=colour, markerfacecolor=PAPER, markeredgewidth=1.6, linestyle="none",
                label=f"{alg}: {len(tr)} chains, {failed} find nothing on a held-out phantom")
    ax.plot([0.85, 0.94], [0.85, 0.94], color=HAIRLINE, linewidth=0.9)
    ax.text(0.857, 0.80, "hairline: test = training; ring sizes: chains with the same scores", color=MUTED, fontsize=10, ha="left")
    ax.set_xlim(0.855, 0.94)
    ax.set_ylim(-0.04, 1.0)
    ax.set_xlabel("Dice on the training phantoms")
    ax.set_ylabel("worst Dice over the three held-out phantoms")
    ax.legend(loc="center left", fontsize=10.5)
    ax.annotate("one scale, 2.35 mm, and a high threshold\nnear 0.144: on phantom 4 nothing is kept", xy=(0.9257, 0.0), xytext=(0.858, 0.22), fontsize=11, color=INK,
                arrowprops=dict(arrowstyle="-", color=HAIRLINE, linewidth=0.8))
    panel_title(ax, "On phantoms not seen in tuning", "each non-dominated chain's training Dice against its worst held-out Dice")
    source(fig, f"heuresis multiobjective::Nsga2 (μ 32, SBX and polynomial mutation), RandomSearch, strategies::Cma; eikon phantoms, features::frangi, segmentation, fusion::dice, cldice.  {VERSION}")
    save(fig, "search-47-vessel-chain-front")


# ---- search-48: what four chains reconstruct on a held-out phantom --------------------------------------------

def projections_figure():
    pr = rows("projections.csv")
    layers = OrderedDict()
    for x in pr:
        layers.setdefault(x["layer"], []).append([float(v) for v in x["values"].split()])
    arr = {k: np.array(v) for k, v in layers.items()}
    held = rows("held_out.csv")
    score = {(x["chain"], x["phantom"]): x for x in held}
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "Four chains on a phantom none of them was tuned on",
          "Held-out phantom 4: a tree of tubes (radii 0.7 to 2 mm) among four balls, noise σ 0.5; projections along y of what each chain keeps")
    gs = fig.add_gridspec(2, 4, left=0.04, right=0.98, top=0.82, bottom=0.08, wspace=0.08, hspace=0.28)
    extent = [-17.75, 17.75, -17.75, 17.75]  # mm, the 72 voxels of 0.5 mm, centred

    def show(ax, img, cmap, vmin, vmax, name, sub):
        ax.imshow(img, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, extent=extent, interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.grid(False)
        ax.set_title(name, loc="left", fontsize=13, color=INK, pad=20)
        ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=10, color=MUTED, ha="left", va="bottom")

    grey = LinearSegmentedColormap.from_list("paper_ink", [PAPER, "#B8ACB4", INK])
    show(fig.add_subplot(gs[0, 0]), arr["slice"], grey, -0.5, 1.5, "The image", "the middle slice, y = 0.25 mm (intensities)")
    show(fig.add_subplot(gs[0, 1]), arr["image"], grey, 0.6, 2.4, "Its maximum along y", "the tree and the balls, through the noise")
    truth_voxels = int(float(score[("hand-set", "4")]["truth_voxels"]))
    show(fig.add_subplot(gs[0, 2]), arr["truth"], ListedColormap([PAPER, INK]), 0, 1, "The truth", f"the tubes only, {truth_voxels:,} voxels")
    # the key
    kx = fig.add_subplot(gs[0, 3])
    kx.axis("off")
    key = [("#6E5A6A", "kept, in the truth"), (ALGORITHM["PSO"], "kept, not in the truth"), ("#7FB5B0", "in the truth, missed")]
    for i, (c, t) in enumerate(key):
        kx.add_patch(plt.Rectangle((0.05, 0.78 - i * 0.14), 0.09, 0.09, color=c, transform=kx.transAxes))
        kx.text(0.18, 0.825 - i * 0.14, t, transform=kx.transAxes, fontsize=12, color=INK, va="center")
    kx.text(0.05, 0.28, "Each pixel: whether any voxel along y\nis kept, or is in the truth.\nDice and clDice are of the volumes.", transform=kx.transAxes,
            fontsize=11, color=MUTED, va="top", linespacing=1.5)
    names = [("hand-set", "hand-set chain"), ("NSGA-II best Dice", "NSGA-II: best Dice"), ("NSGA-II best mean", "NSGA-II: best mean"), ("CMA-ES seed 2", "CMA-ES: seed 2")]
    overlay = ListedColormap([PAPER, "#7FB5B0", ALGORITHM["PSO"], "#6E5A6A"])
    T = arr["truth"] > 0.5
    for j, (layer, label) in enumerate(names):
        P = arr[layer] > 0.5
        code = np.zeros(P.shape)
        code[T & ~P] = 1
        code[P & ~T] = 2
        code[P & T] = 3
        s = score[(layer, "4")]
        sub = f"Dice {float(s['dice']):.3f} · clDice {float(s['cldice']):.3f} · {int(float(s['voxels'])):,} voxels kept"
        if float(s["voxels"]) == 0:
            sub = "nothing kept: Dice 0, clDice 0"
        show(fig.add_subplot(gs[1, j]), code, overlay, -0.5, 3.5, label, sub)
    source(fig, f"eikon phantoms::tree (depth 3, Murray's γ = 3) and four balls, rasterised with partial volumes; the chains of results/held_out.csv.  {VERSION}")
    save(fig, "search-48-vessel-reconstructions")


# ---- search-49: one score, CMA-ES against random search, and the cliff the thresholds make -----------------------

def score_figure():
    conv = rows("convergence.csv")
    ev = rows("evaluations.csv")
    held = rows("held_out.csv")
    fig = plt.figure(figsize=(16, 10), dpi=100)
    title(fig, "One score for the chain: CMA-ES against random search, on a landscape of plateaus",
          "1 − (Dice + clDice)/2 on the two training phantoms: the best found so far, and how the score is distributed over the box")
    gs = fig.add_gridspec(1, 2, left=0.065, right=0.97, top=0.80, bottom=0.12, wspace=0.22)
    ax = fig.add_subplot(gs[0, 0])
    finals = defaultdict(list)
    for alg, colour, ls in (("CMA-ES", ALGORITHM["CMA-ES"], "-"), ("random", ALGORITHM["random"], (0, (2, 2))), ("NSGA-II", ALGORITHM["DE"], (0, (6, 2)))):
        seeds = sorted(set(int(x["seed"]) for x in conv if x["algorithm"] == alg))
        for sd in seeds:
            r = sorted([(float(x["evaluations"]), float(x["best"])) for x in conv if x["algorithm"] == alg and int(x["seed"]) == sd])
            end_ = 800 if alg == "NSGA-II" else 320
            xs = [e for e, _ in r] + [end_]
            ys = [v for _, v in r] + [r[-1][1]]
            finals[alg].append(ys[-1])
            ax.step(xs, ys, where="post", color=colour, linewidth=1.8, linestyle=ls, alpha=0.9,
                    label=(f"{alg}, {len(seeds)} seeds" if alg != "NSGA-II" else "NSGA-II's best mean, one run") if sd == seeds[0] else None)
    ax.set_yscale("log")
    ax.set_ylim(0.03, 1.05)
    ax.set_xlim(0, 820)
    ax.set_xlabel("evaluations (each two phantoms reconstructed and scored)")
    ax.set_ylabel("best score so far, 1 − (Dice + clDice)/2")
    ax.legend(loc="upper right", fontsize=11, handlelength=3)
    lines = ["final best", "CMA-ES   " + ", ".join(f"{v:.4f}" for v in finals["CMA-ES"]), "random   " + ", ".join(f"{v:.4f}" for v in finals["random"]),
             "NSGA-II  " + ", ".join(f"{v:.4f}" for v in finals["NSGA-II"])]
    ax.text(440, 0.30, "\n".join(lines), fontsize=11, color=INK, va="top", ha="left", linespacing=1.6)
    panel_title(ax, "Convergence", "CMA-ES with the box transformation, σ₀ 0.3, λ 10, 320 evaluations; random search uniform in the box")

    ax = fig.add_subplot(gs[0, 1])
    bins = np.logspace(np.log10(0.03), np.log10(1.0001), 46)
    groups = (("random", ("random", "random (score)"), ALGORITHM["random"]), ("CMA-ES", ("CMA-ES",), ALGORITHM["CMA-ES"]))
    for name, algs, colour in groups:
        sc = np.array([1 - 0.5 * (float(x["dice"]) + float(x["cldice"])) for x in ev if x["algorithm"] in algs])
        w = np.full(len(sc), 100.0 / len(sc))
        ax.hist(np.clip(sc, 0.03, 1.0), bins=bins, weights=w, histtype="stepfilled" if name == "random" else "step", color=colour,
                alpha=0.45 if name == "random" else 1.0, linewidth=1.8, label=f"{name}: {len(sc):,} evaluations")
        good = 100.0 * np.mean(sc < 0.05)
        ax.text(0.031, 52 if name == "random" else 46, f"{name}: {good:.1f}% of its evaluations below 0.05", color=colour if name == "CMA-ES" else INK, fontsize=11, ha="left")
    ax.set_xscale("log")
    ax.set_xlim(0.03, 1.08)
    ax.set_ylim(0, 60)
    hand = [x for x in held if x["chain"] == "hand-set" and x["set"] == "training"]
    hscore = 1 - 0.5 * (np.mean([float(x["dice"]) for x in hand]) + np.mean([float(x["cldice"]) for x in hand]))
    ax.annotate(f"the tree with the balls:\nthe hand-set chain, {hscore:.3f}", xy=(hscore, 14), xytext=(0.30, 33), fontsize=11, color=INK,
                arrowprops=dict(arrowstyle="-", color=HAIRLINE, linewidth=0.8))
    ax.annotate("nothing kept: 1", xy=(0.99, 6), xytext=(0.52, 24), fontsize=11, color=INK, arrowprops=dict(arrowstyle="-", color=HAIRLINE, linewidth=0.8))
    ax.annotate("the tree alone", xy=(0.04, 3), xytext=(0.045, 20), fontsize=11, color=INK, arrowprops=dict(arrowstyle="-", color=HAIRLINE, linewidth=0.8))
    ax.set_xlabel("score, 1 − (Dice + clDice)/2 (logarithmic)")
    ax.set_ylabel("share of the algorithm's evaluations (%)")
    ax.legend(loc="upper right", fontsize=11)
    panel_title(ax, "Where the evaluations fell", "random search samples the box uniformly; CMA-ES, three runs, straddles the plateaus")
    source(fig, f"heuresis strategies::Cma (CmaBounds::transform), search::RandomSearch, multiobjective::Nsga2; eikon features::frangi on cached Hessian eigenvalues at 13 scales.  {VERSION}")
    save(fig, "search-49-vessel-chain-score")


if __name__ == "__main__":
    front_figure()
    projections_figure()
    score_figure()
