---
title: Tsunami Travel Times
status: draft
implemented-by:
  - domains/oceanography/js/tsunami.js
  - applications/earth-atlas/plugins/tsunami.js
references: [Lamb1932, Dijkstra1959, Shokin1987, Wessel2009, Watada2014, KongRosenfeld1989, Atwater2005]
---

# Tsunami Travel Times

## 1. The Long Wave

A tsunami is a train of surface gravity waves set off by a sudden displacement of the sea floor or the sea. Its wavelengths, 100 to 500 km, are long beside the depth of the ocean, a few kilometres.

*Hypotheses.* Small-amplitude gravity waves on water of depth $h$, inviscid and incompressible, on a rigid bottom, with gravity $g$.

**Proposition 1.** For wavenumber $k$ the phase speed is $c = \sqrt{(g/k)\tanh(kh)}$; for $kh \to 0$, $c \to \sqrt{g h}$, the same for every wavelength (Lamb 1932, §170).

*Proof.* The dispersion relation is $\omega^2 = g k \tanh(kh)$ ([Wave energy propagation](wave-energy-propagation.md), eq. 1), and $\tanh x = x + O(x^3)$. ∎

Over 4000 m, $c = 198$ m/s, 713 km/h; over 100 m, 31 m/s. A wave of 200 km in 4000 m has $kh = 0.13$, and its speed is within 0.3% of $\sqrt{gh}$. The long wave is non-dispersive, and its front, the first arrival, travels at $\sqrt{gh}$.

## 2. The First Arrival as a Shortest Path

**Definition 1.** The *travel time* $T(x)$ from a source $x_0$ is the least time in which the front can reach $x$ through the sea: the infimum, over paths $\gamma$ from $x_0$ to $x$ that do not cross land, of $\int_\gamma \mathrm{d}s / c$, $s$ the length along the sphere.

$T$ is the viscosity solution of the eikonal equation $|\nabla T| = 1/c$ with $T(x_0) = 0$ (Fermat's principle). It is computed here on a graph.

**Definition 2 (the graph).** The nodes are the wet points of a longitude–latitude grid ($h > 0$). Each node is joined to the nodes of its row's *stencil* (§3) by a straight step whose cost is its great-circle length $\ell$ times the mean slowness of its ends, $\ell\,(1/c_i + 1/c_j)/2$. $T$ is the length of a shortest path from the source, found by Dijkstra's (1959) algorithm with a binary heap, in $O(n \log n)$ for $n$ nodes.

This is the method of the travel-time charts of Shokin et al. (1987) and of Wessel's (2009) TTT program, whose times he compared with observed arrivals in the Pacific and Indian Oceans.

## 3. The Stencil

A shortest path on a graph can only turn through the directions of its steps, so in a uniform sea it is longer than the straight path.

**Proposition 2.** In a uniform medium on a plane, let the stencil's directions leave a largest angular gap $\alpha$ between consecutive ones. A shortest grid path to any point is at most $1/\cos(\alpha/2)$ times the straight distance, and this bound is attained.

*Proof.* A point in the sector between consecutive step directions $e_1, e_2$ is reached by $a e_1 + b e_2$, $a, b \ge 0$, of length $a + b$. Its distance is $\sqrt{a^2 + b^2 + 2ab\cos\alpha'}$, $\alpha' \le \alpha$ the sector's angle, and the ratio $(a + b)/\sqrt{\cdot}$ is greatest at $a = b$, where it is $1/\cos(\alpha'/2)$. ∎

On square cells, the 16 directions $(0, \pm1), (\pm1, 0), (\pm1, \pm1), (\pm1, \pm2), (\pm2, \pm1)$ leave $\alpha = 26.6°$: 2.75%. On a longitude–latitude grid the cells' east–west side is $\cos\varphi$ times the north–south one, and the gaps open toward the poles. The solver therefore multiplies each row's column offsets by $1/\cos\varphi$ and rounds them, keeping the unit steps $(0, \pm1)$. Adding the steps $(\pm1, \pm3), (\pm3, \pm1), (\pm2, \pm3), (\pm3, \pm2)$ gives 32 directions. The worst excess, by Proposition 2:

| Latitude | 16, fixed | 16, scaled | 32, scaled |
|---|---|---|---|
| 0° | 2.75% | 2.75% | 1.31% |
| 30° | 3.53% | 3.53% | 1.71% |
| 45° | 4.93% | 3.41% | 1.46% |
| 50° | 5.72% | 4.22% | 2.08% |
| 60° | 8.24% | 2.75% | 1.31% |
| 75° | 17.05% | 2.92% | 1.40% |

The mean excess over all directions is about a third of the worst. The plug-in uses 32 neighbours by default.

## 4. Land and Corners

**Definition 3.** A step from cell $(0,0)$ to cell $(\Delta r, \Delta c)$ *passes through* the cells its straight segment, between the cells' centres, enters; and *through the corner* between two cells where it goes from a cell to its diagonal neighbour. A step is taken only if every cell it passes through is wet, and one of the two cells beside each corner it passes through.

**Proposition 3.** A coast of land cells that is a closed curve in the 8-neighbourhood sense (cells touching at their corners count as joined) is not crossed by any step.

*Proof.* Replace, at each corner a step passes through, the corner by the wet cell beside it. The cells the step passes through then form a 4-connected chain of wet cells from its start to its end. A closed 8-connected curve separates the plane into 4-connected components (the digital Jordan curve theorem; Kong and Rosenfeld 1989), so the chain lies within one of them. ∎

*Remark.* Testing only the cells at the points of a sampled segment is not enough: a step through a corner between two diagonal land cells passes untested, and the wave leaks through coasts of single cells. The checks include such a coast.

## 5. The Source

The source is a point. It is moved to the nearest wet grid point within 100 km, if it is not on one; then that point and those its stencil reaches are given their straight-line times from the true source, $\ell(x_0, x)\,(1/c_0 + 1/c_x)/2$, so that the grid's offset from the source adds no time.

## 6. Verification

| Test | Result |
|---|---|
| $\sqrt{g h}$ over 4000 m | 198.06 m/s |
| A uniform ocean of 4000 m on a ½° grid, sources at 0.2°, 47.3° and −62.1°, all points 1000–9000 km away: worst excess of the 32-neighbour stencil (mean) | 1.3%, 1.7%, 2.6% (0.34%, 0.30%, 0.28%) |
| The same, 16 neighbours | 2.7%, 3.1%, 4.3% (0.96%, 0.79%, 0.72%) |
| A wall of land from pole to pole with a gap at 30° N: the time behind it against the path through the gap | −1.3%: the gap is 2° wide, and the shortest path passes near its southern edge, not through its middle |
| A closed sea whose coast is a diamond of single cells touching at their corners | No point outside reached |
| A source on an island of 5 × 5 cells: the nearest sea; none if beyond reach | 164.38 km, as the great circle gives; none |
| Plug-in: a point 3000 km away over 4000 m (4.21 h) | $3 \times 10^{-8}$ |
| 22 May 1960, from the USGS epicentre of the Chilean earthquake (38.24° S, 73.05° W, moved 46 km to the sea) over ETOPO 2022 at 20′: Hilo, Sanriku coast of Japan | 14.9 h and 22.4 h; observed about 15 h and 22 h (Atwater et al. 2005) |
| Run time, ETOPO 2022 at 10′ (2160 × 1080, 1.53 million points of sea), from Panama: 16 and 32 neighbours, in Node | 2.7 s; 4.7 s |
| The plug-in in the page, 32 neighbours at 10′ (the 5′ relief averaged in blocks of 2 × 2) | 4.0 s |

*Caveat.*
- The times are those of linear long waves from a point. A great earthquake breaks 100–1000 km of fault, and coasts facing the rupture are reached sooner than from its epicentre.
- Observed first arrivals across an ocean are later than the long-wave times by up to a few per cent, as sea water is compressible, the solid Earth yields under the wave, and the wave's own mass perturbs gravity (Watada et al. 2014).
- Near the coast the times depend on the bathymetry's resolution: a cell's depth is a mean over 18 km at 10′, and shoaling within it is not resolved.
- A travel time says when a wave could first arrive, not whether one was made, nor its height. Warnings are the tsunami warning centres'.
