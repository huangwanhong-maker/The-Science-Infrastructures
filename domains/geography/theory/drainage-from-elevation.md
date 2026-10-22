---
title: Drainage from Elevation
status: draft
implemented-by:
  - domains/geography/js/hydrology.js
  - applications/earth-atlas/plugins/terrain.js
references: [Barnes2014, OCallaghan1984, Zevenbergen1987, Kahn1962, Lindsay2008, WRI2003]
---

# Drainage from Elevation

## 1. Water on a Grid of Cells

Water on land runs downhill until it reaches the sea. Given only the height of the surface, the rivers, their basins and the divides between them follow, once three things are fixed: where water leaves the land, what it does in a hollow, and which way it goes from a cell.

**Definition 1.** On a grid (as in the note on zonal statistics), the *neighbours* of a cell are the eight cells that share an edge or a corner with it, across the antimeridian on a grid round the Earth and never across a pole. A set of *outlets* $O$ is given. A *flow* is a map $d$ from the cells not in $O$ to cells, each to one of its neighbours. It *drains* when every chain $c, d(c), d^2(c), \dots$ reaches $O$.

**Proposition 1.** A flow drains if and only if its graph, the cells with an edge from each $c$ to $d(c)$, has no cycle; it is then a forest whose roots are the outlets.

*Proof.* Every cell outside $O$ has one edge out and the cells are finite, so a chain either reaches $O$, where it stops, or repeats a cell, which is a cycle. A graph in which every vertex but the roots has one edge out and there is no cycle is a forest of trees directed to the roots. $\square$

The tree of an outlet is its drainage basin; the subtree above a cell is that cell's *watershed*. The rest of this note makes a flow that drains, from heights.

## 2. The Outlets Are the Seas

ETOPO 2022 gives the height of the land and the depth of the sea floor as one surface, so that the cells at or below sea level, $z \le 0$, are the sea. They are not all sea: ETOPO also carries the depths of large lakes, whose floors lie below sea level (Superior to −137 m, Baikal to −1177 m on the 10′ grid), and the floors of dry depressions (Qattara, the Faiyum, the Dead Sea). Taking these as outlets would end the St Lawrence in Lake Superior and the Nile in the Faiyum.

**Definition 2.** The outlets are the cells with $z \le 0$ whose *body*, the component of such cells under the relation of neighbours, has an area of at least $10^5$ km², or reaches the edge of a grid that does not go round the Earth; and the cells without a value, and that grid's edges.

On the 10′ grid the bodies at or below sea level are the ocean (360.1 million km²), the Caspian with its shores (685 000 km²) and the Black Sea with the Sea of Azov (464 000 km², cut from the ocean where the Bosporus is narrower than a cell); then Lake Superior's deep part (26 000 km²), Qattara (20 000), Baikal (16 000), Lake Eyre (15 000) and some 540 smaller ones. The threshold falls in the gap between the third and the fourth. The Caspian, which has no outflow, is an outlet, as it is in nature.

## 3. Hollows Are Filled

A cell all of whose neighbours are higher is a pit: water reaching it stands, and in a grid of heights most pits are errors of measurement or of sampling. A hollow fills to the level at which it spills.

**Definition 3.** The *spill level* of a cell $c$ is
$$s(c) = \min_{\gamma : c \to O}\ \max_{x \in \gamma \setminus O} z(x),$$
the least, over the paths of neighbours from $c$ to an outlet, of the greatest height on the path. The *filled* surface is $\max(z, s)$.

**Proposition 2 (Priority-Flood; Barnes et al. 2014, Algorithm 1).** Flood the grid from the outlets inward, always taking next the cell of least height among those reached (a priority queue), and raise each newly reached cell to the height of the cell it was reached from if it lies lower. The heights so given are the filled surface.

*Proof (sketch).* The cells are taken in order of their filled height, as in Dijkstra's algorithm with the path cost $\max$ in place of $+$: when a cell is taken, no unexplored path to an outlet can reach it lower, since every cell left in the queue is at least as high. Barnes et al. (2014) give the proof in full. $\square$

A filled hollow is flat, and on a flat no neighbour is lower: a flow is not yet defined there. Priority-Flood+ε (Algorithm 3 of the same paper) raises a cell reached from a lower or equal one to the next floating-point number above it, $\mathrm{nextUp}(z)$, and takes it next from a plain queue (the *pit*), so that the filled surface rises along the order of the flood by the least possible steps. Here, when the next cell of the priority queue is no higher than the next of the pit, it goes first, so that the cells are taken in order of their filled height.

**Proposition 3.** After Priority-Flood+ε every cell outside $O$ has a neighbour strictly lower than itself.

*Proof.* A cell other than an outlet was reached from a cell $p$ already taken. If it went to the pit its height is $\mathrm{nextUp}(h(p)) > h(p)$; otherwise $h(c) = z(c) > \mathrm{nextUp}(h(p)) > h(p)$. $\square$

In float64 a step of $\mathrm{nextUp}$ at 1000 m is $1.1 \times 10^{-13}$ m; the 10′ grid raises 100 000 cells in all, so that no filled flat departs from a level by more than $10^{-7}$ m.

## 4. Steepest Descent on the Sphere

**Definition 4 (D8; O'Callaghan and Mark 1984).** From a cell $c$ outside $O$, water goes to the neighbour $n$ of greatest descent $(h(c) - h(n)) / \delta(c, n)$, where $h$ is the filled surface and $\delta$ the distance between the centres of the cells.

On a grid of longitude and latitude the distance depends on the direction and the latitude: a step east is $R\cos\varphi\,\Delta\lambda$ and a step north $R\,\Delta\varphi$, half as long at 60°. With the steps of a plane, a flow would turn east and west too readily at high latitudes. Here $\delta$ is the great-circle distance (by the haversine formula, as `centralAngle` of geodesy), computed once for the eight directions of each row.

**Corollary 2.** On the filled surface the flow of Definition 4 drains.

*Proof.* By Proposition 3 the greatest descent is positive, so $h$ decreases strictly along every chain; a chain cannot return to a cell, and by Proposition 1 the flow drains. $\square$

## 5. Upstream Area

**Definition 5.** The *upstream area* of a cell is the area of its watershed, the cell included: $U(c) = \sum_{x \to^* c} |x|$, with $|x|$ the area of cell $x$ on the sphere (zonal statistics, Proposition 1).

$U$ is computed by passing each cell's area to its receiver in a topological order of the forest (Kahn 1962): a cell is passed on once every cell draining into it has been, which needs no sorting and takes $O(n)$ steps. The watershed of a cell is found by the reverse search, from the cell to every neighbour whose receiver it is.

**Proposition 4.** The area of the watershed of $c$, found by search, equals $U(c)$; at the outlet of a basin closed but for that outlet it is the area of the basin.

*Proof.* Both are the sum of $|x|$ over the subtree above $c$ (Proposition 1). $\square$

The *main stem* above a cell climbs to the neighbour draining into it with the greatest upstream area, until a cell into which none drains; the path *downstream* follows the flow to the sea. A place chosen near a river is taken onto it by moving it to the cell of greatest upstream area within a cell or two, a pour point snapped to the stream (Lindsay et al. 2008).

## 6. Slope

**Definition 6.** The slope of the surface is $|\nabla z| = (G^2 + H^2)^{1/2}$ with centred differences in the metric of the sphere,
$$G = \frac{z_E - z_W}{2 R \cos\varphi\, \Delta\lambda}, \qquad H = \frac{z_N - z_S}{2 R\, \Delta\varphi}, \tag{1}$$
the slope of Zevenbergen and Thorne's (1987) quadratic surface, one-sided at the first and last rows.

**Proposition 5.** (1) is exact for $z = a R \varphi$ (slope $a$) and $z = b \lambda$ (slope $b / (R\cos\varphi)$).

*Proof.* A centred difference is exact for a function linear in its variable. $\square$

*Caveat.* On a 10′ grid (18 km) the slope is that of the surface smoothed over 37 km: in mountains several times less than the slopes of the ground. ETOPO's 10′ product takes every fortieth point of its 15″ grid rather than averaging, so it is also aliased.

## 7. Limits

- **Resolution.** On 10′ cells a river is a chain of cells and a valley narrower than a cell may be missed: a river may leave its basin through a pass the cells cannot see. On the flat divide between the Madeira and the Tapajós in Mato Grosso, the 10′ grid sends part of the Tapajós into the Madeira, and the Amazon at Óbidos receives 5.38 million km² against about 4.7 (+15 %), though the basin as a whole is right (6.20 against 6.2).
- **Basins without an outlet.** Filling makes every basin drain to the sea: Lake Chad's basin and the hollows of the Sahara spill into the Niger, which then drains 6.6 million km² (2.1 in nature); Lake Eyre spills into Spencer Gulf. Basins that end inland are a separate question; Barnes et al. (2014) give a variant of Priority-Flood that labels basins, from which such outlets could be chosen.
- **The Caspian** is an outlet, with the Volga's basin and the closed basins of Central Asia filled into it (6.6 million km² on the 10′ grid).

## 8. Verification

| Check | Result |
|---|---|
| A plane rising east (north) drains west (south) from every cell | exact |
| The path downstream on that plane: 25 steps along the parallel, by great circles | $2 \times 10^{-16}$ |
| Slope of $z = aR\varphi$ and $z = b\lambda$ | $2 \times 10^{-8}$ |
| A pit in a rim of 200 m with a notch at 150 m: filled to 150 m, every cell of it and no other raised | $3 \times 10^{-13}$ m |
| A closed basin through one outlet: upstream area at the outlet against $R^2 \Delta\lambda (\sin\varphi_2 - \sin\varphi_1)$ | $5 \times 10^{-14}$ |
| A watershed's area against the upstream area at its outlet | exact |
| The watershed of a cone's apex: the apex alone, outlined by four corners | exact |
| An L of three cells: an outline of six corners enclosing their area | $5 \times 10^{-13}$ |
| ETOPO 10′, 2.33 million cells: every land cell drains | 0 without a way down; filled in 0.55 s, flow in 0.17 s (Node) |
| Basins on ETOPO 10′ against their published areas (WRI 2003; Atlas of Canada; Guerrero et al. 1997): Amazon 6.20 (6.2), Río de la Plata 3.34 (3.17), Mackenzie 1.72 (1.81), St Lawrence with the Great Lakes 1.32 (1.34), Murray–Darling 1.04 (1.05) million km² | within 5.2 % (limit 10 %) |

The checks are in `domains/geography/tests/run-tests.mjs` and, for the plug-in, `applications/earth-atlas/tests/plugins-geography.mjs`.

## References

- Barnes, R., Lehman, C., and Mulla, D. (2014), Priority-flood: an optimal depression-filling and watershed-labeling algorithm for digital elevation models, *Computers & Geosciences* 62, 117–127.
- Guerrero, R. A., Acha, E. M., Framiñan, M. B., and Lasta, C. A. (1997), Physical oceanography of the Río de la Plata Estuary, Argentina, *Continental Shelf Research* 17(7), 727–742.
- Kahn, A. B. (1962), Topological sorting of large networks, *Communications of the ACM* 5(11), 558–562.
- Lindsay, J. B., Rothwell, J. J., and Davies, H. (2008), Mapping outlet points used for watershed delineation onto DEM-derived stream networks, *Hydrological Processes* 22, 1701–1709.
- Natural Resources Canada, *The Atlas of Canada*: Rivers (drainage areas of the Mackenzie and the St Lawrence).
- NOAA National Centers for Environmental Information (2022), *ETOPO 2022 15 Arc-Second Global Relief Model*, doi:10.25921/fd45-gt74.
- O'Callaghan, J. F., and Mark, D. M. (1984), The extraction of drainage networks from digital elevation data, *Computer Vision, Graphics, and Image Processing* 28, 323–344.
- World Resources Institute (2003), *Watersheds of the World*, Water Resources eAtlas.
- Zevenbergen, L. W., and Thorne, C. R. (1987), Quantitative analysis of land surface topography, *Earth Surface Processes and Landforms* 12, 47–56.
