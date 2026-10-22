---
title: Zonal Statistics on the Sphere
status: draft
implemented-by:
  - domains/geography/js/zonal.js
  - data/products/regions.py
  - data/products/population_grid.py
references: [Schiavina2023, Pesaresi2024, GHSL2023, Foley1990, Moritz2000]
---

# Zonal Statistics on the Sphere

## 1. Fields and Regions as Functions on the Sphere

A question of human geography is often a question about a region: how many people live in a country, how warm a province will be, how many people in it the coming storm will reach. The data that answer it come as fields on grids, and the regions come as polygons. Zonal statistics join the two.

**Definition 1 (grid).** A regular grid of longitude and latitude is given by the centre $(\lambda_0, \varphi_0)$ of its first cell, its steps $\Delta\lambda, \Delta\varphi$ and its numbers of columns and rows $n_\lambda, n_\varphi$; rows run from north to south. Its cell $(r, c)$ is the rectangle of longitude and latitude
$$C_{rc} = [\lambda_c - \tfrac12\Delta\lambda,\ \lambda_c + \tfrac12\Delta\lambda] \times [\varphi_r - \tfrac12\Delta\varphi,\ \varphi_r + \tfrac12\Delta\varphi] \cap \{|\varphi| \le \tfrac\pi2\},$$
with $\lambda_c = \lambda_0 + c\,\Delta\lambda$ and $\varphi_r = \varphi_0 - r\,\Delta\varphi$. The grid *goes round the Earth* when $n_\lambda \Delta\lambda = 2\pi$; longitudes are then taken modulo $2\pi$.

**Definition 2 (field).** A field on a grid is a function $f$ on its cells, with values in $\mathbb R \cup \{\bot\}$, $\bot$ meaning no value. It is read as the function on the sphere $S^2_R$ of radius $R$ that is constant on each cell. Its *domain* $D_f$ is the union of the cells with a value.

**Definition 3 (regions).** A regions raster is a function $\rho$ from the cells of a grid to $\{0, 1, \dots, K\}$; it partitions the sphere into the sets $S_k = \bigcup \{C : \rho(C) = k\}$, region $k$ for $k \ge 1$ and none for $k = 0$.

A regions raster is made from polygons by the rule of the centre: a cell belongs to the region whose polygon contains its centre, by the even–odd rule, scanline by scanline (Foley et al. 1990, §3.6). Two polygons that share an edge never both claim a cell, and a boundary is placed to within half a cell of the raster. A region smaller than a cell may hold no cell: at 0.1° (about 11 km), 14 of Natural Earth's 242 countries hold none (the Vatican, Monaco, Macao, the Maldives, Nauru, Tuvalu), and 457 of its 4596 first-order divisions (63 of Malta's local councils, 62 of Slovenia's municipalities, Liechtenstein's communes). Their people are counted with their neighbours, or, for islands far from any other region, not at all (§5).

## 2. The Area of a Cell

**Proposition 1.** The area of the rectangle $[\lambda_1, \lambda_2] \times [\varphi_1, \varphi_2]$ on $S^2_R$ is
$$|A| = R^2 (\lambda_2 - \lambda_1)(\sin\varphi_2 - \sin\varphi_1). \tag{1}$$

*Proof.* The area element of the sphere is $\mathrm dA = R^2 \cos\varphi\, \mathrm d\varphi\, \mathrm d\lambda$, and $\int_{\varphi_1}^{\varphi_2} \cos\varphi\, \mathrm d\varphi = \sin\varphi_2 - \sin\varphi_1$. $\square$

**Corollary 1.** The band between the parallels $\varphi_1 < \varphi_2$ has area $2\pi R^2 (\sin\varphi_2 - \sin\varphi_1)$, and the cells of a grid round the Earth sum to $4\pi R^2$.

The coordinate $y = \sin\varphi$ is the ordinate of Lambert's cylindrical equal-area projection: in the plane $(\lambda, \sin\varphi)$ areas are those of the sphere divided by $R^2$, which is why `featureAreaFraction` of `raster.js` measures polygons there (see the note on equal-area maps).

*Remark.* The data of this domain are on the WGS 84 ellipsoid, and the sphere here has its mean radius $R = 6371.0088$ km (Moritz 2000). Between the same parallels of geodetic latitude, a cell of the ellipsoid is 0.45 % smaller than (1) at the equator, 0.2 % larger at 45° and 0.7 % larger at 70°: the error of an area, and of a density per area, computed on the sphere.

## 3. The Overlap of Two Grids

A field and a regions raster rarely share a grid: the population is at 0.25°, the regions at 0.1°, a forecast at 1°, the relief at 5′.

**Proposition 2.** The intersection of two rectangles of longitude and latitude is a rectangle of longitude and latitude, and its area factorises:
$$|C \cap C'| = R^2\, \ell_\lambda(C, C')\, \ell_{\sin\varphi}(C, C'), \tag{2}$$
where $\ell_\lambda$ is the length of the intersection of their intervals of longitude (modulo $2\pi$) and $\ell_{\sin\varphi}$ that of their intervals of $\sin\varphi$.

*Proof.* A rectangle is a product of an interval of longitude and an interval of latitude, the intersection of products is the product of intersections, and $\sin$ is increasing on $[-\tfrac\pi2, \tfrac\pi2]$, so that the interval of latitude maps to an interval of $\sin\varphi$ whose length gives the area by (1). $\square$

The overlaps of every cell of one grid with every cell of the other are therefore the products of two lists: for each row of the field, the rows of the regions it meets and $\ell_{\sin\varphi}$; for each column, the columns it meets and $\ell_\lambda$ (`overlapTable`). The lists are computed once, in $O(n_\varphi + n_\lambda)$ entries of each grid, and a cell of the field visits $(m_\varphi + 1)(m_\lambda + 1)$ cells of the regions at most, $m$ being the ratio of the steps: 9 for the population on regions at 0.1°, 121 for a forecast at 1°, 4 for the relief at 5′.

## 4. The Statistics

**Definition 4.** For a field $f$ and regions $\rho$, with $w_k(C) = |C \cap S_k|$ for a cell $C$ of the field:

- the *area* of region $k$ is $A_k = \sum_{C \subseteq D_f} w_k(C)$;
- its *mean* is $\bar f_k = A_k^{-1} \sum_C w_k(C)\, f(C)$, the area-weighted mean;
- its *least* and *greatest* values are $\min$ and $\max$ of $f(C)$ over the cells with $w_k(C) > 0$;
- its *cells* are the number of those cells.

**Proposition 3.** $\bar f_k = A_k^{-1} \int_{S_k \cap D_f} f\, \mathrm dA$ exactly; its least and greatest values are the essential infimum and supremum of $f$ on $S_k \cap D_f$. A constant field has its value as the mean, least and greatest value of every region.

*Proof.* $f$ is constant on each cell, so $\int_{S_k \cap D_f} f\, \mathrm dA = \sum_C f(C)\, |C \cap S_k|$. The values taken on a set of positive area are those of the cells that meet it with positive area. $\square$

The least and greatest values are those of the field as it is given: a forecast's cell of 1° that touches a country by a sliver contributes its value, as it should, since the forecast holds that value over the whole cell.

**Counts.** A field of counts per cell, such as people, is extensive: it is not averaged but divided. The count of a cell is shared among the regions it meets in proportion to the areas it shares with them,
$$\Sigma_k = \sum_C f(C)\, \frac{w_k(C)}{W(C)},\qquad W(C) = \sum_{j \ge 1} w_j(C), \tag{3}$$
over the regions only (`share: 'regions'`): a cell half in a country and half in the sea gives the country all its people, who live on land. With `share: 'area'`, $W(C) = |C|$, and the part outside every region goes to region 0.

**Proposition 4 (conservation).** $\sum_{k \ge 0} \Sigma_k = \sum_C f(C)$: every count is counted once.

*Proof.* For a cell meeting some region, $\sum_{k \ge 1} w_k(C) / W(C) = 1$; a cell meeting none gives its count to $\Sigma_0$ (or, §5, to one region). With `share: 'area'`, $\sum_{k \ge 0} w_k(C) = |C|$. $\square$

## 5. Coasts Finer than the Regions

A coast drawn by the cells of a raster is a staircase: a cell of the regions is land when its centre is. A cell of the field may then hold people on a strip of coast narrower than the raster's cells, and meet no cell of any region. Kushiro, on the coast of Hokkaidō, is such a place: its cell of 0.25° holds 74 000 people, and every cell of 0.1° within it has its centre in the sea. Worldwide, 20.9 million people (0.26 %) of GHS-POP's 0.25° grid live in cells that meet no country at 0.1°.

With a *reach* of $q$ cells, the count of a cell that meets no region goes to the region of the nearest cell of the regions, within $q$ cells of those it meets, by distance in the tangent plane at its centre. With $q = 3$ (0.3°, about 33 km), the people outside every country fall from 20.9 to 1.96 million (0.024 %), those of ships, oil platforms and islands beyond 33 km of any region's cell. The reach can misplace a small island near another country's coast to that country; it is used for counts only.

## 6. People per Region

GHS-POP R2023A (Schiavina et al. 2023) gives the residential population of 2025 on a 30″ grid, 43 202 × 21 384 cells in 256 × 256 tiles of float64 values compressed by LZW. The product `fields-ghsl-population` sums them onto 0.25° (`population_grid.py`): each 30″ cell goes whole to the cell containing its centre, which places a boundary within 15″ and keeps the total exactly, 8 191 988 536 people, the release's own total for 2025 (GHSL Data Package 2023, Table 20). Stored as 16-bit integers with a quantum of 254 people, the total becomes 8 191 211 825 (0.01 % fewer, the cells with fewer than half a quantum being read as empty).

By the regions of Natural Earth at 0.1° with a reach of 3, Japan holds 121.92 million (its 47 prefectures 121.94 million), against 126 146 099 at the census of 1 October 2020 (Statistics Bureau of Japan): 3.3 % fewer, of which some 2 % is the decline from 2020 to 2025 in the UN's estimates (the World Bank's 123.98 million for 2024), and the rest the coasts and small islands of §5.

## 7. Verification

| Check | Result |
|---|---|
| A constant field: mean, least and greatest value of every region of a Voronoi partition, fields finer (0.2°) and coarser (0.75°, shifted) than the regions (0.5°) | exact to $5 \times 10^{-12}$ |
| The area of the bands 10°–20.5° N and 30°–5° S against Corollary 1 | $2 \times 10^{-13}$ |
| The regions and region 0 against $4\pi R^2$ | $4 \times 10^{-14}$ |
| Counts kept (Proposition 4), each way of sharing, with and without reach | $4 \times 10^{-16}$ |
| GHS-POP by country and by province: every person counted once | $4 \times 10^{-16}$ |
| Japan against its census of 2020 | −3.3 % (limit 5 %) |
| The countries' raster of `regions.py` against `fillRings` of `raster.js` on the same polygons | $1.5 \times 10^{-7}$ of the cells differ (float32 vertices) |
| The raster's area of each country over 300 000 km² against its polygon's | within 0.63 % (Vietnam) |

The checks are in `domains/geography/tests/run-tests.mjs` and `data/tests/check_geography.py`.

## References

- European Commission, Joint Research Centre (2023), *GHSL Data Package 2023*, Publications Office of the European Union, JRC133256, doi:10.2760/098587.
- Foley, J. D., van Dam, A., Feiner, S. K., and Hughes, J. F. (1990), *Computer Graphics: Principles and Practice*, 2nd ed., Addison-Wesley, §3.6.
- Moritz, H. (2000), Geodetic Reference System 1980, *Journal of Geodesy* 74, 128–133.
- Pesaresi, M., Schiavina, M., Politis, P., Freire, S., Krasnodębska, K., Uhl, J. H., Carioli, A., et al. (2024), Advances on the Global Human Settlement Layer by joint assessment of Earth Observation and population survey data, *International Journal of Digital Earth* 17(1), 2390454, doi:10.1080/17538947.2024.2390454.
- Schiavina, M., Freire, S., Carioli, A., and MacManus, K. (2023), *GHS-POP R2023A: GHS population grid multitemporal (1975–2030)*, European Commission, Joint Research Centre, doi:10.2905/2FF68A52-5B5B-4A22-8F40-C41DA8332CFE.
- Statistics Bureau of Japan (2021), *2020 Population Census*: final count of the population, 126 146 099.
