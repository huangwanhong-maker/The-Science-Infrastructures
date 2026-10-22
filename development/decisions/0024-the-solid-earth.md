# 0024. The Solid Earth: Plates, GPS, Earthquake Mechanisms, the Crust and Sections

- **Status:** Accepted, provisionally: built at the user's request; open to revision on review. Making geophysics a field of `core/domains/` and `applications/` is proposed (§2.7) and awaits the user's decision.
- **Date:** 2026-10-22
- **Decided by:** the user asked which ways of seeing the motion of the solid Earth were feasible from global open data, and then asked for all but the real-time streams (NTRIP, SeedLink, InSAR, LiDAR). Claude chose the sources and the design below.

## 1. Context

The Earth Atlas showed earthquakes as dots of a month (USGS) and disaster alerts (GDACS), but nothing of why the ground moves or what lies beneath it. The open sources, checked on 2026-10-22, are these:

- **Plate kinematics.** NNR-MORVEL56 (Argus, Gordon and DeMets 2011, 56 plates; Table 1, transcribed) and ITRF2020's plate motion model (Altamimi et al. 2023, 13 plates; as transcribed in PROJ). The plates' outlines are Bird's PB2002 (2003), as GeoJSON from fraxen/tectonicplates (ODC-By). EarthScope's Plate Motion Calculator gives velocities for checking.
- **GPS.** The Nevada Geodetic Laboratory publishes the velocities of about 22,000 stations (MIDAS; Blewitt et al. 2016) and each station's daily positions, in IGS20, without CORS headers and without a licence, asking for citation.
- **Earthquake mechanisms.** The Global CMT catalogue: about 70,000 moment tensors since 1976, as NDK text; the request is citation only.
- **The crust.** CRUST1.0 (Laske et al. 2013), on 1-degree cells, 1.2 MB, citation only.
- **Slab geometry.** USGS Slab2 (public domain) is behind a Cloudflare challenge on ScienceBase and cannot be fetched by a script.

## 2. Decision

1. **Plates and their motion.**
   - The plates are PB2002's 52 outlines, with the 5,819 steps of their boundaries. Each step carries its class, its left and right plate, the polarity of subduction, and Bird's own relative velocity (data library: `pb2002`, `plate_model.py`, product kind `plates`).
   - Their motion is NNR-MORVEL56, kept in the core (`src/core/plates.js`) as a model is, with its citation. The velocity is v = ω × r on WGS 84.
   - The motion is offered as a vector field computed in the page, so the existing flow and surface layers draw it. Their pace and the conversion of mm/yr to m/s are properties of the quantity.
   - The place card names a place's plate, its motion and its pole. It gives the nearest boundary and the relative motion across it, resolved by Bird's conventions into opening and right-lateral sliding, and says which plate sinks beneath which.
2. **GPS stations.**
   - **Velocities.** The MIDAS velocities (`ngl-midas`, product kind `stations`) are drawn as arrows of constant length on screen, with dots coloured by vertical rate.
   - **Frames.** They are shown either in IGS20 or relative to the station's plate. The plate's motion subtracted is ITRF2020's model, which is in the stations' own frame, for its 13 plates, and NNR-MORVEL56 for the others. With NNR-MORVEL56 everywhere, stations in stable interiors kept 1.7 to 3.2 mm/yr of frame difference; with ITRF2020's model they keep 0.3 to 0.6 mm/yr.
   - **Daily positions.** A station's record (`ngl-station`, product kind `series`) is requested on demand (ADR 0023). It is drawn as east, north and up against time, with the MIDAS line and the earthquakes near enough to have moved it. Text parameters may now be held to a pattern, since a station's name is put into a URL.
3. **Earthquake mechanisms.**
   - The Global CMT catalogue (`gcmt`) becomes an events product of kind `mechanism`. It has its own layer, not the earthquakes' layer.
   - The core (`src/core/moment.js`) computes the axes by Jacobi rotations, the nodal planes from T and P (Aki and Richards 2002), Frohlich's classes, the radiation pattern, and the equal-area projection.
   - Beach balls are drawn by a shader that evaluates xᵀMx for each pixel's direction, coloured by class and turned to the map's north. The layer shows either the past year or every earthquake of magnitude 6.5 and more since 1976.
   - A chosen earthquake gets its faulting, planes and axes, and a stereonet, in the place card.
4. **The crust.**
   - CRUST1.0 (`crust1`, `crust_model.py`) becomes a `crust` product, holding the nine layers' tops, P and S speeds and densities as scaled 16-bit integers.
   - It also becomes a field product of the Moho's depth and of the thicknesses of the crust and of its sediments, which are offered as maps.
   - The place card draws the column beneath a place to scale, from the cell containing it, as the model's authors advise.
5. **Sections through the Earth.**
   - From the place card, a section runs across the nearest plate boundary within 900 km. At a subduction zone it heads toward the overriding plate, from 400 km before the boundary to 1,200 km beyond it, so that the sinking plate dips to the right.
   - It shows the relief, the crust to the Moho, the boundaries crossed, the 410 and 660 km discontinuities, and the Global CMT earthquakes within 50, 100 or 200 km of the line.
   - The section can be turned by 15°, is drawn on the map, and reads under the pointer. Its geometry is in `src/core/section.js`.
6. **A latent bug fixed.** The stage updated the camera's matrices only at the next render, while names, arrows and beach balls are placed on the screen earlier in the frame. They used the previous view until the camera moved again. The matrices are now updated when the camera moves.
7. **Geophysics, proposed as a field.** CLAUDE.md §1 names geophysics as a future field.
   - Its first theory note is `theory/geophysics/earthquake-mechanisms.md`. Plate kinematics, which belongs to geodesy, is `theory/geodesy/plate-kinematics.md`.
   - The code's comments name `core/domains/geophysics` as the stand-in's place.
   - Adding geophysics to the fields of `core/domains/` and `applications/` (ADR 0010) is a change of structure, and awaits the user's agreement. Until then the notes stay where they are and the code stays in the Earth Atlas.

## 3. Alternatives Considered

- **Real-time streams** (GNSS by NTRIP, seismograms by SeedLink, InSAR, LiDAR). NTRIP carries raw observations that would need a positioning engine and accounts. SeedLink needs a streaming server. InSAR needs accounts and very large files, and LiDAR local surveys. The user asked to leave these out.
- **One plate model for everything.** NNR-MORVEL56 covers all 56 plates and is geological, a 3-million-year average. ITRF2020's model covers 13 plates and is geodetic, from the present motion in GNSS's own frame. The first describes the plates and their boundaries; the second is what a GPS velocity should be compared with. Each is used where it fits, and the atlas says which.
- **Slab2 for sections.** It is blocked to scripts. The earthquakes themselves outline the slabs, and Slab2 can be added if fetched by hand.
- **USGS ComCat for sections.** It is denser, and its hypocentres are often better than centroids. It could be requested on demand by region (ADR 0023). Global CMT's centroids, already at hand, are enough to show the Wadati–Benioff zones, and come with mechanisms.

## 4. Consequences

- **The motion of the solid Earth, its causes and its interior** can be explored from any place: the plate and its pole, the ground's motion at 20,000 stations, the faulting of 70,000 earthquakes, the crust beneath, and a section through the boundary nearby.
- **Verification** (150 checks), among them:
  - plate velocities against the EarthScope calculator, to 0.005 mm/yr;
  - the 52 plates tiling the sphere;
  - Bird's boundary classes: 99.8 % of ridge steps open, 99.6 % of trench steps close, and 99.4 % of transforms slip in the sense he gives;
  - GPS stations in four stable plate interiors within 0.6 mm/yr of their plate (median);
  - a 20-year record against its MIDAS velocity, to 0.3 mm/yr;
  - the nodal planes of all 70,451 mechanisms against the catalogue's, 99 % within 0.76°;
  - constructed double couples;
  - CRUST1.0 beneath Tibet, the Pacific and the Canadian Shield;
  - section geometry, exact to 10⁻¹¹ km.
- **Data:**
  - the plates product, 0.3 MB;
  - the stations product, 1.5 MB;
  - the mechanisms product, 10 MB, read only when the layer or a section needs it;
  - the crust products, 5 MB.
- *Caveat.*
  - PB2002 has no outlines for MORVEL's Capricorn, Lwandle, Macquarie and Sur; there the containing plate's motion is used.
  - NNR-MORVEL56 is a geological average and differs from present GPS by a millimetre or two a year.
  - Global CMT's centroids are less precise in depth than relocated hypocentres.
  - CRUST1.0 averages 1-degree cells.
  - None of the providers but Bird's converter states a licence. The data are fetched into the cache and not redistributed.
