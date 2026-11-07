# 0078. The Whole Earth in the Urban Science Workbench: One 3D GIS, Any City, the World's Development Data

- **Status:** Accepted, provisionally (2026-11-07). Amends ADR 0060 (the workbench's map), and the workbench's colour in CLAUDE.md §7. Amended by [ADR 0079](0079-data-in-space-and-time.md): works are named again beside the lenses.
- **Date:** 2026-11-06 to 2026-11-07
- **Decided by:**
  - **The request.** The user asked to "build a 3D GIS view, largely revise the UI", reusing the Earth Atlas and Cosmopoiesis, to integrate "such tools and data that we can integrate" from a list of the World Bank's and others' tools and data "and visualization on the 3D GIS", in pink-purple, because the page was "too normal, and unclear". Then: "we doesn't want just do several cities... If impossible to get all data, then ignore this requirement".
  - **The design.** Claude designed it, found that every city could be had, and built it.
  - **To be confirmed.** The user has yet to see it. Nothing is committed or pushed until they say.

## 1. Context

The workbench was a three-dimensional map of two cities on a flat page, with plug-ins in a drawer. The user wanted one GIS, the question put first, and the world's development data on it.

Three findings shaped it.
- **Every city can be had.** Overture Maps publishes the world's buildings, ways, water and land use as PMTiles archives, read tile by tile by range requests. Kathmandu's central tile has 19,342 buildings; it comes in about 2 s. GlobalBuildingAtlas's service forbids bulk use, and Overpass was at capacity.
- **Most of those buildings have no height.** In Freetown almost none of 99,000 do, and in Kathmandu's central tile about 13 % have a height or floors. GHS-BUILT-H gives the average of the buildings on a 100 m cell in 2018: an estimate that must look like one.
- **Geometry cannot be built on the page's thread.** A dense tile needs tens of thousands of footprints triangulated, and picking through them cannot go triangle by triangle.

## 2. Decision

### 2.1 One View from the Earth to a Doorstep

The Earth (`visualization/js/gis/earth.js`) is a globe of the shared `globe` layer on a canvas of its own, in double precision relative to the eye. A city scene rides above it on a transparent canvas, driven by one set of gestures. Plug-ins draw in the local frame as before and know nothing of it.

The Earth's imagery is graded into the page's palette, and a sky follows the real Sun. A shader may level the ground under a city's model and feather it back into the terrain. All of this is optional per frame and off by default, so Cosmopoiesis is unchanged.

### 2.2 Any Place on Earth

Coming down anywhere opens that place: its name from the gazetteer, its clock by longitude (and said to be so), and its tiles.

- **Server** (`data/acquisition/world_city.py`). Level-14 tiles are reduced from Overture's archives and kept in the cache, read in worker processes so the page stays answered.
- **Heights.** Measured, else floors at 3.2 m, else GHS-BUILT-H's average (a 10 degree tile fetched behind the scenes on first visit), else a default. The building keeps which, and an estimate is drawn paler and says so when picked.
- **Workers** (`world-worker.js`). Parsing, triangulating (`geometry.js`'s `triangulate`, by ear clipping with holes bridged) and placing in the exact local frame are done off the page's thread, returned as typed arrays. Each building stands on the real ground from a grid of the terrain; beside a model, the ground follows the same levelling the shader applies.
- **Scene** (`addBuiltSolids`, `addBuiltPaint`, `addBuiltGround`). Indexed geometry with flat shading from derivatives, about 400 bytes a building. Roads and open land are painted before what stands on them, so they cannot fight the ground. A depth-only grid of the terrain lets a ridge hide what is behind it, and each tile's own grid takes shadows.
- **Bounds.** At most seven tiles across and 420,000 buildings. What was seen stays on disk. The edge dissolves in a stipple.

*Caveat.* Heights are mostly estimates, footprints have no date, and roof shapes are absent. A tile's first visit costs a few seconds, and a region's first visit 20 s more for its heights.

### 2.3 Lenses Before Plug-ins

The page asks a question first: People, Movement, Fabric, Risk, Change. Each layer and tool names its lens in its manifest, and the scales at which it means something, shown as a small altitude tape. An `anywhere` key says what holds beyond the workbench's own cities.

Around the map: a plate with the place's name in a large serif italic; an altitude tape from Earth to Doorstep; a compass and scale bar; a day ribbon drawn from the real Sun; a scale of years when a layer has years; a probe with a leader line; folded legends, one open at a time; layers lifted apart as sheets.

Pink-purple on a plum night, lilac paper by day. The workbench's sea blue is retired, at the user's word; the street's own data keep the `sea` maps.

### 2.4 What Was Integrated

As data or service, on the map:
- **Space2Stats:** H3 level 6 cells of people, lights, built land, urbanisation and hazards, anywhere.
- **GHS-SMOD** and **WSF Evolution** on the ground, the latter growing with a scale of years.
- **Black Marble**, on the globe by GIBS.
- **Poverty and Inequality Platform:** each economy's latest survey, on the world's H3 cells.
- **Data360:** any indicator by meaning, through the API behind data360-mcp.
- **GloFAS river floods:** depth by return period, anywhere.
- **ThinkHazard!:** levels by administrative unit.
- **The Freetown Flood Risk Assessment:** a plug-in of its own, with its maps, 153,435 buildings, wards and loss tables.
- **Risk Data Library Standard:** the assessment's record.

As methods, written here from their statements, none of their code copied:
- H3 binning (giga-spatial), on an H3 ported to physis and held to h3-py;
- composites (GeoE3), in nomos;
- exposure (RAPIDA), in urban-science;
- exceedance and annual loss (CCDR-tools), in tyche;
- reach over ways (GOSTnets), in urban-science on morphe's shortest paths;
- the Urban Space Usage Index (eca-resilience), in urban-science, on Melbourne's counters where it has GPS pings.

Left aside, each with its reason in the workbench's *Sources* tool: montevideo-brt, GRAS, Metadata Editor, Proof-Carrying Numbers, DevDataBench, REaLTabFormer and the guide to synthetic data, MANAGE-WB, the Mobile Phone Data manual, ai4data.

Under the naming policy (ADR 0069), a hexagon holding any part of an area in dispute takes no economy's value or name on these mosaics. The regions product gives such areas to whoever holds them, so that the mosaics cannot rely on it for this.

### 2.5 What the Foundations Gained

- **physis** H3 (decision 0023), **nomos** composites (decision 0006), **tyche** exceedance (decision 0012), **morphe** shortest paths (decision 0011).
- **urban-science** gained `exposure`, `access` and `usage`.
- **data** gained `pmtiles.py`, `shapefile.py`, `world_city.py`, `services.py`, the windows and service manifests, and products `poverty_latest`, `freetown_flood`, `rdls_records`.

## 3. Results

- **Any place.** Kathmandu, with no preparation, opens as a named place with about 154,000 buildings standing on the terrain within reach, 8 % of them with a recorded height or floors.
- **Freetown's annual loss.** Recomputed from the assessment's 1,943 events it is $808,904.91, its own figure to the cent. The sum over the annual loss gives the catalogue's length, 10,000 years.
- **Its curve.** The published occurrence curve is the events in order at every rank a return period falls on (to 1 part in 10⁷). From its dozen published points alone, the area within them is $762,246, 5.8 % below; holding the last loss beyond 1 in 1,000 years adds $82,630, ending 4.4 % above. The tail is an assumption that decides the sign of the error.
- **What a flood reaches.** A 1-in-100-year river flood puts 7,293 of the 153,435 buildings, $384 million, in over 15 cm.

## 4. Consequences

- The workbench needs network access to reach a new place. What it has seen stays.
- The cache under `data/cache/world-city/` holds ODbL data. It must not be redistributed without its terms.
- The city canvas cannot be occluded by the globe's terrain directly, so that each tile brings a depth grid.
- A model and tiles can disagree beside a model's edge. The tile buildings whose middles fall within a model's box are left out.

## 5. Open

- Roof forms and dates of buildings; levels of detail beyond a tile's reach; a time zone map instead of longitude.
- Pedestrian networks in the reach analysis have no signals or gradients; `by road` no traffic.
- Whether the *Fabric* lens should show what the buildings are used for, from Overture's places.
