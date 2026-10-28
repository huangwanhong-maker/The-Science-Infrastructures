# 0049. The Spacetime of the Earth and the Moon: the Domains' C++ of Time, Frames and Figures, the Navigation Bridge, and the Globe

- **Status:** Accepted, provisionally (the synthetic detail is a texture, not a model of terrain; light time and aberration are not applied to the frames; the GCRS's relativistic scale is not modelled)
- **Date:** 2026-10-28
- **Decided by:** the user asked Cosmopoiesis for its first substantive work: "(1) our spacetime model + space data (or spacetime data, if the data can be obtained and updated in real time) + visualization", with a frontend like the Earth Atlas or Google Earth, its multi-resolution rendering, the Earth's, the ocean's and the Moon's data and something better than interpolation where data are missing; and "(2) coordinate systems and conversion", "like aircraft and spacecraft navigation systems … various systems, and convert among them". They agreed that the domains' C++ be written for metrology, geodesy and the Moon's frames, that the navigation bridge be founded, that the renderer live in the shared visualization layer, and that sub-resolution terrain be synthetic, labelled and off by default. Claude designed and wrote it, with two ports by agents it directed (Karney's geodesics, SGP4), and this record.

## 1. Context

The domains of time, the Earth and the sky were provisional JavaScript (ADR 0011, 0014): good for drawing, not for a world whose frames must agree to the metre and the milliarcsecond. The navigation bridge was planned. physis had no rotations or moving frames of its own beyond SU(2), and its `mechanics` was planned. The data library had no Earth orientation, leap seconds, geoid, terrain tiles or satellites. The visualization layer had flat map tiles and no globe.

Cosmopoiesis (its decision 0001) holds that a world's frames form a groupoid and that observables are invariants; its first slice had been planned as a thermal box. The user moved spacetime and coordinates first.

## 2. Decision

### 2.1 The Foundations

physis gained (its decision 0009) rotations in every form, rigid motions with their rates, and the frame graph: a tree of frames joined by time-dependent motions, a relation closing a loop between two routes, and the loop's holonomy measured by Chasles's angle and slide; and `mechanics`, the two-body problem (elements, Kepler's equation, propagation by universal variables, the frames of an orbit).

### 2.2 The Domains

| Domain | Modules written in C++ |
|---|---|
| metrology | `time`: every time scale after SOFA, leap seconds (IANA's list or built in), TDB − TT by Fairhead and Bretagnon, UT1 from the IERS or approximated and saying so, ISO 8601 with scales; `earth_orientation`: the IERS's finals2000A; `frames`: IAU 2006/2000A in full (nutation and CIO series transcribed from ERFA), the Earth's frames as a frame graph with the equinox route a loop, GMST of 1982 |
| geodesy | `figure`: ellipsoids, geodetic coordinates both ways, local frames and look angles; `geodesic`: Karney 2013 (from GeographicLib, MIT); `projection`: the transverse Mercator by Krüger's series to n⁶, UTM and UPS, MGRS; `geoid`: EGM96's grid |
| astronomy | `ephemerides`: DAF, SPK types 2 and 3, binary PCK; `rotation`: the IAU models and the Moon's DE440 frames |
| navigation (founded) | `frames`: the catalogue of the Solar System's frames as one graph (ICRS, GCRS, the Earth's frames, TEME, the Moon's, the planets', the ecliptic, the galaxy, the Earth–Moon rotating frame, places), with loops between standards; `orbits`: two-line element sets, SGP4/SDP4 (from python-sgp4, MIT, after Vallado 2006), TEME; `vehicle`: the aerospace Euler angles, rates, track and flight path, air data |

Every module is checked against an authority, its reference values in its tests' fixtures with the script that made them: SOFA's test values and ERFA, JPL Horizons, NAIF's worked example, GeographicLib, PROJ, GEOTRANS, Vallado's verification set and skyfield. Where two standards define one frame differently, the catalogue does not pick one silently; it closes a loop, and the holonomy states the difference: the equinox and CIO routes to the ITRS (0.3 mas today, the celestial pole offsets), the IAU's Moon against DE440's (about 12″), GMST of 1982 against IAU 2006's in TEME (tens of milliarcseconds).

### 2.3 The Data

| Dataset | Use |
|---|---|
| `naif-lunar-frames`, `iers-eop`, `iana-leap-seconds` | the Moon's orientation, the Earth's, and UTC; the last two refreshed by fetching again |
| `egm96-geoid` (product `geoid`) | heights above mean sea level |
| `aws-terrain-tiles`, `lola-gdr`, ETOPO 2022 (product `fields`), Natural Earth (product `land-mask`) | the terrain of the Earth and the Moon |
| `usgs-moon-nomenclature`, Natural Earth (product `places`) | the names of places and features |
| `celestrak-gp` | satellites' element sets, by group, refreshed every few hours |

The new module `data/tiles` builds terrain tiles on demand in the globes' geographic pyramid. Each sample carries a provenance code: measured, filled by the source, compiled, finer than the data, synthetic, or the sea's surface. The sea is drawn at its surface (the geoid) with its floor kept beneath. Synthetic detail below the data's resolution is made only on request. It is a self-affine extrapolation of the region's roughness (H = 0.8 on the Earth, 0.9 on the Moon), a function of the point so that tiles agree at their seams, each octave entering over two levels of the pyramid so that the levels do not show as patches, and it is labelled.

### 2.4 The Visualization Layer

The layer gained the namespace `globe`: the pyramid; the figures; terrain meshes ordered by quarter; the level of detail; imagery services; a virtual globe's camera; terrain workers; sampled motions; and three.js renderers of a globe and of the sky. Its level of detail follows Cesium's quadtree, with five choices of its own:

- the screen-space error lessened toward the horizon;
- the horizon test by a tile's culling point;
- a walk through tiles not yet loaded, so that a view brought close asks for its levels at once;
- a parent drawn only in the quarters its children cannot yet cover, with skirts under the lines between its quarters as well as at its edges;
- proxy bounds narrowed by the nearest loaded ancestor's heights.

A close view settles in a few seconds rather than a level per round trip. The renderer draws relative to the camera in double precision, with a logarithmic depth buffer. The colour of a star moved here from the Celestial Atlas, its second user.

### 2.5 The Work

Cosmopoiesis's engine (`engine/cosmopoiesis.cpp`) is kept alive by its server and answers nine operations:

- `status`, `time`, `frames`, `ephemeris`, `transform` and `holonomy`;
- `convert`: a point in every system;
- `geodesic`;
- `satellites`: SGP4 into the GCRS, with the shadow and look angles.

Its page draws both globes and the sky at any instant, refined as one comes near. Each tab does one job:

- **Layers:** imagery and overlays, terrain options, satellites.
- **Time:** the scales and the Earth's rotation.
- **Place:** the chosen point in every system, an observer, the geodesic, and going to coordinates.
- **Frames:** the catalogue, conversions and loops.
- **Data:** sources and the recipe.

The page draws only when something changes.

*Rationale.* The page does not recompute the frames: it asks the engine for windows of samples and interpolates them (slerp with the Earth's spin factored out, Hermite paths), so that there is one implementation of each frame, the C++, and the page draws at the frame rate. A WebAssembly build of the engine for the page was considered and set aside for now. Borrowing CesiumJS for the globe was set aside too: the layer owns its interfaces (CLAUDE.md §2.2), the atlases draw with three.js, and the level of detail is small enough to own.

## 3. Consequences

- The domains of the Earth and the sky have C++ beneath their JavaScript; the JavaScript remains for the pages that use it until they move to the engines.
- Any work can compose the catalogue of frames, the conversions and SGP4 through the domains' C++, and draw a globe through the layer's `globe`.
- The root `CMakeLists.txt` composes the new targets in the domains' order.
- Open: apparent places (light time, aberration) as their own frames in astronomy; a WebAssembly or WebGPU path for the engine's frames in the page; the lake surfaces that ETOPO's coarsest levels draw as basins; EOX's newer cloudless mosaics are licensed for non-commercial use only and are not offered; the age of satellites' elements bounds their accuracy (kilometres a day), which the page says; on the Moon, a population of small craters in equilibrium would be a better synthetic detail than noise.
