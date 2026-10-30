# 0055. The Earth Beneath: CRUST1.0 over PREM, Sections along Geodesics, and the World below the Ground and the Sea

- **Status:** Accepted, provisionally. The Earth only. The Moon's interior, the ocean's water column and the C++ of the rest of the geophysics domain are open (§5).
- **Date:** 2026-10-30
- **Decided by:** the user, who asked for the environment page (ENV) of Cosmopoiesis: "we may even visualize under ground/sea geography". The design is Claude's, made after the precision console of Cosmopoiesis's decision 0009.

## 1. Context

Cosmopoiesis's ENV page showed what the observer's instruments read above the ground: the air, gravity, the magnetic field, the sky, space weather. Nothing showed what lies below the ground or the sea. The world could already take the sea away and show its floor (the LAYERS page; the terrain tiles' second layer, ADR 0049). The infrastructure held the pieces of the rest, but not where Cosmopoiesis could use them:
- **The crust.** The geophysics domain's `structure` (ADR 0024) read CRUST1.0 and cut sections along great circles, in provisional JavaScript for the Earth Atlas. Cosmopoiesis computes in C++ in its engine.
- **The deep Earth.** No model reached below the Moho.
- **The data.** The data library held CRUST1.0 as the product `crust-crust1`, Global CMT's 70,451 earthquake centroids since 1976 (`events-gcmt`), ETOPO 2022 and the Terrarium tiles with their sea floor.

## 2. Decision

### 2.1 What Went Where

| Where | What |
|---|---|
| geophysics `structure` (C++, new beneath the JavaScript) | `crust.hpp`: CRUST1.0's columns, the same as `js/crust.js` (its cells, its layers present, the Moho, the thicknesses); `prem.hpp`: PREM (Dziewonski and Anderson 1981, Table I, isotropic), its shells, discontinuities, the mass within a radius integrated exactly, gravity and the moment of inertia; `interior.hpp`: the column composed of the two, as a profile to the centre and as physis fields of density and of the speeds of P and S waves; `section.hpp`: a section along a geodesic of the ellipsoid (geodesy's geodesics), and where a place lies along and across it |
| Cosmopoiesis's engine | the instrument `seismic_column.crust1` (the encyclopedia's B09, elastic waves); the operation `beneath.section`; the world's fields `earth.interior.density`, `earth.interior.vp` and `earth.interior.vs` (its decision 0010) |
| Cosmopoiesis's server | the terrain along a section (the data library's tiles: ETOPO 2022 or Terrarium, the sea told from land by Natural Earth's land); Global CMT's earthquakes near it; CRUST1.0 and Global CMT in its recipe |
| Cosmopoiesis's page | ENV's section *Beneath*: the column's readings, a section through the observer in three views, the column to the centre, both separable; the section's trace on the world |

Nothing general was written in a domain. The geodesics are geodesy's; the minimiser that finds the foot of a place on a section is calx's (Brent); the fields are physis's.

### 2.2 The Column: CRUST1.0 over PREM

A place takes the column of its 1° cell of CRUST1.0, from the top of its water or ice to its Moho. Below the Moho, the column is PREM's. This is the usual composition of a crustal model over a radial one:
- **Depth.** Depth is measured below mean sea level, along the ellipsoid's normal. PREM's spherical depth is taken as that depth, without the ellipticity correction; the error is of the order of the flattening, at most 21 km near the centre.
- **A shallow Moho.** Where the Moho is shallower than PREM's own (24.4 km, under the oceans), the mantle beneath is PREM's lid, its polynomial continued upward.
- **PREM's crust is not used.** PREM's crust and ocean are global averages; CRUST1.0 replaces them.

The column reads, for each depth: vp, vs, the density, and PREM's gravity G M(r)/r². It also reads the time a P wave takes to rise from the Moho to the surface, and the layer the mount is in when it is below the surface.

### 2.3 The Section

- **The line.** A section is the geodesic of WGS 84 through a place at an azimuth, sampled from −half to +half.
- **Where an earthquake lies.** A place's distance along the section is that of its foot: the point of the line nearest to it. Brent's minimisation finds the foot to √ε |s|, and two steps of the spherical right triangle at the foot take it to rounding. The place's distance across is the length of the geodesic from the foot to it, signed by its side.
- **Two resolutions, composed.** The surface, the sea floor and the sea come from the terrain, a sample every step. The layers beneath come from CRUST1.0, constant in each cell. A layer is drawn where the model puts it and cut by the terrain above it, and the first solid layer of a cell reaches up to ground that stands above the cell's mean. The figure says so. Its layers step at the cells' edges; they are not smoothed into a picture the model does not make.
- **SRTM's water.** Where Terrarium has exactly 0 m at sea (SRTM writes water as 0), the sea floor is ETOPO 2022's.

### 2.4 Verification

The structure module's checks (`structure/tests/check_structure.cpp`, 30 checks):
- **PREM against ObsPy's tabulation** (`tests/fixtures/prem-obspy.txt`, written by `tests/make_prem_fixtures.py` from `obspy/taup/data/prem.nd` at a fixed commit, its SHA-256 checked):
  - both sides of each of its seven discontinuities, to 2·10⁻⁴ km/s and g/cm³;
  - in the lid and the low-velocity zone (24.4 to 220 km), where the table is not linear in the radius as PREM's isotropic polynomials are, the two depart by up to 1.2·10⁻³ km/s, and are compared apart.
- **PREM's mass**, integrated exactly, against WGS 84's GM over CODATA 2018's G, to 1.7·10⁻⁴. Its gravity is greatest at the core–mantle boundary, 10.69 m/s². Its moment of inertia comes to 0.3308 M a².
- **CRUST1.0** passes the JavaScript's checks (a constructed oceanic column, the cells' corners, Tibet, the central Pacific, the Canadian Shield, Pn). All 64,800 cells are contiguous from their tops to their Moho, every layer slower than the mantle beneath.
- **Sections**, against geodesy's direct problem: places built 80 km to the right and 120 km to the left of the line, and on it behind its centre, are found along and across to 10⁻⁹ m; successive samples lie a step apart.
- **The interior field** is refused in the air and beyond the centre, each with its reason.

Cosmopoiesis's checks hold the engine's section, its events, its columns and the interior's fields against the same fixture (its decision 0010).

## 3. Alternatives Considered

- **Use the domain's JavaScript from the page**, as the Earth Atlas does. Rejected: the infrastructure's computation is C++ wherever possible (ADR 0003), the crust is needed by an instrument of the engine, and C++ was already beneath the other domains Cosmopoiesis reads (ADR 0049).
- **Interpolate CRUST1.0 between cells** for a smooth section. Rejected: the model's own readme asks that a cell's values be taken as its average. A smooth picture would show structure the model does not have.
- **Draw PREM's crust and ocean.** Rejected: they are global averages, and CRUST1.0 is local.
- **AK135 or IASP91 in place of PREM.** Deferred: they are tabulated models, fitted to travel times; PREM is polynomial and fitted to the Earth's mass and free oscillations as well, so its mass, gravity and moment of inertia can be checked exactly.
- **A sphere for sections**, as `js/section.js` uses. Rejected for the C++: the geodesic of the ellipsoid costs nothing more here and is exact.

## 4. Consequences

- **The geophysics domain** has its second C++ module, `geophysics::structure`. The suite runs 128 checks. The JavaScript stays for the Earth Atlas until it reads the engine's answers.
- **Cosmopoiesis** reads the Earth beneath a place as an instrument, as a section and as fields (its decision 0010). It has thirty operations and thirteen fields. The encyclopedia's B09 is written.
- **The data recipe** of Cosmopoiesis fetches CRUST1.0 and Global CMT and builds their products.

## 5. Open

- **The Moon's interior.** GRAIL's crustal thickness (Wieczorek et al. 2013) and a seismic model of the Moon's interior (Weber et al. 2011). Until then a lunar section shows LOLA's terrain and says why it shows nothing beneath.
- **The ocean's water column.** Temperature and salinity from the World Ocean Atlas, the speed of sound by TEOS-10 (oceanography's planned `seawater`), and the SOFAR channel.
- **The terrain tiles' sea floor.** The tiles of the data library still take SRTM's 0 m at sea as the floor near coasts. The section and the probe correct it (§2.3), and the tiles should too.
- **The ellipticity correction** of PREM's depths, and a three-dimensional mantle (tomography).
- **The rest of geophysics in C++**: plates, mechanisms and aftershocks.
