# Celestial Atlas

*Status: draft (2026-10-22). Kind: tool. Field: astronomy.* Decision records: [ADR 0011](../../../development/decisions/0011-first-application-celestial-atlas.md), [ADR 0013](../../../development/decisions/0013-galaxies-and-their-relations.md), [ADR 0016](../../../development/decisions/0016-maps-of-bodies-and-iau-rotation.md), [ADR 0020](../../../development/decisions/0020-detail-of-bodies-at-close-range.md), [ADR 0021](../../../development/decisions/0021-the-sky-from-a-place.md).

The Celestial Atlas shows the Sun, Moon, planets, stars and galaxies, at any instant from 1849 to 2150 or live. It is organised as five views, from the Solar System out to the cosmic web.

| View | What it shows |
|---|---|
| Solar System | A 3D heliocentric view in the ecliptic frame of J2000, with DE440 positions and osculating orbits. Bodies carry maps of their surfaces, turned to their orientation at the instant shown, with finer detail of the Moon, Mars and Mercury when one fills the view; the Earth shows its city lights at night; the Milky Way lies behind the stars. Bodies can be enlarged by a stated factor. |
| Neighbourhood | Stars in 3D at their measured distances, in galactic coordinates, from a chosen catalogue: HYG (Hipparcos) or any fetched Gaia DR3 subset, against the Milky Way as seen from the Sun. The radius ranges from 5 to 5000 pc. |
| Galaxies | Galaxies and the relations between them, in supergalactic coordinates:<ul><li>**Local Volume** (UNGC): each galaxy is linked to its main disturber, with tidal groups and a Hubble-flow diagram.</li><li>**Local Universe** (2MRS groups): galaxies are linked to their group centres, with named clusters and a mass–light diagram.</li></ul> |
| Sky Chart | The sky for any place and time, drawn as a star atlas, with constellations named in Latin, English, Chinese or Japanese, the Moon in its phase, and the planets. The place is chosen on a small globe, lit as at the instant shown, or by clicking the Earth in the Solar System view; the chart gives local time there and links to the weather in the Earth Atlas. |
| Orbits | The planets seen from the north ecliptic pole, with lines of sight and their configuration |

## 1. Running

From the repository root:

```
python data/acquisition/fetch.py                                    # DE440, HYG, constellations
python data/products/star_catalog.py                                # HYG star product
python data/acquisition/fetch.py gaia-dr3 --profile small           # optional: Gaia subset
python data/products/star_catalog.py --source gaia-dr3 --profile small
python data/acquisition/fetch.py ungc tully-2mass-groups            # galaxies
python data/products/galaxy_catalog.py --source ungc
python data/products/galaxy_catalog.py --source tully-2mass-groups  # names from SIMBAD, cached
python data/acquisition/fetch.py solar-system-maps                  # optional: maps of the bodies and the sky
python data/products/body_maps.py
python applications/astronomy/celestial-atlas/serve.py --open
```

The page also loads three.js from jsDelivr, and the Source Serif 4 font from Google Fonts. Close to the Moon, Mars or Mercury it fetches map tiles from NASA's Trek services.

To choose which part of a catalogue to fetch, and how much of it, see `python data/acquisition/fetch.py --list` and [data/registry](../../../data/registry/README.md).

## 2. Data

| Data | Source | Licence |
|---|---|---|
| Sun, Moon, planets | JPL DE440 (Park et al. 2021), SPK kernel from NAIF | Freely distributed; cite |
| Stars | HYG v4.1 (Hipparcos, Yale, Gliese); Gaia DR3 with Bailer-Jones distances | CC BY-SA 4.0; CC BY-SA 3.0 IGO |
| Constellations | d3-celestial, O. Frohn | BSD 3-Clause |
| Local Volume galaxies | Updated Nearby Galaxy Catalog (Karachentsev et al. 2013), CDS VizieR | Published catalogue; cite |
| Galaxy groups | Tully (2015), from the 2MASS Redshift Survey (Huchra et al. 2012), CDS VizieR | Published catalogue; cite |
| Galaxy and cluster names | SIMBAD (CDS), resolved at build time | Cite SIMBAD |
| Maps of the Sun, the Moon, the planets, Saturn's rings and the Milky Way | Solar System Scope (INOVE), from NASA, JPL and ESA mission data | CC BY 4.0 |
| The Earth by day and by night | Blue Marble and Black Marble (NASA), via GIBS; shared with the Earth Atlas | NASA; acknowledge |
| Detail of the Moon, Mars and Mercury at close range | LRO WAC mosaic (about 100 m), Viking MDIM 2.1 colour (232 m), MESSENGER MDIS (665 m), as tiles from NASA's Trek services, named by the map products | NASA; acknowledge the missions |

## 3. Method

- **Time.**
  - From 1972, TT = UTC + (TAI − UTC) + 32.184 s.
  - TDB is TT plus its main periodic term.
  - UT1 is taken equal to UTC.
  - Before 1972, ΔT follows Espenak and Meeus (2006).
- **Ephemeris.** DE440's type-2 Chebyshev segments are read directly from the SPK kernel and chained to the barycentre.
- **Apparent places.**
  - **Applied:** light-time, IAU 1976 precession, GMST (IAU 2006), topocentric parallax and refraction.
  - **Not applied:** aberration, nutation and deflection, which together amount to under 1′.
- **Orbits.** Osculating two-body conics, computed from state vectors.
- **Where the observer stands.** In the Solar System view the Earth is turned by GMST and precession, the Sky Chart's frame, so a place clicked on it is the place whose sky is drawn; the IAU model of the Earth differs by 0.05–0.26° in longitude over 2000–2100. Clicked points are converted to geodetic latitude on WGS 84, named after the nearest Natural Earth town within 120 km, and given its time zone within 500 km, or nautical time at sea.
- **Orientation of the bodies.** The IAU rotation models (Archinal et al. 2018) give each body's pole (α₀, δ₀) and prime meridian W in the ICRF, including the Moon's periodic terms and those of Mercury and Mars. The rotation to body-fixed axes is B = R_z(W) R_x(90° − δ₀) R_z(90° + α₀).
- **Detail at close range.** When the focused body fills a third of the view or more, rays from the screen find the region in view in body-fixed longitude and latitude; the tiles of the level whose pixels match the screen's (at most 48) are fetched and drawn over the global map by the body's shader, fading in at their edges.
- **Maps.** A body's map is in east longitude with 0 at its centre; the textured sphere's local axes are related to body-fixed ones by (x, y, z)_body = (x, −z, y)_local. The Milky Way map is in galactic coordinates; its orientation was checked against twelve bright stars when it was built. The Earth's night lights show where the Sun is below the horizon.
- **Galaxies.** The theory is in [theory/astronomy/galaxy-relations.md](../../../theory/astronomy/galaxy-relations.md).
  - Positions use supergalactic coordinates.
  - The tidal index Θ₁ = max log(M/D³) − 10.96 is recomputed, along with each galaxy's main disturber.
  - The Milky Way's missing luminosity is inferred from the published indices of its satellites.
  - Tidal groups are the connected components of the relations with Θ₁ > 0.
  - The Hubble flow is drawn with H₀ = 73 km s⁻¹ Mpc⁻¹.
  - Group distances in Tully (2015) are kinematic, in h⁻¹ Mpc, and are labelled as such.

## 4. Verification

```
node applications/astronomy/celestial-atlas/tests/run-tests.mjs
```

There are 66 checks.

| Quantity | Reference | Agreement |
|---|---|---|
| Positions, 1850–2150 | JPL Horizons (DE441) | < 1 mm for the planets, 18 m for the Moon |
| Geocentric astrometric directions | JPL Horizons | < 0.05″ |
| GMST, precession, horizontal coordinates | Meeus, worked examples | Within the stated tolerances |
| Supergalactic directions, 43 036 galaxies | Tully (2015) | Median 0.00012°, maximum 0.0096° (one listed anomaly) |
| Main disturbers, 770 galaxies | Karachentsev et al. (2013) | 99.2 % identical |
| Tidal index Θ₁ | Karachentsev et al. (2013) | Median difference 0.027 dex |
| Sub-Earth points of the Sun, the Moon and the planets, 2000–2036 | JPL Horizons | Within 0.011° (Mercury); the others within 0.003° |
| Tiles for detail, 2000 random regions | Cover the region; the level needed or limited by the count; at most 48 | All |
| Geodetic and geocentric latitude, WGS 84 | Round trip; greatest difference 0.1924° at 45° (Torge and Müller 2012) | 1 × 10⁻¹⁴ °; 2 × 10⁻⁵ ° |
| The globe's orthographic projection, 2000 random points | Round trip (Snyder 1987) | 2 × 10⁻¹⁴ rad |

## 5. Structure

The computation in `src/core/` is a provisional stand-in for the laboratory's C++ core ([ADR 0011](../../../development/decisions/0011-first-application-celestial-atlas.md)).

| Path | Role | Stands in for |
|---|---|---|
| `src/core/time.js` | Time scales, Earth rotation angle, GMST | `core/metrology` (time) |
| `src/core/frames.js` | Rotations, precession, horizon, galactic and supergalactic frames | `core/metrology` (frames) |
| `src/core/spk.js` | DE440 reader and evaluator | `core/domains/astronomy` |
| `src/core/sky.js` | Apparent places, elongation, phase | `core/domains/astronomy` |
| `src/core/kepler.js` | Osculating orbits | physis (mechanics) |
| `src/core/galaxies.js` | Tidal index, main disturbers, inferred masses, tidal groups | `core/domains/astronomy` |
| `src/core/rotation.js` | IAU rotation models of the Sun, the Moon and the planets | `core/domains/astronomy` |
| `src/core/tiles.js` | Tiles of an equirectangular map service for a region and a resolution | `core/domains/geography` (cartography) |
| `src/core/geodesy.js` | Geodetic and geocentric latitude, distances, the nearest place | `core/domains/geodesy` |
| `src/data/catalog.js` | Loads products: stars, constellations, galaxies | — |
| `src/render/` | Palette, star points, label placement, maps of the bodies and the sky, detail at close range, the globe for choosing a place | `visualization/` |
| `src/views/` | The five views and the galaxy diagrams | `visualization/renderers` |

## 6. Limits

- UT1 − UTC and future leap seconds are not modelled.
- Star positions are drawn at their catalogue epoch, without proper motion.
- Gaia subsets taken brightest-first are limited by magnitude, not complete within a distance.
- The browser loads the 32 MB DE440 kernel, and about 17 MB of maps after the first frame.
- The maps are global images at 4096 pixels (the sky at 8192). Finer detail exists for the Moon, Mars and Mercury only, and the tile mosaics differ a little in tone from the global maps.
- The tile services do not wrap round 180° of longitude: a view across it shows detail on its larger side only.
- The Earth's day map is the Blue Marble, a cloud-free composite, not the Earth on the day shown.
