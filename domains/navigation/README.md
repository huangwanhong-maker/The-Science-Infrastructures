# navigation

*Status: draft. Bridge, height 3.* Finding one's place and way, on the Earth, in its air and above it: the frames in which a place, a vehicle and an orbit are stated, the conversions between them, and the motion of what flies, by the sky, by reckoning, by charts and by orbits.

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `frames` | The catalogue of the Solar System's frames as one physis frame graph, from the barycentre's ICRS through the GCRS and the Earth's frames to the ITRS, TEME, the Moon's mean-Earth frame, the Earth–Moon rotating frame, the ecliptic, the galactic and the galactocentric frames (astropy's, `galactocentric.hpp`), and places on the Earth and the Moon; loops between standards, whose holonomies measure their disagreement. Addresses (`address.hpp`, [ADR 0050](../../development/decisions/0050-instruments-and-their-acquisition.md)): the language in which an event of the Solar System or the Galaxy is written in one line (a body or a frame, a chart, coordinates, steps from a local horizon, an instant), LL(1), printed canonically, resolved to an event and expressed in any frame and chart | draft (C++) |
| `orbits` | Two-line element sets; the SGP4/SDP4 propagator; the TEME frame and its rotation to the Earth | draft (C++) |
| `vehicle` | A vehicle's frames: the aerospace Euler angles, body and Euler rates, ground speed, track and flight path, air data, the body frame on the Earth | draft (C++) |
| `celestial` | Celestial navigation and sight reduction | planned |
| `dead_reckoning` | Dead reckoning | planned |
| `charts` | Charts and routes: great circles and rhumb lines (the geodesic itself is geodesy's) | planned |
| later | GNSS: satellites' broadcast orbits, pseudoranges and positioning | planned |

## 2. Concepts

- **One graph of frames.** Every frame is a node of a tree joined to its parent by a motion that depends on the instant (metrology's `Instant`, which carries every time scale and the Earth's orientation); the transformation between any two is composed along the path between them, as NAIF's SPICE composes its frames. Where two standards define the same frame by different routes (the equinox and CIO routes to the ITRS; the IAU's model of the Moon and DE440's libration; GMST of 1982 and IAU 2006's in TEME), a relation closes a loop, and its holonomy, the motion once round it, says how far they disagree: under a milliarcsecond, some arcseconds, tens of milliarcseconds.
- **TEME** is SGP4's frame: the true equator and mean equinox of date, turned to the Earth by GMST of 1982 alone (Vallado, Crawford, Hujsak and Kelso 2006, AIAA 2006-6753). Its sidereal time is metrology's.
- **The Earth–Moon rotating frame** is the synodic frame of the restricted three-body problem: the RTN frame of the Moon's orbit about the Earth (physis's `orbit_frame`), its origin moved to the barycentre that DE440 gives, its rate from the ephemeris's acceleration.
- **SGP4** propagates the mean elements of a two-line set: a port of python-sgp4 (Brandon Rhodes, MIT), itself Vallado's code of 2006, reproducing his verification set to the digits it prints.
- **A vehicle** carries the body axes of aeronautics (x forward, y to starboard, z down), turned from the local north-east-down frame by the yaw, pitch and roll of the 3-2-1 sequence (ISO 1151; Stevens, Lewis and Johnson 2016).

## 3. Dependencies

| Source | Used for |
|---|---|
| astronomy | the ephemerides, the IAU rotation models, the Moon's frames |
| geodesy | the figure of the Earth, places and their local frames; geodesics and grids for what is built on them |
| metrology (beneath astronomy and geodesy) | time scales, the Earth's orientation and frames, GMST of 1982 |
| geography | cartography (for `charts`, planned) |
| physis | rotations, moving frames and the frame graph (`geometry`); the two-body problem and the frames of an orbit (`mechanics`) |
| calx | the checks' verification |

## 4. Code

The C++ is header-only, targets `navigation::frames`, `navigation::orbits` and `navigation::vehicle` ([ADR 0049](../../development/decisions/0049-the-spacetime-of-the-earth-and-the-moon.md)). Its checks:

- `frames/tests/check_catalogue.cpp`: the catalogue's frames and paths; NAIF's worked example of the Moon's frames through the catalogue; the holonomies of its loops (the IAU Moon, TEME); places on the Earth and the Moon; the Earth–Moon rotating frame (the Moon on its x axis, its origin the barycentre by DE440's mass ratio to 10⁻⁹).
- `frames/tests/check_address.cpp`: the grammar's examples read and printed canonically, `print(parse(w)) = w`; syntax errors with their columns; every chart's round trip; steps from a horizon; the galactocentric frame against astropy 8.0.2 at 300 stars (`tests/fixtures/galactocentric-reference.txt`, made by `tests/make_frames_fixtures.py`), to $6\times10^{-16}$ relative and $4\times10^{-9}$ m/s.
- `orbits/tests/check_sgp4.cpp`: Vallado's verification set (`tests/fixtures/SGP4-VER.TLE` and `tcppver.out`, copied from python-sgp4 unaltered; their source in `tests/fixtures/README.md`), every case and error code, to the 10 µm the file prints; against Vallado's compiled code through python-sgp4 in both operation modes and gravity models; TEME to the ITRS against Vallado's worked example, verified with skyfield.
- `vehicle/tests/check_attitude.cpp`: the Euler angles' round trips and gimbal lock, agreement with physis's intrinsic Z-Y-X rotation, rates, ground speed and track, air data, and the body frame on the Earth.

```
node domains/navigation/tests/run-tests.mjs
```

The catalogue's checks read the ephemeris and the Moon's kernels from the data cache and skip what is missing.

## 5. Used By

Cosmopoiesis (its spacetime, coordinates and satellites); the sight-reduction trainer and the route planner (planned applications).
