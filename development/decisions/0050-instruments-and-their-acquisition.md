# 0050. Instruments and Their Acquisition: Addresses of Events, the Models Instruments Read, and Scheduling with Deadlines

- **Status:** Accepted, provisionally (the timing is soft real time, measured and not guaranteed; the kernels of observation are Cosmopoiesis's until a second work needs them)
- **Date:** 2026-10-28
- **Decided by:** the user asked Cosmopoiesis for one world of the Solar System, "extendable to the galaxy", in which "a global spacetime addressable system" converts between every coordinate system "by our unified framework … algebraic structure and field theory based. May be gauge theory"; and for instruments, an abstract sensor "with its physical property and available spacetime range", physical and virtual, read "when clock ticks" and on events, scheduled with priorities, in parallel and in batches, "no confliction, safety", without starvation. Claude designed it and wrote it, with two modules ported by agents it directed (the World Magnetic Model; the standard atmosphere and normal gravity), and this record. The work's own records are Cosmopoiesis's decisions [0003](https://github.com/huangwanhong-maker/Cosmopoiesis-Lab/blob/main/development/decisions/0003-one-world-its-addresses-and-a-console.md) and [0004](https://github.com/huangwanhong-maker/Cosmopoiesis-Lab/blob/main/development/decisions/0004-instruments-and-their-acquisition.md); calx's is its [decision 0009](https://github.com/huangwanhong-maker/calx/blob/main/development/decisions/0009-real-time-scheduling.md).

## 1. Context

After ADR 0049, any point of the Earth or the Moon could be converted between the catalogue's frames, but only through an operation per system, and only near the two bodies. There was no way to write an arbitrary event of the Solar System in one line and read it back. Instruments had nothing to read beyond the ephemerides and the figure: the infrastructure had no geomagnetic field, no standard atmosphere in C++, and no normal gravity. And scheduling with deadlines had no home: the neuroscience engine has phases, the Genome Workbench has jobs, and neither has a deadline.

## 2. Decision

### 2.1 The Foundations

- **calx** gained `compute/scheduler/realtime.hpp` (its decision 0009). It holds:
  - the periodic task model;
  - the analyses: utilisation, Liu and Layland's bound, rate- and deadline-monotonic orders, response times with blocking, the busy period, the demand bound and EDF's feasibility;
  - constant bandwidth servers and elastic compression;
  - a deterministic simulation on several processors, by EDF, fixed priority with aging, or FIFO;
  - a non-preemptive dispatcher on the pool of threads.

  Its 34 checks hold the analyses to the simulation.
- **calx** also gained `special_functions/legendre.hpp`, the associated Legendre functions in Schmidt's semi-normalisation, with their derivatives and their quotient by cos φ finite at the poles. They were written in the geophysics domain and moved down, since they name no object in nature; they are held to the standard library's `std::assoc_legendre`.

### 2.2 The Domains

| Domain | Modules written in C++ |
|---|---|
| navigation | `frames/address.hpp`: the address language. A body or a frame, a chart (cartesian, spherical, cylindrical, geodetic or planetocentric, UTM, MGRS), coordinates, steps from a local horizon (ENU, NED, AER) and an instant. It is LL(1), with a canonical print in the shortest round-trip decimals, resolution to an event and expression of an event in any frame and chart. `frames/galactocentric.hpp`: the galactocentric frame with astropy's parameters (v4.0), under the ICRS in the catalogue. Its checks hold the frame to astropy at 300 stars, to $6\times10^{-16}$ relative in position and $4\times10^{-9}$ m/s in velocity |
| geophysics | `geomagnetism`: the World Magnetic Model 2025, with the field and its secular variation at a geodetic place, the Earth-fixed vector, the centred dipole and McIlwain's $L$, the model's validity and its error model. It is held to the report's high-precision example to $10^{-10}$ nT and to NOAA's test values. It found three misprints of the report (eq. 15's sign, the north pole's meridian east and not west, eq. 29's δX for δZ) and states them, and NOAA's file of test values departs from the report's equations by $7\times10^{-4}$ nT. The domain now depends on geodesy and metrology |
| meteorology | `atmosphere`, the domain's first C++: the U.S. Standard Atmosphere 1976 from −5 to 86 km, with its temperatures, pressure, density, the speed of sound, viscosities, conductivity, number density, mean free path and collision frequency, and columns for instruments rising through it. It is held to the document's Tables I–III, transcribed from its scan at thirty-eight heights (the script records the scan's SHA-256), to fluids' and to ICAO's; it applies the correction above 80 km that the document prescribes and its own tables omit |
| geodesy | `gravity`: normal gravity of a level ellipsoid (WGS 84, GRS 80, or any four defining constants), with its derived constants as NGA and Moritz print them, Somigliana's formula and the closed field at any height in ellipsoidal-harmonic coordinates |

Every module is checked against its authority, with its reference values in fixtures written by a script beside them, as CLAUDE.md §2.5 asks; the agents' transcriptions say where they were done by eye.

### 2.3 The Data

| Dataset | Use |
|---|---|
| `noaa-wmm` | WMM2025's coefficients, NOAA's test values and the technical report |
| `noaa-swpc` | space weather in real time, from SWPC's JSON feeds: the solar wind and the interplanetary field at L1, the planetary $K$ index, GOES's X-rays and energetic protons |

### 2.4 The Visualization Layer

`globe` gained:

- flights between bodies (`transferCamera`, `slerpVec`, `ease`): the distance interpolated in its logarithm and raised to see both bodies;
- a body's map as the two root tiles of the geographic pyramid, finer close in from a WMTS (`mapLayer`);
- the Sun's lighting with limb darkening;
- the Milky Way drawn in galactic coordinates;
- the camera's range bounded per body.

### 2.5 The Work

Cosmopoiesis draws one world of every body of DE440 in the ICRS, addressed in the navigation domain's language, on a console laid out as a flight deck (its decision 0003). It reads that world through instruments: twelve kernels of observation, nine families of virtual instruments on the domains' models and ten physical sources. It schedules them by a runtime in its engine and its server, which composes calx's scheduling (its decision 0004). Its formalism gains the events and their addresses, the frames as a connection whose holonomies measure the disagreement of standards, the kinds of phenomena on three axes, the kernels, provenance and the determinism of readings, and the timing of acquisition with its proposition against starvation.

*Rationale.* Each piece went where it names its objects:

- The scheduling and the Legendre functions name nothing in nature, so they went to calx.
- The address language names the bodies and frames of the Solar System, so it went to the navigation bridge, beside the catalogue it reads.
- The models instruments read went to the domains whose objects they are.
- The kernels of observation stay in Cosmopoiesis: one work uses them, and the infrastructure writes an abstraction only after its pattern has appeared twice (CLAUDE.md §2.2).

## 3. Consequences

- Any work can write and read an event of the Solar System or the Galaxy in one line, convert it to every frame and chart, and read the catalogue's holonomies as the disagreements of its standards.
- Any work can evaluate the geomagnetic field, the standard atmosphere and normal gravity in C++, and schedule duties with deadlines through calx.
- The root `CMakeLists.txt` composes geophysics and meteorology, with checks (`GEOPHYSICS_CHECKS`, `METEOROLOGY_CHECKS`).
- Open questions:
  - the kernels of observation in a foundation, with their registry in the observation edge, when the Earth Atlas's observations or the Brain Atlas's forward models join them;
  - a provenance for indices derived from observations without a model (the $K$ index);
  - a backend of calx's scheduler on a real-time kernel;
  - the relativistic transformation between the BCRS and the GCRS;
  - frames beyond the Galaxy, which are boosts;
  - fully normalised Legendre functions of very high degree, with extended range, for gravity models such as EGM2008 and GRGM1200A.
