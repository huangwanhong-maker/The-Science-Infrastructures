# domains: Fields and the Bridges Between Them

*Status: planned.* Layer 2 of `core/`.

## 1. Purpose

This directory holds the fields that the laboratory investigates together, and the topics that bridge them.

- **Fields.** A field applies the mathematics and physics of the libraries to particular bodies and systems: the sky, the Earth's surface, the ocean and the atmosphere.
- **Bridges.** A bridge joins two or more fields. It declares the fields it joins, and depends on them.

Every field and bridge has its own directory from the start ([ADR 0010](../../development/decisions/0010-domain-organised-fields-and-applications.md)). The modules inside them are created when their first content is written. `applications/` mirrors this structure.

## 2. Fields

| Field | Scope | Status |
|---|---|---|
| [`astronomy/`](astronomy/) | Positional astronomy, ephemerides, astrometry; later astrophysics and cosmology | planned |
| [`geography/`](geography/) | All branches: physical geography, cartography and GIS, remote sensing, human geography, and others as the work requires | planned |
| [`oceanography/`](oceanography/) | Physical, chemical and biological oceanography; seawater thermodynamics; the seafloor | planned |
| [`meteorology/`](meteorology/) | Meteorology in the broad sense of atmospheric science: dynamics, moist processes, radiation, chemistry, weather prediction | planned |
| [`geodesy/`](geodesy/) | The figure of the Earth, its gravity field and rotation, and reference systems. The other fields rest on it. | planned |

## 3. Bridges

| Bridge | Joins | Scope | Status |
|---|---|---|---|
| [`tides/`](tides/) | astronomy, geodesy, oceanography, meteorology | The tide-generating potential, ocean and atmospheric tides, and tidal analysis and prediction | planned |
| [`navigation/`](navigation/) | astronomy, geodesy, geography | Celestial navigation, charts and routes; later GNSS | planned |
| [`air_sea/`](air_sea/) | meteorology, oceanography | Exchange of momentum, heat, freshwater and gases at the sea surface; wind waves | planned |
| [`climate/`](climate/) | meteorology, oceanography, geography, astronomy | The coupled system, energy balance, climatic zones, variability, and orbital forcing | planned |

## 4. Dependencies

- **Between fields.** A field may depend on another field when the dependency is declared and acyclic. For example, geography depends on geodesy.
- **Bridges.** A bridge depends on the fields it joins.
- **No field depends on a bridge.**
- **Shared physics.** Ocean and atmosphere are treated as one coupled system ([ADR 0007](../../development/decisions/0007-atmosphere-and-ocean-as-coupled-system.md)). The fluid dynamics they share lives in the physis library, in its `continuum` module.

## 5. Reference Models and Systems

Mature, open models and systems show how physics, data and observation are combined in practice.

- **Ocean:** the models MITgcm, MOM6, NEMO, ROMS, FESOM and Oceananigans.jl, and the ECCO state estimate.
- **Meteorology and the coupled system:** WRF, MPAS and ICON; the Earth-system models CESM and E3SM; the ERA5 reanalysis; and machine-learning weather models such as GraphCast and ECMWF's AIFS.
- **Geography:**
  - PROJ and GDAL (permissive licences), which may be linked;
  - QGIS and GRASS GIS (GPL), which we study but do not link.

We study these systems, reimplement their algorithms in this framework where useful, and validate against them.

- Reimplementing an algorithm is always allowed.
- Copying code requires a compatible licence.

See also the [idea on learning from ocean models](../../development/ideas/2026-10-21-learning-from-ocean-models.md).
