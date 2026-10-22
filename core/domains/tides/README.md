# tides

*Status: planned.* This bridge joins astronomy, geodesy, oceanography and meteorology.

## 1. Scope

| Module | Scope |
|---|---|
| `potential` | The tide-generating potential in spherical harmonics, the equilibrium tide, and the Doodson arguments |
| `constituents` | Tidal constituents and nodal modulation |
| `analysis` | Tidal analysis and prediction |
| `dynamic` | Dynamic ocean tides |
| `atmospheric` | Atmospheric tides |

## 2. Dependencies

| Source | What tides uses from it |
|---|---|
| `astronomy` | Ephemerides of the Moon and the Sun |
| `geodesy` | The figure and rotation of the Earth |
| `oceanography` | Bathymetry and ocean dynamics |
| `meteorology` | The atmosphere's response, for atmospheric tides |
| physis | `gravitation` |

## 3. Applications

See [`applications/tides/`](../../../applications/tides/). Related seed idea: [the tidal spectrum as a lattice](../../../development/ideas/2026-10-21-tidal-spectrum-as-lattice.md).
