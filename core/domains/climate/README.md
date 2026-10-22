# climate

*Status: planned.* This bridge joins meteorology, oceanography, geography and astronomy.

## 1. Scope

| Module | Scope |
|---|---|
| `energy_balance` | Energy-balance models |
| `coupled` | The coupled ocean–atmosphere system |
| `variability` | Variability such as ENSO and the monsoons |
| `zones` | Climatic zones and their classification. A provisional implementation of the Köppen–Geiger classification is in the Earth Atlas (`applications/geography/earth-atlas/src/core/climate.js`); the theory is in [theory/climate/koppen-geiger-classification.md](../../../theory/climate/koppen-geiger-classification.md). |
| `orbital_forcing` | Insolation and Milankovitch cycles |
| later | The cryosphere and the land surface |

## 2. Dependencies

| Source | Used for |
|---|---|
| `meteorology` | The atmosphere |
| `oceanography` | The ocean |
| `geography` | The land surface and climatic zones |
| `astronomy` | Orbital elements, insolation |
| `air_sea` | Coupling between ocean and atmosphere |

## 3. Applications

See [`applications/climate/`](../../../applications/climate/).
