# air_sea

*Status: planned.* Bridge joining meteorology and oceanography.

## 1. Scope

| Module | Scope |
|---|---|
| `fluxes` | Fluxes of momentum, heat, freshwater and gases across the sea surface; bulk formulae |
| `waves` | Wind waves |

## 2. Dependencies

| Source | Used for |
|---|---|
| `meteorology` | The atmospheric boundary layer |
| `oceanography` | The ocean mixed layer, seawater properties |
| physis | `continuum`, `thermodynamics` |

## 3. Applications

See [`applications/air_sea/`](../../../applications/air_sea/).
