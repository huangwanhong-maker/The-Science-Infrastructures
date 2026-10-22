# air-sea

*Status: planned. Bridge, height 3.* The exchange of momentum, heat, water and gases at the sea surface, where the atmosphere and the ocean act on each other.

## 1. Scope

| Module | Scope |
|---|---|
| `fluxes` | Fluxes of momentum, heat, freshwater and gases across the sea surface; bulk formulae |
| `waves` | Wind waves: their growth under the wind |

## 2. Dependencies

| Source | Used for |
|---|---|
| meteorology | The atmospheric boundary layer |
| oceanography | The ocean mixed layer, seawater properties, the propagation of waves |
| physis | `continuum`, `thermodynamics` |

## 3. Used By

climate; the air–sea exchange application (planned).
