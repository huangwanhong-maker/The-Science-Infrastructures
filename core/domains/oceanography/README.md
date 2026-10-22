# oceanography

*Status: planned.* Field.

## 1. Scope

| Module | Scope |
|---|---|
| `physical` | Circulation and waves |
| `seawater` | Seawater thermodynamics (TEOS-10) |
| `biogeochemistry` | Chemical and biological oceanography |
| `seafloor` | Bathymetry and the seafloor |

## 2. Dependencies

| Dependency | Used for |
|---|---|
| physis | `continuum`, `thermodynamics`, `chemistry` |
| `geodesy` | The geoid and other reference surfaces |
| calx | Numerical methods and device kernels |

## 3. Applications

See [`applications/oceanography/`](../../../applications/oceanography/).
