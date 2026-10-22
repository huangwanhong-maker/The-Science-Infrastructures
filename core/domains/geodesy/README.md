# geodesy

*Status: planned.* Field. Geodesy supports the other fields: geography, oceanography and meteorology all rest on its reference systems.

## 1. Scope

| Module | Scope |
|---|---|
| `figure` | The figure of the Earth: ellipsoids and the geoid |
| `gravity` | The gravity field |
| `rotation` | The rotation of the Earth |
| `reference_systems` | Coordinate reference systems |

## 2. Dependencies

- `core/metrology`: time scales, frames and Earth orientation.
- physis: `gravitation` and `geometry`.
- calx.

## 3. Applications

See [`applications/geodesy/`](../../../applications/geodesy/).
