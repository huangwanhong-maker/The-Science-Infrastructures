# meteorology

*Status: planned.* Field. Here "meteorology" is taken in the broad sense of atmospheric science ([ADR 0007](../../../development/decisions/0007-atmosphere-and-ocean-as-coupled-system.md), [ADR 0010](../../../development/decisions/0010-domain-organised-fields-and-applications.md)).

## 1. Scope

| Module | Scope |
|---|---|
| `dynamics` | Atmospheric dynamics |
| `moist` | Moist thermodynamics, clouds and precipitation |
| `radiation_budget` | The radiation budget of the atmosphere |
| `boundary_layer` | The atmospheric boundary layer |
| `chemistry` | Atmospheric chemistry |
| `prediction` | Weather analysis and prediction |

## 2. Dependencies

- `geodesy`, for reference surfaces.
- physis: `continuum`, `thermodynamics`, `radiation` and `chemistry`.
- calx.

## 3. Applications

See [`applications/meteorology/`](../../../applications/meteorology/).
