# geography

*Status: planned.* Field. Geography is taken in all its branches, as the work requires ([ADR 0010](../../../development/decisions/0010-domain-organised-fields-and-applications.md)).

## 1. Scope

| Branch | Scope |
|---|---|
| Physical geography | Geomorphology and terrain; rivers, drainage basins and the hydrological cycle on land; soils; biogeography; climatic zones |
| Cartography and GIS | Map design; spatial data models (vector and raster); coordinate reference systems in practice; terrain models; spatial analysis and spatial statistics |
| Remote sensing | Sensors and observing geometry; radiometric and geometric correction; classification of land cover |
| Human geography | Population, settlement, land use, and human–environment interaction |
| Other branches | Added as the work requires, for example historical geography |

## 2. Dependencies

- **Within this repository:** `geodesy`, for reference systems and the figure of the Earth.
- **physis:**
  - `geometry`, for map projections and geodesics;
  - `radiation`, for remote sensing.
- **morphe:** the topology of spatial data, such as planar subdivisions and complexes.
- **calx:** numerical methods.

## 3. Candidate Bridges

These bridges will be added to `core/domains/` when work on them begins:

- `coasts`, joining geography and oceanography;
- `hydrology`, joining geography and meteorology.

## 4. Applications

See [`applications/geography/`](../../../applications/geography/).
