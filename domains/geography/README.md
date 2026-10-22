# geography

*Status: draft. Field, height 2.* Geography in all its branches, as the work requires ([ADR 0010](../../development/decisions/0010-domain-organised-fields-and-applications.md)).

## 1. Scope

| Branch | Scope | Status |
|---|---|---|
| Cartography and GIS | Map design; map projections; spatial data models (vector and raster); polygons on rasters and areas on the sphere; joins of data to places; coordinate reference systems in practice; spatial analysis and spatial statistics | draft (JavaScript: `projections.js`, `raster.js`, `join.js`) |
| Physical geography | Geomorphology and terrain; rivers, drainage basins and the hydrological cycle on land; soils; biogeography | planned |
| Remote sensing | Sensors and observing geometry; radiometric and geometric correction; classification of land cover | planned |
| Human geography | Population, settlement, land use, and human–environment interaction | planned |
| Other branches | Added as the work requires, for example historical geography | planned |

## 2. Dependencies

- **Domains:** geodesy, for reference systems and the figure of the Earth.
- **Foundations:** physis (`geometry` for projections and geodesics, `radiation` for remote sensing); morphe (the topology of spatial data: planar subdivisions and complexes); calx.

## 3. Code

`js/index.js` exports `projections.js`, `raster.js` and `join.js`. The checks (projections against Šavrič et al. 2018, rasters, joins of the World Bank's economies to Natural Earth's countries, the land fraction, water in shaded relief) are run by:

```
node domains/geography/tests/run-tests.mjs
```

## 4. Notes

- [Equal-area maps and data per place](theory/equal-area-maps.md)

## 5. Candidate Bridges

- `coasts`, joining geography and oceanography;
- `hydrology`, joining geography and meteorology.

## 6. Used By

navigation and climate (planned); the Earth Atlas.
