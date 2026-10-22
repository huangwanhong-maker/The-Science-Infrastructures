# geography

*Status: draft. Field, height 2.* Geography in all its branches, as the work requires ([ADR 0010](../../development/decisions/0010-domain-organised-fields-and-applications.md)).

## 1. Scope

| Branch | Scope | Status |
|---|---|---|
| Cartography and GIS | Map design; map projections; spatial data models (vector and raster); polygons on rasters and areas on the sphere; joins of data to places; zonal statistics (fields summarised over regions on the sphere); coordinate reference systems in practice; spatial analysis and spatial statistics | draft (JavaScript: `projections.js`, `raster.js`, `join.js`, `zonal.js`) |
| Physical geography | Geomorphology and terrain; rivers, drainage basins and the hydrological cycle on land; soils; biogeography | draft for terrain and drainage (JavaScript: `hydrology.js`: slope, depressions filled by Priority-Flood, D8 flow on the sphere, upstream area, watersheds and paths); the rest planned |
| Remote sensing | Sensors and observing geometry; radiometric and geometric correction; classification of land cover | planned |
| Human geography | Population, settlement, land use, and human–environment interaction | draft for population (GHS-POP on a grid, `data/products/population_grid.py`, summed over regions by `zonal.js`; exposure to hazards in the Earth Atlas); the rest planned |
| Other branches | Added as the work requires, for example historical geography | planned |

## 2. Dependencies

- **Domains:** geodesy, for reference systems and the figure of the Earth.
- **Foundations:** physis (`geometry` for projections and geodesics, `radiation` for remote sensing); morphe (the topology of spatial data: planar subdivisions and complexes); calx.

## 3. Code

`js/index.js` exports `projections.js`, `raster.js`, `join.js`, `zonal.js` and `hydrology.js`. The checks (projections against Šavrič et al. 2018, rasters, joins of the World Bank's economies to Natural Earth's countries, the land fraction, water in shaded relief; zonal statistics against closed forms on the sphere and GHS-POP against Japan's census; drainage of synthetic surfaces and of five great basins on ETOPO 2022) are run by:

```
node domains/geography/tests/run-tests.mjs
```

## 4. Notes

- [Equal-area maps and data per place](theory/equal-area-maps.md)
- [Zonal statistics on the sphere](theory/zonal-statistics.md)
- [Drainage from elevation](theory/drainage-from-elevation.md)

## 5. Candidate Bridges

- `coasts`, joining geography and oceanography;
- `hydrology`, joining geography and meteorology: the water cycle on land, from rain to runoff, on the drainage network that geography's `hydrology.js` draws from elevation.

## 6. Used By

navigation and climate (planned); the Earth Atlas.
