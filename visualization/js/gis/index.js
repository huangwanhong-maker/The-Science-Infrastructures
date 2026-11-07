// A city-scale GIS in three dimensions: the local frame on WGS 84 (east, north, up from an origin), the planar
// geometry of rings, polylines and ribbons, colours of features from the validated maps, and a three.js scene
// of extruded buildings, streets, areas, points and bars with map controls, picking and the Sun's light.

export * from './frame.js';     // geodetic ↔ local east, north, up, exactly on WGS 84
export * from './geometry.js';  // rings, containment, centroids, extents, polylines and ribbons
export * from './style.js';     // colours of categories, magnitudes and signed values
export * from './renderer-three.js'; // the scene: solids, flats, ribbons, points; picking, hover, the Sun, flights
export * from './prisms.js';     // columns on cells, drawn by one mesh and a table of heights and colours
export * from './world.js';      // any place on Earth: tiles of buildings and ways, asked for as the camera comes near
export * from './earth.js';      // one view from the whole globe to a street: the Earth beneath a city, columns and rasters over it, the sky
