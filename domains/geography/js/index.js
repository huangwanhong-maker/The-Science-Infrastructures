// The public interface of the geography domain's provisional JavaScript.
// Code outside it imports from this file only.

export * from './projections.js'; // map projections: Equal Earth and equirectangular
export * from './raster.js';      // polygons drawn onto rasters, and areas on the sphere
export * from './join.js';        // data about places joined to the geometry of places
export * from './zonal.js';       // fields summarised over regions: areal means, sums, extremes
export * from './hydrology.js';   // terrain and water: slope, depressions, flow, watersheds
