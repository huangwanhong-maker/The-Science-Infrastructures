// The 2D GIS of the visualization layer (ADR 0069): maps of the sphere in a projection, drawn in SVG.
//
// - projection.js: the geography domain's projections in the plane of a drawing; outlines, graticules,
//   polygons, polylines and great-circle arcs as paths (pure);
// - view.js: zoom about a point, fitting, centring (pure);
// - labels.js: names placed by priority without overlap, and the data library's naming policy (pure);
// - widget.js: an interactive map: the base, an overlay for an application's layers, zoom, pan, picking
//   and the readout of the projection, the place and the zoom (DOM);
// - figure.js: a map composed as a figure for publication, SVG with its fonts or PNG (DOM).
//
// The pure parts are checked in Node (visualization/tests); the city-scale GIS in three dimensions is
// ../gis.

export * from './projection.js';
export * from './view.js';
export * from './labels.js';
export * from './widget.js';
export * from './figure.js';
