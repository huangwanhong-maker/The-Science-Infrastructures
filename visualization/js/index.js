// The public interface of the shared visualization layer's provisional JavaScript.
// Code outside it imports from this file only.

export * from './contour.js'; // isolines and extrema of gridded fields
export * from './flow.js';    // particles that make a flow on the sphere visible
export * from './tiles.js';   // tiles of an equirectangular map service
export * from './units.js';   // the units in which quantities are shown, and numbers for reading
export * from './isosurface.js'; // a surface at a level of a volume, by surface nets
export * from './complex.js';    // cell complexes: realised in space, unfolded into nets, as Hasse diagrams; words and facets read
