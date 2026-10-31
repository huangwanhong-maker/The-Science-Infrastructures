// The public interface of the shared visualization layer's provisional JavaScript.
// Code outside it imports from this file only.

export * from './contour.js'; // isolines and extrema of gridded fields
export * from './flow.js';    // particles that make a flow on the sphere visible
export * from './tiles.js';   // tiles of an equirectangular map service
export * from './units.js';   // the units in which quantities are shown, and numbers for reading
export * from './isosurface.js'; // a surface at a level of a volume, by surface nets
export * from './complex.js';    // cell complexes: realised in space, unfolded into nets, as Hasse diagrams; words and facets read
export * from './stars.js';   // the colour of a star from its colour index
export * as globe from './globe/index.js'; // virtual globes: tiles, figures, terrain, level of detail, imagery, controls, a three.js renderer
export * as gis from './gis/index.js';     // a city-scale GIS in three dimensions: the local frame, rings and ribbons, colours, a three.js scene
// The scene model, the styles and the renderers (ADR 0005, 0053), each under its name.
export * as scene from '../scene/index.js';         // samplings, a field's answers, scales, visual expressions and their rules, legends, volumes
export * as styles from '../styles/index.js';       // colour in OKLab, the colour maps, colour-vision deficiencies
export * as renderers from '../renderers/index.js'; // a renderer's arithmetic: colouring, drapes, coloured curves
export { fieldMarks, canvasOf } from '../renderers/three/index.js'; // the marks in three.js, which is passed in
