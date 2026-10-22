// The public interface of the geophysics domain's provisional JavaScript.
// Code outside it imports from this file only.

export * from './plates.js';  // the motion of rigid tectonic plates (NNR-MORVEL56, ITRF2020)
export * from './moment.js';  // earthquake mechanisms from moment tensors
export * from './crust.js';   // the crust beneath a place (CRUST1.0)
export * from './section.js'; // sections through the Earth along great circles
