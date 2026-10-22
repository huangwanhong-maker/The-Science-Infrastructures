// The public interface of the meteorology domain's provisional JavaScript.
// Code outside it imports from this file only.

export * from './wind.js';         // wind speed and the direction from which it blows
export * from './barotropic.js';   // the barotropic vorticity equation, and the Earth as Williamson et al. give it
export * from './shallowwater.js'; // the shallow-water equations on the sphere
export * from './cyclones.js';     // tropical cyclones found and tracked in sea-level pressure
