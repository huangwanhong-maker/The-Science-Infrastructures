// The public interface of the meteorology domain's provisional JavaScript.
// Code outside it imports from this file only.

export * from './wind.js';         // wind speed and the direction from which it blows
export * from './barotropic.js';   // the barotropic vorticity equation, and the Earth as Williamson et al. give it
export * from './shallowwater.js'; // the shallow-water equations on the sphere
export * from './cyclones.js';     // tropical cyclones found and tracked in sea-level pressure
export * from './fireweather.js';  // the Canadian Fire Weather Index System, and daily weather from a forecast
export * from './heat.js';         // wet-bulb temperature, the heat index, daily maxima
export * from './settling.js';     // the standard atmosphere, Stokes settling, falling parcels in the wind
