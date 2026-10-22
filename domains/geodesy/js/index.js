// The public interface of the geodesy domain's provisional JavaScript.
// Code outside it imports from this file only.

export * from './geodesy.js'; // places on the WGS 84 ellipsoid: geodetic and geocentric latitude, distances
export * from './sphere.js';  // the Earth as a sphere of mean radius: positions, directions, distances, motion
export * from './compass.js'; // azimuths named by the points of the compass, and their mean
