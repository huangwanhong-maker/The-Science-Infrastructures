// Directions in the local horizon: azimuths, in degrees clockwise from north,
// named by the points of the compass, and averaged.
//
// Provisional JavaScript, standing in for the C++ planned for this domain (ADR 0014).

const POINTS = ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW'];

// The nearest of the sixteen points of the compass.
export const compassPoint = (degrees) => POINTS[Math.round((((degrees % 360) + 360) % 360) / 22.5) % 16];

// Mean of directions (degrees) with weights, through their unit vectors, since
// the mean of 350 and 10 degrees is 0 and not 180.
export function meanDirection(directions, weights) {
  let x = 0, y = 0;
  directions.forEach((d, i) => {
    const w = weights ? weights[i] : 1;
    x += w * Math.sin((d * Math.PI) / 180);
    y += w * Math.cos((d * Math.PI) / 180);
  });
  return x === 0 && y === 0 ? NaN : ((Math.atan2(x, y) * 180) / Math.PI + 360) % 360;
}
