// The forecast at a place from JavaScript (Node 18 or later): GFS's temperature at
// 2 m every six hours, and the grid over Japan at the first time as a Float32Array.
//
//   node api/examples/js/forecast.mjs [lon lat]          (the service: $SCIENCE_API_URL or http://127.0.0.1:8800)

import { Client } from '../../js/science.mjs';

const [lon = 139.69, lat = 35.69] = process.argv.slice(2).map(Number);
const api = new Client();
const gfs = await api.field('gfs-forecast');
const { times, values } = await gfs.series('tmp2m', lon, lat);
console.log(`GFS at ${lon}, ${lat} (UTC)`);
times.forEach((t, i) => console.log(`  ${t.toISOString().slice(0, 16).replace('T', ' ')}  ${(values[i] - 273.15).toFixed(1)} °C`));

const { values: grid, grid: g, time } = await gfs.grid('tmp2m', { time: 0, bbox: [122, 24, 150, 46] });
let warmest = -Infinity, at = null;
grid.forEach((v, k) => { if (v > warmest) { warmest = v; at = [g.lon0 + (k % g.nlon) * g.dlon, g.lat0 - Math.floor(k / g.nlon) * g.dlat]; } });
console.log(`\nOver Japan at ${time}: ${g.nlon} × ${g.nlat} points; the warmest ${(warmest - 273.15).toFixed(1)} °C at ${at}`);
