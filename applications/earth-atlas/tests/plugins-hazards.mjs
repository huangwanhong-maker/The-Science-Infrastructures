// Checks of the hazard plug-ins (fire weather, heat, heavy rain, tsunami travel
// times, aftershocks, volcanic ash), run without a page from constructed
// products whose answers are known, and, where the cache holds it, the tsunami
// of 1960 over the real relief. Called by run-tests.mjs with its harness.

import { existsSync, readFileSync } from 'node:fs';
import { resolve } from 'node:path';

import { FieldProduct } from '../src/data/fields.js';
import { fireWeatherDay, STANDARD_START, wetBulbStull, heatIndex, stokesSettling } from '../../../domains/meteorology/js/index.js';
import { aftershockProbability, subsurfaceRuptureLengthKm } from '../../../domains/geophysics/js/index.js';
import { longWaveSpeed } from '../../../domains/oceanography/js/index.js';
import { distanceM } from '../../../domains/geodesy/js/index.js';

const HOUR = 3600000, DAY = 24 * HOUR;

export async function checks({ check, assert, skip, PRODUCTS, api }) {
  const grid = { lon0: -180, dlon: 2.5, nlon: 144, lat0: 90, dlat: 2.5, nlat: 73, periodic: true };
  const cells = (fn) => {
    const a = new Float32Array(grid.nlon * grid.nlat);
    for (let r = 0; r < grid.nlat; r++) for (let c = 0; c < grid.nlon; c++) a[r * grid.nlon + c] = fn(grid.lon0 + c * grid.dlon, grid.lat0 - r * grid.dlat);
    return a;
  };
  const cell = (g, lon, lat) => Math.round((g.lat0 - lat) / g.dlat) * g.nlon + ((Math.round((lon - g.lon0) / g.dlon) % g.nlon) + g.nlon) % g.nlon;
  const input = (product, names, level = null) => ({ product, names, variables: names.map((n) => product.meta.variables.find((v) => v.name === n)), meta: product.meta, model: product.meta.title, level });
  const context = ({ time = 0, pinned = null, fields = [], events = {}, geography = null } = {}) => ({
    time, pinned, progress() {},
    fields: () => fields.map((p) => p.meta), field: async (name) => fields.find((p) => p.meta.name === name),
    list: (kind) => (kind === 'events' ? Object.values(events).map((s) => ({ name: s.name, kind: s.kind, title: s.title })) : []),
    events: async (name) => events[name], geography: async () => geography,
  });
  // A forecast every 6 h for 6 days, from 18 UTC, as GFS's products are.
  const forecast = (t0, variables) => {
    const times = [...Array(24).keys()].map((k) => t0 + 6 * k * HOUR);
    return api.makeProduct({ name: 'test-forecast', title: 'Test forecast', grid, times,
      variables: variables.map(([name, long_name, standard_name, units, fn]) => ({ name, long_name, standard_name, units, slices: times.map((t) => cells((lon, lat) => fn(lon, lat, t))) })) });
  };
  // Hours of local mean solar time at a longitude and instant.
  const solarHour = (lon, t) => ((((t + (lon / 15) * HOUR) % DAY) + DAY) % DAY) / HOUR;

  // ---- fire weather ----
  {
    // The weather of the first day of Van Wagner and Pickett's test data (17 °C, 42 %,
    // 25 km/h, no rain) at local noon only, and cold and still at the other times: if the
    // plug-in takes each place's noon, its first day is the published one everywhere north of 30° N.
    const t0 = Date.UTC(2026, 3, 11, 18);
    const noonOnly = (lon, t, atNoon, otherwise) => (Math.abs(solarHour(lon, t) - 12) <= 3 ? atNoon : otherwise);
    const dry = forecast(t0, [
      ['t2', 'air temperature at 2m', 'air_temperature', 'K', (lon, lat, t) => noonOnly(lon, t, 290.15, 263.15)],
      ['rh', 'relative humidity at 2m', 'relative_humidity', '%', (lon, lat, t) => noonOnly(lon, t, 42, 90)],
      ['u', 'eastward wind at 10m', 'eastward_wind', 'm s-1', (lon, lat, t) => noonOnly(lon, t, 25 / 3.6, 0)],
      ['v', 'northward wind at 10m', 'northward_wind', 'm s-1', () => 0],
      ['p', 'rainfall rate', 'precipitation_flux', 'kg m-2 s-1', () => 0],
    ]);
    const fire = (await import('../plugins/fire-weather.js')).create(api);
    const ins = { temperature: input(dry, ['t2'], '2 m'), humidity: input(dry, ['rh'], '2 m'), wind: input(dry, ['u', 'v'], '10 m'), rain: input(dry, ['p']) };
    const [result] = await fire.run(ins, { ffmc: 85, dmc: 6, dc: 15 }, context());
    const p = result.product, published = { fire_weather_index: 10.1, initial_spread_index: 10.9, buildup_index: 8.5, fine_fuel_moisture_code: 87.7, duff_moisture_code: 8.5, drought_code: 19.0 };
    let worst = 0;
    for (const lon of [-177.5, -60, 0, 45, 92.5, 180]) for (const lat of [32.5, 47.5, 70]) {
      for (const [name, v] of Object.entries(published)) worst = Math.max(worst, Math.abs(p.slice(name, 0)[cell(grid, lon, lat)] - v));
    }
    check('plug-in fire weather: the first day at local noon is the standard test day (FWI 10.1) at every longitude', worst, 0.05);
    assert('plug-in fire weather: whole days 13-16 April, at 12 UTC', p.times.length === 4 && api.iso(p.times[0]) === '2026-04-13T12:00:00Z', p.times.map(api.iso).join(' '));
    // Rain of 2.4 mm a day, from a constant rate: the first day is the domain's with 2.4 mm.
    const wet = forecast(t0, [
      ['t2', 'air temperature at 2m', 'air_temperature', 'K', () => 293.15], ['rh', 'relative humidity at 2m', 'relative_humidity', '%', () => 21],
      ['u', 'eastward wind at 10m', 'eastward_wind', 'm s-1', () => 25 / 3.6], ['v', 'northward wind at 10m', 'northward_wind', 'm s-1', () => 0],
      ['p', 'rainfall rate', 'precipitation_flux', 'kg m-2 s-1', () => 2.4 / 86400],
    ]);
    const [rained] = await fire.run({ temperature: input(wet, ['t2'], '2 m'), humidity: input(wet, ['rh'], '2 m'), wind: input(wet, ['u', 'v'], '10 m'), rain: input(wet, ['p']) },
      { ffmc: 85, dmc: 6, dc: 15 }, context());
    const expected = fireWeatherDay(STANDARD_START, { T: 20, H: 21, W: 25, rain: 2.4 }, 4, 47.5), i = cell(grid, 10, 47.5);
    check('plug-in fire weather: 24 h of rain from the rate (2.4 mm a day): DMC and FFMC as the domain gives them',
      Math.abs(rained.product.slice('duff_moisture_code', 0)[i] - expected.dmc) + Math.abs(rained.product.slice('fine_fuel_moisture_code', 0)[i] - expected.ffmc), 1e-4);
  }

  // ---- heat ----
  {
    // 30 °C within three hours of local noon and 20 °C otherwise, at 50 %; a normal daily maximum of 25 °C.
    const t0 = Date.UTC(2026, 6, 10, 18);
    const hot = forecast(t0, [
      ['t2', 'air temperature at 2m', 'air_temperature', 'K', (lon, lat, t) => (Math.abs(solarHour(lon, t) - 12) <= 3 ? 303.15 : 293.15)],
      ['rh', 'relative humidity at 2m', 'relative_humidity', '%', () => 50],
    ]);
    const months = [...Array(12).keys()].map((m) => Date.UTC(2005, m, 15));
    const normals = api.makeProduct({ name: 'test-normals', title: 'Test normals', grid, times: months, variables: [
      { name: 'tmax', long_name: 'Daily maximum air temperature at 2 m', standard_name: 'air_temperature', units: 'degree_C', cell_methods: 'time: maximum within days time: mean over days time: mean over years', slices: months.map(() => cells(() => 25)) }] });
    normals.meta.climatology = normals.climatology = { period: '1991-01-01/2020-12-31', year: 2005, cycle: 'year' };
    const heat = (await import('../plugins/heat.js')).create(api);
    const [result] = await heat.run({ temperature: input(hot, ['t2'], '2 m'), humidity: input(hot, ['rh'], '2 m') }, {}, context({ fields: [normals] }));
    const p = result.product, i = cell(grid, 72.5, 20);
    check('plug-in heat: the daily maximum wet-bulb temperature is Stull’s at the day’s 30 °C', Math.abs(p.slice('wet_bulb_max', 1)[i] - 273.15 - wetBulbStull(30, 50)), 1e-4, '°C');
    check('plug-in heat: the heat index likewise (NWS, 31.4 °C)', Math.abs(p.slice('heat_index_max', 1)[i] - heatIndex(303.15, 50)), 1e-4, 'K');
    let off = 0;
    for (const lon of [-180, -97.5, 0, 37.5, 140]) off = Math.max(off, Math.abs(p.slice('tmax_excess', 2)[cell(grid, lon, -30)] - 5));
    check('plug-in heat: 5 °C above a normal maximum of 25 °C, at every longitude', off, 1e-4, 'K');
  }

  // ---- heavy rain ----
  {
    const t0 = Date.UTC(2026, 9, 20, 18), rate = 1e-4;
    const rainy = forecast(t0, [['p', 'rainfall rate', 'precipitation_flux', 'kg m-2 s-1', () => rate]]);
    const months = [...Array(12).keys()].map((m) => Date.UTC(2005, m, 15));
    const normals = api.makeProduct({ name: 'test-ppt', title: 'Test normals', grid, times: months, variables: [
      { name: 'ppt', long_name: 'Precipitation in the month', standard_name: 'precipitation_amount', units: 'kg m-2', slices: months.map(() => cells((lon) => (lon < 0 ? 100 : 0.5))) }] });
    normals.meta.climatology = normals.climatology = { period: '1991-01-01/2020-12-31', year: 2005, cycle: 'year' };
    const rain = (await import('../plugins/heavy-rain.js')).create(api);
    const [result] = await rain.run({ rain: input(rainy, ['p']) }, { least: 1 }, context({ fields: [normals] }));
    const p = result.product, last = p.times.length - 1, w = cell(grid, -60, 10), e = cell(grid, 60, 10);
    check('plug-in heavy rain: 0.36 mm/h for 24 h is 8.64 mm', Math.abs(p.slice('rain_24h', 0)[w] - 8.64), 1e-4, 'mm');
    check('plug-in heavy rain: and for 72 h 25.92 mm', Math.abs(p.slice('rain_72h', last)[w] - 25.92), 1e-4, 'mm');
    check('plug-in heavy rain: as fractions of a normal month of 100 mm', Math.abs(p.slice('rain_24h_ratio', 0)[w] - 0.0864) + Math.abs(p.slice('rain_72h_ratio', last)[w] - 0.2592), 1e-6);
    assert('plug-in heavy rain: no 72 h sum before the forecast covers 72 h, no fraction where the normal is below 1 mm',
      Number.isNaN(p.slice('rain_72h', 0)[w]) && Number.isNaN(p.slice('rain_24h_ratio', 0)[e]) && api.iso(p.times[0]) === '2026-10-21T12:00:00Z', api.iso(p.times[0]));
  }

  // ---- tsunami travel times ----
  {
    // An ocean of 4000 m on a 1° grid with a continent from 30° W to 10° E between 60° S and
    // 60° N: from the chosen place the wave reaches a point at 3000 km in its great-circle time.
    const g = { lon0: -180, dlon: 1, nlon: 360, lat0: 90, dlat: 1, nlat: 181, periodic: true };
    const z = new Float32Array(g.nlon * g.nlat);
    for (let r = 0; r < g.nlat; r++) for (let c = 0; c < g.nlon; c++) { const lon = -180 + c, lat = 90 - r; z[r * g.nlon + c] = lon >= -30 && lon <= 10 && Math.abs(lat) <= 60 ? 500 : -4000; }
    const relief = api.makeProduct({ name: 'test-relief', title: 'Test relief', grid: g, times: [], variables: [{ name: 'z', long_name: 'Elevation', standard_name: 'altitude', units: 'm', slices: [z] }] });
    const tsunami = (await import('../plugins/tsunami.js')).create(api);
    const [field, sketch] = await tsunami.run({ relief: input(relief, ['z']) }, { spacing: 60, stencil: 32, magnitude: 6.5 }, context({ pinned: [-150, 0] }));
    const h = field.product.slice('tsunami_travel_time', 0), target = cell(g, -123, 0);
    const exact = distanceM(-150, 0, -123, 0) / longWaveSpeed(4000) / 3600;
    check(`plug-in tsunami: a point 3000 km away at its great-circle time over 4000 m (${exact.toFixed(2)} h)`, h[target] / exact - 1, 0.03);
    assert('plug-in tsunami: hours, with hourly contours, the source marked; land has no time',
      field.product.meta.variables[0].units === 'h' && field.product.meta.variables[0].contours.step === 1 && field.show.lines === 'tsunami_travel_time'
      && sketch.markers[0].lon === -150 && Number.isNaN(h[cell(g, -10, 0)]));
    // Behind the continent, the wave comes round it: later than straight through.
    const behind = h[cell(g, 20, 0)], straight = distanceM(-150, 0, 20, 0) / longWaveSpeed(4000) / 3600;
    assert('plug-in tsunami: behind a continent the wave comes round it, later than the straight path', behind > straight, `${behind.toFixed(2)} h against ${straight.toFixed(2)} h`);
    // With no place chosen, the latest shallow earthquake of M ≥ 6.5 is the source.
    const quakes = { name: 'test-quakes', kind: 'earthquake', title: 'Test earthquakes', fields: ['time_ms', 'lon', 'lat', 'depth_km', 'mag', 'mag_type', 'place', 'id'], events: [
      [Date.UTC(2026, 9, 1), 160, -20, 10, 7.1, 'mww', 'early', 'a'], [Date.UTC(2026, 9, 5), 140, 10, 120, 7.4, 'mww', 'deep', 'b'],
      [Date.UTC(2026, 9, 8), -100, 0, 15, 6.6, 'mww', 'latest shallow', 'c'], [Date.UTC(2026, 9, 9), -80, -10, 10, 6.0, 'mww', 'small', 'd']] };
    const [, chosen] = await tsunami.run({ relief: input(relief, ['z']) }, { spacing: 60, stencil: 16, magnitude: 6.5 }, context({ events: { 'test-quakes': quakes } }));
    assert('plug-in tsunami: without a place, from the latest shallow earthquake of M ≥ 6.5', chosen.markers[0].lon === -100 && /latest shallow/.test(chosen.title), chosen.title);
    // The tsunami of 22 May 1960 from the Chilean earthquake (Mw 9.5, 19:11 UTC; USGS epicentre
    // 38.24° S, 73.05° W, inland, so the source moves to the nearest sea) reached Hilo about
    // 15 h later and the Sanriku coast of Japan about 22 h later (Atwater et al. 2005, Surviving
    // a Tsunami: Lessons from Chile, Hawaii, and Japan, USGS Circular 1187).
    const meta = resolve(PRODUCTS, 'fields-etopo-2022-medium.json');
    if (existsSync(meta)) {
      const m = JSON.parse(readFileSync(meta, 'utf8')), bin = readFileSync(resolve(PRODUCTS, m.file));
      const etopo = new FieldProduct(m, bin.buffer.slice(bin.byteOffset, bin.byteOffset + bin.byteLength));
      const begun = performance.now();
      const [chile, , note] = await tsunami.run({ relief: input(etopo, ['z']) }, { spacing: 20, stencil: 32, magnitude: 6.5 }, context({ pinned: [-73.05, -38.24] }));
      const seconds = (performance.now() - begun) / 1000;
      const gr = chile.product.grid, t = chile.product.slice('tsunami_travel_time', 0);
      const nearest = (lon, lat, km) => { let best = Infinity; for (let r = 0; r < gr.nlat; r++) for (let c = 0; c < gr.nlon; c++) { const x = t[r * gr.nlon + c]; if (x === x && x < best && distanceM(lon, lat, gr.lon0 + c * gr.dlon, gr.lat0 - r * gr.dlat) < km * 1000) best = x; } return best; };
      const hilo = nearest(-155.08, 19.72, 40), sanriku = nearest(141.9, 39.0, 60);
      assert(`plug-in tsunami: 1960 from Valdivia over ETOPO 2022 at 20′: Hilo in ${hilo.toFixed(1)} h (about 15), Sanriku in ${sanriku.toFixed(1)} h (about 22); ${seconds.toFixed(1)} s`,
        Math.abs(hilo - 15) < 1 && Math.abs(sanriku - 22) < 1.5, note.text);
    } else {
      skip('plug-in tsunami: fields-etopo-2022-medium missing (the tsunami of 1960)');
    }
  }

  // ---- aftershocks ----
  {
    const tm = Date.UTC(2026, 9, 9, 18), now = tm + 3 * DAY;
    const set = { name: 'test-month', kind: 'earthquake', title: 'Test earthquakes', selection: { min_magnitude: 4.5 }, window: { start: '2026-09-21T00:00:00Z', end: '2026-10-21T00:00:00Z' },
      fields: ['time_ms', 'lon', 'lat', 'depth_km', 'mag', 'mag_type', 'place', 'id'], events: [
        [Date.UTC(2026, 8, 25), 120, 30, 10, 7.5, 'mww', 'far and larger', 'x'], [tm, 0, 0, 10, 7.0, 'mww', 'mainshock', 'm'],
        [tm + HOUR, 0.2, 0.1, 10, 5.1, 'mb', 'near', 'a1'], [tm + 2 * DAY, -0.3, 0.2, 10, 4.6, 'mb', 'near', 'a2'], [tm + DAY, 3, 3, 10, 4.8, 'mb', 'outside the zone', 'o']] };
    const after = (await import('../plugins/aftershocks.js')).create(api);
    const params = { radius: 300, a: -1.67, b: 0.91, p: 1.08, c: 0.05 };
    const [table, zone, observed] = await after.run({}, params, context({ time: now, pinned: [1, 0.5], events: { 'test-month': set } }));
    const want = (M, days) => aftershockProbability(3, 3 + days, M, 7.0);
    assert('plug-in aftershocks: the mainshock is the largest within 300 km of the place, not the largest of the month', /M 7 /.test(table.title), table.title);
    assert('plug-in aftershocks: the chances over the next day and week from the instant shown, as the domain gives them',
      table.rows[0][1] === `${Math.round(100 * want(5, 1))} %` && table.rows[1][2] === `${Math.round(100 * want(6, 7))} %`, `${table.rows.map((r) => r.join(' ')).join('; ')}`);
    const L = subsurfaceRuptureLengthKm(7.0), ring = zone.polylines[0].coords;
    let worst = 0;
    for (let k = 0; k < ring.length; k += 2) worst = Math.max(worst, Math.abs(distanceM(0, 0, ring[k], ring[k + 1]) / 1000 - L));
    check(`plug-in aftershocks: the zone is a circle of the rupture length, ${L.toFixed(1)} km`, worst, 0.01, 'km');
    assert('plug-in aftershocks: two of the catalogue’s earthquakes fall in it since the mainshock', /Observed: 2 earthquakes/.test(observed.text), observed.text);
    const [far] = await after.run({}, params, context({ time: now, events: { 'test-month': set } }));
    assert('plug-in aftershocks: without a place, the largest of the month', /M 7.5/.test(far.title), far.title);
  }

  // ---- volcanic ash ----
  {
    const t0 = Date.UTC(2026, 9, 20, 0), times = [...Array(17).keys()].map((k) => t0 + 6 * k * HOUR);
    const winds = (u) => api.makeProduct({ name: 'test-winds', title: 'Test winds', grid, times, variables: ['10m', '850', '500', '250'].flatMap((l) => [
      { name: `u${l}`, long_name: `eastward wind ${l}`, standard_name: 'eastward_wind', units: 'm s-1', slices: times.map(() => cells(() => u)) },
      { name: `v${l}`, long_name: `northward wind ${l}`, standard_name: 'northward_wind', units: 'm s-1', slices: times.map(() => cells(() => 0)) }]) });
    const ins = (p) => ({ wind10: input(p, ['u10m', 'v10m']), wind850: input(p, ['u850', 'v850']), wind500: input(p, ['u500', 'v500']), wind250: input(p, ['u250', 'v250']) });
    const ash = (await import('../plugins/volcanic-ash.js')).create(api);
    const params = { top: 3, count: 3, diameter: 30, density: 2500, duration: 0, hours: 96, step: 0.25 };
    // Still air: each parcel lands at the vent, after the time to fall through ∫ dz / w(z).
    const [still, stillTable] = await ash.run(ins(winds(0)), params, context({ time: t0, pinned: [20, 0] }));
    const fall = (z) => { let s = 0; for (let h = z; h > 0; h -= 1) s += 1 / stokesSettling(30e-6, 2500, h - 0.5).w; return s; };
    const landing = still.markers[1], lowest = still.moving[0];
    check('plug-in volcanic ash: in still air a parcel lands at the vent', Math.hypot(landing.lon - 20, landing.lat), 1e-9, '°');
    check('plug-in volcanic ash: after falling 1 km at its Stokes speed (relative)', ((lowest.times[lowest.times.length - 1] - t0) / 1000 - fall(1000)) / fall(1000), 1e-3);
    assert('plug-in volcanic ash: one path, one moving point and one landing per release height', still.polylines.length === 3 && still.moving.length === 3 && still.markers.length === 4 && stillTable.rows.length === 3);
    // A uniform wind of 10 m/s along the equator: the lowest parcel lands U t from the vent.
    const [windy] = await ash.run(ins(winds(10)), params, context({ time: t0, pinned: [20, 0] }));
    const m = windy.moving[0], tl = (m.times[m.times.length - 1] - t0) / 1000;
    check('plug-in volcanic ash: in a uniform wind the displacement is U t', distanceM(20, 0, windy.markers[1].lon, windy.markers[1].lat) / (10 * tl) - 1, 1e-6);
    // Without a place, the latest volcano of GDACS's events; with none, a note.
    const gdacs = { name: 'test-gdacs', kind: 'disaster', title: 'Test disasters', fields: ['time_ms', 'end_ms', 'lon', 'lat', 'type', 'alert', 'name', 'severity', 'country', 'id', 'report'],
      events: [[t0 - 5 * DAY, t0, 110.44, -7.54, 'VO', 'Orange', 'Merapi', '', 'Indonesia', 1, ''], [t0 - DAY, t0, 0, 0, 'EQ', 'Green', 'An earthquake', '', '', 2, '']] };
    const [found] = await ash.run(ins(winds(0)), params, context({ time: t0, events: { 'test-gdacs': gdacs } }));
    assert('plug-in volcanic ash: without a place, from the latest volcano in GDACS’s events', found.kind === 'sketch' && found.markers[0].lon === 110.44 && /Merapi/.test(found.title), found.title || found.text);
    const [none] = await ash.run(ins(winds(0)), params, context({ time: t0, events: {} }));
    assert('plug-in volcanic ash: and without either, a note', none.kind === 'note');
  }
}
