// Checks of the geography plug-ins (plugins/exposure.js, plugins/terrain.js),
// run by tests/run-tests.mjs: first on synthetic fields and regions whose answers
// are known, then on the library's products where they are built (GFS, GHS-POP,
// the regions, ETOPO). A plug-in fetches a regions raster from data/products/;
// here that address is served from the products' folder.

import { existsSync, readFileSync } from 'node:fs';
import { resolve } from 'node:path';

const HOUR = 3600000;

// fetch('data/products/<file>') answered from a map of buffers, or from PRODUCTS.
function serve(files, products) {
  const real = globalThis.fetch;
  globalThis.fetch = async (url) => {
    const m = /^data\/products\/(.+)$/.exec(String(url));
    if (!m) return real(url);
    const bytes = files[m[1]] || (products && existsSync(resolve(products, m[1])) ? readFileSync(resolve(products, m[1])) : null);
    if (!bytes) return { ok: false, status: 404, arrayBuffer: async () => new ArrayBuffer(0), json: async () => null };
    const copy = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
    return { ok: true, status: 200, arrayBuffer: async () => copy, json: async () => JSON.parse(Buffer.from(copy).toString('utf8')) };
  };
  return () => { globalThis.fetch = real; };
}

const gridOf = (d) => ({ lon0: -180 + d / 2, dlon: d, nlon: 360 / d, lat0: 90 - d / 2, dlat: d, nlat: 180 / d, periodic: true });
const fill = (g, f) => {
  const a = new Float32Array(g.nlon * g.nlat);
  for (let r = 0; r < g.nlat; r++) for (let c = 0; c < g.nlon; c++) a[r * g.nlon + c] = f(g.lon0 + c * g.dlon, g.lat0 - r * g.dlat);
  return a;
};
const input = (product, names, level = null) => ({ product, names, variables: names.map((n) => product.meta.variables.find((v) => v.name === n)), meta: product.meta, model: product.meta.title, level });

export async function checks({ check, assert, skip, PRODUCTS, api }) {
  const exposure = (await import('../plugins/exposure.js')).create(api);
  const terrain = (await import('../plugins/terrain.js')).create(api);

  // ---- exposure, synthetic: two countries, one hot and one windy ----------------------------
  {
    // Regions on a 1° grid: A (10°×10° at 0–10° E, 0–10° N) in two provinces, north and south; B at 20–30° E.
    const rg = gridOf(1), countries = new Uint16Array(rg.nlon * rg.nlat), provinces = new Uint16Array(rg.nlon * rg.nlat);
    for (let r = 0; r < rg.nlat; r++) for (let c = 0; c < rg.nlon; c++) {
      const lon = rg.lon0 + c, lat = rg.lat0 - r, i = r * rg.nlon + c;
      if (lat > 0 && lat < 10 && lon > 0 && lon < 10) { countries[i] = 1; provinces[i] = lat < 5 ? 1 : 2; }
      if (lat > 0 && lat < 10 && lon > 20 && lon < 30) { countries[i] = 2; provinces[i] = 3; }
    }
    const metas = {
      countries: { name: 'countries', title: 'Test countries', grid: rg, file: 'regions-test-countries.bin', regions: [
        { id: 1, code: 'AAA', name: 'Aland', country: 'AAA' }, { id: 2, code: 'BBB', name: 'Beland', country: 'BBB' }] },
      admin1: { name: 'admin1', title: 'Test provinces', grid: rg, file: 'regions-test-admin1.bin', regions: [
        { id: 1, code: 'AA-S', name: 'South', country: 'AAA' }, { id: 2, code: 'AA-N', name: 'North', country: 'AAA' }, { id: 3, code: 'BB-1', name: 'One', country: 'BBB' }] },
    };
    const restore = serve({ 'regions-test-countries.bin': Buffer.from(countries.buffer), 'regions-test-admin1.bin': Buffer.from(provinces.buffer) });
    // 1000 people in each 0.5° cell of the two countries; the forecast at 2°, three times.
    const pg = gridOf(0.5);
    const inland = (lon, lat) => lat > 0 && lat < 10 && ((lon > 0 && lon < 10) || (lon > 20 && lon < 30));
    const people = api.makeProduct({ name: 'test-people', title: 'People', domain: 'geography', grid: pg, times: [Date.UTC(2025, 6, 1)],
      variables: [{ name: 'population', long_name: 'People', standard_name: 'population', units: '1', slices: [fill(pg, (lon, lat) => (inland(lon, lat) ? 1000 : 0))] }] });
    const fg = gridOf(2), t0 = Date.UTC(2026, 9, 22), times = [t0, t0 + 6 * HOUR, t0 + 12 * HOUR];
    // 40 °C over A (and 6° round it) at the second time only, 20 °C elsewhere; wind 20 m/s over B (and round it) at the third.
    const hot = (lon, lat) => lon > -6 && lon < 16 && lat > -6 && lat < 16, windy = (lon, lat) => lon > 14 && lon < 36 && lat > -6 && lat < 16;
    const forecast = api.makeProduct({ name: 'test-forecast', title: 'Forecast', grid: fg, times, variables: [
      { name: 't2m', long_name: 'Temperature at 2 m', standard_name: 'air_temperature', units: 'K', level: '2 m', slices: times.map((t, k) => fill(fg, (lon, lat) => (k === 1 && hot(lon, lat) ? 313.15 : 293.15))) },
      { name: 'u10', long_name: 'Eastward wind at 10 m', standard_name: 'eastward_wind', units: 'm s-1', level: '10 m', slices: times.map((t, k) => fill(fg, (lon, lat) => (k === 2 && windy(lon, lat) ? 20 : 5))) },
      { name: 'v10', long_name: 'Northward wind at 10 m', standard_name: 'northward_wind', units: 'm s-1', level: '10 m', slices: times.map(() => fill(fg, () => 0)) },
    ] });
    const inputs = { temperature: input(forecast, ['t2m'], '2 m'), wind: input(forecast, ['u10', 'v10'], '10 m'), people: input(people, ['population']) };
    const context = (time) => ({ time, pinned: null, progress() {}, library: { regions: async (name) => metas[name] } });
    const heat = await exposure.run(inputs, { hazard: 'heat', heat: 35, wind: 17 }, context(t0));
    const [fieldResult, byCountry, byProvince, top] = heat;
    const exposed = fieldResult.product.slice('exposed', 0).reduce((a, x) => a + x, 0);
    check('plug-in exposure: heat over Aland exposes its 400 000 people', exposed - 400000, 1e-6, 'people');
    check('plug-in exposure: by country, Aland’s 400 000 and none of Beland', Math.abs(byCountry.values.AAA - 400000) + (byCountry.values.BBB || 0), 1e-6, 'people');
    check('plug-in exposure: by province, 200 000 in each half of Aland', Math.abs(byProvince.values['AA-S'] - 2e5) + Math.abs(byProvince.values['AA-N'] - 2e5), 1e-6, 'people');
    assert('plug-in exposure: results are a field, two regions (countries, provinces), a table and a note',
      heat.map((r) => r.kind).join(' ') === 'field regions regions table note' && byCountry.by === 'countries' && byProvince.by === 'admin1');
    assert('plug-in exposure: the table leads with Aland, all of its people', top.rows[0][0] === 'Aland' && top.rows[0][2] === '100.0 %', JSON.stringify(top.rows[0]));
    const later = (await exposure.run(inputs, { hazard: 'heat', heat: 35 }, context(t0 + 12 * HOUR)))[1];
    assert('plug-in exposure: a heat that has passed by the time shown exposes no one', Object.keys(later.values).length === 0);
    const wind = await exposure.run(inputs, { hazard: 'wind', heat: 35, wind: 17 }, context(t0));
    check('plug-in exposure: a gale over Beland exposes its 400 000', Math.abs(wind[1].values.BBB - 400000) + (wind[1].values.AAA || 0), 1e-6, 'people');
    check('plug-in exposure: the greatest wind sampled at the people’s cells is 20 m/s there', Math.abs(wind[0].product.slice('wspd10m_max', 0)[pg.nlon * Math.round((pg.lat0 - 5.25) / 0.5) + Math.round((25.25 - pg.lon0) / 0.5)] - 20), 1e-4, 'm/s');
    restore();
  }

  // ---- exposure on the library: GFS's heat over GHS-POP ------------------------------------
  {
    const need = ['fields-gfs-forecast.json', 'fields-ghsl-population.json', 'regions-countries.json', 'regions-admin1.json'];
    if (need.every((f) => existsSync(resolve(PRODUCTS, f)))) {
      const { FieldProduct } = await import('../src/data/fields.js');
      const load = (f) => { const m = JSON.parse(readFileSync(resolve(PRODUCTS, f), 'utf8')), b = readFileSync(resolve(PRODUCTS, m.file)); return new FieldProduct(m, b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); };
      const gfs = load('fields-gfs-forecast.json'), pop = load('fields-ghsl-population.json');
      const regions = async (name) => JSON.parse(readFileSync(resolve(PRODUCTS, `regions-${name}.json`), 'utf8'));
      const restore = serve({}, PRODUCTS);
      const inputs = { temperature: input(gfs, ['tmp2m'], '2 m'), wind: input(gfs, ['ugrd10m', 'vgrd10m'], '10 m'), people: input(pop, ['population']) };
      const t = Date.now();
      const out = await exposure.run(inputs, { hazard: 'heat', heat: 30 }, { time: gfs.times[0], pinned: null, progress() {}, library: { regions } });
      const ms = Date.now() - t;
      const total = out[0].product.slice('exposed', 0).reduce((a, x) => a + (x || 0), 0);
      const byCountry = Object.values(out[1].values).reduce((a, x) => a + x, 0);
      assert(`plug-in exposure on GFS and GHS-POP (${ms} ms): ${(total / 1e6).toFixed(1)} million above 30 °C, fewer than the world`, total > 0 && total < 8.2e9);
      check('plug-in exposure on GFS and GHS-POP: the countries hold all but 0.1 % of them', 1 - byCountry / total, 1e-3);
      restore();
    } else {
      skip('plug-in exposure on the library: GFS, GHS-POP or the regions missing');
    }
  }

  // ---- terrain, synthetic: a conical island -------------------------------------------------
  {
    const g = gridOf(0.5), apex = [10.25, 0.25];
    const z = fill(g, (lon, lat) => { const d = Math.hypot(lon - apex[0], lat - apex[1]); return d < 8 ? 2000 - 200 * d : -100; });
    const relief = api.makeProduct({ name: 'test-relief', title: 'Relief', domain: 'geography', grid: g, times: [], variables: [{ name: 'z', standard_name: 'altitude', units: 'm', slices: [z] }] });
    const context = (pinned) => ({ time: Date.now(), pinned, progress() {} });
    const [field, sketch, table] = await terrain.run({ relief: input(relief, ['z']) }, { cells: 0, snap: 0 }, context(apex));
    const up = field.product.slice('upstream_area', 0), i = Math.round((g.lat0 - apex[1]) / g.dlat) * g.nlon + Math.round((apex[0] - g.lon0) / g.dlon);
    const cellArea = (6371008.8 ** 2 * g.dlon * (Math.PI / 180) * (Math.sin(0.5 * Math.PI / 180) - Math.sin(0))) / 1e6;
    check('plug-in terrain: the apex of a cone drains only itself (float32)', up[i] / cellArea - 1, 1e-6);
    assert('plug-in terrain: its watershed is one cell, outlined by four corners', table.rows[1][1] === '1' && sketch.polylines[0].coords.length === 10);
    let sea = 0;
    for (let k = 0; k < up.length; k++) if (z[k] <= 0 && up[k] === up[k]) sea++;
    assert('plug-in terrain: no upstream area at sea', sea === 0);
    const coast = await terrain.run({ relief: input(relief, ['z']) }, { cells: 0, snap: 0 }, context([apex[0] + 6, apex[1]]));
    const path = coast[1].polylines[coast[1].polylines.length - 1].coords;
    assert('plug-in terrain: from the cone’s flank, water runs away from the apex to the sea', path[path.length - 2] > apex[0] + 7, `${path.slice(-2)}`);
  }

  // ---- terrain on ETOPO: the Murray at Murray Bridge, the Amazon at Óbidos -----------------------
  {
    const f = 'fields-etopo-2022-medium.json';
    if (existsSync(resolve(PRODUCTS, f))) {
      const { FieldProduct } = await import('../src/data/fields.js');
      const m = JSON.parse(readFileSync(resolve(PRODUCTS, f), 'utf8')), b = readFileSync(resolve(PRODUCTS, m.file));
      const etopo = new FieldProduct(m, b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength));
      const at = async (place) => {
        const t = Date.now(), out = await terrain.run({ relief: input(etopo, ['z']) }, { cells: 10, snap: 2 }, { time: Date.now(), pinned: place, progress() {} });
        return { ms: Date.now() - t, area: Number(out[2].rows[0][1].replace(/[^0-9]/g, '')) };
      };
      // The Murray–Darling basin, 1.05 million km² (World Resources Institute 2003), nearly all above Murray Bridge.
      const murray = await at([139.27, -35.12]);
      check(`plug-in terrain on ETOPO 10′ (${murray.ms} ms): the Murray at Murray Bridge drains ${(murray.area / 1e6).toFixed(2)} million km², against 1.05`, murray.area / 1.05e6 - 1, 0.1);
      // The Amazon at Óbidos drains about 4.7 million km² (the gauge 17050001 of Brazil's water
      // agency, ANA). On the 10′ grid the flat divide between the Madeira and the Tapajós is
      // misplaced, and Óbidos receives part of the Tapajós (about a sixth too much): within 20 %.
      const obidos = await at([-55.51, -1.92]);
      check(`plug-in terrain on ETOPO 10′ (${obidos.ms} ms): the Amazon at Óbidos drains ${(obidos.area / 1e6).toFixed(2)} million km², against 4.7`, obidos.area / 4.7e6 - 1, 0.2);
    } else {
      skip('plug-in terrain on the library: fields-etopo-2022-medium missing');
    }
  }
}
