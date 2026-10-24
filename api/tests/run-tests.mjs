// Checks of the API (api/README.md): the service is started on a free port with
// the repository's Python, and its endpoints are held to what they stand for —
// samples, series and grids to the Earth Atlas's own reading of the same product
// (FieldProduct, over physis's grid.js), calls into the domains to the same
// functions called here directly, readings posted to the stream that carries them
// within a second — together with the write policy with and without a token, the
// errors, the plug-ins and their jobs, and the clients: the Python client's
// checks (check_client.py, reporting JSON lines as data/tests/checks.py does), the
// C and C++ clients' where a compiler is found (check_cpp.mjs, built by
// api/Makefile), and the JavaScript client's, here; and the plug-ins over
// networks of places (check_networks.py), in Python and C++.
//
// The services started here keep the sensors' readings and the plug-ins' products
// in a directory of their own, deleted at the end, and build no products: the
// cache is read, never written. What the computer cannot do (no Python with numpy,
// no Node for the bridge, no products in the cache) is skipped, and said.
// Run from the repository root:  node api/tests/run-tests.mjs

import { spawn, spawnSync } from 'node:child_process';
import { existsSync, mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { createServer } from 'node:net';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

import { suite, PRODUCTS, ROOT } from '../../development/tools/checks.mjs';
import { Client } from '../js/science.mjs';
import { native } from './check_cpp.mjs';

const { check, assert, skip, done } = suite();
const HERE = import.meta.dirname;
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const pythons = [resolve(ROOT, '.venv/Scripts/python.exe'), resolve(ROOT, '.venv/bin/python'), 'python3', 'python'];
const python = pythons.find((p) => (p.includes('.venv') ? existsSync(p) : true) && spawnSync(p, ['-c', 'import numpy'], { stdio: 'ignore' }).status === 0);
const SOURCE = 'api-check';
const COMPUTE_SECONDS = 3;        // the main service's limit on a call (10 by default)
const HEARTBEAT = 0.5;            // and its heartbeats' interval (15 by default)
const work = mkdtempSync(join(tmpdir(), 'science-api-'));
const realStore = existsSync(resolve(ROOT, 'observation/sensors/store.py'));
const services = [];

async function freePort() {
  const s = createServer();
  await new Promise((r) => s.listen(0, '127.0.0.1', r));
  const { port } = s.address();
  await new Promise((r) => s.close(r));
  return port;
}

// A service on a free port: the sensors' store (or its stand-in) and the
// plug-ins' products in this run's own directory.
async function start(args = []) {
  const port = await freePort();
  const store = realStore ? [] : ['--sensor-store', resolve(HERE, 'stub_store.py')];
  const child = spawn(python, [resolve(ROOT, 'api/service/server.py'), '--port', String(port), '--product-seconds', '0', ...store, ...args], {
    cwd: ROOT, stdio: ['ignore', 'pipe', 'pipe'],
    env: { ...process.env, SCIENCE_OBSERVATIONS: join(work, 'observations'), SCIENCE_PLUGIN_OUTPUT: join(work, 'plugins'), PYTHONIOENCODING: 'utf-8' },
  });
  let log = '';
  child.stdout.on('data', (d) => { log += d; });
  child.stderr.on('data', (d) => { log += d; });
  const service = { child, port, base: `http://127.0.0.1:${port}/api/v1`, log: () => log };
  services.push(service);
  for (let i = 0; i < 150; i++) {
    if (child.exitCode !== null) break;
    try { if ((await fetch(`${service.base}/`)).ok) return service; } catch { /* not yet */ }
    await sleep(100);
  }
  throw new Error(`the service did not start: ${log.slice(-800)}`);
}

function stop(service) {
  // The tree: the server, and the Node workers of its compute bridge.
  if (process.platform === 'win32') spawnSync('taskkill', ['/pid', String(service.child.pid), '/T', '/F'], { stdio: 'ignore' });
  else service.child.kill();
}

const get = async (base, path, init) => {
  const r = await fetch(`${base}${path}`, init);
  const type = r.headers.get('content-type') || '';
  return { status: r.status, headers: r.headers, body: type.includes('json') ? await r.json() : await r.arrayBuffer() };
};
const post = (base, path, body, headers = {}) => get(base, path, { method: 'POST', headers: { 'Content-Type': 'application/json', ...headers }, body: JSON.stringify(body) });

// The page's own reading of a field product, to hold the service to.
async function pageField(name) {
  const meta = resolve(PRODUCTS, `fields-${name}.json`);
  const reader = resolve(ROOT, 'applications/earth-atlas/src/data/fields.js');
  if (!existsSync(meta) || !existsSync(reader)) return null;
  const { FieldProduct } = await import(pathToFileURL(reader).href);
  const m = JSON.parse(readFileSync(meta, 'utf8'));
  const bin = readFileSync(resolve(PRODUCTS, m.file));
  return new FieldProduct(m, bin.buffer.slice(bin.byteOffset, bin.byteOffset + bin.byteLength));
}

// The largest relative difference of two lists, null matching NaN (else infinite).
function relative(a, b) {
  let worst = 0;
  for (let i = 0; i < Math.max(a.length, b.length); i++) {
    const x = a[i] === null ? NaN : a[i], y = b[i];
    if (Number.isNaN(x) || Number.isNaN(y)) { if (Number.isNaN(x) !== Number.isNaN(y)) return Infinity; continue; }
    worst = Math.max(worst, Math.abs(x - y) / Math.max(Math.abs(y), 1e-12));
  }
  return worst;
}

async function fields(base) {
  const gfs = await pageField('gfs-forecast');
  if (!gfs) { skip('api fields: fields-gfs-forecast is not in the cache'); return; }
  const listed = await get(base, '/fields');
  assert('api fields: the list names gfs-forecast', listed.body.fields.some((f) => f.name === 'gfs-forecast'));
  const described = await get(base, '/fields/gfs-forecast');
  assert('api fields: a description without byte offsets, with quanta', described.body.variables.every((v) => !('offsets' in v) && v.quantum > 0) && !('file' in described.body));

  // Places over the globe, the poles and both sides of the antimeridian; instants between
  // stored times, at one, and beyond the span (clamped).
  const lons = [139.69, -0.13, 179.9, -179.95, 0, 359.9, 12.3, -73.97, 151.2, 37.6];
  const lats = [35.69, 51.51, -0.4, 10.2, 89.95, -89.99, -33.9, 40.78, -33.87, 55.75];
  const t = gfs.times;
  const instants = [t[3] + 2 * 3600e3, t[5], t[0] - 86400e3, t[t.length - 1] + 3600e3, t[1] + 1];
  let worst = 0;
  for (const v of gfs.meta.variables) {
    for (const ms of instants) {
      const r = await get(base, `/fields/gfs-forecast/sample?var=${v.name}&lon=${lons}&lat=${lats}&time=${new Date(ms).toISOString()}`);
      worst = Math.max(worst, relative(r.body.values, lons.map((lon, i) => gfs.sample(v.name, ms, lon, lats[i]))));
    }
  }
  check('api fields: samples equal the page\'s FieldProduct.sample (every variable, five instants, ten places)', worst, 1e-4);
  const many = await post(base, '/fields/gfs-forecast/sample', { var: 'tmp2m', lon: lons, lat: lats, time: new Date(instants[0]).toISOString() });
  check('api fields: samples asked by POST equal the page\'s', relative(many.body.values, lons.map((lon, i) => gfs.sample('tmp2m', instants[0], lon, lats[i]))), 1e-4);

  const ww3 = await pageField('ww3-forecast');
  const dir = ww3 && ww3.meta.variables.find((v) => /direction/.test(v.standard_name || '') || v.units === 'degree');
  if (dir) {
    const ms = ww3.times[2] + 3 * 3600e3;
    const pl = [[-150, 20], [140, 30], [-30, 45], [60, -40]];
    const r = await get(base, `/fields/ww3-forecast/sample?var=${dir.name}&lon=${pl.map((p) => p[0])}&lat=${pl.map((p) => p[1])}&time=${new Date(ms).toISOString()}`);
    let d = 0;
    pl.forEach(([lon, lat], i) => {
      const a = r.body.values[i], b = ww3.sampleDirection(dir.name, ms, lon, lat);
      if ((a === null) !== Number.isNaN(b)) d = Infinity;
      else if (a !== null) d = Math.max(d, Math.abs(((a - b + 540) % 360) - 180));
    });
    assert('api fields: a direction is interpolated through its unit vector', r.body.interpolation === 'direction');
    check('api fields: directions equal the page\'s sampleDirection', d, 1e-6, 'degree');
  } else {
    skip('api fields: no wave direction in the cache');
  }

  const series = await get(base, '/fields/gfs-forecast/series?var=tmp2m&lon=139.69&lat=35.69');
  check('api fields: a series equals the page\'s FieldProduct.series', relative(series.body.values, gfs.series('tmp2m', 139.69, 35.69)), 1e-4);
  assert('api fields: a series has every time', series.body.times.length === t.length && Date.parse(series.body.times[0]) === t[0]);

  // A grid across the antimeridian, every second point, at a stored time: the stored values themselves.
  const g = await get(base, '/fields/gfs-forecast/grid?var=tmp2m&time=3&bbox=170,-10,-170,10&stride=2');
  const s = gfs.slice('tmp2m', 3), G = gfs.grid;
  let gridWorst = 0;
  g.body.values.forEach((row, r) => row.forEach((value, c) => {
    const lon = g.body.grid.lon0 + c * g.body.grid.dlon, lat = g.body.grid.lat0 - r * g.body.grid.dlat;
    const col = ((Math.round((lon - G.lon0) / G.dlon) % G.nlon) + G.nlon) % G.nlon, rr = Math.round((G.lat0 - lat) / G.dlat);
    gridWorst = Math.max(gridWorst, Math.abs(value - s[rr * G.nlon + col]));
  }));
  assert('api fields: a grid across the antimeridian has its shape', g.body.grid.nlon === 11 && g.body.grid.nlat === 11 && g.body.values.length === 11, JSON.stringify(g.body.grid));
  check('api fields: its values are the stored ones', gridWorst, 0, 'K');
  const f32 = await get(base, '/fields/gfs-forecast/grid?var=tmp2m&time=3&bbox=170,-10,-170,10&stride=2&format=f32');
  const header = JSON.parse(f32.headers.get('x-grid'));
  const raw = new Float32Array(f32.body);
  check('api fields: float32 grids equal JSON grids', Math.max(...g.body.values.flat().map((v, i) => Math.abs(v - raw[i]))), 0, 'K');
  assert('api fields: the float32 grid is described in X-Grid', header.grid.nlon === 11 && raw.length === 121);
  // At an instant between two stored times, the grid is their linear blend, as the page draws it.
  const between = await get(base, `/fields/gfs-forecast/grid?var=tmp2m&time=${new Date(t[3] + 2 * 3600e3).toISOString()}&format=f32`);
  const blend = new Float32Array(between.body), a = gfs.slice('tmp2m', 3), c = gfs.slice('tmp2m', 4), w = 2 / 6;
  let blendWorst = 0;
  for (let i = 0; i < a.length; i += 97) blendWorst = Math.max(blendWorst, Math.abs(blend[i] - Math.fround(a[i] + (c[i] - a[i]) * w)));
  check('api fields: a grid between stored times is their blend', blendWorst, 1e-4, 'K');

  const exported = await get(base, '/fields/gfs-forecast/export?var=tmp2m&time=0&format=tif');
  if (exported.status === 501) skip('api fields: export answers 501 until data/formats/export.py exists');
  else {
    const bytes = (b) => new Uint8Array(b);
    assert('api fields: a GeoTIFF export begins as a little-endian TIFF does', exported.status === 200 && bytes(exported.body)[0] === 0x49 && bytes(exported.body)[2] === 42, `status ${exported.status}`);
    assert('api fields: the export names its file', /attachment; filename="gfs-forecast-tmp2m-[^"]+\.tif"/.test(exported.headers.get('content-disposition') || ''), exported.headers.get('content-disposition'));
    const nc = await get(base, '/fields/gfs-forecast/export?var=tmp2m&time=all&format=nc&bbox=120,20,150,50');
    assert('api fields: every time as netCDF classic', nc.status === 200 && new TextDecoder().decode(bytes(nc.body).slice(0, 3)) === 'CDF', `status ${nc.status}`);
    const one = await get(base, '/fields/gfs-forecast/export?var=tmp2m&time=all&format=tif');
    assert('api fields: every time is refused for a GeoTIFF (400)', one.status === 400, `${one.status}`);
  }

  // Errors.
  const errors = [
    ['/fields/no-such-field', 404], ['/fields/gfs-forecast/sample?var=nope&lon=0&lat=0', 404], ['/fields/gfs-forecast/sample?lon=0&lat=0', 400],
    ['/fields/gfs-forecast/sample?var=tmp2m&lon=east&lat=0', 400], ['/fields/gfs-forecast/sample?var=tmp2m&lon=0&lat=95', 400],
    ['/fields/gfs-forecast/sample?var=tmp2m&lon=0&lat=0&time=yesterday', 400], ['/fields/gfs-forecast/grid?var=tmp2m&format=png', 400],
    ['/fields/gfs-forecast/grid?var=tmp2m&time=999', 400], ['/fields/gfs-forecast/sample?var=tmp2m&lon=1,2,3&lat=1,2', 400],
  ];
  for (const [path, status] of errors) {
    const r = await get(base, path);
    assert(`api errors: ${path} answers ${status}`, r.status === status && typeof r.body.error === 'string', `${r.status} ${JSON.stringify(r.body)}`);
  }
}

async function compute(base) {
  const listed = await get(base, '/compute');
  if (listed.status === 503) { skip(`api compute: ${listed.body.error}`); return; }
  const units = listed.body.units;
  const astronomy = units.find((u) => u.unit === 'astronomy');
  assert('api compute: the units are the domains and foundations with a JavaScript interface', astronomy && units.some((u) => u.unit === 'physis' && u.kind === 'foundation'));
  const sun = astronomy.functions.find((f) => f.name === 'sunPosition');
  assert('api compute: a function is listed with its parameters and the comment above it', sun && sun.parameters[0] === 'jd' && /right ascension/.test(sun.doc), JSON.stringify(sun));
  const { sunPosition } = await import(pathToFileURL(resolve(ROOT, 'domains/astronomy/js/index.js')).href);
  const { windFromDirection } = await import(pathToFileURL(resolve(ROOT, 'domains/meteorology/js/index.js')).href);
  const { distanceKm } = await import(pathToFileURL(resolve(ROOT, 'domains/geodesy/js/index.js')).href);
  const jd = 2461000.25;
  const r = await post(base, '/compute/astronomy/sunPosition', { args: [jd] });
  const direct = sunPosition(jd);
  check('api compute: the Sun\'s position equals a direct call of astronomy\'s sunPosition',
    Math.max(...Object.keys(direct).map((k) => Math.abs(r.body.result[k] - direct[k]))), 0, 'degree');
  const w = await post(base, '/compute/meteorology/windFromDirection', [-3, 4]);
  check('api compute: wind direction equals meteorology\'s windFromDirection', Math.abs(w.body.result - windFromDirection(-3, 4)), 0, 'degree');
  const d = await post(base, '/compute/geodesy/distanceKm', { args: [35.69, 139.69, 51.51, -0.13] });
  check('api compute: a distance equals geodesy\'s distanceKm', Math.abs(d.body.result - distanceKm(35.69, 139.69, 51.51, -0.13)), 0, 'km');
  const grid = { lon0: 0, dlon: 1, nlon: 2, lat0: 1, dlat: 1, nlat: 2, periodic: false };
  const b = await post(base, '/compute/physis/bilinear', { args: [{ $typed: 'Float32Array', values: [1, 2, 3, null] }, grid, 0.25, 0.75] });
  // At (0.25, 0.75) the weights are 9/16, 3/16, 3/16 and 1/16, the last on the missing corner.
  check('api compute: a typed array argument, with a missing value renormalised', Math.abs(b.body.result - (9 * 1 + 3 * 2 + 3 * 3) / 15), 1e-12);
  const someClass = units.flatMap((u) => u.classes.map((c) => [u.unit, c]))[0];
  const refusals = [['/compute/astronomy/noSuchFunction', [], 404], ['/compute/no-such-unit/f', [], 404], ['/compute/astronomy/J2000', [], 404],
    ['/compute/astronomy/sunPosition', { args: 'jd' }, 400]];
  if (someClass) refusals.push([`/compute/${someClass[0]}/${someClass[1]}`, [], 404]);
  for (const [path, body, status] of refusals) {
    const x = await post(base, path, body);
    assert(`api compute: ${path} answers ${status}`, x.status === status, `${x.status} ${JSON.stringify(x.body)}`);
  }
  const failed = await post(base, '/compute/physis/bilinear', { args: [null, null, 0, 0] });
  assert('api compute: a function that throws answers 422 with its error', failed.status === 422 && /TypeError/.test(failed.body.error), JSON.stringify(failed.body));
  // Gauss-Legendre nodes of degree 200 000 take minutes: the call is stopped at the
  // service's limit (COMPUTE_SECONDS here), and its worker replaced.
  const t0 = performance.now();
  const slow = await post(base, '/compute/physis/gaussLegendre', [200000]);
  const seconds = (performance.now() - t0) / 1000;
  assert('api compute: a call past its time limit answers 504', slow.status === 504, `${slow.status} ${JSON.stringify(slow.body)}`);
  check('api compute: ... and is stopped at the limit', seconds - COMPUTE_SECONDS, 2, 's');
  const after = await post(base, '/compute/meteorology/windSpeed', [3, 4]);
  assert('api compute: the next call is answered by a fresh worker', after.status === 200 && after.body.result === 5, `${after.status}`);
}

// A stream with no readings sends a heartbeat every HEARTBEAT seconds.
async function heartbeats(base) {
  const control = new AbortController();
  const response = await fetch(`${base}/stream`, { signal: control.signal });
  const reader = response.body.getReader();
  const decoder = new TextDecoder();
  let text = '';
  const until = performance.now() + 1000 * (2.6 * HEARTBEAT);
  while (performance.now() < until) {
    const chunk = await Promise.race([reader.read(), sleep(until - performance.now()).then(() => null)]);
    if (!chunk || chunk.done) break;
    text += decoder.decode(chunk.value, { stream: true });
  }
  control.abort();
  const beats = (text.match(/^: heartbeat /gm) || []).length;
  assert('api observations: a quiet stream sends its heartbeats', response.headers.get('content-type').startsWith('text/event-stream') && beats >= 2, `${beats} in ${2.6 * HEARTBEAT} s`);
}

async function observations(base) {
  const define = { title: 'API check', variables: ['air_temperature', 'relative_humidity'] };
  const at = () => new Date(Date.now() - 1000).toISOString();
  // A stream opened before the readings are posted carries them within a second.
  const api = new Client(base);
  const control = new AbortController();
  const received = [];
  const listening = (async () => {
    for await (const event of api.stream(SOURCE, { signal: control.signal })) { received.push([performance.now(), event]); control.abort(); }
  })().catch((error) => received.push([performance.now(), { error: String(error) }]));
  await sleep(500);
  const t0 = performance.now();
  const first = await post(base, `/observations/${SOURCE}`, {
    define, readings: [{ station: 'js1', lon: 139.76, lat: 35.68, time: at(), values: { air_temperature: 290.5, relative_humidity: 61 } },
      { station: 'js1', lon: 139.76, lat: 35.68, time: at(), values: { air_temperature: 1000 } }],
  });
  await Promise.race([listening, sleep(5000)]);
  control.abort();
  assert('api observations: a source is defined and its readings appended', first.status === 200 && first.body.accepted === 1, JSON.stringify(first.body));
  assert('api observations: an implausible reading is refused, with its reason', first.body.rejected?.length === 1 && first.body.rejected[0].index === 1, JSON.stringify(first.body.rejected));
  if (received.length && received[0][1].readings) {
    check('api observations: a posted reading reaches the stream', (received[0][0] - t0) / 1000, 1, 's');
    assert('api observations: the stream carries the source and the reading', received[0][1].source === SOURCE && received[0][1].readings[0].station === 'js1');
    // A client that reconnects with the last id it saw receives what it missed.
    const missedControl = new AbortController();
    const second = post(base, `/observations/${SOURCE}`, [{ station: 'js2', lon: 139.7, lat: 35.6, time: at(), values: { air_temperature: 289 } }]);
    await second;
    let missed = null;
    for await (const event of api.stream(null, { lastEventId: received[0][1].id, signal: missedControl.signal })) { missed = event; break; }
    assert('api observations: a reconnecting stream receives what it missed (Last-Event-ID)', missed && missed.readings[0].station === 'js2', JSON.stringify(missed));
  } else {
    assert('api observations: a posted reading reaches the stream', false, JSON.stringify(received));
  }
  const back = await get(base, `/observations/${SOURCE}?limit=1`);
  assert('api observations: the latest reading reads back', back.body.readings?.length === 1 && back.body.readings[0].station === 'js2', JSON.stringify(back.body).slice(0, 300));
  const listed = await get(base, '/observations');
  assert('api observations: the list names the sensor source', listed.body.sensors.some((s) => s.name === SOURCE));
  const product = listed.body.products[0];
  if (product) {
    const rows = await get(base, `/observations/${product.name}?limit=3`);
    assert(`api observations: the set ${product.name}, its latest three rows`, rows.status === 200 && rows.body.rows.length === Math.min(3, rows.body.of), `${rows.status}`);
  } else {
    skip('api observations: no observation set in the cache');
  }
  const cross = await post(base, `/observations/${SOURCE}`, { readings: [] }, { Origin: 'http://elsewhere.example' });
  assert('api policy: without a token, a page of another origin may not write', cross.status === 403, `${cross.status}`);
  const unknown = await post(base, '/observations/no-such-source', { readings: [] });
  assert('api observations: posting to a source never defined answers 404', unknown.status === 404, `${unknown.status} ${JSON.stringify(unknown.body)}`);
  const bad = await get(base, `/observations/${SOURCE}`, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: '{not json' });
  assert('api observations: a body that is not JSON answers 400', bad.status === 400, `${bad.status}`);
}

async function products(base) {
  for (const kind of ['events', 'tracks', 'regions', 'verification']) {
    const listed = await get(base, `/${kind}`);
    const first = listed.body[kind]?.[0];
    if (!first) { skip(`api products: no ${kind} in the cache`); continue; }
    const one = await get(base, `/${kind}/${first.name}`);
    const file = JSON.parse(readFileSync(resolve(PRODUCTS, first.file), 'utf8'));
    assert(`api products: ${kind}/${first.name} is served as the data library wrote it`, one.status === 200 && JSON.stringify(one.body) === JSON.stringify(file));
  }
  const quakes = (await get(base, '/events')).body.events.find((e) => e.name.startsWith('usgs-earthquakes'));
  if (quakes) {
    const r = await get(base, `/events/${quakes.name}?bbox=120,20,150,50&limit=5`);
    const ok = r.body.events.every((row) => row[1] >= 120 && row[1] <= 150 && row[2] >= 20 && row[2] <= 50) && r.body.events.length <= 5;
    assert('api products: events filtered by a box and a limit', r.status === 200 && ok, `${r.body.count} of ${r.body.of}`);
    const geo = await get(base, `/events/${quakes.name}?format=geojson&bbox=120,20,150,50&limit=5`);
    if (geo.status === 501) skip('api products: GeoJSON answers 501 until data/formats/export.py exists');
    else {
      const fc = geo.body;
      assert('api products: the same events as GeoJSON points', geo.status === 200 && /geo\+json/.test(geo.headers.get('content-type'))
        && fc.type === 'FeatureCollection' && fc.features.length === r.body.events.length
        && fc.features.every((f) => f.geometry.type === 'Point'), `${geo.status} ${fc.features?.length}`);
    }
  }
  const regions = (await get(base, '/regions')).body.regions.find((x) => x.name === 'countries');
  if (regions) {
    const r = await get(base, '/regions/countries/at?lon=139.69,2.35,-150&lat=35.69,48.86,-30');
    assert('api products: the regions at places (Tokyo in Japan, Paris in France, none at sea)',
      r.body.at[0]?.code === 'JPN' && r.body.at[1]?.code === 'FRA' && r.body.at[2] === null, JSON.stringify(r.body.at.map((x) => x && x.code)));
  }
  for (const [path, status] of [['/events/no-such-set', 404], ['/no-such-endpoint', 404], ['/plugins/run', 405]]) {
    const r = await get(base, path);
    assert(`api errors: ${path} answers ${status}`, r.status === status, `${r.status}`);
  }
  const index = await get(base, '/');
  const spec = JSON.parse(readFileSync(resolve(ROOT, 'api/openapi.json'), 'utf8'));
  const missing = index.body.endpoints.filter((e) => !spec.paths[e.path.replace('/api/v1', '') || '/']?.[e.method.toLowerCase()]);
  assert('api: openapi.json describes every endpoint the service lists', missing.length === 0, missing.map((e) => `${e.method} ${e.path}`).join(', '));
}

async function plugins(base) {
  const listed = await get(base, '/plugins');
  const atlas = listed.body.applications?.['earth-atlas'] || [];
  assert('api plug-ins: the Earth Atlas\'s Python plug-ins are listed', atlas.some((p) => p.file === 'plugins/vorticity.py' && p.runtime === 'python'), JSON.stringify(listed.body).slice(0, 200));
  const d = await get(base, '/plugins/describe?app=earth-atlas&file=plugins/vorticity.py');
  assert('api plug-ins: a description with its needs', d.status === 200 && d.body.id === 'vorticity-py' && d.body.needs[0].name === 'wind');
  const outside = await get(base, '/plugins/describe?app=earth-atlas&file=serve.py');
  assert('api plug-ins: a file outside plugins/ is refused', outside.status === 400, `${outside.status}`);
  if (!d.body.needs?.[0]?.chosen) { skip('api plug-ins: no wind in the cache to run the plug-in on'); return; }
  const api = new Client(base);
  const run = api.runPlugin('earth-atlas', 'plugins/vorticity.py', { poll: 300 });
  await sleep(200);
  const second = await post(base, '/plugins/run', { app: 'earth-atlas', file: 'plugins/vorticity.py' });
  assert('api plug-ins: one run at a time (a second answers 409)', second.status === 409, `${second.status}`);
  let job;
  try { job = await run; } catch (error) { job = error.job || { state: 'failed', message: String(error) }; }
  assert('api plug-ins: a run finishes, with a field among its results', job.state === 'done' && job.results.some((r) => r.kind === 'field' && r.api), job.message || '');
  const field = job.results?.find((r) => r.kind === 'field');
  if (field) {
    const made = await get(base, field.api.replace('/api/v1', ''));
    assert('api plug-ins: the run\'s field is served among the fields', made.status === 200 && made.body.variables.some((v) => v.name === 'zeta'));
  }
}

async function policy() {
  const token = 'check-token-2f6c';
  const guarded = await start(['--token', token]);
  try {
    const reading = { define: { title: 'API check', variables: ['air_temperature'] }, readings: [] };
    const without = await post(guarded.base, `/observations/${SOURCE}-t`, reading);
    const wrong = await post(guarded.base, `/observations/${SOURCE}-t`, reading, { Authorization: 'Bearer wrong' });
    const right = await post(guarded.base, `/observations/${SOURCE}-t`, reading, { Authorization: `Bearer ${token}` });
    const read = await get(guarded.base, '/fields');
    assert('api policy: with a token, a write without it answers 401', without.status === 401, `${without.status}`);
    assert('api policy: with a token, a wrong token answers 401', wrong.status === 401, `${wrong.status}`);
    assert('api policy: with a token, a write with it is accepted', right.status === 200, `${right.status} ${JSON.stringify(right.body)}`);
    assert('api policy: with a token, reads need none', read.status === 200);
    const index = await get(guarded.base, '/');
    assert('api policy: the index states the policy', index.body.writes.policy === 'token' && index.body.writes.writable === false);
  } finally {
    stop(guarded);
  }
  const closed = await start(['--read-only']);
  try {
    const r = await post(closed.base, `/observations/${SOURCE}`, { readings: [] });
    const c = await post(closed.base, '/compute/meteorology/windSpeed', [3, 4]);
    assert('api policy: a read-only service refuses writes (403)', r.status === 403, `${r.status}`);
    assert('api policy: a read-only service still computes', c.status === 200 && c.body.result === 5, `${c.status}`);
  } finally {
    stop(closed);
  }
}

if (!python) {
  skip('api: no Python with numpy found, so no service to check');
} else {
  let service;
  try {
    service = await start(['--compute-seconds', String(COMPUTE_SECONDS), '--heartbeat-seconds', String(HEARTBEAT)]);
    const base = service.base;
    if (!existsSync(resolve(PRODUCTS, 'index.json'))) skip('api: the data cache has no products; the checks that read them are skipped');
    await fields(base);
    await compute(base);
    await observations(base);
    await heartbeats(base);
    await products(base);
    await plugins(base);

    // The Python client, against the same service.
    const r = spawnSync(python, [resolve(HERE, 'check_client.py'), `http://127.0.0.1:${service.port}`, `${SOURCE}-py`],
      { cwd: ROOT, encoding: 'utf8', timeout: 300000, maxBuffer: 1 << 26, env: { ...process.env, PYTHONIOENCODING: 'utf-8' } });
    let seen = 0;
    for (const line of (r.stdout || '').split('\n')) {
      if (!line.startsWith('{')) continue;
      let rec;
      try { rec = JSON.parse(line); } catch { continue; }
      seen++;
      if (rec.kind === 'check') check(rec.name, rec.value, rec.limit, rec.unit || '');
      else if (rec.kind === 'assert') assert(rec.name, rec.ok, rec.detail || '');
      else if (rec.kind === 'skip') skip(rec.name);
    }
    assert('api: check_client.py ran to the end', r.status === 0 && seen > 0, (r.stderr || '').trim().split('\n').slice(-4).join(' | '));

    // The C and C++ clients, where a compiler is found (check_cpp.mjs says how, and where they run).
    await native({ service, python, freePort, check, assert, skip, work });

    // The JavaScript client's own surface.
    const api = new Client(`http://127.0.0.1:${service.port}`);
    const info = await api.info();
    assert('js client: the service\'s index', info.api === 'v1' && info.writes.policy === 'open');
    if (existsSync(resolve(PRODUCTS, 'fields-gfs-forecast.json'))) {
      const gfs = await api.field('gfs-forecast');
      const one = await gfs.sample('tmp2m', 139.69, 35.69, { time: gfs.times[1] });
      const page = await pageField('gfs-forecast');
      check('js client: a sample equals the page\'s', Math.abs(one - page.sample('tmp2m', gfs.times[1].getTime(), 139.69, 35.69)) / one, 1e-4);
      const { values, grid } = await gfs.grid('tmp2m', { time: 0, bbox: [120, 20, 150, 50] });
      assert('js client: a grid is a Float32Array of its shape', values instanceof Float32Array && values.length === grid.nlon * grid.nlat);
      try {
        await api.field('no-such-field');
        assert('js client: a missing field throws', false);
      } catch (error) {
        assert('js client: a missing field throws an ApiError of status 404', error.status === 404, error.message);
      }
    }
    const units = await api.units().catch(() => null);
    if (units) check('js client: compute', Math.abs((await api.compute('meteorology', 'windSpeed', 3, 4)) - 5), 0, 'm/s');
  } catch (error) {
    assert('api: the checks ran', false, `${error.stack || error}\n${service ? service.log().slice(-1500) : ''}`);
  } finally {
    if (service) stop(service);
  }
  try {
    await policy();
  } catch (error) {
    assert('api policy: the services with a token and read-only ran', false, String(error.stack || error));
  }
  for (const s of services) stop(s);

  // The plug-ins over networks of places (check_networks.py), which need no service.
  const r = spawnSync(python, [resolve(HERE, 'check_networks.py')],
    { cwd: ROOT, encoding: 'utf8', timeout: 600000, maxBuffer: 1 << 24, env: { ...process.env, PYTHONIOENCODING: 'utf-8' } });
  let seen = 0;
  for (const line of (r.stdout || '').split('\n')) {
    if (!line.startsWith('{')) continue;
    let rec;
    try { rec = JSON.parse(line); } catch { continue; }
    seen++;
    if (rec.kind === 'check') check(rec.name, rec.value, rec.limit, rec.unit || '');
    else if (rec.kind === 'assert') assert(rec.name, rec.ok, rec.detail || '');
    else if (rec.kind === 'skip') skip(rec.name);
  }
  assert('api: check_networks.py ran to the end', r.status === 0 && seen > 0, (r.stderr || '').trim().split('\n').slice(-4).join(' | '));
}
rmSync(work, { recursive: true, force: true });
done();
