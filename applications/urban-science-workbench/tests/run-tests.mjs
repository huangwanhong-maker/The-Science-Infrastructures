// Checks of the Urban Science Workbench: its plug-ins' manifests (every layer's and tool's module present;
// planned plug-ins with what they will hold), its areas and recipe, the syntax of every page module; and, when the engine of
// Street Oceanography is built from the present sources, its answers: the paper's Figure 1 read right by the
// amendment and wrongly by the paper in the fourth case; the catalytic pathway's coupling length against the
// killed walk's closed form; the kinematic street conserving its crowd; the district's detailed balance without
// circulation and its affinities untouched by attraction; and, when the data library's products are present,
// a sensor's rhythm (the day forced, harmonic constants resolved), the civil share of its clock, a street's
// hours and mixing (M = between + within); the tides (high water on a weekday, a day's prediction near its
// counts), the current meters (each way summing to the total, the handedness test), a street's crowd second
// by second (conserved, a red crossing passing nobody) and a visitor's plume (its masses summing to one).
// When the engine of Human flow is built and the records are, on Melbourne's district as the page sends it: the walkway
// and its counters, the departures calibrated each way, the walkers' counters reading what the calibrated flow expects,
// a quarter's frames inside the district, a street closed to walkers passing through, an event's crowd and its record
// (every event simulated), the permits, and one decision of POMCP on the walkers.
// When the engine of Landscape Dynamics is built, on a grid district with a hole: its cells, the lane's rewrites
// (refinements, then the split) undone and planned again, a day's run (walkers kept, each way's hops), the event
// skeleton, the scales (harmonic flows as many as holes, memory keeping every walker) and a building's events.
// The mathematics itself is checked in the foundations and in the urban-science domain.
// Run from the repository root:  node applications/urban-science-workbench/tests/run-tests.mjs

import { existsSync, mkdtempSync, readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { spawnSync } from 'node:child_process';

import { suite, ROOT } from '../../../development/tools/checks.mjs';

const { assert, skip, done } = suite();
const APP = resolve(ROOT, 'applications/urban-science-workbench');

// ---- plug-ins ---------------------------------------------------------------------------------------
const plugins = readdirSync(join(APP, 'plugins')).filter((p) => existsSync(join(APP, 'plugins', p, 'plugin.json')));
assert('at least one plug-in with layers or tools, and three planned', plugins.length >= 4, plugins.join(', '));
for (const p of plugins) {
  const m = JSON.parse(readFileSync(join(APP, 'plugins', p, 'plugin.json'), 'utf8'));
  if (m.status === 'planned') {
    assert(`plug-in ${p}: planned, with what it will hold and no views`, !!m.title && (m.planned || []).length > 0 && (m.views || []).length === 0);
    continue;
  }
  const parts = [...(m.layers || []), ...(m.tools || [])];
  const missing = parts.filter((v) => !existsSync(join(APP, 'plugins', p, v.module)));
  assert(`plug-in ${p}: a title, layers or tools, and every module present`, !!m.title && parts.length > 0 && missing.length === 0, missing.map((v) => v.id).join(', '));
  if (m.engine) assert(`plug-in ${p}: its engine's source present`, existsSync(join(APP, 'plugins', p, m.engine)));
}
{
  const areas = JSON.parse(readFileSync(join(APP, 'areas.json'), 'utf8')).areas;
  assert('the areas: each with an origin and a time zone', areas.length >= 1 && areas.every((a) => Number.isFinite(a.lat) && Number.isFinite(a.lon) && a.timezone), areas.map((a) => a.id).join(', '));
}
{
  const recipe = JSON.parse(readFileSync(join(APP, 'recipe.json'), 'utf8'));
  assert('the recipe: every step a script that exists', recipe.steps.every((s) => existsSync(join(ROOT, s.run[0]))), recipe.steps.map((s) => s.run[0]).join(', '));
}

// ---- Human flow's colours: the walkers' patterns, every two apart in every vision -------------------------
{
  const src = readFileSync(join(APP, 'plugins/human-flow/shared.js'), 'utf8');
  const list = (name) => JSON.parse(src.match(new RegExp(`export const ${name} = (\\[[^\\]]*\\])`))[1].replace(/'/g, '"'));
  const S = await import(new URL(`file:///${join(ROOT, 'visualization/styles/index.js').replace(/\\/g, '/')}`));
  const least = (hexes) => {
    const rgb = hexes.map(S.hexRgb), out = {};
    for (const kind of ['normal', ...S.DEFICIENCIES]) {
      const seen = rgb.map((c) => (kind === 'normal' ? c : S.simulate(c, kind).map((v) => Math.min(1, Math.max(0, v)))));
      let m = Infinity;
      for (let i = 0; i < seen.length; i++) for (let j = i + 1; j < seen.length; j++) m = Math.min(m, S.deltaE(seen[i], seen[j]));
      out[kind] = m;
    }
    return out;
  };
  const modes = least(list('MODE_COLOURS')), status = least(list('STATUS_COLOURS').slice(1, 5));
  assert(`human flow: the walkers' six patterns apart in every vision (least OKLab ΔE ${Math.min(...Object.values(modes)).toFixed(3)}), and what they do (${Math.min(...Object.values(status)).toFixed(3)})`,
    Object.values(modes).every((d) => d >= 0.09) && Object.values(status).every((d) => d >= 0.06), JSON.stringify(modes));
}

// ---- the page's modules -----------------------------------------------------------------------------
{
  // node --check on a .js ES module may pass syntax errors; each module is checked as .mjs.
  const dir = mkdtempSync(join(tmpdir(), 'urban-check-'));
  const files = [];
  const walk = (d) => { for (const f of readdirSync(d, { withFileTypes: true })) { const p = join(d, f.name); if (f.isDirectory()) walk(p); else if (f.name.endsWith('.js')) files.push(p); } };
  walk(join(APP, 'src'));
  for (const p of plugins) for (const d of ['tools', 'layers']) if (existsSync(join(APP, 'plugins', p, d))) walk(join(APP, 'plugins', p, d));
  const bad = [];
  files.forEach((f, i) => {
    const copy = join(dir, `m${i}.mjs`);
    writeFileSync(copy, readFileSync(f, 'utf8'));
    const r = spawnSync(process.execPath, ['--check', copy], { encoding: 'utf8' });
    if (r.status !== 0) bad.push(`${f.slice(APP.length + 1)}: ${r.stderr.split('\n').find((l) => l.includes('Error')) || r.stderr}`);
  });
  assert(`the page's ${files.length} modules parse`, bad.length === 0, bad.join('; '));
}

// ---- the engine, when it is built -------------------------------------------------------------------
{
  const status = spawnSync('python', [join(APP, 'engine.py'), 'status'], { cwd: ROOT, encoding: 'utf8' });
  let state = null;
  try { state = JSON.parse(status.stdout).find((s) => s.plugin === 'street-oceanography'); } catch { /* not built */ }
  if (!state || !state.current) {
    skip('the engine of Street Oceanography is not built from the present sources (python applications/urban-science-workbench/engine.py build)');
  } else {
    const ask = (request) => {
      const r = spawnSync('python', [join(APP, 'engine.py'), 'ask', 'street-oceanography', JSON.stringify(request)],
        { cwd: ROOT, encoding: 'utf8', maxBuffer: 1 << 28, env: { ...process.env, PYTHONIOENCODING: 'utf-8' } });
      const answer = JSON.parse(r.stdout);
      if (!answer.ok) throw new Error(answer.error);
      return answer.result;
    };
    const cases = ask({ op: 'constructed' });
    assert('Figure 1: four cases; the amendment reads the fourth (linear, non-separable) as proportional, the paper as disproportionate',
      cases.length === 4 && cases[3].amended === 'proportional response' && cases[3].paper === 'disproportionate response' && cases[1].amended === 'disproportionate response',
      cases.map((c) => `${c.paper} / ${c.amended}`).join('; '));
    const pw = ask({ op: 'pathway', search: 0, length: 600, bridal: 300 });
    assert('pathway: the wandering coupling length is the killed walk\'s closed form (2%)', Math.abs(pw.walk.length / pw.scales.walkClosedForm - 1) < 0.02,
      `${pw.walk.length} against ${pw.scales.walkClosedForm}`);
    const sd = ask({ op: 'streetday', peak: 8, blocks: 3, from: 11, to: 13 });
    assert('street day: the crowd conserved on both days (to a millionth of those who entered)',
      Math.abs(sd.ordinary.budgetResidual) <= 1e-6 * Math.max(1, sd.ordinary.entered.at(-1)) && Math.abs(sd.rule.budgetResidual) <= 1e-6 * Math.max(1, sd.rule.entered.at(-1)));
    const flat = ask({ op: 'district', circulation: 0, holeCirculation: 0, holes: [] });
    assert('district: an attraction alone satisfies detailed balance (σ ≈ 0)', flat.detailedBalance && flat.entropyProduction < 1e-10, String(flat.entropyProduction));
    const round = ask({ op: 'district', circulation: 0.8 });
    assert('district: circulation produces entropy, and the attraction changes no cycle affinity', round.entropyProduction > 1e-6 && round.affinityGapWithoutAttraction < 1e-9);

    const products = ask({ op: 'status' });
    if (!products['street-sensors.json']) skip('the data library\'s street products are not built (the recipe)');
    else {
      const sensors = ask({ op: 'sensors' });
      const mall = sensors.find((s) => /Bourke Street Mall \(North\)/.test(s.name)) || sensors[0];
      const rh = ask({ op: 'rhythm', sensor: mall.id, era: 'recent' });
      const day = rh.peaks.find((p) => p.name === 'D1');
      assert(`rhythm (${mall.name}): the daily line is forced and far above the integration null`, day && day.verdict === 'forced' && day.ratio > 10, day ? `× ${day.ratio}` : 'no D1');
      assert('rhythm: harmonic constants resolved and the selectivity in [0, 1]', rh.harmonics.constituents > 20 && rh.selectivity.S >= 0 && rh.selectivity.S <= 1,
        `${rh.harmonics.constituents} constituents, S = ${rh.selectivity.S}`);
      const profiles = ask({ op: 'profiles', era: 'recent' });
      const pm = profiles.find((x) => x.id === mall.id);
      assert("profiles: every outdoor sensor's week, 168 hours, busier at 13:00 on a weekday than at 04:00", pm && pm.profile.length === 168 && pm.profile[13] > pm.profile[4],
        pm ? `${pm.profile[13]} against ${pm.profile[4]}` : 'none');
      const clock = ask({ op: 'clock', sensor: mall.id });
      assert('clock: the changes of daylight saving inside the record are found', clock.shifts.length >= 2, String(clock.shifts.length));
      // The street read as a sea.
      const tides = ask({ op: 'tides' });
      const tm = tides.find((t) => t.id === mall.id);
      assert(`tides: ${tides.length} streets; ${mall.name}'s weekday high water between 11:00 and 15:00, its mean week 168 hours`,
        tides.length >= 80 && tm && tm.weekday.highAt > 11 && tm.weekday.highAt < 15 && tm.week.length === 168, tm ? `${tm.weekday.highAt} h` : 'none');
      const sea = ask({ op: 'sea', date: '2025-03-12' });
      let pred = 0, obs = 0;
      for (const x of sea.sensors) for (let hr = 0; hr < 24; hr++) if (Number.isFinite(x.observed[hr])) pred += x.level[4 * hr + 2], obs += x.observed[hr];
      assert('sea: on an ordinary Wednesday the city\'s counts lie within 15% of its tide', sea.observed > 50 && Math.abs(obs / pred - 1) < 0.15, `${sea.observed} streets, observed/tide ${(obs / pred).toFixed(3)}`);
      const both = sea.sensors.filter((x) => x.way1.some(Number.isFinite));
      const sums = both.every((x) => x.observed.every((v, hr) => !Number.isFinite(v) || !Number.isFinite(x.way1[hr]) || x.way1[hr] + x.way2[hr] === v));
      assert('sea: each way sums to the count, hour by hour', both.length > 50 && sums, String(both.length));
      const cur = ask({ op: 'currents' });
      assert(`currents: ${cur.sensors.length} current meters, the residual the mean of the net current; handedness ${cur.handedness.left} of ${cur.handedness.footpaths}`,
        cur.sensors.length >= 60 && cur.sensors.every((m) => Math.abs(m.net.reduce((a, v) => a + v, 0) / 168 - m.residual) < 0.06) && cur.handedness.footpaths >= 30 && cur.handedness.pLeft > 0 && cur.handedness.pLeft <= 1);
      const wave = ask({ op: 'streetwave', hour: 12, blockLengths: [140, 90, 160], signalled: [1, 0, 1], demand: Array.from({ length: 24 }, (_, hr) => (hr >= 7 && hr < 20 ? 0.4 : 0)), window: 120, step: 2 });
      const at = wave.interfaces.findIndex((x) => Math.abs(x - wave.crossings[0]) < 1e-6);
      assert('streetwave: 61 snapshots of a street of three blocks, conserved, a red crossing passing nobody',
        wave.times.length === 61 && Math.abs(wave.budgetResidual) < 1e-6 * Math.max(1, wave.entered) && wave.green.every((g, k) => g[0] || Math.abs(wave.flow[k][at]) < 1e-5));
      const pl = ask({ op: 'plume', length: 300, bridal: 100, ring: 140, horizon: 3600, steps: 30 });
      const mass = pl.times.map((_, k) => pl.bought[k] + pl.left[k] + pl.errand[k].reduce((a, v) => a + v, 0) + pl.need[k].reduce((a, v) => a + v, 0));
      assert('plume: the visitor\'s masses sum to one at every time, and P(bought) stays below the yield', mass.every((m) => Math.abs(m - 1) < 1e-5) && pl.bought.at(-1) <= pl.yield + 1e-9);
      if (products['street-holdings.json']) {
        const streets = ask({ op: 'streets' });
        const g = streets.find((s) => s.id === 'melbourne-swanston') || streets[0];
        const st = ask({ op: 'street', street: g.id });
        assert(`street (${g.name}): holdings, and M = between + within`, st.holdings.length > 0 && (!st.mixing || Math.abs(st.mixing.between + st.mixing.within - st.mixing.departure) < 1e-9));
      } else skip('OpenStreetMap\'s holdings are not built');
      // New York: the Brooklyn Bridge's counter, an area of its own, with its own weather and holidays.
      if (existsSync(join(ROOT, 'data/cache/products/street-areas.json')) && existsSync(join(ROOT, 'data/cache/products/street-sensors-new-york.json'))) {
        const ny = ask({ op: 'tides', area: 'new-york' });
        const bridge = ny[0];
        assert(`New York: one tide, the Brooklyn Bridge's, one crest a day (form ${bridge && bridge.form.toFixed(2)}), weekday high water at ${bridge && bridge.weekday.highAt} h`,
          ny.length === 1 && bridge.formClass === 'diurnal' && bridge.weekday.highAt > 13 && bridge.weekday.highAt < 18 && bridge.week.length === 168);
        assert('the areas kept apart: Melbourne\'s tides hold no New York sensor, and New York\'s none of Melbourne\'s', tides.every((t) => t.id !== bridge.id) && tides.length > 50);
        const cur = ask({ op: 'currents', area: 'new-york' });
        assert('New York: the bridge a current meter each way, toward Manhattan and toward Brooklyn, the residual the mean of the net current',
          cur.sensors.length === 1 && cur.sensors[0].towardNames.join() === 'toward Manhattan,toward Brooklyn' && Math.abs(cur.sensors[0].net.reduce((a, v) => a + v, 0) / 168 - cur.sensors[0].residual) < 0.06);
        const july = ask({ op: 'sea', area: 'new-york', date: '2018-07-04' });
        const s = july.sensors[0];
        assert('New York: on 4 July 2018 its holiday and its weather, each way summing to the count', july.holiday === 'Independence Day' && july.weather && july.weather.tmax > 20
          && s.observed.every((v, hr) => !Number.isFinite(v) || s.way1[hr] + s.way2[hr] === v) && s.observed.some(Number.isFinite), `${july.holiday}, ${JSON.stringify(july.weather)}`);
      } else skip('New York\'s street records are not built (python data/products/street_records.py --only new-york)');
    }
  }
}

// ---- the engine of Event Spacetime, when it is built ------------------------------------------------
{
  const status = spawnSync('python', [join(APP, 'engine.py'), 'status'], { cwd: ROOT, encoding: 'utf8' });
  let state = null;
  try { state = JSON.parse(status.stdout).find((s) => s.plugin === 'event-spacetime'); } catch { /* not built */ }
  if (!state || !state.current) {
    skip('the engine of Event Spacetime is not built from the present sources (python applications/urban-science-workbench/engine.py build event-spacetime)');
  } else {
    const ask = (request) => {
      const r = spawnSync('python', [join(APP, 'engine.py'), 'ask', 'event-spacetime', JSON.stringify(request)],
        { cwd: ROOT, encoding: 'utf8', maxBuffer: 1 << 28, env: { ...process.env, PYTHONIOENCODING: 'utf-8' } });
      const answer = JSON.parse(r.stdout);
      if (!answer.ok) throw new Error(answer.error);
      return answer.result;
    };
    // One street between two junctions, two buildings side by side north of it and one south; a record of 2002–2024
    // in which the third is built in 2010, the first refurbished in 2015, the street has an event in 2016 and every
    // building is counted at the censuses of 2004 and 2012.
    const junctions = [[-40, 0], [40, 0]];
    const segments = [{ u: 0, v: 1, line: [[-40, 0], [40, 0]], kind: 'street', name: 'Main Street' }];
    const buildings = [
      { key: 'a', name: 'A', ring: [[-20, 5], [0, 5], [0, 25], [-20, 25]], base: 0, top: 12 },
      { key: 'b', name: 'B', ring: [[0, 5], [20, 5], [20, 25], [0, 25]], base: 0, top: 18 },
      { key: 'c', name: 'C', ring: [[-10, -25], [10, -25], [10, -6], [-10, -6]], base: 0, top: 9 },
    ];
    const table = (columns, rows) => ({ columns, rows });
    const record = {
      window: { from: '2002-01-01', to: '2025-01-01' }, sources: [{ id: 'a check' }],
      properties: table(['id', 'building'], [['base:1', 'a'], ['base:2', 'b'], ['base:3', 'c']]),
      censuses: table(['year', 'property'], [[2004, 'base:1'], [2004, 'base:2'], [2004, 'base:3'], [2012, 'base:1'], [2012, 'base:2'], [2012, 'base:3']]),
      events: table(['id', 'type', 'kind', 'from', 'to', 'match', 'source', 'provenance', 'detail'], [
        ['c built', 'birth', 'constructed', '2010-01-01', '2011-01-01', { property: 'base:3' }, 'a check', 'recorded', {}],
        ['a refurbished', 'alteration', 'refurbished', '2015-01-01', '2016-01-01', { property: 'base:1' }, 'a check', 'recorded', {}],
        ['a parade', 'regime', 'Public Event - Parade', '2016-03-05', '2016-03-06', { streets: ['Main Street'] }, 'a check', 'recorded', {}],
        ['a closure', 'regime', 'closed', '2017-07-01', '2017-07-02', { streets: ['Main Street'] }, 'a check', 'recorded', { closure: true }],
      ]),
    };
    const R = ask({ op: 'record', area: 'check', junctions, segments, buildings, record });
    const byId = Object.fromEntries(R.couplings.map((c, i) => [c.id, i]));
    assert('event spacetime: three buildings, a street space and two junction spaces; a party wall, three frontages, two cross-sections',
      R.network.buildings === 3 && R.network.streets === 1 && R.network.junctions === 2 && R.network.partyWalls === 1 && R.network.frontages === 3 && R.network.crossSections === 2,
      JSON.stringify(R.network));
    const C = R.couplings[byId['building:c']];
    assert('event spacetime: C begins at its construction in 2010, A and B stand from the window\'s start', C.life[0] === 55197 && R.couplings[byId['building:a']].life[0] === R.window.from);
    const frontC = R.primitives.find((p) => p.type === 'frontage' && R.couplings[p.a].id === 'building:c');
    assert('event spacetime: C\'s frontage begins with C', frontC && frontC.life[0] === C.life[0]);
    const kinds = Object.keys(R.kinds).sort().join(', ');
    assert('event spacetime: the events kept — the birth, the refurbishment, the parade and its end, the closure and its reopening, six censuses less the one before C stood',
      R.kinds['birth/constructed'] === 1 && R.kinds['alteration/refurbished'] === 1 && R.kinds['regime/Public Event - Parade'] === 1 && R.kinds['regime/ended'] === 1
        && R.kinds['regime/closed'] === 1 && R.kinds['regime/reopened'] === 1 && R.kinds['observation/census'] === 5 && R.dropped["a census before the building's construction"] === 1, kinds);
    const parade = R.events.find((e) => e.kind === 'Public Event - Parade'), closure = R.events.find((e) => e.kind === 'closed');
    assert('event spacetime: the parade is a state of the street space (it meets no interface); the closure meets the street\'s two cross-sections, its junctions attending',
      parade.meets.length === 0 && parade.attend.length === 1 && closure.meets.length === 2 && closure.attend.length === 3);
    const K = R.topology.K;
    assert('event spacetime: the event complex valid, connected, its Euler characteristic V − E + F, its Betti numbers by elimination',
      K.b0 === 1 && K.euler === K.vertices - K.edges + K.faces && K.b1 !== undefined && K.b2 !== undefined && K.b0 - K.b1 + K.b2 === K.euler, JSON.stringify(K));
    const named = { area: 'check', junctions, segments, buildings, record };
    const before = ask({ op: 'cut', t: 54000, ...named }), after = ask({ op: 'cut', t: 58000, ...named });
    assert('event spacetime: the network of the moment in 2006 holds five volumes and five interfaces, in 2017 six and six',
      before.couplings.length === 5 && before.primitives.length === 5 && after.couplings.length === 6 && after.primitives.length === 6,
      `${before.couplings.length}/${before.primitives.length} → ${after.couplings.length}/${after.primitives.length}`);
  }
}

// ---- the engine of Human flow, when it is built --------------------------------------------------------
{
  const status = spawnSync('python', [join(APP, 'engine.py'), 'status'], { cwd: ROOT, encoding: 'utf8' });
  let state = null;
  try { state = JSON.parse(status.stdout).find((s) => s.plugin === 'human-flow'); } catch { /* not built */ }
  const products = ['city-model-melbourne.json', 'street-sensors.json', 'street-holidays-victoria.json', 'street-terrain-melbourne.json', 'city-events-melbourne.json'];
  if (!state || !state.current) {
    skip('the engine of Human flow is not built from the present sources (python applications/urban-science-workbench/engine.py build human-flow)');
  } else if (!products.every((p) => existsSync(join(ROOT, 'data/cache/products', p)))) {
    skip('Human flow: the data library\'s city model and street records are not built (the recipe)');
  } else {
    // The district the page sends (src/district.js on the city model, the sensors placed in the area's frame with their
    // nearest way), for Melbourne.
    const { districtOf, boxOf, networkOf, footprintsIn } = await import(new URL(`file:///${join(APP, 'src/district.js').replace(/\\/g, '/')}`));
    const { localFrameAt } = await import(new URL(`file:///${join(ROOT, 'visualization/js/gis/frame.js').replace(/\\/g, '/')}`));
    const area = JSON.parse(readFileSync(join(APP, 'areas.json'), 'utf8')).areas.find((a) => a.id === 'melbourne');
    const model = JSON.parse(readFileSync(join(ROOT, 'data/cache/products/city-model-melbourne.json'), 'utf8'));
    const frame = localFrameAt(model.origin.lat, model.origin.lon);
    const B = boxOf(frame, districtOf(area));
    const net = networkOf(model, B);
    const WALK = new Set(['primary', 'secondary', 'tertiary', 'residential', 'unclassified', 'living_street', 'pedestrian', 'service', 'primary_link', 'secondary_link', 'tertiary_link']);
    const roads = model.roads.filter((r) => WALK.has(r.kind) && r.line && r.line.length >= 2);
    const sensors = [];
    for (const s of JSON.parse(readFileSync(join(ROOT, 'data/cache/products/street-sensors.json'), 'utf8')).sensors) {
      if (!Number.isFinite(s.lat) || !s.lat) continue;
      const p = frame.toLocal(s.lat, s.lon), at = [+p[0].toFixed(2), +p[1].toFixed(2)];
      if (!(at[0] > B.x0 - 500 && at[0] < B.x1 + 500 && at[1] > B.y0 - 500 && at[1] < B.y1 + 500)) continue;
      let best = [Infinity, null, ''];
      for (const r of roads) for (let i = 0; i + 1 < r.line.length; i++) {
        const a = r.line[i], b = r.line[i + 1], ex = b[0] - a[0], ey = b[1] - a[1], L2 = ex * ex + ey * ey;
        const t = L2 ? Math.max(0, Math.min(1, ((at[0] - a[0]) * ex + (at[1] - a[1]) * ey) / L2)) : 0;
        const d = Math.hypot(at[0] - a[0] - t * ex, at[1] - a[1] - t * ey);
        if (d < best[0]) best = [d, r.kind, r.name || ''];
      }
      sensors.push({ id: s.id, at, way: { name: best[2], kind: best[1] === 'service' ? 'laneway' : best[1] === 'pedestrian' || best[1] === 'living_street' ? 'mall' : 'street' } });
    }
    const district = { box: { x0: B.x0, y0: B.y0, x1: B.x1, y1: B.y1 }, ...net, buildings: footprintsIn(model, B).map((b) => ({ key: b.key, ring: b.ring, base: b.base, top: b.top })), sensors };
    // One engine process for every request (its runs are jobs it keeps), polled until each is done.
    const dir = mkdtempSync(join(tmpdir(), 'urban-hf-'));
    writeFileSync(join(dir, 'district.json'), JSON.stringify(district));
    const driver = [
      'import json, sys, time', `sys.path.insert(0, ${JSON.stringify(APP)})`, 'import engine',
      'D = json.load(open(sys.argv[1], encoding="utf-8"))', 'w = engine.Worker("human-flow")', 'out = {}',
      'def ask(op, **kw):', '    a = w.ask({"op": op, **D, **kw})', '    if not a.get("ok"): raise SystemExit(op + ": " + str(a.get("error")))', '    return a["result"]',
      'def wait(op, **kw):', '    while True:', '        r = ask(op, **kw)', '        if not r.get("pending"): return r', '        time.sleep(1)',
      'base = {"weekday": 2, "from": 12, "to": 14}',
      'opened = {**base, "regime": {"street": "Little Collins Street", "kind": "open", "from": 13, "to": 14}}',
      'closed = {**base, "regime": {"street": "Little Collins Street", "kind": "closed", "from": 13, "to": 14}}',
      'event = {**base, "event": {"street": "Bourke Street Mall", "crowd": 3000, "from": 13, "to": 14, "dwell": 30}}',
      'out["prepare"] = ask("prepare")', 'out["base"] = wait("result", **base)', 'out["frame"] = ask("frame", hour=13.5, **base)',
      'out["closed"] = wait("result", **closed)', 'out["opened"] = wait("result", **opened)', 'out["compare"] = ask("compare", a=opened, b=closed, over=13, until=14)',
      'out["event"] = wait("result", **event)', 'out["record"] = wait("record", mjd0=61000.5, **event)', 'out["permits"] = ask("permits")',
      'out["advice"] = wait("advice", epochs=1, epochMinutes=5, simulations=4, truth=1.5, **event)',
      'w.stop()', 'print(json.dumps(out))',
    ].join('\n');
    const r = spawnSync('python', ['-c', driver, join(dir, 'district.json')], { cwd: ROOT, encoding: 'utf8', maxBuffer: 1 << 29, env: { ...process.env, PYTHONIOENCODING: 'utf-8' } });
    let O = null;
    try { O = JSON.parse(r.stdout); } catch { assert('human flow: the engine answers', false, (r.stderr || r.stdout || '').slice(-400)); }
    if (O) {
      const P = O.prepare;
      assert(`human flow: the district's walkway (${P.junctions.length} junctions, ${P.segments.length} segments), its walk sums converging (radius ${P.radius})`,
        P.segments.length >= 300 && P.radius < 1 && P.signals.length >= 10, `${P.signals.length} signalled junctions`);
      assert(`human flow: ${P.counters.length} counters each way and ${P.entrances.length} entrances; the departures calibrated to them each way (R² ${P.calibration.r2Way}), visits a tenth of the trips`,
        P.counters.length >= 10 && P.entrances.length >= 15 && P.calibration.r2Way > 0.6 && Math.abs(P.calibration.visitShare - 0.09) < 0.02);
      const S = O.base;
      let sim = 0, exp = 0;
      for (const c of S.counters) for (let w = 0; w < 2; w++) for (let hh = 1; hh < c.simulated[w].length; hh++) sim += c.simulated[w][hh], exp += c.expected[w][hh];
      assert(`human flow: the walkers' counters read what the calibrated flow expects (Σ ${Math.round(sim)} against ${Math.round(exp)}, after the warm-up hour), and every walker arrives`,
        Math.abs(sim / exp - 1) < 0.1 && S.totals.failed === 0 && S.totals.arrived > 0.9 * S.totals.departed, `${S.totals.arrived} of ${S.totals.departed}`);
      const F = O.frame, inBox = F.frames.every((fr) => { for (let i = 0; i < fr.length; i += 7) { const x = fr[i + 4] / 10, y = fr[i + 5] / 10; if (x < B.x0 - 30 || x > B.x1 + 30 || y < B.y0 - 30 || y > B.y1 + 30 || fr[i + 3] < 1 || fr[i + 3] > 4) return false; } return true; });
      assert(`human flow: a quarter's frames, seven a minute apart by ten seconds, ${F.frames[0].length / 7} walkers at 13:30 inside the district, each walking, waiting or standing`,
        F.frames.length === 7 && F.frames[0].length / 7 > 1000 && inBox && F.density.length === P.segments.length);
      const lc = P.streets.find((s) => s.name === 'Little Collins Street');
      const inner = lc.edges.filter((e) => { const g = P.segments[e]; return [g.u, g.v].every((j) => P.segments.filter((x) => x.u === j || x.v === j).length > 1); });
      const onIt = (v) => inner.reduce((a, e) => a + v[e], 0);
      // Against the street open with the same walkers' draws (its own walkers kept apart alike): those passing through
      // leave it; those beginning or ending on it, or on a lane reached only through it, keep it, and nobody is stranded.
      assert(`human flow: Little Collins Street closed to walkers passing through carries fewer walkers (${Math.round(onIt(O.compare.entriesA))} → ${Math.round(onIt(O.compare.entriesB))} entries in the hour, the street open with the same draws against closed), none left with no way`,
        onIt(O.compare.entriesB) < onIt(O.compare.entriesA) && O.closed.totals.failed === 0 && O.opened.totals.failed === 0);
      // The other streams draw alike with the crowd and without it, so the crowd's departures are the difference.
      const crowd = O.event.totals.departed - O.base.totals.departed;
      assert(`human flow: an event's crowd of 3,000 an hour for an hour comes (${crowd} departures), and its run's record is sound with every event simulated`,
        Math.abs(crowd - 3000) < 4 * Math.sqrt(3000) && O.record.simulated && O.record.problems.length === 0 && O.record.events.length > 0);
      assert(`human flow: ${O.permits.permits.length} of the City's permits on the district's streets, with their days' counts`, O.permits.permits.length > 100 && O.permits.permits.some((p) => p.counts.length > 0));
      assert('human flow: POMCP decides the street\'s regime on the walkers (one decision), beside each regime held throughout',
        O.advice.steps.length === 1 && ['open', 'one-way', 'closed'].includes(O.advice.steps[0].action) && O.advice.fixed.length === 3);
    }
  }
}

// ---- the engine of Landscape Dynamics, when it is built ----------------------------------------------
{
  const status = spawnSync('python', [join(APP, 'engine.py'), 'status'], { cwd: ROOT, encoding: 'utf8' });
  let state = null;
  try { state = JSON.parse(status.stdout).find((s) => s.plugin === 'landscape-dynamics'); } catch { /* not built */ }
  if (!state || !state.current) {
    skip('the engine of Landscape Dynamics is not built from the present sources (python applications/urban-science-workbench/engine.py build landscape-dynamics)');
  } else {
    const ask = (request) => {
      const r = spawnSync('python', [join(APP, 'engine.py'), 'ask', 'landscape-dynamics', JSON.stringify(request)],
        { cwd: ROOT, encoding: 'utf8', maxBuffer: 1 << 28, env: { ...process.env, PYTHONIOENCODING: 'utf-8' } });
      const answer = JSON.parse(r.stdout);
      if (!answer.ok) throw new Error(answer.error);
      return answer.result;
    };
    // A district of three by three blocks 100 m apart, a dead end in one of them (a hole), two stores facing the
    // middle block from either side.
    const junctions = [], segments = [], at = (i, j) => 4 * j + i;
    for (let j = 0; j < 4; j++) for (let i = 0; i < 4; i++) junctions.push([100 * i, 100 * j]);
    for (let j = 0; j < 4; j++) for (let i = 0; i < 3; i++) segments.push({ u: at(i, j), v: at(i + 1, j), line: [junctions[at(i, j)], junctions[at(i + 1, j)]], kind: 'street', name: `row ${j}` });
    for (let i = 0; i < 4; i++) for (let j = 0; j < 3; j++) segments.push({ u: at(i, j), v: at(i, j + 1), line: [junctions[at(i, j)], junctions[at(i, j + 1)]], kind: 'street', name: `column ${i}` });
    junctions.push([250, 250]);
    segments.push({ u: at(2, 2), v: 16, line: [junctions[at(2, 2)], [250, 250]], kind: 'laneway', name: 'a dead end' });
    const net = { junctions, segments, planarise: true, stores: [{ name: 'A', use: 'other-retail', at: [150, 104] }, { name: 'B', use: 'food-and-drink', at: [150, 196] }] };
    const D = ask({ op: 'district', ...net });
    assert('landscape: the district a regular complex of 17 junctions, 25 segments and 8 blocks, one hole (b = 1, 1, 0)',
      D.topology.regular && D.topology.counts.join() === '17,25,8' && D.topology.betti.join() === '1,1,0', JSON.stringify(D.topology.counts));
    const L = ask({ op: 'lane', ...net });
    assert('landscape: the lane in three rewrites, two refinements then the split; closed, the same district, and planned again the same',
      L.steps.length === 3 && L.steps[0].refinement && L.steps[1].refinement && !L.steps[2].refinement && L.undo.isomorphic && L.undo.again.length === L.lane.length);
    const edgeAt = (a, b) => D.cells.edges.find((e) => { const p = e.line[0], q = e.line[e.line.length - 1]; return (p[0] === a[0] && p[1] === a[1] && q[0] === b[0] && q[1] === b[1]) || (p[0] === b[0] && p[1] === b[1] && q[0] === a[0] && q[1] === a[1]); });
    const X = ask({ op: 'complex', ...net, edits: [{ kind: 'cut', segment: edgeAt([100, 0], [100, 100]).id, at: [100, 50] }, { kind: 'join', segment: edgeAt([200, 0], [200, 100]).id },
      { kind: 'join', segment: edgeAt([0, 0], [100, 0]).id }] });
    assert('landscape: free edits — a street cut (a refinement), two blocks joined (an intervention), a join across the edge refused with its reason; the complex regular',
      X.edits.length === 2 && X.edits[0].refinement && !X.edits[1].refinement && X.refusedAt === 2 && /edge/.test(X.refused) && X.topology.regular && X.topology.counts.join() === '18,25,7',
      `${X.topology.counts} ${X.refused}`);
    const entries = Array.from({ length: 24 }, (_, hr) => 200 + 1400 * Math.exp(-0.5 * ((hr - 13) / 4) ** 2));
    const S0 = ask({ op: 'run', ...net, edits: [{ kind: 'cut', segment: edgeAt([100, 0], [100, 100]).id, at: [100, 50] }], params: { seed: 3, start: 12, end: 15 },
      initial: { walkers: { [D.cells.vertices[5].id]: 300 }, openA: true }, entries, altitude: Array(24).fill(0.5) });
    assert('landscape: a run on the edited complex from noon, from the state set: frames from 12:00, the first store open, the walkers set there and kept',
      S0.frames[0].time === 12 && S0.frames[0].a && Math.abs(S0.balance.start - 300) < 1e-6 && Math.abs(S0.balance.residual) < 1e-6);
    const R = ask({ op: 'run', ...net, params: { seed: 3 }, entries, altitude: Array(24).fill(0.5) });
    const eachWay = R.frames.every((F) => Object.keys(F.flow).every((id) => Math.abs(F.forth[id] - F.back[id] - F.flow[id]) <= 0.02 + 1e-6 * (F.forth[id] + F.back[id])));
    assert(`landscape: a day of ${R.events.length} events, walkers kept, every configuration regular, each way's hops differing by the net flow`,
      Math.abs(R.balance.residual) < 1e-6 && R.strata.every((s) => s.regular) && R.frames.length === 96 && eachWay, String(R.balance.residual));
    const S = ask({ op: 'skeleton' });
    const count = S.views.find((v) => /how many/.test(v.name));
    assert('landscape: the event skeleton, 16 states composed either way alike, driven; counting the open stores exact as a whole, not change by change',
      S.states.length === 16 && S.coherent && !S.detailedBalance && S.entropyProduction > 0 && count.whole && !count.eventwise);
    const shops = [{ at: [150, 104], stay: 0.3 }, { at: [157, 104], stay: 0.3 }, { at: [104, 150], stay: 0.8 }];
    const C = ask({ op: 'scales', ...net, shops });
    assert('landscape: as many flows round the holes as holes; walkers kept with the shops; with the stays forgotten, all on the streets',
      C.harmonic.length === C.betti[1] && C.memory.kept === 0 && C.memory.onStreetMarkov.every((v) => Math.abs(v - 100) < 1e-6) && C.memory.inside > 0 && C.operators.deviation === 0);
    const sun = { altitude: Array.from({ length: 24 }, (_, hr) => Math.max(-0.3, Math.sin((Math.PI * (hr - 6)) / 14))), azimuth: Array.from({ length: 24 }, (_, hr) => (Math.PI * (hr - 6)) / 12) };
    const Ms = ask({ op: 'massing', buildings: [{ key: 'a', ring: [[0, 0], [10, 0], [10, 10], [0, 10]], base: 0, top: 10 }, { key: 'b', ring: [[10, 0], [20, 0], [20, 10], [10, 10]], base: 0, top: 20 },
      { key: 'c', ring: [[40, 0], [50, 0], [50, 10], [40, 10]], base: 0, top: 10 }], ...sun, air: Array(24).fill(18), start: 0, end: 6, initial: { a: 30 } });
    const T6 = Ms.frames[Ms.frames.length - 1].temperature;
    assert('landscape: the buildings as one 3-complex — three 3-cells, one wall shared, two terraces; heat from the warmed building into its neighbour, not into the one apart',
      Ms.regular && Ms.counts[3] === 3 && Ms.walls.length === 1 && Ms.terraces === 2 && T6.b > T6.c && T6.a < 30, JSON.stringify(T6));
    const B = ask({ op: 'building', x: 0, y: 0, heading: 0, width: 20, depth: 12, storey: 3.5, nx: 4, ny: 2, floors: 3, ...sun, events: [{ hour: 9, kind: 'doors' }, { hour: 12, kind: 'merge' }, { hour: 16, kind: 'partition' }] });
    assert('landscape: a building of 24 rooms a ball; doors, a merge and a parting each keep its heat and its complex regular',
      B.counts[3] === 24 && B.betti.join() === '1,0,0,0' && B.events.length === 3 && B.events.every((e) => e.regular && Math.abs(e.heatChange) <= 1e-9 * Math.abs(B.hours[0].heat)));
  }
}

// The territory: a hexagon over an area in dispute takes no economy; one beside it keeps its own.
{
  const driver = [
    'import json, sys', `sys.path.insert(0, ${JSON.stringify(APP)})`, 'import engine',
    'if not engine.status("territory")["current"]: print(json.dumps({"skip": "the engine of The territory is not built from the present sources"})); raise SystemExit',
    'w = engine.Worker("territory")',
    'places = {"Srinagar": [34.08, 74.80], "Simferopol": [44.95, 34.10], "Laayoune": [27.15, -13.20], "Delhi": [28.61, 77.21], "Paris": [48.85, 2.35], "Tokyo": [35.68, 139.69]}',
    'out = {}',
    'for res in (2, 3):',
    '    a = w.ask({"op": "earth_people", "res": res})',
    '    if not a.get("ok"): print(json.dumps({"skip": str(a.get("error"))})); raise SystemExit',
    '    a = a["result"]',
    '    at = dict(zip(a["cells"], a["code"]))',
    '    for name, (lat, lon) in places.items():',
    '        c = w.ask({"op": "bin_points", "res": res, "points": [[lat, lon, 1]]})["result"]["cells"][0]',
    '        out[f"{name}:{res}"] = at.get(c, "no cell")',
    'w.stop()', 'print(json.dumps(out))',
  ].join('\n');
  const r = spawnSync('python', ['-c', driver], { cwd: ROOT, encoding: 'utf8', maxBuffer: 1 << 29, env: { ...process.env, PYTHONIOENCODING: 'utf-8' } });
  let O = null;
  try { O = JSON.parse(r.stdout); } catch { assert('the territory: the engine answers', false, (r.stderr || r.stdout || '').slice(-400)); }
  if (O && O.skip) skip(`the territory: ${O.skip}`);
  else if (O) {
    assert('the territory: no economy on a cell over an area in dispute', [2, 3].every((k) => ['Srinagar', 'Simferopol', 'Laayoune'].every((n) => O[`${n}:${k}`] === null)), JSON.stringify(O));
    assert('the territory: an economy on the cells of others, at their most people', O['Paris:3'] === 'FRA' && O['Tokyo:3'] === 'JPN' && O['Delhi:3'] === 'IND', JSON.stringify(O));
  }
}

done();
