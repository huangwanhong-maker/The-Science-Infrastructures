// Checks of the Urban Building Workbench: its plug-ins' manifests (every layer's and tool's module present;
// planned plug-ins with what they will hold), its areas and recipe, the syntax of every page module; and, when the engine of
// Street Oceanography is built from the present sources, its answers: the paper's Figure 1 read right by the
// amendment and wrongly by the paper in the fourth case; the catalytic pathway's coupling length against the
// killed walk's closed form; the kinematic street conserving its crowd; the district's detailed balance without
// circulation and its affinities untouched by attraction; and, when the data library's products are present,
// a sensor's rhythm (the day forced, harmonic constants resolved), the civil share of its clock, a street's
// hours and mixing (M = between + within); the tides (high water on a weekday, a day's prediction near its
// counts), the current meters (each way summing to the total, the handedness test), a street's crowd second
// by second (conserved, a red crossing passing nobody) and a visitor's plume (its masses summing to one).
// The mathematics itself is checked in the foundations and in the urban-science domain.
// Run from the repository root:  node applications/urban-building-workbench/tests/run-tests.mjs

import { existsSync, mkdtempSync, readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { spawnSync } from 'node:child_process';

import { suite, ROOT } from '../../../development/tools/checks.mjs';

const { assert, skip, done } = suite();
const APP = resolve(ROOT, 'applications/urban-building-workbench');

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
    skip('the engine of Street Oceanography is not built from the present sources (python applications/urban-building-workbench/engine.py build)');
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
    }
  }
}

done();
