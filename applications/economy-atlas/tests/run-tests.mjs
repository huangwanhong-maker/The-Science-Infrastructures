// Checks of what is the Economy Atlas's own: its examples (every leg between positions it has, every
// transport positive, every table square), its colours, its recipe (each step a manifest or script
// that exists); and, when the engine is built, its answers held to the protocol of the README (§3)
// and to the numbers of the user's paper and the textbook: the triangle's holonomy 0.60 and deficit
// 0.511, the five positions' chord holonomies 0.975 and 1.0725, the fee's deficit −ln 0.9, the market
// with prices flat with its log prices recovered, Miller and Blair's Leontief inverse; and, where
// the products are built, the world's trade in a year (exchange and net adding up to the whole, the
// shares of the net adding to 1, the attraction summing to zero, every circulating share between 0
// and 1) and the currencies' spacetime (one lens for each currency and month, its curvature measured
// against the exchange rates' alone). The science is checked with the economics domain and the
// foundations.
// Run from the repository root:  node applications/economy-atlas/tests/run-tests.mjs

import { existsSync, readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { spawnSync } from 'node:child_process';

import { suite, ROOT } from '../../../development/tools/checks.mjs';
import { NETWORKS, TABLES } from '../src/presets.js';
import { diverging, seriesColour } from '../src/palette.js';

const { assert, check, skip, done } = suite();
const APP = resolve(ROOT, 'applications/economy-atlas');

// ---- the examples ------------------------------------------------------------------------------------
for (const p of NETWORKS) {
  const n = p.network.positions.length;
  assert(`example ${p.id}: every leg joins positions it has, with a positive transport`,
    p.network.legs.every((l) => l.from >= 0 && l.from < n && l.to >= 0 && l.to < n && l.from !== l.to && (l.ratio ?? 1) > 0));
}
for (const t of TABLES) {
  const n = t.table.sectors.length;
  assert(`table ${t.id}: square deliveries and a final demand for each sector`,
    t.table.z.length === n && t.table.z.every((r) => r.length === n) && t.table.final_demand.length === n);
}
assert('colours: the diverging ramp runs from slate through paper to terracotta',
  diverging(-1) === 'rgb(63, 107, 124)' && diverging(1) === 'rgb(185, 97, 46)' && diverging(0) === 'rgb(243, 236, 223)', `${diverging(-1)} ${diverging(0)} ${diverging(1)}`);
assert('colours: the series cycles', seriesColour(0) === seriesColour(20) && seriesColour(-1) === seriesColour(19));

// ---- the recipe -------------------------------------------------------------------------------------
{
  const recipe = JSON.parse(readFileSync(join(APP, 'recipe.json'), 'utf8'));
  const missing = recipe.steps.filter((s) => {
    const [script, arg] = s.run;
    if (!existsSync(join(ROOT, script))) return true;
    return script.endsWith('fetch.py') && !existsSync(join(ROOT, 'data/registry', `${arg}.json`));
  });
  assert('the recipe: every step a script and, for fetching, a manifest that exist', missing.length === 0, missing.map((s) => s.id).join(', '));
}

// ---- the engine, when it is built --------------------------------------------------------------------
{
  const status = spawnSync('python', [join(APP, 'engine.py'), 'status'], { cwd: ROOT, encoding: 'utf8' });
  let state = null;
  try { state = JSON.parse(status.stdout); } catch { /* not built */ }
  const sidePath = (p) => (state?.toolchain === 'wsl' ? p.replace(/^([A-Za-z]):[\\/]/, (_, c) => `/mnt/${c.toLowerCase()}/`).replace(/\\/g, '/') : p);
  const ask = (request) => {
    const r = spawnSync('python', [join(APP, 'engine.py'), 'run', '-'], { cwd: ROOT, encoding: 'utf8', input: JSON.stringify(request), maxBuffer: 1 << 28 });
    try { return JSON.parse(r.stdout); } catch { return { error: (r.stdout || r.stderr || 'no answer').replace(/\0/g, '').trim().slice(0, 200) }; }
  };
  if (!state?.current) skip('the engine is not built from the present sources: its answers are not checked');
  else {
    const net = (id) => NETWORKS.find((p) => p.id === id).network;
    const tri = ask({ op: 'circuits', network: net('triangle') });
    check('engine: the paper\'s triangle has holonomy 0.60', (tri.cycles?.[0]?.holonomy ?? NaN) - 0.6, 1e-12);
    check('engine: ... and deficit −ln 0.6', (tri.cycles?.[0]?.deficit ?? NaN) + Math.log(0.6), 1e-12);
    check('engine: ... its retentions multiply to it', (tri.cycles?.[0]?.retention ?? [NaN]).reduce((s, r) => s * r, 1) - 0.6, 1e-12);
    const five = ask({ op: 'circuits', network: net('five') });
    const h = (five.cycles || []).map((c) => c.holonomy).sort((a, b) => a - b);
    check('engine: the paper\'s five positions, chord holonomies 0.975 and 1.0725', Math.abs(h[0] - 0.975) + Math.abs(h[1] - 1.0725), 1e-12);
    const fee = ask({ op: 'circuits', network: net('fee') });
    check('engine: a fee of a tenth round a circuit is a deficit of −ln 0.9', (fee.cycles?.[0]?.deficit ?? NaN) + Math.log(0.9), 1e-12);
    const prices = ask({ op: 'circuits', network: net('prices') });
    const p = [1, 2, 5, 8];
    const phi = prices.log_price || [];
    const recovered = p.every((x, i) => Math.abs((phi[i] - phi[0]) - (-Math.log(x) + Math.log(p[0]))) < 1e-12);
    assert('engine: a market with prices is flat, its potential minus the log prices', prices.defect < 1e-12 && recovered, JSON.stringify(phi));
    const mb = TABLES.find((t) => t.id === 'miller-blair').table;
    const x = mb.z.map((r, i) => r.reduce((s, v) => s + v, 0) + mb.final_demand[i]);
    const va = x.map((xj, j) => xj - mb.z.reduce((s, r) => s + r[j], 0));
    const pr = ask({ op: 'production', table: { ...mb, value_added: va } });
    const d = 0.7575, L = [[0.95 / d, 0.25 / d], [0.2 / d, 0.85 / d]];
    check('engine: Miller and Blair\'s Leontief inverse', pr.L ? Math.max(...L.flatMap((r, i) => r.map((v, j) => Math.abs(pr.L[i][j] - v)))) : NaN, 1e-12);
    // The products, where built.
    const trade = join(ROOT, 'data/cache/products/trade-flows.json');
    if (existsSync(trade)) {
      const m = ask({ op: 'map', flows: sidePath(trade), year: 2022 });
      const T = m.totals || {};
      assert('engine: the map of 2022, its economies, net legs and triangles', !m.error && m.economies?.length > 150 && m.legs?.length > 1000 && m.triangles?.length > 100, m.error || '');
      check('engine: ... exchange and net transfer add up to the whole', (T.reciprocal + T.net - T.total) / T.total, 1e-9);
      check('engine: ... the shares of the net flows add to 1', T.driven + T.local + T.global - 1, 1e-9);
      check('engine: ... the attraction sums to zero', m.economies ? m.economies.reduce((s, e) => s + e.attraction, 0) : NaN, 1e-4);
      assert("engine: ... every economy's circulating share between 0 and 1", m.economies?.every((e) => e.circulating >= 0 && e.circulating <= 1));
    } else skip('the trade flows are not built: the map is not checked');
    const cur = join(ROOT, 'data/cache/products/currency-spacetime.json');
    if (existsSync(cur)) {
      const s = ask({ op: 'spacetime', record: sidePath(cur), first: '2019-01', last: '2026-06' });
      assert('engine: the currencies\' spacetime, one lens for each currency and month',
        !s.error && s.lenses?.length === (s.currencies.length - 1) * (s.months.length - 1), s.error || `${s.lenses?.length}`);
      const K = s.curvature || {};
      assert("engine: ... its curvature against the exchange rates' alone, 0 < rms and the share removed below 1",
        K.rms > 0 && K.rms_exchange > 0 && K.removed < 1 && Math.abs(K.removed - (1 - (K.rms / K.rms_exchange) ** 2)) < 1e-9, JSON.stringify(K));
    } else skip('the currencies\' record is not built: the spacetime is not checked');
  }
}

done();
