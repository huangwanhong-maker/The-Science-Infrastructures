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
import { diverging, sequential, seriesColour } from '../src/palette.js';

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
assert("colours: the diverging ramp runs from sky blue through the neutral to rose (the day's, outside a page)",
  diverging(-1) === 'rgb(74, 134, 197)' && diverging(1) === 'rgb(224, 96, 126)' && diverging(0) === 'rgb(241, 240, 234)', `${diverging(-1)} ${diverging(0)} ${diverging(1)}`);
assert('colours: the sequential ramp runs from the neutral through gold to deep amber', sequential(0) === 'rgb(241, 240, 234)' && sequential(Math.pow(0.5, 1.25)) === 'rgb(221, 182, 90)'
  && sequential(1) === 'rgb(168, 112, 26)', `${sequential(0)} ${sequential(Math.pow(0.5, 1.25))} ${sequential(1)}`);
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
    toolbox(ask, sidePath);
  }
}

// ---- the toolbox (ADR 0048): each tool's answer against a closed form or an identity ------------------------
function toolbox(ask, sidePath) {
  const figure = (answer, title, label) => answer.blocks?.find((b) => b.kind === 'figures' && b.title.startsWith(title))?.items.find(([k]) => k === label)?.[1] ?? NaN;
  const solow = ask({ op: 'tool', tool: 'solow', params: { saving: 0.2, labour_growth: 0.01, technology_growth: 0.02, depreciation: 0.05, alpha: 1 / 3, sigma: 1 } });
  check('toolbox: Solow, k* = (s/(n+g+δ))^(1/(1−α))', figure(solow, 'The steady state', 'capital k*') - Math.pow(0.2 / 0.08, 1.5), 1e-9);
  check('toolbox: ... the golden rule saves α', figure(solow, 'The steady state', 'golden-rule saving') - 1 / 3, 1e-9);
  const islm = ask({ op: 'tool', tool: 'islm', params: { c0: 200, c1: 0.75, taxes: 100, i0: 200, i1: 25, government: 100, money: 1000, price: 2, k: 0.5, h: 50, delta_government: 50 } });
  check('toolbox: IS–LM, Y* = 1350 and r* = 3.5', Math.abs(figure(islm, 'The equilibrium', 'output Y*') - 1350) + Math.abs(figure(islm, 'The equilibrium', 'interest rate r*') - 3.5), 1e-9);
  check('toolbox: ... the fiscal multiplier h/(h(1−c₁)+i₁k) = 2, and fifty more spent raise Y by 100', Math.abs(figure(islm, 'The equilibrium', 'fiscal multiplier dY/dG') - 2)
    + Math.abs(figure(islm, 'After the policy', 'change in output') - 100), 1e-9);
  const t = [2015, 2016, 2017, 2018, 2019, 2020], Y = [100, 104, 107, 112, 111, 118], K = [300, 310, 322, 330, 341, 350], L = [50, 51, 51.5, 52, 52.2, 53];
  const ga = ask({ op: 'tool', tool: 'growth-accounting', columns: { year: t, output: Y, capital: K, labour: L }, roles: { time: 'year', output: 'output', capital: 'capital', labour: 'labour' },
    params: { capital_share: 0.3 } });
  const rows = ga.blocks?.find((b) => b.kind === 'table')?.rows || [];
  check('toolbox: growth accounts add up, output = capital + labour + TFP, period by period',
    rows.length === 5 ? Math.max(...rows.map((r) => Math.abs(r[1] - r[2] - r[3] - r[4]))) : NaN, 1e-9);

  // tyche and nomos through the engine.
  const tableOf = (answer, title) => answer.blocks?.find((b) => b.kind === 'table' && b.title.startsWith(title));
  const xs = [1, 2, 3, 4, 5, 6, 7, 8], zs = [3, 1, 4, 1, 5, 9, 2, 6];
  const ols = ask({ op: 'tool', tool: 'ols', columns: { y: xs.map((x, i) => 1 + 2 * x - 0.5 * zs[i]), x: xs, z: zs }, roles: { y: 'y', x: ['x', 'z'] }, params: { covariance: 'HC1', logs: 'none', constant: true } });
  const beta = (tableOf(ols, 'Coefficients')?.rows || []).map((r) => r[1]);
  check('toolbox: least squares on exact data returns (1, 2, −0.5)', beta.length === 3 ? Math.abs(beta[0] - 1) + Math.abs(beta[1] - 2) + Math.abs(beta[2] + 0.5) : NaN, 1e-10);
  const y4 = [10, 11, 12, 16, 9, 10.5, 11.5, 15.5];
  const did = ask({ op: 'tool', tool: 'did', columns: { y: y4, g: [0, 0, 1, 1, 0, 0, 1, 1], t: [0, 1, 0, 1, 0, 1, 0, 1] }, roles: { y: 'y', treated: 'g', post: 't' }, params: {} });
  check('toolbox: difference in differences is the difference of the four means: (15.75 − 11.75) − (10.75 − 9.5) = 2.75',
    figure(did, 'The effect', 'δ̂, difference in differences') - 2.75, 1e-12);
  const bos = ask({ op: 'tool', tool: 'game', params: { payoffs: '2,1  0,0\n0,0  1,2' } });
  const eq = tableOf(bos, 'Nash equilibria')?.rows || [];
  assert("toolbox: the Battle of the Sexes has two pure equilibria and the mixed one (⅔, ⅓), (⅓, ⅔)",
    eq.length === 3 && eq.some((r) => r[0] === 'R1 0.667, R2 0.333' && r[1] === 'C1 0.333, C2 0.667'), JSON.stringify(eq));
  const vote = ask({ op: 'tool', tool: 'voting', params: { ballots: '8: A > B > C > D\n6: B > C > D > A\n5: C > D > B > A\n4: D > C > B > A' } });
  const winners = Object.fromEntries(vote.blocks?.find((b) => b.kind === 'figures')?.items || []);
  assert('toolbox: on one profile, plurality elects A, instant runoff C, and B is the Condorcet winner',
    winners.plurality === 'A' && winners['instant runoff'] === 'C' && winners['the Condorcet winner'] === 'B', JSON.stringify(winners));
  const pi = ask({ op: 'tool', tool: 'price-index', columns: { period: ['1', '1', '2', '2'], item: ['a', 'b', 'a', 'b'], price: [1, 2, 1.5, 2.2], quantity: [10, 5, 8, 6] },
    roles: { period: 'period', item: 'item', price: 'price', quantity: 'quantity' }, params: { chained: false } });
  const last = tableOf(pi, 'The indices')?.rows?.[1] || [];
  // P_F · Q_F = Σp₁q₁ / Σp₀q₀ (the factor reversal test): (12 + 13.2) / (10 + 10).
  check("toolbox: Fisher's price index times his quantity index is the ratio of values", (last[3] * last[7]) / 1e4 - 25.2 / 20, 1e-12);
  const hp = ask({ op: 'tool', tool: 'hp-filter', columns: { t: [...Array(30).keys()], v: [...Array(30).keys()].map((k) => Math.exp(0.02 * k + 0.1 * Math.sin(k))) }, roles: { time: 't', value: 'v' }, params: { lambda: 1600, log: true } });
  const [trendBlock, cycleBlock] = hp.blocks?.filter((b) => b.kind === 'lines') || [];
  const [ser, trend] = trendBlock?.series || [];
  check('toolbox: the Hodrick–Prescott trend and cycle add up to the series', ser ? Math.max(...ser.y.map((v, k) => Math.abs(v - trend.y[k] - cycleBlock.series[0].y[k]))) : NaN, 1e-9);
  const trade = join(ROOT, 'data/cache/products/trade-flows.json');
  if (existsSync(trade)) {
    const g = ask({ op: 'tool', tool: 'trade-gravity', flows: sidePath(trade), params: { year: 2022, method: 'ppml' } });
    const elasticity = figure(g, 'The gravity', 'distance elasticity');
    assert('toolbox: the gravity of 2022 by PPML with effects, a distance elasticity between −1.2 and −0.6 (the literature: about −0.9)',
      elasticity > -1.2 && elasticity < -0.6, String(elasticity));
  }
}

// ---- one's own data, and the example plug-ins (ADR 0046): their own checks, one JSON object a line ----------
const pythons = [resolve(ROOT, '.venv/Scripts/python.exe'), resolve(ROOT, '.venv/bin/python'), 'python3', 'python'];
const python = pythons.find((p) => (p.includes('.venv') ? existsSync(p) : true) && spawnSync(p, ['-c', 'import numpy'], { stdio: 'ignore' }).status === 0);
for (const script of ['check_datasets.py', 'check_plugins.py']) {
  if (!python) { skip(`${script}: no Python with numpy`); continue; }
  const r = spawnSync(python, [join(APP, 'tests', script)], { cwd: ROOT, encoding: 'utf8', timeout: 900000, maxBuffer: 1 << 24, env: { ...process.env, PYTHONIOENCODING: 'utf-8' } });
  let seen = 0;
  for (const line of (r.stdout || '').split(/\r?\n/)) {
    if (!line.startsWith('{')) continue;
    let rec;
    try { rec = JSON.parse(line); } catch { continue; }
    seen++;
    if (rec.kind === 'check') check(rec.name, rec.value, rec.limit, rec.unit || '');
    else if (rec.kind === 'assert') assert(rec.name, rec.ok, rec.detail || '');
    else if (rec.kind === 'skip') skip(rec.name);
  }
  assert(`${script} ran to the end`, r.status === 0 && seen > 0, (r.stderr || '').trim().split(/\r?\n/).slice(-4).join(' | '));
}

done();
