// Checks of what is the Spiking Simulator's own: the quantities it shows and writes, the
// specification it edits (renaming and removing keep every reference whole), its examples
// (every model, parameter and receptor they name exists in the engine's catalogue), and, when
// the engine is built, every example run by the engine itself. The science is checked with
// the neuroscience domain.
// Run from the repository root:  node applications/spiking-simulator/tests/run-tests.mjs

import { spawnSync } from 'node:child_process';
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';

import { suite, ROOT, CACHE } from '../../../development/tools/checks.mjs';
import { format, parseField, wellFormed, seconds, hint } from '../src/units.js';
import { empty, addPopulation, renamePopulation, removePopulation, removeProjection, problems, forEngine, setRule, ruleOf, clone, normalise } from '../src/spec.js';
import { PRESETS } from '../src/presets.js';
import { POPULATION_COLOURS } from '../src/palette.js';
import { ticks, decades, scale, extent } from '../src/charts.js';
import { toSI } from '../src/views.js';

const { assert, skip, done } = suite();
const APP = resolve(ROOT, 'applications/spiking-simulator');

// ---- quantities ----------------------------------------------------------------------------------
{
  const cases = [[2.5e-10, 'F', '250 pF'], [0.02, 's', '20 ms'], [20000, 'Hz', '20 kHz'], [-0.07, 'V', '-70 mV'], [1.5, 's', '1.5 s'],
    [0.7e-6, 'kg^-2 m^-4 s^6 A^3', '0.7 nS/mV'], [3e-6, 's', '3 µs'], [0, 'A', '0 A'], [0.08, '1', '0.08'], [120, 's', '120 s']];
  const wrong = cases.filter(([v, u, want]) => format(v, u) !== want).map(([v, u, want]) => `${v} ${u}: ${format(v, u)}, not ${want}`);
  assert('simulator: values are shown with the prefix that suits them', !wrong.length, wrong.join('; '));
  assert('simulator: a prefix is written in ASCII for the engine (µ as u)', format(3e-6, 's', { ascii: true }) === '3 us');
  const parsed = [parseField('20 ms'), parseField('3'), parseField('{"uniform": ["0 mV", "20 mV"]}'), parseField(''), parseField('5 µs')];
  assert('simulator: a field is read as a quantity, a number, a distribution or nothing',
    parsed[0] === '20 ms' && parsed[1] === 3 && parsed[2].uniform[1] === '20 mV' && parsed[3] === null && parsed[4] === '5 us', JSON.stringify(parsed));
  assert('simulator: a bare number is refused where a unit is needed, and a word refused anywhere',
    !wellFormed('3', 's') && wellFormed('3', '1') && !wellFormed('fast', 's') && wellFormed('2.5e-3 s', 's') && !wellFormed('{bad', 's'));
  assert('simulator: durations are read in seconds', seconds('500 ms') === 0.5 && seconds('2 s') === 2 && seconds(1.5) === 1.5 && Number.isNaN(seconds('soon')));
  assert('simulator: the engine\'s compound units read as neuroscience writes them', hint('kg^-2 m^-4 s^6 A^3') === 'nS/mV' && hint('1') === 'a number');
}

// ---- the charts ------------------------------------------------------------------------------------------
{
  const same = (a, b) => a.length === b.length && a.every((v, i) => Math.abs(v - b[i]) < 1e-12);
  assert('simulator: ticks at 1, 2 or 5 times a power of ten', same(ticks(0, 1, 5), [0, 0.2, 0.4, 0.6, 0.8, 1]) && same(ticks(-3, 7, 5), [-2, 0, 2, 4, 6]) && same(ticks(0, 500, 5), [0, 100, 200, 300, 400, 500]), JSON.stringify(ticks(0, 1, 5)));
  assert('simulator: decades between two positive values', same(decades(1, 1000), [1, 10, 100, 1000]) && same(decades(0.003, 0.2), [0.01, 0.1]));
  const lin = scale([0, 10], [100, 300]), log = scale([1, 100], [0, 200], true);
  assert('simulator: scales map and invert, linearly and logarithmically', lin(5) === 200 && Math.abs(lin.invert(250) - 7.5) < 1e-12 && Math.abs(log(10) - 100) < 1e-12 && Math.abs(log.invert(200) - 100) < 1e-9);
  const e = extent([3, 1, NaN, 2], { zero: true });
  assert('simulator: an extent covers the finite values, and zero when asked', e[0] === 0 && e[1] > 3 && e[1] < 3.2, JSON.stringify(e));
  assert('simulator: quantities as specifications write them, in SI units', toSI('40 pA') === 4e-11 && Math.abs(toSI('3 nA') - 3e-9) < 1e-24 && toSI(5) === 5 && toSI('1.5 ms') === 1.5e-3 && Number.isNaN(toSI('fast')));
}

// ---- the specification ------------------------------------------------------------------------------
{
  const s = empty();
  assert('simulator: an empty network is refused before it reaches the engine', problems(s).some((p) => p.includes('no population')));
  addPopulation(s, 'lif_psc_exp');
  addPopulation(s, 'lif_psc_exp');
  assert('simulator: a new population takes a name not yet used', s.populations[0].name !== s.populations[1].name, s.populations.map((p) => p.name).join(', '));
  const [a, b] = s.populations.map((p) => p.name);
  s.projections.push({ source: a, target: b, receptor: 'excitatory', rule: { all_to_all: true } }, { source: b, target: a, receptor: 'excitatory', rule: { one_to_one: true } });
  s.devices.push({ model: 'dc', target: a });
  s.recorders.push({ model: 'spikes', population: b }, { model: 'weights', projection: 1 });
  renamePopulation(s, a, 'E');
  assert('simulator: renaming a population renames it in every projection, device and recorder',
    s.projections[0].source === 'E' && s.projections[1].target === 'E' && s.devices[0].target === 'E' && !problems(s).length, problems(s).join('; '));
  const t = clone(s);
  removeProjection(t, 0);
  assert('simulator: removing a projection renumbers the weight recorders after it', t.recorders[1].projection === 0 && !problems(t).length);
  removePopulation(s, 'E');
  assert('simulator: removing a population removes what names it, and what recorded its projections',
    s.projections.length === 0 && s.devices.length === 0 && s.recorders.length === 1 && !problems(s).length, JSON.stringify(s.recorders));
  const p = { rule: { fixed_indegree: 100, autapses: false } };
  setRule(p, 'fixed_probability', 0.2);
  assert('simulator: changing a rule keeps its options', ruleOf(p).kind === 'fixed_probability' && ruleOf(p).value === 0.2 && p.rule.autapses === false);
  const stripped = forEngine({ duration: '1 s', populations: [{ name: 'n', model: 'parrot', size: 1, parameters: {}, bias: null }], projections: [] });
  assert('simulator: empty fields are not sent to the engine', !('parameters' in stripped.populations[0]) && !('bias' in stripped.populations[0]));
  assert('simulator: eight population colours, Dusty Orchid first', POPULATION_COLOURS.length === 8 && POPULATION_COLOURS[0] === '#bf68ae');
}

// ---- the examples, against the engine's catalogue ------------------------------------------------------
const catalogueFile = join(CACHE, 'engines/snn/models.json');
const catalogue = existsSync(catalogueFile) ? new Map(JSON.parse(readFileSync(catalogueFile, 'utf8')).map((m) => [m.name, m])) : null;
{
  const ids = PRESETS.map((p) => p.id);
  assert('simulator: every example has its own name, a title and an account of what it shows',
    new Set(ids).size === ids.length && PRESETS.every((p) => p.title && p.about && p.about.length > 80), ids.join(', '));
  const wrong = PRESETS.flatMap((p) => problems(p.spec()).map((x) => `${p.id}: ${x}`));
  assert('simulator: every example is whole (its references, sizes and names)', !wrong.length, wrong.join('; '));
  if (!catalogue) skip('simulator: the examples against the catalogue: the engine is not built');
  else {
    const missing = [];
    for (const p of PRESETS) {
      const s = normalise(p.spec());
      const model = (name) => catalogue.get(s.populations.find((q) => q.name === name)?.model);
      for (const q of s.populations) {
        const m = catalogue.get(q.model);
        if (!m) { missing.push(`${p.id}: no model ${q.model}`); continue; }
        for (const k of Object.keys(q.parameters ?? {})) if (!m.parameters.some((x) => x.name === k)) missing.push(`${p.id}: ${q.model} has no parameter ${k}`);
        for (const k of Object.keys(q.initial ?? {})) if (!m.components.some((x) => x.name === k)) missing.push(`${p.id}: ${q.model} has no component ${k}`);
        if (Array.isArray(q.bias) && q.bias.length !== q.size) missing.push(`${p.id}: ${q.name}'s bias has ${q.bias.length} values for ${q.size} neurons`);
      }
      for (const x of [...s.projections, ...s.devices.filter((d) => d.receptor)]) {
        const m = model(x.target);
        if (m && !m.receptors.some((r) => r.name === x.receptor)) missing.push(`${p.id}: ${m.name} has no receptor ${x.receptor}`);
      }
      for (const r of s.recorders.filter((x) => x.components)) {
        const m = model(r.population);
        for (const c of r.components) if (m && !m.components.some((x) => x.name === c)) missing.push(`${p.id}: ${m.name} has no component ${c}`);
      }
    }
    assert('simulator: every model, parameter, component and receptor the examples name is in the engine\'s catalogue', !missing.length, missing.join('; '));
  }
}

// ---- the examples, run by the engine ---------------------------------------------------------------------------------
{
  const python = process.platform === 'win32' ? 'python' : 'python3';
  const status = spawnSync(python, [join(APP, 'engine.py'), 'status'], { encoding: 'utf8' });
  let current = false;
  try { current = JSON.parse(status.stdout).current; } catch { /* no engine */ }
  if (!current) skip('simulator: the examples run by the engine: the engine is not built from the present sources (python data/acquisition/recipe.py applications/spiking-simulator/recipe.json)');
  else {
    const dir = mkdtempSync(join(tmpdir(), 'spiking-simulator-'));
    try {
      let described = 0, profiled = 0;
      for (const p of PRESETS) {
        // A short run of each: enough for every population to fire.
        const s = forEngine(p.spec());
        s.duration = '300 ms';
        const file = join(dir, `${p.id}.json`);
        writeFileSync(file, JSON.stringify(s));
        const r = spawnSync(python, [join(APP, 'engine.py'), 'run', file, '--threads', '4'], { encoding: 'utf8', maxBuffer: 1 << 28 });
        let res = null;
        try { res = JSON.parse(r.stdout); } catch { /* reported below */ }
        const recorded = res ? res.spikes.reduce((n, x) => n + x.times.length, 0) : 0;
        assert(`simulator: the example "${p.title}" runs in the engine and its neurons fire`,
          r.status === 0 && res && Math.abs(res.time - 0.3) < 1e-9 && recorded > 0 && res.populations.length === s.populations.length,
          r.status === 0 ? `${recorded} recorded spikes in ${res?.statistics.run_seconds.toFixed(2)} s` : (r.stderr || '').trim().split('\n').pop());
        if (res) {
          const prof = res.statistics.profile;
          if (res.analysis.length === res.spikes.length && res.analysis.every((a) => a.rates.length === a.count && Number.isFinite(a.rate_mean))) described++;
          if (prof && prof.threads.length === 4 && prof.windows > 0 && prof.update + prof.deliver > 0) profiled++;
        }
      }
      assert('simulator: every run describes each spike recorder\'s activity, neuron by neuron', described === PRESETS.length, `${described} of ${PRESETS.length}`);
      assert('simulator: every run gives the engine\'s profile, thread by thread', profiled === PRESETS.length, `${profiled} of ${PRESETS.length}`);
    } finally {
      rmSync(dir, { recursive: true, force: true });
    }
  }
}

done();
