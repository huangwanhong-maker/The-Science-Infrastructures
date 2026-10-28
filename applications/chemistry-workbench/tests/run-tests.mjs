// Checks of the Chemistry Workbench: its plug-ins' manifests (every view's module there), its
// recipe, the syntax of every page module, the drawing kit's reading of labels; and, when the engine
// of Symmetry and generation is built from the present sources, its answers: every check of the
// manuscript's claims passing; the CIP descriptors of the two configurations of CHBrClF
// against RDKit's (tests/fixtures/cip.json, made by make_fixtures.py); the point groups of the
// catalogue; the grammars' derivation counts; the complexes' homology; the Wick sums against the
// Gaussian moments; the field's values along the characters; Frank's network.
// The mathematics itself is checked in the foundations (morphe's graph rewriting, physis's point
// groups, stranded field theories and kinetics).
// Run from the repository root:  node applications/chemistry-workbench/tests/run-tests.mjs

import { existsSync, mkdtempSync, readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { spawnSync } from 'node:child_process';

import { suite, ROOT } from '../../../development/tools/checks.mjs';
import { parseLabel } from '../src/draw/graph.js';

const { assert, skip, done } = suite();
const APP = resolve(ROOT, 'applications/chemistry-workbench');

// ---- plug-ins ---------------------------------------------------------------------------------------
const plugins = readdirSync(join(APP, 'plugins')).filter((p) => existsSync(join(APP, 'plugins', p, 'plugin.json')));
assert('at least one plug-in', plugins.length > 0, plugins.join(', '));
for (const p of plugins) {
  const m = JSON.parse(readFileSync(join(APP, 'plugins', p, 'plugin.json'), 'utf8'));
  const missing = (m.views || []).filter((v) => !existsSync(join(APP, 'plugins', p, v.module)));
  assert(`plug-in ${p}: a title, views, and every view's module present`, !!m.title && (m.views || []).length > 0 && missing.length === 0, missing.map((v) => v.id).join(', '));
  if (m.engine) assert(`plug-in ${p}: its engine's source present`, existsSync(join(APP, 'plugins', p, m.engine)));
}
// ---- the recipe -------------------------------------------------------------------------------------
{
  const recipe = JSON.parse(readFileSync(join(APP, 'recipe.json'), 'utf8'));
  assert('the recipe: every step a script that exists', recipe.steps.every((s) => existsSync(join(ROOT, s.run[0]))));
}

// ---- the page's modules -----------------------------------------------------------------------------
{
  // node --check on a .js ES module may pass syntax errors; each module is checked as .mjs.
  const dir = mkdtempSync(join(tmpdir(), 'chemistry-check-'));
  const files = [];
  const walk = (d) => { for (const f of readdirSync(d, { withFileTypes: true })) { const p = join(d, f.name); if (f.isDirectory()) walk(p); else if (f.name.endsWith('.js')) files.push(p); } };
  walk(join(APP, 'src'));
  for (const p of plugins) if (existsSync(join(APP, 'plugins', p, 'views'))) walk(join(APP, 'plugins', p, 'views'));
  const bad = [];
  files.forEach((f, i) => {
    const copy = join(dir, `m${i}.mjs`);
    writeFileSync(copy, readFileSync(f, 'utf8'));
    const r = spawnSync(process.execPath, ['--check', copy], { encoding: 'utf8' });
    if (r.status !== 0) bad.push(`${f.slice(APP.length + 1)}: ${r.stderr.split('\n').find((l) => l.includes('Error')) || r.stderr}`);
  });
  assert(`the page's ${files.length} modules parse`, bad.length === 0, bad.join('; '));
}
{
  const c = parseLabel('C@+'), o = parseLabel('O-'), n = parseLabel('A_tet'), f = parseLabel('C•');
  assert('labels: an oriented carbon, a charged oxygen, a nonterminal, a free terminus',
    c.element === 'C' && c.tau === '+' && o.element === 'O' && o.charge === '-' && n.nonterminal && !f.element && f.free);
}

// ---- the engine, when it is built -------------------------------------------------------------------
{
  const status = spawnSync('python', [join(APP, 'engine.py'), 'status'], { cwd: ROOT, encoding: 'utf8' });
  let state = null;
  try { state = JSON.parse(status.stdout).find((s) => s.plugin === 'symmetry-generation'); } catch { /* not built */ }
  if (!state || !state.current) {
    skip('the engine of Symmetry and generation is not built from the present sources (python applications/chemistry-workbench/engine.py build)');
  } else {
    const ask = (request) => {
      const r = spawnSync('python', [join(APP, 'engine.py'), 'ask', 'symmetry-generation', JSON.stringify(request)],
        { cwd: ROOT, encoding: 'utf8', maxBuffer: 1 << 28, env: { ...process.env, PYTHONIOENCODING: 'utf-8' } });
      const answer = JSON.parse(r.stdout);
      if (!answer.ok) throw new Error(answer.error);
      return answer.result;
    };
    const results = ask({ op: 'checks' });
    const all = Object.values(results).flat();
    const failed = Object.entries(results).flatMap(([id, cs]) => cs.filter((c) => !c.ok).map((c) => `${id}: ${c.name} (${c.detail})`));
    assert(`the manuscript's claims: all ${all.length} checks of ${Object.keys(results).length} claims pass`, failed.length === 0, failed.join('; '));

    const fixture = JSON.parse(readFileSync(join(APP, 'tests/fixtures/cip.json'), 'utf8'));
    const geo = ask({ op: 'grammar', example: 'stereocentre' }).geometry;
    assert(`CHBrClF: τ and CIP against ${fixture.source.split(':')[0]} (X₊ ${fixture.plus.cip}, X₋ ${fixture.minus.cip})`,
      geo.tauPlus === fixture.plus.tau && geo.cipPlus === fixture.plus.cip && geo.tauMinus === fixture.minus.tau && geo.cipMinus === fixture.minus.cip,
      `engine: ${geo.tauPlus} ${geo.cipPlus}, ${geo.tauMinus} ${geo.cipMinus}`);

    const expected = { 'chbrclf-plus': 'C1', 'chbrclf-minus': 'C1', 'dichloroethene-z': 'C2v', 'dichloroethene-e': 'C2h', benzene: 'D6h', water: 'C2v', ammonia: 'C3v',
      methane: 'Td', ethylene: 'D2h', 'boron-trifluoride': 'D3h', 'sulfur-hexafluoride': 'Oh', ethane: 'D3d', 'hydrogen-peroxide': 'C2', allene: 'D2d', 'carbon-dioxide': 'D∞h' };
    const cat = ask({ op: 'catalogue' });
    const wrong = cat.filter((m) => expected[m.id] !== m.pointGroup.schoenflies).map((m) => `${m.id}: ${m.pointGroup.schoenflies}`);
    assert(`the catalogue's ${cat.length} idealized molecules have their textbook point groups (Cotton 1990)`, wrong.length === 0 && cat.length === Object.keys(expected).length, wrong.join(', '));
    const tables = cat.filter((m) => m.pointGroup.characters).every((m) => {
      const C = m.pointGroup.characters;
      return C.rows.length === C.classes.length && C.classes.reduce((s, c) => s + c.size, 0) === m.pointGroup.order;
    });
    assert('every character table is square and its classes partition the group', tables);

    const gly = ask({ op: 'grammar', example: 'glycine', units: 4 });
    assert('four labelled glycines: 24 backbones, 144 derivations, one orbit of S₄', gly.labelled.terminalClasses === 24 && gly.labelled.terminalDerivations === 144 && gly.labelled.orbits === 1);

    const cx = ask({ op: 'complexes' });
    assert('complexes: the six-cycle, filled, and the two-vertex complex have H = (ℤ, ℤ), (ℤ, 0, 0), (ℤ, 0, ℤ³)',
      cx.ring.homology.map((h) => h.group).join() === 'ℤ,ℤ' && cx.filled.homology.map((h) => h.group).join() === 'ℤ,0,0' && cx.dipole.homology.map((h) => h.group).join() === 'ℤ,0,ℤ^3');

    for (const [n, g] of [[1, 'Z2'], [2, 'Z3'], [3, 'Z2']]) {
      const w = ask({ op: 'wick', n, group: g });
      assert(`Wick, order ${n} over ${g}: Σ over pairings of |G|^F = brute-force moment (${w.sum})`, w.sum === w.moment);
    }
    const f = ask({ op: 'field', group: 'Z2' });
    assert('the field over ℤ₂: V along the characters is 64·[s₁ = s₃]', f.modes.every((m) => m.V === (m.s[0] === m.s[2] ? 64 : 0)));

    const fr = ask({ op: 'frank', precursor: 16 });
    assert('Frank: the generator commutes with L ↔ R and the law of outcomes is symmetric', fr.chain.commutes && fr.chain.asymmetry < 1e-12);
    assert('Frank: ṁ = (kᵢs/2) m (1 − m²) along the solution', fr.ode.identityError < 1e-9, String(fr.ode.identityError));
  }
}

done();
