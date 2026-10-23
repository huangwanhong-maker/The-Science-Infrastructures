// Checks of what is the Complex Atlas's own: how it writes groups, words and presentations; what
// it reads as a complex; its examples; the fixtures it can show without the engine (each a whole
// complex, each analysis obeying the Euler-Poincaré formula), and those of the generation of
// spacetime (graphs weighted by N to their faces, a density of mass 1, a condensate classical at low
// density, spin foams whose vertices carry their faces' spins); its page and server; and, when the
// engine is built, every family of the catalogue built and analysed by the engine itself. The
// drawing is checked with the visualization layer; the mathematics with morphe and physis.
// Run from the repository root:  node applications/complex-atlas/tests/run-tests.mjs

import { spawnSync } from 'node:child_process';
import { existsSync, readdirSync, readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';

import { suite, ROOT } from '../../../development/tools/checks.mjs';
import { checkComplex, eulerCharacteristic, cellCounts } from '../../../visualization/js/index.js';
import { groupText, wordText, presentationText, normalFormNames, chainText, countsText, number } from '../src/format.js';
import { kindOf, readComplex } from '../src/input.js';
import { EXAMPLES } from '../src/examples.js';
import { selectionSets } from '../src/selection.js';

const { assert, check, skip, done } = suite();
const APP = resolve(ROOT, 'applications/complex-atlas');
const FIX = join(APP, 'tests/fixtures');
const read = (p) => JSON.parse(readFileSync(p, 'utf8'));

// ---- writing --------------------------------------------------------------------------------------
{
  const cases = [[[1, []], 'ℤ'], [[2, []], 'ℤ²'], [[0, [2]], 'ℤ/2'], [[1, [2, 3, 2]], 'ℤ ⊕ (ℤ/2)² ⊕ ℤ/3'], [[0, []], '0'], [[12, []], 'ℤ¹²']];
  const wrong = cases.filter(([[r, t], want]) => groupText(r, t) !== want).map(([[r, t], want]) => `${r} ${t}: ${groupText(r, t)}, not ${want}`);
  assert('atlas: abelian groups written as ℤ^r ⊕ ℤ/d', !wrong.length, wrong.join('; '));
  assert('atlas: words written with powers and inverses', wordText([1, 1, 2, -1], ['a', 'b']) === 'a² b a⁻¹' && wordText([-2, -2, -2], ['a', 'b']) === 'b⁻³' && wordText([], ['a']) === '1');
  assert('atlas: presentations as ⟨ generators | relators ⟩, the trivial group as ⟨ | ⟩',
    presentationText(['a', 'b'], [[1, 2, -1, -2]]) === '⟨ a, b | a b a⁻¹ b⁻¹ ⟩' && presentationText([], []) === '⟨ | ⟩' && presentationText(['a'], []) === '⟨ a | ⟩');
  assert('atlas: the letters of normal forms (a₁ b₁ a₂ b₂; a a for one crosscap)',
    normalFormNames({ orientable: true, normal_form: [1, 2, -1, -2, 3, 4, -3, -4] }).join() === 'a₁,b₁,a₂,b₂' && normalFormNames({ orientable: false, normal_form: [1, 1] }).join() === 'a');
  assert('atlas: chains written as sums of cells', chainText([[0, 1], [1, -2]], 1, ['a', 'b']) === 'a − 2 b' && chainText([[2, -1]], 2, []) === '−f3');
  assert('atlas: counts of cells in words', countsText([1, 2, 1]) === '1 vertex · 2 edges · 1 face' && countsText([5, 10, 10, 5]) === '5 vertices · 10 edges · 10 faces · 5 3-cells');
  assert('atlas: numbers for reading', number(0.69777) === '0.6978' && number(-3) === '−3' && number(null) === '—');
}

// ---- reading ---------------------------------------------------------------------------------------
{
  assert('atlas: text is read as JSON, facets or words by its look', kindOf('{"vertices": 1}') === 'json' && kindOf('0 1 2; 1 2 3') === 'facets' && kindOf('abAB') === 'words');
  const t = readComplex('abAB'), f = readComplex('0 1 2\n0 2 3', 'auto'), j = readComplex(JSON.stringify({ complex: { vertices: 1, edges: [[0, 0]], faces: [] } }));
  assert('atlas: a complex is read from words, facets, or JSON (bare or wrapped as the engine answers)',
    t.complex.edges.length === 2 && f.complex.faces.length === 2 && j.complex.edges.length === 1 && j.kind === 'json');
  let message = '';
  try { readComplex('{"vertices": 2, "edges": [[0, 1]], "faces": [{"base": 0, "boundary": [1]}]}'); } catch (error) { message = error.message; }
  assert('atlas: a complex that is not whole is refused with the reason', message.includes('does not close'), message);
  try { readComplex('ab#'); message = ''; } catch (error) { message = error.message; }
  assert('atlas: words that are not words are refused', message.includes('not a letter'), message);
}

// ---- the examples -----------------------------------------------------------------------------------
{
  const chi = { torus: 0, 'genus-2': -2, 'projective-plane': 1, klein: 0, mobius: 0, 'dunce-hat': 1, sphere: 2, 'two-faces': 0, tetrahedron: 2, octahedron: 2, 'rp2-6': 1, 'torus-7': 0, 'simplex-4': 0 };
  const wrong = EXAMPLES.flatMap((ex) => {
    try {
      const { complex } = readComplex(ex.text, ex.kind);
      return eulerCharacteristic(complex) === chi[ex.id] ? [] : [`${ex.id}: χ = ${eulerCharacteristic(complex)}`];
    } catch (error) { return [`${ex.id}: ${error.message}`]; }
  });
  assert('atlas: every example reads whole, with the Euler characteristic of its space', !wrong.length && EXAMPLES.length === Object.keys(chi).length, wrong.join('; '));
  assert('atlas: every example has its own name, a title and an account', new Set(EXAMPLES.map((e) => e.id)).size === EXAMPLES.length && EXAMPLES.every((e) => e.title && e.about.length > 30));
}

// ---- the fixtures ------------------------------------------------------------------------------------
{
  const files = readdirSync(FIX);
  const builds = files.filter((f) => f.startsWith('build-')).map((f) => [f, read(join(FIX, f)).complex]);
  const bad = builds.flatMap(([f, c]) => checkComplex(c).map((p) => `${f}: ${p}`));
  assert('atlas: every fixture complex is whole', builds.length >= 4 && !bad.length, bad.join('; '));
  const problems = [];
  for (const [f, c] of builds) {
    const file = join(FIX, f.replace('build-', 'analyse-'));
    if (!existsSync(file)) continue;
    const a = read(file).analysis;
    const ranks = a.homology.reduce((s, h) => s + (h.dimension % 2 ? -h.rank : h.rank), 0);
    if (a.euler !== eulerCharacteristic(c) || ranks !== a.euler || cellCounts(c).join() !== a.counts.join()) problems.push(`${f}: χ ${a.euler}, ranks ${ranks}, counts ${a.counts}`);
    for (const [p, b] of Object.entries(a.betti || {})) if (b.reduce((s, x, k) => s + (k % 2 ? -x : x), 0) !== a.euler) problems.push(`${f}: Betti numbers over F${p}`);
    // Each generator of the fundamental group is a closed path at the base vertex 0.
    for (const loop of a.fundamental_group?.loops || []) {
      let at = 0;
      for (const k of loop) { const [t, h] = c.edges[Math.abs(k) - 1]; if ((k > 0 ? t : h) !== at) { problems.push(`${f}: a loop is not a path`); break; } at = k > 0 ? h : t; }
      if (at !== 0) problems.push(`${f}: a loop does not close`);
    }
  }
  assert('atlas: every fixture analysis obeys Euler-Poincaré (Σ(−1)^k rank H_k = χ, over ℤ and each field), and its loops close', !problems.length, problems.join('; '));
  const lattice = read(join(FIX, 'build-lattice-torus.json')).complex, g = read(join(FIX, 'gauge-lattice-torus.json'));
  let worst = 0;
  lattice.faces.forEach((f, i) => {
    const a = f.boundary.reduce((s, k) => s + Math.sign(k) * g.link_angle[Math.abs(k) - 1], 0);
    worst = Math.max(worst, Math.abs(Math.cos(a) - g.plaquettes[i]));
  });
  assert('atlas: the fixture gauge field is a field on the lattice: each plaquette the holonomy of its links', g.link_angle.length === lattice.edges.length && g.plaquettes.length === lattice.faces.length && worst < 1e-5, `${worst}`);
  const marks = selectionSets({ kind: 'loop', word: [1, 2, -1, -2] }, lattice);
  assert('atlas: a loop marks its edges and their ends', marks.edges.size === 2 && marks.edges.get(0) === 1 && marks.vertices.size >= 2);
}

// ---- the generation of spacetime: theories, expansions, large N, condensates, geometries, spin foams ----
{
  const { colouredGraphComplex, colouredGraphLayout } = await import('../../../visualization/js/index.js');
  const files = new Set(readdirSync(FIX));
  const has = (name) => files.has(`${name}.json`);
  // The faces of colours i, j of a coloured graph: the cycles of sigma_j⁻¹ sigma_i.
  const faces = (cg, i, j) => {
    const inv = []; cg.sigma[j].forEach((w, b) => { inv[w] = b; });
    const seen = new Uint8Array(cg.p); let n = 0;
    for (let b0 = 0; b0 < cg.p; b0++) { if (seen[b0]) continue; n++; for (let b = b0; !seen[b]; b = inv[cg.sigma[i][b]]) seen[b] = 1; }
    return n;
  };
  const permutation = (s, p) => s.length === p && new Set(s).size === p && s.every((x) => Number.isInteger(x) && x >= 0 && x < p);
  const wholeGraph = (cg) => cg && Number.isInteger(cg.colours) && cg.sigma.length === cg.colours && cg.sigma.every((s) => permutation(s, cg.p));
  if (!has('theories')) skip('atlas: the generation fixtures: none yet');
  else {
    const T = read(join(FIX, 'theories.json')).theories;
    const ids = ['scalar', 'vector', 'matrix', 'tensor-melonic', 'tensor-wheel', 'coloured-simplicial', 'gft-boulatov', 'gft-ooguri', 'gft-melonic'];
    assert('atlas: the theories are the nine of the protocol, each with its action, family, rank and bubbles',
      ids.every((id) => T.some((t) => t.id === id)) && T.every((t) => t.title && t.action && ['QFT', 'GFT'].includes(t.family) && Number.isInteger(t.rank) && Array.isArray(t.params) && Array.isArray(t.bubbles) && t.bubbles.every((b) => b.name && wholeGraph(b.coloured_graph))),
      T.map((t) => t.id).join(', '));
    const byId = Object.fromEntries(T.map((t) => [t.id, t]));
    // Expansions: every graph a coloured graph whose N-power is its count of faces: of colours
    // (0, i) when the theory has bubbles (colour 0 the propagators), of all pairs when it has none.
    const problems = [];
    let expansions = 0;
    for (const f of files) {
      const m = f.match(/^expand-(.+)-(\d+)\.json$/);
      if (!m) continue;
      expansions++;
      const t = byId[m[1]], x = read(join(FIX, f));
      if (!t || !Array.isArray(x.graphs) || !x.series?.Z || !x.series?.F) { problems.push(`${f}: shape`); continue; }
      for (const gr of x.graphs) {
        const cg = gr.coloured_graph;
        if (!wholeGraph(cg)) { problems.push(`${f} ${gr.id}: not a coloured graph`); continue; }
        let F = 0;
        if (t.bubbles.length) for (let i = 1; i < cg.colours; i++) F += faces(cg, 0, i);
        else for (let i = 0; i < cg.colours; i++) for (let j = i + 1; j < cg.colours; j++) F += faces(cg, i, j);
        if (F !== gr.n_power) problems.push(`${f} ${gr.id}: N^${gr.n_power}, faces ${F}`);
        if (typeof gr.weight !== 'string' || !/^-?\d+(\/\d+)?$/.test(gr.weight)) problems.push(`${f} ${gr.id}: weight ${gr.weight}`);
        if (cg.colours >= 2 && checkComplex(colouredGraphComplex(cg)).length) problems.push(`${f} ${gr.id}: its 2-complex is not whole`);
        const L = colouredGraphLayout(cg);
        if (!L.nodes.every((n) => Number.isFinite(n.x) && Number.isFinite(n.y))) problems.push(`${f} ${gr.id}: not drawn`);
      }
      if (x.series.Z.length !== Number(m[2]) || !x.series.Z.every((s) => Array.isArray(s.terms))) problems.push(`${f}: series to order ${x.series.Z.length}`);
    }
    assert('atlas: every Feynman graph of the fixtures is a coloured graph weighted by N to its faces, drawn, its 2-complex whole', expansions >= 4 && !problems.length, problems.slice(0, 6).join('; '));
    if (has('expand-matrix-1')) {
      const z = read(join(FIX, 'expand-matrix-1.json')).series.Z[0].terms;
      assert('atlas: the matrix model at first order: Z = 1 − λ⟨tr (MM†)²⟩ = 1 − 2N³ λ', z.length === 1 && z[0].power === 3 && z[0].coeff === '-2', JSON.stringify(z));
    }
    if (has('expand-vector-2')) {
      // ⟨(φ̄·φ)⁴⟩ = N(N+1)(N+2)(N+3): Z at λ² = ½ (N⁴ + 6N³ + 11N² + 6N).
      const z = read(join(FIX, 'expand-vector-2.json')).series.Z[1].terms.map((t) => `${t.power}:${t.coeff}`).join(' ');
      assert('atlas: the vector model at second order is the moment N(N+1)(N+2)(N+3)/2 of a complex Gaussian vector', z === '4:1/2 3:3 2:11/2 1:3', z);
    }
  }
  if (has('large-n-matrix')) {
    const L = read(join(FIX, 'large-n-matrix.json'));
    const { x, rho } = L.density;
    let mass = 0;
    for (let i = 1; i < x.length; i++) mass += 0.5 * (rho[i] + rho[i - 1]) * (x[i] - x[i - 1]);
    check('atlas: the eigenvalue density of the matrix model has total mass 1', mass - 1, 2e-3);
    assert('atlas: the planar matrix model has a critical coupling below zero and γ = −1/2', L.critical.g_c < 0 && Math.abs(L.critical.gamma + 0.5) < 1e-12);
  }
  if (has('condensate')) {
    const C = read(join(FIX, 'condensate.json'));
    const n = C.phi.length, f = C.friedmann;
    const lo = f.rho.indexOf(Math.min(...f.rho));
    assert('atlas: the condensate\'s curves are one per time and its effective Friedmann equation becomes classical at low density',
      [C.V, C.H, C.rho, C.classical.V].every((a) => a.length === n) && Math.abs(f.H2[lo] - f.classical_H2[lo]) < 0.02 * f.classical_H2[lo], `${f.H2[lo]} vs ${f.classical_H2[lo]}`);
    if (C.bounce) assert('atlas: at the bounce the volume is least and the expansion stops (H = V′/3V = 0)', Math.min(...C.V) >= C.bounce.V * (1 - 1e-9) && Math.min(...C.H.map(Math.abs)) < 0.05 * Math.max(...C.H.map(Math.abs)));
  }
  for (const kind of ['quadrangulation', 'melonic']) {
    if (!has(`random-geometry-${kind}`)) continue;
    const c = read(join(FIX, `random-geometry-${kind}.json`)).complex;
    const ok = checkComplex(c).length === 0 && (kind !== 'quadrangulation' || (c.faces.every((f) => f.boundary.length === 4) && eulerCharacteristic(c) === 2));
    assert(`atlas: the random ${kind} of the fixtures is whole${kind === 'quadrangulation' ? ', of squares, a sphere' : ''}`, ok, checkComplex(c).join('; '));
    if (has(`dimensions-${c.name}`)) {
      const d = read(join(FIX, `dimensions-${c.name}.json`));
      assert(`atlas: the dimensions of the random ${kind}: as many values as times and radii, all finite`,
        d.spectral.t.length === d.spectral.ds.length && d.hausdorff.r.length === d.hausdorff.dh.length && [...d.spectral.ds, ...d.hausdorff.dh].every(Number.isFinite));
    }
  }
  const foams = [...files].filter((f) => f.startsWith('spin-foam-'));
  if (!foams.length) skip('atlas: spin foams of the fixtures: none yet');
  else {
    const wrong = [];
    for (const f of foams) {
      const r = read(join(FIX, f));
      if (checkComplex(r.foam).length) { wrong.push(`${f}: the foam is not whole`); continue; }
      if (r.face_spins.length !== r.foam.faces.length) wrong.push(`${f}: ${r.face_spins.length} spins, ${r.foam.faces.length} faces`);
      for (const [i, v] of r.vertices.entries()) {
        if (v.faces.length !== 6 || v.spins.some((s, k) => s !== r.face_spins[v.faces[k]])) wrong.push(`${f} v${i}: its spins are not its faces'`);
        if (v.tetrahedron?.lengths && v.tetrahedron.lengths.some((l, k) => Math.abs(l - (v.spins[k] / 2 + 0.5)) > 1e-12)) wrong.push(`${f} v${i}: lengths are not j + ½`);
      }
    }
    assert('atlas: every spin foam of the fixtures is a whole 2-complex, a spin per face, each vertex six of them, edges of length j + ½', !wrong.length, wrong.slice(0, 5).join('; '));
  }
}

// ---- the page and the server ------------------------------------------------------------------------
{
  const html = readFileSync(join(APP, 'index.html'), 'utf8');
  assert('atlas: three.js is pinned to 0.170.0 in the import map', html.includes('three@0.170.0/build/three.module.js') && html.includes('three@0.170.0/examples/jsm/'));
  const sources = readdirSync(join(APP, 'src'), { recursive: true }).filter((f) => String(f).endsWith('.js')).map(String);
  const missing = [];
  for (const f of sources) {
    const text = readFileSync(join(APP, 'src', f), 'utf8');
    for (const m of text.matchAll(/from '(\.[^']+)'/g)) {
      const target = resolve(join(APP, 'src', f, '..'), m[1]);
      if (!existsSync(target)) missing.push(`${f} -> ${m[1]}`);
    }
  }
  assert('atlas: every module the page imports exists', !missing.length && sources.length >= 12, missing.join('; '));
  const recipe = read(join(APP, 'recipe.json'));
  assert('atlas: the recipe builds the engine by engine.py', recipe.steps.length === 1 && recipe.steps[0].run[0] === 'applications/complex-atlas/engine.py');
  const py = spawnSync(process.platform === 'win32' ? 'python' : 'python3', ['-c', 'import ast,sys; [ast.parse(open(f, encoding="utf-8").read()) for f in sys.argv[1:]]', join(APP, 'serve.py'), join(APP, 'engine.py')], { encoding: 'utf8' });
  if (py.error) skip('atlas: the server\'s Python: no Python here');
  else assert('atlas: the server and the engine\'s builder are Python that parses', py.status === 0, py.stderr.trim().split('\n').pop());
}

// ---- the engine, when it is built ----------------------------------------------------------------------
{
  const python = process.platform === 'win32' ? 'python' : 'python3';
  const st = spawnSync(python, [join(APP, 'engine.py'), 'status'], { encoding: 'utf8' });
  let current = false;
  try { current = JSON.parse(st.stdout).current; } catch { current = false; }
  if (!current) skip('atlas: the catalogue built and analysed by the engine: the engine is not built (python applications/complex-atlas/engine.py build)');
  else {
    const ask = (request) => {
      const r = spawnSync(python, [join(APP, 'engine.py'), 'run', '-'], { input: JSON.stringify(request), encoding: 'utf8', maxBuffer: 1 << 28 });
      try { return JSON.parse(r.stdout); } catch { return { error: (r.stderr || r.stdout || '').trim().split('\n').pop() }; }
    };
    const cat = ask({ op: 'catalogue' });
    assert('atlas: the engine gives its catalogue, every family with a name, a title and a group',
      Array.isArray(cat.families) && cat.families.length > 0 && cat.families.every((f) => f.id && f.title && f.group), cat.error || '');
    const wrong = [];
    for (const f of cat.families || []) {
      const params = Object.fromEntries((f.params || []).map((p) => [p.name, p.default]));
      const b = ask({ op: 'build', family: f.id, params });
      if (!b.complex) { wrong.push(`${f.id}: ${b.error}`); continue; }
      const p = checkComplex(b.complex);
      if (p.length) { wrong.push(`${f.id}: ${p[0]}`); continue; }
      const a = ask({ op: 'analyse', complex: b.complex }).analysis;
      if (!a) { wrong.push(`${f.id}: no analysis`); continue; }
      const ranks = a.homology.reduce((s, h) => s + (h.dimension % 2 ? -h.rank : h.rank), 0);
      if (a.euler !== eulerCharacteristic(b.complex) || ranks !== a.euler) wrong.push(`${f.id}: χ ${a.euler}, Σ(−1)^k rank ${ranks}`);
    }
    assert('atlas: every family built by the engine is a whole complex, and its analysis obeys Euler-Poincaré', !wrong.length, wrong.join('; '));
    // The generation of spacetime, where the engine has it.
    const th = ask({ op: 'theories' });
    if (!th.theories) skip(`atlas: the engine's theories: ${th.error || 'none'}`);
    else {
      assert('atlas: the engine gives its theories, each with an action, a family and bubbles', th.theories.length >= 4 && th.theories.every((t) => t.id && t.action && t.family && Array.isArray(t.bubbles)));
      const x = ask({ op: 'expand', theory: 'matrix', params: {}, order: 2 });
      if (!x.graphs) skip(`atlas: the engine's expansion of the matrix model: ${x.error || 'none'}`);
      else {
        const { colouredGraphComplex } = await import('../../../visualization/js/index.js');
        const bad = x.graphs.filter((gr) => checkComplex(colouredGraphComplex(gr.coloured_graph)).length);
        assert('atlas: the engine expands the matrix model to order 2 into coloured graphs whose 2-complexes are whole', x.graphs.length > 0 && !bad.length, `${bad.length} not whole`);
      }
      const b = ask({ op: 'build', family: 'boundary-simplex', params: { n: 4 } });
      const sf = b.complex ? ask({ op: 'spin_foam', complex: b.complex, r: 5, mode: 'uniform', two_j: 2, sweeps: 10, seed: 1 }) : { error: 'no ∂Δ⁴' };
      if (!sf.foam) skip(`atlas: the engine's spin foam of ∂Δ⁴: ${sf.error || 'none'}`);
      else assert('atlas: the engine lays a spin foam on ∂Δ⁴: 5 vertices of 6 faces each, the spins of a vertex those of its faces',
        sf.vertices.length === 5 && sf.vertices.every((v) => v.faces.length === 6 && v.spins.every((s, k) => s === sf.face_spins[v.faces[k]])) && checkComplex(sf.foam).length === 0);
    }
  }
}

done();
