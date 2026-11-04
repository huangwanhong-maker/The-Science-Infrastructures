// Checks of Fabrica: what is the atlas's own. The science (the frames, the parts and their systems, the
// measures, the sections) is checked with the anatomy domain; these check that the page reads the product as
// it is written, maps LPS to its scene as it says, names what the tours name, and that the engine answers as
// the page expects. Those that need the data or the engine skip, and say so, when they are missing.
// Run from the repository root:  node applications/fabrica/tests/run-tests.mjs

import { spawn } from 'node:child_process';
import { existsSync, readFileSync, statSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite } from '../../../development/tools/checks.mjs';

const { assert, skip, done } = suite();
const APP = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const ROOT = resolve(APP, '../..');
const PRODUCTS = join(ROOT, 'data/cache/products');

const { capital, search, searchIndex } = await import(toUrl('src/tree.js'));
const { SYSTEM_STYLE } = await import(toUrl('src/palette.js'));
const { TOURS } = await import(toUrl('src/tours.js'));
const data = await import(toUrl('src/data.js'));

function toUrl(rel, base = APP) {
  const p = join(base, rel).replace(/\\/g, '/');
  return /^[A-Za-z]:/.test(p) ? `file:///${p}` : `file://${p}`;
}
function rootUrl(rel) { return toUrl(rel, ROOT); }

// ---- without the data ------------------------------------------------------------------------------------------

assert('names read as a reader expects: a capital first, a closing roman numeral in capitals',
  capital('hepatovenous segment iv') === 'Hepatovenous segment IV' && capital('left kidney') === 'Left kidney' && capital('vein of the axis') === 'Vein of the axis');
{
  data.setCentre(800);
  const p = [12, -34, 1000];
  const w = data.lpsToWorld(p);
  const back = data.worldToLps(w);
  const det = 1 * (0 * 0 - 1 * -1); // the map's matrix [[1,0,0],[0,0,1],[0,−1,0]]
  assert('LPS to the scene and back is the identity, and the map is a rotation (determinant +1)',
    back.every((x, k) => Math.abs(x - p[k]) < 1e-12) && det === 1);
  const left = data.lpsToWorld([1, 0, 0]), up = data.lpsToWorld([0, 0, 801]), front = data.lpsToWorld([0, -1, 0]);
  assert('in the scene the subject\'s left is +x, superior +y and anterior +z', left[0] === 1 && up[1] === 1 && front[2] === 1);
}
{
  const keys = ['skin', 'skeleton', 'muscles', 'fascia', 'heart', 'arteries', 'veins', 'brain', 'nerves', 'eye', 'ear', 'respiratory', 'digestive', 'urinary', 'genital', 'endocrine', 'lymphoid'];
  assert('every system of the domain has a colour and a first style', keys.every((k) => SYSTEM_STYLE[k] && /^#[0-9a-f]{6}$/.test(SYSTEM_STYLE[k].colour) && ['solid', 'glass', 'off'].includes(SYSTEM_STYLE[k].show)));
  const { hexRgb, lightness } = await import(rootUrl('visualization/styles/colour.js'));
  const { separation } = await import(rootUrl('visualization/styles/cvd.js'));
  const pair = [hexRgb(SYSTEM_STYLE.arteries.colour), hexRgb(SYSTEM_STYLE.veins.colour)];
  const apart = separation(pair);
  assert('arteries and veins stay apart in normal vision and under each deficiency (OKLab ΔE > 0.1)',
    Object.values(apart).every((d) => d > 0.1), JSON.stringify(Object.fromEntries(Object.entries(apart).map(([k, d]) => [k, Number(d.toFixed(3))]))));
  assert('and the veins are the darker', lightness(pair[1]) < lightness(pair[0]) - 0.05);
}
{
  const steps = TOURS.flatMap((t) => t.steps);
  assert('every tour has steps with a title and a text', TOURS.length >= 4 && steps.every((s) => s.title && s.text && s.text.length > 40));
}

// ---- with the product --------------------------------------------------------------------------------------------

const productPath = join(PRODUCTS, 'body-bp3d.json');
if (!existsSync(productPath)) {
  skip('the product is not built: python data/acquisition/fetch.py bodyparts3d; python data/products/body_anatomy.py');
} else {
  const product = JSON.parse(readFileSync(productPath, 'utf8'));
  assert('the product holds 2,234 elements in four packs, each pack as long as it says',
    product.elements.length === 2234 && product.packs.length === 4 && product.packs.every((p) => existsSync(join(PRODUCTS, p.raw)) && statSync(join(PRODUCTS, p.raw)).size === p.bytes && existsSync(join(PRODUCTS, p.file))));
  const packs = product.packs.map((p) => readFileSync(join(PRODUCTS, p.raw)));
  let worst = 0, inside = true;
  data.setCentre(0);
  for (const e of product.elements.filter((x, i) => i % 97 === 0)) {
    const d = data.decodeElement(e, packs[e.pack]);
    for (let v = 0; v < e.vertices; v++) {
      const lps = data.worldToLps([d.position[3 * v], d.position[3 * v + 1], d.position[3 * v + 2]]);
      for (let k = 0; k < 3; k++) if (lps[k] < e.box[k] - 1e-3 || lps[k] > e.box[k + 3] + 1e-3) inside = false;
      const n = Math.hypot(d.normal[3 * v], d.normal[3 * v + 1], d.normal[3 * v + 2]) / 127;
      worst = Math.max(worst, Math.abs(n - 1));
    }
    if (Math.max(...d.index) >= e.vertices) inside = false;
  }
  assert('decoded meshes lie in their boxes, their corners among their points', inside);
  assert('their normals are unit vectors to the quantisation of a byte', worst < 0.02, `worst ${worst.toFixed(4)}`);
}

// ---- with the engine -------------------------------------------------------------------------------------------------

const build = (() => { try { return JSON.parse(readFileSync(join(ROOT, 'data/cache/engines/fabrica/build.json'), 'utf8')); } catch (_) { return null; } })();
const kept = existsSync(join(ROOT, 'data/cache/fabrica')) ? readFileSync : null;
let catalogue = null;
if (kept) {
  const { readdirSync } = await import('node:fs');
  const f = readdirSync(join(ROOT, 'data/cache/fabrica')).find((x) => /^catalogue-.*\.json$/.test(x));
  if (f) catalogue = JSON.parse(readFileSync(join(ROOT, 'data/cache/fabrica', f), 'utf8'));
}
if (!catalogue) {
  skip('the catalogue has not been computed: start the atlas once (python applications/fabrica/serve.py)');
} else {
  const concepts = new Set(catalogue.concepts.map((c) => c.id));
  const named = TOURS.flatMap((t) => t.steps.flatMap((s) => [s.select, ...(s.isolate || [])])).filter(Boolean);
  const missing = named.filter((id) => !concepts.has(id));
  assert('every part a tour names is a concept of the catalogue', missing.length === 0, missing.join(', '));
  const levels = new Set(catalogue.levels.map((l) => l.label));
  const vertebrae = TOURS.flatMap((t) => t.steps.map((s) => s.plate?.vertebra)).filter(Boolean);
  assert('every vertebra a tour cuts at is a level of the body', vertebrae.every((v) => levels.has(v)));
  assert('every system of the catalogue has a style', catalogue.systems.every((s) => SYSTEM_STYLE[s.key]));
  const atlas = { concepts: new Map(catalogue.concepts.map((c) => [c.id, c])) };
  const found = search(searchIndex(atlas), 'kidney');
  assert('finding "kidney" puts the kidney first', found[0]?.name === 'kidney', found.slice(0, 3).map((f) => f.name).join(', '));
  const gsv = search(searchIndex(atlas), 'great saphenous');
  assert('finding "great saphenous" finds the vein of either side', gsv.some((f) => f.name === 'right great saphenous vein') && gsv.some((f) => f.name === 'left great saphenous vein'));

  if (!build || !existsSync(join(ROOT, 'data/cache/engines/fabrica/fabrica'))) {
    skip('the engine is not built: python applications/fabrica/engine.py build');
  } else {
    // One engine process, asked as the server asks it.
    const wsl = build.toolchain === 'wsl';
    const toWsl = (p) => p.replace(/\\/g, '/').replace(/^([A-Za-z]):/, (m, d) => `/mnt/${d.toLowerCase()}`);
    const program = join(ROOT, 'data/cache/engines/fabrica/fabrica');
    const args = [wsl ? toWsl(program) : program, '--serve', '--data', wsl ? toWsl(PRODUCTS) : PRODUCTS];
    const proc = spawn(wsl ? 'wsl' : args[0], wsl ? ['--exec', ...args] : args.slice(1), { stdio: ['pipe', 'pipe', 'ignore'] });
    let buffer = '';
    const waiting = [];
    proc.stdout.on('data', (b) => {
      buffer += b.toString('utf8');
      let k;
      while ((k = buffer.indexOf('\n')) >= 0) { const line = buffer.slice(0, k); buffer = buffer.slice(k + 1); waiting.shift()?.(JSON.parse(line)); }
    });
    const ask = (req) => new Promise((r) => { waiting.push(r); proc.stdin.write(`${JSON.stringify(req)}\n`); });
    const t10 = catalogue.levels.find((l) => l.label === 'T10').z;
    const sec = await ask({ op: 'section', plane: 'transverse', level: t10, convention: 'radiological' });
    assert('the engine cuts the body at T10 and names the level', sec.ok && sec.result.vertebra === 'T10' && sec.result.pieces.length > 50);
    const pieces = sec.ok ? sec.result.pieces : [];
    assert('every piece\'s loops are closed lists of (u, v) and its area positive', pieces.every((p) => p.loops.every((l) => l.length >= 6 && l.length % 2 === 0) && (p.loops.length === 0 || p.area > 0)));
    const meanU = (name) => {
      const ps = pieces.filter((p) => catalogue.elements[p.e].name === name);
      let s = 0, n = 0;
      for (const p of ps) for (const l of p.loops) for (let k = 0; k < l.length; k += 2) { s += l[k]; n++; }
      return n ? s / n : NaN;
    };
    assert('seen as in imaging, the aorta lies on the image\'s right of the inferior vena cava', meanU('Descending thoracic aorta') > meanU('Inferior vena cava'));
    const neuro = await ask({ op: 'section', plane: 'transverse', level: t10, convention: 'neurological' });
    const nU = (name) => {
      const ps = neuro.result.pieces.filter((p) => catalogue.elements[p.e].name === name);
      let s = 0, n = 0;
      for (const p of ps) for (const l of p.loops) for (let k = 0; k < l.length; k += 2) { s += l[k]; n++; }
      return s / n;
    };
    assert('seen from above, on its left', nU('Descending thoracic aorta') < nU('Inferior vena cava') && neuro.result.seenFrom === 'superior');
    const kidney = catalogue.elements.find((e) => e.name === 'Left kidney');
    const prof = await ask({ op: 'profile', element: kidney.i, axis: 'long', samples: 200 });
    assert('a profile along the kidney\'s long axis integrates to its volume within 1 % (Cavalieri)', prof.ok && Math.abs(prof.result.integral / prof.result.volume - 1) < 0.01,
      prof.ok ? `${prof.result.integral} against ${prof.result.volume}` : prof.error);
    const liver = await ask({ op: 'part', concept: 'FMA7197', relation: 'partof' });
    assert('the liver, a compound of its segments and vessels, is measured from its parts: the solids among them over a litre',
      liver.ok && liver.result.elements.length > 10 && liver.result.solidVolume > 1e6 && liver.result.solids <= liver.result.elements.length,
      liver.ok ? `${liver.result.solids} of ${liver.result.elements.length} parts solid, ${(liver.result.solidVolume / 1000).toFixed(0)} cm³` : liver.error);
    const bad = await ask({ op: 'section', plane: 'oblique', level: 0 });
    assert('a plane the domain does not name is refused with a reason', !bad.ok && /sagittal, coronal or transverse/.test(bad.error));
    proc.stdin.end();
    proc.kill();
  }
}

done();
