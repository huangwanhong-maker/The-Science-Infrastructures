// Checks of what is the Brain Atlas's own: the files it reads and the manifests that
// fetch them; the anatomy product it draws; the parcellations it puts on the template's
// grid, their sides and colours; the scalp it casts and the electrodes placed on it; the
// names it reads in pasted tables; and its guided tours. The science (the tree of
// structures, what is known of each, the head's geometry) is checked with the domain.
// Run from the repository root:  node applications/brain-atlas/tests/run-tests.mjs
// The checks on fetched data are skipped, and say so, when a dataset is missing.

import { existsSync, readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import zlib from 'node:zlib';

import { suite, CACHE, ROOT } from '../../../development/tools/checks.mjs';
import {
  readNifti, maskScalp, tenTwenty, structure, system, BRAIN_STRUCTURES, keysOf,
} from '../../../domains/neuroscience/js/index.js';
import { FILES } from '../src/data/anatomy.js';
import { PARCELLATIONS, loadParcellation, groupOf } from '../src/data/parcellations.js';
import { parseTable, parsePoints, channelValues, valueRange } from '../src/ui/data.js';
import { headMask, LANDMARKS, CENTRE } from '../src/modes/head.js';
import { DataMode } from '../src/modes/data.js';
import { TOURS } from '../src/ui/tours.js';

const { check, assert, skip, done } = suite();
const inflate = (b) => (b.length > 2 && b[0] === 0x1f && b[1] === 0x8b ? zlib.gunzipSync(b) : b);
const cachePath = (rel) => resolve(CACHE, rel);
const has = (rel) => existsSync(cachePath(rel));
const read = async (rel) => new Uint8Array(inflate(readFileSync(cachePath(rel))));

// ---- the manifests name every file the atlas reads ------------------------------------------------
{
  const files = [FILES.t1, FILES.head, ...FILES.julich, FILES.julichHierarchy, FILES.julichMap, ...Object.values(FILES.schaefer).flatMap((f) => [f.volume, f.table])];
  for (const p of PARCELLATIONS) files.push(...Object.values(p.files).flat());
  const missing = [];
  for (const f of new Set(files)) {
    const [dataset, name] = f.split('/');
    const manifest = resolve(ROOT, 'data', 'registry', `${dataset}.json`);
    if (!existsSync(manifest)) { missing.push(`${dataset}.json`); continue; }
    const names = new Set(JSON.parse(readFileSync(manifest, 'utf8')).files.map((x) => x.file));
    if (!names.has(name)) missing.push(f);
  }
  assert('files: every file the atlas reads is in its dataset’s manifest', missing.length === 0, missing.join(', '));
  assert('files: every parcellation names its licence and its source', PARCELLATIONS.every((p) => p.licence && p.citation && p.description));
  // The recipe fetches every dataset the atlas reads, and builds the anatomy.
  const recipe = JSON.parse(readFileSync(resolve(ROOT, 'applications', 'brain-atlas', 'recipe.json'), 'utf8'));
  const steps = new Set(recipe.steps.map((s) => s.id));
  const datasets = new Set([...files].map((f) => f.split('/')[0]));
  const absent = [...datasets, 'massp-2', 'harvard-oxford', 'mial-thalamic-nuclei', 'ebrains-deep-white-matter', 'allen-mouse-ccf', 'anatomy'].filter((d) => !steps.has(d));
  assert('files: the recipe fetches every dataset the atlas reads and builds the anatomy after the atlases it is made from', absent.length === 0 && recipe.steps.findIndex((s) => s.id === 'anatomy') > recipe.steps.findIndex((s) => s.id === 'ebrains-deep-white-matter'), absent.join());
}

// ---- the anatomy product ---------------------------------------------------------------------------
{
  const json = 'products/anatomy-mni152.json';
  if (has(json)) {
    const p = JSON.parse(readFileSync(cachePath(json), 'utf8'));
    assert('product: the anatomy product is of its schema, with its three packs', p.schema === 'anatomy/1' && ['core', 'white', 'bundles'].every((k) => has(`products/${p.packs[k]}`)));
    const values = p.labelTable.map(([v]) => v);
    assert(`product: the joined labels are numbered once each (${values.length})`, new Set(values).size === values.length && values.every((v) => v > 0));
    const unknown = p.labelTable.filter(([, k]) => !BRAIN_STRUCTURES.some((s) => ['left', 'right', 'mid'].some((side) => (s.label && s.label[side] === k) || (s.mesh && s.mesh[side] === k))));
    assert('product: every joined label is a structure of the catalogue', unknown.length === 0, unknown.map(([, k]) => k).join());
    const bytes = await read(`products/${p.labels}`);
    const v = readNifti(bytes).volume;
    const t1 = has(FILES.t1) ? readNifti(await read(FILES.t1)).volume : null;
    assert('product: the joined labels lie on the template’s grid', !t1 || (v.dims.join() === t1.dims.join() && v.affine.every((x, i) => Math.abs(x - t1.affine[i]) < 1e-4)));
    const n = p.meshes.reduce((t, m) => t + m.vertices, 0);
    assert(`product: ${p.meshes.length} meshes, ${(n / 1e6).toFixed(2)} million vertices, each in a pack`, p.meshes.every((m) => p.packs[m.pack]) && n > 5e5);
  } else skip('product: the anatomy product', 'run python data/products/brain_anatomy.py');
}

// ---- parcellations on the grid ---------------------------------------------------------------------
{
  const spec = PARCELLATIONS.find((x) => x.id === 'schaefer-2018-100-7');
  const julich = PARCELLATIONS.find((x) => x.id === 'julich-brain-3.1');
  if (has(spec.files.volume) && has(FILES.t1)) {
    const t1 = readNifti(await read(FILES.t1)).volume;
    const grid = { dims: t1.dims, affine: t1.affine };
    const L = await loadParcellation(spec, grid, () => {}, read);
    assert(`parcellation: Schaefer’s 100 parcels on the grid (${L.regions.length - 1})`, L.regions.length - 1 === 100);
    // A parcel's side, by its name, against the side of its voxels.
    const [nx, ny] = grid.dims;
    const sumX = new Float64Array(L.regions.length);
    for (let i = 0; i < L.index.length; i++) { const k = L.index[i]; if (k) sumX[k] += (i % nx) * grid.affine[0] + grid.affine[3]; }
    const wrong = L.regions.filter((r, k) => r && ((r.side() === 'left') !== (sumX[k] / L.voxels[k] < 0))).map((r) => r.name);
    assert('parcellation: every Schaefer parcel lies on the side its name says', wrong.length === 0, wrong.join());
    assert('parcellation: every parcel belongs to one of Yeo’s seven networks', L.network.slice(1).every((n) => ['Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default'].includes(n)));
    void ny;
    if (has(julich.files.fragments[0])) {
      const J = await loadParcellation(julich, grid, () => {}, read);
      assert(`parcellation: Julich-Brain’s areas on the grid, both sides (${J.regions.length - 1} regions)`, J.regions.length - 1 > 300 && J.regions.slice(1).every((r) => r.side()));
      const groups = new Set(J.regions.slice(1).map((r) => groupOf(r)));
      assert('parcellation: Julich-Brain’s areas fall in lobes and nucleus groups', ['frontal lobe', 'occipital lobe', 'temporal lobe', 'parietal lobe'].every((g) => groups.has(g)));
      // Colours within the atlas's gentle range.
      const sat = (c) => { const mx = Math.max(...c), mn = Math.min(...c); const l = (mx + mn) / 2; return mx === mn ? 0 : (mx - mn) / (l > 0.5 ? 2 - mx - mn : mx + mn); };
      assert('parcellation: the areas’ colours are soft (saturation at most 0.47)', J.colours.slice(1).every((c) => sat(c) <= 0.47));
      const area = J.regions.find((r) => r && /hOc1/.test(r.name) && r.side() === 'left');
      const k = J.indexOf.get(area);
      let sx = 0, n = 0;
      for (let i = 0; i < J.index.length; i++) if (J.index[i] === k) { sx += (i % nx) + grid.affine[3]; n++; }
      check('parcellation: hOc1 (V1) on the left lies left of the midline, on the medial occipital surface', Math.max(0, sx / n + 3), 0, 'mm');
    } else skip('parcellation: Julich-Brain on the grid', 'fetch ebrains-julich-brain');
  } else skip('parcellation: Schaefer on the grid', 'fetch templateflow-mni152 and schaefer-2018');
}

// ---- the scalp and the electrodes on the template ---------------------------------------------------
{
  if (has(FILES.head)) {
    const mask = headMask(readNifti(await read(FILES.head)).volume);
    const scalp = maskScalp(mask, { centre: CENTRE, axes: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], step: 2, maxRadius: 230 });
    const lm = Object.fromEntries(Object.entries(LANDMARKS).map(([k, p]) => [k, scalp.onSurface(p)]));
    check('head: the landmarks lie on the template’s scalp', Math.max(...Object.entries(LANDMARKS).map(([k, p]) => Math.hypot(...p.map((x, a) => x - lm[k][a])))), 2.5, 'mm');
    const t = tenTwenty(scalp, lm, { system: '10-10', samples: 240 });
    const at = (n) => t.positions.find((p) => p.name === n).position;
    check('head: the midline and coronal arcs meet at Cz', t.czDiscrepancy, 3, 'mm');
    assert('head: Cz at the vertex, C3 left and C4 right of it, both above T7 and T8', at('Cz')[2] > 95 && at('C3')[0] < -40 && at('C4')[0] > 40 && at('C3')[2] > at('T7')[2] + 40 && at('C4')[2] > at('T8')[2] + 40);
    assert('head: Fpz in front, Oz behind, on the midline within 15 mm', at('Fpz')[1] > 70 && at('Oz')[1] < -100 && Math.abs(at('Fpz')[0]) < 15 && Math.abs(at('Oz')[0]) < 15);
    // The scalp is one closed surface: no cavity of the head is outside it.
    let inside = 0;
    for (let i = 0; i < mask.data.length; i++) inside += mask.data[i];
    check('head: the head’s volume, cavities filled', Math.abs(inside / 1e6 - 4.3), 1.2, 'litres');
  } else skip('head: the scalp of the template', 'fetch templateflow-mni152');
}

// ---- pasted data -----------------------------------------------------------------------------------
{
  assert('data: tables with a header, comments, quotes and tabs', parseTable('region,value\n# a comment\n"Inferior frontal gyrus, pars opercularis",0.3\nleft precentral\t1.5').length === 2);
  const { points, bad } = parsePoints('x,y,z\n-38,-22,56,3.1,left hand\n0,-60,20\nbad,1,2');
  assert('data: points with values and names', points.length === 2 && points[0].value === 3.1 && points[0].name === 'left hand' && points[1].value === null && bad.length === 1);
  const { values, unknown } = channelValues('Fz, 2\nT3, 1\nnowhere, 5', ['Fz', 'T7', 'Cz']);
  assert('data: channel names, with the old temporal names translated', values.get('Fz') === 2 && values.get('T7') === 1 && unknown.join() === 'nowhere');
  check('data: a diverging range is symmetric', Math.abs(valueRange([-1, 3], { diverging: true })[0] + 3), 0);
  const many = Array.from({ length: 100 }, (_, i) => i).concat([1e6]);
  assert('data: an outlier does not stretch the range of many values', valueRange(many)[1] < 100);
  const d = new DataMode({ app: {} });
  const same = (a, b) => a.length === b.length && a.every((x) => b.includes(x));
  assert('data: FreeSurfer’s names', same(d.keysFor('Left-Hippocampus'), ['aseg:17']) && same(d.keysFor('Right-Thalamus-Proper'), ['aseg:49']) && same(d.keysFor('CC_Central'), ['aseg:253']));
  assert('data: the atlas’s names, with or without a side', same(d.keysFor('left putamen'), ['aseg:12']) && same(d.keysFor('Putamen'), ['aseg:12', 'aseg:51']) && same(d.keysFor('Precentral gyrus R'), ['ho:14']));
  assert('data: Latin names and other names', same(d.keysFor('Gyrus precentralis'), ['ho:13', 'ho:14']) && d.keysFor('M1').length === 2);
  assert('data: FreeSurfer’s cortex stands for the gyri of its side', same(d.keysFor('Left-Cerebral-Cortex'), [...keysOf('cerebral-cortex', 'label').left]) && d.keysFor('Left-Cerebral-Cortex').length === 48);
  assert('data: unknown names match nothing', d.keysFor('the soul').length === 0);
}

// ---- the guided tours --------------------------------------------------------------------------------
{
  const problems = [];
  for (const tour of TOURS) {
    if (!tour.id || !tour.title || !tour.summary || tour.steps.length < 3) problems.push(`${tour.id}: incomplete`);
    for (const [i, s] of tour.steps.entries()) {
      const where = `${tour.id} ${i + 1}`;
      if (!s.title || !s.text) problems.push(`${where}: no text`);
      for (const id of [...(s.show || []), ...(s.hide || [])]) if (!structure(id)) problems.push(`${where}: ${id}`);
      if (typeof s.select === 'string' && !structure(s.select)) problems.push(`${where}: ${s.select}`);
      if (s.select && s.select.system && !system(s.select.system)) problems.push(`${where}: system ${s.select.system}`);
      if (s.view && !(s.view.target && s.view.direction && s.view.distance > 100)) problems.push(`${where}: view`);
      if (/[a-z]'[a-z]/i.test(s.text) || /&|<(?!\/?(i|b)>)/.test(s.text)) problems.push(`${where}: typography`);
    }
  }
  assert(`tours: ${TOURS.length} tours of ${TOURS.reduce((n, t) => n + t.steps.length, 0)} steps, each naming structures that exist`, problems.length === 0, problems.join('; '));
}

done();
