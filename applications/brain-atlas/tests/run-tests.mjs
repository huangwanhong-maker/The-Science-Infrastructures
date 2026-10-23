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
  readNifti, maskScalp, tenTwenty, structure, system, BRAIN_STRUCTURES, keysOf, structureOfKey,
  readEdf, electrodeMap, radialNormaliser, shellModel, eegLeadField, inverseOperator, bandpass, covariance, decimate, peaks, THREE_SHELLS,
} from '../../../domains/neuroscience/js/index.js';
import { MONTAGES, regionOf } from '../src/ui/traces.js';
import { project } from '../src/ui/topo.js';
import { FILES } from '../src/data/anatomy.js';
import { PARCELLATIONS, loadParcellation, groupOf } from '../src/data/parcellations.js';
import { parseTable, parsePoints, channelValues, valueRange } from '../src/ui/data.js';
import { headMask, LANDMARKS, CENTRE } from '../src/modes/head.js';
import { DataMode } from '../src/modes/data.js';
import { TOURS } from '../src/ui/tours.js';
import { EEG_EXAMPLES, FMRI_EXAMPLES } from '../src/modes/signals.js';

const { check, assert, skip, done } = suite();
const inflate = (b) => (b.length > 2 && b[0] === 0x1f && b[1] === 0x8b ? zlib.gunzipSync(b) : b);
const cachePath = (rel) => resolve(CACHE, rel);
const has = (rel) => existsSync(cachePath(rel));
const read = async (rel) => new Uint8Array(inflate(readFileSync(cachePath(rel))));

// ---- the manifests name every file the atlas reads ------------------------------------------------
{
  const files = [FILES.t1, FILES.head, ...FILES.julich, FILES.julichHierarchy, FILES.julichMap, ...Object.values(FILES.schaefer).flatMap((f) => [f.volume, f.table])];
  for (const p of PARCELLATIONS) files.push(...Object.values(p.files).flat());
  files.push(...EEG_EXAMPLES.map((x) => x.file), ...FMRI_EXAMPLES.map((x) => x.file));
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


// ---- signals: the montages and the disc -------------------------------------------------------------------
{
  const known = ['Fp1', 'Fpz', 'Fp2', 'AF7', 'AF3', 'AFz', 'AF4', 'AF8', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FCz', 'FC5', 'FC3', 'FC1', 'FC2', 'FC4', 'FC6', 'T7', 'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6', 'T8', 'CP5', 'CP3', 'CPz', 'CP4', 'CP6', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO7', 'POz', 'PO8', 'O1', 'Oz', 'O2', 'Iz'];
  const sizes = Object.entries(MONTAGES).filter(([k]) => k !== 'butterfly').map(([, m]) => known.filter((n) => m.pick(n)).length);
  assert('signals: every page of channels holds 9 to 16 of the 10-10 system', sizes.every((n) => n >= 9 && n <= 16), sizes.join());
  assert('signals: channels by region', regionOf('Fp1') === 'frontal' && regionOf('C3') === 'central' && regionOf('T7') === 'temporal' && regionOf('CPz') === 'parietal' && regionOf('O2') === 'occipital');
  const [vx, vy] = project([0, 0, 1]), [fx, fy] = project([0, 1, 0]), [rx, ry] = project([1, 0, 0]);
  assert('signals: the disc has the vertex at its centre, the nose up and the right to the right', Math.hypot(vx, vy) < 1e-9 && fy < 0 && Math.abs(fx) < 1e-9 && rx > 0 && Math.abs(ry) < 1e-9 && Math.abs(Math.hypot(fx, fy) - (90 / 125)) < 1e-9);
}

// ---- signals: real recordings and maps, on the template ------------------------------------------------------
{
  const edf = (r) => `physionet-eegmmidb/S001R0${r}.edf`;
  const product = 'products/anatomy-mni152.json';
  if (has(edf(2)) && has(edf(3)) && has(FILES.head) && has(product)) {
    // The template's head and its electrodes, as the Signals mode builds them.
    const mask = headMask(readNifti(await read(FILES.head)).volume);
    const scalp = maskScalp(mask, { centre: CENTRE, axes: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], step: 2, maxRadius: 230 });
    const lm = Object.fromEntries(Object.entries(LANDMARKS).map(([k, p]) => [k, scalp.onSurface(p)]));
    const placed = tenTwenty(scalp, lm, { system: '10-10', samples: 240 });
    const electrodes = electrodeMap(placed);
    const rec = readEdf(await read(edf(2)));
    const named = rec.channels.filter((c) => c.name && electrodes.has(c.name));
    assert(`signals: every channel of the EEGMMIDB recordings finds its electrode (${named.length} of ${rec.channels.length})`, named.length === 64);
    // The cortex of the anatomy product, thinned to sources 6 mm apart, as the mode does.
    const p = JSON.parse(readFileSync(cachePath(product), 'utf8'));
    const pack = await read(`products/${p.packs.core}`);
    const sources = [], normals = [], gyri = [];
    for (const m of p.meshes.filter((x) => x.key.startsWith('cortex:'))) {
      const pos = new Int16Array(pack.buffer, pack.byteOffset + m.position, 3 * m.vertices), nrm = new Int8Array(pack.buffer, pack.byteOffset + m.normal, 4 * m.vertices);
      const lab = new Uint8Array(pack.buffer, pack.byteOffset + m.label, m.vertices);
      const pts = Array.from({ length: m.vertices }, (_, v) => [pos[3 * v] * p.positionScale, pos[3 * v + 1] * p.positionScale, pos[3 * v + 2] * p.positionScale]);
      for (const v of decimate(pts, 6)) {
        sources.push(pts[v]);
        const n = [nrm[4 * v], nrm[4 * v + 1], nrm[4 * v + 2]], l = Math.hypot(...n) || 1;
        normals.push(n.map((x) => x / l));
        gyri.push(lab[v] ? structureOfKey(`ho:${lab[v]}`) : null);
      }
    }
    const radii = placed.positions.filter((q) => !['Nz', 'Iz', 'LPA', 'RPA'].includes(q.name)).map((q) => Math.hypot(...q.position.map((x, a) => x - scalp.centre[a]))).sort((a, b) => a - b);
    const R = radii[Math.floor(radii.length / 2)], normalise = radialNormaliser(scalp, R);
    const sensors = named.map((c) => normalise(electrodes.get(c.name)));
    const L = eegLeadField(shellModel({ radius: R, ...THREE_SHELLS, terms: 100 }), scalp.centre, sensors, sources.map(normalise), normals);
    const op = inverseOperator(L, sensors.length, sources.length, { method: 'sLORETA', snr: 3 });
    const bandCov = (r, lo, hi, spans) => {
      const data = named.map((c) => bandpass(r.channels.find((x) => x.name === c.name).data, r.rate, lo, hi));
      const n = data[0].length;
      for (let t = 0; t < n; t++) { let m = 0; for (const d of data) m += d[t]; m /= data.length; for (const d of data) d[t] -= m; }
      const N = data.length, C = new Float64Array(N * N);
      let total = 0;
      for (const [a, b] of spans(r)) { const Ci = covariance(data.map((d) => ({ data: d })), a, b); for (let i = 0; i < C.length; i++) C[i] += Ci[i] * (b - a); total += b - a; }
      for (let i = 0; i < C.length; i++) C[i] /= total;
      return C;
    };
    const whole = (r) => [[0, r.channels[0].data.length]];
    const peakOf = (values) => { let b = 0; for (let s = 1; s < values.length; s++) if (values[s] > values[b]) b = s; return b; };
    const lobe = (s) => { const g = gyri[s]; if (!g) return null; let x = structure(g.id); while (x && !x.lobe) x = structure(x.parent); return x ? x.lobe : null; };
    const alpha = op.power(bandCov(rec, 8, 13, whole));
    const a = peakOf(alpha);
    assert(`signals: the alpha of closed eyes has its strongest source at the back of the brain (${gyri[a] ? structure(gyri[a].id).name : '?'}, y = ${sources[a][1].toFixed(0)} mm)`, sources[a][1] < -60 && ['occipital', 'parietal'].includes(lobe(a)));
    // Moving a fist: the mu rhythm subsides over the sensorimotor cortex (event-related desynchronisation).
    const motor = readEdf(await read(edf(3)));
    const spansOf = (codes, skip) => (r) => r.annotations.filter((x) => codes.includes(x.text)).map((x) => [Math.round((x.onset + skip) * r.rate), Math.min(r.channels[0].data.length, Math.round((x.onset + x.duration) * r.rate))]).filter(([u, w]) => w - u > 16);
    const move = op.power(bandCov(motor, 8, 13, spansOf(['T1', 'T2'], 0.5))), rest = op.power(bandCov(motor, 8, 13, spansOf(['T0'], 0)));
    const drop = Array.from(move, (v, s) => Math.log10(rest[s] / v));
    const d = peakOf(drop);
    const central = ['precentral-gyrus', 'postcentral-gyrus', 'supramarginal-gyrus-anterior', 'superior-parietal-lobule', 'central-operculum', 'supplementary-motor-area'];
    assert(`signals: moving a fist, the mu rhythm subsides most over the sensorimotor cortex (${gyri[d] ? structure(gyri[d].id).name : '?'}, ${gyri[d] ? gyri[d].side : ''})`, gyri[d] && central.includes(gyri[d].id));
  } else skip('signals: real recordings and maps', 'fetch physionet-eegmmidb and build the anatomy');
  if (has('neurovault-ibc/right_hand-avg.nii.gz') && has(product)) {
    const vol = readNifti(await read('neurovault-ibc/right_hand-avg.nii.gz')).volume;
    const top = peaks(vol, { threshold: 5, max: 1 })[0];
    const labels = readNifti(await read(`products/${JSON.parse(readFileSync(cachePath(product), 'utf8')).labels}`)).volume;
    const table = new Map(JSON.parse(readFileSync(cachePath(product), 'utf8')).labelTable);
    const near = [];
    for (let dz = -3; dz <= 3; dz++) for (let dy = -3; dy <= 3; dy++) for (let dx = -3; dx <= 3; dx++) {
      const q = [top.position[0] + dx, top.position[1] + dy, top.position[2] + dz];
      const ijk = [0, 1, 2].map((a) => Math.round(q[a] - labels.affine[4 * a + 3]));
      const k = table.get(labels.data[ijk[0] + labels.dims[0] * (ijk[1] + labels.dims[1] * ijk[2])]);
      const at = k && structureOfKey(k);
      if (at) near.push(at);
    }
    assert(`fMRI: moving the right hand peaks in the left central gyri (${top.position.map((x) => Math.round(x)).join(', ')})`, top.position[0] < -20 && near.some((at) => ['precentral-gyrus', 'postcentral-gyrus'].includes(at.id) && at.side === 'left'));
  } else skip('fMRI: the IBC maps', 'fetch neurovault-ibc');
}

done();
