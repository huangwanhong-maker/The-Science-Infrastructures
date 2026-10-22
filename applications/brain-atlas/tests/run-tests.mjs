// Checks of what is the Brain Atlas's own: its reading of the fetched atlases
// through the neuroscience domain, the label index it builds on the template's
// grid, the sides and the anatomy coming out right on real data, and the
// manifests naming every file it reads. The science is checked with the domain.
// Run from the repository root:  node applications/brain-atlas/tests/run-tests.mjs
// The checks need the fetched data (see the README); they are skipped, and say
// so, when a dataset is missing.

import { existsSync, readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import zlib from 'node:zlib';

import { suite, CACHE, ROOT } from '../../../development/tools/checks.mjs';
import {
  readNifti, readNrrd, fromSchaefer, fromSiibraParcellation, fromAllenStructureGraph, labelledMap, nearlyEqual, apply, ccfVoxelAffine, volume as makeVolume,
  SPACES, ATLASES, hemisphereOf,
} from '../../../domains/neuroscience/js/index.js';
import { FILES, labelIndex, regionIndices } from '../src/data/atlases.js';
import { parseTable, valuesByRegion, parsePoints, channelValues, valueRange } from '../src/ui/data.js';

const { check, assert, skip, done } = suite();
const inflate = (b) => zlib.gunzipSync(b);
const cached = (dataset, file) => { const p = resolve(CACHE, dataset, file); return existsSync(p) ? readFileSync(p) : null; };

// ---- the manifests name every file the atlas reads ------------------------------------------------
{
  const missing = [];
  for (const [atlasId, files] of Object.entries(FILES)) {
    const specs = [files.template, ...Object.values(files.parcellations)];
    for (const spec of specs) {
      const manifest = resolve(ROOT, 'data', 'registry', `${spec.dataset}.json`);
      if (!existsSync(manifest)) { missing.push(`${atlasId}: ${spec.dataset}.json`); continue; }
      const names = new Set(JSON.parse(readFileSync(manifest, 'utf8')).files.map((f) => f.file));
      for (const [key, value] of Object.entries(spec)) {
        if (key === 'dataset' || key === 'reader') continue;
        for (const f of [].concat(value)) if (typeof f === 'string' && !names.has(f)) missing.push(`${spec.dataset}: ${f}`);
      }
    }
  }
  assert('atlas: every file of every atlas is in its dataset’s manifest', missing.length === 0, missing.join(', '));
  assert('atlas: the catalogue and the files agree on the parcellations', Object.keys(FILES.human.parcellations).every((id) => ATLASES.human.parcellations.some((p) => p.id === id)) && Object.keys(FILES.mouse.parcellations).every((id) => ATLASES.mouse.parcellations.some((p) => p.id === id)));
}

// ---- pasted data -----------------------------------------------------------------------------------
{
  assert('data: tables with a header, comments, quotes and tabs', parseTable('region,value\n# a comment\n"Area 44 (IFG) - left hemisphere",0.3\nleft precentral\t1.5').length === 2);
  const { points, bad } = parsePoints('x,y,z\n-38,-22,56,3.1,left hand\n0,-60,20\nbad,1,2');
  assert('data: points with values and names', points.length === 2 && points[0].value === 3.1 && points[0].name === 'left hand' && points[1].value === null && bad.length === 1);
  const { values, unknown } = channelValues('Fz, 2\nT3, 1\nnowhere, 5', ['Fz', 'T7', 'Cz']);
  assert('data: channel names, with the old temporal names translated', values.get('Fz') === 2 && values.get('T7') === 1 && unknown.join() === 'nowhere');
  check('data: a diverging range is symmetric', Math.abs(valueRange([-1, 3], { diverging: true })[0] + 3), 0);
}

// ---- the human template and Schaefer ----------------------------------------------------------------
{
  const t1 = cached('templateflow-mni152', FILES.human.template.t1);
  const dseg = cached('schaefer-2018', FILES.human.parcellations['schaefer-2018-100-7'].volume);
  const tsv = cached('schaefer-2018', FILES.human.parcellations['schaefer-2018-100-7'].table);
  if (t1 && dseg && tsv) {
    const template = readNifti(t1, { inflate }).volume;
    const grid = SPACES['mni152-nlin-2009c-asym'].templates['1mm'];
    assert('human: the template is on the catalogue’s 1 mm grid', template.dims.join() === grid.dims.join() && nearlyEqual(template.affine, grid.affine, 1e-4));
    const parcellation = fromSchaefer(tsv.toString('utf8'));
    const vol = readNifti(dseg, { inflate }).volume;
    assert('human: Schaefer’s map is on the template’s grid', vol.dims.join() === template.dims.join() && nearlyEqual(vol.affine, template.affine, 1e-4));
    const map = labelledMap({ parcellation, volumes: [vol] });
    const index = labelIndex(map);
    assert('human: the label index joins every label to a parcel, 100 of them', index.unknownLabels === 0 && index.regions.length === 101 && index.index.length === vol.data.length);
    // The sides: every parcel named left has its centroid at negative x; the visual network is the most posterior.
    let wrongSide = 0;
    const byNetwork = new Map();
    for (const r of parcellation.leaves()) {
      const s = map.statistics(r);
      if (!s) { wrongSide++; continue; }
      if ((r.side() === 'left') !== (s.centroid[0] < 0)) wrongSide++;
      const net = r.parent.parent.acronym;
      if (!byNetwork.has(net)) byNetwork.set(net, []);
      byNetwork.get(net).push(s);
    }
    check('human: parcels named left lie at negative x, right at positive (100 parcels)', wrongSide, 0, 'parcels');
    const meanY = (list) => list.reduce((a, s) => a + s.centroid[1] * s.voxels, 0) / list.reduce((a, s) => a + s.voxels, 0);
    const meanZ = (list) => list.reduce((a, s) => a + s.centroid[2] * s.voxels, 0) / list.reduce((a, s) => a + s.voxels, 0);
    const ys = Object.fromEntries([...byNetwork].map(([k, v]) => [k, meanY(v)]));
    const zs = Object.fromEntries([...byNetwork].map(([k, v]) => [k, meanZ(v)]));
    assert('human: the visual network is the most posterior, the limbic network the lowest, the somatomotor network high', Object.entries(ys).every(([k, y]) => k === 'Vis' || y > ys.Vis) && Object.entries(zs).every(([k, z]) => k === 'Limbic' || z > zs.Limbic) && zs.SomMot > 35,
      `y: ${JSON.stringify(Object.fromEntries(Object.entries(ys).map(([k, y]) => [k, +y.toFixed(1)])))}; z: ${JSON.stringify(Object.fromEntries(Object.entries(zs).map(([k, z]) => [k, +z.toFixed(1)])))}`);
    // Region indices of a network cover its parcels.
    const vis = parcellation.get('Vis');
    assert('human: a network’s indices are those of its parcels', regionIndices(index, vis).size === vis.leaves().length);
    // The cortex in the parcellation occupies a plausible share of the brain mask.
    const maskBytes = cached('templateflow-mni152', FILES.human.template.mask);
    if (maskBytes) {
      const mask = readNifti(maskBytes, { inflate }).volume;
      let inMask = 0, labelled = 0;
      for (let i = 0; i < vol.data.length; i++) { if (vol.data[i]) { labelled++; if (mask.data[i]) inMask++; } }
      check('human: Schaefer’s parcels lie within the brain mask, but for a few per cent at its edge', 1 - inMask / labelled, 0.05);
    }
  } else {
    skip('human: the MNI template or the Schaefer parcellation is missing (python data/acquisition/fetch.py templateflow-mni152 schaefer-2018)');
  }
}

// ---- Julich-Brain ------------------------------------------------------------------------------------
{
  const spec = FILES.human.parcellations['julich-brain-3.1'];
  const hierarchy = cached('ebrains-julich-brain', spec.hierarchy), mapJson = cached('ebrains-julich-brain', spec.map);
  const lh = cached('ebrains-julich-brain', spec.fragments[0]), rh = cached('ebrains-julich-brain', spec.fragments[1]);
  if (hierarchy && mapJson && lh && rh) {
    const parcellation = fromSiibraParcellation(JSON.parse(hierarchy.toString('utf8')), { map: JSON.parse(mapJson.toString('utf8')) });
    const volumes = [lh, rh].map((b, i) => ({ volume: readNifti(b, { inflate }).volume, fragment: i, name: parcellation.fragments[i].fragment }));
    const map = labelledMap({ parcellation, volumes });
    const index = labelIndex(map);
    assert('julich: every label of both hemispheres is an area of the configuration', index.unknownLabels === 0 && index.regions.length === 415, `${index.unknownLabels} unknown, ${index.regions.length - 1} regions`);
    let wrong = 0;
    for (let k = 1; k < index.regions.length; k++) {
      const r = index.regions[k];
      const fragment = r.labels[0].fragment;
      if ((fragment === 0) !== (r.hemisphere === 'left')) wrong++;
    }
    check('julich: the left file holds the left-hemisphere areas, the right the right', wrong, 0, 'areas');
    const stats = map.labelStatistics();
    let sides = 0;
    for (const [, s] of stats) { const r = parcellation.byKey.get(s.fragment * 1048576 + s.label); if (r && ((r.hemisphere === 'left') !== (s.centroid[0] < 0))) sides++; }
    check('julich: every area’s centroid is on its own side', sides, 0, 'areas');
    const v1 = parcellation.get('Area hOc1 (V1, 17, CalcS)');
    const s = map.statistics(v1);
    assert('julich: the primary visual cortex lies at the back, near the midline, below the AC-PC plane', s && s.centroid[1] < -70 && Math.abs(s.centroid[0]) < 15 && s.centroid[2] < 15 && s.volume > 5000, s ? s.centroid.map((x) => x.toFixed(1)).join(', ') : 'none');
    const m1 = parcellation.get('Area 4p (PreCG)');
    const sm = map.statistics(m1);
    assert('julich: the primary motor area 4p lies high and near the central sulcus', sm && sm.centroid[2] > 35 && sm.centroid[1] > -40 && sm.centroid[1] < 0, sm ? sm.centroid.map((x) => x.toFixed(1)).join(', ') : 'none');
  } else {
    skip('julich: Julich-Brain is missing (python data/acquisition/fetch.py ebrains-julich-brain)');
  }
}

// ---- the mouse ---------------------------------------------------------------------------------------
{
  const spec = FILES.mouse.parcellations['allen-ccf-2017'];
  const graph = cached('allen-mouse-ccf', spec.hierarchy), annotation = cached('allen-mouse-ccf', spec.volume), template = cached('allen-mouse-ccf', FILES.mouse.template.volume);
  if (graph && annotation && template) {
    const parcellation = fromAllenStructureGraph(JSON.parse(graph.toString('utf8')));
    const raw = readNrrd(annotation, { inflate }).volume;
    const vol = makeVolume({ dims: raw.dims, data: raw.data, affine: ccfVoxelAffine(raw.headerAffine[0]) });
    const grid = SPACES['allen-ccf-v3'].templates['50um'];
    assert('mouse: the 50 um annotation is on the catalogue’s grid, and the template too', vol.dims.join() === grid.dims.join() && readNrrd(template, { inflate }).volume.dims.join() === grid.dims.join());
    const map = labelledMap({ parcellation, volumes: [vol] });
    const index = labelIndex(map);
    assert('mouse: every label of the annotation is a structure of the ontology', index.unknownLabels === 0 && index.regions.length > 600, `${index.unknownLabels} unknown, ${index.regions.length - 1} structures`);
    // The hemispheres: the brain's extent mirrors about the midline (the third axis runs left to right).
    const [nx, ny, nz] = vol.dims;
    let same = 0, total = 0;
    for (let k = 0; k < nz / 2; k += 3) for (let j = 0; j < ny; j += 3) for (let i = 0; i < nx; i += 3) {
      const a = vol.data[i + nx * (j + ny * k)] > 0, b = vol.data[i + nx * (j + ny * (nz - 1 - k))] > 0;
      if (a || b) { total++; if (a === b) same++; }
    }
    check('mouse: the brain’s extent mirrors about the midline', 1 - same / total, 0.02);
    // The isocortex at the dorsal surface: at a few columns, the first labelled voxel from the top is in the isocortex.
    const iso = parcellation.get('Isocortex');
    let dorsal = 0, tried = 0;
    for (const [i, k] of [[110, 80], [130, 70], [150, 90], [120, 150], [140, 140]]) {
      for (let j = 0; j < ny; j++) {
        const label = vol.data[i + nx * (j + ny * k)];
        if (!label) continue;
        tried++;
        const r = parcellation.byLabel(label);
        if (r && (r === iso || r.isDescendantOf(iso))) dorsal++;
        break;
      }
    }
    assert('mouse: the first structure beneath the dorsal surface is isocortex (5 columns)', tried === 5 && dorsal === 5, `${dorsal} of ${tried}`);
    // The hippocampus beneath the cortex, posterior of the middle, in each hemisphere.
    // The framework's origin is the volume's corner, so the midline is at x = half the width (5.7 mm).
    const hpf = parcellation.get('HPF');
    const s = map.statistics(hpf);
    const mid = apply(vol.affine, [nx / 2, ny / 2, nz / 2]);
    assert('mouse: the hippocampal formation spans both sides of the midline and lies behind the middle of the brain', s && Math.abs(s.centroid[0] - mid[0]) < 0.3 && s.bounds.min[0] < mid[0] - 1.5 && s.bounds.max[0] > mid[0] + 1.5 && s.centroid[1] < mid[1], s ? `${s.centroid.map((x) => x.toFixed(2))} about ${mid.map((x) => x.toFixed(2))}` : 'none');
  } else {
    skip('mouse: the Allen atlas is missing (python data/acquisition/fetch.py allen-mouse-ccf)');
  }
}

void hemisphereOf; void valuesByRegion;
done();
