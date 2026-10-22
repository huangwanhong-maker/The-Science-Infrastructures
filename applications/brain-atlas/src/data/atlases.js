// The atlas data: which files of the data library's cache make each atlas, how
// they are fetched and read, and how a parcellation becomes a label volume on
// the template's grid, ready for the stage.
//
// The science is the neuroscience domain's: it reads the files from bytes and
// joins the parcellations to their maps. This module only fetches the bytes
// (inflating gzip with the browser's DecompressionStream) and composes the
// results for drawing.

import {
  readNifti, readNrrd, fromSchaefer, fromSiibraParcellation, fromAllenStructureGraph, labelledMap, resample,
  ccfVoxelAffine, SPACES, ATLASES, toUint8, robustRange, volume as makeVolume, nearlyEqual,
} from '../../../../domains/neuroscience/js/index.js';

const DATA = 'data';

// The files of each atlas, under the data library's cache, by the manifests
// that fetch them (data/registry/*.json).
export const FILES = {
  human: {
    template: { dataset: 'templateflow-mni152', t1: 'tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz', mask: 'tpl-MNI152NLin2009cAsym_res-01_desc-brain_mask.nii.gz', head: 'tpl-MNI152NLin2009cAsym_res-01_T1w.nii.gz' },
    parcellations: {
      'julich-brain-3.1': { dataset: 'ebrains-julich-brain', reader: 'siibra', fragments: ['JulichBrainAtlas_3.1_207areas_MPM_lh_MNI152.nii.gz', 'JulichBrainAtlas_3.1_207areas_MPM_rh_MNI152.nii.gz'], hierarchy: 'julichbrain_v3_1.json', map: 'mni152-jba31-labelled.json' },
      'schaefer-2018-100-7': { dataset: 'schaefer-2018', reader: 'schaefer', volume: 'tpl-MNI152NLin2009cAsym_res-01_atlas-Schaefer2018_desc-100Parcels7Networks_dseg.nii.gz', table: 'tpl-MNI152NLin2009cAsym_atlas-Schaefer2018_desc-100Parcels7Networks_dseg.tsv' },
      'schaefer-2018-400-7': { dataset: 'schaefer-2018', reader: 'schaefer', volume: 'tpl-MNI152NLin2009cAsym_res-01_atlas-Schaefer2018_desc-400Parcels7Networks_dseg.nii.gz', table: 'tpl-MNI152NLin2009cAsym_atlas-Schaefer2018_desc-400Parcels7Networks_dseg.tsv' },
    },
    // Landmarks on the template head: scalp points in the direction of the usual positions,
    // found by casting rays from (0, -18, 0) through the head of the 2009c template on
    // 2026-10-22. Provisional; digitised positions should replace them.
    landmarks: { nasion: [0, 88.5, -45], inion: [0, -121.5, -33], lpa: [-84.5, -21, -50], rpa: [85.5, -21, -50], vertex: [0, -18, 105.5] },
    headCentre: [0, -18, 0],
  },
  mouse: {
    template: { dataset: 'allen-mouse-ccf', volume: 'average_template_50.nrrd', fine: 'average_template_25.nrrd' },
    parcellations: {
      'allen-ccf-2017': { dataset: 'allen-mouse-ccf', reader: 'allen', volume: 'annotation_50.nrrd', fine: 'annotation_25.nrrd', hierarchy: 'structure_graph_1.json' },
    },
    landmarks: null,
  },
};

export const MISSING_DATA = `The atlas data are not fetched yet. From the repository root, run:<br><br>
<code>python data/acquisition/fetch.py templateflow-mni152 schaefer-2018 ebrains-julich-brain allen-mouse-ccf</code><br><br>
then reload this page. See the README for what each dataset is.`;

// ---- fetching ---------------------------------------------------------------------------------

async function inflateStream(bytes) {
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

const isGzip = (b) => b.length > 2 && b[0] === 0x1f && b[1] === 0x8b;

// The bytes of a file of the cache, inflated when gzipped, with progress.
export async function fetchBytes(path, onProgress) {
  const url = `${DATA}/${path}`;
  const response = await fetch(url);
  if (!response.ok) throw new Error(`${url}: ${response.status}`);
  const total = Number(response.headers.get('content-length')) || 0;
  const reader = response.body.getReader();
  const chunks = [];
  let received = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    received += value.length;
    if (onProgress && total) onProgress(received / total);
  }
  const bytes = new Uint8Array(received);
  let offset = 0;
  for (const c of chunks) { bytes.set(c, offset); offset += c.length; }
  return isGzip(bytes) ? inflateStream(bytes) : bytes;
}

const fetchText = async (path) => new TextDecoder().decode(await fetchBytes(path));
const fetchJson = async (path) => JSON.parse(await fetchText(path));

// ---- the template ------------------------------------------------------------------------------

// The template of an atlas: its intensity volume (as 8-bit values for the
// stage, with the range used), the brain mask (human), and the head (human,
// loaded on demand for the scalp).
export async function loadTemplate(atlasId, progress) {
  const files = FILES[atlasId];
  if (atlasId === 'human') {
    const t1 = readNifti(await fetchBytes(`${files.template.dataset}/${files.template.t1}`, (f) => progress(0.1 + 0.5 * f, 'Reading the MNI 152 template…'))).volume;
    const mask = readNifti(await fetchBytes(`${files.template.dataset}/${files.template.mask}`, (f) => progress(0.6 + 0.1 * f, 'Reading the brain mask…'))).volume;
    const range = robustRange(t1.data, 0.01, 0.995);
    return { atlas: ATLASES.human, space: SPACES[ATLASES.human.space], volume: t1, display: toUint8(t1.data, range), range, mask, landmarks: files.landmarks, headCentre: files.headCentre, credit: 'Template: MNI 152 ICBM 2009c nonlinear asymmetric (Fonov et al. 2011), via TemplateFlow' };
  }
  if (atlasId === 'mouse') {
    const { volume: raw } = await readNrrdAsync(await fetchBytes(`${files.template.dataset}/${files.template.volume}`, (f) => progress(0.1 + 0.5 * f, 'Reading the Allen mouse brain template…')));
    // The file's header says LPS; the atlas is PIR in micrometres from the corner (spaces.js).
    const spacing = raw.headerAffine ? raw.headerAffine[0] : 50;
    const vol = makeVolume({ dims: raw.dims, data: raw.data, affine: ccfVoxelAffine(spacing) });
    const range = robustRange(vol.data, 0.01, 0.995);
    const mask = makeVolume({ dims: vol.dims, data: new Uint8Array(vol.data.length).map((_, i) => (vol.data[i] > range[0] + 0.08 * (range[1] - range[0]) ? 1 : 0)), affine: vol.affine });
    return { atlas: ATLASES.mouse, space: SPACES[ATLASES.mouse.space], volume: vol, display: toUint8(vol.data, range), range, mask, landmarks: null, headCentre: null, credit: 'Template: Allen Mouse Brain Common Coordinate Framework v3 (Wang et al. 2020), Allen Institute' };
  }
  throw new Error(`unknown atlas ${atlasId}`);
}

// A gzip-aware inflate for the domain's readers (the NRRD reader inflates the
// data part itself; it needs a synchronous function, so the data are inflated
// here beforehand and the reader given an identity).
async function readNrrdAsync(bytes) {
  // The NRRD header is text; the data after the blank line are gzipped. Inflate them first.
  let end = -1;
  for (let i = 0; i < bytes.length - 1; i++) if (bytes[i] === 0x0a && (bytes[i + 1] === 0x0a || (bytes[i + 1] === 0x0d && bytes[i + 2] === 0x0a))) { end = i; break; }
  const dataOffset = bytes[end + 1] === 0x0d ? end + 3 : end + 2;
  const header = new TextDecoder('latin1').decode(bytes.subarray(0, dataOffset));
  if (/encoding:\s*(gzip|gz)/i.test(header)) {
    const raw = await inflateStream(bytes.subarray(dataOffset));
    const text = header.replace(/encoding:\s*(gzip|gz)/i, 'encoding: raw');
    const head = new TextEncoder().encode(text);
    const joined = new Uint8Array(head.length + raw.length);
    joined.set(head, 0); joined.set(raw, head.length);
    return readNrrd(joined);
  }
  return readNrrd(bytes);
}

// ---- parcellations ------------------------------------------------------------------------------

// A parcellation of an atlas with its labelled map, on the template's grid.
export async function loadParcellation(atlasId, parcellationId, template, progress) {
  const spec = FILES[atlasId].parcellations[parcellationId];
  const meta = ATLASES[atlasId].parcellations.find((p) => p.id === parcellationId);
  if (!spec) throw new Error(`unknown parcellation ${parcellationId}`);
  let parcellation, volumes;
  if (spec.reader === 'schaefer') {
    const table = await fetchText(`${spec.dataset}/${spec.table}`);
    parcellation = fromSchaefer(table, { name: meta.name });
    const v = readNifti(await fetchBytes(`${spec.dataset}/${spec.volume}`, (f) => progress(0.3 + 0.6 * f, `Reading ${meta.name}…`))).volume;
    volumes = [v];
  } else if (spec.reader === 'siibra') {
    const [hierarchy, map] = await Promise.all([fetchJson(`${spec.dataset}/${spec.hierarchy}`), fetchJson(`${spec.dataset}/${spec.map}`)]);
    parcellation = fromSiibraParcellation(hierarchy, { map });
    volumes = [];
    for (let i = 0; i < spec.fragments.length; i++) {
      const v = readNifti(await fetchBytes(`${spec.dataset}/${spec.fragments[i]}`, (f) => progress(0.3 + 0.3 * (i + f), `Reading Julich-Brain, ${parcellation.fragments[i].fragment}…`))).volume;
      volumes.push({ volume: v, fragment: i, name: parcellation.fragments[i].fragment });
    }
  } else if (spec.reader === 'allen') {
    const graph = await fetchJson(`${spec.dataset}/${spec.hierarchy}`);
    parcellation = fromAllenStructureGraph(graph, { name: meta.name });
    const raw = (await readNrrdAsync(await fetchBytes(`${spec.dataset}/${spec.volume}`, (f) => progress(0.3 + 0.6 * f, 'Reading the Allen annotation…')))).volume;
    const spacing = raw.headerAffine ? raw.headerAffine[0] : 50;
    volumes = [makeVolume({ dims: raw.dims, data: raw.data, affine: ccfVoxelAffine(spacing) })];
  } else {
    throw new Error(`no reader for ${parcellationId}`);
  }
  // Every fragment on the template's grid, so that one label texture serves the stage.
  const grid = template.volume;
  const fragments = volumes.map((entry, i) => {
    const v = entry.volume || entry;
    const same = v.dims.join() === grid.dims.join() && nearlyEqual(v.affine, grid.affine, 1e-4);
    return { ...(entry.volume ? entry : { fragment: i, name: null }), volume: same ? v : resample(v, { dims: grid.dims, affine: grid.affine }, { method: 'nearest', outside: 0 }) };
  });
  const map = labelledMap({ parcellation, volumes: fragments, name: meta.name, space: template.space.id });
  return { id: parcellationId, meta, parcellation, map, ...labelIndex(map) };
}

// The compact label index of a map: one uint16 per voxel of the grid, 0 for
// none, else the index of a leaf region in `regions` (1-based), joining the
// fragments; and each region's index.
export function labelIndex(map) {
  const grid = map.fragments[0].volume;
  const n = grid.dims[0] * grid.dims[1] * grid.dims[2];
  const index = new Uint16Array(n);
  const regions = [null];
  const indexOf = new Map();
  const keyIndex = new Map();
  const unknown = new Map();
  for (const { volume: v, fragment } of map.fragments) {
    const d = v.data;
    for (let i = 0; i < n; i++) {
      const label = d[i];
      if (!label || index[i]) continue;
      const key = fragment * 1048576 + label;
      let k = keyIndex.get(key);
      if (k === undefined) {
        const region = map.parcellation.byKey.get(key) || null;
        if (region) {
          k = indexOf.get(region);
          if (k === undefined) { k = regions.length; regions.push(region); indexOf.set(region, k); }
        } else {
          k = 0;
          unknown.set(key, (unknown.get(key) || 0) + 1);
        }
        keyIndex.set(key, k);
      }
      index[i] = k;
    }
  }
  return { index, regions, indexOf, unknownLabels: unknown.size };
}

// The index of the region at a voxel, up the tree: a region's own index, or
// the index of any leaf beneath it.
export function regionIndices(entry, region) {
  const out = new Set();
  const own = entry.indexOf.get(region);
  if (own !== undefined) out.add(own);
  for (const r of region.descendants()) { const k = entry.indexOf.get(r); if (k !== undefined) out.add(k); }
  return out;
}
