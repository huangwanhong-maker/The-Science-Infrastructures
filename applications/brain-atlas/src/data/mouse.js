// The mouse brain: the Allen Institute's Common Coordinate Framework (CCFv3), its average
// template and its 2017 annotation at 50 µm, with the ontology that names every label,
// read by the neuroscience domain from the files of the data library's cache and brought
// to the RAS orientation the sections and the viewer use.

import {
  readNrrd, fromAllenStructureGraph, ccfVoxelAffine, reorient, robustRange, toUint8,
} from '../../../../domains/neuroscience/js/index.js';
import { fetchBytes } from './anatomy.js';

const CACHE = 'data';
export const MOUSE_FILES = {
  template: 'allen-mouse-ccf/average_template_50.nrrd',
  annotation: 'allen-mouse-ccf/annotation_50.nrrd',
  ontology: 'allen-mouse-ccf/structure_graph_1.json',
};

// A NRRD file whose data are gzipped: the data inflated here (the browser inflates
// asynchronously, the domain's reader synchronously), then read by the domain.
async function readNrrdInflated(bytes, inflate) {
  let end = -1;
  for (let i = 0; i < bytes.length - 1; i++) if (bytes[i] === 0x0a && (bytes[i + 1] === 0x0a || (bytes[i + 1] === 0x0d && bytes[i + 2] === 0x0a))) { end = i; break; }
  const dataOffset = bytes[end + 1] === 0x0d ? end + 3 : end + 2;
  const header = new TextDecoder('latin1').decode(bytes.subarray(0, dataOffset));
  if (!/encoding:\s*(gzip|gz)/i.test(header)) return readNrrd(bytes);
  const raw = await inflate(bytes.subarray(dataOffset));
  const head = new TextEncoder().encode(header.replace(/encoding:\s*(gzip|gz)/i, 'encoding: raw'));
  const joined = new Uint8Array(head.length + raw.length);
  joined.set(head, 0); joined.set(raw, head.length);
  return readNrrd(joined);
}

const browserInflate = async (bytes) => {
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
};

// The mouse atlas: {ontology, t1 (display bytes, RAS), labels (compact index, RAS),
// regions (index -> Region), indexOf (Region -> index), colours, voxels, stats}.
// `read` gives the raw bytes of a file of the cache; `inflate` inflates gzip.
export async function loadMouse({ onProgress = () => {}, read = (path, p) => fetchBytes(`${CACHE}/${path}`, p), inflate = browserInflate } = {}) {
  const [tBytes, aBytes, graphBytes] = await Promise.all([
    read(MOUSE_FILES.template, (f) => onProgress(0.5 * f)),
    read(MOUSE_FILES.annotation),
    read(MOUSE_FILES.ontology),
  ]);
  const ontology = fromAllenStructureGraph(JSON.parse(new TextDecoder().decode(graphBytes)));
  const tRaw = (await readNrrdInflated(tBytes, inflate)).volume;
  const aRaw = (await readNrrdInflated(aBytes, inflate)).volume;
  // The files declare their axes as left-posterior-superior; the atlas's documentation and
  // the Allen tools give posterior-inferior-right, which the domain's catalogue follows.
  const affine = ccfVoxelAffine(50);
  const t1 = reorient({ dims: tRaw.dims, data: tRaw.data, affine }, 'RAS');
  const ann = reorient({ dims: aRaw.dims, data: aRaw.data, affine }, 'RAS');
  const range = robustRange(t1.data, 0.01, 0.995);
  // A compact index of the labels present, so that the sections' palette stays small.
  const n = ann.data.length;
  const labels = new Uint16Array(n);
  const regions = [null];
  const indexOf = new Map();
  const byLabel = new Map();
  for (let i = 0; i < n; i++) {
    const id = ann.data[i];
    if (!id) continue;
    let k = byLabel.get(id);
    if (k === undefined) {
      const r = ontology.byLabel(id);
      k = r ? regions.length : 0;
      if (r) { regions.push(r); indexOf.set(r, k); }
      byLabel.set(id, k);
    }
    labels[i] = k;
  }
  const colours = regions.map((r) => (r && r.colour ? r.colour.map((c) => c / 255) : [0.8, 0.75, 0.8]));
  // Each index's voxels, centre and bounds.
  const [nx, ny, nz] = ann.dims, A = ann.affine;
  const m = regions.length;
  const voxels = new Uint32Array(m), sum = new Float64Array(3 * m);
  let v = 0;
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++, v++) {
    const r = labels[v];
    if (!r) continue;
    voxels[r]++;
    sum[3 * r] += A[0] * i + A[3]; sum[3 * r + 1] += A[5] * j + A[7]; sum[3 * r + 2] += A[10] * k + A[11];
  }
  return {
    ontology, regions, indexOf, colours, voxels, sum,
    t1: { dims: t1.dims, affine: t1.affine, data: toUint8(t1.data, range), range },
    labels: { dims: ann.dims, affine: ann.affine, data: labels },
    voxelVolume: 0.05 ** 3,
  };
}

// The indices of a region and of everything beneath it that has voxels.
export function indicesOf(mouse, region) {
  const out = [];
  for (const r of [region, ...region.descendants()]) { const k = mouse.indexOf.get(r); if (k) out.push(k); }
  return out;
}

// The voxels, centre and bounds of a region with everything beneath it.
export function measureRegion(mouse, region) {
  const ks = indicesOf(mouse, region);
  let n = 0;
  const c = [0, 0, 0];
  for (const k of ks) { n += mouse.voxels[k]; for (let a = 0; a < 3; a++) c[a] += mouse.sum[3 * k + a]; }
  return n ? { voxels: n, volume: n * mouse.voxelVolume, centroid: c.map((x) => x / n), indices: ks } : { voxels: 0, volume: 0, centroid: null, indices: ks };
}
