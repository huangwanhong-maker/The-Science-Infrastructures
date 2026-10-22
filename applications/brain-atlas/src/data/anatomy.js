// The atlas's data: the anatomy product of the data library (its catalogue, its
// packs of meshes and its joined labels), the template's image, and the volumes of
// the atlases a user may colour the cortex by. Everything is fetched from the
// data library's cache, under /data/, and read by the neuroscience domain.

import { readNifti, toUint8, robustRange } from '../../../../domains/neuroscience/js/index.js';

export const PRODUCT = 'data/products/anatomy-mni152.json';
const CACHE = 'data';

export const FILES = {
  t1: 'templateflow-mni152/tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz',
  // The whole head, for the scalp: the template's T1 image before the brain was extracted.
  // (TemplateFlow's desc-head_mask is on another, larger grid, and is not used.)
  head: 'templateflow-mni152/tpl-MNI152NLin2009cAsym_res-01_T1w.nii.gz',
  julich: ['ebrains-julich-brain/JulichBrainAtlas_3.1_207areas_MPM_lh_MNI152.nii.gz', 'ebrains-julich-brain/JulichBrainAtlas_3.1_207areas_MPM_rh_MNI152.nii.gz'],
  julichHierarchy: 'ebrains-julich-brain/julichbrain_v3_1.json',
  julichMap: 'ebrains-julich-brain/mni152-jba31-labelled.json',
  schaefer: {
    100: { volume: 'schaefer-2018/tpl-MNI152NLin2009cAsym_res-01_atlas-Schaefer2018_desc-100Parcels7Networks_dseg.nii.gz', table: 'schaefer-2018/tpl-MNI152NLin2009cAsym_atlas-Schaefer2018_desc-100Parcels7Networks_dseg.tsv' },
    400: { volume: 'schaefer-2018/tpl-MNI152NLin2009cAsym_res-01_atlas-Schaefer2018_desc-400Parcels7Networks_dseg.nii.gz', table: 'schaefer-2018/tpl-MNI152NLin2009cAsym_atlas-Schaefer2018_desc-400Parcels7Networks_dseg.tsv' },
  },
};

export const MISSING = `The atlas's data are not fetched and built yet. From the repository root, run:<br><br>
<code>python data/acquisition/recipe.py applications/brain-atlas/recipe.json</code><br><br>
then reload this page. It fetches about 110 MB and builds the anatomy in a minute.`;

const isGzip = (b) => b.length > 2 && b[0] === 0x1f && b[1] === 0x8b;

async function inflate(bytes) {
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

// The bytes of a file, inflated when gzipped, reporting progress as a fraction.
export async function fetchBytes(url, onProgress) {
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
  let o = 0;
  for (const c of chunks) { bytes.set(c, o); o += c.length; }
  return isGzip(bytes) ? inflate(bytes) : bytes;
}

export const fetchJson = async (url) => JSON.parse(new TextDecoder().decode(await fetchBytes(url)));
export const fetchText = async (url) => new TextDecoder().decode(await fetchBytes(url));

export async function loadProduct() { return fetchJson(PRODUCT); }

export async function loadPack(product, pack, onProgress) {
  return fetchBytes(`data/products/${product.packs[pack]}`, onProgress);
}

// A NIfTI volume of the cache, read by the domain.
export async function loadVolume(path, onProgress) {
  const bytes = await fetchBytes(path.startsWith('data/') ? path : `${CACHE}/${path}`, onProgress);
  return readNifti(bytes).volume;
}

// The template for the sections: its values brought to 0..255 over a robust range.
export async function loadTemplate(onProgress) {
  const t1 = await loadVolume(FILES.t1, onProgress);
  const range = robustRange(t1.data, 0.005, 0.997);
  return { dims: t1.dims, affine: t1.affine, data: toUint8(t1.data, range), range };
}

export async function loadLabels(product, onProgress) {
  const v = await loadVolume(`data/products/${product.labels}`, onProgress);
  return { dims: v.dims, affine: v.affine, data: v.data instanceof Uint16Array ? v.data : Uint16Array.from(v.data) };
}
