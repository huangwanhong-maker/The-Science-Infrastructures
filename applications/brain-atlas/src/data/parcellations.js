// The parcellations the Atlases mode draws on the brain: Julich-Brain's
// cytoarchitectonic areas and Schaefer's parcels of functional connectivity, each
// read by the neuroscience domain from the files the data library fetched, and
// brought to the template's grid as one compact index of regions.

import {
  readNifti, fromSchaefer, fromSiibraParcellation, labelledMap, resample, nearlyEqual, YEO_NETWORKS, LOBE_COLOURS,
} from '../../../../domains/neuroscience/js/index.js';
import { fetchBytes } from './anatomy.js';

const CACHE = 'data';

// The parcellations offered, with what the panel says of them.
export const PARCELLATIONS = [
  {
    id: 'julich-brain-3.1', name: 'Julich-Brain 3.1', short: 'Cytoarchitectonic areas',
    kind: 'cytoarchitecture',
    description: 'Areas told apart by the arrangement of their cells, mapped in ten post-mortem brains, cut in 20 µm sections and stained; each area’s probability map, made from the ten, gives way to its neighbours in a maximum-probability map.',
    licence: 'CC BY-NC-SA 4.0', citation: 'Amunts et al. 2020', reader: 'siibra',
    files: {
      fragments: ['ebrains-julich-brain/JulichBrainAtlas_3.1_207areas_MPM_lh_MNI152.nii.gz', 'ebrains-julich-brain/JulichBrainAtlas_3.1_207areas_MPM_rh_MNI152.nii.gz'],
      hierarchy: 'ebrains-julich-brain/julichbrain_v3_1.json',
      map: 'ebrains-julich-brain/mni152-jba31-labelled.json',
    },
  },
  {
    id: 'schaefer-2018-100-7', name: 'Schaefer 2018, 100 parcels', short: 'Parcels of resting-state networks',
    kind: 'function',
    description: 'Parcels of the cortex whose activity at rest rises and falls together, found in the functional MRI of 1489 people; each belongs to one of the seven networks of Yeo et al. (2011).',
    licence: 'MIT', citation: 'Schaefer et al. 2018', reader: 'schaefer',
    files: {
      volume: 'schaefer-2018/tpl-MNI152NLin2009cAsym_res-01_atlas-Schaefer2018_desc-100Parcels7Networks_dseg.nii.gz',
      table: 'schaefer-2018/tpl-MNI152NLin2009cAsym_atlas-Schaefer2018_desc-100Parcels7Networks_dseg.tsv',
    },
  },
  {
    id: 'schaefer-2018-400-7', name: 'Schaefer 2018, 400 parcels', short: 'Finer parcels of the same networks',
    kind: 'function',
    description: 'The same networks divided four times as finely: 400 parcels, 200 on each side.',
    licence: 'MIT', citation: 'Schaefer et al. 2018', reader: 'schaefer',
    files: {
      volume: 'schaefer-2018/tpl-MNI152NLin2009cAsym_res-01_atlas-Schaefer2018_desc-400Parcels7Networks_dseg.nii.gz',
      table: 'schaefer-2018/tpl-MNI152NLin2009cAsym_atlas-Schaefer2018_desc-400Parcels7Networks_dseg.tsv',
    },
  },
];

// What the seven networks of Yeo et al. (2011) are: the card of a Schaefer parcel says it.
export const NETWORKS = {
  Vis: { name: 'Visual network', text: 'The occipital cortex and its neighbours on the ventral temporal surface: primary visual cortex around the calcarine sulcus and the extrastriate areas that analyse form, colour and motion.' },
  SomMot: { name: 'Somatomotor network', text: 'The cortex around the central sulcus, the supplementary motor area, the opercula and the auditory cortex of Heschl’s gyrus: the primary areas of movement, touch and hearing, which rise and fall together at rest.' },
  DorsAttn: { name: 'Dorsal attention network', text: 'The intraparietal sulcus and superior parietal lobule, the frontal eye fields and the motion area of the temporo-occipital junction: attention directed voluntarily to places and features.' },
  SalVentAttn: { name: 'Salience and ventral attention network', text: 'The anterior insula, the anterior and mid-cingulate cortex, the temporoparietal junction and the frontal operculum: attention captured by what is unexpected or matters to the body.' },
  Limbic: { name: 'Limbic network', text: 'The orbitofrontal cortex and the temporal pole. These regions lie where the MRI signal is weakest, near the sinuses and the ear canals, so this network is the least certain of the seven.' },
  Cont: { name: 'Frontoparietal control network', text: 'The lateral prefrontal cortex, the inferior parietal lobule and the precuneus’s anterior part: holding and switching goals, the flexible control of what the other networks do.' },
  Default: { name: 'Default mode network', text: 'The medial prefrontal cortex, the posterior cingulate and precuneus, the angular gyrus and the lateral temporal cortex: active when attention turns from the world to memory, the self and others.' },
};

const isFile = async (path) => {
  try { const r = await fetch(`${CACHE}/${path}`, { method: 'HEAD' }); return r.ok; } catch { return false; }
};

// The parcellations whose files are in the cache.
export async function availableParcellations() {
  const out = [];
  for (const p of PARCELLATIONS) {
    const first = p.files.fragments ? p.files.fragments[0] : p.files.volume;
    if (await isFile(first)) out.push(p);
  }
  return out;
}

const hsl = (h, s, l) => {
  const k = (n) => (n + h / 30) % 12, a = s * Math.min(l, 1 - l);
  const f = (n) => l - a * Math.max(-1, Math.min(k(n) - 3, 9 - k(n), 1));
  return [f(0), f(8), f(4)];
};
function toHsl([r, g, b]) {
  const max = Math.max(r, g, b), min = Math.min(r, g, b), l = (max + min) / 2;
  if (max === min) return [0, 0, l];
  const d = max - min, s = l > 0.5 ? d / (2 - max - min) : d / (max + min);
  const h = max === r ? (g - b) / d + (g < b ? 6 : 0) : max === g ? (b - r) / d + 2 : (r - g) / d + 4;
  return [h * 60, s, l];
}
const hexRgb = (hex) => [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16) / 255);

// Colours brought into the atlas's gentle range: saturation at most 0.46, lightness 0.56 to 0.8.
export function soften(rgb, shift = 0) {
  const [h, s, l] = toHsl(rgb);
  return hsl((h + 360) % 360, Math.min(0.46, s * 0.8), Math.min(0.8, Math.max(0.56, l + shift)));
}

// Julich-Brain's configuration has no colours. Each area takes the hue of its lobe or
// of the nucleus group it belongs to, as the anatomy draws them, varied a little from
// area to area by the golden ratio, so that neighbours stand apart and the two sides
// of an area share one colour.
const GROUP_COLOURS = {
  'frontal lobe': LOBE_COLOURS.frontal, 'parietal lobe': LOBE_COLOURS.parietal, 'temporal lobe': LOBE_COLOURS.temporal,
  'occipital lobe': LOBE_COLOURS.occipital, 'limbic lobe': LOBE_COLOURS.limbic, insula: LOBE_COLOURS.insula,
  'basal ganglia': '#79b4cf', 'basal forebrain': '#f0d08a', amygdala: '#e27f7a', 'dorsal thalamus': '#8fc39f',
  'ventral thalamus': '#a6d4a8', metathalamus: '#9fc6ef', 'cerebellar nuclei': '#c98b5b', midbrain: '#cdb6ab',
  'STN (Subthalamus)': '#f0cf63',
};
// The group of a Julich-Brain region: its lobe, or the nucleus group it belongs to.
export function groupOf(region) {
  const a = region.ancestors().reverse();
  return a[2] ? a[2].name : a[a.length - 1] ? a[a.length - 1].name : '';
}
function julichColour(region, i) {
  const group = groupOf(region);
  const base = GROUP_COLOURS[group] || (/GapMap/.test(group) ? '#c9bcc6' : '#cdbfc9');
  const [h, s, l] = toHsl(hexRgb(base));
  const t = ((i * 0.618034) % 1) - 0.5;
  return hsl((h + 30 * t + 360) % 360, Math.min(0.44, s + 0.04 * Math.sign(t)), Math.min(0.8, Math.max(0.56, l + 0.12 * (((i * 0.381966) % 1) - 0.5))));
}

// A parcellation read and put on the grid: {spec, parcellation, index (Uint16Array on the
// grid; 0 none, else 1 + the leaf's place in `regions`), regions, colours (rgb 0..1 by
// index), network (by index, for Schaefer), voxels (count by index)}. `read` gives the
// inflated bytes of a file of the cache (in the page, by fetching it).
export async function loadParcellation(spec, grid, onProgress = () => {}, read = (path, progress) => fetchBytes(`${CACHE}/${path}`, progress)) {
  const text = async (path) => new TextDecoder().decode(await read(path));
  let parcellation, volumes;
  if (spec.reader === 'schaefer') {
    parcellation = fromSchaefer(await text(spec.files.table), { name: spec.name });
    volumes = [{ volume: readNifti(await read(spec.files.volume, onProgress)).volume, fragment: 0 }];
  } else {
    const [hierarchy, map] = await Promise.all([text(spec.files.hierarchy).then(JSON.parse), text(spec.files.map).then(JSON.parse)]);
    parcellation = fromSiibraParcellation(hierarchy, { map });
    volumes = [];
    for (let i = 0; i < spec.files.fragments.length; i++) {
      const v = readNifti(await read(spec.files.fragments[i], (f) => onProgress((i + f) / spec.files.fragments.length))).volume;
      volumes.push({ volume: v, fragment: i, name: parcellation.fragments[i]?.fragment || null });
    }
  }
  const fragments = volumes.map((e) => {
    const same = e.volume.dims.join() === grid.dims.join() && nearlyEqual(e.volume.affine, grid.affine, 1e-4);
    return { ...e, volume: same ? e.volume : resample(e.volume, { dims: grid.dims, affine: grid.affine }, { method: 'nearest', outside: 0 }) };
  });
  const map = labelledMap({ parcellation, volumes: fragments, name: spec.name, space: 'MNI152NLin2009cAsym' });
  const n = grid.dims[0] * grid.dims[1] * grid.dims[2];
  const index = new Uint16Array(n);
  const regions = [null];
  const indexOf = new Map();
  const keyIndex = new Map();
  for (const { volume: v, fragment } of map.fragments) {
    const d = v.data;
    for (let i = 0; i < n; i++) {
      const label = d[i];
      if (!label || index[i]) continue;
      const key = fragment * 1048576 + label;
      let k = keyIndex.get(key);
      if (k === undefined) {
        const region = parcellation.byKey.get(key) || null;
        if (region) {
          k = indexOf.get(region);
          if (k === undefined) { k = regions.length; regions.push(region); indexOf.set(region, k); }
        } else k = 0;
        keyIndex.set(key, k);
      }
      index[i] = k;
    }
  }
  const voxels = new Uint32Array(regions.length);
  for (let i = 0; i < n; i++) voxels[index[i]]++;
  const colours = regions.map((r) => {
    if (!r) return [0.8, 0.75, 0.78];
    for (let x = r; x; x = x.parent) if (x.colour) return soften(Array.isArray(x.colour) ? x.colour.map((c) => c / 255) : hexRgb(x.colour));
    const area = r.hemisphere && r.parent ? r.parent : r;
    return julichColour(area, area.index);
  });
  // A Schaefer parcel's network: the acronym of its top region.
  const network = regions.map((r) => (r && spec.reader === 'schaefer' ? r.root().acronym : null));
  return { spec, parcellation, index, regions, indexOf, colours, network, voxels };
}


// The colours of Yeo's seven networks as the Schaefer tables give them, softened.
export function networkColours(loaded) {
  const out = new Map();
  for (let k = 1; k < loaded.regions.length; k++) {
    const net = loaded.network[k];
    if (net && !out.has(net)) out.set(net, loaded.parcellation.get(net)?.colour || null);
  }
  return out;
}

export { YEO_NETWORKS };
