// The Signals mode: what the brain does, mapped onto it.
//
//   EEG   a recording (EDF, BDF, BrainVision or a table) read by the domain; its channels
//         placed at the electrodes of the 10-10 system on the template's scalp; a band
//         chosen; a map at the cursor, of band power, or of the contrast between two
//         conditions, drawn on the scalp (a topographic disc and the electrodes) and on the
//         cortex as an estimate of its sources (MNE, dSPM or sLORETA through a three-shell
//         spherical head); the traces below, played in real time or slowed down; a probe
//         that reads the estimate at any point of the cortex.
//   fMRI  a statistical map (fMRI, PET or any z, t or F map in MNI space) read through the
//         cortical ribbon onto the cortex and over the sections, with its peaks named by
//         the anatomy beneath them.
//   MEG   a dipole placed anywhere on the cortex, and what EEG and MEG see of it: the
//         potential on the scalp, the radial field over an idealised helmet, and the
//         estimate each would make of where it lies.
//
// The forward and inverse models are the domain's (forward.js, inverse.js). They are a
// template's: a standard head, spheres for its tissues, the sources on the outer edge of
// the template's cortex. They show what such data imply, not a subject's own anatomy.

import { loadHeadModel, ScalpSpline } from '../signals/head-model.js';
import { buildSourceSpace, spread, maxima } from '../signals/sources.js';
import { TraceView, regionColour } from '../ui/traces.js';
import { drawTopo } from '../ui/topo.js';
import { fetchBytes, loadVolume } from '../data/anatomy.js';
import { esc } from '../ui/card.js';
import {
  readEdf, readBrainVision, readSignalTable, readNifti, bandpass, welch, bandPower, peakFrequency, BANDS, globalFieldPower, covariance, intervals,
  inverseOperator, shellModel, shellPotential, sarvas, idealHelmet, THREE_SHELLS, projectToSurface, peaks, smooth as smoothVolume, standardise, invert, mniToTalairach, structureOfKey,
} from '../../../../domains/neuroscience/js/index.js';

// ---- the examples ---------------------------------------------------------------------------------

const MOTOR = { T0: { name: 'rest', rest: true }, T1: { name: 'left fist', colour: 'rgba(159, 199, 232, 0.17)' }, T2: { name: 'right fist', colour: 'rgba(232, 160, 180, 0.17)' } };
export const EEG_EXAMPLES = [
  { id: 'R02', file: 'physionet-eegmmidb/S001R02.edf', name: 'At rest, eyes closed', note: 'One minute of rest with the eyes closed: the alpha rhythm over the back of the head.', meaning: { T0: { name: 'rest', rest: true } }, band: 'alpha', view: 'power', companion: 'R01' },
  { id: 'R01', file: 'physionet-eegmmidb/S001R01.edf', name: 'At rest, eyes open', note: 'One minute of rest looking at a point: the alpha rhythm subsides.', meaning: { T0: { name: 'rest', rest: true } }, band: 'alpha', view: 'power', companion: 'R02' },
  { id: 'R03', file: 'physionet-eegmmidb/S001R03.edf', name: 'Opening and closing a fist', note: 'Two minutes: the left or the right fist, four seconds at a time, as a target appears on that side. The mu rhythm of the motor cortex, of alpha’s frequencies, subsides as the hand moves.', meaning: MOTOR, band: 'alpha', view: 'contrast', contrast: ['T1+T2', 'T0'] },
  { id: 'R04', file: 'physionet-eegmmidb/S001R04.edf', name: 'Imagining it', note: 'The same task, imagined without moving: the mu rhythm subsides too, if less.', meaning: MOTOR, band: 'alpha', view: 'contrast', contrast: ['T1+T2', 'T0'] },
];
export const FMRI_EXAMPLES = [
  { id: 'right_hand', file: 'neurovault-ibc/right_hand-avg.nii.gz', fwhm: 6, threshold: 5, name: 'Moving the right hand', note: 'Moving the right hand, against the left hand, the feet and the tongue: the hand area of the left motor and somatosensory cortex.' },
  { id: 'tongue', file: 'neurovault-ibc/tongue-avg.nii.gz', fwhm: 6, threshold: 5, look: [-54, -11, 30], name: 'Moving the tongue', note: 'Moving the tongue, against the hands and the feet: the face area, low on both precentral gyri.' },
  { id: 'story', file: 'neurovault-ibc/story-math.nii.gz', fwhm: 6, threshold: 3, look: [-58, -34, -3], name: 'Listening to stories', note: 'Listening to stories, against mental arithmetic: the language network of the temporal lobes and the frontal cortex.' },
  { id: 'faces', file: 'neurovault-ibc/face-avg.nii.gz', fwhm: 6, threshold: 4, look: [48, -56, -24], direction: [0.8, -0.48, -0.36], name: 'Seeing faces', note: 'Pictures of faces, against bodies, places and tools: the fusiform and occipital face areas, on the underside of the occipital and temporal lobes.' },
  { id: 'wm', file: 'neurovault-ibc/2back-0back.nii.gz', fwhm: 6, threshold: 3, look: [51, 8, 36], name: 'Working memory', note: 'Holding the last two items in mind (2-back), against the current one (0-back): the frontoparietal network.' },
  { id: 'tom', file: 'neurovault-ibc/mental-random.nii.gz', fwhm: 6, threshold: 3, look: [63, -46, 22], name: 'Reading minds', note: 'Watching shapes that seem to interact with intent, against random motion: the theory-of-mind network of the temporoparietal junction and the medial prefrontal cortex.' },
];
const SIGNAL_BANDS = { broad: { name: 'Broad, 1–45 Hz', lo: 1, hi: 45 }, ...BANDS };

// ---- colours ----------------------------------------------------------------------------------------

const RAMPS = {
  diverging: [[0, [46, 82, 140]], [0.25, [118, 160, 206]], [0.5, [244, 236, 232]], [0.75, [234, 150, 162]], [1, [170, 46, 92]]],
  ember: [[0, [212, 150, 172]], [0.4, [234, 140, 124]], [0.72, [244, 190, 110]], [1, [255, 244, 216]]],
  heat: [[0, [214, 112, 140]], [0.5, [244, 170, 96]], [1, [255, 243, 206]]],
  cool: [[0, [112, 116, 206]], [0.5, [104, 176, 222]], [1, [206, 242, 250]]],
};
function ramp(name, t) {
  const s = RAMPS[name];
  t = Math.max(0, Math.min(1, t));
  let k = 0;
  while (k < s.length - 2 && t > s[k + 1][0]) k++;
  const [t0, a] = s[k], [t1, b] = s[k + 1], f = (t - t0) / (t1 - t0);
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}
const TABLES = Object.fromEntries(Object.keys(RAMPS).map((k) => [k, Array.from({ length: 256 }, (_, i) => ramp(k, i / 255))]));
const css = (c) => `rgb(${c.map(Math.round).join(',')})`;
const gradient = (name) => `linear-gradient(90deg, ${Array.from({ length: 9 }, (_, i) => `${css(ramp(name, i / 8))} ${i * 12.5}%`).join(', ')})`;
const BASE = [0.86, 0.79, 0.77];
const smooth = (x) => (x <= 0 ? 0 : x >= 1 ? 1 : x * x * (3 - 2 * x));
const tidy = (x) => (!Number.isFinite(x) ? '—' : Math.abs(x) >= 100 ? x.toFixed(0) : Math.abs(x) >= 10 ? x.toFixed(1) : Math.abs(x) >= 1 ? x.toFixed(2) : x.toPrecision(2));
const signed = (x) => (x < 0 ? `−${Math.abs(Math.round(x))}` : `${Math.round(x)}`);

let lfWorker = null, lfRequests = 0;
const lfPending = new Map();
function leadField(message) {
  if (!lfWorker) {
    lfWorker = new Worker(new URL('../workers/leadfield.worker.js', import.meta.url), { type: 'module' });
    lfWorker.onmessage = (e) => { const p = lfPending.get(e.data.id); if (!p) return; lfPending.delete(e.data.id); if (e.data.error) p.reject(new Error(e.data.error)); else p.resolve(e.data.L); };
  }
  const id = ++lfRequests;
  return new Promise((resolve, reject) => { lfPending.set(id, { resolve, reject }); lfWorker.postMessage({ id, ...message }); });
}

export class SignalsMode {
  constructor(ctx) {
    this.ctx = ctx;
    this.modality = 'eeg';
    this.eeg = { band: 'alpha', view: 'power', method: 'sLORETA', snr: 3, skull: 80, threshold: 0.35, speed: 1, time: 0, playing: false, window: 2, recordings: new Map(), contrast: null, probe: null };
    // The maps of one person, unsmoothed, are smoothed as a group analysis would be (FWHM 6 mm)
    // and restandardised; each example brings its threshold: z = 5 for the strong motor maps,
    // about 3 for the weaker cognitive contrasts.
    this.fmri = { threshold: 4, sign: 'both', fwhm: 6, example: null, raw: null, volume: null };
    this.meg = { orientation: 'normal', reconstruct: 'eeg', strength: 20, source: null };
    this.scalpStyle = 'glass';
    this.names = false;
  }

  get app() { return this.ctx.app; }
  get viewer() { return this.app.viewer; }

  // ---- entering and leaving ---------------------------------------------------------------------------

  async enter() {
    this.active = true;
    const app = this.app;
    this.savedCortex = this.ctx.state.cortexStyle;
    if (this.ctx.state.cortexStyle !== 'solid') { this.ctx.state.cortexStyle = 'solid'; this.ctx.applyStyles(); }
    app.plateOff = true;
    this.ctx.schedulePlate(0);
    this.showDock(this.modality === 'eeg');
    app.cardOverride = () => this.cardHtml();
    this.ctx.renderCard();
    if (!this.head) {
      this.ctx.toast('Placing the electrodes on the template’s scalp…', 8000);
      this.head = await loadHeadModel();
    }
    if (!this.active) return;
    if (!this.space) {
      this.ctx.toast('Placing the sources on the cortex…', 6000);
      await new Promise((r) => setTimeout(r, 30));
      this.space = buildSourceSpace(app.instances.get('cerebral-cortex') || []);
      this.buffers = this.space.instances.map((inst) => ({ values: new Float32Array(inst.shade.length), rgb: new Float32Array(3 * inst.shade.length) }));
    }
    this.ensureScalp();
    app.cortexPainter = (inst) => this.paintCortexInst(inst);
    app.meshPainter = (inst) => this.paintMuted(inst);
    this.labelsEl = this.labelsEl || this.makeLabels();
    if (!this.viewHook) { this.viewHook = true; this.ctx.onViewHook(() => { if (this.active) this.placeLabels(); }); }
    const v = this.viewer;
    v.flyTo({ target: [0, -14, 10], direction: [-0.55, -0.2, 0.81], distance: v.fitDistance(110, 1.3), duration: 900 });
    await this.switchTo(this.modality, true);
  }

  leave() {
    this.active = false;
    const app = this.app;
    this.stop();
    app.plateOff = false;
    app.cortexPainter = null;
    app.meshPainter = null;
    app.cardOverride = null;
    if (app.sections) { app.sections.stat = null; for (const view of Object.values(app.views)) view.redraw(); }
    if (this.savedCortex && this.savedCortex !== this.ctx.state.cortexStyle) { this.ctx.state.cortexStyle = this.savedCortex; }
    const scalp = this.viewer.instances.get('scalp');
    if (scalp) this.viewer.setVisible(scalp, false);
    const helmet = this.viewer.instances.get('helmet');
    if (helmet) this.viewer.setVisible(helmet, false);
    this.viewer.setMarkers([]);
    this.viewer.setArrow('dipole', null);
    if (this.labelsEl) this.labelsEl.innerHTML = '';
    this.showDock(false);
    document.getElementById('legend').innerHTML = '';
    this.ctx.applyStyles();
    this.ctx.paintAll();
    this.ctx.schedulePlate();
  }

  showDock(on) {
    const stage = document.getElementById('stage');
    stage.classList.toggle('with-dock', on);
    const dock = document.getElementById('dock');
    dock.hidden = !on;
    if (on && !this.traces) {
      this.traces = new TraceView(dock, {
        onTime: (t) => { this.setTime(t); },
        onControl: (what, value) => { if (what === 'play') this.toggle(); if (what === 'speed') this.eeg.speed = value; },
      });
    }
  }

  async switchTo(modality, entering = false) {
    if (this.modality !== modality || entering) {
      this.stop();
      this.modality = modality;
      const v = this.viewer;
      v.setMarkers([]);
      v.setArrow('dipole', null);
      const helmet = v.instances.get('helmet');
      if (helmet) v.setVisible(helmet, false);
      if (this.app.sections) { this.app.sections.stat = null; }
      this.app.meshPainter = modality === 'fmri' ? null : (inst) => this.paintMuted(inst);
      this.ctx.paintAll();
      this.ctx.state.show = this.modality === 'fmri' ? this.ctx.state.show : '3d';
    }
    this.showDock(modality === 'eeg' && this.active);
    this.ctx.renderLeft();
    this.ctx.renderToolbar();
    if (modality === 'eeg') {
      if (!this.eeg.rec) await this.loadExample(EEG_EXAMPLES[0]);
      else this.refreshEeg();
    } else if (modality === 'fmri') {
      this.ensureScalp();
      if (!this.fmri.volume) await this.loadFmriExample(FMRI_EXAMPLES[0]);
      else this.refreshFmri();
    } else if (modality === 'meg') {
      this.ensureScalp();
      this.refreshMeg();
    }
    this.ctx.renderCard();
  }

  // ---- the scalp and its electrodes -------------------------------------------------------------------

  ensureScalp() {
    const v = this.viewer, h = this.head;
    if (!h) return;
    let scalp = v.instances.get('scalp');
    if (!scalp) scalp = v.addFreeMesh('scalp', { positions: h.mesh.positions, normals: h.mesh.normals, indices: h.mesh.indices }, { kind: 'skin', colour: '#efd3c4', pickable: false, render: 4 });
    if (!this.scalpDirs || this.scalpDirs.length !== scalp.shade.length * 3) {
      const p = scalp.mesh.geometry.attributes.position.array, c = h.centre;
      this.scalpDirs = new Float32Array(p.length);
      for (let i = 0; i < p.length; i += 3) { this.scalpDirs[i] = p[i] - c[0]; this.scalpDirs[i + 1] = p[i + 1] - c[1]; this.scalpDirs[i + 2] = p[i + 2] - c[2]; }
      this.scalpRGB = new Float32Array(p.length);
    }
    const show = this.modality !== 'fmri' && this.scalpStyle !== 'hidden';
    v.setVisible(scalp, show);
    if (show) v.setStyle(scalp, this.scalpStyle === 'map' ? 'solid' : 'glass');
    return scalp;
  }

  makeLabels() {
    const el = document.createElement('div');
    el.id = 'signal-labels';
    el.className = 'electrode-layer';
    document.getElementById('viewport').appendChild(el);
    return el;
  }

  // The electrodes in use, as markers coloured by their values (or plain).
  showElectrodes(values = null, colourOf = null) {
    const list = this.electrodes || [];
    this.viewer.setMarkers(list.map((e, i) => ({ position: e.position, radius: 2.2, colour: values && colourOf ? css(colourOf(values[i])) : '#f6dfb8' })));
    this.markerValues = values;
    this.placeLabels();
  }

  placeLabels() {
    const el = this.labelsEl;
    if (!el) return;
    if (!this.names || !this.electrodes || this.modality === 'fmri') { el.innerHTML = ''; return; }
    const v = this.viewer, cam = v.camera.position, c = this.head.centre;
    const out = [];
    for (const e of this.electrodes) {
      const p = e.position, w = v.root.localToWorld(cam.clone().set(...p)), cw = v.root.localToWorld(cam.clone().set(...c));
      if ((w.x - cw.x) * (cam.x - w.x) + (w.y - cw.y) * (cam.y - w.y) + (w.z - cw.z) * (cam.z - w.z) < 0) continue;
      const q = v.project(p);
      if (q) out.push(`<div class="electrode-label" style="left:${q.x.toFixed(1)}px;top:${q.y.toFixed(1)}px">${esc(e.name)}</div>`);
    }
    el.innerHTML = out.join('');
  }

  // ---- painting ----------------------------------------------------------------------------------

  // The cortex: the base colour, and the map where it is beyond the threshold.
  paintCortexInst(inst) {
    const i = this.space ? this.space.instances.indexOf(inst) : -1;
    if (i < 0) return;
    const { values, rgb } = this.buffers[i];
    const n = inst.shade.length;
    const map = this.cortexMap;          // {kind: 'sources' | 'vertices', values, scale, signed, threshold, rampPos, rampNeg}
    if (map && map.kind === 'sources') spread(this.space, i, map.values, values);
    const vertexValues = map && map.kind === 'vertices' ? map.values[i] : values;
    const thr = map ? map.threshold : 1, posT = map && TABLES[map.rampPos], negT = map && TABLES[map.rampNeg || map.rampPos];
    for (let v = 0; v < n; v++) {
      let r = BASE[0], g = BASE[1], b = BASE[2];
      if (map) {
        const x = vertexValues[v] / map.scale;
        const a = Math.abs(x);
        if (a >= thr && Number.isFinite(a)) {
          const t = Math.min(1, (a - thr) / Math.max(1e-6, 1 - thr));
          const c = (x < 0 && map.signed ? negT : posT)[Math.min(255, (t * 255) | 0)];
          const k = 0.55 + 0.45 * smooth(Math.min(1, (a - thr) / Math.max(1e-6, (1 - thr) * 0.25)));
          r = r * (1 - k) + (c[0] / 255) * k; g = g * (1 - k) + (c[1] / 255) * k; b = b * (1 - k) + (c[2] / 255) * k;
        }
      }
      rgb[3 * v] = r; rgb[3 * v + 1] = g; rgb[3 * v + 2] = b;
    }
    this.viewer.paintRGB(inst, rgb);
  }

  paintCortex() { for (const inst of this.space ? this.space.instances : []) this.paintCortexInst(inst); }

  // A structure other than the cortex, in a quiet tone of its own colour, so that the maps stand out.
  paintMuted(inst) {
    if (this.modality === 'fmri') return false;
    const hex = this.ctx.colourOf(inst.id) || '#d9c3cf';
    const c = [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16) / 255).map((x, i) => 0.35 * x + 0.65 * BASE[i] * 0.92);
    const n = inst.shade.length, rgb = new Float32Array(3 * n);
    for (let v = 0; v < n; v++) { rgb[3 * v] = c[0]; rgb[3 * v + 1] = c[1]; rgb[3 * v + 2] = c[2]; }
    this.viewer.paintRGB(inst, rgb);
    return true;
  }

  // The scalp, as a map of a spline (or plain).
  paintScalp(spline = null, colourOf = null) {
    const scalp = this.viewer.instances.get('scalp');
    if (!scalp || !scalp.visible) return;
    if (!spline || this.scalpStyle !== 'map') { this.viewer.paint(scalp); return; }
    const d = this.scalpDirs, rgb = this.scalpRGB;
    for (let i = 0; i < d.length; i += 3) {
      const c = colourOf(spline.at([d[i], d[i + 1], d[i + 2]]));
      rgb[i] = c[0] / 255; rgb[i + 1] = c[1] / 255; rgb[i + 2] = c[2] / 255;
    }
    this.viewer.paintRGB(scalp, rgb);
  }

  // ---- EEG: recordings --------------------------------------------------------------------------------

  async readExample(ex) {
    const e = this.eeg;
    let rec = e.recordings.get(ex.id);
    if (rec) return rec;
    this.ctx.toast(`Reading ${ex.name.toLowerCase()}…`, 6000);
    try { rec = readEdf(await fetchBytes(`data/${ex.file}`)); } catch (err) {
      this.ctx.toast('The example recordings are not fetched: run the Brain Atlas’s recipe.', 6000);
      console.error(err);
      return null;
    }
    rec.title = ex.name; rec.note = ex.note; rec.meaning = ex.meaning; rec.source = 'PhysioNet EEGMMIDB, subject 1';
    rec.id = ex.id; rec.example = ex;
    e.recordings.set(ex.id, rec);
    return rec;
  }

  async loadExample(ex) {
    const e = this.eeg;
    const rec = await this.readExample(ex);
    if (!rec) return;
    e.band = ex.band || e.band;
    e.view = ex.view || e.view;
    e.contrast = ex.contrast ? { a: `${ex.id}|${ex.contrast[0]}`, b: `${ex.id}|${ex.contrast[1]}` } : null;
    e.probe = null;
    this.turnToPeak = true;
    await this.useRecording(rec);
  }

  async openFiles(files) {
    const list = [...files];
    const byExt = (re) => list.find((f) => re.test(f.name));
    try {
      if (this.modality === 'fmri' || byExt(/\.nii(\.gz)?$/i)) { await this.openNifti(byExt(/\.nii(\.gz)?$/i) || list[0]); return; }
      let rec;
      const edf = byExt(/\.(edf|bdf|rec)$/i);
      if (edf) rec = readEdf(new Uint8Array(await edf.arrayBuffer()));
      else if (byExt(/\.vhdr$/i)) {
        const header = await byExt(/\.vhdr$/i).text();
        const dataName = (header.match(/DataFile=(.+)/) || [])[1]?.trim();
        const dataFile = list.find((f) => f.name === dataName) || byExt(/\.(eeg|dat)$/i);
        if (!dataFile) throw new Error('choose the .vhdr together with its .eeg data file (and the .vmrk markers)');
        const mk = byExt(/\.vmrk$/i);
        rec = readBrainVision({ header, data: new Uint8Array(await dataFile.arrayBuffer()), markers: mk ? await mk.text() : null });
      } else if (byExt(/\.(csv|tsv|txt)$/i)) {
        const text = await byExt(/\.(csv|tsv|txt)$/i).text();
        rec = readSignalTable(text, { rate: Number(prompt('The sampling rate, in hertz (none is needed if the first column is time):', '250')) || null });
      } else throw new Error('a recording is an EDF, BDF, BrainVision (.vhdr with .eeg) or a table');
      rec.title = list.map((f) => f.name).join(', ');
      rec.note = `${rec.format} recording`;
      rec.id = `file:${rec.title}`;
      rec.meaning = {};
      this.eeg.recordings.set(rec.id, rec);
      this.eeg.contrast = null;
      await this.useRecording(rec);
    } catch (err) {
      console.error(err);
      this.ctx.toast(`Not read: ${err.message}`, 7000);
    }
  }

  // A recording made ready: its channels placed, referenced to their average, the lead field
  // and the inverse operator for its electrodes.
  async useRecording(rec) {
    const e = this.eeg, h = this.head;
    const used = rec.channels.filter((c) => c.name && h.electrodes.has(c.name));
    if (used.length < 8) { this.ctx.toast(`Only ${used.length} channels are named by the 10-10 system; at least eight are needed.`, 7000); return; }
    // Duplicate names: the first is kept.
    const seen = new Set();
    const chans = used.filter((c) => (seen.has(c.name) ? false : (seen.add(c.name), true)));
    rec.used = chans;
    rec.unused = rec.channels.filter((c) => !chans.includes(c)).map((c) => c.label);
    e.rec = rec;
    e.filtered = new Map();
    e.time = Math.min(e.time, rec.duration);
    this.electrodes = chans.map((c) => ({ name: c.name, position: h.electrodes.get(c.name) }));
    this.spline = new ScalpSpline(this.electrodes.map((x) => x.position.map((v, a) => v - h.centre[a])));
    const key = `${chans.map((c) => c.name).join(',')}|${e.skull}`;
    if (e.lfKey !== key) {
      this.ctx.toast(`Computing the lead field: ${chans.length} electrodes × ${this.space.count} sources…`, 15000);
      const sources = this.space.positions.map((p) => h.normalise(p));
      const sensors = this.electrodes.map((x) => h.normalise(x.position));
      e.L = await leadField({ kind: 'eeg', model: { radius: h.radius, radii: THREE_SHELLS.radii, conductivities: [0.33, 0.33 / e.skull, 0.33], terms: 100 }, centre: h.centre, sensors, sources, normals: this.space.normals });
      e.lfKey = key;
      e.op = null;
    }
    this.ctx.renderLeft();
    this.refreshEeg();
  }

  operator() {
    const e = this.eeg;
    const key = `${e.lfKey}|${e.method}|${e.snr}`;
    if (!e.op || e.opKey !== key) {
      e.op = inverseOperator(e.L, this.electrodes.length, this.space.count, { method: e.method, snr: e.snr, depth: 0.8, reference: 'average' });
      e.opKey = key;
      e.scaleJ = null;
    }
    return e.op;
  }

  // The channels of a recording in the band, referenced to their average.
  bandData(rec = this.eeg.rec, band = this.eeg.band) {
    const e = this.eeg;
    rec.filtered = rec.filtered || new Map();
    if (!rec.filtered.has(band)) {
      const b = SIGNAL_BANDS[band];
      const chans = rec.used || rec.channels;
      const data = chans.map((c) => bandpass(c.data, rec.rate, b.lo, Math.min(b.hi, rec.rate / 2 - 1)));
      const n = data[0].length, m = data.length;
      for (let t = 0; t < n; t++) { let s = 0; for (let i = 0; i < m; i++) s += data[i][t]; s /= m; for (let i = 0; i < m; i++) data[i][t] -= s; }
      rec.filtered.set(band, data);
    }
    void e;
    return rec.filtered.get(band);
  }

  refreshEeg() {
    const e = this.eeg;
    if (!e.rec || !this.active || this.modality !== 'eeg') return;
    const data = this.bandData();
    const chans = e.rec.used.map((c, i) => ({ name: c.name, data: data[i], colour: regionColour(c.name) }));
    if (e.probe !== null && e.probe !== undefined) chans.unshift({ name: 'probe', data: this.probeSeries(data), colour: '#f6c47e' });
    const gfp = globalFieldPower(chans.filter((c) => c.name !== 'probe'));
    this.traces.setRecording({ rate: e.rec.rate, duration: e.rec.duration, channels: chans, annotations: e.rec.annotations, meaning: e.rec.meaning }, gfp, { band: SIGNAL_BANDS[e.band].name, note: e.view === 'potential' ? 'the map follows the cursor' : e.view === 'power' ? `power over ${e.window} s about the cursor` : 'a contrast of two conditions' });
    this.traces.setTime(e.time);
    this.spectrum = null;
    this.updateMap();
    this.ctx.renderCard();
  }

  // The estimate at the probe, through the band, over the recording.
  probeSeries(data) {
    const e = this.eeg, op = this.operator(), N = this.electrodes.length, s = e.probe, n = data[0].length;
    const out = new Float32Array(n);
    const row = op.K.subarray(s * N, s * N + N), sc = op.scale ? op.scale[s] : 1;
    for (let t = 0; t < n; t++) { let v = 0; for (let i = 0; i < N; i++) v += row[i] * data[i][t]; out[t] = v * sc; }
    // In the traces' units: scaled to the channels' typical size.
    let a = 0, b = 0;
    for (let t = 0; t < n; t += 5) { a = Math.max(a, Math.abs(out[t])); b = Math.max(b, Math.abs(data[0][t])); }
    const k = a > 0 ? (0.6 * b) / a : 1;
    for (let t = 0; t < n; t++) out[t] *= k;
    return out;
  }

  // ---- EEG: maps --------------------------------------------------------------------------------------

  updateMap() {
    const e = this.eeg;
    if (!e.rec) return;
    if (e.view === 'potential') this.mapPotential();
    else if (e.view === 'power') this.mapPower();
    else this.mapContrast();
  }

  mapPotential() {
    const e = this.eeg, data = this.bandData(), op = this.operator();
    const k = Math.min(data[0].length - 1, Math.round(e.time * e.rec.rate));
    const phi = data.map((d) => d[k]);
    // Scales fixed for the recording and band, from the 98th percentiles.
    if (!e.scaleJ || e.scaleBand !== e.band) {
      const vs = [], js = [];
      const step = Math.max(1, Math.floor(data[0].length / 160));
      for (let t = 0; t < data[0].length; t += step) {
        const x = data.map((d) => d[t]);
        for (const v of x) vs.push(Math.abs(v));
        const J = op.apply(x);
        let m = 0;
        for (let s = 0; s < J.length; s++) m = Math.max(m, Math.abs(J[s]));
        js.push(m);
      }
      vs.sort((a, b) => a - b); js.sort((a, b) => a - b);
      e.scaleV = vs[Math.floor(0.98 * (vs.length - 1))] || 1;
      // The typical moment: the median of the largest estimate at each instant.
      e.scaleJ = js[Math.floor(0.5 * (js.length - 1))] || 1;
      e.scaleBand = e.band;
    }
    // Signed: the current flows out of the cortex (warm) or into it (cool), and an
    // oscillation shows as the colours taking turns.
    const J = op.apply(phi);
    const mag = J.map(Math.abs);
    this.cortexMap = { kind: 'sources', values: J, scale: e.scaleJ, signed: true, threshold: e.threshold, rampPos: 'heat', rampNeg: 'cool' };
    this.scalpMap = { values: phi, range: [-e.scaleV, e.scaleV], signed: true, unit: 'µV', title: 'Potential' };
    this.paintMaps(mag);
  }

  mapPower() {
    const e = this.eeg, data = this.bandData(), op = this.operator(), rate = e.rec.rate;
    const n = data[0].length;
    const half = Math.max(0.25, e.window / 2);
    const t0 = Math.max(0, Math.round((e.time - half) * rate)), t1 = Math.min(n, Math.max(t0 + 8, Math.round((e.time + half) * rate)));
    const C = covariance(data.map((d) => ({ data: d })), t0, t1);
    const N = data.length;
    const pch = Array.from({ length: N }, (_, i) => C[i * N + i]);
    const P = op.power(C);
    const amp = P.map((v) => Math.sqrt(Math.max(0, v)));
    let m = 0;
    for (const v of amp) m = Math.max(m, v);
    const pmax = Math.max(...pch);
    this.cortexMap = { kind: 'sources', values: amp, scale: m || 1, signed: false, threshold: e.threshold, rampPos: 'ember' };
    this.scalpMap = { values: pch, range: [0, pmax || 1], signed: false, unit: 'µV²', title: `${SIGNAL_BANDS[e.band].name.split(',')[0]} power` };
    this.paintMaps(amp);
  }

  // The conditions a contrast may take: the annotations of the recording, and the other
  // recordings read.
  conditions() {
    const e = this.eeg, out = [];
    if (!e.rec) return out;
    const codes = [...new Set(e.rec.annotations.map((a) => a.text))];
    // A recording of one condition is named by its title; the conditions of a task by their meaning.
    for (const c of codes) out.push({ key: `${e.rec.id}|${c}`, name: codes.length === 1 ? e.rec.title : `${(e.rec.meaning[c] && e.rec.meaning[c].name) || c}`, rec: e.rec, code: c });
    // The two conditions of a task together: either fist.
    if (codes.includes('T1') && codes.includes('T2') && e.rec.meaning.T1) out.unshift({ key: `${e.rec.id}|T1+T2`, name: `either ${e.rec.meaning.T1.name.split(' ').pop()}`, rec: e.rec, code: ['T1', 'T2'] });
    for (const [id, r] of e.recordings) if (id !== e.rec.id) out.push({ key: `${id}|*`, name: r.title, rec: r, code: null });
    return out;
  }

  // The covariance of the band's data over a condition's intervals (the whole of another recording).
  conditionCovariance(cond) {
    const rec = cond.rec;
    if (!rec.used) {
      const h = this.head, names = this.electrodes.map((x) => x.name);
      rec.used = names.map((n) => rec.channels.find((c) => c.name === n)).filter(Boolean);
      if (rec.used.length !== names.length) throw new Error(`${rec.title} lacks some of the channels`);
      void h;
    }
    const data = this.bandData(rec);
    // A task's spans begin half a second after their cue, when the response has begun.
    const codes = cond.code === null ? null : [].concat(cond.code);
    const task = (c) => rec.meaning && rec.meaning[c] && !rec.meaning[c].rest;
    const spans = codes === null ? [{ from: 0, to: rec.duration }] : codes.flatMap((c) => intervals(rec, c).map((s) => (task(c) ? { ...s, from: s.from + 0.5 } : s)));
    const N = data.length, C = new Float64Array(N * N);
    let total = 0;
    for (const sp of spans) {
      const a = Math.round(sp.from * rec.rate), b = Math.min(data[0].length, Math.round(sp.to * rec.rate));
      if (b - a < 16) continue;
      const Ci = covariance(data.map((d) => ({ data: d })), a, b);
      for (let i = 0; i < C.length; i++) C[i] += Ci[i] * (b - a);
      total += b - a;
    }
    for (let i = 0; i < C.length; i++) C[i] /= total || 1;
    return { C, seconds: total / rec.rate, count: spans.length };
  }

  mapContrast() {
    const e = this.eeg, op = this.operator();
    const conds = this.conditions();
    if (conds.length < 2) {
      const companion = e.rec.example && e.rec.example.companion && EEG_EXAMPLES.find((x) => x.id === e.rec.example.companion);
      if (companion && !e.recordings.has(companion.id)) {
        this.readExample(companion).then((r) => { if (r) { e.contrast = { a: `${e.rec.id}|T0`, b: `${companion.id}|*` }; this.ctx.renderLeft(); this.refreshEeg(); } });
      } else { e.view = 'power'; this.ctx.renderLeft(); }
      this.mapPower();
      return;
    }
    if (!e.contrast || !conds.find((c) => c.key === e.contrast.a) || !conds.find((c) => c.key === e.contrast.b)) e.contrast = { a: conds[0].key, b: conds[1].key };
    const A = conds.find((c) => c.key === e.contrast.a), B = conds.find((c) => c.key === e.contrast.b);
    let ca, cb;
    try { ca = this.conditionCovariance(A); cb = this.conditionCovariance(B); } catch (err) { this.ctx.toast(err.message, 6000); return; }
    const N = this.electrodes.length;
    const chan = Array.from({ length: N }, (_, i) => Math.log10(ca.C[i * N + i] / cb.C[i * N + i]));
    const pa = op.power(ca.C), pb = op.power(cb.C);
    const src = Array.from(pa, (v, s) => Math.log10(Math.max(v, 1e-30) / Math.max(pb[s], 1e-30)));
    let m = 0;
    for (const v of src) m = Math.max(m, Math.abs(v));
    const mc = Math.max(...chan.map(Math.abs)) || 1;
    e.contrastInfo = { a: A.name, b: B.name, sa: ca.seconds, sb: cb.seconds, na: ca.count, nb: cb.count };
    this.cortexMap = { kind: 'sources', values: Float64Array.from(src), scale: m || 1, signed: true, threshold: e.threshold, rampPos: 'heat', rampNeg: 'cool' };
    this.scalpMap = { values: chan, range: [-mc, mc], signed: true, unit: 'log₁₀ ratio', title: `${A.name} ÷ ${B.name}` };
    this.paintMaps(Float64Array.from(src, Math.abs));
  }

  paintMaps(strength) {
    const sm = this.scalpMap;
    this.spline.fit(sm.values);
    const [lo, hi] = sm.range;
    const colourOf = (v) => (sm.signed ? TABLES.diverging[Math.max(0, Math.min(255, Math.round(((v - lo) / (hi - lo || 1)) * 255)))] : TABLES.ember[Math.max(0, Math.min(255, Math.round(((v - lo) / (hi - lo || 1)) * 255)))]);
    this.scalpColour = colourOf;
    this.paintScalp(this.spline, colourOf);
    this.showElectrodes(sm.values, colourOf);
    this.paintCortex();
    this.peaks = maxima(this.space, this.cortexMap.values, { count: 4, floor: this.cortexMap.scale * this.eeg.threshold });
    this.peakStrength = strength;
    this.drawTopoCard();
    this.legend();
    if (!this.playing) this.renderPeaks();
    if (this.turnToPeak && this.peaks.length) {
      this.turnToPeak = false;
      const p = this.space.positions[this.peaks[0]], c = this.head.centre;
      const d = p.map((x, a) => x - c[a]);
      const n = Math.hypot(...d) || 1;
      const dir = [d[0] / n, d[1] / n, d[2] / n + 0.55];
      const v = this.viewer;
      v.flyTo({ target: [0, -16, 12], direction: dir, distance: v.fitDistance(112, 1.25), duration: 1100 });
    }
  }

  drawTopoCard() {
    const canvas = document.querySelector('#card canvas.topo');
    const sm = this.scalpMap;
    if (!canvas || !sm) return;
    const [lo, hi] = sm.range;
    const span = hi - lo || 1;
    drawTopo(canvas, {
      field: (d) => this.spline.at(d),
      colour: (v) => (sm.signed ? TABLES.diverging : TABLES.ember)[Math.max(0, Math.min(255, Math.round(((v - lo) / span) * 255)))],
      electrodes: this.electrodes.map((x) => ({ name: x.name, direction: x.position.map((v, a) => v - this.head.centre[a]) })),
      contours: span / 10, size: 236, title: sm.title, unit: sm.unit,
    });
  }

  legend() {
    const el = document.getElementById('legend');
    const map = this.cortexMap;
    if (!map) { el.innerHTML = ''; return; }
    const e = this.eeg;
    const title = this.modality === 'meg' ? `${this.meg.reconstruct.toUpperCase()} estimate, ${e.method}` : e.view === 'contrast' ? 'Sources, log₁₀ of the power ratio' : e.view === 'potential' ? `Sources, ${e.method}: out, in` : `Sources, ${e.method}`;
    const bar = map.signed ? `linear-gradient(90deg, ${css(ramp('cool', 1))}, ${css(ramp('cool', 0))} 42%, #ddd0cc 50%, ${css(ramp('heat', 0))} 58%, ${css(ramp('heat', 1))})` : gradient('ember');
    el.innerHTML = `<div><div class="title"><span>${esc(title)}</span><span class="unit">beyond ${Math.round(map.threshold * 100)} %</span></div><div class="bar" style="background:${bar}"></div>
<div class="ticks"><span style="left:0%">${map.signed ? `−${tidy(map.scale)}` : '0'}</span><span style="left:100%">${tidy(map.scale)}</span></div>
<div class="note">On the template’s cortex, through three spheres; an estimate, not an image.</div></div>`;
  }

  renderPeaks() {
    const el = document.querySelector('#card .peaks');
    if (!el || !this.peaks) return;
    el.innerHTML = this.peaks.map((s, i) => {
      const at = this.space.structure[s];
      const st = at && this.ctx.structure(at.id);
      const v = this.cortexMap.values[s] / this.cortexMap.scale;
      return `<button type="button" data-peak="${s}"${st ? ` data-select="${esc(at.id)}" data-side="${esc(at.side)}"` : ''}><i style="background:${esc(st ? this.ctx.colourOf(at.id) : '#ccc')}"></i>${i + 1}. ${esc(st ? st.short || st.name : 'the cortex')}${at ? ` <small>${esc(at.side)}</small>` : ''} <small>${tidy(v * 100)} %</small></button>`;
    }).join('') || '<span class="note">Nothing beyond the threshold.</span>';
  }

  // ---- time --------------------------------------------------------------------------------------

  setTime(t) {
    const e = this.eeg;
    if (!e.rec) return;
    e.time = Math.max(0, Math.min(e.rec.duration, t));
    this.traces.setTime(e.time);
    if (e.view !== 'contrast') this.updateMap();
  }

  toggle() { if (this.playing) this.stop(); else this.play(); }

  play() {
    const e = this.eeg;
    if (!e.rec || this.playing) return;
    if (e.view === 'contrast') { e.view = 'potential'; this.ctx.renderLeft(); }
    this.playing = true;
    this.traces.setPlaying(true);
    let last = performance.now(), lastMap = 0;
    const frame = (now) => {
      if (!this.playing) return;
      const dt = Math.min(0.1, (now - last) / 1000);
      last = now;
      e.time += dt * e.speed;
      if (e.time >= e.rec.duration) e.time = 0;
      this.traces.setTime(e.time);
      const interval = e.view === 'potential' ? 50 : 220;
      if (now - lastMap > interval) { lastMap = now; this.updateMap(); }
      requestAnimationFrame(frame);
    };
    requestAnimationFrame(frame);
  }

  stop() {
    if (!this.playing) return;
    this.playing = false;
    if (this.traces) this.traces.setPlaying(false);
    this.renderPeaks();
  }

  // ---- the spectrum ---------------------------------------------------------------------------------

  spectrumOf() {
    const e = this.eeg;
    if (this.spectrum && this.spectrum.rec === e.rec) return this.spectrum;
    const chans = e.rec.used;
    let sum = null, freqs = null;
    const per = new Map();
    for (const c of chans) {
      const p = welch(c.data, e.rec.rate, { segment: 2 });
      per.set(c.name, p);
      if (!sum) { sum = new Float64Array(p.power.length); freqs = p.freqs; }
      for (let k = 0; k < sum.length; k++) sum[k] += p.power[k] / chans.length;
    }
    const mean = { freqs, power: sum };
    this.spectrum = { rec: e.rec, mean, per, alpha: peakFrequency(mean, 7, 14) };
    return this.spectrum;
  }

  spectrumSvg() {
    const sp = this.spectrumOf(), { freqs, power } = sp.mean;
    const W = 330, H = 132, x0 = 30, y0 = 8, w = W - x0 - 8, h = H - y0 - 24;
    const fmax = Math.min(45, freqs[freqs.length - 1]);
    let lo = Infinity, hi = -Infinity;
    for (let k = 1; k < freqs.length; k++) if (freqs[k] <= fmax && freqs[k] >= 1) { const v = Math.log10(power[k] || 1e-12); lo = Math.min(lo, v); hi = Math.max(hi, v); }
    const X = (f) => x0 + (f / fmax) * w, Y = (v) => y0 + h - ((Math.log10(v || 1e-12) - lo) / (hi - lo || 1)) * h;
    const bands = Object.entries(BANDS).map(([k, b]) => `<rect x="${X(b.lo).toFixed(1)}" y="${y0}" width="${(X(Math.min(b.hi, fmax)) - X(b.lo)).toFixed(1)}" height="${h}" class="band${k === this.eeg.band ? ' on' : ''}"/><text x="${((X(b.lo) + X(Math.min(b.hi, fmax))) / 2).toFixed(1)}" y="${y0 + h + 19}" class="bl">${b.name.slice(0, 1) === 'D' ? 'δ' : b.name.slice(0, 1) === 'T' ? 'θ' : b.name.slice(0, 1) === 'A' ? 'α' : b.name.slice(0, 1) === 'B' ? 'β' : 'γ'}</text>`).join('');
    let d = '';
    for (let k = 1; k < freqs.length && freqs[k] <= fmax; k++) d += `${k === 1 ? 'M' : 'L'}${X(freqs[k]).toFixed(1)},${Y(power[k]).toFixed(1)}`;
    const ticks = [0, 10, 20, 30, 40].filter((f) => f <= fmax).map((f) => `<text x="${X(f).toFixed(1)}" y="${y0 + h + 9}" class="tk">${f}</text>`).join('');
    const alpha = sp.alpha ? `<line x1="${X(sp.alpha).toFixed(1)}" x2="${X(sp.alpha).toFixed(1)}" y1="${y0}" y2="${y0 + h}" class="peak"/><text x="${(X(sp.alpha) + 4).toFixed(1)}" y="${y0 + 10}" class="pk">${sp.alpha.toFixed(1)} Hz</text>` : '';
    return `<svg class="spectrum" viewBox="0 0 ${W} ${H}" width="100%">${bands}<path d="${d}" class="psd"/>${alpha}${ticks}<text x="${x0 - 4}" y="${y0 + 8}" class="ax" text-anchor="end">log</text><text x="${W - 8}" y="${y0 + h + 9}" class="tk" text-anchor="end">Hz</text></svg>`;
  }

  // ---- fMRI -----------------------------------------------------------------------------------------

  async loadFmriExample(ex) {
    this.ctx.toast(`Reading ${ex.name.toLowerCase()}…`, 6000);
    let vol;
    try { vol = await loadVolume(ex.file); } catch (err) {
      console.error(err);
      this.ctx.toast('The example maps are not fetched: run the Brain Atlas’s recipe.', 6000);
      return;
    }
    await this.useVolume(vol, ex);
  }

  async openNifti(file) {
    let bytes = new Uint8Array(await file.arrayBuffer());
    if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
      bytes = new Uint8Array(await new Response(stream).arrayBuffer());
    }
    const vol = readNifti(bytes).volume;
    await this.useVolume(vol, { id: `file:${file.name}`, name: file.name, note: 'A statistical map read from a file. It is drawn as it is: its coordinates are taken to be those of MNI space.' });
  }

  useVolume(vol, ex) {
    const f = this.fmri;
    if (vol.dims.length > 3 && vol.dims[3] > 1) this.ctx.toast(`A 4-D volume: its first of ${vol.dims[3]} volumes is shown.`, 6000);
    f.raw = vol;
    f.example = ex;
    if (ex.threshold) f.threshold = ex.threshold;
    // A file is drawn as it is, unless its reader asks for smoothing.
    f.fwhm = ex.fwhm ?? 0;
    this.turnToPeak = true;
    return this.smoothFmri();
  }

  // The map as drawn: the volume as read, or smoothed by a Gaussian of FWHM f.fwhm and brought
  // back to unit spread by its median and scaled MAD, so that its thresholds read as z again.
  async smoothFmri() {
    const f = this.fmri;
    let vol = f.raw;
    if (f.fwhm > 0) {
      this.ctx.toast(`Smoothing the map, ${f.fwhm} mm…`, 4000);
      await new Promise((r) => setTimeout(r, 30));
      vol = standardise(smoothVolume(f.raw, f.fwhm, { edge: 'zero' })).volume;
    }
    f.volume = vol;
    f.inv = invert(vol.affine);
    // The cortex reads the map through its ribbon, from 1 mm outside to 6 mm beneath its
    // outer surface: the template's surface is an average's, smoother than a person's, and a
    // person's activity lies within a few millimetres of it.
    f.vertex = this.space.instances.map((inst) => projectToSurface(vol, inst.mesh.geometry.attributes.position.array, inst.mesh.geometry.attributes.normal.array, { depths: [-1, 0, 1, 2, 3, 4, 5, 6] }));
    // The deep structures, at their surfaces and a little within.
    f.deep = new Map();
    for (const [id, list] of this.app.instances) {
      if (id === 'cerebral-cortex' || id === 'cerebral-white-matter' || id === 'cerebellar-cortex' || id === 'brainstem') continue;
      for (const inst of list) f.deep.set(inst.key, projectToSurface(vol, inst.mesh.geometry.attributes.position.array, inst.mesh.geometry.attributes.normal.array, { depths: [0.5, 1.5] }));
    }
    let mx = 0;
    for (let i = 0; i < vol.data.length; i++) { const a = Math.abs(vol.data[i]); if (a > mx && Number.isFinite(a)) mx = a; }
    f.max = mx;
    this.refreshFmri();
  }

  fmriColour(v) {
    const f = this.fmri, thr = f.threshold, top = Math.max(thr + 1, f.max * 0.9);
    if (v >= thr && f.sign !== 'negative') return TABLES.heat[Math.min(255, Math.round(((v - thr) / (top - thr)) * 255))];
    if (v <= -thr && f.sign !== 'positive') return TABLES.cool[Math.min(255, Math.round(((-v - thr) / (top - thr)) * 255))];
    return null;
  }

  refreshFmri() {
    const f = this.fmri, app = this.app;
    if (!f.volume || this.modality !== 'fmri') return;
    this.viewer.setMarkers([]);
    this.viewer.setArrow('dipole', null);
    if (this.labelsEl) this.labelsEl.innerHTML = '';
    const top = Math.max(f.threshold + 1, f.max * 0.9);
    // The cortex: values are in z units; the map's scale is its top, the threshold the z threshold.
    const keep = (v) => (f.sign === 'positive' ? (v > 0 ? v : 0) : f.sign === 'negative' ? (v < 0 ? v : 0) : v);
    this.cortexMap = { kind: 'vertices', values: f.vertex.map((arr) => Float32Array.from(arr, (v) => (Number.isFinite(v) ? keep(v) : 0))), scale: top, signed: true, threshold: f.threshold / top, rampPos: 'heat', rampNeg: 'cool' };
    this.paintCortex();
    // The deep structures, when seen through a glass cortex.
    app.meshPainter = (inst) => {
      const vals = f.deep.get(inst.key);
      if (!vals) return false;
      const n = inst.shade.length, rgb = new Float32Array(3 * n), base = this.ctx.colourOf(inst.id) || '#d9c3cf';
      const b = [1, 3, 5].map((i) => parseInt(base.slice(i, i + 2), 16) / 255).map((x) => 0.55 * x + 0.35);
      for (let v = 0; v < n; v++) {
        const c = Number.isFinite(vals[v]) ? this.fmriColour(vals[v]) : null;
        if (c) { rgb[3 * v] = c[0] / 255; rgb[3 * v + 1] = c[1] / 255; rgb[3 * v + 2] = c[2] / 255; } else { rgb[3 * v] = b[0]; rgb[3 * v + 1] = b[1]; rgb[3 * v + 2] = b[2]; }
      }
      this.viewer.paintRGB(inst, rgb);
      return true;
    };
    for (const [id, list] of app.instances) if (id !== 'cerebral-cortex') for (const inst of list) app.meshPainter(inst);
    // The sections.
    if (app.sections) {
      app.sections.stat = { vol: f.volume, inv: f.inv, colour: (v) => this.fmriColour(v) };
      for (const view of Object.values(app.views)) view.redraw();
    }
    // The peaks, positive and negative, each named by the anatomy at it.
    const pos = f.sign !== 'negative' ? peaks(f.volume, { threshold: f.threshold, sign: 1, minDistance: 12, max: 12 }) : [];
    const neg = f.sign !== 'positive' ? peaks(f.volume, { threshold: f.threshold, sign: -1, minDistance: 12, max: 6 }) : [];
    f.peaks = [...pos, ...neg].sort((a, b) => Math.abs(b.value) - Math.abs(a.value)).slice(0, 12).map((p) => ({ ...p, where: this.anatomyAt(p.position) }));
    // A new map: the sections pass through its strongest peak, and the camera turns to it.
    if (this.turnToPeak && f.peaks.length) {
      this.turnToPeak = false;
      // The site an example names as its network's; else the strongest activation that lies
      // on the lateral surface, where it can be seen; failing that, the strongest peak.
      const posPeaks = f.peaks.filter((q) => q.value > 0);
      const lateral = posPeaks.find((q) => Math.abs(q.position[0]) >= 25 && q.value >= 0.75 * (posPeaks[0] ? posPeaks[0].value : 0));
      const p = (f.example && f.example.look) || (lateral || posPeaks[0] || f.peaks[0]).position, c = this.head.centre;
      this.ctx.setPoint(p);
      const d = p.map((x, a) => x - c[a]), n = Math.hypot(...d) || 1;
      const v = this.viewer;
      v.flyTo({ target: [0, -16, 12], direction: (f.example && f.example.direction) || [d[0] / n, d[1] / n, d[2] / n + 0.35], distance: v.fitDistance(105, 1.25), duration: 1100 });
    }
    const el = document.getElementById('legend');
    el.innerHTML = `<div><div class="title"><span>${esc(f.example ? f.example.name : 'Map')}</span><span class="unit">z</span></div>
<div class="bar" style="background:linear-gradient(90deg, ${css(ramp('cool', 1))}, ${css(ramp('cool', 0))} 38%, #3a2e3a 38%, #3a2e3a 62%, ${css(ramp('heat', 0))} 62%, ${css(ramp('heat', 1))})"></div>
<div class="ticks"><span style="left:0%">−${tidy(top)}</span><span style="left:38%">−${tidy(f.threshold)}</span><span style="left:62%">${tidy(f.threshold)}</span><span style="left:100%">${tidy(top)}</span></div>
<div class="note">Read through the cortex onto its surface; drawn where |z| ≥ ${tidy(f.threshold)}.</div></div>`;
    this.ctx.renderCard();
  }

  // The structure of the anatomy at a point, or the nearest within 5 mm.
  anatomyAt(p) {
    const s = this.app.sections;
    if (!s) return null;
    const human = s.labels;
    void human;
    const ignore = new Set(['cerebral-white-matter', 'csf', 'vessel', 'white-matter-hypointensities']);
    const at = (q) => { const key = this.app.labelValues.get(s.labelAt(q)); const st = key && structureOfKey(key); return st && !ignore.has(st.id) ? st : null; };
    const here = at(p);
    if (here) return { ...here, distance: 0 };
    let best = null;
    for (let r = 1; r <= 5 && !best; r++) for (let dx = -r; dx <= r; dx++) for (let dy = -r; dy <= r; dy++) for (let dz = -r; dz <= r; dz++) {
      if (Math.max(Math.abs(dx), Math.abs(dy), Math.abs(dz)) !== r) continue;
      const st = at([p[0] + dx, p[1] + dy, p[2] + dz]);
      if (st) { const d = Math.hypot(dx, dy, dz); if (!best || d < best.distance) best = { ...st, distance: d }; }
    }
    return best;
  }

  // ---- MEG and EEG of one dipole ----------------------------------------------------------------------

  refreshMeg() {
    const m = this.meg, h = this.head, v = this.viewer;
    if (this.modality !== 'meg' || !h) return;
    if (!this.helmet) this.buildHelmet();
    const helmet = v.instances.get('helmet');
    if (helmet) { v.setVisible(helmet, true); v.setStyle(helmet, 'solid'); }
    if (!m.source) {
      // A first source: the left hand area, normal to the cortex there.
      const target = [-38, -24, 58];
      let best = 0, bd = Infinity;
      this.space.positions.forEach((p, s) => { const d = Math.hypot(...p.map((x, a) => x - target[a])); if (d < bd) { bd = d; best = s; } });
      m.source = { s: best, position: this.space.positions[best], normal: this.space.normals[best] };
    }
    const q = this.dipoleOrientation();
    const p = m.source.position;
    v.setArrow('dipole', p, q, { length: 30 });
    const qAm = q.map((x) => x * m.strength * 1e-9);
    // EEG: the potential at every point of the scalp, through three spheres after the radial normalisation.
    const model = this.megModel = this.megModel || shellModel({ radius: h.radius, ...THREE_SHELLS, terms: 100 });
    const c = h.centre;
    const toM = (x) => x.map((u, a) => (u - c[a]) / 1000);
    const r0n = toM(h.normalise(p));
    const inner = 0.97 * model.radii[0] * model.radius / 1000, rr = Math.hypot(...r0n);
    const r0 = rr > inner ? r0n.map((x) => (x * inner) / rr) : r0n;
    const eegAt = (dir) => { const n = Math.hypot(...dir); return shellPotential(model, dir.map((x) => (x / n) * model.radius / 1000), r0, qAm) * 1e6; };
    const megAt = (pos) => { const rS = toM(pos), b = sarvas(rS, toM(p), qAm), n = Math.hypot(...rS); return ((b[0] * rS[0] + b[1] * rS[1] + b[2] * rS[2]) / n) * 1e15; };
    // The electrodes and the helmet's sensors.
    const names = [...h.electrodes.keys()].filter((n) => !['Nz', 'LPA', 'RPA', 'A1', 'A2', 'M1', 'M2'].includes(n));
    this.electrodes = names.map((n) => ({ name: n, position: h.electrodes.get(n) }));
    const ev = this.electrodes.map((e) => eegAt(e.position.map((x, a) => x - c[a])));
    const sensors = this.helmetSensors;
    const mv = sensors.positions.map((s) => megAt(s));
    const eMax = Math.max(...ev.map(Math.abs)) || 1, mRaw = Math.max(...mv.map(Math.abs));
    const silent = mRaw < 1e-3;
    const mMax = silent ? 1 : mRaw;
    if (silent) for (let i = 0; i < mv.length; i++) mv[i] = 0;
    m.eegMax = eMax; m.megMax = silent ? 0 : mRaw; m.silent = silent;
    // The scalp map of EEG, on the 3-D scalp when asked, and on the disc.
    this.spline = new ScalpSpline(this.electrodes.map((e) => e.position.map((x, a) => x - c[a]))).fit(ev);
    const colourE = (x) => TABLES.diverging[Math.max(0, Math.min(255, Math.round(((x + eMax) / (2 * eMax)) * 255)))];
    this.paintScalp(this.spline, colourE);
    this.showElectrodes(ev, colourE);
    // The helmet's field, at its vertices.
    if (helmet) {
      const pos = helmet.mesh.geometry.attributes.position.array, rgb = new Float32Array(pos.length);
      for (let i = 0; i < pos.length; i += 3) {
        const b = silent ? 0 : megAt([pos[i], pos[i + 1], pos[i + 2]]);
        const t = Math.max(0, Math.min(255, Math.round(((b + mMax) / (2 * mMax)) * 255)));
        const col = TABLES.diverging[t];
        rgb[i] = col[0] / 255; rgb[i + 1] = col[1] / 255; rgb[i + 2] = col[2] / 255;
      }
      v.paintRGB(helmet, rgb);
    }
    // The estimates each modality would make of the source: sLORETA, noiseless.
    this.megEstimate(ev, mv);
    this.scalpMap = { values: ev, range: [-eMax, eMax], signed: true, unit: 'µV', title: 'EEG potential' };
    this.megDisc = { values: mv, positions: sensors.positions, range: [-mMax, mMax] };
    this.ctx.renderCard();
    this.legend();
  }

  dipoleOrientation() {
    const m = this.meg, n = m.source.normal, p = m.source.position, c = this.head.centre;
    const radial = p.map((x, a) => x - c[a]);
    const rl = Math.hypot(...radial);
    const r = radial.map((x) => x / rl);
    if (m.orientation === 'radial') return r;
    if (m.orientation === 'tangential') {
      // The part of the normal across the radius, or else any direction across it.
      let t = n.map((x, a) => x - (n[0] * r[0] + n[1] * r[1] + n[2] * r[2]) * r[a]);
      if (Math.hypot(...t) < 0.2) t = [r[1], -r[0], 0];
      const tl = Math.hypot(...t);
      return t.map((x) => x / tl);
    }
    return n;
  }

  buildHelmet() {
    const h = this.head, c = h.centre, R = h.radius + 30;
    // A cap of a sphere: rings of latitude down to 20 degrees below the equator.
    const rings = 40, segs = 96, pos = [], nrm = [], idx = [];
    const zmin = Math.sin(-0.35);
    for (let i = 0; i <= rings; i++) {
      const z = 1 - (i / rings) * (1 - zmin), r = Math.sqrt(Math.max(0, 1 - z * z));
      for (let j = 0; j <= segs; j++) {
        const t = (j / segs) * 2 * Math.PI, d = [r * Math.cos(t), r * Math.sin(t), z];
        pos.push(c[0] + R * d[0], c[1] + R * d[1], c[2] + R * d[2]);
        nrm.push(...d);
      }
    }
    for (let i = 0; i < rings; i++) for (let j = 0; j < segs; j++) {
      const a = i * (segs + 1) + j, b = a + segs + 1;
      idx.push(a, b, a + 1, a + 1, b, b + 1);
    }
    const helmet = this.viewer.addFreeMesh('helmet', { positions: Float32Array.from(pos), normals: Float32Array.from(nrm), indices: Uint32Array.from(idx) }, { kind: 'skin', colour: '#e8dce8', pickable: false, opacity: 0.34, render: 5 });
    helmet.material.side = 2;
    this.helmet = true;
    this.helmetSensors = idealHelmet(c, R, { count: 248, lowest: -0.35 });
  }

  megEstimate(ev, mv) {
    const m = this.meg, h = this.head;
    const key = 'meg-ops';
    if (!this.megOps) {
      this.megOps = 'building';
      this.ctx.toast('Computing the lead fields of EEG and MEG…', 12000);
      const sources = this.space.positions;
      const elec = this.electrodes.map((e) => h.normalise(e.position));
      Promise.all([
        leadField({ kind: 'eeg', model: { radius: h.radius, radii: THREE_SHELLS.radii, conductivities: THREE_SHELLS.conductivities, terms: 100 }, centre: h.centre, sensors: elec, sources: sources.map((p) => h.normalise(p)), normals: this.space.normals }),
        leadField({ kind: 'meg', centre: h.centre, sensors: this.helmetSensors.positions, orientations: this.helmetSensors.orientations, sources, normals: this.space.normals }),
      ]).then(([Le, Lm]) => {
        this.megOps = {
          eeg: inverseOperator(Le, elec.length, this.space.count, { method: 'sLORETA', snr: 10, reference: 'average' }),
          meg: inverseOperator(Lm, this.helmetSensors.positions.length, this.space.count, { method: 'sLORETA', snr: 10, reference: null }),
        };
        if (this.modality === 'meg') this.refreshMeg();
      });
      void key;
      return;
    }
    if (this.megOps === 'building') return;
    const op = this.megOps[m.reconstruct];
    const data = m.reconstruct === 'eeg' ? ev : mv;
    const J = op.apply(data).map(Math.abs);
    let best = 0;
    for (let s = 1; s < J.length; s++) if (J[s] > J[best]) best = s;
    let mx = J[best] || 1;
    m.estimate = { peak: best, error: Math.hypot(...this.space.positions[best].map((x, a) => x - m.source.position[a])) };
    // A dipole the radial way is invisible to MEG: then the estimate is noise.
    if (m.reconstruct === 'meg' && Math.max(...mv.map(Math.abs)) < 1e-3) mx = Infinity;
    this.cortexMap = { kind: 'sources', values: J, scale: mx, signed: false, threshold: 0.45, rampPos: 'ember' };
    this.paintCortex();
  }

  // ---- picking and hovering ------------------------------------------------------------------------

  vertexAt(hit) {
    const pos = hit.inst.mesh.geometry.attributes.position.array;
    let best = hit.face.a, bd = Infinity;
    for (const v of [hit.face.a, hit.face.b, hit.face.c]) {
      const d = (pos[3 * v] - hit.point[0]) ** 2 + (pos[3 * v + 1] - hit.point[1]) ** 2 + (pos[3 * v + 2] - hit.point[2]) ** 2;
      if (d < bd) { bd = d; best = v; }
    }
    return best;
  }

  // The source nearest a picked point of the cortex.
  sourceAt(hit) {
    const p = hit.point;
    let best = -1, bd = Infinity;
    this.space.positions.forEach((q, s) => { if (this.space.side[s] !== hit.inst.side) return; const d = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2; if (d < bd) { bd = d; best = s; } });
    return best;
  }

  pick(event) {
    const v = this.viewer;
    const mi = v.pickMarker(event.clientX, event.clientY);
    if (mi >= 0 && this.electrodes && this.electrodes[mi]) {
      this.chosenChannel = this.electrodes[mi].name;
      this.ctx.renderCard();
      return true;
    }
    const hit = v.pick(event.clientX, event.clientY);
    if (!hit || hit.inst.id !== 'cerebral-cortex' || !hit.face) return false;
    if (this.modality === 'meg') {
      const s = this.sourceAt(hit);
      if (s >= 0) { this.meg.source = { s, position: this.space.positions[s], normal: this.space.normals[s] }; this.refreshMeg(); }
      return true;
    }
    if (this.modality === 'eeg' && this.eeg.rec) {
      const s = this.sourceAt(hit);
      this.eeg.probe = s;
      this.refreshEeg();
      const at = this.space.structure[s];
      this.ctx.toast(`The probe reads the estimate at the ${at ? this.ctx.structure(at.id).name.toLowerCase() : 'cortex'}${at ? `, ${at.side}` : ''}: its trace is the first, in gold.`, 5000);
      return true;
    }
    if (this.modality === 'fmri' && this.fmri.volume) {
      this.ctx.setPoint(hit.point);
      const vv = this.vertexAt(hit), i = this.space.instances.indexOf(hit.inst);
      const z = this.fmri.vertex[i][vv];
      const at = this.anatomyAt(hit.point);
      this.ctx.toast(`z = ${tidy(z)} at ${hit.point.map((x) => signed(x)).join(', ')} mm${at ? `, ${this.ctx.structure(at.id).name.toLowerCase()}` : ''}`, 5000);
      return true;
    }
    return false;
  }

  hoverInfo(pointer) {
    const v = this.viewer;
    const mi = v.pickMarker(pointer.x, pointer.y);
    if (mi >= 0 && this.electrodes && this.electrodes[mi]) {
      const val = this.markerValues ? this.markerValues[mi] : null;
      return { name: this.electrodes[mi].name, sub: val !== null && val !== undefined ? `${tidy(val)} ${this.scalpMap ? this.scalpMap.unit : ''}` : 'an electrode', colour: '#f6dfb8' };
    }
    const hit = v.pick(pointer.x, pointer.y);
    if (!hit || hit.inst.id !== 'cerebral-cortex' || !hit.face) return null;
    const i = this.space.instances.indexOf(hit.inst), vv = this.vertexAt(hit);
    const label = hit.inst.labels ? hit.inst.labels[vv] : 0;
    const at = label ? structureOfKey(`ho:${label}`) : null;
    const name = at ? this.ctx.structure(at.id).name : 'Cerebral cortex';
    let sub = at ? `${at.side} hemisphere` : '';
    const map = this.cortexMap;
    if (map && i >= 0) {
      const val = map.kind === 'vertices' ? map.values[i][vv] : this.buffers[i].values[vv];
      sub += map.kind === 'vertices' ? ` · z ${tidy(val)}` : ` · ${tidy((Math.abs(val) / map.scale) * 100)} % of the scale`;
    }
    return { name, sub, colour: at ? this.ctx.colourOf(at.id) : '#e3cfc9' };
  }

  // ---- the card --------------------------------------------------------------------------------------

  cardHtml() {
    if (this.modality === 'fmri') return this.fmriCard();
    if (this.modality === 'meg') return this.megCard();
    return this.eegCard();
  }

  eegCard() {
    const e = this.eeg, rec = e.rec;
    setTimeout(() => this.bindCard(), 0);
    if (!rec) return '<div class="card-empty"><h2>Signals</h2><p>Reading the recording…</p></div>';
    const facts = [['Channels', `${rec.used.length} <small>of ${rec.channels.length}</small>`], ['Rate', `${tidy(rec.rate)} <small>Hz</small>`], ['Length', `${tidy(rec.duration)} <small>s</small>`]];
    const sp = this.spectrumOf();
    const chosen = this.chosenChannel && sp.per.get(this.chosenChannel);
    const bandPow = chosen ? bandPower(chosen, SIGNAL_BANDS[e.band].lo, SIGNAL_BANDS[e.band].hi) : null;
    const cinfo = e.view === 'contrast' && e.contrastInfo ? `<p class="note">${esc(e.contrastInfo.a)}: ${tidy(e.contrastInfo.sa)} s in ${e.contrastInfo.na} ${e.contrastInfo.na === 1 ? 'span' : 'spans'}; ${esc(e.contrastInfo.b)}: ${tidy(e.contrastInfo.sb)} s in ${e.contrastInfo.nb}. Warm where the first has more power, cool where the second has.</p>` : '';
    return `<div class="card-head"><div class="over">EEG · ${esc(rec.source || rec.format)}</div>
<h2><span class="sw" style="background:#f6c47e"></span><span>${esc(rec.title)}</span></h2>
<div class="latin">${esc(rec.note || '')}</div></div>
<div class="card-body">
<div class="facts">${facts.map(([k, v]) => `<div><div class="k">${esc(k)}</div><div class="v">${v}</div></div>`).join('')}</div>
<div class="sect"><h4>On the scalp</h4><div class="topo-wrap"><canvas class="topo"></canvas></div>${cinfo}</div>
<div class="sect"><h4>Strongest sources</h4><div class="chips peaks"></div><p class="note">Local maxima of the estimate on the cortex. Choose one to read its anatomy; click the cortex to probe the estimate’s time course anywhere.</p></div>
<div class="sect"><h4>Spectrum${chosen ? ` at ${esc(this.chosenChannel)}` : ', the mean of the channels'}</h4>${this.spectrumSvg()}
<p class="note">${sp.alpha ? `The alpha peak lies at ${sp.alpha.toFixed(1)} Hz. ` : ''}${chosen ? `${esc(SIGNAL_BANDS[e.band].name.split(',')[0])} power at ${esc(this.chosenChannel)}: ${tidy(bandPow)} µV². ` : 'Click an electrode for its own power. '}Welch’s method, segments of 2 s.</p></div>
<div class="sect"><h4>How the map is made</h4><p>The channels are referenced to their average and filtered to the band. On the scalp, their values are interpolated by spherical splines (Perrin et al. 1989). On the cortex, ${esc(e.method)} estimates the currents at ${this.space.count} places normal to the cortex, through a head of three spheres: brain, skull of 1/${e.skull} of the brain’s conductivity, scalp (Rush and Driscoll 1968), with the template’s head mapped onto them. It is an estimate from a template, blurred by the skull and by the few electrodes, not an image of the brain.</p></div>
${rec.unused.length ? `<p class="note">Not placed: ${esc(rec.unused.slice(0, 12).join(', '))}${rec.unused.length > 12 ? '…' : ''}.</p>` : ''}
<div class="card-refs"><p>Schalk, G. et al. (2004), BCI2000: a general-purpose brain-computer interface (BCI) system, <i>IEEE Transactions on Biomedical Engineering</i> 51, 1034–1043; Goldberger, A. L. et al. (2000), PhysioBank, PhysioToolkit, and PhysioNet, <i>Circulation</i> 101, e215–e220: the example recordings, ODC-By 1.0.</p>
<p>Pascual-Marqui, R. D. (2002), Standardized low-resolution brain electromagnetic tomography (sLORETA): technical details, <i>Methods and Findings in Experimental and Clinical Pharmacology</i> 24D, 5–12. Hämäläinen, M. S. and Ilmoniemi, R. J. (1994), Interpreting magnetic fields of the brain: minimum norm estimates, <i>Medical and Biological Engineering and Computing</i> 32, 35–42. Dale, A. M. et al. (2000), Dynamic statistical parametric mapping, <i>Neuron</i> 26, 55–67.</p></div></div>`;
  }

  fmriCard() {
    const f = this.fmri;
    setTimeout(() => this.bindCard(), 0);
    if (!f.volume) return '<div class="card-empty"><h2>Maps of activity</h2><p>Reading the map…</p></div>';
    const ex = f.example || {};
    const rows = (f.peaks || []).map((p, i) => {
      const st = p.where && this.ctx.structure(p.where.id);
      return `<tr data-peak-pos="${p.position.map((x) => x.toFixed(1)).join(',')}"><td>${i + 1}</td><td class="${p.value < 0 ? 'neg' : 'pos'}">${tidy(p.value)}</td><td>${p.position.map((x) => signed(x)).join(', ')}</td><td>${st ? `<button type="button" data-select="${esc(p.where.id)}" data-side="${esc(p.where.side === 'mid' ? '' : p.where.side)}">${esc(st.short || st.name)}${p.where.side !== 'mid' ? ` <small>${esc(p.where.side[0].toUpperCase())}</small>` : ''}</button>` : '<span class="note">—</span>'}</td></tr>`;
    }).join('');
    const vox = Math.abs(f.volume.affine[0]).toFixed(1);
    return `<div class="card-head"><div class="over">fMRI · ${esc(ex.id && ex.id.startsWith('file:') ? 'a file' : 'Individual Brain Charting, subject 1')}</div>
<h2><span class="sw" style="background:${css(ramp('heat', 0.6))}"></span><span>${esc(ex.name || 'A map')}</span></h2>
<div class="latin">${esc(ex.note || '')}</div></div>
<div class="card-body">
<div class="facts"><div><div class="k">Voxel</div><div class="v">${vox} <small>mm</small></div></div><div><div class="k">Smoothing</div><div class="v">${f.fwhm ? `${f.fwhm} <small>mm</small>` : 'none'}</div></div><div><div class="k">Largest</div><div class="v">${tidy(f.max)} <small>z</small></div></div></div>
<div class="sect"><h4>Peaks</h4><table class="peak-table"><tr><th></th><th>z</th><th>MNI, mm</th><th>Lies in</th></tr>${rows || '<tr><td colspan="4" class="note">Nothing beyond the threshold.</td></tr>'}</table>
<p class="note">Local maxima at least 12 mm apart, named by the structure of the anatomy at them, or the nearest within 5 mm. Click a row to go there in the sections.</p></div>
<div class="sect"><h4>How it is drawn</h4>${f.fwhm ? `<p>A map of one person, unsmoothed, is mostly specks: noise passes any threshold at single voxels. It is smoothed here as a group analysis would be, by a Gaussian of ${f.fwhm} mm at half its height, which averages each voxel with its neighbours within a few millimetres; the noise shrinks more than the activity, which spans centimetres. The smoothed map is then divided by its spread over the brain (the median absolute deviation, scaled to a standard deviation), so that z again counts standard deviations of the noise.</p>` : ''}<p>The map is read at every point of the cortex along its inward normal, from 1 mm outside to 6 mm within, keeping the strongest value; the deep structures read it at their surfaces. It is drawn where |z| reaches the threshold: warm for activity above the baseline, cool for below. The template’s cortex is not the subject’s, so the map is placed by its coordinates, as on any standard brain.</p></div>
<div class="card-refs"><p>Pinho, A. L. et al. (2018), Individual Brain Charting, a high-resolution fMRI dataset for cognitive mapping, <i>Scientific Data</i> 5, 180105; the maps via NeuroVault (Gorgolewski et al. 2015), CC0.</p></div></div>`;
  }

  megCard() {
    const m = this.meg;
    setTimeout(() => this.bindCard(), 0);
    if (!m.source) return '<div class="card-empty"><h2>EEG and MEG</h2><p>Placing a source…</p></div>';
    const at = this.space.structure[m.source.s], st = at && this.ctx.structure(at.id);
    const est = m.estimate;
    const radial = m.orientation === 'radial';
    return `<div class="card-head"><div class="over">One source, two ways of seeing it</div>
<h2><span class="sw" style="background:#f6c47e"></span><span>A dipole in the ${esc(st ? st.name.toLowerCase() : 'cortex')}</span></h2>
<div class="latin">${esc(at ? `${at.side} hemisphere · ` : '')}${m.strength} nA·m, ${esc({ normal: 'normal to the cortex', tangential: 'tangential to the head', radial: 'radial to the head' }[m.orientation])}</div></div>
<div class="card-body">
<div class="topo-pair"><div><canvas class="topo"></canvas></div><div><canvas class="topo meg"></canvas></div></div>
<div class="facts"><div><div class="k">EEG, largest</div><div class="v">${tidy(m.eegMax)} <small>µV</small></div></div><div><div class="k">MEG, largest</div><div class="v">${m.silent ? 'none' : `${tidy(m.megMax)} <small>fT</small>`}</div></div><div><div class="k">${esc(m.reconstruct.toUpperCase())} misplaces it by</div><div class="v">${est && this.megOps !== 'building' ? `${tidy(est.error)} <small>mm</small>` : '…'}</div></div></div>
<div class="sect"><h4>What each sees</h4><p>A current in the cortex makes a potential on the scalp and a magnetic field outside the head. The potential spreads through the skull, which conducts poorly and smears it; the field passes the skull untouched. Outside a sphere, a current pointing along the radius makes no field at all (Sarvas 1987): MEG sees the currents of the walls of the sulci, tangential to the head, and is blind to those of the crowns of the gyri. EEG sees both. Turn the source to see it${radial ? ': radial now, its field is nothing, while its potential is a single round pole.' : '.'}</p>
<p>The cortex shows where sLORETA would place the source from the ${esc(m.reconstruct.toUpperCase())} alone, without noise: ${est ? `${tidy(est.error)} mm from where it is` : '…'}${radial && m.reconstruct === 'meg' ? ', which, with no field to go by, is nowhere in particular' : ''}. Click the cortex to move the source.</p></div>
<div class="card-refs"><p>Sarvas, J. (1987), Basic mathematical and electromagnetic concepts of the biomagnetic inverse problem, <i>Physics in Medicine and Biology</i> 32, 11–22. Hämäläinen, M. et al. (1993), Magnetoencephalography: theory, instrumentation, and applications to noninvasive studies of the working human brain, <i>Reviews of Modern Physics</i> 65, 413–497. Rush, S. and Driscoll, D. A. (1968), Current distribution in the brain from surface electrodes, <i>Anesthesia and Analgesia</i> 47, 717–723.</p></div></div>`;
  }

  bindCard() {
    const el = document.getElementById('card');
    this.renderPeaks();
    if (this.modality === 'eeg') this.drawTopoCard();
    if (this.modality === 'meg' && this.scalpMap) {
      const c = this.head.centre;
      const canvas = el.querySelector('canvas.topo.meg');
      if (canvas && this.megDisc) {
        const d = this.megDisc, [lo, hi] = d.range;
        const sp = new ScalpSpline(d.positions.map((p) => p.map((x, a) => x - c[a]))).fit(d.values);
        drawTopo(canvas, { field: (dir) => sp.at(dir), colour: (v) => TABLES.diverging[Math.max(0, Math.min(255, Math.round(((v - lo) / (hi - lo || 1)) * 255)))], electrodes: [], contours: (hi - lo) / 10 || null, size: 160, title: 'MEG field', unit: 'fT' });
      }
      const eegCanvas = el.querySelector('canvas.topo:not(.meg)');
      if (eegCanvas) { const sm = this.scalpMap, [lo, hi] = sm.range; drawTopo(eegCanvas, { field: (dir) => this.spline.at(dir), colour: (v) => TABLES.diverging[Math.max(0, Math.min(255, Math.round(((v - lo) / (hi - lo || 1)) * 255)))], electrodes: this.electrodes.map((x) => ({ name: x.name, direction: x.position.map((v, a) => v - c[a]) })), contours: (hi - lo) / 10 || null, size: 160, title: 'EEG potential', unit: 'µV' }); }
    }
    el.querySelectorAll('[data-peak]').forEach((b) => b.addEventListener('click', (ev) => {
      ev.stopPropagation();
      const s = Number(b.dataset.peak);
      const p = this.space.positions[s];
      this.ctx.focusBox(p.map((x) => x - 14), p.map((x) => x + 14), 2.2);
      this.eeg.probe = s;
      this.refreshEeg();
    }));
    el.querySelectorAll('tr[data-peak-pos]').forEach((tr) => tr.addEventListener('click', (ev) => {
      if (ev.target.closest('button')) return;
      const p = tr.dataset.peakPos.split(',').map(Number);
      this.ctx.setPoint(p);
      this.ctx.state.show = 'sections';
      this.ctx.renderToolbar();
      document.querySelector('[data-show=sections]')?.click();
    }));
  }

  // ---- the panel and the toolbar -----------------------------------------------------------------------

  panel() {
    const tabs = `<div class="pad seg modality">${[['eeg', 'EEG'], ['meg', 'EEG and MEG'], ['fmri', 'fMRI']].map(([k, n]) => `<button type="button" class="opt" data-modality="${k}" aria-pressed="${this.modality === k}">${n}</button>`).join('')}</div>`;
    const drop = '<p class="pad note drop-note">Or drop a file on the brain.</p>';
    if (this.modality === 'fmri') {
      const f = this.fmri;
      return `<h3>Signals</h3>${tabs}
<h3>A map of activity</h3><div class="pad atlas-choices">${FMRI_EXAMPLES.map((x) => `<button type="button" class="atlas-choice" data-fmri="${x.id}" aria-pressed="${f.example && f.example.id === x.id}"><b>${esc(x.name)}</b><span class="what">${esc(x.note)}</span></button>`).join('')}</div>
<div class="pad row-btns"><label class="btn file">Open a map… <input type="file" accept=".nii,.gz" data-open hidden></label></div>${drop}
<h3>Smoothing</h3><div class="pad seg">${[[0, 'None'], [4, '4 mm'], [6, '6 mm'], [8, '8 mm']].map(([k, n]) => `<button type="button" class="opt" data-fwhm="${k}" aria-pressed="${f.fwhm === k}">${n}</button>`).join('')}</div>
<h3>Threshold</h3><div class="pad"><div class="slider"><span>|z| ≥</span><input type="range" min="1.5" max="10" step="0.1" value="${f.threshold}" data-fthr><span class="value">${f.threshold.toFixed(1)}</span></div>
<div class="seg">${[['both', 'Both'], ['positive', 'Above'], ['negative', 'Below']].map(([k, n]) => `<button type="button" class="opt" data-fsign="${k}" aria-pressed="${f.sign === k}">${n}</button>`).join('')}</div></div>`;
    }
    if (this.modality === 'meg') {
      const m = this.meg;
      return `<h3>Signals</h3>${tabs}
<p class="pad note">One current in the cortex, of ${m.strength} nA·m, about what a few square centimetres of cortex make when they act together. Click the cortex to place it.</p>
<h3>Its direction</h3><div class="pad seg">${[['normal', 'Normal to the cortex'], ['tangential', 'Tangential'], ['radial', 'Radial']].map(([k, n]) => `<button type="button" class="opt" data-orient="${k}" aria-pressed="${m.orientation === k}">${n}</button>`).join('')}</div>
<h3>Estimated from</h3><div class="pad seg">${[['eeg', 'EEG'], ['meg', 'MEG']].map(([k, n]) => `<button type="button" class="opt" data-recon="${k}" aria-pressed="${m.reconstruct === k}">${n}</button>`).join('')}</div>
<h3>Its strength</h3><div class="pad"><div class="slider"><span>nA·m</span><input type="range" min="5" max="100" step="5" value="${m.strength}" data-strength><span class="value">${m.strength}</span></div></div>
<p class="pad note">The helmet is an idealised array of 248 radial magnetometers, 3 cm beyond the scalp.</p>`;
    }
    const e = this.eeg;
    const conds = this.conditions();
    const ex = EEG_EXAMPLES.map((x) => `<button type="button" class="atlas-choice" data-eeg="${x.id}" aria-pressed="${e.rec && e.rec.id === x.id}"><b>${esc(x.name)}</b><span class="what">${esc(x.note)}</span></button>`).join('');
    const opt = (sel, list) => list.map((c) => `<option value="${esc(c.key)}"${c.key === sel ? ' selected' : ''}>${esc(c.name)}</option>`).join('');
    return `<h3>Signals</h3>${tabs}
<h3>A recording</h3><div class="pad atlas-choices">${ex}</div>
<div class="pad row-btns"><label class="btn file">Open a recording… <input type="file" accept=".edf,.bdf,.rec,.vhdr,.vmrk,.eeg,.dat,.csv,.tsv,.txt" multiple data-open hidden></label></div>${drop}
<h3>The band</h3><div class="pad seg">${Object.entries(SIGNAL_BANDS).map(([k, b]) => `<button type="button" class="opt" data-band="${k}" aria-pressed="${e.band === k}">${esc(b.name.split(',')[0])}${k === 'broad' ? '' : ` <small>${b.lo}–${b.hi}</small>`}</button>`).join('')}</div>
<h3>The map</h3><div class="pad seg">${[['potential', 'At the cursor'], ['power', 'Band power'], ['contrast', 'Contrast']].map(([k, n]) => `<button type="button" class="opt" data-view="${k}" aria-pressed="${e.view === k}">${n}</button>`).join('')}</div>
${e.view === 'power' ? `<div class="pad"><div class="slider"><span>over</span><input type="range" min="0.5" max="20" step="0.5" value="${e.window}" data-window><span class="value">${e.window} s</span></div></div>` : ''}
${e.view === 'contrast' && conds.length >= 2 ? `<div class="pad contrast"><label>Of <select data-ca>${opt(e.contrast?.a, conds)}</select></label><label>over <select data-cb>${opt(e.contrast?.b, conds)}</select></label><p class="note">Read another example to compare with it, as eyes closed over eyes open.</p></div>` : ''}
<h3>The sources</h3><div class="pad seg">${['sLORETA', 'dSPM', 'MNE'].map((k) => `<button type="button" class="opt" data-method="${k}" aria-pressed="${e.method === k}">${k}</button>`).join('')}</div>
<div class="pad"><div class="slider"><span>shown above</span><input type="range" min="0" max="0.9" step="0.05" value="${e.threshold}" data-thr><span class="value">${Math.round(e.threshold * 100)} %</span></div>
<label class="field"><span>Skull, a fraction of the brain’s conductivity</span><select data-skull>${[15, 25, 40, 80].map((k) => `<option value="${k}"${e.skull === k ? ' selected' : ''}>1/${k}</option>`).join('')}</select></label></div>`;
  }

  bind(left) {
    left.querySelectorAll('[data-modality]').forEach((b) => b.addEventListener('click', () => this.switchTo(b.dataset.modality)));
    left.querySelectorAll('[data-eeg]').forEach((b) => b.addEventListener('click', () => this.loadExample(EEG_EXAMPLES.find((x) => x.id === b.dataset.eeg))));
    left.querySelectorAll('[data-fmri]').forEach((b) => b.addEventListener('click', () => this.loadFmriExample(FMRI_EXAMPLES.find((x) => x.id === b.dataset.fmri)).then(() => this.ctx.renderLeft())));
    left.querySelectorAll('[data-open]').forEach((inp) => inp.addEventListener('change', () => { if (inp.files.length) this.openFiles(inp.files); }));
    left.querySelectorAll('[data-band]').forEach((b) => b.addEventListener('click', () => { this.eeg.band = b.dataset.band; this.eeg.scaleJ = null; this.turnToPeak = true; this.ctx.renderLeft(); this.refreshEeg(); }));
    left.querySelectorAll('[data-view]').forEach((b) => b.addEventListener('click', () => { this.eeg.view = b.dataset.view; this.turnToPeak = true; this.ctx.renderLeft(); this.refreshEeg(); }));
    left.querySelectorAll('[data-method]').forEach((b) => b.addEventListener('click', () => { this.eeg.method = b.dataset.method; this.ctx.renderLeft(); this.refreshEeg(); }));
    const thr = left.querySelector('[data-thr]');
    if (thr) thr.addEventListener('input', () => { this.eeg.threshold = Number(thr.value); thr.nextElementSibling.textContent = `${Math.round(this.eeg.threshold * 100)} %`; if (this.cortexMap) { this.cortexMap.threshold = this.eeg.threshold; this.paintCortex(); this.legend(); } });
    const win = left.querySelector('[data-window]');
    if (win) win.addEventListener('input', () => { this.eeg.window = Number(win.value); win.nextElementSibling.textContent = `${this.eeg.window} s`; this.updateMap(); });
    const skull = left.querySelector('[data-skull]');
    if (skull) skull.addEventListener('change', () => { this.eeg.skull = Number(skull.value); this.useRecording(this.eeg.rec); });
    const ca = left.querySelector('[data-ca]'), cb = left.querySelector('[data-cb]');
    if (ca) ca.addEventListener('change', () => { this.eeg.contrast.a = ca.value; this.refreshEeg(); });
    if (cb) cb.addEventListener('change', () => { this.eeg.contrast.b = cb.value; this.refreshEeg(); });
    const fthr = left.querySelector('[data-fthr]');
    if (fthr) fthr.addEventListener('change', () => { this.fmri.threshold = Number(fthr.value); this.refreshFmri(); });
    if (fthr) fthr.addEventListener('input', () => { fthr.nextElementSibling.textContent = Number(fthr.value).toFixed(1); });
    left.querySelectorAll('[data-fwhm]').forEach((b) => b.addEventListener('click', () => { this.fmri.fwhm = Number(b.dataset.fwhm); this.ctx.renderLeft(); if (this.fmri.raw) this.smoothFmri(); }));
    left.querySelectorAll('[data-fsign]').forEach((b) => b.addEventListener('click', () => { this.fmri.sign = b.dataset.fsign; this.ctx.renderLeft(); this.refreshFmri(); }));
    left.querySelectorAll('[data-orient]').forEach((b) => b.addEventListener('click', () => { this.meg.orientation = b.dataset.orient; this.ctx.renderLeft(); this.refreshMeg(); }));
    left.querySelectorAll('[data-recon]').forEach((b) => b.addEventListener('click', () => { this.meg.reconstruct = b.dataset.recon; this.ctx.renderLeft(); this.refreshMeg(); }));
    const strength = left.querySelector('[data-strength]');
    if (strength) strength.addEventListener('change', () => { this.meg.strength = Number(strength.value); strength.nextElementSibling.textContent = strength.value; this.refreshMeg(); });
  }

  toolbar() {
    const opt = (action, label, pressed) => `<button type="button" class="opt" data-mode-action="${action}" aria-pressed="${pressed}">${label}</button>`;
    const state = this.ctx.state;
    if (this.modality === 'fmri') {
      const show = `<div class="group"><span class="caption">Show</span><div class="row"><button type="button" class="opt" data-show="3d" aria-pressed="${state.show !== 'sections'}">In the round</button><button type="button" class="opt" data-show="sections" aria-pressed="${state.show === 'sections'}">In sections</button></div></div>`;
      if (state.show === 'sections') return `${show}<div class="group"><span class="caption">Main section</span><div class="row">${['axial', 'coronal', 'sagittal'].map((p) => `<button type="button" class="opt" data-main="${p}" aria-pressed="${state.mainPlane === p}">${p[0].toUpperCase() + p.slice(1)}</button>`).join('')}</div></div>`;
      return `${show}<div class="group"><span class="caption">Cortex</span><div class="row">${opt('cortex:solid', 'Solid', state.cortexStyle === 'solid')}${opt('cortex:glass', 'Glass', state.cortexStyle === 'glass')}</div></div>`;
    }
    return `<div class="group"><span class="caption">Scalp</span><div class="row">${opt('scalp:glass', 'Glass', this.scalpStyle === 'glass')}${opt('scalp:map', 'Map', this.scalpStyle === 'map')}${opt('scalp:hidden', 'Hidden', this.scalpStyle === 'hidden')}</div></div>
<div class="group"><span class="caption">Electrodes</span><div class="row">${opt('names', 'Names', this.names)}</div></div>`;
  }

  toolbarAction(action) {
    const [what, value] = action.split(':');
    if (what === 'scalp') {
      this.scalpStyle = value;
      this.ensureScalp();
      if (this.modality === 'meg') this.refreshMeg(); else this.updateMap();
    }
    if (what === 'names') { this.names = !this.names; this.placeLabels(); }
    if (what === 'cortex') { this.ctx.state.cortexStyle = value; this.ctx.applyStyles(); this.refreshFmri(); }
    this.ctx.renderToolbar();
  }
}
