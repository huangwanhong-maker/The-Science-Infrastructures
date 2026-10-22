// The Head mode: the template's scalp around the brain, the electrodes of the 10–20
// or 10–10 system placed on it by the arcs between the nasion, the inion and the
// ears, the brain and head scaled to a subject's measurements, and values measured
// at the electrodes interpolated over the scalp by spherical splines and carried
// radially onto the cortex beneath.
//
// The geometry is the neuroscience domain's (head.js): maskScalp, tenTwenty,
// scaleToHead, sphericalSpline (Jasper 1958; Oostenveld and Praamstra 2001; Perrin
// et al. 1989). The scalp's surface is extracted from the template's head mask in a
// worker, by the visualization layer's isosurface.

import { loadVolume, FILES } from '../data/anatomy.js';
import { channelValues, valueRange, tidy } from '../ui/data.js';
import { colourAt, cssGradient } from '../render/colormaps.js';
import { esc } from '../ui/card.js';
import {
  maskScalp, tenTwenty, scaleToHead, sphericalSpline, apply, invert, structureOfKey, TEN_TWENTY_NAMES_1985, robustRange,
} from '../../../../domains/neuroscience/js/index.js';

// Landmarks on the template's head: scalp points in the directions of the usual positions,
// found by casting rays from (0, −18, 0) through the head of the MNI 152 2009c template on
// 2026-10-22. Provisional: digitised positions should replace them.
const LANDMARKS = { nasion: [0, 88.5, -45], inion: [0, -121.5, -33], lpa: [-84.5, -21, -50], rpa: [85.5, -21, -50], vertex: [0, -18, 105.5] };
const CENTRE = [0, -18, 0];
const LANDMARK_NAMES = { Nz: 'Nasion', Iz: 'Inion', LPA: 'Left preauricular point', RPA: 'Right preauricular point' };
const REGION_LETTERS = { Fp: 'frontopolar', AF: 'anterior frontal', F: 'frontal', FC: 'frontocentral', FT: 'frontotemporal', C: 'central', T: 'temporal', TP: 'temporoparietal', CP: 'centroparietal', P: 'parietal', PO: 'parieto-occipital', O: 'occipital', I: 'inion', N: 'nasion' };
const EXAMPLE = `# An invented pattern, for trying the map: strongest at the back of the head,
# as the alpha rhythm of a person resting with closed eyes is.
Fp1, 0.8\nFp2, 0.9\nF7, 1.2\nF3, 1.6\nFz, 1.7\nF4, 1.6\nF8, 1.1\nT7, 1.9\nC3, 2.6\nCz, 2.8\nC4, 2.7\nT8, 2.0\nP7, 4.1\nP3, 4.9\nPz, 5.3\nP4, 5.0\nP8, 4.2\nO1, 6.4\nO2, 6.6`;

const signed = (x) => (x < 0 ? `−${-x}` : `${x}`);
const coords = (c) => c.map((x) => signed(Math.round(x))).join(', ');
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const unit = (v) => { const n = Math.hypot(...v) || 1; return [v[0] / n, v[1] / n, v[2] / n]; };

let worker = null, requests = 0;
const pending = new Map();
function surface(message) {
  if (!worker) {
    worker = new Worker(new URL('../workers/surface.worker.js', import.meta.url), { type: 'module' });
    worker.onmessage = (e) => { const p = pending.get(e.data.id); if (!p) return; pending.delete(e.data.id); if (e.data.error) p.reject(new Error(e.data.error)); else p.resolve(e.data); };
  }
  const id = ++requests;
  return new Promise((resolve, reject) => { pending.set(id, { resolve, reject }); worker.postMessage({ id, ...message }); });
}

export { LANDMARKS, CENTRE };

export class HeadMode {
  constructor(ctx) {
    this.ctx = ctx;
    this.system = '10-20';
    this.scalpStyle = 'glass';
    this.names = true;
    this.template = scaleToHead(LANDMARKS, {}, CENTRE).template;
    this.measurements = { ...this.template };
    this.affine = null;            // template -> subject, when scaled
    this.scales = [1, 1, 1];
    this.channelsText = '';
    this.report = '';
    this.spline = null;
    this.positions = [];           // electrodes in subject space
    this.chosen = null;
  }

  get app() { return this.ctx.app; }

  async enter() {
    this.active = true;
    const app = this.app;
    app.plateOff = true;
    // The brain is seen whole beneath the scalp.
    this.savedCortex = this.ctx.state.cortexStyle;
    if (this.ctx.state.cortexStyle !== 'solid') { this.ctx.state.cortexStyle = 'solid'; this.ctx.applyStyles(); }
    this.ctx.schedulePlate(0);
    if (!this.labelsEl) {
      this.labelsEl = document.createElement('div');
      this.labelsEl.id = 'electrode-labels';
      document.getElementById('viewport').appendChild(this.labelsEl);
      this.ctx.onViewHook(() => { if (this.active) this.placeLabels(); });
    }
    app.viewer.setScale(this.scales, CENTRE);
    const v = app.viewer;
    v.flyTo({ target: [0, -14, 8], direction: [-0.72, 0.52, 0.46], distance: v.fitDistance(122, 1.18), duration: 1000 });
    if (!this.mask) {
      this.ctx.toast('Reading the template’s head…', 8000);
      this.mask = headMask(await loadVolume(FILES.head));
      if (!this.active) return;
    }
    if (!this.mesh) {
      this.ctx.toast('Shaping the scalp…', 8000);
      this.mesh = await surface({ kind: 'mask', dims: this.mask.dims, affine: this.mask.affine, data: this.mask.data, step: 2, smoothing: 4 });
      if (!this.active) return;
    }
    this.showScalp();
    this.place();
    this.install();
  }

  leave() {
    this.active = false;
    const app = this.app;
    app.plateOff = false;
    if (this.savedCortex) { this.ctx.state.cortexStyle = this.savedCortex; this.ctx.applyStyles(); }
    app.viewer.setScale([1, 1, 1], CENTRE);
    if (this.scalp) app.viewer.setVisible(this.scalp, false);
    app.viewer.setMarkers([]);
    app.cortexColour = null; app.cardOverride = null;
    if (this.labelsEl) this.labelsEl.innerHTML = '';
    document.getElementById('legend').innerHTML = '';
    this.ctx.paintCortex();
    this.ctx.schedulePlate();
  }

  showScalp() {
    const app = this.app;
    if (!this.scalp) {
      this.scalp = app.viewer.addFreeMesh('scalp', { positions: this.mesh.positions, normals: this.mesh.normals, indices: this.mesh.indices }, { kind: 'skin', colour: '#efd3c4', pickable: false, render: 4 });
    }
    app.viewer.setVisible(this.scalp, this.scalpStyle !== 'hidden');
    app.viewer.setStyle(this.scalp, this.scalpStyle === 'glass' ? 'glass' : 'solid');
    this.paintScalp();
  }

  // The scalp cast through the head (scaled to the subject when measurements are set), the
  // landmarks put on it, and the electrodes placed by arcs.
  place() {
    const A = this.affine;
    const mask = A ? { ...this.mask, affine: compose(A, this.mask.affine) } : this.mask;
    const centre = A ? apply(A, CENTRE) : CENTRE;
    this.scalpModel = maskScalp(mask, { centre, axes: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], step: 2, maxRadius: 230 });
    const landmarks = {};
    for (const [k, p] of Object.entries(LANDMARKS)) landmarks[k] = this.scalpModel.onSurface(A ? apply(A, p) : p);
    this.landmarks = landmarks;
    const placed = tenTwenty(this.scalpModel, landmarks, { system: this.system === '10-20' ? '10-10' : '10-10', samples: 240 });
    const keep = this.system === '10-20' ? new Set(['Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'T7', 'C3', 'Cz', 'C4', 'T8', 'P7', 'P3', 'Pz', 'P4', 'P8', 'O1', 'O2', 'Fpz', 'Oz']) : null;
    this.positions = placed.positions.filter((p) => !LANDMARK_NAMES[p.name] && (!keep || keep.has(p.name)));
    this.czDiscrepancy = placed.czDiscrepancy;
    this.landmarkPoints = placed.positions.filter((p) => LANDMARK_NAMES[p.name]);
    if (this.channels) this.mapChannels(false);
  }

  toTemplate(p) { return this.affine ? apply(invert(this.affine), p) : p; }

  install() {
    const app = this.app;
    this.markers();
    this.placeLabels();
    app.cortexColour = this.spline ? (inst, i) => this.cortexValueColour(inst, i) : null;
    for (const inst of app.instances.get('cerebral-cortex') || []) inst.borderFlags = null;
    this.ctx.paintCortex();
    this.legend();
    this.ctx.renderLeft();
    this.ctx.renderToolbar();
    if (this.chosen) this.chooseElectrode(this.chosen);
  }

  markers() {
    const list = [];
    for (const e of this.positions) {
      let colour = '#f6dfb8';
      if (this.spline && this.channels.has(e.name)) colour = this.css(this.channels.get(e.name));
      list.push({ position: this.toTemplate(e.position), colour, radius: e.name === this.chosen ? 4.2 : 3.1 });
    }
    for (const l of this.landmarkPoints || []) list.push({ position: this.toTemplate(l.position), colour: '#e8c3dc', radius: 2.2 });
    this.app.viewer.setMarkers(list);
  }

  // Names beside the electrodes that face the viewer.
  placeLabels() {
    const el = this.labelsEl;
    if (!el) return;
    if (!this.names || !this.active) { el.innerHTML = ''; return; }
    const v = this.app.viewer;
    const cam = v.camera.position;
    const out = [];
    const centre = CENTRE;
    for (const e of this.positions) {
      const p = this.toTemplate(e.position);
      // Facing: the outward direction from the centre against the line of sight (in the world).
      const world = v.root.localToWorld(v.camera.position.clone().set(...p));
      const cw = v.root.localToWorld(v.camera.position.clone().set(...centre));
      const out1 = [world.x - cw.x, world.y - cw.y, world.z - cw.z];
      const toCam = [cam.x - world.x, cam.y - world.y, cam.z - world.z];
      if (out1[0] * toCam[0] + out1[1] * toCam[1] + out1[2] * toCam[2] < 0) continue;
      const q = v.project(p);
      if (!q) continue;
      const value = this.spline && this.channels.has(e.name) ? ` <small>${tidy(this.channels.get(e.name))}</small>` : '';
      out.push(`<div class="electrode-label${e.name === this.chosen ? ' chosen' : ''}" style="left:${q.x.toFixed(1)}px;top:${q.y.toFixed(1)}px">${esc(e.name)}${value}</div>`);
    }
    el.innerHTML = out.join('');
  }

  // ---- values at the electrodes ------------------------------------------------------------

  mapChannels(report = true) {
    const { values, unknown } = channelValues(this.channelsText, this.positions.map((p) => p.name));
    if (values.size < 3) {
      this.spline = null; this.channels = null;
      if (report) this.report = '<span class="bad">At least three channels of the system shown are needed.</span>';
      this.install();
      return;
    }
    const centre = this.scalpModel.centre;
    const dirs = [], vals = [];
    for (const [name, value] of values) { const e = this.positions.find((p) => p.name === name); dirs.push(sub(e.position, centre)); vals.push(value); }
    this.channels = values;
    this.spline = sphericalSpline(dirs, vals, { regularisation: 1e-9 });
    this.directions = dirs.map(unit);
    this.range = valueRange(vals, { diverging: vals.some((x) => x < 0) && vals.some((x) => x > 0) });
    this.colormap = vals.some((x) => x < 0) && vals.some((x) => x > 0) ? 'balance' : 'thermal';
    // g(t) tabulated, so that the map can be evaluated at every vertex quickly.
    const N = 2048;
    this.gTable = new Float64Array(N + 1);
    for (let i = 0; i <= N; i++) this.gTable[i] = this.spline.g(-1 + (2 * i) / N);
    this.bins = Array.from({ length: 256 }, (_, i) => colourAt(this.colormap, i / 255).map((c) => c / 255));
    if (report) this.report = `<span class="good">${values.size} channels mapped</span>${unknown.length ? `; <span class="bad">not in the ${esc(this.system)} system: ${esc(unknown.slice(0, 8).join(', '))}</span>` : ''}.`;
    this.paintScalp();
    if (report) this.install();
  }

  // The interpolated value in a direction from the centre (subject space).
  valueAt(d) {
    const u = unit(d), c = this.spline.coefficients, e = this.directions, g = this.gTable, N = g.length - 1;
    let s = this.spline.constant;
    for (let i = 0; i < e.length; i++) {
      const t = u[0] * e[i][0] + u[1] * e[i][1] + u[2] * e[i][2];
      const x = (Math.max(-1, Math.min(1, t)) + 1) * 0.5 * N, k = Math.min(N - 1, Math.floor(x)), f = x - k;
      s += c[i] * (g[k] + (g[k + 1] - g[k]) * f);
    }
    return s;
  }

  bin(v) { const [lo, hi] = this.range; return this.bins[Math.max(0, Math.min(255, Math.round((255 * (v - lo)) / (hi - lo || 1))))]; }
  css(v) { const c = this.bin(v); return `rgb(${c.map((x) => Math.round(255 * x)).join(',')})`; }

  cortexValueColour(inst, i) {
    const pos = inst.mesh.geometry.attributes.position.array;
    const p = [pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]];
    const q = this.affine ? apply(this.affine, p) : p;
    return this.bin(this.valueAt(sub(q, this.scalpModel.centre)));
  }

  paintScalp() {
    if (!this.scalp) return;
    const v = this.app.viewer;
    if (!this.spline) { v.paint(this.scalp); return; }
    const pos = this.scalp.mesh.geometry.attributes.position.array;
    v.paint(this.scalp, (i) => {
      const p = [pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]];
      const q = this.affine ? apply(this.affine, p) : p;
      return this.bin(this.valueAt(sub(q, this.scalpModel.centre)));
    });
  }

  legend() {
    const el = document.getElementById('legend');
    if (!this.spline) { el.innerHTML = ''; return; }
    const [lo, hi] = this.range;
    el.innerHTML = `<div><div class="title"><span>Over the scalp</span><span class="unit">${this.channels.size} channels</span></div>
<div class="bar" style="background:${cssGradient(this.colormap)}"></div>
<div class="ticks"><span style="left:0%">${tidy(lo)}</span><span style="left:50%">${tidy((lo + hi) / 2)}</span><span style="left:100%">${tidy(hi)}</span></div>
<div class="note">Spherical splines (Perrin et al. 1989), carried radially onto the cortex: a picture of the scalp map, not of its sources.</div></div>`;
  }

  // ---- the head's size ------------------------------------------------------------------------

  applySize(reset = false) {
    if (reset) this.measurements = { ...this.template };
    const s = scaleToHead(LANDMARKS, this.measurements, CENTRE);
    const identity = s.scales.every((x) => Math.abs(x - 1) < 1e-6);
    this.scales = s.scales;
    this.affine = identity ? null : s.affine;
    this.app.viewer.setScale(this.scales, CENTRE);
    this.place();
    this.paintScalp();
    this.install();
    this.ctx.toast(identity ? 'The template’s size' : `Scaled by ${this.scales.map((x) => x.toFixed(3)).join(' × ')}`);
  }

  // ---- choosing an electrode --------------------------------------------------------------------

  pick(event) {
    const i = this.app.viewer.pickMarker(event.clientX, event.clientY);
    if (i < 0 || i >= this.positions.length) return i >= 0;
    this.chooseElectrode(this.positions[i].name);
    return true;
  }

  hoverInfo(pointer) {
    const i = this.app.viewer.pickMarker(pointer.x, pointer.y);
    if (i < 0) return null;
    if (i >= this.positions.length) { const l = this.landmarkPoints[i - this.positions.length]; return { name: LANDMARK_NAMES[l.name], sub: 'A landmark of the head', colour: '#e8c3dc' }; }
    const e = this.positions[i];
    const under = this.beneath(e);
    return { name: e.name, sub: under ? `over the ${under.name.toLowerCase()}${under.side ? `, ${under.side}` : ''}` : 'an electrode', colour: '#f6dfb8' };
  }

  // The cortex beneath an electrode: going from it toward the centre of the head, the first
  // structure of the anatomy met, and the depth at which the brain begins.
  beneath(e) {
    const s = this.app.sections;
    if (!s) return null;
    const p0 = this.toTemplate(e.position);
    const d = unit(sub(CENTRE, p0));
    let brainAt = null;
    for (let r = 0; r < 90; r += 0.5) {
      const p = [p0[0] + d[0] * r, p0[1] + d[1] * r, p0[2] + d[2] * r];
      const key = this.app.labelValues.get(s.labelAt(p));
      if (!key) continue;
      const at = structureOfKey(key);
      if (!at) continue;
      if (at.id === 'csf' || at.id === 'vessel') continue;
      if (brainAt === null) brainAt = r;
      if (at.id === 'cerebral-white-matter' || at.id === 'cerebral-cortex') continue;
      const st = this.ctx.structure(at.id);
      return { id: at.id, side: at.side === 'mid' ? null : at.side, name: st.name, depth: brainAt, at: p };
    }
    return null;
  }

  chooseElectrode(name) {
    this.chosen = name;
    const e = this.positions.find((p) => p.name === name);
    if (!e) return;
    this.markers();
    this.placeLabels();
    const under = this.beneath(e);
    if (under) this.ctx.setPoint(under.at);
    const app = this.app;
    app.cardOverride = () => this.electrodeCard(e, under);
    this.ctx.renderCard();
  }

  electrodeCard(e, under) {
    const letters = e.name.replace(/[0-9z]+$/i, '');
    const region = REGION_LETTERS[letters] || '';
    const num = e.name.match(/(\d+)$/);
    const side = /z$/i.test(e.name) ? 'on the midline' : num ? (Number(num[1]) % 2 ? 'over the left hemisphere' : 'over the right hemisphere') : '';
    const older = TEN_TWENTY_NAMES_1985[e.name];
    const tpl = this.toTemplate(e.position);
    const facts = [['Position, MNI', `${coords(tpl)} <small>mm</small>`]];
    if (this.affine) facts.push(['On the subject', `${coords(e.position)} <small>mm</small>`]);
    if (under && under.depth !== null) facts.push(['Scalp to brain', `${Math.round(under.depth)} <small>mm</small>`]);
    if (this.spline && this.channels.has(e.name)) facts.push(['Value', tidy(this.channels.get(e.name))]);
    const chip = under ? `<button type="button" data-select="${esc(under.id)}" data-side="${esc(under.side || '')}"><i style="background:${esc(this.ctx.colourOf(under.id) || '#ccc')}"></i>${esc(under.name)}${under.side ? ` <small>${esc(under.side)}</small>` : ''}</button>` : '';
    return `<div class="card-head"><div class="over">An electrode of the ${esc(this.system)} system</div>
<h2><span class="sw" style="background:#f6dfb8"></span><span>${esc(e.name)}</span></h2>
${older ? `<div class="latin">formerly ${esc(older)}</div>` : ''}
<button type="button" class="close" data-action="electrode-close" aria-label="Close"><svg><use href="#i-close"/></svg></button></div>
<div class="card-body"><p class="summary">The ${esc(region)} position ${esc(side)}${e.arc ? `, placed on the ${esc(e.arc)} ${e.arc === 'ring' ? 'through Fpz, T7, Oz and T8' : 'arc'} at ${Math.round(100 * e.fraction)} % of its length` : ''}.</p>
<div class="facts">${facts.map(([k, v]) => `<div><div class="k">${esc(k)}</div><div class="v">${v}</div></div>`).join('')}</div>
${chip ? `<div class="sect"><h4>The cortex beneath</h4><p>Going straight in, toward the centre of the head, the first part of the brain met is:</p><div class="chips">${chip}</div><p class="note">What an electrode records comes from a wide patch of cortex, mostly from gyri whose surfaces face it; this is where it sits, not all it hears.</p></div>` : ''}
<div class="sect"><h4>How the names are made</h4><p>The letters name the region beneath (Fp frontopolar, F frontal, C central, T temporal, P parietal, O occipital, and their combinations); <i>z</i> marks the midline; odd numbers lie on the left, even on the right, rising away from the midline. The positions divide the arcs from the nasion to the inion and between the ears into tenths and fifths, so that they fall on every head at the same proportional places.</p></div>
<div class="card-refs"><p>Jasper, H. H. (1958), The ten-twenty electrode system of the International Federation, <i>Electroencephalography and Clinical Neurophysiology</i> 10, 371–375.</p><p>Oostenveld, R. and Praamstra, P. (2001), The five percent electrode system for high-resolution EEG and ERP measurements, <i>Clinical Neurophysiology</i> 112, 713–719.</p></div></div>`;
  }

  cardAction(a) {
    if (a === 'electrode-close') { this.chosen = null; this.app.cardOverride = null; this.markers(); this.placeLabels(); this.ctx.renderCard(); }
  }

  // ---- the toolbar and the panel ------------------------------------------------------------------

  toolbar() {
    const opt = (action, label, pressed) => `<button type="button" class="opt" data-mode-action="${action}" aria-pressed="${pressed}">${label}</button>`;
    return `<div class="group"><span class="caption">Scalp</span><div class="row">${opt('scalp:glass', 'Glass', this.scalpStyle === 'glass')}${opt('scalp:skin', 'Skin', this.scalpStyle === 'skin')}${opt('scalp:hidden', 'Hidden', this.scalpStyle === 'hidden')}</div></div>
<div class="group"><span class="caption">Electrodes</span><div class="row">${opt('system:10-20', '10–20', this.system === '10-20')}${opt('system:10-10', '10–10', this.system === '10-10')}</div></div>
<div class="group"><span class="caption">Names</span><div class="row">${opt('names', 'Names', this.names)}</div></div>`;
  }

  toolbarAction(action) {
    const [what, value] = action.split(':');
    if (what === 'scalp') { this.scalpStyle = value; this.showScalp(); }
    if (what === 'system') { this.system = value; this.chosen = null; this.place(); this.install(); }
    if (what === 'names') { this.names = !this.names; this.placeLabels(); }
    this.ctx.renderToolbar();
  }

  panel() {
    const field = (key, label) => `<label class="field"><span>${label} <span class="was">template ${Math.round(this.template[key])} mm</span></span><input type="number" step="1" min="60" max="400" data-measure="${key}" value="${Math.round(this.measurements[key])}"></label>`;
    return `<h3>The head</h3>
<p class="pad note">The template’s scalp, with the electrodes of the ${esc(this.system)} system placed on it as a technician places them: by tenths and fifths of the arcs from the nasion over the vertex to the inion, and from ear to ear.</p>
<h3>A subject’s head</h3>
<div class="pad">
${field('length', 'Length, nasion to inion')}
${field('breadth', 'Breadth, ear to ear')}
${field('height', 'Height, ears to vertex')}
<div class="row-btns"><button type="button" class="btn primary" data-head="scale">Scale to this head</button><button type="button" class="btn" data-head="reset">Template size</button></div>
<p class="note">Straight distances, not tape measures. The brain and the scalp are scaled about the centre of the head; the electrodes are placed again on the scaled scalp. Scales now ${this.scales.map((s) => s.toFixed(3)).join(' × ')}.</p>
</div>
<h3>Values at the electrodes</h3>
<p class="pad note">One channel a line, as <b>Cz, 12.5</b>: a power, an amplitude, a potential. The values are interpolated over the scalp and shown on the brain beneath.</p>
<div class="pad data-form">
<textarea id="head-channels" spellcheck="false" placeholder="Fz, 2.1&#10;Cz, 4.0&#10;Pz, 1.2&#10;C3, -1.0&#10;C4, 3.5">${esc(this.channelsText)}</textarea>
<div class="row-btns"><button type="button" class="btn primary" data-head="map">Map the channels</button><button type="button" class="btn" data-head="clear">Clear</button><button type="button" class="btn" data-head="example">An example</button></div>
<div class="report">${this.report}</div>
</div>
${this.positions.length ? `<p class="pad note">${this.positions.length} electrodes; the coronal and midline arcs meet within ${this.czDiscrepancy.toFixed(1)} mm at Cz. The landmarks of the template are provisional.</p>` : ''}`;
  }

  bind(left) {
    const t = left.querySelector('#head-channels');
    t.addEventListener('input', () => { this.channelsText = t.value; });
    left.querySelectorAll('[data-measure]').forEach((inp) => inp.addEventListener('change', () => { const v = Number(inp.value); if (v > 0) this.measurements[inp.dataset.measure] = v; }));
    left.addEventListener('click', (e) => {
      const b = e.target.closest('[data-head]');
      if (!b) return;
      const a = b.dataset.head;
      if (a === 'scale') this.applySize(false);
      if (a === 'reset') this.applySize(true);
      if (a === 'map') this.mapChannels();
      if (a === 'clear') { this.channelsText = ''; this.spline = null; this.channels = null; this.report = ''; this.paintScalp(); this.install(); }
      if (a === 'example') { this.channelsText = EXAMPLE; this.mapChannels(); }
    });
  }
}

// The head as a mask: the T1 image above a tenth of its robust range, and everything that
// cannot be reached from the edges of the volume without crossing the head (the sinuses,
// the ear canals, the dark bone) filled in, so that the scalp is one closed surface.
export function headMask(t1) {
  const [nx, ny, nz] = t1.dims;
  const [lo, hi] = robustRange(t1.data, 0.02, 0.98);
  const thr = lo + 0.12 * (hi - lo);
  const n = nx * ny * nz;
  const data = new Uint8Array(n);
  for (let i = 0; i < n; i++) data[i] = t1.data[i] > thr ? 1 : 0;
  // Flood the outside from the faces of the volume (6-connected), but not from the lowest
  // slice, which cuts through the neck and the mouth: the airways there are inside the head.
  const outside = new Uint8Array(n);
  const queue = new Int32Array(n);
  let head = 0, tail = 0;
  const push = (i) => { if (!outside[i] && !data[i]) { outside[i] = 1; queue[tail++] = i; } };
  for (let k = 1; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    if (i === 0 || j === 0 || i === nx - 1 || j === ny - 1 || k === nz - 1) push(i + nx * (j + ny * k));
  }
  const sxy = nx * ny;
  while (head < tail) {
    const v = queue[head++];
    const i = v % nx, j = Math.floor(v / nx) % ny, k = Math.floor(v / sxy);
    if (i > 0) push(v - 1); if (i < nx - 1) push(v + 1);
    if (j > 0) push(v - nx); if (j < ny - 1) push(v + nx);
    if (k > 1) push(v - sxy); if (k < nz - 1) push(v + sxy);
  }
  for (let i = 0; i < sxy; i++) if (outside[i + sxy]) outside[i] = 1;
  for (let i = 0; i < n; i++) data[i] = outside[i] ? 0 : 1;
  return { dims: t1.dims, affine: t1.affine, data };
}

// The product of two 4 × 4 affines given as 16 numbers in rows.
function compose(A, B) {
  const m = new Array(16).fill(0);
  for (let r = 0; r < 4; r++) for (let c = 0; c < 4; c++) { let s = 0; for (let k = 0; k < 4; k++) s += A[4 * r + k] * B[4 * k + c]; m[4 * r + c] = s; }
  return m;
}
