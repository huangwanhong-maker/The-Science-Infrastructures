// The mouse brain in the Atlases mode: the Allen Mouse Brain Common Coordinate Framework
// (CCFv3), its major divisions drawn in the round in the atlas's own colours (softened),
// its ontology as a tree, its template and annotation in the sections, and a card for
// every structure. The human brain is put away while the mouse is shown.
//
// The mouse brain is some 10 mm long; it is drawn twelve times enlarged, about its centre,
// so that it fills the stage as the human brain does. Coordinates stay in millimetres.

import { loadMouse, indicesOf, measureRegion } from '../data/mouse.js';
import { soften } from '../data/parcellations.js';
import { Sections } from '../sections/sections.js';
import { esc, volumeText } from '../ui/card.js';
import { rasToCcf, CITATIONS } from '../../../../domains/neuroscience/js/index.js';

const ENLARGE = 12;
// The major divisions drawn in the round, by acronym (Allen ontology).
const DIVISIONS = ['Isocortex', 'OLF', 'HPF', 'CTXsp', 'STR', 'PAL', 'TH', 'HY', 'MB', 'P', 'MY', 'CB'];
const signed = (x, d = 2) => (x < 0 ? `−${(-x).toFixed(d)}` : x.toFixed(d));

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

const css = (rgb) => `rgb(${rgb.map((x) => Math.round(255 * x)).join(',')})`;
const hex = (rgb) => `#${rgb.map((x) => Math.round(255 * x).toString(16).padStart(2, '0')).join('')}`;

export class MouseView {
  constructor(ctx) {
    this.ctx = ctx;
    this.mouse = null;
    this.selected = null;          // a Region of the ontology
    this.expanded = new Set();
    this.filter = '';
    this.style = 'glass';          // the whole brain's surface: glass, or hidden
    this.instances = [];           // the free meshes added
  }

  get app() { return this.ctx.app; }

  async enter() {
    this.active = true;
    const app = this.app;
    app.otherSpecies = true;
    app.plateOff = true;
    this.ctx.applyStyles();
    this.ctx.schedulePlate(0);
    if (!this.mouse) {
      this.ctx.toast('Reading the Allen mouse brain…', 9000);
      this.mouse = await loadMouse();
      const root = this.mouse.ontology.get('Basic cell groups and regions');
      this.expanded = new Set([this.mouse.ontology.roots[0], root, this.mouse.ontology.get('Cerebrum'), this.mouse.ontology.get('Brain stem')].filter(Boolean));
      this.sections = new Sections({ t1: this.mouse.t1, labels: this.mouse.labels, palette: this.palette() });
      this.sections.overlay = 0.6;
      const [nx, ny, nz] = this.mouse.t1.dims, A = this.mouse.t1.affine;
      this.centre = [A[3] + (nx / 2) * A[0], A[7] + (ny / 2) * A[5], A[11] + (nz / 2) * A[10]];
    }
    if (!this.active) return;
    // The sections show the mouse.
    this.humanSections = app.sections;
    this.humanPoint = this.ctx.state.point;
    app.sections = this.sections;
    for (const v of Object.values(app.views)) v.sections = this.sections;
    this.ctx.setPoint(this.centre);
    // The round view: the mouse enlarged about its centre.
    app.viewer.setScale([ENLARGE, ENLARGE, ENLARGE], this.centre);
    app.viewer.flyTo({ target: this.centre, direction: [-0.62, 0.5, 0.6], distance: app.viewer.fitDistance(8 * ENLARGE, 1.35), duration: 900 });
    if (!this.built) await this.build();
    for (const inst of this.instances) app.viewer.setVisible(inst, true);
    this.showStyles();
    this.legend();
  }

  leave() {
    this.active = false;
    const app = this.app;
    app.otherSpecies = false;
    app.plateOff = false;
    for (const inst of this.instances) app.viewer.setVisible(inst, false);
    app.viewer.setScale([1, 1, 1], [0, 0, 0]);
    if (this.humanSections) {
      app.sections = this.humanSections;
      for (const v of Object.values(app.views)) v.sections = this.humanSections;
      this.ctx.setPoint(this.humanPoint || [0, -18, 10]);
    }
    app.cardOverride = null;
    document.getElementById('legend').innerHTML = '';
    this.ctx.applyStyles();
    const v = app.viewer;
    v.flyTo({ target: [0, -16, 12], direction: [-0.78, 0.38, 0.5], distance: v.fitDistance(98, 1.5), duration: 900 });
  }

  palette() {
    const pal = new Uint8Array(4 * this.mouse.regions.length);
    this.mouse.regions.forEach((r, k) => { if (k) { const c = soften(this.mouse.colours[k], 0.02); pal.set([Math.round(255 * c[0]), Math.round(255 * c[1]), Math.round(255 * c[2]), 255], 4 * k); } });
    return pal;
  }

  colourOf(region) {
    for (let r = region; r; r = r.parent) if (r.colour) return soften(r.colour.map((c) => c / 255), 0.02);
    return [0.82, 0.76, 0.8];
  }

  // The whole brain as glass, and its major divisions as solid parts.
  async build() {
    const app = this.app, m = this.mouse;
    this.built = true;
    this.ctx.toast('Shaping the mouse brain…', 12000);
    const dims = m.labels.dims, affine = m.labels.affine;
    const whole = await surface({ kind: 'mask', dims, affine, data: m.labels.data, step: 2, smoothing: 3 });
    this.brain = app.viewer.addFreeMesh('mouse:brain', whole, { kind: 'tissue', colour: '#ead2cc', pickable: false, render: 3 });
    this.instances.push(this.brain);
    this.parts = new Map();
    for (const acr of DIVISIONS) {
      const region = m.ontology.get(acr);
      if (!region) continue;
      const indices = indicesOf(m, region);
      if (!indices.length) continue;
      const mesh = await surface({ kind: 'region', dims, affine, data: m.labels.data, indices, step: 2, smoothing: 3 });
      if (!this.active) return;
      const inst = app.viewer.addFreeMesh(`mouse:${acr}`, mesh, { kind: 'tissue', colour: hex(this.colourOf(region)) });
      inst.mouseRegion = region;
      this.parts.set(region, inst);
      this.instances.push(inst);
    }
    this.ctx.toast('The mouse brain: Allen CCFv3, 2017 annotation');
  }

  showStyles() {
    const v = this.app.viewer;
    if (this.brain) { v.setVisible(this.brain, this.style !== 'hidden'); if (this.style !== 'hidden') v.setStyle(this.brain, 'glass'); }
    const sel = this.selected;
    for (const [region, inst] of this.parts || []) {
      const inSel = sel && (region === sel || region.isAncestorOf(sel) || sel.isAncestorOf(region));
      v.setStyle(inst, sel && !inSel ? 'ghost' : 'solid');
      v.setEmphasis(inst, sel && (region === sel) ? 0.18 : 0, '#f6d39a');
    }
    if (this.regionInst) v.setStyle(this.regionInst, 'solid');
  }

  // A structure chosen: its surface, its card, and its voxels lit in the sections.
  async select(region) {
    const app = this.app, m = this.mouse, v = app.viewer;
    this.selected = region;
    if (this.regionInst) { v.removeInstance(this.regionInst); this.instances = this.instances.filter((i) => i !== this.regionInst); this.regionInst = null; }
    this.showStyles();
    const meas = region ? measureRegion(m, region) : null;
    this.sections.selected = meas && meas.indices.length ? new Set(meas.indices) : null;
    this.sections.dimOthers = Boolean(meas);
    if (meas && meas.centroid) this.ctx.setPoint(meas.centroid);
    for (const view of Object.values(app.views)) view.redraw();
    app.cardOverride = () => { setTimeout(() => this.bindCard(), 0); return this.selected ? this.cardHtml() : this.introHtml(); };
    this.ctx.renderCard();
    const tree = document.querySelector('.mouse-tree');
    if (tree) tree.innerHTML = this.treeHtml();
    if (!region || !meas.indices.length || this.parts.has(region)) return;
    const small = meas.voxels < 20000;
    const mesh = await surface({ kind: 'region', dims: m.labels.dims, affine: m.labels.affine, data: m.labels.data, indices: meas.indices, step: small ? 1 : 2, smoothing: small ? 2 : 3 });
    if (this.selected !== region || !this.active) return;
    this.regionInst = v.addFreeMesh('mouse:region', mesh, { kind: 'tissue', colour: hex(this.colourOf(region)), render: 1 });
    v.setEmphasis(this.regionInst, 0.16, '#f6d39a');
    this.instances.push(this.regionInst);
    this.showStyles();
  }

  // ---- picking, sections ----------------------------------------------------------------------

  pick(event) {
    const hit = this.app.viewer.pick(event.clientX, event.clientY);
    const region = hit && hit.inst.mouseRegion;
    if (!region) { if (this.selected) this.select(null); return true; }
    this.select(region);
    return true;
  }
  hoverInfo(pointer) {
    const hit = this.app.viewer.pick(pointer.x, pointer.y);
    const region = hit && hit.inst.mouseRegion;
    if (!region) return null;
    return { name: region.name, sub: `${region.acronym} · Allen Mouse Brain Atlas`, colour: css(this.colourOf(region)) };
  }
  sectionAt(w) {
    const k = this.sections.labelAt(w);
    const r = k && this.mouse.regions[k];
    if (!r) return null;
    const path = r.ancestors().reverse().slice(2, 4).map((a) => a.acronym).join(' › ');
    return { name: `${r.name} (${r.acronym})`, sub: path, colour: css(this.colourOf(r)), choose: () => this.select(r) };
  }
  sectionSelection() { return this.sections.selected; }
  sectionPalette() { return this.palette(); }

  // ---- the card -------------------------------------------------------------------------------

  introHtml() {
    return `<div class="card-empty"><h2>The mouse brain</h2>
<p>The Allen Institute’s Common Coordinate Framework: the average of 1,675 young adult mice imaged by serial two-photon tomography, every voxel labelled by a structure of the Allen ontology (CCFv3, 2017 annotation), drawn here at 50 µm.</p>
<p>The mouse brain is about a thousandth of the human brain’s volume, yet built on the same plan: a cortex, though smooth and without gyri, over the striatum and the thalamus, a hippocampus, a cerebellum and a brainstem. Its olfactory bulbs, by comparison, are enormous.</p>
<p class="hint">Click a part, a section or an entry of the tree to read it. The brain is drawn twelve times enlarged.</p></div>`;
  }

  cardHtml() {
    const r = this.selected, m = this.mouse;
    const meas = measureRegion(m, r);
    const path = r.ancestors().reverse().slice(1);
    const facts = [];
    if (meas.voxels) {
      facts.push(['Volume', volumeText(meas.volume)]);
      facts.push(['Centre, RAS', `${meas.centroid.map((x) => signed(x)).join(', ')} <small>mm</small>`]);
      const ccf = rasToCcf(meas.centroid);
      facts.push(['Centre, CCF', `${ccf.map((x) => Math.round(x)).join(', ')} <small>µm</small>`]);
    }
    const kids = r.children.filter((c) => measureRegion(m, c).voxels).slice(0, 18);
    const chip = (x) => `<button type="button" data-mouse="${x.id}"><i style="background:${css(this.colourOf(x))}"></i>${esc(x.name)} <small>${esc(x.acronym)}</small></button>`;
    const share = meas.voxels && r.parent ? (() => { const pm = measureRegion(m, r.parent); return pm.voxels ? Math.round((100 * meas.voxels) / pm.voxels) : null; })() : null;
    return `<div class="card-head"><div class="over">${path.slice(-3).map((a) => `<button type="button" data-mouse="${a.id}">${esc(a.acronym || a.name)}</button>`).join('<span class="sep">›</span>')}</div>
<h2><span class="sw" style="background:${css(this.colourOf(r))}"></span><span>${esc(r.name)}</span></h2>
<div class="latin">${esc(r.acronym)} · Allen Mouse Brain Atlas, structure ${r.id}</div>
<button type="button" class="close" data-mouse-close aria-label="Close"><svg><use href="#i-close"/></svg></button></div>
<div class="card-body">
${facts.length ? `<div class="facts">${facts.map(([k, v]) => `<div><div class="k">${esc(k)}</div><div class="v">${v}</div></div>`).join('')}</div>` : '<p class="note">This structure has no voxels of its own in the 50 µm annotation: it is too small, or it is labelled only through its parts.</p>'}
${share !== null ? `<p class="note">${share} % of the ${esc(r.parent.name.toLowerCase())}.</p>` : ''}
${kids.length ? `<div class="sect"><h4>Contains</h4><div class="chips">${kids.map(chip).join('')}</div></div>` : ''}
<div class="sect"><h4>In the sections</h4><div class="mini-sections"><div><canvas data-mini="axial"></canvas><span>axial</span></div><div><canvas data-mini="coronal"></canvas><span>coronal</span></div><div><canvas data-mini="sagittal"></canvas><span>sagittal</span></div></div></div>
<div class="sect"><h4>The coordinates</h4><p>The framework’s own axes run from the anterior, superior, left corner of the volume: posterior, inferior and to the right, in micrometres. Here they are also given as right, anterior and superior millimetres from the same corner. Stereotaxic coordinates from bregma are a separate fit.</p></div>
<div class="card-refs"><p>${esc(CITATIONS.Wang2020 || 'Wang, Q. et al. (2020), The Allen Mouse Brain Common Coordinate Framework: a 3D reference atlas, Cell 181, 936–953.')}</p><p>Data: Allen Institute for Brain Science, under its terms of use.</p></div></div>`;
  }

  bindCard() {
    const el = document.getElementById('card');
    el.querySelectorAll('[data-mouse]').forEach((b) => b.addEventListener('click', (e) => { e.stopPropagation(); this.select(this.mouse.ontology.byId.get(Number(b.dataset.mouse))); }));
    const close = el.querySelector('[data-mouse-close]');
    if (close) close.addEventListener('click', (e) => { e.stopPropagation(); this.select(null); });
    const meas = this.selected && measureRegion(this.mouse, this.selected);
    if (meas && meas.centroid) {
      const axes = { axial: 2, coronal: 1, sagittal: 0 };
      for (const c of el.querySelectorAll('canvas[data-mini]')) this.sections.draw(c, c.dataset.mini, meas.centroid[axes[c.dataset.mini]], { scale: 2, flipU: c.dataset.mini === 'sagittal' });
    }
  }

  legend() {
    const keys = DIVISIONS.map((a) => this.mouse.ontology.get(a)).filter(Boolean).map((r) => `<span><i style="background:${css(this.colourOf(r))}"></i>${esc(r.acronym)}</span>`).join('');
    document.getElementById('legend').innerHTML = `<div><div class="title"><span>Major divisions</span><span class="unit">Allen ontology</span></div><div class="keys">${keys}</div><div class="note">Colours of the Allen atlas, softened.</div></div>`;
  }

  // ---- the panel and the toolbar ---------------------------------------------------------------

  panelHtml() {
    if (!this.mouse) return '<p class="pad note">Reading the Allen mouse brain…</p>';
    return `<h3>The Allen ontology</h3>
<div class="pad"><input type="search" class="region-filter mouse-filter" placeholder="Filter by name or acronym…" value="${esc(this.filter)}" autocomplete="off" spellcheck="false"></div>
<div class="tree mouse-tree">${this.treeHtml()}</div>`;
  }

  treeHtml() {
    const m = this.mouse;
    const q = this.filter.trim().toLowerCase();
    const has = (r) => measureRegion(m, r).voxels > 0;
    const row = (r, depth, twist) => `<div class="row${r === this.selected ? ' selected' : ''}" data-mouse="${r.id}" style="padding-left:${4 + depth * 13}px"><span class="twist ${twist ? 'has' : 'leaf'}" data-twist="${r.id}"><svg><use href="#i-chev"/></svg></span><span class="sw" style="background:${css(this.colourOf(r))}"></span><span class="name">${esc(r.name)}</span><span class="count">${esc(r.acronym)}</span></div>`;
    if (q) {
      const hits = m.ontology.all.filter((r) => (r.name.toLowerCase().includes(q) || (r.acronym || '').toLowerCase() === q || (r.acronym || '').toLowerCase().startsWith(q)) && has(r)).slice(0, 160);
      return hits.map((r) => row(r, 0, false)).join('') || '<p class="note" style="padding:4px 8px">Nothing by that name.</p>';
    }
    const out = [];
    // The selected structure's ancestors are opened, so that it can be seen in the tree.
    if (this.selected) for (const a of this.selected.ancestors()) this.expanded.add(a);
    const walk = (r, depth) => {
      if (!has(r)) return;
      const kids = r.children.filter(has);
      out.push(row(r, depth, kids.length).replace('class="twist has"', `class="twist has${this.expanded.has(r) ? ' open' : ''}"`));
      if (this.expanded.has(r)) for (const c of kids) walk(c, depth + 1);
    };
    for (const root of m.ontology.roots) walk(root, 0);
    return out.join('');
  }

  bindPanel(left) {
    const filter = left.querySelector('.mouse-filter');
    const tree = left.querySelector('.mouse-tree');
    if (filter) filter.addEventListener('input', () => { this.filter = filter.value; tree.innerHTML = this.treeHtml(); });
    if (tree) tree.addEventListener('click', (e) => {
      const tw = e.target.closest('[data-twist]');
      const id = Number((tw || e.target.closest('[data-mouse]'))?.dataset.twist || e.target.closest('[data-mouse]')?.dataset.mouse);
      const r = this.mouse.ontology.byId.get(id);
      if (!r) return;
      if (tw && tw.classList.contains('has')) { if (this.expanded.has(r)) this.expanded.delete(r); else this.expanded.add(r); tree.innerHTML = this.treeHtml(); return; }
      this.select(r);
    });
  }

  toolbar(state) {
    const opt = (attr, value, label, pressed) => `<button type="button" class="opt" ${attr}="${value}" aria-pressed="${pressed}">${label}</button>`;
    const g = [`<div class="group"><span class="caption">Show</span><div class="row">${opt('data-show', '3d', 'In the round', state.show !== 'sections')}${opt('data-show', 'sections', 'In sections', state.show === 'sections')}</div></div>`];
    if (state.show === 'sections') g.push(`<div class="group"><span class="caption">Main section</span><div class="row">${['axial', 'coronal', 'sagittal'].map((p) => opt('data-main', p, p[0].toUpperCase() + p.slice(1), state.mainPlane === p)).join('')}</div></div>`);
    else g.push(`<div class="group"><span class="caption">Outline</span><div class="row">${opt('data-mode-action', 'mouse-outline:glass', 'Glass', this.style === 'glass')}${opt('data-mode-action', 'mouse-outline:hidden', 'Hidden', this.style === 'hidden')}</div></div>`);
    return g.join('');
  }

  toolbarAction(action) {
    const [what, value] = action.split(':');
    if (what === 'mouse-outline') { this.style = value; this.showStyles(); }
    this.ctx.renderToolbar();
  }
}
