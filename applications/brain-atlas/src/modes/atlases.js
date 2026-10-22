// The Atlases mode: the cortex coloured by another division of it, Julich-Brain's
// cytoarchitectonic areas or Schaefer's parcels of the resting-state networks, in
// three dimensions and in the sections; each area or parcel with a card that says
// where it lies in the gross anatomy.
//
// An area is assigned to each vertex of the cortex by the label of the voxel at the
// vertex or, the surface lying at the outer edge of the grey matter, the first label
// met going inward along the normal (up to 4 mm); then one pass of a majority filter
// over the mesh's neighbours removes stray vertices, and unlabelled vertices take
// their neighbours' label where they have one.

import { availableParcellations, loadParcellation, NETWORKS, soften, groupOf } from '../data/parcellations.js';
import { esc, volumeText } from '../ui/card.js';
import { structureOfKey, invert, apply, CITATIONS, mniToTalairach } from '../../../../domains/neuroscience/js/index.js';
import { MouseView } from './mouse.js';
import { MOUSE_FILES } from '../data/mouse.js';
import { toolbarHtml } from '../ui/controls.js';

const NEUTRAL = [0.83, 0.78, 0.8];
const STEPS = [0, -1, -2, -3, -4, 1];
const signed = (x) => (x < 0 ? `−${-x}` : `${x}`);
const coords = (c) => c.map((x) => signed(Math.round(x))).join(', ');

export class AtlasesMode {
  constructor(ctx) {
    this.ctx = ctx;
    this.available = null;
    this.loaded = new Map();
    this.current = null;        // the loaded parcellation shown
    this.colourBy = 'parcel';   // or 'network' (Schaefer)
    this.selected = 0;          // the index of the chosen region
    this.filter = '';
    this.busy = false;
    this.species = 'human';     // or 'mouse'
    this.mouseView = new MouseView(ctx);
    this.mouseAvailable = false;
  }

  get app() { return this.ctx.app; }

  async enter() {
    this.active = true;
    // The divisions are drawn on the cortex, which is shown solid while the mode lasts.
    this.savedCortex = this.ctx.state.cortexStyle;
    if (this.ctx.state.cortexStyle !== 'solid') { this.ctx.state.cortexStyle = 'solid'; this.ctx.applyStyles(); this.ctx.renderToolbar(); }
    this.ctx.state.selected = null;
    this.app.cardOverride = () => this.introHtml();
    this.ctx.renderCard();
    if (!this.available) {
      this.available = await availableParcellations();
      try { this.mouseAvailable = (await fetch(`data/${MOUSE_FILES.annotation}`, { method: 'HEAD' })).ok; } catch { this.mouseAvailable = false; }
      this.ctx.renderLeft();
    }
    if (!this.active) return;
    if (this.species === 'mouse') { await this.showMouse(); return; }
    if (!this.current && this.available.length) await this.choose(this.available[0].id);
    else if (this.current) this.install();
  }

  leave() {
    this.active = false;
    if (this.species === 'mouse') this.mouseView.leave();
    this.leaveHuman();
  }

  // The human brain's divisions put away (on leaving the mode, or turning to the mouse).
  leaveHuman() {
    const app = this.app;
    if (this.savedCortex && this.savedCortex !== this.ctx.state.cortexStyle) { this.ctx.state.cortexStyle = this.savedCortex; this.ctx.applyStyles(); }
    app.cortexColour = null; app.cortexKey = null; app.cortexRegion = null; app.regionLit = null; app.cardOverride = null;
    for (const inst of app.instances.get('cerebral-cortex') || []) inst.borderFlags = null;
    if (this.saved && app.sections) { app.sections.labels = this.saved.labels; app.sections.palette = this.saved.palette; app.sections.selected = null; }
    this.saved = null;
    document.getElementById('legend').innerHTML = '';
    this.ctx.paintCortex();
  }

  async setSpecies(species) {
    if (species === this.species) return;
    if (species === 'mouse') {
      this.leaveHuman();
      this.species = 'mouse';
      this.ctx.renderLeft();
      await this.showMouse();
    } else {
      this.mouseView.leave();
      this.species = 'human';
      if (this.ctx.state.cortexStyle !== 'solid') { this.ctx.state.cortexStyle = 'solid'; this.ctx.applyStyles(); }
      this.app.cardOverride = () => this.introHtml();
      this.ctx.renderToolbar();
      if (this.current) this.install(); else if (this.available?.length) await this.choose(this.available[0].id);
    }
  }

  async showMouse() {
    await this.mouseView.enter();
    const view = this.mouseView;
    this.app.cardOverride = () => { setTimeout(() => view.bindCard(), 0); return view.selected ? view.cardHtml() : view.introHtml(); };
    this.ctx.renderCard();
    this.ctx.renderToolbar();
    this.ctx.renderLeft();
  }

  // Hooks of the app, given to the mouse when it is shown.
  pick(event) { return this.species === 'mouse' && this.active ? this.mouseView.pick(event) : false; }
  hoverInfo(pointer) { return this.species === 'mouse' && this.active ? this.mouseView.hoverInfo(pointer) : null; }
  toolbar() { return this.species === 'mouse' ? this.mouseView.toolbar(this.ctx.state) : toolbarHtml(this.ctx.state); }
  toolbarAction(action, b) { if (this.species === 'mouse') this.mouseView.toolbarAction(action, b); }

  async choose(id) {
    const spec = this.available.find((p) => p.id === id);
    if (!spec) return;
    if (!this.loaded.has(id)) {
      this.busy = true;
      this.ctx.renderLeft();
      this.ctx.toast(`Reading ${spec.name}…`, 6000);
      const app = this.app;
      while (!app.sections) await new Promise((r) => setTimeout(r, 200));
      const grid = app.sections.labels.dims ? (this.saved ? this.saved.labels : app.sections.labels) : null;
      const loaded = await loadParcellation(spec, grid);
      this.measure(loaded, grid);
      this.loaded.set(id, loaded);
      this.busy = false;
      this.ctx.toast(`${spec.name}: ${loaded.regions.length - 1} regions`);
    }
    this.current = this.loaded.get(id);
    this.selected = 0;
    this.colourBy = this.current.spec.reader === 'schaefer' ? this.colourBy : 'parcel';
    if (this.active) this.install();
  }

  // Each region's voxels, centre, bounds, and the anatomical structures it overlaps.
  measure(loaded, grid) {
    const n = loaded.regions.length;
    const [nx, ny, nz] = grid.dims;
    const A = grid.affine;
    const sum = new Float64Array(3 * n), lo = new Float32Array(3 * n).fill(Infinity), hi = new Float32Array(3 * n).fill(-Infinity);
    const overlap = new Map();
    const idx = loaded.index, anat = grid.data;
    let v = 0;
    for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++, v++) {
      const r = idx[v];
      if (!r) continue;
      const x = A[0] * i + A[1] * j + A[2] * k + A[3], y = A[4] * i + A[5] * j + A[6] * k + A[7], z = A[8] * i + A[9] * j + A[10] * k + A[11];
      sum[3 * r] += x; sum[3 * r + 1] += y; sum[3 * r + 2] += z;
      if (x < lo[3 * r]) lo[3 * r] = x; if (x > hi[3 * r]) hi[3 * r] = x;
      if (y < lo[3 * r + 1]) lo[3 * r + 1] = y; if (y > hi[3 * r + 1]) hi[3 * r + 1] = y;
      if (z < lo[3 * r + 2]) lo[3 * r + 2] = z; if (z > hi[3 * r + 2]) hi[3 * r + 2] = z;
      const a = anat[v];
      if (a) { const key = r * 65536 + a; overlap.set(key, (overlap.get(key) || 0) + 1); }
    }
    const byRegion = Array.from({ length: n }, () => []);
    for (const [key, count] of overlap) byRegion[Math.floor(key / 65536)].push([key % 65536, count]);
    loaded.stats = Array.from({ length: n }, (_, r) => {
      const c = loaded.voxels[r];
      return c ? {
        voxels: c,
        centroid: [sum[3 * r] / c, sum[3 * r + 1] / c, sum[3 * r + 2] / c],
        bounds: { min: [lo[3 * r], lo[3 * r + 1], lo[3 * r + 2]], max: [hi[3 * r], hi[3 * r + 1], hi[3 * r + 2]] },
        overlap: byRegion[r].sort((p, q) => q[1] - p[1]),
      } : null;
    });
  }

  // The parcellation's region at each vertex of a cortex instance.
  vertexRegions(inst) {
    const L = this.current;
    inst.parcels = inst.parcels || new Map();
    if (inst.parcels.has(L.spec.id)) return inst.parcels.get(L.spec.id);
    const app = this.app;
    const grid = this.saved ? this.saved.labels : app.sections.labels;
    const [nx, ny, nz] = grid.dims;
    const inv = invert(grid.affine);
    const pos = inst.mesh.geometry.attributes.position.array, nrm = inst.mesh.geometry.attributes.normal.array;
    const nv = pos.length / 3;
    const side = inst.side;
    const ok = new Uint8Array(L.regions.length);
    for (let k = 1; k < L.regions.length; k++) { const s = L.regions[k].side(); ok[k] = !s || !side || s === side ? 1 : 0; }
    let out = new Uint16Array(nv);
    for (let v = 0; v < nv; v++) {
      for (const d of STEPS) {
        const p = apply(inv, [pos[3 * v] + d * nrm[3 * v], pos[3 * v + 1] + d * nrm[3 * v + 1], pos[3 * v + 2] + d * nrm[3 * v + 2]]);
        const i = Math.round(p[0]), j = Math.round(p[1]), k = Math.round(p[2]);
        if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz) continue;
        const r = L.index[i + nx * (j + ny * k)];
        if (r && ok[r]) { out[v] = r; break; }
      }
    }
    // Neighbours of each vertex, from the triangles.
    const tri = inst.mesh.geometry.index.array;
    const deg = new Uint32Array(nv + 1);
    for (let t = 0; t < tri.length; t++) deg[tri[t] + 1] += 2;
    for (let v = 0; v < nv; v++) deg[v + 1] += deg[v];
    const nb = new Uint32Array(deg[nv]);
    const fill = deg.slice(0, nv);
    for (let t = 0; t < tri.length; t += 3) {
      const a = tri[t], b = tri[t + 1], c = tri[t + 2];
      nb[fill[a]++] = b; nb[fill[a]++] = c; nb[fill[b]++] = a; nb[fill[b]++] = c; nb[fill[c]++] = a; nb[fill[c]++] = b;
    }
    const vote = (src, v) => {
      const counts = new Map();
      for (let e = deg[v]; e < deg[v + 1]; e++) { const r = src[nb[e]]; if (r) counts.set(r, (counts.get(r) || 0) + 1); }
      let best = src[v], bc = 0;
      for (const [r, c] of counts) if (c > bc) { best = r; bc = c; }
      return { best, bc, total: deg[v + 1] - deg[v] };
    };
    // A majority filter, then the gaps filled from their neighbours.
    const next = out.slice();
    for (let v = 0; v < nv; v++) { const { best, bc, total } = vote(out, v); if (bc > total / 2) next[v] = best; }
    out = next;
    for (let pass = 0; pass < 6; pass++) {
      let changed = 0;
      const n2 = out.slice();
      for (let v = 0; v < nv; v++) if (!out[v]) { const { best } = vote(out, v); if (best) { n2[v] = best; changed++; } }
      out = n2;
      if (!changed) break;
    }
    inst.parcels.set(L.spec.id, out);
    return out;
  }

  // The colours of the regions, as shared arrays (the viewer caches by identity).
  palette() {
    const L = this.current;
    if (this.colourBy === 'network' && L.spec.reader === 'schaefer') {
      const byNet = new Map();
      return L.regions.map((r, k) => {
        if (!k) return NEUTRAL;
        const net = L.network[k];
        if (!byNet.has(net)) {
          const root = r.root();
          byNet.set(net, root.colour && Array.isArray(root.colour) ? soften(root.colour.map((c) => c / 255)) : L.colours[k]);
        }
        return byNet.get(net);
      });
    }
    return L.colours.map((c, k) => (k ? c : NEUTRAL));
  }

  install() {
    const app = this.app, L = this.current;
    if (!L || !app.sections) return;
    const rgb = this.palette();
    this.rgb = rgb;
    for (const inst of app.instances.get('cerebral-cortex') || []) { inst.borderFlags = null; this.vertexRegions(inst); }
    app.cortexColour = (inst, i) => rgb[inst.parcels.get(L.spec.id)[i]];
    app.cortexKey = (inst, i) => (this.colourBy === 'network' ? L.network[inst.parcels.get(L.spec.id)[i]] || 0 : inst.parcels.get(L.spec.id)[i]);
    app.cortexRegion = (inst, hit) => {
      const k = this.regionAtHit(inst, hit);
      if (!k) return null;
      return { region: k, key: `parc:${k}`, id: null, name: this.regionName(k), sub: this.regionSub(k), colour: this.css(k) };
    };
    // The sections show the parcellation.
    if (!this.saved) this.saved = { labels: app.sections.labels, palette: app.sections.palette };
    app.sections.labels = { dims: this.saved.labels.dims, affine: this.saved.labels.affine, data: L.index };
    app.sections.palette = this.sectionPalette();
    this.light();
    this.legend();
    app.cardOverride = this.selected ? () => this.cardHtml() : () => this.introHtml();
    this.ctx.paintCortex();
    this.ctx.renderCard();
    this.ctx.renderLeft();
    for (const v of Object.values(app.views)) v.redraw();
    this.ctx.schedulePlate();
  }

  sectionPalette() {
    if (this.species === 'mouse') return this.mouseView.mouse ? this.mouseView.sectionPalette() : null;
    if (!this.current || !this.rgb) return null;
    const rgb = this.rgb;
    const pal = new Uint8Array(4 * rgb.length);
    rgb.forEach((c, k) => { if (k) pal.set([Math.round(255 * c[0]), Math.round(255 * c[1]), Math.round(255 * c[2]), 235], 4 * k); });
    return pal;
  }

  regionAtHit(inst, hit) {
    const L = this.current;
    const regions = inst.parcels && inst.parcels.get(L.spec.id);
    if (!regions || !hit.face) return 0;
    const pos = inst.mesh.geometry.attributes.position.array;
    let best = hit.face.a, bd = Infinity;
    for (const v of [hit.face.a, hit.face.b, hit.face.c]) {
      const d = (pos[3 * v] - hit.point[0]) ** 2 + (pos[3 * v + 1] - hit.point[1]) ** 2 + (pos[3 * v + 2] - hit.point[2]) ** 2;
      if (d < bd) { bd = d; best = v; }
    }
    return regions[best];
  }

  regionName(k) {
    const r = this.current.regions[k];
    if (!r) return '';
    if (this.current.spec.reader === 'schaefer') {
      const net = NETWORKS[this.current.network[k]];
      const side = r.side() === 'left' ? 'left' : 'right';
      return `${net ? net.name.replace(/ network$/, '') : this.current.network[k]}, parcel ${r.attributes.number || ''}${r.attributes.area ? ` (${r.attributes.area})` : ''}, ${side}`;
    }
    return r.name.replace(/ - (left|right) hemisphere$/, '');
  }

  regionSub(k) {
    const r = this.current.regions[k];
    const path = r.ancestors().map((a) => a.name).filter((n) => !/hemisphere$/.test(n));
    return [this.current.spec.name, path[0]].filter(Boolean).join(' · ');
  }

  css(k) { const c = this.rgb[k] || NEUTRAL; return `rgb(${c.map((x) => Math.round(255 * x)).join(',')})`; }

  // Light the chosen region on the cortex and in the sections.
  light() {
    const app = this.app, L = this.current, k = this.selected;
    if (!k) { app.regionLit = null; if (app.sections) app.sections.selected = null; return; }
    const lit = new Map();
    for (const inst of app.instances.get('cerebral-cortex') || []) {
      const regions = inst.parcels.get(L.spec.id);
      const flags = new Uint8Array(regions.length);
      for (let i = 0; i < regions.length; i++) flags[i] = regions[i] === k ? 1 : 0;
      lit.set(inst.key, flags);
    }
    app.regionLit = lit;
    app.sections.selected = new Set([k]);
    app.sections.dimOthers = false;
  }

  selectRegion(k, w = null) {
    const app = this.app;
    this.selected = k;
    this.ctx.state.selected = null;
    this.light();
    const st = this.current.stats[k];
    if (w && w.point) this.ctx.setPoint(w.point);
    else if (st) this.ctx.setPoint(st.centroid);
    app.cardOverride = k ? () => this.cardHtml() : () => this.introHtml();
    this.ctx.paintCortex();
    this.ctx.renderCard();
    for (const v of Object.values(app.views)) v.redraw();
    const list = document.querySelector('.region-list');
    if (list) list.innerHTML = this.listHtml();
    this.drawMinis();
  }

  // What the sections say at a point, and what a click there chooses.
  sectionAt(w) {
    if (this.species === 'mouse') return this.mouseView.mouse ? this.mouseView.sectionAt(w) : null;
    const k = this.app.sections.labelAt(w);
    if (!k) return null;
    return { name: this.regionName(k), sub: this.regionSub(k), colour: this.css(k), choose: () => this.selectRegion(k, { point: w }) };
  }
  sectionSelection() {
    if (this.species === 'mouse') return this.mouseView.mouse ? this.mouseView.sectionSelection() : null;
    return this.selected ? new Set([this.selected]) : null;
  }

  cardAction(a) {
    const app = this.app;
    if (a === 'region-focus') {
      const st = this.current.stats[this.selected];
      if (!st) return;
      import('three').then(({ Box3, Vector3 }) => app.viewer.frame(new Box3(new Vector3(...st.bounds.min), new Vector3(...st.bounds.max)), { margin: 1.9 }));
    }
    if (a === 'region-close') this.selectRegion(0);
  }

  // The card before an area is chosen: what the division shown is.
  introHtml() {
    const L = this.current;
    if (!L) return `<div class="card-empty"><h2>Divisions of the cortex</h2><p>The cortex divided otherwise than by its gyri: by the arrangement of its cells, or by the networks whose activity rises and falls together at rest.</p><p class="hint">Reading the atlas…</p></div>`;
    const isSchaefer = L.spec.reader === 'schaefer';
    return `<div class="card-empty"><h2>${esc(L.spec.name)}</h2>
<p>${esc(L.spec.description)}</p>
<p>${isSchaefer
    ? 'Where the gyri are the cortex’s geography, the networks are its politics: regions far apart that work together. A parcel is a patch of cortex whose activity at rest follows one course.'
    : 'Where the gyri are drawn by the folds of the surface, the cytoarchitectonic areas are drawn by the tissue itself: the size, density and layering of the neurons, which change at the borders between areas, often away from any sulcus.'}</p>
<p class="hint">Click the brain or an entry of the list to read ${isSchaefer ? 'a parcel' : 'an area'}: where it lies in the gross anatomy, how large it is, and what the atlas says of it.</p></div>`;
  }

  cardHtml() {
    const L = this.current, k = this.selected, r = L.regions[k];
    if (!r) return '';
    const st = L.stats[k];
    const app = this.app;
    const isSchaefer = L.spec.reader === 'schaefer';
    const path = r.ancestors().reverse().map((x) => x.name).filter((n) => !/hemisphere$/.test(n) && n !== this.regionName(k)).slice(1);
    const facts = [];
    if (st) {
      facts.push(['Volume', volumeText(st.voxels)]);
      facts.push(['Centre, MNI', `${coords(st.centroid)} <small>mm</small>`]);
      facts.push(['Talairach', coords(mniToTalairach(st.centroid))]);
    }
    // The gross anatomy beneath: the structures of the anatomy that share its voxels.
    const within = [];
    if (st) {
      let rest = 0;
      for (const [value, count] of st.overlap) {
        const key = app.labelValues.get(value);
        const at = key && structureOfKey(key);
        if (!at) { rest += count; continue; }
        if (['cerebral-white-matter', 'csf', 'vessel', 'white-matter-hypointensities', 'cerebral-cortex'].includes(at.id)) { rest += count; continue; }
        const prev = within.find((x) => x.id === at.id);
        if (prev) prev.count += count; else within.push({ id: at.id, side: at.side, count });
      }
      within.sort((a, b) => b.count - a.count);
    }
    const total = st ? st.voxels : 1;
    const chips = within.filter((x) => x.count / total >= 0.04).slice(0, 7).map((x) => {
      const s = this.ctx.structure(x.id);
      return `<button type="button" data-anatomy="${esc(x.id)}" data-side="${esc(x.side)}"><i style="background:${esc(this.ctx.colourOf(x.id) || '#ccc')}"></i>${esc(s.name)} <small>${Math.round((100 * x.count) / total)} %</small></button>`;
    }).join('');
    const net = isSchaefer ? NETWORKS[L.network[k]] : null;
    const html = [];
    html.push('<div class="card-head">');
    html.push(`<div class="over">${esc(L.spec.name)}${path.length ? ` <span class="sep">›</span> ${path.map(esc).join(' <span class="sep">›</span> ')}` : ''}</div>`);
    html.push(`<h2><span class="sw" style="background:${this.css(k)}"></span><span>${esc(this.regionName(k))}</span></h2>`);
    if (isSchaefer) html.push(`<div class="latin">${esc(r.name)}</div>`);
    if (r.side()) html.push(`<div class="side">${r.side() === 'left' ? 'Left' : 'Right'} hemisphere</div>`);
    html.push('<button type="button" class="close" data-action="region-close" aria-label="Close"><svg><use href="#i-close"/></svg></button>');
    html.push('<div class="card-actions"><button type="button" data-action="region-focus"><svg><use href="#i-focus"/></svg>Focus</button></div>');
    html.push('</div><div class="card-body">');
    if (net) html.push(`<p class="summary"><b>${esc(net.name)}.</b> ${esc(net.text)}</p>`);
    else html.push(`<p class="summary">${julichGloss(r.name)}</p>`);
    if (facts.length) html.push(`<div class="facts">${facts.map(([a, b]) => `<div><div class="k">${esc(a)}</div><div class="v">${b}</div></div>`).join('')}</div>`);
    if (chips) html.push(`<div class="sect"><h4>Lies in, by the anatomy</h4><div class="chips">${chips}</div><p class="note">The share of the ${isSchaefer ? 'parcel' : 'area'}’s voxels in each structure of the anatomy. Choose one to read about it.</p></div>`);
    html.push('<div class="sect"><h4>In the sections</h4><div class="mini-sections"><div><canvas data-mini="axial"></canvas><span>axial</span></div><div><canvas data-mini="coronal"></canvas><span>coronal</span></div><div><canvas data-mini="sagittal"></canvas><span>sagittal</span></div></div></div>');
    html.push(`<div class="sect"><h4>The atlas</h4><p>${esc(L.spec.description)}</p></div>`);
    const cites = isSchaefer ? ['Schaefer2018', 'Yeo2011'] : ['Amunts2020', 'Amunts2015'];
    html.push(`<div class="card-refs">${cites.map((c) => `<p>${esc(CITATIONS[c] || c)}</p>`).join('')}<p>Licence of the data: ${esc(L.spec.licence)}.</p></div>`);
    html.push('</div>');
    setTimeout(() => this.bindCard(), 0);
    return html.join('');
  }

  bindCard() {
    const el = document.getElementById('card');
    el.querySelectorAll('[data-anatomy]').forEach((b) => b.addEventListener('click', (e) => {
      e.stopPropagation();
      const id = b.dataset.anatomy, side = b.dataset.side;
      this.ctx.setMode('anatomy');
      this.ctx.select({ id, side: side && side !== 'mid' ? side : null });
    }));
    this.drawMinis();
  }

  drawMinis() {
    const st = this.current && this.current.stats[this.selected];
    if (!st || !this.app.sections) return;
    const axes = { axial: 2, coronal: 1, sagittal: 0 };
    for (const c of document.querySelectorAll('#card canvas[data-mini]')) {
      const plane = c.dataset.mini;
      this.app.sections.draw(c, plane, st.centroid[axes[plane]], { scale: 2, flipU: plane === 'sagittal' });
    }
  }

  legend() {
    const el = document.getElementById('legend');
    const L = this.current;
    if (!L) { el.innerHTML = ''; return; }
    if (L.spec.reader === 'schaefer' && this.colourBy === 'network') {
      const keys = Object.entries(NETWORKS).map(([acr, n]) => {
        const k = L.network.findIndex((x) => x === acr);
        return k > 0 ? `<span><i style="background:${this.css(k)}"></i>${esc(n.name.replace(/ network$/, ''))}</span>` : '';
      }).join('');
      el.innerHTML = `<div><div class="title"><span>Yeo’s seven networks</span></div><div class="keys">${keys}</div></div>`;
    } else {
      el.innerHTML = `<div><div class="title"><span>${esc(L.spec.name)}</span><span class="unit">${L.regions.length - 1} regions</span></div><div class="note">${L.spec.reader === 'schaefer' ? 'Each parcel in a shade of its network’s colour.' : 'Colours tell neighbours apart; the two sides of an area share one.'}</div></div>`;
    }
  }

  // ---- the panel ------------------------------------------------------------------------

  panel() {
    if (!this.available) return '<h3>Atlases</h3><p class="pad note">Looking for the atlases in the data library…</p>';
    if (!this.available.length) return `<h3>Atlases</h3><p class="pad note">No atlas is fetched yet. From the repository root, run <code>python data/acquisition/fetch.py ebrains-julich-brain schaefer-2018</code>, then reload.</p>`;
    const species = this.mouseAvailable ? `<h3>Species</h3><div class="pad seg species">${['human', 'mouse'].map((s) => `<button type="button" class="opt" data-species="${s}" aria-pressed="${this.species === s}">${s === 'human' ? 'Human' : 'Mouse'}</button>`).join('')}</div>` : '';
    if (this.species === 'mouse') return species + this.mouseView.panelHtml();
    const L = this.current;
    const cards = this.available.map((p) => `<button type="button" class="atlas-choice" data-parcellation="${esc(p.id)}" aria-pressed="${Boolean(L && L.spec.id === p.id)}">
<b>${esc(p.name)}</b><span class="what">${esc(p.short)}</span><span class="meta">${esc(p.citation)} · ${esc(p.licence)}</span></button>`).join('');
    const colourBy = L && L.spec.reader === 'schaefer' ? `<h3>Coloured by</h3><div class="pad seg"><button type="button" class="opt" data-colour-by="parcel" aria-pressed="${this.colourBy === 'parcel'}">Parcel</button><button type="button" class="opt" data-colour-by="network" aria-pressed="${this.colourBy === 'network'}">Network</button></div>` : '';
    return `${species}<h3>Divisions of the cortex</h3><div class="pad atlas-choices">${cards}</div>
${this.busy ? '<p class="pad note">Reading…</p>' : ''}${colourBy}
<h3>${L ? (L.spec.reader === 'schaefer' ? 'Parcels' : 'Areas') : 'Regions'}</h3>
<div class="pad"><input type="search" class="region-filter" placeholder="Filter…" value="${esc(this.filter)}" autocomplete="off" spellcheck="false"></div>
<div class="region-list">${L ? this.listHtml() : ''}</div>`;
  }

  // The list of regions: Julich-Brain's areas by lobe or nucleus group, one row an area with
  // its two sides; Schaefer's parcels by network, the left before the right.
  listHtml() {
    const L = this.current;
    const q = this.filter.trim().toLowerCase();
    const isSchaefer = L.spec.reader === 'schaefer';
    const groups = new Map();
    const rows = new Map();
    for (let k = 1; k < L.regions.length; k++) {
      const r = L.regions[k];
      const name = this.regionName(k);
      if (q && !name.toLowerCase().includes(q) && !r.name.toLowerCase().includes(q)) continue;
      if (isSchaefer) {
        const g = NETWORKS[L.network[k]]?.name || L.network[k];
        if (!groups.has(g)) groups.set(g, []);
        groups.get(g).push({ name: `${r.side() === 'left' ? 'Left' : 'Right'} ${r.attributes.number || ''}${r.attributes.area ? ` · ${r.attributes.area}` : ''}`, sides: { one: k }, order: (r.side() === 'left' ? 0 : 1000) + (r.attributes.number || 0) });
        continue;
      }
      const area = r.side() && r.parent ? r.parent : r;
      const g = groupOf(r);
      const label = g ? g[0].toUpperCase() + g.slice(1) : 'Other';
      if (!rows.has(area)) {
        const row = { name, sides: {}, order: 0 };
        rows.set(area, row);
        if (!groups.has(label)) groups.set(label, []);
        groups.get(label).push(row);
      }
      rows.get(area).sides[r.side() || 'one'] = k;
    }
    const out = [];
    const sel = this.selected;
    // Lobes from front to back, then the limbic lobe and the insula, then the nuclei; the gaps of the map last.
    const ORDER = ['Visual network', 'Somatomotor network', 'Dorsal attention network', 'Salience and ventral attention network', 'Limbic network', 'Frontoparietal control network', 'Default mode network', 'Frontal lobe', 'Parietal lobe', 'Temporal lobe', 'Occipital lobe', 'Limbic lobe', 'Insula', 'Basal ganglia', 'Basal forebrain', 'Amygdala', 'Dorsal thalamus', 'Ventral thalamus', 'Metathalamus', 'STN (Subthalamus)', 'Midbrain', 'Cerebellar nuclei'];
    const rank = (g) => (/GapMap/.test(g) ? 100 : ORDER.includes(g) ? ORDER.indexOf(g) : 50);
    const ordered = [...groups].sort((x, y) => rank(x[0]) - rank(y[0]));
    for (const [g, list] of ordered) {
      list.sort((x, y) => x.order - y.order || x.name.localeCompare(y.name, 'en', { numeric: true }));
      out.push(`<div class="group">${esc(g)}</div>`);
      for (const row of list) {
        const first = row.sides.left || row.sides.one || row.sides.right;
        const on = Object.values(row.sides).includes(sel);
        const sides = isSchaefer ? '' : ['left', 'right'].filter((sd) => row.sides[sd]).map((sd) => `<span class="lr${row.sides[sd] === sel ? ' on' : ''}" data-region="${row.sides[sd]}" title="${sd} hemisphere">${sd[0].toUpperCase()}</span>`).join('');
        out.push(`<button type="button" data-region="${first}" class="${on ? 'selected' : ''}"><i style="background:${this.css(first)}"></i><span class="nm">${esc(row.name)}</span>${sides}</button>`);
      }
    }
    return out.join('') || '<p class="note" style="padding:4px 8px">Nothing by that name.</p>';
  }

  bind(left) {
    left.querySelectorAll('[data-species]').forEach((b) => b.addEventListener('click', () => this.setSpecies(b.dataset.species)));
    if (this.species === 'mouse') { this.mouseView.bindPanel(left); return; }
    left.querySelectorAll('[data-parcellation]').forEach((b) => b.addEventListener('click', () => this.choose(b.dataset.parcellation)));
    left.querySelectorAll('[data-colour-by]').forEach((b) => b.addEventListener('click', () => { this.colourBy = b.dataset.colourBy; this.install(); }));
    const filter = left.querySelector('.region-filter');
    if (filter) filter.addEventListener('input', () => { this.filter = filter.value; left.querySelector('.region-list').innerHTML = this.listHtml(); });
    const list = left.querySelector('.region-list');
    if (list) list.addEventListener('click', (e) => { const b = e.target.closest('[data-region]'); if (b) this.selectRegion(Number(b.dataset.region)); });
  }
}

// A Julich-Brain area's name read aloud: the letters of its short name say where it lies,
// the brackets its older names. (Amunts et al. 2020; the EBRAINS area descriptions.)
const PREFIXES = [
  [/GapMap/, 'a part of the cortex whose areas are not yet mapped; Julich-Brain fills it provisionally, so that the map covers the whole cortex'],
  [/^Area hOc/, 'an area of the occipital cortex (hOc: human occipital cortex); hOc1 is the primary visual cortex, V1'],
  [/^Area FG/, 'an area of the fusiform gyrus (FG), in the ventral visual stream'],
  [/^Area hIP/, 'an area of the intraparietal sulcus (hIP: human intraparietal)'],
  [/^Area hPO/, 'an area of the parieto-occipital sulcus (hPO)'],
  [/^Area P[FG]/, 'an area of the inferior parietal lobule, named in von Economo’s letters (PF, PG and their parts)'],
  [/^Area [57][A-Za-z]* /, 'an area of the superior parietal lobule, numbered after Brodmann'],
  [/^Area Op/, 'an area of the parietal operculum (Op), the second somatosensory cortex (Eickhoff et al. 2006)'],
  [/^Area (Id|Ig|Ia)/, 'an area of the insula: agranular (Ia), dysgranular (Id) or granular (Ig), by whether it has an inner granular layer'],
  [/^Area Fo/, 'an area of the orbitofrontal cortex (Fo: frontal orbital)'],
  [/^Area Fp/, 'an area of the frontal pole (Fp)'],
  [/^Area Te /, 'an area of the auditory cortex: Te 1 on Heschl’s gyrus is the primary auditory cortex (Morosan et al. 2001), Te 2 and Te 3 lie around it on the superior temporal gyrus'],
  [/^Area (TE|TI|TeI|STS)/, 'an area of the temporal cortex, around the auditory cortex or in the superior temporal sulcus'],
  [/^Area (Ph|CoS|OTS)/, 'an area of the ventral temporal cortex: of the parahippocampal gyrus (Ph), the collateral sulcus (CoS) or the occipitotemporal sulcus (OTS)'],
  [/^Area TPJ/, 'an area of the temporoparietal junction, where the superior temporal and supramarginal gyri meet'],
  [/^Area 4[ap] /, 'part of the primary motor cortex, Brodmann’s area 4: 4a on the crown of the precentral gyrus, 4p in the depth of the central sulcus (Geyer et al. 1996)'],
  [/^Area 6[a-z0-9]* /, 'part of the premotor cortex, Brodmann’s area 6, or of the supplementary motor areas on its medial side'],
  [/^Area (3a|3b|1|2) /, 'part of the primary somatosensory cortex, Brodmann’s areas 3a, 3b, 1 and 2'],
  [/^Area (44|45) /, 'part of Broca’s region, Brodmann’s areas 44 and 45'],
  [/^Area [aps]?(24|25|32|33)[a-z]* /, 'an area of the cingulate cortex, numbered after Brodmann'],
  [/^Area (8|9|46)[a-z0-9]* /, 'an area of the dorsolateral prefrontal cortex'],
  [/^Area (IFS|IFJ|MFG|SFS|SFG)/, 'an area of the lateral prefrontal cortex, named by the sulcus or gyrus where it lies: IFS the inferior frontal sulcus, IFJ the inferior frontal junction, MFG, SFS and SFG the middle frontal gyrus, the superior frontal sulcus and gyrus'],
  [/^(CA|DG|Subiculum|Entorhinal|HATA|Area EC|Area TR|TrS)|\(Hippocampus/, 'a part of the hippocampal formation or of the cortex beside it'],
  [/Amygdala/i, 'a nucleus group of the amygdala'],
  [/^Ch /, 'a cholinergic nucleus of the basal forebrain (Ch1–Ch4 of Mesulam et al. 1983)'],
  [/Ventral Striatum|Ventral Pallidum|Basal Forebrain/, 'a part of the ventral striatum, the ventral pallidum or the basal forebrain'],
  [/\((Thalamus|Metathalamus)/, 'a nucleus of the thalamus'],
  [/^STN/, 'the subthalamic nucleus'],
  [/\(Midbrain/, 'a nucleus of the midbrain'],
  [/\(Cerebellum/, 'a deep nucleus of the cerebellum'],
];
function julichGloss(name) {
  const hit = PREFIXES.find(([re]) => re.test(name));
  if (/GapMap/.test(name)) return `Not an area but ${esc(hit[1])}. Its name says which mapped regions it lies between.`;
  const what = hit ? hit[1] : 'a part of the brain told apart by the size, density and layering of its cells';
  const naming = /^Area /.test(name)
    ? 'Its name gives, in brackets, the gyrus or sulcus where it lies and, where there is one, its classical name or Brodmann’s number.'
    : 'Its name gives its abbreviation and, in brackets, the region and the full name.';
  return `A cytoarchitectonic area of Julich-Brain: ${esc(what)}. ${naming} Its border was drawn where the arrangement of the cells changes, in ten brains, and the map gives each voxel to the area most likely there.`;
}
