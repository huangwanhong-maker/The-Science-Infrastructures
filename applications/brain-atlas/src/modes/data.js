// The Data mode: numbers mapped onto the brain. Values by structure, read from a
// table someone pastes (by name, by FreeSurfer's name, or by an area of the atlas
// shown in the Atlases mode), colour the structures in three dimensions and in the
// sections; points given by coordinates appear as markers, each with what lies there.
// Two measures of the atlas itself serve as examples: the volume of each structure
// and the asymmetry between its two sides.

import { parseTable, parsePoints, valueRange, tidy } from '../ui/data.js';
import { colourAt, cssGradient, COLORMAP_NAMES } from '../render/colormaps.js';
import { esc } from '../ui/card.js';
import {
  BRAIN_STRUCTURES, structureOfKey, keysOf, knowledgeOf, talairachToMni, mniToTalairach,
} from '../../../../domains/neuroscience/js/index.js';

const NONE = [0.42, 0.39, 0.43];
const signed = (x) => (x < 0 ? `−${-x}` : `${x}`);
const coords = (c) => c.map((x) => signed(Math.round(x))).join(', ');

// FreeSurfer's names of the segmentation's labels (FreeSurferColorLUT.txt), for tables
// such as aseg.stats.
const FREESURFER = {
  2: 'Left-Cerebral-White-Matter', 3: 'Left-Cerebral-Cortex', 4: 'Left-Lateral-Ventricle', 5: 'Left-Inf-Lat-Vent',
  7: 'Left-Cerebellum-White-Matter', 8: 'Left-Cerebellum-Cortex', 10: 'Left-Thalamus', 11: 'Left-Caudate', 12: 'Left-Putamen',
  13: 'Left-Pallidum', 14: '3rd-Ventricle', 15: '4th-Ventricle', 16: 'Brain-Stem', 17: 'Left-Hippocampus', 18: 'Left-Amygdala',
  24: 'CSF', 26: 'Left-Accumbens-area', 28: 'Left-VentralDC', 30: 'Left-vessel', 31: 'Left-choroid-plexus',
  41: 'Right-Cerebral-White-Matter', 42: 'Right-Cerebral-Cortex', 43: 'Right-Lateral-Ventricle', 44: 'Right-Inf-Lat-Vent',
  46: 'Right-Cerebellum-White-Matter', 47: 'Right-Cerebellum-Cortex', 49: 'Right-Thalamus', 50: 'Right-Caudate', 51: 'Right-Putamen',
  52: 'Right-Pallidum', 53: 'Right-Hippocampus', 54: 'Right-Amygdala', 58: 'Right-Accumbens-area', 60: 'Right-VentralDC',
  62: 'Right-vessel', 63: 'Right-choroid-plexus', 77: 'WM-hypointensities', 85: 'Optic-Chiasm',
  251: 'CC_Posterior', 252: 'CC_Mid_Posterior', 253: 'CC_Central', 254: 'CC_Mid_Anterior', 255: 'CC_Anterior',
};
const FS_ALIASES = { 'left-thalamus-proper': 10, 'right-thalamus-proper': 49 };

const norm = (s) => String(s).toLowerCase().replace(/[’']/g, '').replace(/[^a-z0-9]+/g, ' ').trim();

// The side a name states, and the name without it.
function sideOf(text) {
  const t = ` ${norm(text)} `;
  const m = t.match(/ (left|right|lh|rh|l|r) /);
  if (!m) return { side: null, base: norm(text) };
  const side = ['left', 'lh', 'l'].includes(m[1]) ? 'left' : 'right';
  return { side, base: t.replace(m[0], ' ').trim().replace(/\s+/g, ' ') };
}

export class DataMode {
  constructor(ctx) {
    this.ctx = ctx;
    this.text = '';
    this.pointsText = '';
    this.space = 'mni';
    this.colormap = 'thermal';
    this.values = null;          // Map key ('aseg:12', 'ho:13', or 'parc:k') -> value
    this.byAtlas = false;        // whether the keys are regions of the Atlases mode's parcellation
    this.range = null;
    this.title = '';
    this.unit = '';
    this.report = '';
    this.points = [];
    this.pointsReport = '';
    this.chosenPoint = -1;
    this.names = null;
  }

  get app() { return this.ctx.app; }

  // Every name a structure is known by, normalised, with its keys by side.
  index() {
    if (this.names) return this.names;
    const names = new Map();
    const put = (name, entry) => { const n = norm(name); if (n && !names.has(n)) names.set(n, entry); };
    for (const s of BRAIN_STRUCTURES) {
      const k = knowledgeOf(s.id).entry;
      const entry = { id: s.id };
      for (const n of [s.name, s.short, s.id.replace(/-/g, ' '), k?.latin, ...(k?.also || [])]) if (n) put(n, entry);
    }
    for (const [v, name] of Object.entries(FREESURFER)) put(name, { key: `aseg:${v}` });
    for (const [name, v] of Object.entries(FS_ALIASES)) put(name, { key: `aseg:${v}` });
    this.names = names;
    return names;
  }

  // The label keys a name in a table stands for: [] when nothing matches.
  keysFor(text) {
    const names = this.index();
    const direct = names.get(norm(text));
    const fromEntry = (entry, side) => {
      // FreeSurfer's cortex labels stand for every gyrus of their side.
      if (entry.key === 'aseg:3' || entry.key === 'aseg:42') return [...keysOf('cerebral-cortex', 'label')[entry.key === 'aseg:3' ? 'left' : 'right']];
      if (entry.key) return [entry.key];
      const k = keysOf(entry.id, 'label');
      const sides = side ? [side, 'mid'] : ['left', 'right', 'mid'];
      return sides.flatMap((sd) => [...k[sd]]);
    };
    if (direct) return fromEntry(direct, null);
    const { side, base } = sideOf(text);
    const e = names.get(base);
    return e ? fromEntry(e, side) : [];
  }

  // Values from the table: by the anatomy's structures, or by the regions of the parcellation
  // the Atlases mode has loaded, whichever matches more of the lines.
  mapValues(text, { title = 'Values', unit = '' } = {}) {
    const rows = parseTable(text).filter((r) => r.length >= 2);
    const anatomy = new Map(), unmatched = [];
    const atlas = this.app.modes.atlases?.current || null;
    const byRegion = new Map(), unmatchedAtlas = [];
    for (const row of rows) {
      const name = row[0], value = Number(row[row.length - 1]);
      if (!Number.isFinite(value)) { unmatched.push(name); unmatchedAtlas.push(name); continue; }
      const keys = this.keysFor(name);
      if (keys.length) for (const k of keys) anatomy.set(k, value); else unmatched.push(name);
      if (atlas) {
        const r = atlas.parcellation.match(name);
        const ks = r ? r.leaves().map((leaf) => atlas.indexOf.get(leaf)).filter(Boolean) : [];
        if (ks.length) for (const k of ks) byRegion.set(`parc:${k}`, value); else unmatchedAtlas.push(name);
      }
    }
    const useAtlas = atlas && byRegion.size && unmatchedAtlas.length < unmatched.length;
    this.byAtlas = Boolean(useAtlas);
    this.values = useAtlas ? byRegion : anatomy;
    const miss = useAtlas ? unmatchedAtlas : unmatched;
    const matched = rows.length - miss.length;
    this.range = valueRange(this.values.values());
    this.title = title; this.unit = unit;
    this.report = `${matched ? `<span class="good">${matched} ${matched === 1 ? 'line' : 'lines'} matched</span>${useAtlas ? ` to ${esc(atlas.spec.name)}` : ''}` : '<span class="bad">Nothing matched</span>'}${miss.length ? `; <span class="bad">not found: ${esc(miss.slice(0, 8).join(', '))}${miss.length > 8 ? '…' : ''}</span>` : ''}.`;
    this.install();
  }

  // The examples: the volume of each structure, and the asymmetry of each pair.
  builtIn(which) {
    const regions = this.app.product.regions;
    const lines = [];
    const seen = new Set();
    for (const s of BRAIN_STRUCTURES) {
      if (!s.label || ['csf', 'vessel', 'white-matter-hypointensities', 'cerebral-white-matter', 'cerebral-cortex'].includes(s.id)) continue;
      const vol = (key) => { const r = key && regions[key]; return r ? (r.joined ? r.joined.volume : r.volume) : null; };
      if (which === 'volume') {
        for (const side of ['left', 'right', 'mid']) {
          const v = vol(s.label[side]);
          if (v && !seen.has(s.label[side])) { seen.add(s.label[side]); lines.push(`"${side === 'mid' ? '' : `${side} `}${s.short}", ${(v / 1000).toFixed(2)}`); }
        }
      } else if (which === 'asymmetry' && s.label.left && s.label.right) {
        const l = vol(s.label.left), r = vol(s.label.right);
        if (l && r) lines.push(`"${s.short}", ${((200 * (l - r)) / (l + r)).toFixed(1)}`);
      }
    }
    this.text = lines.join('\n');
    if (which === 'volume') this.mapValues(this.text, { title: 'Volume', unit: 'cm³' });
    else { this.colormap = 'balance'; this.mapValues(this.text, { title: 'Asymmetry, 2(L − R)/(L + R)', unit: '%' }); this.range = valueRange(this.values.values(), { diverging: true }); this.install(); }
    this.ctx.renderLeft();
  }

  colourOfValue(v) {
    const [lo, hi] = this.range;
    return colourAt(this.colormap, (v - lo) / (hi - lo)).map((c) => c / 255);
  }

  enter() { this.install(); }

  leave() {
    const app = this.app;
    app.cortexColour = null; app.cortexKey = null; app.meshColour = null; app.cardFacts = null;
    for (const inst of app.instances.get('cerebral-cortex') || []) inst.borderFlags = null;
    app.viewer.setMarkers([]);
    document.getElementById('legend').innerHTML = '';
    this.ctx.paintAll();
  }

  install() {
    const app = this.app;
    if (this.ctx.state.mode !== 'data') return;
    const values = this.values;
    const cache = new Map();
    const rgb = (v) => { if (!cache.has(v)) cache.set(v, this.colourOfValue(v)); return cache.get(v); };
    if (values && this.byAtlas) {
      const atlases = app.modes.atlases;
      const L = atlases.current;
      for (const inst of app.instances.get('cerebral-cortex') || []) atlases.vertexRegions(inst);
      app.cortexColour = (inst, i) => { const v = values.get(`parc:${inst.parcels.get(L.spec.id)[i]}`); return v === undefined ? NONE : rgb(v); };
      app.cortexKey = (inst, i) => inst.parcels.get(L.spec.id)[i];
      app.meshColour = () => '#6b636d';
    } else if (values) {
      app.cortexColour = (inst, i) => { const v = values.get(`ho:${inst.labels[i]}`); return v === undefined ? NONE : rgb(v); };
      app.cortexKey = null;
      app.meshColour = (inst) => {
        const keys = [inst.key, inst.region?.label].filter(Boolean);
        const at = structureOfKey(inst.key, 'mesh');
        const labelKey = at && this.ctx.structure(at.id).label?.[at.side];
        const v = values.get(labelKey) ?? values.get(keys[0]);
        if (v === undefined) return '#6b636d';
        const c = rgb(v);
        return `#${c.map((x) => Math.round(255 * x).toString(16).padStart(2, '0')).join('')}`;
      };
    } else {
      app.cortexColour = null; app.cortexKey = null; app.meshColour = null;
    }
    app.cardFacts = (sel) => this.facts(sel);
    for (const inst of app.instances.get('cerebral-cortex') || []) inst.borderFlags = null;
    this.ctx.paintAll();
    this.ctx.applyStyles();
    this.markers();
    this.legend();
    if (app.sections) { app.sections.palette = this.ctx.labelPalette(); for (const v of Object.values(app.views)) v.redraw(); }
    this.ctx.renderCard();
  }

  // The palette of the sections: each label coloured by its value; the rest faint.
  sectionPalette() {
    if (!this.values || this.byAtlas) return null;
    const table = this.app.product.labelTable;
    const pal = new Uint8Array(4 * (table.length + 1));
    for (const [value, key] of table) {
      const v = this.values.get(key);
      if (v === undefined) { pal.set([110, 100, 112, 0], 4 * value); continue; }
      const c = this.colourOfValue(v).map((x) => Math.round(255 * x));
      pal.set([c[0], c[1], c[2], 255], 4 * value);
    }
    return pal;
  }

  // The value of a selected structure, for its card.
  facts(sel) {
    if (!this.values || this.byAtlas || !sel || !sel.id) return [];
    const k = keysOf(sel.id, 'label');
    const out = [];
    for (const side of sel.side ? [sel.side, 'mid'] : ['left', 'right', 'mid']) {
      for (const key of k[side]) {
        const v = this.values.get(key);
        if (v !== undefined) out.push([`${this.title}${side !== 'mid' && !sel.side ? `, ${side}` : ''}`, `${tidy(v)} <small>${esc(this.unit)}</small>`]);
      }
    }
    return out.slice(0, 3);
  }

  legend() {
    const el = document.getElementById('legend');
    if (!this.values || !this.values.size) { el.innerHTML = this.points.length ? '' : ''; return; }
    const [lo, hi] = this.range;
    el.innerHTML = `<div><div class="title"><span>${esc(this.title)}</span><span class="unit">${esc(this.unit)}</span></div>
<div class="bar" style="background:${cssGradient(this.colormap)}"></div>
<div class="ticks"><span style="left:0%">${tidy(lo)}</span><span style="left:50%">${tidy((lo + hi) / 2)}</span><span style="left:100%">${tidy(hi)}</span></div>
<div class="note">${this.values.size} ${this.byAtlas ? 'regions' : 'labels'} coloured; the rest grey.</div></div>`;
  }

  // ---- points ---------------------------------------------------------------------------

  showPoints() {
    const { points, bad } = parsePoints(this.pointsText);
    const toMni = this.space === 'mni' ? (p) => p : (p) => talairachToMni(p, { method: this.space === 'brett' ? 'brett' : 'lancaster', software: this.space === 'spm' ? 'spm' : this.space === 'fsl' ? 'fsl' : 'other' });
    this.points = points.map((p) => ({ ...p, mni: toMni(p.position) }));
    const vals = this.points.map((p) => p.value).filter((v) => v !== null);
    this.pointRange = vals.length ? valueRange(vals) : null;
    this.pointsReport = `${this.points.length ? `<span class="good">${this.points.length} ${this.points.length === 1 ? 'point' : 'points'}</span>` : '<span class="bad">No points</span>'}${bad.length ? `; <span class="bad">not read: ${esc(bad.slice(0, 4).join(' | '))}</span>` : ''}.`;
    if (this.points.length && this.ctx.state.cortexStyle === 'solid') { this.ctx.state.cortexStyle = 'glass'; this.ctx.applyStyles(); this.ctx.renderToolbar(); }
    this.markers();
    this.ctx.renderLeft();
  }

  markers() {
    if (this.ctx.state.mode !== 'data') return;
    const list = this.points.map((p, i) => {
      let colour = '#f6c47e';
      if (p.value !== null && this.pointRange) { const c = colourAt('orchid', (p.value - this.pointRange[0]) / (this.pointRange[1] - this.pointRange[0] || 1)); colour = `rgb(${c.map(Math.round).join(',')})`; }
      return { position: p.mni, colour, radius: i === this.chosenPoint ? 3.6 : 2.6 };
    });
    this.app.viewer.setMarkers(list);
  }

  // A marker under the pointer: chosen on click, named on hover.
  pick(event) {
    const i = this.app.viewer.pickMarker(event.clientX, event.clientY);
    if (i < 0) return false;
    this.choosePoint(i);
    return true;
  }
  hoverInfo(pointer) {
    const i = this.app.viewer.pickMarker(pointer.x, pointer.y);
    if (i < 0) return null;
    const p = this.points[i];
    return { name: p.name || `Point ${i + 1}`, sub: `MNI ${coords(p.mni)} mm${p.value !== null ? ` · ${tidy(p.value)}` : ''}`, colour: '#f6c47e' };
  }

  // What lies at a point: the structure of the joined labels there, or the nearest within 6 mm.
  whatIsAt(p) {
    const s = this.app.sections;
    if (!s) return null;
    const ignore = new Set(['cerebral-white-matter', 'csf', 'vessel', 'white-matter-hypointensities', 'cerebral-cortex']);
    const at = (q) => { const key = this.app.labelValues.get(s.labelAt(q)); const st = key && structureOfKey(key); return st && !ignore.has(st.id) ? st : null; };
    const here = at(p);
    if (here) return { ...here, distance: 0 };
    let best = null;
    for (let r = 1; r <= 6 && !best; r++) {
      for (let dx = -r; dx <= r; dx++) for (let dy = -r; dy <= r; dy++) for (let dz = -r; dz <= r; dz++) {
        if (Math.max(Math.abs(dx), Math.abs(dy), Math.abs(dz)) !== r) continue;
        const st = at([p[0] + dx, p[1] + dy, p[2] + dz]);
        if (st) { const d = Math.hypot(dx, dy, dz); if (!best || d < best.distance) best = { ...st, distance: d }; }
      }
    }
    const raw = this.app.labelValues.get(s.labelAt(p));
    const tissue = raw && structureOfKey(raw);
    return best ? { ...best, tissue: tissue ? tissue.id : null } : (tissue ? { ...tissue, distance: 0 } : null);
  }

  choosePoint(i) {
    this.chosenPoint = i;
    const p = this.points[i];
    this.ctx.setPoint(p.mni);
    this.markers();
    const app = this.app;
    app.cardOverride = () => this.pointCard(i);
    this.ctx.renderCard();
    app.cardOverride = null;
  }

  pointCard(i) {
    const p = this.points[i];
    const w = this.whatIsAt(p.mni);
    const st = w && this.ctx.structure(w.id);
    const tal = mniToTalairach(p.mni);
    const facts = [['MNI', `${coords(p.mni)} <small>mm</small>`], ['Talairach', coords(tal)]];
    if (p.value !== null) facts.push(['Value', tidy(p.value)]);
    const where = st
      ? `<button type="button" data-select="${esc(st.id)}" data-side="${esc(w.side === 'mid' ? '' : w.side)}"><i style="background:${esc(this.ctx.colourOf(st.id) || '#ccc')}"></i>${esc(st.name)}${w.side && w.side !== 'mid' ? ` <small>${esc(w.side)}</small>` : ''}</button>`
      : '';
    const how = !w ? 'Outside the labelled brain.' : w.distance ? `The point lies in ${w.tissue === 'cerebral-white-matter' ? 'white matter' : 'unlabelled tissue'}; the nearest structure is ${Math.round(w.distance)} mm away:` : 'The point lies in:';
    return `<div class="card-head"><div class="over">A point of the data</div><h2><span class="sw" style="background:#f6c47e"></span><span>${esc(p.name || `Point ${i + 1}`)}</span></h2>
<button type="button" class="close" data-action="deselect" aria-label="Close"><svg><use href="#i-close"/></svg></button></div>
<div class="card-body"><div class="facts">${facts.map(([k, v]) => `<div><div class="k">${esc(k)}</div><div class="v">${v}</div></div>`).join('')}</div>
<div class="sect"><h4>Where it lies</h4><p>${how}</p><div class="chips">${where}</div></div>
<p class="note">Talairach coordinates by Brett’s transform; given in ${esc({ mni: 'MNI', brett: 'Talairach (Brett)', spm: 'Talairach (Lancaster, SPM)', fsl: 'Talairach (Lancaster, FSL)', other: 'Talairach (Lancaster)' }[this.space])} coordinates.</p></div>`;
  }

  // ---- the panel ------------------------------------------------------------------------

  panel() {
    const maps = COLORMAP_NAMES.map((n) => `<option value="${n}"${n === this.colormap ? ' selected' : ''}>${n}</option>`).join('');
    const atlas = this.app.modes.atlases?.current;
    return `<h3>Values by structure</h3>
<p class="pad note">A table pasted here, a structure and a value on each line, colours the brain. Names may be the atlas’s (“left putamen”, “Precentral gyrus”), FreeSurfer’s (“Left-Hippocampus”, as in aseg.stats)${atlas ? `, or the regions of <b>${esc(atlas.spec.name)}</b>` : ''}.</p>
<div class="pad data-form">
<textarea id="data-values" spellcheck="false" placeholder="structure, value&#10;left hippocampus, 3.9&#10;Right-Amygdala, 1.6&#10;Precentral gyrus, 2.2">${esc(this.text)}</textarea>
<div class="row-btns"><button type="button" class="btn primary" data-data="map">Colour the brain</button><button type="button" class="btn" data-data="clear">Clear</button>
<label class="inline">Colours <select id="data-colormap">${maps}</select></label></div>
<div class="report">${this.report}</div>
<div class="row-btns examples"><span class="note">Examples from the atlas itself:</span><button type="button" class="btn" data-data="volume">Volumes</button><button type="button" class="btn" data-data="asymmetry">Left–right asymmetry</button></div>
</div>
<h3>Points</h3>
<p class="pad note">Coordinates, one point a line: x, y, z, then a value and a name if you like. Each point is shown with what lies there.</p>
<div class="pad data-form">
<textarea id="data-points" spellcheck="false" placeholder="x, y, z, value, name&#10;-38, -22, 56, 3.1, left hand area&#10;0, -52, 26, 2.0, posterior cingulate">${esc(this.pointsText)}</textarea>
<div class="row-btns"><button type="button" class="btn primary" data-data="points">Show the points</button><button type="button" class="btn" data-data="clear-points">Clear</button>
<label class="inline">In <select id="data-space">${[['mni', 'MNI'], ['brett', 'Talairach (Brett)'], ['spm', 'Talairach (Lancaster, SPM)'], ['fsl', 'Talairach (Lancaster, FSL)']].map(([v, t]) => `<option value="${v}"${v === this.space ? ' selected' : ''}>${t}</option>`).join('')}</select></label></div>
<div class="report">${this.pointsReport}</div>
${this.points.length ? `<div class="region-list">${this.points.map((p, i) => `<button type="button" data-point="${i}" class="${i === this.chosenPoint ? 'selected' : ''}"><i style="background:#f6c47e"></i>${esc(p.name || `Point ${i + 1}`)} <small style="margin-left:auto;color:var(--text-faint)">${coords(p.mni)}</small></button>`).join('')}</div>` : ''}
</div>`;
  }

  bind(left) {
    const values = left.querySelector('#data-values'), points = left.querySelector('#data-points');
    values.addEventListener('input', () => { this.text = values.value; });
    points.addEventListener('input', () => { this.pointsText = points.value; });
    left.querySelector('#data-colormap').addEventListener('change', (e) => { this.colormap = e.target.value; if (this.values) { this.install(); } });
    left.querySelector('#data-space').addEventListener('change', (e) => { this.space = e.target.value; });
    left.addEventListener('click', (e) => {
      const b = e.target.closest('[data-data], [data-point]');
      if (!b) return;
      if (b.dataset.point) { this.choosePoint(Number(b.dataset.point)); this.ctx.renderLeft(); return; }
      const a = b.dataset.data;
      if (a === 'map') { this.mapValues(this.text, { title: 'Values', unit: '' }); this.ctx.renderLeft(); }
      if (a === 'clear') { this.values = null; this.text = ''; this.report = ''; this.install(); this.ctx.renderLeft(); }
      if (a === 'volume' || a === 'asymmetry') { if (a === 'volume') this.colormap = 'thermal'; this.builtIn(a); }
      if (a === 'points') this.showPoints();
      if (a === 'clear-points') { this.points = []; this.pointsText = ''; this.pointsReport = ''; this.chosenPoint = -1; this.markers(); this.ctx.renderLeft(); }
    });
  }
}
