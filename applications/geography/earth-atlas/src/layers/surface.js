// The surface layer: one quantity over the whole surface, as a field, as the
// difference of two models' fields, or as a value per country, coloured by a
// colour map; or a class per cell (climate classes), alone or as the changes
// from another period.

import { lookupTable, colourAt, cssGradient } from '../render/colormaps.js';
import { FieldTextures, setField, setLut, setCountryColours, classTexture, setClasses } from '../render/surface.js';
import { classColours, soften, cssColour } from '../render/classes.js';
import { displayUnit, toShown, toStored } from '../core/units.js';
import { joinCountries, valueForYear, latestWellCoveredYear } from '../core/join.js';

// Rounded ends for an automatic range.
function niceRange(lo, hi) {
  const span = hi - lo || Math.abs(hi) || 1;
  const step = 10 ** Math.floor(Math.log10(span / 5));
  return [Math.floor(lo / step) * step, Math.ceil(hi / step) * step];
}

function percentile(sorted, p) {
  if (!sorted.length) return NaN;
  const x = (sorted.length - 1) * p, i = Math.floor(x);
  return sorted[i] + (sorted[Math.min(i + 1, sorted.length - 1)] - sorted[i]) * (x - i);
}

// Ticks for a legend: about n values, nicely rounded, within [lo, hi].
export function ticks(lo, hi, scale = 'linear', n = 5) {
  if (scale === 'log' && lo > 0) {
    const out = [];
    for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) {
      for (const m of [1, 2, 5]) {
        const v = m * 10 ** e;
        if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v);
      }
    }
    const decade = (v) => Math.abs(Math.log10(v) - Math.round(Math.log10(v))) < 1e-9;
    return out.length > 7 ? out.filter(decade) : out;
  }
  const raw = (hi - lo) / n;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw) || raw;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-6; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return out;
}

export function position(v, lo, hi, scale) {
  if (scale === 'log') return (Math.log(v) - Math.log(lo)) / (Math.log(hi) - Math.log(lo));
  if (scale === 'sqrt') return (Math.sqrt(Math.max(0, v)) - Math.sqrt(Math.max(0, lo))) / (Math.sqrt(hi) - Math.sqrt(Math.max(0, lo)));
  return (v - lo) / (hi - lo);
}

export class SurfaceLayer {
  constructor(stage, library, countries) {
    this.stage = stage;
    this.library = library;
    this.countries = countries;        // the geography product's countries
    this.textures = new FieldTextures();
    this.classTextures = new Map();
    this.def = null;
    this.product = null;
    this.info = null;
  }

  // def: the layer; compare: another layer with the same key (another model),
  // to show def minus compare.
  async select(def, compare = null, time = Date.now()) {
    this.def = def;
    this.compare = null;
    this.product = null;
    const u = this.stage.material.uniforms;
    if (!def) { u.uMode.value = 0; return; }
    if (def.kind === 'classes') {
      const other = compare && compare.kind === 'classes';
      const [a, b] = await Promise.all([this.library.classes(def.product), other ? this.library.classes(compare.product) : null]);
      if (this.def !== def) return;
      this.product = a;
      this.compare = b ? { def: compare, product: b } : null;
      setClasses(this.stage.material, this.classTexture(a), a, classColours(a.meta.classes),
        b ? { raster: b, texture: this.classTexture(b) } : null);
      return;
    }
    if (def.kind === 'field' && compare && compare.kind === 'field' && !def.variables[1]) {
      const [a, b] = await Promise.all([this.library.field(def.product), this.library.field(compare.product)]);
      if (this.def !== def) return;
      this.product = a;
      this.compare = { def: compare, product: b };
      const unit = displayUnit(def.units);
      // A difference: the offset of the units cancels.
      this.unit = def.percent ? { unit: '%', a: 100, b: 0 } : { unit: unit.unit, a: unit.a, b: 0 };
      const r = this.differenceRange(time);
      this.range = [-r, r];
      this.scale = 'linear';
      setLut(this.stage.material, lookupTable('balance'));
      u.uRange.value.set(-r / this.unit.a, r / this.unit.a, 0);
      return;
    }
    if (def.kind === 'field') {
      const product = await this.library.field(def.product);
      if (this.def !== def) return;
      this.product = product;
      this.unit = displayUnit(def.units);
      // A spread, like a difference, is a size: the offset of the units does not apply.
      if (def.spread) this.unit = { ...this.unit, b: 0 };
      const st = def.style;
      let range = st.range;
      if (!range) range = niceRange(toShown(st.auto[0], this.unit), toShown(st.auto[1], this.unit));
      if (def.percent) this.unit = { unit: '%', a: 100, b: 0 };
      this.range = def.percent ? [range[0] * 100, range[1] * 100] : range;
      this.scale = st.scale || 'linear';
      setLut(this.stage.material, lookupTable(st.colormap, { reverse: st.reverse, clip: st.clip, fade: st.fade }));
      const stored = this.range.map((x) => toStored(x, this.unit));
      u.uRange.value.set(Math.min(...stored), Math.max(...stored), this.scale === 'log' ? 1 : this.scale === 'sqrt' ? 2 : 0);
    } else {
      const product = await this.library.indicator(def.product);
      if (this.def !== def) return;
      this.product = product;
      this.join = joinCountries(this.countries, product.economies);
      this.unit = { unit: def.units || '', a: 1, b: 0 };
      // One range for all years, so that years can be compared.
      const all = [];
      for (const e of Object.values(product.economies)) for (const v of e.values) if (v !== null) all.push(v);
      all.sort((a, b) => a - b);
      const lo = percentile(all, 0.02), hi = percentile(all, 0.98);
      this.scale = lo > 0 && hi / lo > 30 ? 'log' : 'linear';
      this.range = this.scale === 'log'
        ? [10 ** Math.floor(Math.log10(lo)), 10 ** Math.ceil(Math.log10(hi))]
        : niceRange(lo, hi);
      if (/^%/.test(def.units || '')) this.range = [Math.max(0, this.range[0]), Math.min(100, this.range[1])];
      this.defaultYear = latestWellCoveredYear(product.economies, product.years);
      this.year = this.year && product.years.includes(this.year) ? this.year : this.defaultYear;
      this.paintCountries();
    }
  }

  // A symmetric range for a difference: the 98th percentile of its size over a
  // lattice of places and over the times the two products share, rounded up,
  // so that the scale holds while time runs.
  differenceRange(time = Date.now()) {
    const a = this.product, b = this.compare.product, va = this.def.variables[0], vb = this.compare.def.variables[0];
    const lo = Math.max(a.times[0] ?? time, b.times[0] ?? time);
    const hi = Math.min(a.times[a.times.length - 1] ?? time, b.times[b.times.length - 1] ?? time);
    const times = hi > lo ? [0, 0.25, 0.5, 0.75, 1].map((f) => lo + f * (hi - lo)) : [time];
    const sizes = [];
    for (const t of times) for (let lat = -80; lat <= 80; lat += 5) for (let lon = -180; lon < 180; lon += 5) {
      const d = a.sample(va, t, lon, lat) - b.sample(vb, t, lon, lat);
      if (Number.isFinite(d)) sizes.push(Math.abs(d) * this.unit.a);
    }
    sizes.sort((x, y) => x - y);
    const p = sizes.length ? sizes[Math.floor(0.98 * (sizes.length - 1))] : 1;
    const mag = 10 ** Math.floor(Math.log10(p || 1));
    return [1, 2, 2.5, 5, 10].map((m) => m * mag).find((x) => x >= p) || p;
  }

  paintCountries() {
    const p = this.product, st = this.def.style;
    const colours = new Uint8Array(4 * (this.countries.length + 1));
    this.values = new Array(this.countries.length).fill(null);
    this.countries.forEach((c, i) => {
      const code = this.join.codes[i];
      const r = code && valueForYear(p.economies[code], p.years, this.year, 4);
      const k = 4 * (i + 1);
      if (!r) {
        colours.set([96, 84, 104, code ? 150 : 90], k);       // no value: a quiet grey-plum
        return;
      }
      this.values[i] = { ...r, code, name: p.economies[code].name };
      let t = Math.max(0, Math.min(1, position(r.value, this.range[0], this.range[1], this.scale)));
      t = st.clip[0] + (st.clip[1] - st.clip[0]) * (st.reverse ? 1 - t : t);
      const [cr, cg, cb] = colourAt(st.colormap, t);
      colours.set([cr, cg, cb, 236], k);
    });
    setCountryColours(this.stage.material, colours);
    this.stage.material.uniforms.uMode.value = 3;
  }

  setYear(year) {
    if (!this.def || this.def.kind !== 'indicator') return;
    this.year = year;
    this.paintCountries();
  }

  // An integer texture for a class map, kept for a few maps.
  classTexture(raster) {
    const name = raster.meta.name;
    let t = this.classTextures.get(name);
    if (!t) {
      t = classTexture(raster);
      this.classTextures.set(name, t);
      if (this.classTextures.size > 4) {
        const [k, old] = this.classTextures.entries().next().value;
        if (old !== this.stage.material.uniforms.uK0.value && old !== this.stage.material.uniforms.uK1.value) {
          this.classTextures.delete(k);
          old.dispose();
        }
      }
    }
    return t;
  }

  samePeriod(other) {
    return other.period && this.def.period && other.period[0] === this.def.period[0] && other.period[1] === this.def.period[1];
  }

  // The class at a place, and the other map's class there when comparing.
  classesAt(lon, lat) {
    if (!this.def || this.def.kind !== 'classes' || !this.product) return null;
    return { now: this.product.classAt(lon, lat), then: this.compare ? this.compare.product.classAt(lon, lat) : null };
  }

  // Points the GPU at the data for an instant; returns what is shown.
  update(time) {
    if (!this.def || !this.product) return null;
    if (this.def.kind === 'indicator') return { year: this.year };
    if (this.def.kind === 'classes') return null;
    const b = this.product.bracket(time, this.def.variables[0]);
    const second = this.compare && {
      product: this.compare.product, name: this.compare.def.variables[0],
      bracket: this.compare.product.bracket(time, this.compare.def.variables[0]),
    };
    setField(this.stage.material, this.textures, this.product, this.def.variables, b, second);
    this.shown = b;
    return b;
  }

  // The value at a place (stored units), or NaN.
  sample(lon, lat, time) {
    if (!this.product || !this.def || this.def.kind !== 'field') return NaN;
    const [a, b] = this.def.variables;
    const x = this.product.sample(a, time, lon, lat);
    if (this.compare) return x - this.compare.product.sample(this.compare.def.variables[0], time, lon, lat);
    return b ? Math.hypot(x, this.product.sample(b, time, lon, lat)) : x;
  }

  legend() {
    if (!this.def || !this.product) return null;
    if (this.def.kind === 'classes') {
      const c = this.compare;
      return {
        kind: 'classes',
        title: c ? `${this.def.label}: ${this.def.model}, where it differs from ${this.samePeriod(c.def) ? c.def.sourceName : c.def.model}` : `${this.def.label}, ${this.def.model}`,
        unit: '',
        classes: this.product.meta.classes.map((k) => ({ code: k.code, name: k.name, colour: cssColour(soften(k.colour)) })),
      };
    }
    const st = this.compare ? { colormap: 'balance' } : this.def.style;
    return {
      title: this.compare ? `${this.def.label}: ${this.def.model} minus ${this.compare.def.model}` : this.def.label,
      unit: this.unit.unit,
      gradient: cssGradient(st.colormap, { reverse: st.reverse, clip: st.clip, fade: st.fade }),
      ticks: ticks(this.range[0], this.range[1], this.scale).map((v) => ({ at: position(v, this.range[0], this.range[1], this.scale), value: v })),
      range: this.range,
      scale: this.scale,
    };
  }
}
