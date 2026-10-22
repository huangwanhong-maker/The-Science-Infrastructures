// Loading data products from the laboratory's data library (data/products/).
//
// The application reads only products, never raw downloads. The index lists
// what has been built; each product's metadata is small and read at once, and
// its values (the larger binary files) are read when a layer first needs them.

import { FieldProduct, ClassRaster } from '../core/grid.js';

async function json(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  return r.json();
}

async function binary(url, onProgress) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  const total = Number(r.headers.get('Content-Length')) || 0;
  if (!onProgress || !total || !r.body) return r.arrayBuffer();
  const reader = r.body.getReader();
  const out = new Uint8Array(total);
  let done = 0;
  for (;;) {
    const { value, done: end } = await reader.read();
    if (end) break;
    out.set(value, done);
    done += value.length;
    onProgress(done / total);
  }
  return out.buffer;
}

export class Library {
  constructor(base = 'data/products') {
    this.base = base;
    this.index = null;
    this.fieldMeta = new Map();
    this.fields = new Map();
    this.indicators = new Map();
    this.events = new Map();
    this.classMeta = new Map();
    this.classMaps = new Map();
  }

  async open() {
    this.index = await json(`${this.base}/index.json`);
    const fields = this.index.fields || [];
    const metas = await Promise.all(fields.map((e) => json(`${this.base}/${e.file}`).catch(() => null)));
    fields.forEach((e, i) => { if (metas[i]) this.fieldMeta.set(e.name, metas[i]); });
    const classes = this.index.classes || [];
    const cmetas = await Promise.all(classes.map((e) => json(`${this.base}/${e.file}`).catch(() => null)));
    classes.forEach((e, i) => { if (cmetas[i]) this.classMeta.set(e.name, cmetas[i]); });
    return this;
  }

  list(kind) { return (this.index && this.index[kind]) || []; }

  async geography(name) {
    const entry = this.list('geography').find((e) => !name || e.name === name);
    if (!entry) return null;
    const meta = await json(`${this.base}/${entry.file}`);
    const buffer = await binary(`${this.base}/${meta.coordinates.file}`);
    return { meta, coords: new Float32Array(buffer) };
  }

  // A field product with its values, read once.
  field(name, onProgress) {
    if (!this.fields.has(name)) {
      const meta = this.fieldMeta.get(name);
      this.fields.set(name, binary(`${this.base}/${meta.file}`, onProgress).then((b) => new FieldProduct(meta, b)));
    }
    return this.fields.get(name);
  }

  // A class product (a category per cell), read once.
  classes(name) {
    if (!this.classMaps.has(name)) {
      const meta = this.classMeta.get(name);
      this.classMaps.set(name, binary(`${this.base}/${meta.file}`).then((b) => new ClassRaster(meta, b)));
    }
    return this.classMaps.get(name);
  }

  // A field product computed in the application, offered like one from the library.
  addField(product) {
    this.fieldMeta.set(product.meta.name, product.meta);
    this.fields.set(product.meta.name, Promise.resolve(product));
  }

  // A plates product: the plates' outlines (rings of points) and their boundary steps, read once.
  plates(name) {
    if (!this.plateSets) this.plateSets = new Map();
    const entry = this.list('plates').find((e) => !name || e.name === name);
    if (!entry) return Promise.resolve(null);
    if (!this.plateSets.has(entry.name)) {
      this.plateSets.set(entry.name, json(`${this.base}/${entry.file}`)
        .then(async (meta) => ({ meta, coords: new Float32Array(await binary(`${this.base}/${meta.file}`)) })));
    }
    return this.plateSets.get(entry.name);
  }

  // A crust product (the layers of CRUST1.0): its description and its values, read once.
  crust() {
    if (!this.crustSet) {
      const entry = this.list('crust')[0];
      this.crustSet = entry ? json(`${this.base}/${entry.file}`)
        .then(async (meta) => ({ meta, values: new Int16Array(await binary(`${this.base}/${meta.file}`)) })) : Promise.resolve(null);
    }
    return this.crustSet;
  }

  // A stations product (GPS velocities), read once.
  stationSet(name) {
    if (!this.stationSets) this.stationSets = new Map();
    if (!this.stationSets.has(name)) {
      const entry = this.list('stations').find((e) => e.name === name);
      this.stationSets.set(name, json(`${this.base}/${entry.file}`));
    }
    return this.stationSets.get(name);
  }

  // A class product computed in the application, offered like one from the library.
  addClasses(raster) {
    this.classMeta.set(raster.meta.name, raster.meta);
    this.classMaps.set(raster.meta.name, Promise.resolve(raster));
  }

  indicator(name) {
    if (!this.indicators.has(name)) {
      const entry = this.list('indicators').find((e) => e.name === name);
      this.indicators.set(name, json(`${this.base}/${entry.file}`));
    }
    return this.indicators.get(name);
  }

  trackSet(name) {
    if (!this.trackSets) this.trackSets = new Map();
    if (!this.trackSets.has(name)) {
      const entry = this.list('tracks').find((e) => e.name === name);
      this.trackSets.set(name, json(`${this.base}/${entry.file}`));
    }
    return this.trackSets.get(name);
  }

  eventSet(name) {
    if (!this.events.has(name)) {
      const entry = this.list('events').find((e) => e.name === name);
      this.events.set(name, json(`${this.base}/${entry.file}`));
    }
    return this.events.get(name);
  }
}
