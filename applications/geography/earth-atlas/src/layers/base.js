// The base map: relief, coastlines, boundaries, rivers, the graticule and the
// names of places, from the geography product and the relief field.

import { SegmentLayer, polylineSegments } from '../render/segments.js';
import { baseImage } from '../render/relief.js';
import { setBase, setCountryRaster } from '../render/surface.js';
import { indexRaster } from '../core/raster.js';

const RASTER = [4096, 2048];      // countries and seas, for looking up places
const BASE = [8192, 4096];        // the painted base map: about 5 km to a pixel
const OBLIQUITY = 23.4362;          // degrees, IAU 2006 value for 2026 (23.4393 - 0.013 T)

export class BaseLayer {
  constructor(stage, geography) {
    this.stage = stage;
    this.geo = geography.meta;
    this.coords = geography.coords;
    const L = this.geo.layers;
    this.countries = L.countries;
    this.layers = {
      graticule: new SegmentLayer({ color: '#efe3ee', width: 0.8, opacity: 0.13, order: 1 }),
      parallels: new SegmentLayer({ color: '#CA8DBE', width: 1.0, opacity: 0.45, dash: 0.012, order: 1 }),
      rivers: new SegmentLayer({ color: '#a8c0e0', width: 0.9, opacity: 0.38, order: 1 }),
      coast: new SegmentLayer({ color: '#f6ece8', width: 1.0, opacity: 0.55, order: 2 }),
      borders: new SegmentLayer({ color: '#f6ece8', width: 0.8, opacity: 0.32, order: 2 }),
      disputed: new SegmentLayer({ color: '#f6ece8', width: 0.8, opacity: 0.32, dash: 0.006, order: 2 }),
      outline: new SegmentLayer({ color: '#cdb6c8', width: 1.2, opacity: 0.5, order: 2 }),
    };
    for (const l of Object.values(this.layers)) stage.overlay.add(l.mesh);
    const [W, H] = RASTER;
    this.countryRaster = indexRaster(this.coords, L.countries, W, H);
    this.marineRaster = indexRaster(this.coords, L.marine, W / 2, H / 2);
    setCountryRaster(stage.material, this.countryRaster, W, H);
    this.rebuild();
  }

  // Paints the base map in a worker (src/workers/relief.worker.js), or here if
  // workers are not available; resolves when the map is drawn.
  paint(elevation, style) {
    const [W, H] = BASE;
    const ticket = (this.painting = {});
    const L = this.geo.layers;
    const job = {
      layers: { land: L.land, lakes: L.lakes, glaciers: L.glaciers, ice_shelves: L.ice_shelves },
      coords: this.coords.slice(),
      elevation: elevation ? { meta: elevation.meta, codes: elevation.codes.slice() } : null,
      width: W, height: H, style,
    };
    const done = (rgba) => { if (this.painting === ticket) setBase(this.stage.material, rgba, W, H); };
    // Here, on the main thread: when workers are not available, or one fails.
    const here = () => {
      const { data } = baseImage({ elevation, geography: this.geo, coords: this.coords, width: W, height: H, style });
      done(new Uint8Array(data.buffer));
    };
    return new Promise((resolve) => {
      let worker;
      try {
        worker = new Worker(new URL('../workers/relief.worker.js', import.meta.url), { type: 'module' });
      } catch {
        here();
        resolve();
        return;
      }
      worker.onmessage = (e) => {
        worker.terminate();
        if (e.data.rgba) done(e.data.rgba);
        else { console.error('relief:', e.data.error); here(); }
        resolve();
      };
      worker.onerror = (e) => {
        worker.terminate();
        console.error('relief worker:', e.message || 'failed');
        here();
        resolve();
      };
      worker.postMessage(job, [job.coords.buffer, ...(job.elevation ? [job.elevation.codes.buffer] : [])]);
    });
  }

  // Lines in world coordinates for the current projection.
  rebuild() {
    const proj = this.stage.projector();
    const L = this.geo.layers;
    const runs = (features, key = 'lines', test = () => true) => features.filter(test).flatMap((f) => f[key]);
    const set = (layer, segs) => layer.set(segs.a, segs.b, segs.info);
    const lift = (r) => (lon, lat, out) => proj.project(lon, lat, out, r);
    set(this.layers.coast, polylineSegments(this.coords, runs(L.coastline), lift(1.0010), proj.breaks));
    set(this.layers.borders, polylineSegments(this.coords, runs(L.boundaries, 'lines', (f) => f.class === 'international'), lift(1.0011), proj.breaks));
    set(this.layers.disputed, polylineSegments(this.coords, runs(L.boundaries, 'lines', (f) => f.class !== 'international'), lift(1.0011), proj.breaks));
    set(this.layers.rivers, polylineSegments(this.coords, runs(L.rivers), lift(1.0008), proj.breaks));

    // Meridians and parallels every 15 degrees; the tropics and polar circles,
    // where the Sun can stand at the zenith or stay up all day, dashed.
    const pts = [], grat = [];
    const line = (coords) => { grat.push([pts.length / 2, coords.length / 2]); pts.push(...coords); };
    for (let lon = -180; lon < 180; lon += 15) {
      const c = [];
      for (let lat = -90; lat <= 90; lat += 1) c.push(lon, lat);
      line(c);
    }
    for (let lat = -75; lat <= 75; lat += 15) {
      const c = [];
      for (let lon = -180; lon <= 180; lon += 1) c.push(lon, lat);
      line(c);
    }
    const g = new Float32Array(pts);
    set(this.layers.graticule, polylineSegments(g, grat, lift(1.0005), proj.breaks));
    const special = [], sruns = [];
    for (const lat of [OBLIQUITY, -OBLIQUITY, 90 - OBLIQUITY, -(90 - OBLIQUITY)]) {
      sruns.push([special.length / 2, 361]);
      for (let lon = -180; lon <= 180; lon += 1) special.push(lon, lat);
    }
    set(this.layers.parallels, polylineSegments(new Float32Array(special), sruns, lift(1.0006), proj.breaks));

    const outline = this.stage.mapOutline();
    if (outline) {
      const a = [], b = [], info = [];
      for (let i = 0; i + 1 < outline.length; i++) {
        a.push(outline[i][0], outline[i][1], 0.003); b.push(outline[i + 1][0], outline[i + 1][1], 0.003); info.push(0, 1, 0, 0);
      }
      this.layers.outline.set(new Float32Array(a), new Float32Array(b), new Float32Array(info));
    } else {
      this.layers.outline.clear();
    }
  }

  show(options) {
    const l = this.layers;
    l.coast.visible = options.coast;
    l.borders.visible = l.disputed.visible = options.borders;
    l.rivers.visible = options.rivers;
    l.graticule.visible = l.parallels.visible = options.graticule;
    l.outline.visible = this.stage.mode !== 'globe';
  }

  // The country (index into the geography's countries) or sea at a place.
  countryAt(lon, lat) {
    const [W, H] = RASTER;
    const i = Math.min(W - 1, Math.max(0, Math.floor(((lon + 180) / 360) * W)));
    const j = Math.min(H - 1, Math.max(0, Math.floor(((90 - lat) / 180) * H)));
    const id = this.countryRaster[j * W + i];
    return id ? id - 1 : -1;
  }

  seaAt(lon, lat) {
    const W = RASTER[0] / 2, H = RASTER[1] / 2;
    const i = Math.min(W - 1, Math.max(0, Math.floor(((lon + 180) / 360) * W)));
    const j = Math.min(H - 1, Math.max(0, Math.floor(((90 - lat) / 180) * H)));
    const id = this.marineRaster[j * W + i];
    return id ? this.geo.layers.marine[id - 1] : null;
  }

  // Names to place, in order of priority, for a zoom level and a language.
  labelItems(zoom, lang) {
    const L = this.geo.layers, items = [];
    // In English, Natural Earth's short name (China, not People's Republic of China).
    const name = (names) => (lang === 'en' ? names.default : names[lang]) || names.en || names.default;
    for (const m of L.marine) {
      if (m.class === 'ocean' && zoom < 6) {
        items.push({ key: `m${m.names.default}`, text: name(m.names), lon: m.label[0], lat: m.label[1], className: 'ocean', font: 'italic 13px "Source Serif 4", serif', tracking: 2.9, size: 14 });
      }
    }
    const upper = /^(en|fr|de|es|pt|it|nl|pl|sv|ru|uk|el|tr|id|vi|hu)$/.test(lang);
    for (const c of [...L.countries].sort((a, b) => a.min_label - b.min_label || a.label_rank - b.label_rank)) {
      const sovereign = c.type === 'Sovereign country' || c.type === 'Country';
      if (c.min_label > zoom + (sovereign ? 0.7 : -0.6)) continue;
      items.push({ key: `c${c.id}`, text: name(c.names), lon: c.label[0], lat: c.label[1], className: 'country', font: '600 11px "Source Serif 4", serif', upper, tracking: 1.3, size: 12 });
    }
    for (const m of L.marine) {
      if (m.class === 'ocean' || m.min_label > zoom + 0.6) continue;
      items.push({ key: `m${m.names.default}`, text: name(m.names), lon: m.label[0], lat: m.label[1], className: 'sea', font: 'italic 11.5px "Source Serif 4", serif', size: 12 });
    }
    const rankLimit = zoom < 2.8 ? -1 : zoom < 3.4 ? 0 : zoom < 4 ? 1 : zoom < 4.6 ? 2 : zoom < 5.4 ? 4 : 10;
    for (const p of L.places) {
      if (p.rank > rankLimit && !(p.capital && zoom >= 3.6)) continue;
      items.push({ key: `p${p.names.default}${p.lon}`, text: name(p.names), lon: p.lon, lat: p.lat, className: p.capital ? 'place capital' : 'place', font: '12px "Source Serif 4", serif', anchor: 'right', size: 13 });
    }
    return items;
  }

  frame(resolution) { for (const l of Object.values(this.layers)) l.frame(resolution, 0); }
}
