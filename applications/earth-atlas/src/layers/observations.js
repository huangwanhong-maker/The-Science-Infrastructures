// Observations at places: weather stations (METAR), buoys (NDBC) and our own
// sensors (observation/sensors), from products of kind `observations`
// (data/README.md §4). Each station is a dot coloured by one variable: its
// report nearest the time shown, within a window, or none. A set of sensors is
// read again when it changes, so a reading posted to the API appears at once.

import { MarkerLayer } from '../render/markers.js';
import { colourAt, cssGradient } from '../render/colormaps.js';
import { ticks, position } from './surface.js';
import { displayUnit, toShown, format } from '../../../../visualization/js/index.js';

const HOUR = 3600000;
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

// Colours and ranges by variable, in the units shown, the same as the catalogue
// gives the quantity (src/layers/catalogue.js), so that a station's dot and the
// model's field beneath it compare by eye; others take their observed range.
export const OBSERVED = {
  air_temperature: { colormap: 'thermal', range: [-35, 40] },
  dew_point_temperature: { colormap: 'thermal', range: [-35, 30] },
  sea_surface_temperature: { colormap: 'thermal', range: [-2, 32] },
  relative_humidity: { colormap: 'tempo', reverse: true, range: [0, 100] },
  air_pressure_at_sea_level: { colormap: 'dense', reverse: true, range: [975, 1040] },
  wind_speed: { colormap: 'orchid', range: [0, 25] },
  wind_speed_of_gust: { colormap: 'orchid', range: [0, 35] },
  sea_surface_wave_significant_height: { colormap: 'deep', reverse: true, range: [0, 9] },
  visibility_in_air: { colormap: 'deep', range: [0, 10000] },
};

export class ObservationsLayer {
  constructor(stage, library) {
    this.stage = stage;
    this.library = library;
    this.markers = new MarkerLayer({ capacity: 8192, fill: 0.85, order: 5 });
    stage.overlay.add(this.markers.mesh);
    this.name = null;
    this.variable = null;
    this.data = null;
    this.shown = [];
    this.window = 3 * HOUR;          // a report counts within this of the time shown
    this.time = 0;
  }

  // Shows a set and one of its variables (the first if none is given); null hides.
  async select(name, variable = null) {
    this.name = name;
    if (!name) { this.data = null; this.shown = []; this.markers.clear(); return; }
    const data = await this.library.observationSet(name);
    if (this.name !== name) return;
    this.data = data;
    this.index = Object.fromEntries(data.fields.map((f, i) => [f, i]));
    const names = data.variables.map((v) => v.name);
    this.variable = names.includes(variable) ? variable : names[0];
    // The sensors' own window is wider: a home station may report every few minutes or every hour.
    this.window = /^sensor-/.test(name) ? 12 * HOUR : 3 * HOUR;
    this.style();
    this.update(this.time, true);
  }

  // Reads the set again (new readings), keeping what is shown.
  async refresh() { if (this.name) await this.select(this.name, this.variable); }

  // Readings as they arrive from the API's stream ({station, lon, lat, elevation, time, values}),
  // added to the rows held, so that a sensor's reading is drawn the moment it is posted.
  addReadings(readings) {
    if (!this.data) return;
    const names = this.data.variables.map((v) => v.name);
    for (const r of readings) {
      const time = typeof r.time === 'number' ? new Date(r.time).toISOString().replace('.000Z', 'Z') : r.time;
      const values = r.values || {};
      this.data.rows.push([r.station, r.name || r.station, r.lon, r.lat, r.elevation ?? null, time, ...names.map((n) => (Number.isFinite(values[n]) ? values[n] : null))]);
      const t = Date.parse(time);
      if (this.data.window && t > Date.parse(this.data.window.end)) this.data.window.end = time;
    }
    this.update(this.time, true);
  }

  setVariable(variable) {
    if (!this.data || !this.data.variables.some((v) => v.name === variable)) return;
    this.variable = variable;
    this.style();
    this.update(this.time, true);
  }

  // The colour scale of the variable shown, in the units shown.
  style() {
    const v = this.meta();
    this.unit = displayUnit(v.units);
    const known = OBSERVED[v.standard_name || v.name];
    let range = known && known.range;
    if (!range) {
      const values = this.data.rows.map((r) => r[this.index[v.name]]).filter((x) => x !== null && Number.isFinite(x)).map((x) => toShown(x, this.unit)).sort((a, b) => a - b);
      range = values.length ? [values[Math.floor(values.length * 0.02)], values[Math.floor(values.length * 0.98)]] : [0, 1];
      if (range[0] === range[1]) range = [range[0] - 1, range[1] + 1];
    }
    this.range = range;
    this.colormap = (known && known.colormap) || 'matter';
    this.reverse = Boolean(known && known.reverse);
  }

  meta() { return this.data.variables.find((v) => v.name === this.variable); }

  // Each station's report nearest the time shown, within the window.
  update(time, force = false) {
    if (!force && time === this.time) return;
    this.time = time;
    if (!this.data) return;
    const ix = this.index, col = ix[this.variable], best = new Map();
    for (const row of this.data.rows) {
      const value = row[col];
      if (value === null || !Number.isFinite(value)) continue;
      const dt = Math.abs(Date.parse(row[ix.time]) - time);
      if (dt > this.window) continue;
      const prior = best.get(row[ix.id]);
      if (!prior || dt < prior.dt) best.set(row[ix.id], { row, dt, value });
    }
    this.shown = [...best.values()];
    const proj = this.stage.projector(), p = [0, 0, 0], pos = [], colours = [], sizes = [];
    const [lo, hi] = this.range;
    for (const s of this.shown) {
      if (!proj.project(s.row[ix.lon], s.row[ix.lat], p, 1.0028)) continue;
      const t = Math.max(0, Math.min(1, (toShown(s.value, this.unit) - lo) / (hi - lo)));
      const [r, g, b] = colourAt(this.colormap, this.reverse ? 1 - t : t);
      pos.push(p[0], p[1], p[2]); colours.push(r / 255, g / 255, b / 255, 1); sizes.push(/^sensor-/.test(this.name) ? 5 : 3.2);
    }
    this.markers.set(new Float32Array(pos), new Float32Array(colours), new Float32Array(sizes));
  }

  frame(resolution) { this.markers.frame(resolution); }

  // The legend: as a surface's, with the count of stations reporting.
  legend() {
    if (!this.data) return null;
    const v = this.meta(), [lo, hi] = this.range;
    return {
      title: `${this.data.title}: ${(v.long_name || v.name).replace(/^./, (c) => c.toUpperCase())}`,
      unit: this.unit.unit,
      gradient: cssGradient(this.colormap, { reverse: this.reverse }),
      ticks: ticks(lo, hi).map((value) => ({ at: position(value, lo, hi, 'linear'), value })),
      note: `${this.shown.length} reporting within ${Math.round(this.window / HOUR)} h of the time shown · ${(this.data.source && this.data.source.credit) || ''}`,
    };
  }

  // The nearest station within 150 km of a place, for the place card: [term, html, kind].
  rows(lon, lat) {
    if (!this.data || !this.shown.length) return [];
    const R = Math.PI / 180, ix = this.index;
    let near = null, km = 150;
    for (const s of this.shown) {
      const a = s.row[ix.lat] * R, b = lat * R, dl = (s.row[ix.lon] - lon) * R;
      const d = 6371 * Math.acos(Math.min(1, Math.sin(a) * Math.sin(b) + Math.cos(a) * Math.cos(b) * Math.cos(dl)));
      if (d < km) { km = d; near = s; }
    }
    if (!near) return [];
    const row = near.row, out = [[`${row[ix.name] || row[ix.id]}`, `<span class="dim">${Math.round(km)} km away, ${esc(row[ix.time].slice(11, 16))} UTC</span>`, 'gloss']];
    for (const v of this.data.variables) {
      const x = row[ix[v.name]];
      if (x === null || !Number.isFinite(x)) continue;
      const u = displayUnit(v.units);
      out.push([(v.long_name || v.name).replace(/^./, (c) => c.toUpperCase()), `${format(toShown(x, u))} ${esc(u.unit)}`]);
    }
    return out;
  }

  // The time span of the set, for the timeline.
  windowSpan() {
    return this.data && this.data.window ? [Date.parse(this.data.window.start), Date.parse(this.data.window.end)] : null;
  }

  set visible(v) { this.markers.visible = v; }
}
