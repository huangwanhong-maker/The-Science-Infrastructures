// The ocean beneath a chosen place: asked of the ocean model on demand (through
// the local server and the data library), summarised in the place card at the
// time shown, and drawn as a profiler's record (src/ui/column.js).

import { Requests } from '../data/requests.js';
import { ColumnFigure } from './column.js';
import { columnSampler, mixedLayerDepth, flowToward } from '../core/ocean.js';
import { compassPoint } from '../core/derived.js';

const DATASET = 'hycom-column';
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const clock = (s) => `${Math.floor(s / 60)}:${String(Math.floor(s % 60)).padStart(2, '0')}`;

export class OceanBeneath {
  constructor(host, { onChange }) {
    this.requests = new Requests();
    this.figure = new ColumnFigure(host);
    this.onChange = onChange;
    this.timer = 0;
  }

  // Places are asked for to 0.01 degree (about 1 km), finer than the model's grid.
  key(lon, lat) { return `${lat.toFixed(2)},${lon.toFixed(2)}`; }

  // The model covers the sea from 80 S northward; depth in metres, positive down.
  offered(lat, depth) { return Number.isFinite(depth) && depth > 5 && lat > -80; }

  ask(lon, lat) {
    this.requests.ask(this.key(lon, lat), DATASET, { lat: lat.toFixed(2), lon: lon.toFixed(2) }, () => this.changed());
    this.changed();
  }

  // While a request runs, the card is refreshed every second for its clock.
  changed() {
    this.onChange();
    const running = [...this.requests.jobs.values()].some((e) => !e.product && !e.error);
    if (running && !this.timer) this.timer = setInterval(() => this.changed(), 1000);
    if (!running && this.timer) { clearInterval(this.timer); this.timer = 0; }
  }

  show(lon, lat, t, seaFloor) {
    const entry = this.requests.get(this.key(lon, lat));
    if (entry && entry.product) this.figure.show(entry.product, { lon, lat, now: t, seaFloor });
  }

  // Rows for the place card: [term, html, kind].
  rows(lon, lat, t, depth) {
    if (!this.offered(lat, depth)) return [];
    const entry = this.requests.get(this.key(lon, lat));
    const button = (action, text) => `<button type="button" class="action" data-action="${action}">${esc(text)}</button>`;
    if (!entry) {
      return [['Beneath the surface', button('ocean-column', 'Currents and temperature by depth')],
        ['', 'From the US Navy’s ocean model (HYCOM), two days back and four ahead: about four minutes the first time.', 'gloss']];
    }
    if (entry.error) {
      return [['Beneath the surface', button('ocean-column', 'Ask again')], ['', esc(entry.error), 'gloss']];
    }
    if (!entry.product) {
      const job = entry.job || {};
      const elapsed = (Date.now() - entry.asked) / 1000;
      const expected = job.seconds ? ` <span class="dir">of about ${clock(job.seconds)}</span>` : '';
      return [['Beneath the surface', `${job.state === 'queued' ? 'waiting' : 'asking HYCOM'} ${clock(elapsed)}${expected}`]];
    }
    const p = entry.product;
    const times = p.times.map((s) => Date.parse(s)), depths = p.depths;
    const tt = Math.max(times[0], Math.min(times[times.length - 1], t));
    const field = (std) => {
      const v = p.variables.find((x) => x.standard_name === std);
      return v ? columnSampler(times, depths, v.values.map((row) => row.map((x) => (x === null ? NaN : x)))) : null;
    };
    const u = field('eastward_sea_water_velocity'), v = field('northward_sea_water_velocity'), T = field('sea_water_temperature');
    const rows = [];
    const current = (z) => {
      const a = u && u(tt, z), b = v && v(tt, z);
      if (!Number.isFinite(a) || !Number.isFinite(b)) return null;
      return `${Math.hypot(a, b).toFixed(2)} m s⁻¹ <span class="dir">toward ${compassPoint(flowToward(a, b))}</span>`;
    };
    const surface = current(0);
    if (surface) rows.push(['Surface current', surface]);
    const deep = current(500);
    if (deep) rows.push(['Current at 500 m', deep]);
    if (T) {
      const profile = depths.map((z) => T(tt, z));
      const m = mixedLayerDepth(profile, depths);
      if (Number.isFinite(m.depth)) rows.push(['Mixed layer', `${m.reachesBottom ? 'to the floor, ' : ''}${Math.round(m.depth)} m`]);
      const t500 = T(tt, 500);
      if (Number.isFinite(profile[0])) rows.push(['Water temperature', `${profile[0].toFixed(1)} °C${Number.isFinite(t500) ? ` <span class="dir">${t500.toFixed(1)} °C at 500 m</span>` : ''}`]);
    }
    if (tt !== t) rows.push(['', `Values for ${new Date(tt).toISOString().slice(0, 16).replace('T', ' ')} UTC, the nearest time the column covers.`, 'gloss']);
    rows.push(['By depth and time', button('ocean-show', 'Show the column')]);
    return rows;
  }
}
