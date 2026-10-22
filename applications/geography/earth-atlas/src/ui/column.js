// The ocean beneath a place, drawn as a current profiler's record: the
// eastward and northward current and the temperature against depth and time.
// Along the top runs the height of the Sun at the place (day light, night
// dark); the present and the start of the forecast are marked; isotherms and
// the depth of the mixed layer are drawn on the temperature; the sea floor is
// grey. Pointing at the record reads every value there.

import { lookupTable } from '../render/colormaps.js';
import { columnSampler, symmetricRange, mixedLayerDepth, flowToward, inertialPeriodHours } from '../core/ocean.js';
import { isolines } from '../core/contour.js';
import { subsolarPoint, sunAltAz } from '../core/sun.js';
import { compassPoint } from '../core/derived.js';
import { arcDeg } from '../core/cyclones.js';

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const DAY = 86400000, HOUR = 3600000;
const FONT = '"Source Serif 4", Georgia, serif';
const INK = '#f6eef5', DIM = '#cdb6c8', FAINT = '#9a8496', GOLD = '#f6c47e', FLOOR = '#2a2229';

// The panels, recognised by standard name.
const PANELS = [
  { standard: 'eastward_sea_water_velocity', title: 'Eastward current', note: 'red flows east, blue west', map: 'balance', signed: true },
  { standard: 'northward_sea_water_velocity', title: 'Northward current', note: 'red flows north, blue south', map: 'balance', signed: true },
  { standard: 'sea_water_temperature', title: 'Temperature', note: 'isotherms; the mixed layer dashed', map: 'thermal' },
];

const fmtDay = (t) => new Date(t).toLocaleDateString('en-GB', { weekday: 'short', day: 'numeric', month: 'short', timeZone: 'UTC' });
const fmtTime = (t) => `${fmtDay(t)} ${new Date(t).toISOString().slice(11, 16)} UTC`;
const coords = (lon, lat) => `${Math.abs(lat).toFixed(2)}° ${lat >= 0 ? 'N' : 'S'}, ${Math.abs(lon).toFixed(2)}° ${lon >= 0 ? 'E' : 'W'}`;

// Round numbers for an axis from lo to hi, about n of them.
function ticks(lo, hi, n = 5) {
  const raw = (hi - lo) / n, p = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * p).find((s) => s >= raw) || 10 * p;
  const out = [];
  for (let x = Math.ceil(lo / step - 1e-9) * step; x <= hi + 1e-9; x += step) out.push(+x.toFixed(10));
  return out;
}

// The colour of the sky strip for the Sun at altitude h (degrees): night plum
// through twilight to a pale day.
function skyColour(h) {
  const s = Math.max(0, Math.min(1, (h + 18) / 24));
  const e = s * s * (3 - 2 * s);
  const night = [36, 22, 42], day = [243, 222, 186];
  return night.map((c, i) => Math.round(c + (day[i] - c) * e));
}

export class ColumnFigure {
  constructor(host, { onClose } = {}) {
    this.host = host;
    this.onClose = onClose;
    this.depthMode = 'upper';
    host.innerHTML = `
      <div class="figure-head">
        <div class="figure-titles"><h3></h3><p class="sub"></p></div>
        <div class="figure-controls">
          <div class="segmented" role="group" aria-label="Depth shown">
            <button type="button" data-depth="upper">0–1000 m</button>
            <button type="button" data-depth="all">To the sea floor</button>
          </div>
          <button type="button" class="close" aria-label="Close">✕</button>
        </div>
      </div>
      <div class="figure-body"><canvas class="plot"></canvas><canvas class="overlay"></canvas><div class="readout" hidden></div></div>
      <p class="caption"></p>`;
    this.body = host.querySelector('.figure-body');
    this.plot = host.querySelector('.plot');
    this.overlay = host.querySelector('.overlay');
    this.readout = host.querySelector('.readout');
    host.querySelectorAll('[data-depth]').forEach((b) => b.addEventListener('click', () => { this.depthMode = b.dataset.depth; this.draw(); }));
    host.querySelector('.close').addEventListener('click', () => { this.hide(); if (this.onClose) this.onClose(); });
    this.overlay.addEventListener('pointermove', (e) => this.point(e));
    this.overlay.addEventListener('pointerleave', () => { this.cursor = null; this.drawOverlay(); this.readout.hidden = true; });
    new ResizeObserver(() => { if (!this.host.hidden && this.product) this.draw(); }).observe(this.body);
  }

  get open() { return !this.host.hidden && Boolean(this.product); }

  // product: a columns/1 product; info: { lon, lat, now, seaFloor (m, positive down, or NaN) }
  show(product, info) {
    this.product = product;
    this.info = info;
    this.prepare();
    this.host.hidden = false;
    this.draw();
  }

  hide() {
    this.host.hidden = true;
    this.product = null;
    this.readout.hidden = true;
  }

  setNow(t) {
    if (!this.open) return;
    this.info.now = t;
    this.drawOverlay();
  }

  prepare() {
    const p = this.product;
    this.times = p.times.map((s) => Date.parse(s));
    this.depths = p.depths;
    this.latestRun = p.latest_run ? Date.parse(p.latest_run) : null;
    this.panels = [];
    for (const def of PANELS) {
      const v = p.variables.find((x) => x.standard_name === def.standard);
      if (!v) continue;
      const values = v.values.map((row) => row.map((x) => (x === null ? NaN : x)));
      this.panels.push({ ...def, variable: v, values, sample: columnSampler(this.times, this.depths, values) });
    }
    const temp = this.panels.find((d) => d.standard === 'sea_water_temperature');
    this.mld = temp ? temp.values.map((profile) => mixedLayerDepth(profile, this.depths)) : null;
    const k = p.deepest_level;
    this.modelFloor = k >= 0 ? (k + 1 < this.depths.length ? (this.depths[k] + this.depths[k + 1]) / 2 : this.depths[k]) : 0;
    const g = p.place;
    this.gridLon = g.grid_longitude;
    this.gridLat = g.grid_latitude;
    this.offsetKm = arcDeg(this.info.lon, this.info.lat, this.gridLon, this.gridLat) * 111.195;
    const head = this.host.querySelector('h3');
    head.textContent = `The ocean beneath ${coords(this.gridLon, this.gridLat)}`;
    const stepH = this.times.length > 1 ? Math.round((this.times[1] - this.times[0]) / HOUR) : 0;
    this.host.querySelector('.sub').textContent = [
      `${p.model || 'Ocean model'}, at its grid point ${this.offsetKm < 1 ? 'at' : `${this.offsetKm.toFixed(0)} km from`} the place chosen`,
      this.latestRun ? `its latest run ${fmtTime(this.latestRun)}; forecast after that` : '',
      stepH ? `every ${stepH} h` : '',
    ].filter(Boolean).join(' · ');
    const inertial = inertialPeriodHours(this.gridLat);
    const floor = Number.isFinite(this.info.seaFloor) && this.info.seaFloor > 0
      ? `Sea floor about ${Math.round(this.info.seaFloor).toLocaleString('en-GB')} m (ETOPO 2022); the model's deepest water here ${Math.round(this.depths[Math.max(0, k)]).toLocaleString('en-GB')} m.` : '';
    this.host.querySelector('.caption').textContent = [
      Number.isFinite(inertial) && Math.abs(this.gridLat) > 2
        ? `Inertial period here ${inertial.toFixed(1)} h: a current set going and then left alone turns full circle in that time (2π/f), which is why the currents sway.`
        : 'Near the equator the Coriolis force vanishes, and currents do not turn in inertial circles.',
      'A current profiler on a ship or a mooring records such a column directly, and its echoes show the plankton that rise at dusk and sink at dawn; a model has currents and temperature, but no plankton.',
      floor,
      `Source: ${p.source.credit}; mixed layer where the temperature first differs from that at 10 m by 0.2 °C (de Boyer Montégut et al. 2004).`,
    ].filter(Boolean).join(' ');
  }

  // Panel rectangles and scales for the body's size.
  layout() {
    const W = this.body.clientWidth, H = this.body.clientHeight;
    const left = 62, right = 104, top = 30, bottom = 30, gap = 12;
    const n = this.panels.length;
    const ph = Math.max(40, (H - top - bottom - gap * (n - 1)) / n);
    // The whole depth: down to the sea floor (the relief's or the model's, the
    // deeper), rounded up to a round number.
    const deepest = Math.max(this.modelFloor, Number.isFinite(this.info.seaFloor) ? this.info.seaFloor : 0, 50);
    const step = ticks(0, deepest, 5);
    const zMax = this.depthMode === 'upper' ? 1000 : Math.ceil((deepest * 1.02) / (step[1] - step[0])) * (step[1] - step[0]);
    const t0 = this.times[0], t1 = this.times[this.times.length - 1];
    const rects = this.panels.map((d, i) => ({ x: left, y: top + i * (ph + gap), w: W - left - right, h: ph }));
    return { W, H, left, right, top, bottom, rects, zMax, t0, t1,
      x: (t) => left + ((t - t0) / (t1 - t0)) * (W - left - right),
      t: (x) => t0 + ((x - left) / (W - left - right)) * (t1 - t0),
      y: (r, z) => r.y + (z / zMax) * r.h,
      z: (r, y) => ((y - r.y) / r.h) * zMax };
  }

  draw() {
    if (!this.product) return;
    // A new scale: the pointer's reading no longer applies.
    this.cursor = null;
    this.readout.hidden = true;
    this.host.querySelectorAll('[data-depth]').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.depth === this.depthMode)));
    const dpr = window.devicePixelRatio || 1;
    const L = (this.L = this.layout());
    for (const c of [this.plot, this.overlay]) {
      c.width = Math.round(L.W * dpr); c.height = Math.round(L.H * dpr);
      c.style.width = `${L.W}px`; c.style.height = `${L.H}px`;
    }
    const ctx = this.plot.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, L.W, L.H);
    // Ranges: currents symmetric about zero; temperature over the depths shown.
    for (const d of this.panels) {
      const shown = [];
      d.values.forEach((row) => row.forEach((v, k) => { if (this.depths[k] <= L.zMax) shown.push(v); }));
      if (d.signed) { const r = symmetricRange(shown, { step: 0.1, floor: 0.2 }); d.range = [-r, r]; }
      else {
        const ok = shown.filter(Number.isFinite);
        d.range = ok.length ? [Math.floor(Math.min(...ok)), Math.ceil(Math.max(...ok))] : [0, 1];
        if (d.range[1] - d.range[0] < 2) d.range[1] = d.range[0] + 2;
      }
      d.lut = lookupTable(d.map);
    }
    this.panels.forEach((d, i) => this.drawPanel(ctx, d, L.rects[i], L, dpr));
    this.drawSky(ctx, L);
    this.drawTimeAxis(ctx, L);
    this.drawOverlay();
  }

  drawPanel(ctx, d, r, L, dpr) {
    // The field, sampled at every device pixel: linear in time and depth.
    const w = Math.max(1, Math.round(r.w * dpr)), h = Math.max(1, Math.round(r.h * dpr));
    const image = ctx.createImageData(w, h);
    const [lo, hi] = d.range;
    const floorRGB = [42, 34, 41];
    const ts = new Float64Array(w), zs = new Float64Array(h);
    for (let i = 0; i < w; i++) ts[i] = L.t0 + ((i + 0.5) / w) * (L.t1 - L.t0);
    for (let j = 0; j < h; j++) zs[j] = ((j + 0.5) / h) * L.zMax;
    for (let j = 0; j < h; j++) {
      for (let i = 0; i < w; i++) {
        const v = d.sample(ts[i], zs[j]);
        const o = 4 * (j * w + i);
        if (Number.isFinite(v)) {
          const k = 4 * Math.max(0, Math.min(255, Math.round(((v - lo) / (hi - lo)) * 255)));
          image.data[o] = d.lut[k]; image.data[o + 1] = d.lut[k + 1]; image.data[o + 2] = d.lut[k + 2];
        } else {
          image.data[o] = floorRGB[0]; image.data[o + 1] = floorRGB[1]; image.data[o + 2] = floorRGB[2];
        }
        image.data[o + 3] = 255;
      }
    }
    ctx.putImageData(image, Math.round(r.x * dpr), Math.round(r.y * dpr));
    ctx.save();
    ctx.beginPath(); ctx.rect(r.x, r.y, r.w, r.h); ctx.clip();
    // The sea floor from the relief, where it lies within the depths shown.
    const floor = this.info.seaFloor;
    if (Number.isFinite(floor) && floor > 0 && floor < L.zMax) {
      ctx.strokeStyle = 'rgba(205, 182, 200, 0.75)'; ctx.setLineDash([5, 4]); ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(r.x, L.y(r, floor)); ctx.lineTo(r.x + r.w, L.y(r, floor)); ctx.stroke();
      ctx.setLineDash([]);
    }
    // Midnights, faintly.
    ctx.strokeStyle = 'rgba(246, 238, 245, 0.10)'; ctx.lineWidth = 1;
    for (let t = Math.ceil(L.t0 / DAY) * DAY; t <= L.t1; t += DAY) {
      ctx.beginPath(); ctx.moveTo(Math.round(L.x(t)) + 0.5, r.y); ctx.lineTo(Math.round(L.x(t)) + 0.5, r.y + r.h); ctx.stroke();
    }
    if (d.standard === 'sea_water_temperature') this.drawTemperatureLines(ctx, d, r, L);
    ctx.restore();
    // Frame, title, depth axis and colour bar.
    ctx.strokeStyle = 'rgba(240, 226, 238, 0.22)'; ctx.lineWidth = 1;
    ctx.strokeRect(r.x + 0.5, r.y + 0.5, r.w - 1, r.h - 1);
    // The title on a small dark label, legible over any colour.
    ctx.textBaseline = 'top'; ctx.textAlign = 'left';
    ctx.font = `600 13px ${FONT}`;
    const tw = ctx.measureText(d.title).width;
    ctx.font = `italic 12px ${FONT}`;
    const nw = ctx.measureText(`· ${d.note}`).width;
    ctx.fillStyle = 'rgba(20, 14, 22, 0.78)';
    ctx.beginPath(); ctx.roundRect(r.x + 5, r.y + 5, tw + nw + 20, 21, 3); ctx.fill();
    ctx.font = `600 13px ${FONT}`; ctx.fillStyle = INK; ctx.fillText(d.title, r.x + 11, r.y + 8);
    ctx.font = `italic 12px ${FONT}`; ctx.fillStyle = DIM;
    ctx.fillText(`· ${d.note}`, r.x + 16 + tw, r.y + 9);
    ctx.font = `11.5px ${FONT}`; ctx.fillStyle = FAINT; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const z of ticks(0, L.zMax, Math.max(2, Math.round(r.h / 45)))) {
      const y = L.y(r, z);
      ctx.fillText(z === 0 ? '0' : z.toLocaleString('en-GB'), r.x - 8, Math.min(r.y + r.h - 6, Math.max(r.y + 6, y)));
      ctx.fillRect(r.x - 4, Math.round(y), 3, 1);
    }
    if (d === this.panels[0]) {
      ctx.save(); ctx.translate(14, r.y + r.h / 2); ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center'; ctx.font = `italic 12px ${FONT}`; ctx.fillStyle = FAINT;
      ctx.fillText('depth (m)', 0, 0); ctx.restore();
    }
    this.drawColourBar(ctx, d, r);
  }

  drawTemperatureLines(ctx, d, r, L) {
    const nt = this.times.length, nz = this.depths.length;
    const grid = { nlon: nt, nlat: nz, lon0: 0, dlon: 1, lat0: 0, dlat: -1, periodic: false };
    const values = new Float32Array(nt * nz);
    for (let i = 0; i < nt; i++) for (let k = 0; k < nz; k++) values[k * nt + i] = d.values[i][k];
    // Fractional indices to time and depth.
    const tAt = (x) => { const i = Math.min(nt - 2, Math.floor(x)); return this.times[i] + (this.times[i + 1] - this.times[i]) * (x - i); };
    const zAt = (y) => { const k = Math.min(nz - 2, Math.floor(y)); return this.depths[k] + (this.depths[k + 1] - this.depths[k]) * (y - k); };
    const [lo, hi] = d.range;
    const step = ticks(lo, hi, 8)[1] - ticks(lo, hi, 8)[0] || 1;
    ctx.strokeStyle = 'rgba(255, 246, 238, 0.32)'; ctx.lineWidth = 0.8;
    ctx.beginPath();
    for (let level = Math.ceil(lo / step) * step; level < hi; level += step) {
      const seg = isolines(values, grid, level);
      for (let s = 0; s < seg.length; s += 4) {
        ctx.moveTo(L.x(tAt(seg[s])), L.y(r, zAt(seg[s + 1])));
        ctx.lineTo(L.x(tAt(seg[s + 2])), L.y(r, zAt(seg[s + 3])));
      }
    }
    ctx.stroke();
    // The mixed layer.
    if (this.mld) {
      ctx.strokeStyle = 'rgba(253, 241, 236, 0.95)'; ctx.lineWidth = 1.4; ctx.setLineDash([6, 4]);
      ctx.beginPath();
      let pen = false;
      this.mld.forEach((m, i) => {
        if (!Number.isFinite(m.depth)) { pen = false; return; }
        const x = L.x(this.times[i]), y = L.y(r, m.depth);
        if (pen) ctx.lineTo(x, y); else ctx.moveTo(x, y);
        pen = true;
      });
      ctx.stroke(); ctx.setLineDash([]);
      const first = this.mld.findIndex((m) => Number.isFinite(m.depth));
      if (first >= 0 && this.mld[first].depth < L.zMax * 0.8) {
        ctx.font = `italic 11.5px ${FONT}`; ctx.fillStyle = INK; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
        ctx.shadowColor = 'rgba(18, 13, 20, 0.95)'; ctx.shadowBlur = 4;
        ctx.fillText('mixed layer', L.x(this.times[first]) + 6, L.y(r, this.mld[first].depth) + 4);
        ctx.shadowBlur = 0;
      }
    }
  }

  drawColourBar(ctx, d, r) {
    const x = r.x + r.w + 14, w = 10, y0 = r.y + 4, h = r.h - 8;
    const [lo, hi] = d.range;
    for (let j = 0; j < h; j++) {
      const k = 4 * Math.round((1 - j / (h - 1)) * 255);
      ctx.fillStyle = `rgb(${d.lut[k]},${d.lut[k + 1]},${d.lut[k + 2]})`;
      ctx.fillRect(x, y0 + j, w, 1);
    }
    ctx.strokeStyle = 'rgba(240, 226, 238, 0.25)'; ctx.strokeRect(x + 0.5, y0 + 0.5, w - 1, h - 1);
    ctx.font = `11px ${FONT}`; ctx.fillStyle = DIM; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    const marks = d.signed ? [lo, lo / 2, 0, hi / 2, hi] : ticks(lo, hi, 4);
    for (const v of marks) {
      const y = y0 + (1 - (v - lo) / (hi - lo)) * (h - 1);
      ctx.fillRect(x + w, Math.round(y), 3, 1);
      ctx.fillText(`${v > 0 && d.signed ? '+' : ''}${(+v.toFixed(2)).toString().replace('-', '−')}`, x + w + 6, y);
    }
    ctx.font = `italic 11px ${FONT}`; ctx.fillStyle = FAINT; ctx.textBaseline = 'bottom';
    ctx.fillText(d.variable.units === 'degC' ? '°C' : d.variable.units.replace('m/s', 'm s⁻¹'), x - 2, r.y - 1);
  }

  // The Sun's height at the place, as a strip above the first panel.
  drawSky(ctx, L) {
    const r = L.rects[0], y = r.y - 12, h = 7;
    for (let x = Math.floor(r.x); x < r.x + r.w; x++) {
      const t = L.t(x + 0.5);
      const alt = sunAltAz(this.gridLon, this.gridLat, subsolarPoint(t)).altitude;
      const [cr, cg, cb] = skyColour(alt);
      ctx.fillStyle = `rgb(${cr},${cg},${cb})`;
      ctx.fillRect(x, y, 1, h);
    }
    ctx.font = `italic 11px ${FONT}`; ctx.fillStyle = FAINT; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    ctx.fillText('Sun', r.x - 8, y + h / 2);
    if (this.latestRun && this.latestRun > L.t0 && this.latestRun < L.t1) {
      const x = L.x(this.latestRun);
      ctx.font = `italic 11px ${FONT}`; ctx.textBaseline = 'bottom';
      ctx.fillStyle = DIM; ctx.textAlign = 'right'; ctx.fillText('analysis', x - 6, y - 2);
      ctx.textAlign = 'left'; ctx.fillText('forecast', x + 6, y - 2);
    }
  }

  drawTimeAxis(ctx, L) {
    const last = L.rects[L.rects.length - 1], y = last.y + last.h + 6;
    ctx.font = `11.5px ${FONT}`; ctx.fillStyle = DIM; ctx.textBaseline = 'top';
    for (let t = Math.ceil(L.t0 / DAY) * DAY; t <= L.t1; t += DAY) {
      const x = L.x(t);
      ctx.fillStyle = FAINT; ctx.fillRect(Math.round(x), y - 6, 1, 4);
      ctx.fillStyle = DIM; ctx.textAlign = 'center';
      if (x > last.x + 24 && x < last.x + last.w - 24) ctx.fillText(fmtDay(t), x, y);
    }
    for (let t = Math.ceil(L.t0 / (6 * HOUR)) * 6 * HOUR; t <= L.t1; t += 6 * HOUR) {
      if (t % DAY) { ctx.fillStyle = FAINT; ctx.fillRect(Math.round(L.x(t)), y - 6, 1, 2); }
    }
    ctx.fillStyle = FAINT; ctx.textAlign = 'right'; ctx.font = `italic 11px ${FONT}`;
    ctx.fillText('00:00 UTC', last.x - 8, y);
  }

  // The present, the start of the forecast, and the pointer's crosshair.
  drawOverlay() {
    if (!this.L) return;
    const L = this.L, dpr = window.devicePixelRatio || 1;
    const ctx = this.overlay.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, L.W, L.H);
    const first = L.rects[0], last = L.rects[L.rects.length - 1];
    const vline = (x, colour, dash) => {
      ctx.strokeStyle = colour; ctx.setLineDash(dash); ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.moveTo(Math.round(x) + 0.5, first.y - 14); ctx.lineTo(Math.round(x) + 0.5, last.y + last.h); ctx.stroke();
      ctx.setLineDash([]);
    };
    if (this.latestRun && this.latestRun > L.t0 && this.latestRun < L.t1) vline(L.x(this.latestRun), 'rgba(246, 238, 245, 0.55)', [3, 3]);
    const now = this.info.now;
    if (now >= L.t0 && now <= L.t1) {
      vline(L.x(now), GOLD, []);
      ctx.font = `600 11px ${FONT}`; ctx.fillStyle = GOLD; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
      ctx.fillText('time shown', L.x(now), first.y - 16);
    }
    if (this.cursor) {
      const { x, y, r } = this.cursor;
      ctx.strokeStyle = 'rgba(246, 238, 245, 0.7)'; ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(Math.round(x) + 0.5, first.y); ctx.lineTo(Math.round(x) + 0.5, last.y + last.h);
      for (const q of L.rects) { ctx.moveTo(q.x, Math.round(y - r.y + q.y) + 0.5); ctx.lineTo(q.x + q.w, Math.round(y - r.y + q.y) + 0.5); }
      ctx.stroke();
    }
  }

  point(e) {
    if (!this.L) return;
    const box = this.overlay.getBoundingClientRect();
    const x = e.clientX - box.left, y = e.clientY - box.top;
    const r = this.L.rects.find((q) => x >= q.x && x <= q.x + q.w && y >= q.y && y <= q.y + q.h);
    if (!r) { this.cursor = null; this.readout.hidden = true; this.drawOverlay(); return; }
    this.cursor = { x, y, r };
    this.drawOverlay();
    const t = this.L.t(x), z = this.L.z(r, y);
    const value = (std) => { const d = this.panels.find((p) => p.standard === std); return d ? d.sample(t, z) : NaN; };
    const u = value('eastward_sea_water_velocity'), v = value('northward_sea_water_velocity'), T = value('sea_water_temperature');
    const lines = [`<b>${esc(fmtTime(t))}</b> · ${Math.round(z).toLocaleString('en-GB')} m`];
    if (Number.isFinite(u) && Number.isFinite(v)) {
      lines.push(`current ${Math.hypot(u, v).toFixed(2)} m s⁻¹ toward ${compassPoint(flowToward(u, v))} <span>(east ${u >= 0 ? '+' : '−'}${Math.abs(u).toFixed(2)}, north ${v >= 0 ? '+' : '−'}${Math.abs(v).toFixed(2)})</span>`);
    }
    if (Number.isFinite(T)) lines.push(`temperature ${T.toFixed(1)} °C`);
    if (!Number.isFinite(u) && !Number.isFinite(T)) lines.push('<span>beneath the sea floor</span>');
    const alt = sunAltAz(this.gridLon, this.gridLat, subsolarPoint(t)).altitude;
    lines.push(`<span>Sun ${alt >= 0 ? `${alt.toFixed(0)}° up` : `${(-alt).toFixed(0)}° below the horizon`}${this.latestRun && t > this.latestRun ? ' · forecast' : ''}</span>`);
    this.readout.innerHTML = lines.join('<br>');
    this.readout.hidden = false;
    const W = this.L.W;
    this.readout.style.left = `${Math.min(W - 290, x + 14)}px`;
    this.readout.style.top = `${Math.max(4, y - 70)}px`;
  }
}
