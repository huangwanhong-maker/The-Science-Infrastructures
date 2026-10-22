// A section through the Earth: depth against distance along a great circle.
// Above, the relief (exaggerated); below, the crust to the Moho (CRUST1.0),
// and the earthquakes near the section (Global CMT centroids), coloured by
// their faulting. Across a subduction zone they trace the sinking plate, the
// Wadati-Benioff zone, down to about 660 km, where the mantle's minerals change
// phase and deep earthquakes end.

import { sectionLine, relativeTo } from '../core/section.js';
import { KIND_COLOURS } from '../layers/mechanisms.js';

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const FONT = '"Source Serif 4", Georgia, serif';
const INK = '#f6eef5', DIM = '#cdb6c8', FAINT = '#9a8496', GOLD = '#f6c47e';
const coords = (lon, lat) => `${Math.abs(lat).toFixed(1)}° ${lat >= 0 ? 'N' : 'S'} ${Math.abs(lon).toFixed(1)}° ${lon >= 0 ? 'E' : 'W'}`;
const compass = (az) => ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW'][Math.round((((az % 360) + 360) % 360) / 22.5) % 16];

export class SectionFigure {
  // sources: { quakes(): [{ lon, lat, depth, mw, kind, t, name }], elevation(lon, lat) m, crust(lon, lat) column, crossings(points) [{ s, label }] }
  constructor(host, sources, { onChange } = {}) {
    this.host = host;
    this.sources = sources;
    this.onChange = onChange;
    this.swath = 100;
    host.innerHTML = `
      <div class="figure-head">
        <div class="figure-titles"><h3></h3><p class="sub"></p></div>
        <div class="figure-controls">
          <div class="segmented" role="group" aria-label="Turn the section">
            <button type="button" data-turn="-15">↺ 15°</button><button type="button" data-turn="15">↻ 15°</button>
          </div>
          <div class="segmented" role="group" aria-label="Width of the band">
            ${[50, 100, 200].map((w) => `<button type="button" data-swath="${w}">±${w} km</button>`).join('')}
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
    host.querySelectorAll('[data-turn]').forEach((b) => b.addEventListener('click', () => { this.az = (this.az + Number(b.dataset.turn) + 360) % 360; this.prepare(); this.draw(); }));
    host.querySelectorAll('[data-swath]').forEach((b) => b.addEventListener('click', () => { this.swath = Number(b.dataset.swath); this.prepare(); this.draw(); }));
    host.querySelector('.close').addEventListener('click', () => this.hide());
    this.overlay.addEventListener('pointermove', (e) => this.point(e));
    this.overlay.addEventListener('pointerleave', () => { this.readout.hidden = true; this.cursor = null; this.drawOverlay(); });
    new ResizeObserver(() => { if (this.open) this.draw(); }).observe(this.body);
  }

  get open() { return !this.host.hidden && Boolean(this.line); }

  // centre, azimuth (degrees) and half-length (km); why: a sentence on how the section was chosen.
  // mark: where along the section (km) the place chosen lies.
  show({ lon, lat, az, half, why, mark = 0 }) {
    Object.assign(this, { lon, lat, az, half, why, centreS: mark });
    this.prepare();
    this.host.hidden = false;
    this.draw();
  }

  hide() { this.host.hidden = true; this.line = null; this.readout.hidden = true; if (this.onChange) this.onChange(); }

  prepare() {
    this.line = sectionLine(this.lon, this.lat, this.az, this.half, 5);
    for (const p of this.line) {
      p.z = this.sources.elevation(p.lon, p.lat);
      p.crust = this.sources.crust(p.lon, p.lat);
    }
    this.crossings = this.sources.crossings(this.line);
    this.quakes = [];
    for (const q of this.sources.quakes()) {
      const r = relativeTo(this.lon, this.lat, this.az, q.lon, q.lat);
      if (Math.abs(r.across) <= this.swath && Math.abs(r.along) <= this.half) this.quakes.push({ ...q, s: r.along, across: r.across });
    }
    this.quakes.sort((a, b) => a.mw - b.mw);
    const deepest = this.quakes.reduce((m, q) => Math.max(m, q.depth), 0);
    this.zMax = Math.max(200, Math.min(750, Math.ceil((deepest + 40) / 100) * 100));
    const a = this.line[0], b = this.line[this.line.length - 1];
    this.host.querySelector('h3').textContent = `A section through the Earth, toward ${compass(this.az)} (${Math.round(this.az)}°)`;
    this.host.querySelector('.sub').textContent = `from ${coords(a.lon, a.lat)} to ${coords(b.lon, b.lat)} · ${(2 * this.half).toLocaleString('en-GB')} km · ${this.quakes.length.toLocaleString('en-GB')} earthquakes within ${this.swath} km of the line · ${this.why}`;
    this.host.querySelectorAll('[data-swath]').forEach((x) => x.setAttribute('aria-pressed', String(Number(x.dataset.swath) === this.swath)));
    if (this.onChange) this.onChange();
  }

  layout() {
    const W = this.body.clientWidth, H = this.body.clientHeight;
    const left = 64, right = 22, top = 8, bottom = 26, strip = Math.max(46, Math.round(H * 0.18)), gap = 8;
    const x = (s) => left + ((s + this.half) / (2 * this.half)) * (W - left - right);
    const relief = { x: left, y: top, w: W - left - right, h: strip };
    const depth = { x: left, y: top + strip + gap, w: W - left - right, h: H - top - strip - gap - bottom };
    const kmPerPx = (2 * this.half) / depth.w, kmPerPxV = this.zMax / depth.h;
    return { W, H, left, right, relief, depth, x, s: (px) => ((px - left) / (W - left - right)) * 2 * this.half - this.half,
      zy: (z) => depth.y + (z / this.zMax) * depth.h, exaggeration: kmPerPx / kmPerPxV };
  }

  draw() {
    if (!this.line) return;
    this.cursor = null; this.readout.hidden = true;
    const dpr = window.devicePixelRatio || 1, L = (this.L = this.layout());
    for (const c of [this.plot, this.overlay]) {
      c.width = Math.round(L.W * dpr); c.height = Math.round(L.H * dpr);
      c.style.width = `${L.W}px`; c.style.height = `${L.H}px`;
    }
    const ctx = this.plot.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, L.W, L.H);
    this.drawRelief(ctx, L);
    this.drawDepth(ctx, L);
    this.drawOverlay();
    this.host.querySelector('.caption').textContent = `Earthquakes within ${this.swath} km of the section: Global CMT centroids since 1976, sized by magnitude and coloured by faulting (thrust, normal, strike-slip, oblique). Beneath a trench they deepen along the sinking plate, the Wadati–Benioff zone, and end near 660 km, where the mantle's minerals change phase. Crust to the Moho: CRUST1.0. Relief above: ETOPO 2022, ${Math.round(this.reliefExaggeration)}× exaggerated; below, depth is ${L.exaggeration.toFixed(1)}× exaggerated against distance.`;
  }

  drawRelief(ctx, L) {
    const r = L.relief, lo = -8000, hi = 6000;
    const y = (z) => r.y + ((hi - Math.max(lo, Math.min(hi, z))) / (hi - lo)) * r.h;
    this.reliefExaggeration = ((hi - lo) / 1000 / r.h) ** -1 * ((2 * this.half) / r.w);
    ctx.fillStyle = 'rgba(169, 189, 224, 0.10)'; ctx.fillRect(r.x, y(0), r.w, r.y + r.h - y(0));
    // The ground, and the sea above it.
    ctx.beginPath(); ctx.moveTo(r.x, r.y + r.h);
    for (const p of this.line) ctx.lineTo(L.x(p.s), y(Number.isFinite(p.z) ? p.z : 0));
    ctx.lineTo(r.x + r.w, r.y + r.h); ctx.closePath();
    ctx.fillStyle = '#5a4651'; ctx.fill();
    ctx.beginPath();
    for (const [i, p] of this.line.entries()) { const yy = y(Number.isFinite(p.z) ? p.z : 0); if (i) ctx.lineTo(L.x(p.s), yy); else ctx.moveTo(L.x(p.s), yy); }
    ctx.strokeStyle = '#e8c3dc'; ctx.lineWidth = 1.2; ctx.stroke();
    ctx.strokeStyle = 'rgba(169, 189, 224, 0.7)'; ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(r.x, y(0)); ctx.lineTo(r.x + r.w, y(0)); ctx.stroke(); ctx.setLineDash([]);
    ctx.font = `11px ${FONT}`; ctx.fillStyle = FAINT; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const z of [-6000, -3000, 0, 3000]) ctx.fillText(`${z === 0 ? '0' : (z / 1000).toFixed(0).replace('-', '−')}${z === 3000 ? ' km' : ''}`, r.x - 8, y(z));
    // Plate boundaries where the section crosses them.
    ctx.font = `italic 11px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (const c of this.crossings) {
      const xx = L.x(c.s);
      ctx.fillStyle = c.colour; ctx.beginPath(); ctx.moveTo(xx, r.y + 1); ctx.lineTo(xx - 5, r.y - 6 + 1); ctx.lineTo(xx + 5, r.y - 6 + 1); ctx.closePath();
      ctx.fillRect(Math.round(xx) - 0.5, r.y, 1, r.h);
      ctx.fillStyle = INK; ctx.fillText(c.label, Math.max(r.x + 50, Math.min(r.x + r.w - 50, xx)), r.y + 3);
    }
  }

  drawDepth(ctx, L) {
    const d = L.depth;
    // The mantle, darkening with depth.
    const g = ctx.createLinearGradient(0, d.y, 0, d.y + d.h);
    g.addColorStop(0, '#3a2a36'); g.addColorStop(1, '#1d1420');
    ctx.fillStyle = g; ctx.fillRect(d.x, d.y, d.w, d.h);
    // The crust: sediments and crystalline crust, from the solid surface to the Moho.
    const band = (top, bottom, colour) => {
      ctx.beginPath();
      this.line.forEach((p, i) => { if (!p.crust) return; const v = L.zy(Math.max(0, top(p))); if (i) ctx.lineTo(L.x(p.s), v); else ctx.moveTo(L.x(p.s), v); });
      for (let i = this.line.length - 1; i >= 0; i--) { const p = this.line[i]; if (p.crust) ctx.lineTo(L.x(p.s), L.zy(Math.max(0, bottom(p)))); }
      ctx.closePath(); ctx.fillStyle = colour; ctx.fill();
    };
    band((p) => -p.crust.solid, (p) => -p.crust.solid + p.crust.sediments, 'rgba(217, 191, 140, 0.8)');
    band((p) => -p.crust.solid + p.crust.sediments, (p) => p.crust.moho, 'rgba(201, 154, 166, 0.75)');
    ctx.beginPath();
    this.line.forEach((p, i) => { if (!p.crust) return; const v = L.zy(p.crust.moho); if (i) ctx.lineTo(L.x(p.s), v); else ctx.moveTo(L.x(p.s), v); });
    ctx.strokeStyle = GOLD; ctx.lineWidth = 1.2; ctx.stroke();
    // The transition zone's two discontinuities.
    ctx.font = `italic 11px ${FONT}`; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom';
    for (const z of [410, 660]) {
      if (z > this.zMax) continue;
      ctx.strokeStyle = 'rgba(246, 238, 245, 0.22)'; ctx.setLineDash([5, 5]);
      ctx.beginPath(); ctx.moveTo(d.x, L.zy(z)); ctx.lineTo(d.x + d.w, L.zy(z)); ctx.stroke(); ctx.setLineDash([]);
      ctx.fillStyle = FAINT; ctx.fillText(`${z} km`, d.x + d.w - 6, L.zy(z) - 2);
    }
    // Earthquakes.
    for (const q of this.quakes) {
      const rr = Math.max(1.6, 1.4 + (q.mw - 5) * 2.4);
      ctx.beginPath(); ctx.arc(L.x(q.s), L.zy(q.depth), rr, 0, 2 * Math.PI);
      ctx.fillStyle = KIND_COLOURS[q.kind] || INK; ctx.globalAlpha = 0.78; ctx.fill(); ctx.globalAlpha = 1;
      if (q.mw >= 7) { ctx.strokeStyle = '#2b1d2d'; ctx.lineWidth = 1; ctx.stroke(); }
    }
    // Frame, axes, the place chosen.
    ctx.strokeStyle = 'rgba(240, 226, 238, 0.2)'; ctx.strokeRect(d.x + 0.5, d.y + 0.5, d.w - 1, d.h - 1);
    ctx.font = `11.5px ${FONT}`; ctx.fillStyle = FAINT; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (let z = 0; z <= this.zMax; z += this.zMax > 300 ? 100 : 50) { ctx.fillText(String(z), d.x - 8, L.zy(z)); ctx.fillRect(d.x - 4, Math.round(L.zy(z)), 3, 1); }
    ctx.save(); ctx.translate(14, d.y + d.h / 2); ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center'; ctx.font = `italic 12px ${FONT}`; ctx.fillText('depth (km)', 0, 0); ctx.restore();
    ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillStyle = DIM;
    const step = this.half > 600 ? 200 : 100, last = Math.floor(this.half / step) * step;
    for (let s = -last; s <= this.half; s += step) {
      ctx.fillRect(Math.round(L.x(s)), d.y + d.h, 1, 4);
      ctx.textAlign = s === last ? 'right' : 'center';
      ctx.fillText(`${s === 0 ? '0' : `${s > 0 ? '' : '−'}${Math.abs(s)}`}${s === last ? ' km' : ''}`, L.x(s) + (s === last ? 4 : 0), d.y + d.h + 6);
    }
    ctx.fillStyle = GOLD; ctx.beginPath();
    ctx.moveTo(L.x(this.centreS || 0), d.y - 1); ctx.lineTo(L.x(this.centreS || 0) - 5, d.y - 8); ctx.lineTo(L.x(this.centreS || 0) + 5, d.y - 8); ctx.closePath(); ctx.fill();
  }

  drawOverlay() {
    if (!this.L) return;
    const L = this.L, dpr = window.devicePixelRatio || 1, ctx = this.overlay.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, L.W, L.H);
    if (!this.cursor) return;
    ctx.strokeStyle = 'rgba(246, 238, 245, 0.6)';
    ctx.beginPath(); ctx.moveTo(Math.round(this.cursor.x) + 0.5, L.relief.y); ctx.lineTo(Math.round(this.cursor.x) + 0.5, L.depth.y + L.depth.h); ctx.stroke();
    if (this.cursor.quake) { ctx.strokeStyle = INK; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(this.cursor.qx, this.cursor.qy, 7, 0, 2 * Math.PI); ctx.stroke(); }
  }

  point(e) {
    if (!this.L) return;
    const L = this.L, box = this.overlay.getBoundingClientRect(), x = e.clientX - box.left, y = e.clientY - box.top;
    if (x < L.left || x > L.W - L.right) { this.cursor = null; this.readout.hidden = true; this.drawOverlay(); return; }
    let best = null, bestD = 9;
    for (const q of this.quakes) {
      const dd = Math.hypot(L.x(q.s) - x, L.zy(q.depth) - y);
      if (dd < bestD) { bestD = dd; best = q; }
    }
    const s = L.s(x);
    const p = this.line[Math.max(0, Math.min(this.line.length - 1, Math.round((s + this.half) / 5)))];
    this.cursor = { x, y, quake: best, qx: best && L.x(best.s), qy: best && L.zy(best.depth) };
    this.drawOverlay();
    const lines = [];
    if (best) {
      lines.push(`<b>M ${best.mw.toFixed(1)}</b> · ${esc(best.kind)} · ${best.depth.toFixed(0)} km deep`);
      lines.push(`<span>${new Date(best.t).toISOString().slice(0, 10)}, ${esc(best.name || '')}, ${Math.abs(best.across).toFixed(0)} km ${best.across >= 0 ? 'to the right of' : 'to the left of'} the line</span>`);
    } else {
      lines.push(`<b>${s >= 0 ? '' : '−'}${Math.abs(s).toFixed(0)} km</b> · ${coords(p.lon, p.lat)}`);
      if (Number.isFinite(p.z)) lines.push(`${p.z >= 0 ? 'elevation' : 'depth of the sea'} ${Math.abs(p.z).toLocaleString('en-GB', { maximumFractionDigits: 0 })} m`);
      if (p.crust) lines.push(`<span>crust ${p.crust.crust.toFixed(1)} km thick, Moho at ${p.crust.moho.toFixed(1)} km</span>`);
    }
    this.readout.innerHTML = lines.join('<br>');
    this.readout.hidden = false;
    this.readout.style.left = `${Math.min(L.W - 300, x + 14)}px`;
    this.readout.style.top = `${Math.max(4, y - 50)}px`;
  }
}
