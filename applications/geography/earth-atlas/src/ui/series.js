// The daily positions of a GPS station, drawn as east, north and up against
// time: how far the ground there has moved since the station began. A straight
// line is steady motion (the plate, or slow deformation); a jump is an
// earthquake or a change of equipment; a wave with a period of a year is the
// ground breathing with water, snow and heat. The line of the station's MIDAS
// velocity is drawn over the days. With the plate's rigid motion removed, what
// is left is the deformation near the station.

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const FONT = '"Source Serif 4", Georgia, serif';
const INK = '#f6eef5', DIM = '#cdb6c8', FAINT = '#9a8496', GOLD = '#f6c47e';
const YEAR_MS = 365.25 * 86400000;
const COLOURS = { east: '#e3b5d3', north: '#9fb6e8', up: '#f6c47e' };

function ticks(lo, hi, n = 4) {
  const raw = (hi - lo) / n, p = 10 ** Math.floor(Math.log10(raw || 1));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * p).find((s) => s >= raw) || 10 * p;
  const out = [];
  for (let x = Math.ceil(lo / step - 1e-9) * step; x <= hi + 1e-9; x += step) out.push(+x.toFixed(10));
  return out;
}

const quantile = (sorted, q) => sorted[Math.max(0, Math.min(sorted.length - 1, Math.floor(q * (sorted.length - 1))))];

export class SeriesFigure {
  constructor(host) {
    this.host = host;
    this.relative = false;
    host.innerHTML = `
      <div class="figure-head">
        <div class="figure-titles"><h3></h3><p class="sub"></p></div>
        <div class="figure-controls">
          <div class="segmented" role="group" aria-label="Frame">
            <button type="button" data-frame="igs">IGS20</button>
            <button type="button" data-frame="plate">Without the plate’s motion</button>
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
    host.querySelectorAll('[data-frame]').forEach((b) => b.addEventListener('click', () => { this.relative = b.dataset.frame === 'plate'; this.prepare(); this.draw(); }));
    host.querySelector('.close').addEventListener('click', () => this.hide());
    this.overlay.addEventListener('pointermove', (e) => this.point(e));
    this.overlay.addEventListener('pointerleave', () => { this.cursor = null; this.readout.hidden = true; this.drawOverlay(); });
    new ResizeObserver(() => { if (!this.host.hidden && this.product) this.draw(); }).observe(this.body);
  }

  get open() { return !this.host.hidden && Boolean(this.product); }

  // product: series/1; station: { id, ve, vn, vu, residual: { plate } }; events: [{ t, mw, km, name }]
  show(product, station, events = []) {
    this.product = product;
    this.station = station;
    this.events = events;
    this.prepare();
    this.host.hidden = false;
    this.draw();
  }

  hide() { this.host.hidden = true; this.product = null; this.readout.hidden = true; }

  prepare() {
    const p = this.product, st = this.station;
    this.times = p.days.map((d) => Date.parse(`${d}T12:00:00Z`));
    const t0 = this.times[0];
    // The plate's rigid motion at the station, in mm/yr, to remove when asked.
    const plate = st.residual ? st.residual.motion : null;
    this.canRemove = Boolean(plate);
    if (!plate) this.relative = false;
    const remove = this.relative ? { east: plate.east, north: plate.north, up: 0 } : { east: 0, north: 0, up: 0 };
    const rate = { east: st.ve, north: st.vn, up: st.vu };
    this.panels = p.components.map((c) => {
      const values = c.values.map((v, i) => v - (remove[c.name] * (this.times[i] - t0)) / YEAR_MS);
      const sorted = values.filter(Number.isFinite).sort((a, b) => a - b);
      // The MIDAS line, through the median of the values at the median time.
      const slope = rate[c.name] - remove[c.name];
      const tm = this.times[this.times.length >> 1];
      const residuals = values.map((v, i) => v - (slope * (this.times[i] - tm)) / YEAR_MS).filter(Number.isFinite).sort((a, b) => a - b);
      const intercept = quantile(residuals, 0.5);
      let lo = quantile(sorted, 0.003), hi = quantile(sorted, 0.997);
      const pad = Math.max(2, (hi - lo) * 0.06);
      lo -= pad; hi += pad;
      return { name: c.name, values, slope, intercept, tm, lo, hi, colour: COLOURS[c.name] };
    });
    this.host.querySelectorAll('[data-frame]').forEach((b) => {
      b.setAttribute('aria-pressed', String((b.dataset.frame === 'plate') === this.relative));
      b.disabled = b.dataset.frame === 'plate' && !this.canRemove;
    });
    const years = (this.times[this.times.length - 1] - t0) / YEAR_MS;
    this.host.querySelector('h3').textContent = `GPS station ${p.station}: where the ground has moved`;
    this.host.querySelector('.sub').textContent = [
      `${p.days[0]} to ${p.days[p.days.length - 1]}, ${p.days.length.toLocaleString('en-GB')} days over ${years.toFixed(1)} years`,
      this.relative ? `without the ${st.residual.plate.name} plate’s rigid motion (${st.residual.model})` : 'in IGS20 (ITRF2020)',
      'displacement from the first day',
    ].join(' · ');
    const where = typeof this.events.radius === 'number' ? `within ${this.events.radius} km` : this.events.radius;
    const quakes = this.events.length ? ` Dashed: earthquakes of magnitude ${this.events.minMw ?? 6.5} and more ${where} (Global CMT).` : '';
    this.host.querySelector('.caption').textContent =
      `Dots: daily positions; line: the station’s MIDAS velocity (Blewitt et al. 2016), ${Math.hypot(st.ve, st.vn).toFixed(1)} mm/yr horizontally and ${st.vu.toFixed(1)} mm/yr up. A step is an earthquake or new equipment; a yearly wave, the load of water and snow.${quakes} Source: ${p.source.credit}.`;
  }

  layout() {
    const W = this.body.clientWidth, H = this.body.clientHeight;
    const left = 70, right = 24, top = 10, bottom = 26, gap = 10;
    const n = this.panels.length, ph = (H - top - bottom - gap * (n - 1)) / n;
    const t0 = this.times[0], t1 = this.times[this.times.length - 1];
    return { W, H, left, right, rects: this.panels.map((_, i) => ({ x: left, y: top + i * (ph + gap), w: W - left - right, h: ph })), t0, t1,
      x: (t) => left + ((t - t0) / (t1 - t0)) * (W - left - right), t: (x) => t0 + ((x - left) / (W - left - right)) * (t1 - t0) };
  }

  draw() {
    if (!this.product) return;
    this.cursor = null; this.readout.hidden = true;
    const dpr = window.devicePixelRatio || 1, L = (this.L = this.layout());
    for (const c of [this.plot, this.overlay]) {
      c.width = Math.round(L.W * dpr); c.height = Math.round(L.H * dpr);
      c.style.width = `${L.W}px`; c.style.height = `${L.H}px`;
    }
    const ctx = this.plot.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, L.W, L.H);
    this.panels.forEach((d, i) => {
      const r = L.rects[i], y = (v) => r.y + r.h - ((v - d.lo) / (d.hi - d.lo)) * r.h;
      ctx.fillStyle = 'rgba(255, 255, 255, 0.025)'; ctx.fillRect(r.x, r.y, r.w, r.h);
      ctx.strokeStyle = 'rgba(240, 226, 238, 0.18)'; ctx.strokeRect(r.x + 0.5, r.y + 0.5, r.w - 1, r.h - 1);
      // Years, faintly.
      const y0 = new Date(L.t0).getUTCFullYear() + 1, y1 = new Date(L.t1).getUTCFullYear();
      const every = Math.max(1, Math.ceil((y1 - y0 + 1) / 14));
      ctx.fillStyle = 'rgba(246, 238, 245, 0.08)';
      for (let yr = y0; yr <= y1; yr++) if ((yr - y0) % every === 0) ctx.fillRect(Math.round(L.x(Date.UTC(yr, 0, 1))), r.y, 1, r.h);
      // Earthquakes nearby.
      ctx.save(); ctx.beginPath(); ctx.rect(r.x, r.y, r.w, r.h); ctx.clip();
      ctx.strokeStyle = 'rgba(239, 138, 152, 0.75)'; ctx.setLineDash([3, 3]); ctx.lineWidth = 1;
      for (const e of this.events) { if (e.t < L.t0 || e.t > L.t1) continue; ctx.beginPath(); ctx.moveTo(Math.round(L.x(e.t)) + 0.5, r.y); ctx.lineTo(Math.round(L.x(e.t)) + 0.5, r.y + r.h); ctx.stroke(); }
      ctx.setLineDash([]);
      // The days.
      ctx.fillStyle = d.colour; ctx.globalAlpha = Math.min(0.85, Math.max(0.25, 2400 / d.values.length));
      d.values.forEach((v, k) => { if (Number.isFinite(v)) ctx.fillRect(L.x(this.times[k]) - 0.8, y(v) - 0.8, 1.6, 1.6); });
      ctx.globalAlpha = 1;
      // The MIDAS line.
      ctx.strokeStyle = INK; ctx.lineWidth = 1.2;
      ctx.beginPath();
      ctx.moveTo(L.x(L.t0), y(d.intercept + (d.slope * (L.t0 - d.tm)) / YEAR_MS));
      ctx.lineTo(L.x(L.t1), y(d.intercept + (d.slope * (L.t1 - d.tm)) / YEAR_MS));
      ctx.stroke();
      ctx.restore();
      // Axis and title.
      ctx.font = `11px ${FONT}`; ctx.fillStyle = FAINT; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
      for (const v of ticks(d.lo, d.hi, Math.max(2, Math.round(r.h / 40)))) {
        const yy = y(v);
        if (yy < r.y + 6 || yy > r.y + r.h - 6) continue;
        ctx.fillText(`${v > 0 ? '+' : ''}${v.toLocaleString('en-GB').replace('-', '−')}`, r.x - 8, yy);
        ctx.fillRect(r.x - 4, Math.round(yy), 3, 1);
      }
      ctx.fillStyle = 'rgba(20, 14, 22, 0.78)';
      ctx.font = `600 13px ${FONT}`;
      const title = d.name === 'up' ? 'Up' : d.name === 'east' ? 'East' : 'North';
      const rate = `${d.slope >= 0 ? '+' : '−'}${Math.abs(d.slope).toFixed(1)} mm/yr`;
      const tw = ctx.measureText(title).width;
      ctx.font = `italic 12px ${FONT}`;
      const rw = ctx.measureText(rate).width;
      ctx.fillStyle = 'rgba(20, 14, 22, 0.78)';
      ctx.beginPath(); ctx.roundRect(r.x + 5, r.y + 5, tw + rw + 22, 20, 3); ctx.fill();
      ctx.textAlign = 'left'; ctx.textBaseline = 'top';
      ctx.font = `600 13px ${FONT}`; ctx.fillStyle = d.colour; ctx.fillText(title, r.x + 11, r.y + 8);
      ctx.font = `italic 12px ${FONT}`; ctx.fillStyle = DIM; ctx.fillText(rate, r.x + 17 + tw, r.y + 9);
      if (i === 0) {
        ctx.save(); ctx.translate(14, r.y + r.h * 1.5); ctx.rotate(-Math.PI / 2);
        ctx.textAlign = 'center'; ctx.font = `italic 12px ${FONT}`; ctx.fillStyle = FAINT; ctx.fillText('millimetres', 0, 0); ctx.restore();
      }
    });
    // Years along the bottom.
    const last = L.rects[L.rects.length - 1], y = last.y + last.h + 6;
    const y0 = new Date(L.t0).getUTCFullYear() + 1, y1 = new Date(L.t1).getUTCFullYear();
    const every = Math.max(1, Math.ceil((y1 - y0 + 1) / 14));
    ctx.font = `11.5px ${FONT}`; ctx.fillStyle = DIM; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (let yr = y0; yr <= y1; yr++) if ((yr - y0) % every === 0) ctx.fillText(String(yr), L.x(Date.UTC(yr, 0, 1)), y);
    this.drawOverlay();
  }

  drawOverlay() {
    if (!this.L) return;
    const L = this.L, dpr = window.devicePixelRatio || 1, ctx = this.overlay.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, L.W, L.H);
    if (!this.cursor) return;
    const first = L.rects[0], last = L.rects[L.rects.length - 1];
    ctx.strokeStyle = 'rgba(246, 238, 245, 0.7)';
    ctx.beginPath(); ctx.moveTo(Math.round(this.cursor.x) + 0.5, first.y); ctx.lineTo(Math.round(this.cursor.x) + 0.5, last.y + last.h); ctx.stroke();
  }

  point(e) {
    if (!this.L) return;
    const box = this.overlay.getBoundingClientRect(), x = e.clientX - box.left, y = e.clientY - box.top;
    if (x < this.L.left || x > this.L.W - this.L.right) { this.cursor = null; this.readout.hidden = true; this.drawOverlay(); return; }
    this.cursor = { x, y };
    this.drawOverlay();
    const t = this.L.t(x);
    // The nearest day.
    let lo = 0, hi = this.times.length - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (this.times[m] <= t) lo = m; else hi = m; }
    const k = Math.abs(this.times[hi] - t) < Math.abs(this.times[lo] - t) ? hi : lo;
    const lines = [`<b>${esc(this.product.days[k])}</b>`];
    for (const d of this.panels) lines.push(`<span style="color:${d.colour}">${d.name}</span> ${d.values[k] >= 0 ? '+' : '−'}${Math.abs(d.values[k]).toFixed(1)} mm`);
    const near = this.events.filter((ev) => Math.abs(ev.t - this.times[k]) < 20 * 86400000);
    for (const ev of near) lines.push(`<span>M ${ev.mw.toFixed(1)}, ${Math.round(ev.km)} km away, ${new Date(ev.t).toISOString().slice(0, 10)}</span>`);
    this.readout.innerHTML = lines.join('<br>');
    this.readout.hidden = false;
    this.readout.style.left = `${Math.min(this.L.W - 220, x + 14)}px`;
    this.readout.style.top = `${Math.max(4, y - 60)}px`;
  }
}
