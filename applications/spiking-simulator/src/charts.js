// Charts on canvas: one figure per chart, its title and legend in HTML, its marks on a canvas
// sized to the device's pixels, and a second canvas over it for the crosshair. Axes and grids
// are hairlines one step off the surface; text wears text colours, never a series' colour;
// lines are 2 px; a tooltip on every chart gives the values under the pointer, the value
// first. A chart is drawn again when its box changes size.

import { el } from './dom.js';
import { short } from './units.js';

const INK = { text: '#cfb9ca', faint: '#9a8496', grid: 'rgba(240, 226, 238, 0.07)', axis: 'rgba(240, 226, 238, 0.22)', cross: 'rgba(246, 238, 245, 0.55)', surface: '#1d1520' };
const FONT = "12px 'Source Serif 4', Georgia, serif";

// Ticks at 1, 2 or 5 times a power of ten, about `count` of them.
export function ticks(lo, hi, count = 6) {
  if (!(hi > lo)) return [lo];
  const raw = (hi - lo) / count, p = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((m) => m * p).find((s) => s >= raw * 0.999);
  const out = [];
  for (let k = Math.ceil(lo / step - 1e-9); k * step <= hi + step * 1e-9; k++) out.push(k === 0 ? 0 : Number((k * step).toPrecision(12)));
  return out;
}

// Decades between lo and hi (both positive).
export function decades(lo, hi) {
  const out = [];
  for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) { const v = 10 ** e; if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v); }
  return out;
}

export function scale(domain, range, log = false) {
  const [d0, d1] = domain, [r0, r1] = range;
  if (log) {
    const a = Math.log(d0), b = Math.log(d1);
    const f = (v) => r0 + ((Math.log(v) - a) / (b - a)) * (r1 - r0);
    f.invert = (p) => Math.exp(a + ((p - r0) / (r1 - r0)) * (b - a));
    return f;
  }
  const f = (v) => r0 + ((v - d0) / (d1 - d0)) * (r1 - r0);
  f.invert = (p) => d0 + ((p - r0) / (r1 - r0)) * (d1 - d0);
  return f;
}

// The range of finite values, padded a little; [0, 1] if none.
export function extent(values, { zero = false, pad = 0.05 } = {}) {
  let lo = Infinity, hi = -Infinity;
  for (const v of values) if (Number.isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; }
  if (lo === Infinity) return [0, 1];
  if (zero) { lo = Math.min(lo, 0); hi = Math.max(hi, 0); }
  if (hi === lo) { hi += Math.abs(hi) * 0.1 || 1; lo -= zero && lo === 0 ? 0 : Math.abs(lo) * 0.1 || 1; }
  const p = (hi - lo) * pad;
  return [zero && lo === 0 ? 0 : lo - p, hi + p];
}

const fmt = (v) => short(v);

export class Figure {
  // opts: title, subtitle, height, legend [{label, colour, kind: 'line'|'dot'|'box'}], x {label, log, format},
  // y {label, log, format}, margin, hover(x, y) -> {title, rows: [{colour, label, value}]} | null,
  // onSelect([x0, x1]) (drag to choose a range), onReset() (double-click).
  constructor(parent, opts) {
    this.opts = opts;
    this.legend = el('div', { class: 'legend' });
    this.canvas = el('canvas');
    this.overlay = el('canvas', { class: 'overlay' });
    this.tip = el('div', { class: 'tip', hidden: true });
    this.box = el('div', { class: 'plot', style: `height:${opts.height ?? 220}px` }, this.canvas, this.overlay, this.tip);
    this.node = el('figure', { class: 'fig' },
      el('figcaption', {}, el('span', { class: 'fig-title', text: opts.title ?? '' }), opts.subtitle ? el('span', { class: 'fig-sub', text: opts.subtitle }) : null),
      this.legend, this.box);
    parent.append(this.node);
    if (opts.legend?.length) for (const l of opts.legend) this.legend.append(el('span', { class: 'key' }, el('i', { class: `k-${l.kind ?? 'line'}`, style: `background:${l.colour}` }), l.label));
    this.margin = { left: 58, right: 14, top: 8, bottom: 34, ...opts.margin };
    this.observer = new ResizeObserver(() => this.render());
    this.observer.observe(this.box);
    this.wire();
  }

  // draw(ctx, x, y, area): the marks; x and y are scales from data to pixels.
  set(domainX, domainY, draw) {
    this.dx = domainX;
    this.dy = domainY;
    this.draw = draw;
    this.render();
  }

  render() {
    if (!this.draw) return;
    const w = this.box.clientWidth, h = this.box.clientHeight;
    if (!w || !h) return;
    const r = window.devicePixelRatio || 1;
    for (const c of [this.canvas, this.overlay]) { c.width = Math.round(w * r); c.height = Math.round(h * r); c.style.width = `${w}px`; c.style.height = `${h}px`; }
    const ctx = this.canvas.getContext('2d');
    ctx.setTransform(r, 0, 0, r, 0, 0);
    ctx.clearRect(0, 0, w, h);
    const m = this.margin, a = { left: m.left, right: w - m.right, top: m.top, bottom: h - m.bottom };
    this.area = a;
    this.x = scale(this.dx, [a.left, a.right], this.opts.x?.log);
    this.y = scale(this.dy, [a.bottom, a.top], this.opts.y?.log);
    this.axes(ctx, a);
    ctx.save();
    ctx.beginPath();
    ctx.rect(a.left, a.top - 1, a.right - a.left, a.bottom - a.top + 2);
    ctx.clip();
    this.draw(ctx, this.x, this.y, a);
    ctx.restore();
  }

  axes(ctx, a) {
    const { x: ox = {}, y: oy = {} } = this.opts;
    ctx.font = FONT;
    ctx.lineWidth = 1;
    const xt = ox.ticks?.(this.dx) ?? (ox.log ? decades(...this.dx) : ticks(...this.dx, Math.max(2, Math.floor((a.right - a.left) / 90))));
    const yt = oy.ticks?.(this.dy) ?? (oy.log ? decades(...this.dy) : ticks(...this.dy, Math.max(2, Math.floor((a.bottom - a.top) / 42))));
    ctx.strokeStyle = INK.grid;
    if (oy.grid !== false) for (const v of yt) { const p = Math.round(this.y(v)) + 0.5; if (p < a.top - 1 || p > a.bottom + 1) continue; ctx.beginPath(); ctx.moveTo(a.left, p); ctx.lineTo(a.right, p); ctx.stroke(); }
    if (ox.grid) for (const v of xt) { const p = Math.round(this.x(v)) + 0.5; ctx.beginPath(); ctx.moveTo(p, a.top); ctx.lineTo(p, a.bottom); ctx.stroke(); }
    ctx.strokeStyle = INK.axis;
    ctx.beginPath(); ctx.moveTo(a.left, Math.round(a.bottom) + 0.5); ctx.lineTo(a.right, Math.round(a.bottom) + 0.5); ctx.stroke();
    ctx.fillStyle = INK.faint;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (const v of xt) { const p = this.x(v); if (p < a.left - 1 || p > a.right + 1) continue; ctx.fillText((ox.format ?? fmt)(v), p, a.bottom + 6); }
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    if (!oy.hideTicks) for (const v of yt) { const p = this.y(v); if (p < a.top - 1 || p > a.bottom + 1) continue; ctx.fillText((oy.format ?? fmt)(v), a.left - 8, p); }
    ctx.textAlign = 'right';
    ctx.textBaseline = 'bottom';
    if (ox.label) ctx.fillText(ox.label, a.right, a.bottom + 31);
    if (oy.label) { ctx.save(); ctx.translate(13, a.top); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(oy.label, 0, 0); ctx.restore(); }
  }

  wire() {
    let drag = null;
    const local = (e) => { const r = this.overlay.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
    const inside = (px, py) => this.area && px >= this.area.left && px <= this.area.right && py >= this.area.top && py <= this.area.bottom;
    this.overlay.addEventListener('pointermove', (e) => {
      const [px, py] = local(e);
      const ctx = this.overlay.getContext('2d');
      const r = window.devicePixelRatio || 1;
      ctx.setTransform(r, 0, 0, r, 0, 0);
      ctx.clearRect(0, 0, this.overlay.width, this.overlay.height);
      if (!inside(px, py)) { this.tip.hidden = true; return; }
      if (drag) {
        ctx.fillStyle = 'rgba(202, 141, 190, 0.14)';
        ctx.fillRect(Math.min(drag, px), this.area.top, Math.abs(px - drag), this.area.bottom - this.area.top);
      }
      ctx.strokeStyle = INK.cross;
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(Math.round(px) + 0.5, this.area.top); ctx.lineTo(Math.round(px) + 0.5, this.area.bottom); ctx.stroke();
      const info = this.opts.hover?.(this.x.invert(px), this.y.invert(py));
      if (!info) { this.tip.hidden = true; return; }
      this.tip.replaceChildren(el('div', { class: 'tip-title', text: info.title ?? '' }),
        ...(info.rows ?? []).map((row) => el('div', { class: 'tip-row' }, row.colour ? el('i', { style: `background:${row.colour}` }) : null, el('b', { text: row.value }), el('span', { text: row.label ?? '' }))));
      this.tip.hidden = false;
      const w = this.box.clientWidth;
      this.tip.style.left = `${px + 14 + 220 > w ? px - 14 - this.tip.offsetWidth : px + 14}px`;
      this.tip.style.top = `${Math.max(0, py - 10)}px`;
    });
    this.overlay.addEventListener('pointerleave', () => {
      this.tip.hidden = true;
      const ctx = this.overlay.getContext('2d');
      ctx.clearRect(0, 0, this.overlay.width, this.overlay.height);
    });
    if (this.opts.onSelect) {
      this.overlay.style.cursor = 'crosshair';
      this.overlay.addEventListener('pointerdown', (e) => { const [px, py] = local(e); if (inside(px, py)) { drag = px; this.overlay.setPointerCapture(e.pointerId); } });
      this.overlay.addEventListener('pointerup', (e) => {
        if (drag === null) return;
        const [px] = local(e);
        const a = Math.min(drag, px), b = Math.max(drag, px);
        drag = null;
        if (b - a > 6) this.opts.onSelect([this.x.invert(a), this.x.invert(b)].sort((u, v) => u - v));
      });
      this.overlay.addEventListener('dblclick', () => this.opts.onReset?.());
    }
  }

  destroy() { this.observer.disconnect(); }
}

// ---- marks ----
export function line(ctx, x, y, xs, ys, colour, width = 2, alpha = 1) {
  ctx.save();
  ctx.globalAlpha = alpha;
  ctx.strokeStyle = colour;
  ctx.lineWidth = width;
  ctx.lineJoin = 'round';
  ctx.lineCap = 'round';
  ctx.beginPath();
  let on = false;
  for (let i = 0; i < xs.length; i++) {
    const v = ys[i];
    if (!Number.isFinite(v)) { on = false; continue; }
    const px = x(xs[i]), py = y(v);
    if (!on) { ctx.moveTo(px, py); on = true; } else ctx.lineTo(px, py);
  }
  ctx.stroke();
  ctx.restore();
}

export function band(ctx, x, y, xs, lo, hi, colour, alpha = 0.12) {
  ctx.save();
  ctx.globalAlpha = alpha;
  ctx.fillStyle = colour;
  ctx.beginPath();
  for (let i = 0; i < xs.length; i++) ctx.lineTo(x(xs[i]), y(hi[i]));
  for (let i = xs.length - 1; i >= 0; i--) ctx.lineTo(x(xs[i]), y(lo[i]));
  ctx.closePath();
  ctx.fill();
  ctx.restore();
}

// Histogram bars from edges and counts, touching bars separated by a 2 px gap of surface.
export function bars(ctx, x, y, edges, counts, colour) {
  ctx.fillStyle = colour;
  for (let k = 0; k < counts.length; k++) {
    if (!(counts[k] > 0)) continue;
    const a = x(edges[k]) + 1, b = x(edges[k + 1]) - 1, top = y(counts[k]), base = y(y.domain0 ?? 0);
    if (b - a < 0.5) continue;
    ctx.fillRect(a, top, b - a, Math.max(1, base - top));
  }
}

export function dots(ctx, x, y, xs, ys, colour, r = 3) {
  ctx.fillStyle = colour;
  for (let i = 0; i < xs.length; i++) {
    if (!Number.isFinite(ys[i])) continue;
    ctx.beginPath();
    ctx.arc(x(xs[i]), y(ys[i]), r, 0, 2 * Math.PI);
    ctx.fill();
  }
}
