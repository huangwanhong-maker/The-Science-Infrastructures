// The timeline: a span of time in UTC, the instants at which the shown data
// exist, the present, and the time shown, which can be dragged.

const DAY = 86400000;
const WEEKDAY = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'];
const MONTH = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];

export class Timeline {
  constructor(canvas, onChange) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.onChange = onChange;
    this.span = [Date.now() - DAY, Date.now() + DAY];
    this.steps = [];
    this.coverage = [];
    this.time = Date.now();
    this.now = Date.now();
    new ResizeObserver(() => this.resize()).observe(canvas);
    const at = (e) => {
      const r = canvas.getBoundingClientRect();
      let t = this.timeAt(e.clientX - r.left);
      // Snap to a data time within 5 pixels.
      for (const s of this.steps) if (Math.abs(this.x(s) - (e.clientX - r.left)) < 5) { t = s; break; }
      this.onChange(Math.max(this.span[0], Math.min(this.span[1], t)));
    };
    let dragging = false;
    canvas.addEventListener('pointerdown', (e) => { dragging = true; canvas.setPointerCapture(e.pointerId); at(e); });
    canvas.addEventListener('pointermove', (e) => { if (dragging) at(e); });
    canvas.addEventListener('pointerup', () => { dragging = false; });
    canvas.addEventListener('pointercancel', () => { dragging = false; });
    this.resize();
  }

  resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const r = this.canvas.getBoundingClientRect();
    this.w = Math.max(10, r.width); this.h = Math.max(10, r.height);
    this.canvas.width = Math.round(this.w * dpr); this.canvas.height = Math.round(this.h * dpr);
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.draw();
  }

  pad() { return 14; }
  x(t) { const p = this.pad(); return p + ((t - this.span[0]) / (this.span[1] - this.span[0])) * (this.w - 2 * p); }
  timeAt(x) { const p = this.pad(); return this.span[0] + ((x - p) / (this.w - 2 * p)) * (this.span[1] - this.span[0]); }

  // span: [start, end] ms; steps: data times; coverage: [[start, end], ...] where data exist.
  set({ span, steps, coverage, time, now }) {
    if (span) this.span = span;
    if (steps) this.steps = steps;
    if (coverage) this.coverage = coverage;
    if (time !== undefined) this.time = time;
    if (now !== undefined) this.now = now;
    this.draw();
  }

  draw() {
    const c = this.ctx, w = this.w, h = this.h;
    if (!w) return;
    c.clearRect(0, 0, w, h);
    const y = 30;
    const [t0, t1] = this.span;
    const length = t1 - t0;
    c.font = '11px "Source Serif 4", Georgia, serif';
    c.textBaseline = 'alphabetic';
    // The track, and where the shown data exist.
    c.strokeStyle = 'rgba(240, 226, 238, 0.16)';
    c.lineWidth = 1;
    c.beginPath(); c.moveTo(this.pad(), y + 0.5); c.lineTo(w - this.pad(), y + 0.5); c.stroke();
    c.strokeStyle = 'rgba(202, 141, 190, 0.55)';
    c.lineWidth = 3;
    c.lineCap = 'round';
    for (const [a, b] of this.coverage) {
      const xa = this.x(Math.max(t0, a)), xb = this.x(Math.min(t1, b));
      if (xb >= xa) { c.beginPath(); c.moveTo(xa, y + 0.5); c.lineTo(Math.max(xb, xa + 0.5), y + 0.5); c.stroke(); }
    }
    c.lineCap = 'butt';
    // Data times.
    c.fillStyle = 'rgba(246, 238, 245, 0.55)';
    if (this.steps.length < (w - 2 * this.pad()) / 3) {
      for (const s of this.steps) if (s >= t0 && s <= t1) c.fillRect(Math.round(this.x(s)) - 0.5, y + 4, 1, 4);
    }
    // Calendar: days, or months for long spans.
    c.fillStyle = 'rgba(205, 182, 200, 0.85)';
    c.strokeStyle = 'rgba(240, 226, 238, 0.22)';
    c.lineWidth = 1;
    if (length <= 40 * DAY) {
      const pxPerDay = ((w - 2 * this.pad()) * DAY) / length;
      for (let d = Math.ceil(t0 / DAY) * DAY; d <= t1; d += DAY) {
        const x = Math.round(this.x(d)) + 0.5;
        c.beginPath(); c.moveTo(x, y - 6); c.lineTo(x, y + 6); c.stroke();
        const date = new Date(d);
        const text = pxPerDay > 70 ? `${WEEKDAY[date.getUTCDay()]} ${date.getUTCDate()}` : `${date.getUTCDate()}`;
        if (pxPerDay > 18) c.fillText(text, x + 4, y - 9);
      }
      if (pxPerDay > 220) {
        c.fillStyle = 'rgba(154, 132, 150, 0.9)';
        for (let t = Math.ceil(t0 / (6 * 3600000)) * 6 * 3600000; t <= t1; t += 6 * 3600000) {
          if (t % DAY === 0) continue;
          c.fillText(`${String(new Date(t).getUTCHours()).padStart(2, '0')}`, Math.round(this.x(t)) - 6, y + 17);
        }
      }
    } else {
      const start = new Date(t0);
      for (let m = Date.UTC(start.getUTCFullYear(), start.getUTCMonth() + 1, 1); m <= t1;) {
        const d = new Date(m), x = Math.round(this.x(m)) + 0.5;
        c.beginPath(); c.moveTo(x, y - 6); c.lineTo(x, y + 6); c.stroke();
        c.fillText(d.getUTCMonth() === 0 ? `${d.getUTCFullYear()}` : MONTH[d.getUTCMonth()], x + 4, y - 9);
        m = Date.UTC(d.getUTCFullYear(), d.getUTCMonth() + 1, 1);
      }
    }
    // The present.
    if (this.now >= t0 && this.now <= t1) {
      const x = this.x(this.now);
      c.fillStyle = '#CA8DBE';
      c.beginPath(); c.moveTo(x - 4, y + 15); c.lineTo(x + 4, y + 15); c.lineTo(x, y + 10); c.closePath(); c.fill();
    }
    // The time shown.
    const x = this.x(Math.max(t0, Math.min(t1, this.time)));
    c.strokeStyle = '#f4dbe7';
    c.lineWidth = 1.5;
    c.beginPath(); c.moveTo(x, 4); c.lineTo(x, h - 4); c.stroke();
    c.fillStyle = '#f4dbe7';
    c.beginPath(); c.arc(x, y + 0.5, 4.5, 0, 2 * Math.PI); c.fill();
    c.fillStyle = '#1f1721';
    c.beginPath(); c.arc(x, y + 0.5, 2, 0, 2 * Math.PI); c.fill();
  }
}

export function formatUTC(t, withSeconds = false) {
  const d = new Date(t);
  const hh = String(d.getUTCHours()).padStart(2, '0'), mm = String(d.getUTCMinutes()).padStart(2, '0');
  const ss = withSeconds ? `:${String(d.getUTCSeconds()).padStart(2, '0')}` : '';
  return `${WEEKDAY[d.getUTCDay()]} ${d.getUTCDate()} ${MONTH[d.getUTCMonth()]} ${d.getUTCFullYear()}, ${hh}:${mm}${ss} UTC`;
}

export function formatShortUTC(t) {
  const d = new Date(t);
  return `${d.getUTCDate()} ${MONTH[d.getUTCMonth()]} ${String(d.getUTCHours()).padStart(2, '0')}:${String(d.getUTCMinutes()).padStart(2, '0')} UTC`;
}

export function formatDay(t) {
  const d = new Date(t);
  return `${d.getUTCDate()} ${MONTH[d.getUTCMonth()]} ${d.getUTCFullYear()}`;
}

// "2 days before the time shown", "5 h after"
export function formatOffset(ms) {
  const a = Math.abs(ms), h = a / 3600000;
  const text = h >= 36 ? `${Math.round(h / 24)} days` : h >= 1.5 ? `${Math.round(h)} h` : `${Math.round(a / 60000)} min`;
  return `${text} ${ms > 0 ? 'before' : 'after'} the time shown`;
}
