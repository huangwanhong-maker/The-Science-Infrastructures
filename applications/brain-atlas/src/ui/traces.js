// The trace viewer of the Signals mode: a recording's channels drawn as an
// electroencephalographer reads them, a page of seconds at a time, stacked or laid over
// one another (the butterfly), with the global field power of the whole recording above
// as a map of where to look, its annotations, the page and the cursor.

const REGION_COLOURS = { frontal: '#e7a4b4', central: '#c3a9e6', temporal: '#9fc7e8', parietal: '#9fd6b4', occipital: '#f2c27c', other: '#d9c8d4' };

// The region of a 10-10 name, by its letters.
export function regionOf(name) {
  if (!name) return 'other';
  const l = name.replace(/[0-9z]+$/i, '').toUpperCase();
  if (['FP', 'AF', 'F', 'NZ'].includes(l)) return 'frontal';
  if (['FC', 'C'].includes(l)) return 'central';
  if (['FT', 'T', 'TP', 'A', 'M'].includes(l)) return 'temporal';
  if (['CP', 'P'].includes(l)) return 'parietal';
  if (['PO', 'O', 'I', 'IZ'].includes(l)) return 'occipital';
  return 'other';
}
export const regionColour = (name) => REGION_COLOURS[regionOf(name)];

// The pages of channels shown stacked: at most sixteen, as a reader of EEG takes them; a
// recording without these names shows its first sixteen.
const list = (names) => { const set = new Set(names); return (n) => set.has(n); };
export const MONTAGES = {
  butterfly: { name: 'Butterfly', pick: () => true },
  midline: { name: 'Midline', pick: list(['Fpz', 'AFz', 'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'POz', 'Oz', 'Iz']) },
  frontal: { name: 'Frontal', pick: list(['Fp1', 'Fpz', 'Fp2', 'AF7', 'AF3', 'AFz', 'AF4', 'AF8', 'F7', 'F3', 'Fz', 'F4', 'F8']) },
  central: { name: 'Central', pick: list(['FC5', 'FC3', 'FC1', 'FC2', 'FC4', 'FC6', 'T7', 'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6', 'T8']) },
  posterior: { name: 'Posterior', pick: list(['CP5', 'CP3', 'CPz', 'CP4', 'CP6', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO7', 'POz', 'PO8', 'O1', 'Oz', 'O2']) },
};

const fmt = (t) => `${t.toFixed(2)} s`;

export class TraceView {
  constructor(el, { onTime = () => {}, onControl = () => {} } = {}) {
    this.el = el;
    this.onTime = onTime;
    this.onControl = onControl;
    this.window = 10;
    this.start = 0;
    this.time = 0;
    this.montage = 'posterior';
    this.gain = 1;
    el.innerHTML = `<div class="dock-bar">
<button type="button" class="play" data-c="play" aria-label="Play"><svg viewBox="0 0 16 16"><path d="M4 2.5v11l9-5.5z" fill="currentColor"/></svg></button>
<span class="clock">0.00 s</span>
<span class="sep"></span>
<label>Speed <select data-c="speed"><option value="0.1">slow motion</option><option value="0.25">¼</option><option value="0.5">½</option><option value="1" selected>real time</option><option value="2">2×</option></select></label>
<label>Page <select data-c="window"><option value="2">2 s</option><option value="5">5 s</option><option value="10" selected>10 s</option><option value="20">20 s</option></select></label>
<label>Show <select data-c="montage">${Object.entries(MONTAGES).map(([k, m]) => `<option value="${k}"${k === this.montage ? ' selected' : ''}>${m.name}</option>`).join('')}</select></label>
<span class="gain"><button type="button" data-c="gain-" aria-label="Smaller">−</button><button type="button" data-c="gain+" aria-label="Larger">+</button></span>
<span class="band-name"></span>
<span class="dock-note"></span>
</div>
<canvas class="overview"></canvas>
<canvas class="traces"></canvas>`;
    this.overview = el.querySelector('canvas.overview');
    this.canvas = el.querySelector('canvas.traces');
    this.bind();
    new ResizeObserver(() => this.render()).observe(el);
  }

  bind() {
    this.el.querySelector('.dock-bar').addEventListener('click', (e) => {
      const b = e.target.closest('[data-c]');
      if (!b || b.tagName === 'SELECT') return;
      if (b.dataset.c === 'play') this.onControl('play');
      if (b.dataset.c === 'gain+') { this.gain *= 1.4; this.render(); }
      if (b.dataset.c === 'gain-') { this.gain /= 1.4; this.render(); }
    });
    this.el.querySelector('.dock-bar').addEventListener('change', (e) => {
      const s = e.target.closest('select');
      if (!s) return;
      if (s.dataset.c === 'speed') this.onControl('speed', Number(s.value));
      if (s.dataset.c === 'window') { this.window = Number(s.value); this.follow(this.time, true); this.render(); }
      if (s.dataset.c === 'montage') { this.montage = s.value; this.render(); }
    });
    const timeAt = (canvas, e, whole) => {
      const r = canvas.getBoundingClientRect();
      const f = Math.max(0, Math.min(1, (e.clientX - r.left - this.inset) / (r.width - this.inset - 6)));
      return whole ? f * this.rec.duration : this.start + f * this.window;
    };
    let dragging = null;
    const down = (canvas, whole) => (e) => {
      if (!this.rec) return;
      dragging = { canvas, whole };
      canvas.setPointerCapture(e.pointerId);
      this.onTime(timeAt(canvas, e, whole));
    };
    this.overview.addEventListener('pointerdown', down(this.overview, true));
    this.canvas.addEventListener('pointerdown', down(this.canvas, false));
    const move = (e) => { if (dragging) this.onTime(timeAt(dragging.canvas, e, dragging.whole)); };
    this.overview.addEventListener('pointermove', move);
    this.canvas.addEventListener('pointermove', move);
    const up = () => { dragging = null; };
    this.overview.addEventListener('pointerup', up);
    this.canvas.addEventListener('pointerup', up);
  }

  // rec: {rate, duration, channels: [{name, data}], annotations, meaning}; gfp: Float32Array.
  setRecording(rec, gfp, { band = '', note = '' } = {}) {
    this.rec = rec;
    this.gfp = gfp;
    // A robust scale: the 98th percentile of |x| over every channel, sampled.
    const sample = [];
    for (const c of rec.channels) for (let t = 0; t < c.data.length; t += 7) sample.push(Math.abs(c.data[t]));
    sample.sort((a, b) => a - b);
    this.scale = sample[Math.floor(0.98 * (sample.length - 1))] || 1;
    let g = 0;
    for (let t = 0; t < gfp.length; t++) g = Math.max(g, gfp[t]);
    this.gfpMax = g || 1;
    this.el.querySelector('.band-name').textContent = band;
    this.el.querySelector('.dock-note').textContent = note;
    this.follow(this.time, true);
    this.render();
  }

  setPlaying(on) {
    this.el.querySelector('.play').innerHTML = on
      ? '<svg viewBox="0 0 16 16"><path d="M4 3h3v10H4zM9 3h3v10H9z" fill="currentColor"/></svg>'
      : '<svg viewBox="0 0 16 16"><path d="M4 2.5v11l9-5.5z" fill="currentColor"/></svg>';
  }

  // Keep the cursor on the page: turn the page when it runs off.
  follow(t, centre = false) {
    if (!this.rec) return;
    const d = this.rec.duration, w = Math.min(this.window, d);
    if (centre) this.start = Math.max(0, Math.min(d - w, t - w / 2));
    else if (t < this.start || t > this.start + w) this.start = Math.max(0, Math.min(d - w, t - 0.1 * w));
  }

  setTime(t) {
    this.time = t;
    this.follow(t);
    this.el.querySelector('.clock').textContent = `${fmt(t)} / ${fmt(this.rec ? this.rec.duration : 0)}`;
    this.render();
  }

  fit(canvas) {
    const r = canvas.getBoundingClientRect(), dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = Math.max(1, Math.round(r.width * dpr)), h = Math.max(1, Math.round(r.height * dpr));
    if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, W: r.width, H: r.height };
  }

  render() {
    if (!this.rec) return;
    this.inset = 46;
    this.renderOverview();
    this.renderTraces();
  }

  annotationColour(text) {
    const m = this.rec.meaning && this.rec.meaning[text];
    if (!m || m.rest) return null;
    return m.colour || 'rgba(242, 194, 124, 0.16)';
  }

  renderOverview() {
    const { ctx, W, H } = this.fit(this.overview);
    const rec = this.rec, x0 = this.inset, w = W - x0 - 6, d = rec.duration;
    ctx.clearRect(0, 0, W, H);
    ctx.font = 'italic 10.5px "Source Serif 4", Georgia, serif';
    ctx.fillStyle = 'rgba(246, 238, 245, 0.45)';
    ctx.textBaseline = 'middle';
    ctx.fillText('field', 6, H / 2);
    for (const a of rec.annotations) {
      const c = this.annotationColour(a.text);
      if (!c) continue;
      ctx.fillStyle = c;
      ctx.fillRect(x0 + (a.onset / d) * w, 2, Math.max(1, ((a.duration || 0.5) / d) * w), H - 4);
    }
    // The global field power, as a soft ridge.
    const g = this.gfp, n = g.length;
    ctx.beginPath();
    ctx.moveTo(x0, H - 3);
    for (let px = 0; px <= w; px++) {
      const i0 = Math.floor((px / w) * (n - 1)), i1 = Math.max(i0 + 1, Math.floor(((px + 1) / w) * (n - 1)));
      let m = 0;
      for (let i = i0; i < Math.min(n, i1); i++) m = Math.max(m, g[i]);
      ctx.lineTo(x0 + px, H - 3 - (m / this.gfpMax) * (H - 8));
    }
    ctx.lineTo(x0 + w, H - 3);
    ctx.closePath();
    const grad = ctx.createLinearGradient(0, 0, 0, H);
    grad.addColorStop(0, 'rgba(242, 194, 124, 0.75)');
    grad.addColorStop(1, 'rgba(202, 141, 190, 0.25)');
    ctx.fillStyle = grad;
    ctx.fill();
    // The page and the cursor.
    ctx.strokeStyle = 'rgba(246, 238, 245, 0.55)';
    ctx.lineWidth = 1;
    ctx.strokeRect(x0 + (this.start / d) * w + 0.5, 1.5, (Math.min(this.window, d) / d) * w, H - 3);
    ctx.fillStyle = '#f6c47e';
    ctx.fillRect(x0 + (this.time / d) * w - 1, 0, 2, H);
  }

  renderTraces() {
    const { ctx, W, H } = this.fit(this.canvas);
    const rec = this.rec, x0 = this.inset, w = W - x0 - 6;
    ctx.clearRect(0, 0, W, H);
    const t0 = this.start, t1 = Math.min(rec.duration, this.start + this.window);
    const i0 = Math.floor(t0 * rec.rate), i1 = Math.min(Math.ceil(t1 * rec.rate), rec.channels[0].data.length - 1);
    const xs = (i) => x0 + ((i / rec.rate - t0) / (t1 - t0)) * w;
    // Annotations on the page, labelled with their meaning.
    ctx.font = 'italic 11px "Source Serif 4", Georgia, serif';
    for (const a of rec.annotations) {
      const end = a.onset + (a.duration || 0);
      if (end < t0 || a.onset > t1) continue;
      const c = this.annotationColour(a.text);
      const xa = x0 + ((Math.max(a.onset, t0) - t0) / (t1 - t0)) * w, xb = x0 + ((Math.min(end, t1) - t0) / (t1 - t0)) * w;
      if (c) { ctx.fillStyle = c; ctx.fillRect(xa, 0, Math.max(1, xb - xa), H); }
      const m = rec.meaning && rec.meaning[a.text];
      ctx.fillStyle = 'rgba(246, 238, 245, 0.6)';
      ctx.textBaseline = 'top';
      if (a.onset >= t0) ctx.fillText(m ? m.name : a.text, xa + 4, 3);
    }
    // Seconds.
    ctx.strokeStyle = 'rgba(240, 226, 238, 0.07)';
    ctx.fillStyle = 'rgba(246, 238, 245, 0.35)';
    ctx.font = '10px "Source Serif 4", Georgia, serif';
    ctx.textBaseline = 'bottom';
    for (let s = Math.ceil(t0); s <= t1; s++) {
      const x = x0 + ((s - t0) / (t1 - t0)) * w;
      ctx.beginPath(); ctx.moveTo(x + 0.5, 14); ctx.lineTo(x + 0.5, H); ctx.stroke();
      ctx.fillText(`${s}`, x + 3, H - 2);
    }
    const pick = MONTAGES[this.montage].pick;
    let shown = rec.channels.filter((c) => c.name === 'probe' || pick(c.name));
    if (!shown.filter((c) => c.name !== 'probe').length) shown = rec.channels.slice(0, 16);
    const butterfly = this.montage === 'butterfly';
    const rows = butterfly ? 1 : Math.max(1, shown.length);
    const rowH = (H - 16) / rows;
    const amp = (butterfly ? (H - 24) / 2 : rowH * 0.9) / (this.scale * (butterfly ? 1.1 : 1.2)) * this.gain;
    ctx.lineWidth = 1;
    const step = Math.max(1, Math.floor((i1 - i0) / (w * 1.5)));
    shown.forEach((c, k) => {
      const y = 14 + (butterfly ? (H - 16) / 2 : rowH * (k + 0.5));
      ctx.strokeStyle = c.colour || regionColour(c.name);
      ctx.globalAlpha = butterfly ? 0.45 : 0.92;
      ctx.beginPath();
      for (let i = i0; i <= i1; i += step) {
        const x = xs(i), v = c.data[i];
        if (i === i0) ctx.moveTo(x, y - v * amp); else ctx.lineTo(x, y - v * amp);
      }
      ctx.stroke();
      ctx.globalAlpha = 1;
      if (!butterfly) {
        ctx.fillStyle = c.colour || regionColour(c.name);
        ctx.font = '600 10.5px "Source Serif 4", Georgia, serif';
        ctx.textBaseline = 'middle';
        ctx.fillText(c.name, 6, y);
      }
    });
    if (butterfly) {
      ctx.fillStyle = 'rgba(246, 238, 245, 0.55)';
      ctx.font = 'italic 10.5px "Source Serif 4", Georgia, serif';
      ctx.textBaseline = 'middle';
      ctx.fillText(`${shown.length} ch`, 6, H / 2);
    }
    // The cursor.
    if (this.time >= t0 && this.time <= t1) {
      const x = x0 + ((this.time - t0) / (t1 - t0)) * w;
      ctx.fillStyle = 'rgba(246, 196, 126, 0.95)';
      ctx.fillRect(x - 0.75, 0, 1.5, H);
    }
  }
}
