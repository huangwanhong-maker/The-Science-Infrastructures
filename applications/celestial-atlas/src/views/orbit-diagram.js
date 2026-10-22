// The planets' orbits seen from the north ecliptic pole, with the vernal
// equinox to the right; planets move anticlockwise. Two scales: the inner
// planets on a linear scale, or all planets with a square-root radial scale.

import { PLANETS, BY_KEY } from '../core/bodies.js';
import { AU_KM, ICRS_TO_ECLIPTIC, apply } from '../core/frames.js';
import { osculatingOrbit } from '../core/sky.js';
import { PALETTE, hexToRgb, rgba } from '../render/colour.js';
import { LabelPlacer } from '../render/declutter.js';

const SERIF = "'Source Serif 4', Georgia, 'Times New Roman', serif";
const ORBIT_REFRESH_S = 6 * 3600;
const RANGES = {
  inner: { extent: 1.75, rings: [0.5, 1, 1.5], keys: ['mercury', 'venus', 'earth', 'mars'] },
  all: { extent: 50, rings: [1, 5, 10, 20, 30, 40, 50], keys: PLANETS.map((b) => b.key) },
};

const ecliptic = (v) => {
  const e = apply(ICRS_TO_ECLIPTIC, v);
  return [e[0] / AU_KM, e[1] / AU_KM, e[2] / AU_KM];
};

// Body colours darkened towards ink, so that they read on paper.
const onPaper = (hex) => {
  const [r, g, b] = hexToRgb(hex);
  const ink = hexToRgb(PALETTE.ink);
  const t = 0.35;
  return `rgb(${[r, g, b].map((x, i) => Math.round(255 * (x + (ink[i] - x) * t))).join(',')})`;
};

export class OrbitDiagramView {
  constructor({ canvas }) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.options = { range: 'inner', sightLines: true };
    this.orbits = {};
    this.lastOrbitEt = -Infinity;
    this.hits = [];
    this.hover = null;
    this.selectedKey = null;
  }

  resize(width, height) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    this.canvas.width = Math.round(width * dpr);
    this.canvas.height = Math.round(height * dpr);
    this.dpr = dpr;
    this.width = width;
    this.height = height;
  }

  setOptions(options) {
    Object.assign(this.options, options);
  }

  // Screen position of an ecliptic point (AU), with the current radial scale.
  toScreen(x, y) {
    const range = RANGES[this.options.range];
    const r = Math.hypot(x, y);
    if (r === 0) return [this.cx, this.cy];
    const rs = this.options.range === 'all'
      ? this.R * Math.sqrt(r / range.extent)
      : this.R * (r / range.extent);
    return [this.cx + (rs * x) / r, this.cy - (rs * y) / r];
  }

  radiusOnScreen(r) {
    const range = RANGES[this.options.range];
    return this.options.range === 'all' ? this.R * Math.sqrt(r / range.extent) : this.R * (r / range.extent);
  }

  refreshOrbits(snap) {
    for (const b of PLANETS) {
      const pts = osculatingOrbit(snap, b.key, 360);
      const out = new Float64Array(2 * 360);
      for (let k = 0; k < 360; k++) {
        const e = ecliptic([pts[3 * k], pts[3 * k + 1], pts[3 * k + 2]]);
        out[2 * k] = e[0];
        out[2 * k + 1] = e[1];
      }
      this.orbits[b.key] = out;
    }
    this.lastOrbitEt = snap.et;
  }

  draw(snap) {
    const { ctx } = this;
    const w = this.width, h = this.height;
    if (Math.abs(snap.et - this.lastOrbitEt) > ORBIT_REFRESH_S) this.refreshOrbits(snap);
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.fillStyle = PALETTE.paper;
    ctx.fillRect(0, 0, w, h);
    this.cx = w / 2;
    this.cy = h / 2;
    this.R = Math.max(40, Math.min(w, h) / 2 - 40);
    this.hits = [];
    const range = RANGES[this.options.range];

    // Distance rings and the direction of the vernal equinox.
    ctx.save();
    ctx.strokeStyle = rgba(PALETTE.inkSoft, 0.3);
    ctx.fillStyle = rgba(PALETTE.inkSoft, 0.85);
    ctx.lineWidth = 0.8;
    ctx.setLineDash([2, 4]);
    ctx.font = `11px ${SERIF}`;
    for (const r of range.rings) {
      const rs = this.radiusOnScreen(r);
      ctx.beginPath();
      ctx.arc(this.cx, this.cy, rs, 0, 2 * Math.PI);
      ctx.stroke();
      ctx.fillText(`${r} AU`, this.cx + rs * 0.71 + 4, this.cy + rs * 0.71 + 12);
    }
    ctx.setLineDash([]);
    ctx.strokeStyle = rgba(PALETTE.roseDeep, 0.45);
    ctx.beginPath();
    ctx.moveTo(this.cx, this.cy);
    ctx.lineTo(this.cx + this.R + 18, this.cy);
    ctx.stroke();
    ctx.fillStyle = PALETTE.accentDeep;
    ctx.font = `13px ${SERIF}`;
    ctx.fillText('♈︎', this.cx + this.R + 22, this.cy + 4);
    ctx.restore();

    // Orbits.
    for (const key of range.keys) {
      const pts = this.orbits[key];
      ctx.beginPath();
      for (let k = 0; k <= 360; k++) {
        const j = k % 360;
        const [x, y] = this.toScreen(pts[2 * j], pts[2 * j + 1]);
        if (k === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.strokeStyle = onPaper(BY_KEY[key].colour);
      ctx.globalAlpha = 0.55;
      ctx.lineWidth = 1.1;
      ctx.stroke();
      ctx.globalAlpha = 1;
    }

    // Lines of sight from the Earth.
    const earth = ecliptic(snap.bodies.earth.helio.p);
    const earthScreen = this.toScreen(earth[0], earth[1]);
    if (this.options.sightLines) {
      ctx.save();
      ctx.setLineDash([3, 4]);
      ctx.lineWidth = 0.9;
      ctx.strokeStyle = rgba(PALETTE.ink, 0.35);
      for (const key of range.keys) {
        if (key === 'earth') continue;
        const p = ecliptic(snap.bodies[key].helio.p);
        const s = this.toScreen(p[0], p[1]);
        ctx.beginPath();
        ctx.moveTo(earthScreen[0], earthScreen[1]);
        ctx.lineTo(s[0], s[1]);
        ctx.stroke();
      }
      ctx.beginPath();
      ctx.moveTo(earthScreen[0], earthScreen[1]);
      ctx.lineTo(this.cx, this.cy);
      ctx.strokeStyle = rgba(PALETTE.goldDeep, 0.6);
      ctx.stroke();
      ctx.restore();
    }

    // The Sun.
    ctx.beginPath();
    ctx.arc(this.cx, this.cy, 8, 0, 2 * Math.PI);
    ctx.fillStyle = PALETTE.gold;
    ctx.fill();
    ctx.strokeStyle = PALETTE.goldDeep;
    ctx.lineWidth = 1.2;
    ctx.stroke();
    this.hits.push({ kind: 'body', key: 'sun', x: this.cx, y: this.cy });

    // Planets inside the range, and markers at the edge for those beyond it.
    const names = [];
    for (const b of PLANETS) {
      const p = ecliptic(snap.bodies[b.key].helio.p);
      const inside = range.keys.includes(b.key);
      let [x, y] = this.toScreen(p[0], p[1]);
      if (!inside) {
        const a = Math.atan2(p[1], p[0]);
        x = this.cx + (this.R + 6) * Math.cos(a);
        y = this.cy - (this.R + 6) * Math.sin(a);
      }
      ctx.beginPath();
      ctx.arc(x, y, inside ? 5.5 : 3, 0, 2 * Math.PI);
      ctx.fillStyle = b.colour;
      ctx.fill();
      ctx.lineWidth = 1.1;
      ctx.strokeStyle = PALETTE.ink;
      ctx.stroke();
      const selected = this.selectedKey === b.key;
      if (selected) {
        ctx.beginPath();
        ctx.arc(x, y, 11, 0, 2 * Math.PI);
        ctx.strokeStyle = PALETTE.accentDeep;
        ctx.lineWidth = 1.4;
        ctx.stroke();
      }
      names.push({ text: inside ? b.name : `${b.name} →`, x, y, r: selected ? 11 : inside ? 6 : 4, inside, selected });
      this.hits.push({ kind: 'body', key: b.key, x, y });
    }

    // Names without overlaps: the selected planet first, then those in range.
    const placer = new LabelPlacer(ctx);
    placer.reserve(this.cx, this.cy, 9);
    for (const n of names) placer.reserve(n.x, n.y, n.r);
    names.sort((a, b) => (b.selected - a.selected) || (b.inside - a.inside));
    for (const n of names) {
      ctx.font = `${n.selected ? '600 ' : ''}${n.inside ? '' : 'italic '}13px ${SERIF}`;
      ctx.fillStyle = n.inside ? PALETTE.ink : PALETTE.inkSoft;
      placer.place(n.text, n.x, n.y, n.r, { force: true });
    }

    ctx.fillStyle = PALETTE.inkSoft;
    ctx.font = `13px ${SERIF}`;
    ctx.textBaseline = 'bottom';
    ctx.fillText(this.options.range === 'all'
      ? 'Seen from the north ecliptic pole · radial scale ∝ √distance'
      : 'Seen from the north ecliptic pole · linear scale', 16, h - 14);
    ctx.textBaseline = 'alphabetic';
  }

  hitTest(x, y) {
    let best = null, bestD = 14;
    for (const h of this.hits) {
      const d = Math.hypot(h.x - x, h.y - y);
      if (d < bestD) { best = h; bestD = d; }
    }
    return best;
  }
}
