// Two-dimensional diagrams of relations between galaxies, on chart paper.
//
// HubbleDiagram (Local Volume): distance against velocity in the Local Group
//   frame. Galaxies that follow cosmic expansion lie near the line v = H0 d;
//   members of bound groups depart from it.
// MassLightDiagram (groups): virial mass against K-band luminosity of groups
//   with five or more members, with lines of constant mass-to-light ratio.

import { PALETTE, rgba, ramp, TIDAL_STOPS, GROUP_STOPS, hexToRgb } from '../render/colour.js';
import { LabelPlacer } from '../render/declutter.js';

const SERIF = "'Source Serif 4', Georgia, 'Times New Roman', serif";
export const HUBBLE_CONSTANT = 73;   // km/s/Mpc

// A colour from a ramp, darkened towards the ink so that it reads on paper.
function onPaper(rgb, t = 0.35) {
  const ink = hexToRgb(PALETTE.ink);
  return `rgb(${rgb.map((x, i) => Math.round(255 * (x + (ink[i] - x) * t))).join(',')})`;
}

const niceStep = (span, target) => {
  const raw = span / target, p = 10 ** Math.floor(Math.log10(raw)), f = raw / p;
  return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p;
};

class Chart {
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.hits = [];
    this.hover = null;
    this.selected = null;
  }

  resize(width, height) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    this.canvas.width = Math.round(width * dpr);
    this.canvas.height = Math.round(height * dpr);
    this.dpr = dpr;
    this.width = width;
    this.height = height;
  }

  frame(xRange, yRange) {
    const { ctx } = this;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.fillStyle = PALETTE.paper;
    ctx.fillRect(0, 0, this.width, this.height);
    this.box = { left: 84, right: this.width - 36, top: 40, bottom: this.height - 70 };
    this.xRange = xRange;
    this.yRange = yRange;
    this.hits = [];
  }

  x(v) { const b = this.box, [a, z] = this.xRange; return b.left + ((v - a) / (z - a)) * (b.right - b.left); }
  y(v) { const b = this.box, [a, z] = this.yRange; return b.bottom - ((v - a) / (z - a)) * (b.bottom - b.top); }

  axes(xLabel, yLabel, xFormat = (v) => v, yFormat = (v) => v) {
    const { ctx, box } = this;
    ctx.save();
    ctx.strokeStyle = rgba(PALETTE.inkSoft, 0.25);
    ctx.fillStyle = PALETTE.inkSoft;
    ctx.font = `12px ${SERIF}`;
    ctx.lineWidth = 0.8;
    const xs = niceStep(this.xRange[1] - this.xRange[0], 8), ys = niceStep(this.yRange[1] - this.yRange[0], 7);
    ctx.textAlign = 'center';
    for (let v = Math.ceil(this.xRange[0] / xs) * xs; v <= this.xRange[1] + 1e-9; v += xs) {
      ctx.beginPath(); ctx.moveTo(this.x(v), box.top); ctx.lineTo(this.x(v), box.bottom); ctx.stroke();
      ctx.fillText(String(xFormat(+v.toFixed(6))), this.x(v), box.bottom + 18);
    }
    ctx.textAlign = 'right';
    for (let v = Math.ceil(this.yRange[0] / ys) * ys; v <= this.yRange[1] + 1e-9; v += ys) {
      ctx.beginPath(); ctx.moveTo(box.left, this.y(v)); ctx.lineTo(box.right, this.y(v)); ctx.stroke();
      ctx.fillText(String(yFormat(+v.toFixed(6))), box.left - 8, this.y(v) + 4);
    }
    ctx.strokeStyle = PALETTE.ink;
    ctx.lineWidth = 1.2;
    ctx.strokeRect(box.left, box.top, box.right - box.left, box.bottom - box.top);
    ctx.fillStyle = PALETTE.ink;
    ctx.font = `13px ${SERIF}`;
    ctx.textAlign = 'center';
    ctx.fillText(xLabel, (box.left + box.right) / 2, box.bottom + 42);
    ctx.translate(box.left - 58, (box.top + box.bottom) / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(yLabel, 0, 0);
    ctx.restore();
  }

  dot(x, y, r, fill, hit) {
    const { ctx } = this;
    ctx.beginPath();
    ctx.arc(x, y, r, 0, 2 * Math.PI);
    ctx.fillStyle = fill;
    ctx.fill();
    ctx.lineWidth = 0.6;
    ctx.strokeStyle = rgba(PALETTE.ink, 0.55);
    ctx.stroke();
    this.hits.push({ ...hit, x, y });
  }

  ring(x, y) {
    const { ctx } = this;
    ctx.beginPath();
    ctx.arc(x, y, 10, 0, 2 * Math.PI);
    ctx.strokeStyle = PALETTE.accentDeep;
    ctx.lineWidth = 1.5;
    ctx.stroke();
  }

  caption(text) {
    const { ctx } = this;
    ctx.save();
    ctx.fillStyle = PALETTE.inkSoft;
    ctx.font = `13px ${SERIF}`;
    ctx.textBaseline = 'bottom';
    ctx.fillText(text, 16, this.height - 12);
    ctx.restore();
  }

  tooltip(text, x, y) {
    const { ctx } = this;
    ctx.save();
    ctx.font = `13px ${SERIF}`;
    const w = ctx.measureText(text).width + 16;
    const tx = Math.min(x + 12, this.width - w - 6), ty = y + 12;
    ctx.fillStyle = rgba(PALETTE.ink, 0.94);
    ctx.fillRect(tx, ty, w, 24);
    ctx.fillStyle = PALETTE.paper;
    ctx.textBaseline = 'middle';
    ctx.fillText(text, tx + 8, ty + 12);
    ctx.restore();
  }

  hitTest(px, py) {
    let best = null, bestD = 10;
    for (const h of this.hits) {
      const d = Math.hypot(h.x - px, h.y - py);
      if (d < bestD) { best = h; bestD = d; }
    }
    return best;
  }
}

export class HubbleDiagram extends Chart {
  draw(cat, relation) {
    const maxDist = 12;
    const shown = [];
    let lo = 0, hi = 0, beyond = 0;
    for (let i = 0; i < cat.count; i++) {
      const d = cat.dist[i], v = cat.vLG[i];
      if (!Number.isFinite(d) || !Number.isFinite(v)) continue;
      if (d > maxDist) { beyond++; continue; }
      shown.push(i);
      lo = Math.min(lo, v);
      hi = Math.max(hi, v);
    }
    this.frame([0, maxDist], [Math.floor((lo - 50) / 100) * 100, Math.ceil((hi + 50) / 100) * 100]);
    this.axes('Distance from the Milky Way (Mpc)', 'Velocity relative to the Local Group (km/s)');
    const { ctx } = this;
    ctx.save();
    ctx.beginPath();
    ctx.rect(this.box.left, this.box.top, this.box.right - this.box.left, this.box.bottom - this.box.top);
    ctx.clip();
    ctx.strokeStyle = rgba(PALETTE.roseDeep, 0.7);
    ctx.setLineDash([6, 5]);
    ctx.lineWidth = 1.3;
    ctx.beginPath();
    ctx.moveTo(this.x(0), this.y(0));
    ctx.lineTo(this.x(maxDist), this.y(HUBBLE_CONSTANT * maxDist));
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.strokeStyle = rgba(PALETTE.ink, 0.35);
    ctx.beginPath();
    ctx.moveTo(this.x(0), this.y(0));
    ctx.lineTo(this.x(maxDist), this.y(0));
    ctx.stroke();
    // Field galaxies first, so that group members are drawn on top.
    shown.sort((a, b) => (relation.theta1[a] || -9) - (relation.theta1[b] || -9));
    for (const i of shown) {
      const L = cat.logL[i];
      const r = Number.isFinite(L) ? Math.min(5.5, Math.max(1.8, 1.8 + (L - 7) * 0.8)) : 2;
      this.dot(this.x(cat.dist[i]), this.y(cat.vLG[i]), r, onPaper(ramp(TIDAL_STOPS, relation.theta1[i])),
        { kind: 'galaxy', index: i, text: cat.names[i] });
    }
    ctx.restore();
    this.legendTidal();
    if (this.selected != null) {
      const h = this.hits.find((x) => x.index === this.selected);
      if (h) this.ring(h.x, h.y);
    }
    this.caption(`${shown.length} galaxies with velocities${beyond ? ` · ${beyond} beyond ${maxDist} Mpc not shown` : ''}`
      + ' · colour: tidal index Θ₁');
    if (this.hover) this.tooltip(this.hover.text, this.hover.x, this.hover.y);
  }

  legendTidal() {
    const { ctx } = this;
    const x0 = this.box.left + 16, y0 = this.box.top + 18;
    ctx.save();
    ctx.fillStyle = `rgba(${[251, 246, 250].join(',')},0.88)`;
    ctx.fillRect(x0 - 8, y0 - 12, 250, 4 * 18 + 6);
    ctx.font = `12px ${SERIF}`;
    ctx.textBaseline = 'middle';
    [[-1.5, 'field galaxy, Θ₁ < 0'], [1, 'group member'], [3, 'close companion, Θ₁ > 2']].forEach(([t, text], k) => {
      ctx.beginPath();
      ctx.arc(x0 + 5, y0 + k * 18, 4.5, 0, 2 * Math.PI);
      ctx.fillStyle = onPaper(ramp(TIDAL_STOPS, t));
      ctx.fill();
      ctx.fillStyle = PALETTE.inkSoft;
      ctx.fillText(text, x0 + 16, y0 + k * 18);
    });
    const yl = y0 + 3 * 18;
    ctx.strokeStyle = rgba(PALETTE.roseDeep, 0.8);
    ctx.setLineDash([6, 5]);
    ctx.lineWidth = 1.3;
    ctx.beginPath();
    ctx.moveTo(x0 - 2, yl);
    ctx.lineTo(x0 + 12, yl);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillText(`Hubble flow, v = H₀d, H₀ = ${HUBBLE_CONSTANT} km s⁻¹ Mpc⁻¹`, x0 + 16, yl);
    ctx.restore();
  }
}

export class MassLightDiagram extends Chart {
  draw(cat) {
    const gr = cat.groups;
    const shown = [];
    for (let j = 0; j < gr.count; j++) {
      if (gr.nHere[j] >= 5 && gr.mVir[j] > 0 && Number.isFinite(gr.logLK[j])) shown.push(j);
    }
    const xs = shown.map((j) => gr.logLK[j]), ys = shown.map((j) => Math.log10(gr.mVir[j] * 1e12));
    const xr = [Math.floor(Math.min(...xs) * 2) / 2, Math.ceil(Math.max(...xs) * 2) / 2];
    const yr = [Math.floor(Math.min(...ys) * 2) / 2, Math.ceil(Math.max(...ys) * 2) / 2];
    this.frame(xr, yr);
    this.axes('Group luminosity, log L_K (solar)', 'Virial mass, log M (solar)');
    const { ctx } = this;
    ctx.save();
    ctx.beginPath();
    ctx.rect(this.box.left, this.box.top, this.box.right - this.box.left, this.box.bottom - this.box.top);
    ctx.clip();
    ctx.setLineDash([5, 5]);
    ctx.lineWidth = 1;
    ctx.font = `italic 12px ${SERIF}`;
    ctx.textAlign = 'right';
    for (const ml of [10, 30, 100, 300, 1000]) {
      const k = Math.log10(ml);
      ctx.strokeStyle = rgba(PALETTE.roseDeep, 0.45);
      ctx.beginPath();
      ctx.moveTo(this.x(xr[0]), this.y(xr[0] + k));
      ctx.lineTo(this.x(xr[1]), this.y(xr[1] + k));
      ctx.stroke();
      // Label where the line leaves the frame, through the right edge or the top.
      const xExit = Math.min(xr[1], yr[1] - k);
      if (xExit > xr[0] + 0.3) {
        ctx.fillStyle = PALETTE.roseDeep;
        ctx.fillText(`M/L = ${ml}`, this.x(xExit) - 6, this.y(xExit + k) + 16);
      }
    }
    ctx.textAlign = 'left';
    ctx.setLineDash([]);
    shown.sort((a, b) => gr.nHere[a] - gr.nHere[b]);
    for (const j of shown) {
      const r = Math.min(7, 1.6 + Math.log10(gr.nHere[j]) * 2.6);
      this.dot(this.x(gr.logLK[j]), this.y(Math.log10(gr.mVir[j] * 1e12)), r,
        onPaper(ramp(GROUP_STOPS, Math.log10(gr.nHere[j]))),
        { kind: 'group', group: j, index: gr.brightest[j] >= 0 ? gr.brightest[j] : null,
          text: `${cat.groupName(j) || `Group ${gr.nest[j]}`} · ${gr.nHere[j]} galaxies` });
    }
    ctx.restore();
    // Cluster names, richest first, without overlapping one another; a paper
    // halo keeps them legible over the points.
    const placer = new LabelPlacer(ctx);
    ctx.save();
    ctx.font = `italic 13px ${SERIF}`;
    ctx.lineJoin = 'round';
    const clusters = this.hits.filter((h) => (gr.names.get(h.group) || '').endsWith('Cluster'))
      .sort((p, q) => gr.nHere[q.group] - gr.nHere[p.group]);
    for (const h of clusters) placer.reserve(h.x, h.y, 5);
    const halo = { fillText: (text, x, y) => {
      ctx.strokeStyle = 'rgba(251,246,250,0.92)';
      ctx.lineWidth = 3.5;
      ctx.strokeText(text, x, y);
      ctx.fillStyle = PALETTE.accentDeep;
      ctx.fillText(text, x, y);
    }, measureText: (t) => ctx.measureText(t) };
    placer.ctx = halo;
    for (const h of clusters) placer.place(gr.names.get(h.group), h.x, h.y, 6);
    ctx.restore();
    if (this.selected != null) {
      const h = this.hits.find((x) => x.index === this.selected);
      if (h) this.ring(h.x, h.y);
    }
    this.caption(`${shown.length} groups with five or more members (Tully 2015) · masses from velocity dispersions`);
    if (this.hover) this.tooltip(this.hover.text, this.hover.x, this.hover.y);
  }
}
