// A chart of the sky for an observer at an instant, in the manner of a printed
// star atlas: stereographic projection of the visible hemisphere with the
// zenith at the centre, north at the top and east at the left, which is how the
// sky looks when the chart is held overhead.

import { BODIES, BY_KEY } from '../core/bodies.js';
import { DEG, ECLIPTIC_TO_ICRS, apply, fromSpherical, refraction, dot, sub, scale, add, unit } from '../core/frames.js';
import { PALETTE, rgba } from '../render/colour.js';
import { LabelPlacer } from '../render/declutter.js';

const SERIF = "'Source Serif 4', Georgia, 'Times New Roman', serif";
const NIGHT_SKY = [255, 250, 253];
const DAY_SKY = [252, 226, 216];

const mix = (a, b, t) => a.map((x, i) => Math.round(x + (b[i] - x) * t));
const clamp = (x, lo, hi) => Math.min(hi, Math.max(lo, x));

export function skyState(sunAltDeg) {
  if (sunAltDeg > 0) return 'daylight';
  if (sunAltDeg > -6) return 'civil twilight';
  if (sunAltDeg > -12) return 'nautical twilight';
  if (sunAltDeg > -18) return 'astronomical twilight';
  return 'night';
}

export class SkyChartView {
  constructor({ canvas, stars, constellations }) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.stars = stars;
    this.constellations = constellations;
    this.options = {
      magLimit: 5.5, lines: true, names: true, nameLang: 'la',
      ecliptic: true, grid: true, refraction: true, starNames: true,
    };
    this.hits = [];
    this.hover = null;
    this.selectedKey = null;
    this.width = this.height = 0;
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

  // Stereographic projection of (alt, az); null below the horizon unless raw.
  project(alt, az, raw = false) {
    if (!raw && alt < -0.002) return null;
    const r = this.R * Math.tan((Math.PI / 2 - alt) / 2);
    return [this.cx - r * Math.sin(az), this.cy - r * Math.cos(az)];
  }

  apparent(alt) {
    return this.options.refraction ? alt + refraction(alt) : alt;
  }

  // Horizontal coordinates of an ICRS unit vector, without allocation.
  horizonOf(m, x, y, z) {
    const hx = m[0] * x + m[1] * y + m[2] * z;
    const hy = m[3] * x + m[4] * y + m[5] * z;
    const hz = m[6] * x + m[7] * y + m[8] * z;
    let az = Math.atan2(hy, -hx);
    if (az < 0) az += 2 * Math.PI;
    return [Math.asin(clamp(hz, -1, 1)), az];
  }

  draw(snap) {
    const { ctx, options } = this;
    const w = this.width, h = this.height;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.fillStyle = PALETTE.paper;
    ctx.fillRect(0, 0, w, h);
    this.R = Math.max(40, Math.min(w, h) / 2 - 46);
    this.cx = w / 2;
    this.cy = h / 2;
    this.hits = [];
    const M = snap.icrsToHorizon;
    const sunAlt = snap.bodies.sun.alt / DEG;
    const day = clamp((sunAlt + 12) / 18, 0, 1);

    const [r, g, b] = mix(NIGHT_SKY, DAY_SKY, day);
    ctx.beginPath();
    ctx.arc(this.cx, this.cy, this.R, 0, 2 * Math.PI);
    ctx.fillStyle = `rgb(${r},${g},${b})`;
    ctx.fill();

    ctx.save();
    ctx.beginPath();
    ctx.arc(this.cx, this.cy, this.R, 0, 2 * Math.PI);
    ctx.clip();
    if (options.grid) this.drawGrid();
    if (options.ecliptic) this.drawEcliptic(M);
    if (options.lines) this.drawConstellationLines(M);
    const starNames = this.drawStars(M, day);
    const bodyNames = this.drawBodies(snap);
    ctx.restore();

    // Labels, in order of importance, without overlaps; they may cross the rim.
    ctx.save();
    const placer = new LabelPlacer(ctx);
    for (const n of bodyNames) placer.reserve(n.x, n.y, n.r);
    ctx.fillStyle = PALETTE.accentDeep;
    for (const n of bodyNames) {
      ctx.font = `${n.selected ? '600 ' : ''}italic 13px ${SERIF}`;
      placer.place(n.text, n.x, n.y, n.r, { force: true });
    }
    ctx.font = `italic 12px ${SERIF}`;
    ctx.fillStyle = rgba(PALETTE.ink, 0.8);
    for (const n of starNames) placer.place(n.text, n.x, n.y, n.r);
    if (options.names) this.drawConstellationNames(M, placer);
    ctx.restore();

    this.drawRim();
    this.drawCaption(snap, sunAlt);
    if (this.hover) this.drawTooltip(this.hover);
  }

  drawGrid() {
    const { ctx } = this;
    ctx.save();
    ctx.strokeStyle = rgba(PALETTE.inkSoft, 0.3);
    ctx.lineWidth = 0.8;
    ctx.setLineDash([2, 4]);
    for (const alt of [30, 60]) {
      ctx.beginPath();
      ctx.arc(this.cx, this.cy, this.R * Math.tan((90 - alt) * DEG / 2), 0, 2 * Math.PI);
      ctx.stroke();
    }
    for (let az = 0; az < 360; az += 45) {
      const a = this.project(0, az * DEG), b = this.project(80 * DEG, az * DEG);
      ctx.beginPath();
      ctx.moveTo(a[0], a[1]);
      ctx.lineTo(b[0], b[1]);
      ctx.stroke();
    }
    ctx.setLineDash([]);
    ctx.beginPath();
    ctx.moveTo(this.cx - 5, this.cy); ctx.lineTo(this.cx + 5, this.cy);
    ctx.moveTo(this.cx, this.cy - 5); ctx.lineTo(this.cx, this.cy + 5);
    ctx.stroke();
    ctx.restore();
  }

  polyline(points) {
    const { ctx } = this;
    let drawing = false;
    ctx.beginPath();
    for (const p of points) {
      if (!p) { drawing = false; continue; }
      if (drawing) ctx.lineTo(p[0], p[1]); else ctx.moveTo(p[0], p[1]);
      drawing = true;
    }
    ctx.stroke();
  }

  drawEcliptic(M) {
    const points = [];
    for (let lon = 0; lon <= 360; lon += 1) {
      const v = apply(ECLIPTIC_TO_ICRS, fromSpherical(lon * DEG, 0));
      const [alt, az] = this.horizonOf(M, v[0], v[1], v[2]);
      points.push(this.project(this.apparent(alt), az));
    }
    const { ctx } = this;
    ctx.save();
    ctx.strokeStyle = rgba(PALETTE.roseDeep, 0.55);
    ctx.lineWidth = 1.1;
    ctx.setLineDash([6, 5]);
    this.polyline(points);
    ctx.restore();
  }

  drawConstellationLines(M) {
    const { ctx } = this;
    const s = this.constellations.segments;
    ctx.save();
    ctx.strokeStyle = rgba(PALETTE.ink, 0.28);
    ctx.lineWidth = 0.8;
    ctx.beginPath();
    for (let k = 0; k < s.length; k += 6) {
      const [a1, z1] = this.horizonOf(M, s[k], s[k + 1], s[k + 2]);
      const [a2, z2] = this.horizonOf(M, s[k + 3], s[k + 4], s[k + 5]);
      if (a1 < 0 && a2 < 0) continue;
      const p = this.project(this.apparent(a1), z1, true);
      const q = this.project(this.apparent(a2), z2, true);
      ctx.moveTo(p[0], p[1]);
      ctx.lineTo(q[0], q[1]);
    }
    ctx.stroke();
    ctx.restore();
  }

  // Draws the stars; returns the names of bright named stars for labelling.
  drawStars(M, day) {
    const { ctx, stars, options } = this;
    const limit = options.magLimit;
    const size = clamp(this.R / 330, 0.55, 1.15);   // smaller dots on smaller charts
    const ink = rgba(PALETTE.ink, (1 - 0.45 * day).toFixed(3));
    const halo = `rgba(${NIGHT_SKY.join(',')},0.9)`;
    const names = [];
    for (const i of stars.byBrightness) {
      const m = stars.mag[i];
      if (m > limit) break;
      const [alt0, az] = this.horizonOf(M, stars.unit[3 * i], stars.unit[3 * i + 1], stars.unit[3 * i + 2]);
      if (alt0 < 0) continue;
      const p = this.project(this.apparent(alt0), az);
      if (!p) continue;
      const radius = clamp(0.55 + (limit - m) * 0.44, 0.55, 4.4) * size;
      if (radius > 1.4) {
        ctx.beginPath();
        ctx.arc(p[0], p[1], radius + 1.1, 0, 2 * Math.PI);
        ctx.fillStyle = halo;
        ctx.fill();
      }
      ctx.beginPath();
      ctx.arc(p[0], p[1], radius, 0, 2 * Math.PI);
      ctx.fillStyle = ink;
      ctx.fill();
      if (m < 4.5) this.hits.push({ kind: 'star', index: i, x: p[0], y: p[1], r: radius });
      if (options.starNames && m < 1.6 && stars.starName(i)) names.push({ text: stars.starName(i), x: p[0], y: p[1], r: radius });
    }
    return names;
  }

  drawConstellationNames(M, placer) {
    const { ctx, options } = this;
    const lang = options.nameLang;
    const cjk = lang === 'zh' || lang === 'ja';
    ctx.save();
    ctx.fillStyle = rgba(PALETTE.inkSoft, 0.85);
    const visible = [];
    for (const c of this.constellations.labels) {
      const [alt, az] = this.horizonOf(M, c.unit[0], c.unit[1], c.unit[2]);
      if (alt < 0.05) continue;
      visible.push([c, this.project(this.apparent(alt), az)]);
    }
    visible.sort((a, b) => a[0].rank - b[0].rank);
    for (const [c, p] of visible) {
      const name = c.names[lang] || c.names.la;
      ctx.font = cjk ? `12px ${SERIF}` : `${c.rank === 1 ? 10.5 : 9.5}px ${SERIF}`;
      ctx.letterSpacing = cjk ? '0px' : '1px';
      placer.place(cjk ? name : name.toUpperCase(), p[0], p[1], 0, { centred: true });
    }
    ctx.restore();
  }

  // Draws the Sun, Moon and planets; returns their names for labelling.
  drawBodies(snap) {
    const { ctx } = this;
    const names = [];
    for (const b of BODIES) {
      if (b.key === 'earth') continue;
      const e = snap.bodies[b.key];
      const p = this.project(this.apparent(e.alt), e.az);
      if (!p) continue;
      if (b.key === 'sun') this.drawSun(p);
      else if (b.key === 'moon') this.drawMoon(snap, p);
      else {
        ctx.beginPath();
        ctx.arc(p[0], p[1], 4.6, 0, 2 * Math.PI);
        ctx.fillStyle = b.colour;
        ctx.fill();
        ctx.lineWidth = 1.2;
        ctx.strokeStyle = PALETTE.ink;
        ctx.stroke();
      }
      const selected = this.selectedKey === b.key;
      if (selected) {
        ctx.beginPath();
        ctx.arc(p[0], p[1], 12, 0, 2 * Math.PI);
        ctx.strokeStyle = PALETTE.accentDeep;
        ctx.lineWidth = 1.4;
        ctx.stroke();
      }
      const r = b.key === 'sun' ? 13 : b.key === 'moon' ? 8.5 : selected ? 12 : 5.5;
      names.push({ text: b.name, x: p[0], y: p[1], r, selected });
      this.hits.push({ kind: 'body', key: b.key, x: p[0], y: p[1], r: 9 });
    }
    return names;
  }

  drawSun(p) {
    const { ctx } = this;
    ctx.save();
    ctx.strokeStyle = PALETTE.goldDeep;
    ctx.lineWidth = 1.2;
    for (let k = 0; k < 12; k++) {
      const a = (k * Math.PI) / 6;
      ctx.beginPath();
      ctx.moveTo(p[0] + 9 * Math.cos(a), p[1] + 9 * Math.sin(a));
      ctx.lineTo(p[0] + 13 * Math.cos(a), p[1] + 13 * Math.sin(a));
      ctx.stroke();
    }
    ctx.beginPath();
    ctx.arc(p[0], p[1], 7, 0, 2 * Math.PI);
    ctx.fillStyle = PALETTE.gold;
    ctx.fill();
    ctx.stroke();
    ctx.restore();
  }

  // The Moon's phase: the bright limb faces the Sun along the great circle
  // through both; the terminator is a half-ellipse of semi-axis |cos i|.
  drawMoon(snap, p) {
    const { ctx } = this;
    const moon = snap.bodies.moon, sun = snap.bodies.sun;
    const m = unit(moon.geo), s = unit(sun.geo);
    const towardSun = unit(sub(s, scale(m, dot(m, s))));
    const q = add(m, scale(towardSun, 0.02));
    const [qa, qz] = this.horizonOf(snap.icrsToHorizon, q[0], q[1], q[2]);
    const pq = this.project(this.apparent(qa), qz, true);
    const angle = Math.atan2(pq[1] - p[1], pq[0] - p[0]);
    const radius = 7.5;
    ctx.save();
    ctx.translate(p[0], p[1]);
    ctx.rotate(angle);
    ctx.beginPath();
    ctx.arc(0, 0, radius, 0, 2 * Math.PI);
    ctx.fillStyle = rgba(PALETTE.ink, 0.85);
    ctx.fill();
    const ci = Math.cos(moon.phaseAngle);
    ctx.beginPath();
    for (let k = 0; k <= 24; k++) {
      const t = -Math.PI / 2 + (Math.PI * k) / 24;
      const x = radius * Math.cos(t), y = radius * Math.sin(t);
      if (k === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    for (let k = 24; k >= 0; k--) {
      const t = -Math.PI / 2 + (Math.PI * k) / 24;
      ctx.lineTo(-ci * radius * Math.cos(t), radius * Math.sin(t));
    }
    ctx.closePath();
    ctx.fillStyle = '#fbeacb';
    ctx.fill();
    ctx.beginPath();
    ctx.arc(0, 0, radius, 0, 2 * Math.PI);
    ctx.lineWidth = 1;
    ctx.strokeStyle = PALETTE.ink;
    ctx.stroke();
    ctx.restore();
  }

  drawRim() {
    const { ctx } = this;
    ctx.save();
    ctx.beginPath();
    ctx.arc(this.cx, this.cy, this.R, 0, 2 * Math.PI);
    ctx.lineWidth = 1.6;
    ctx.strokeStyle = PALETTE.ink;
    ctx.stroke();
    ctx.lineWidth = 0.8;
    for (let az = 0; az < 360; az += 5) {
      const len = az % 30 === 0 ? 8 : az % 10 === 0 ? 5 : 2.5;
      const a = az * DEG;
      const x = -Math.sin(a), y = -Math.cos(a);
      ctx.beginPath();
      ctx.moveTo(this.cx + x * this.R, this.cy + y * this.R);
      ctx.lineTo(this.cx + x * (this.R + len), this.cy + y * (this.R + len));
      ctx.stroke();
    }
    ctx.fillStyle = PALETTE.ink;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.font = `600 15px ${SERIF}`;
    for (const [letter, az] of [['N', 0], ['E', 90], ['S', 180], ['W', 270]]) {
      const a = az * DEG;
      ctx.fillText(letter, this.cx - Math.sin(a) * (this.R + 22), this.cy - Math.cos(a) * (this.R + 22));
    }
    ctx.restore();
  }

  drawCaption(snap, sunAlt) {
    const { ctx } = this;
    const o = snap.observer;
    const lat = `${Math.abs(o.lat / DEG).toFixed(2)}° ${o.lat >= 0 ? 'N' : 'S'}`;
    const lon = `${Math.abs(o.lon / DEG).toFixed(2)}° ${o.lon >= 0 ? 'E' : 'W'}`;
    ctx.save();
    ctx.fillStyle = PALETTE.inkSoft;
    ctx.font = `13px ${SERIF}`;
    ctx.textBaseline = 'top';
    ctx.fillText(`${o.name ? `${o.name} · ` : ''}${lat}, ${lon}`, 16, 14);
    if (o.timezone) {
      // Local civil time at the place, from its time zone (IANA), as clocks there show it.
      try {
        const local = new Intl.DateTimeFormat('en-GB', { timeZone: o.timezone, weekday: 'short', hour: '2-digit', minute: '2-digit', timeZoneName: 'short' })
          .format(new Date(snap.scales.unixMs));
        ctx.fillText(`${local} local time`, 16, 34);
      } catch { /* an unknown zone: no local time */ }
    }
    ctx.textBaseline = 'bottom';
    const where = sunAlt >= 0 ? 'above' : 'below';
    ctx.fillText(`Sun ${Math.abs(sunAlt).toFixed(1)}° ${where} the horizon · ${skyState(sunAlt)}`, 16, this.height - 14);
    ctx.restore();
  }

  drawTooltip(hit) {
    const { ctx } = this;
    const text = hit.kind === 'body' ? BY_KEY[hit.key].name
      : `${this.stars.label(hit.index) || 'Star'} · mag ${this.stars.mag[hit.index].toFixed(1)}`;
    ctx.save();
    ctx.font = `13px ${SERIF}`;
    const width = ctx.measureText(text).width + 16;
    const x = Math.min(hit.x + 12, this.width - width - 6), y = hit.y + 14;
    ctx.fillStyle = rgba(PALETTE.ink, 0.94);
    ctx.fillRect(x, y, width, 24);
    ctx.fillStyle = PALETTE.paper;
    ctx.textBaseline = 'middle';
    ctx.fillText(text, x + 8, y + 12);
    ctx.restore();
  }

  hitTest(x, y) {
    let best = null, bestD = 14;
    for (const h of this.hits) {
      const d = Math.hypot(h.x - x, h.y - y) - (h.kind === 'body' ? 4 : 0);
      if (d < bestD) { best = h; bestD = d; }
    }
    return best;
  }
}
