// The circular map of a genome, as an SVG string: positions run clockwise from the top.
//
//   r 482–505   the scale, in kb or Mb
//   r 452–472   features on the forward strand;  r 428–448  on the reverse strand
//   r 350–420   GC content against its mean: ochre above, blue below
//   r 290–346   GC skew (G − C)/(G + C): mulberry where G leads, sage where C does
//   r 160–270   the cumulative skew, a closed curve whose lowest point is the origin of
//               replication and highest the terminus (Grigoriev 1998)
//
// Functions of geometry are exported for the view's hit-testing: the angle of a position and
// the position of a point.

import { esc } from '../dom.js';
import { featureColour, token } from '../palette.js';
import { tickStep, tickLabel, length as lengthText, percent } from '../format.js';

const TAU = 2 * Math.PI;
export const angleOf = (p, L) => (TAU * p) / L - Math.PI / 2;
export function positionAt(x, y, L) {
  let a = Math.atan2(y, x) + Math.PI / 2;
  if (a < 0) a += TAU;
  return Math.min(L - 1, Math.max(0, Math.floor((a / TAU) * L)));
}
const pt = (r, a) => `${(r * Math.cos(a)).toFixed(2)},${(r * Math.sin(a)).toFixed(2)}`;

// An annular sector from angle a0 to a1 between radii r0 < r1.
export function sector(r0, r1, a0, a1) {
  if (a1 - a0 >= TAU - 1e-9) {
    // A whole ring: two halves.
    return `M${pt(r1, a0)}A${r1},${r1} 0 1 1 ${pt(r1, a0 + Math.PI)}A${r1},${r1} 0 1 1 ${pt(r1, a0)}M${pt(r0, a0)}A${r0},${r0} 0 1 0 ${pt(r0, a0 + Math.PI)}A${r0},${r0} 0 1 0 ${pt(r0, a0)}Z`;
  }
  const large = a1 - a0 > Math.PI ? 1 : 0;
  return `M${pt(r1, a0)}A${r1},${r1} 0 ${large} 1 ${pt(r1, a1)}L${pt(r0, a1)}A${r0},${r0} 0 ${large} 0 ${pt(r0, a0)}Z`;
}

export function circleSvg({ L, circular, tracks, features, window, selected, title, subtitle, markers = ['origin', 'terminus'], links = [] }) {
  const out = [];
  const ink = token('--ink'), faint = token('--ink-faint'), hair = token('--hair-2');
  // The scale.
  out.push(`<circle r="482" fill="none" stroke="${hair}" />`);
  const step = tickStep(L, 14);
  for (let p = 0; p < L; p += step) {
    if (p > 0 && L - p < step * 0.45) break;  // the last tick would crowd the first
    const a = angleOf(p, L);
    out.push(`<line x1="${(482 * Math.cos(a)).toFixed(1)}" y1="${(482 * Math.sin(a)).toFixed(1)}" x2="${(492 * Math.cos(a)).toFixed(1)}" y2="${(492 * Math.sin(a)).toFixed(1)}" stroke="${faint}" />`);
    const lx = 502 * Math.cos(a), ly = 502 * Math.sin(a);
    const anchor = Math.abs(Math.cos(a)) < 0.2 ? 'middle' : Math.cos(a) > 0 ? 'start' : 'end';
    out.push(`<text class="tick-label" x="${lx.toFixed(1)}" y="${(ly + 4).toFixed(1)}" text-anchor="${anchor}">${esc(tickLabel(p, step, L))}</text>`);
  }
  if (!circular) {
    // A linear genome's two ends meet at the top.
    out.push(`<line x1="0" y1="-478" x2="0" y2="-160" stroke="${ink}" stroke-dasharray="2 3" stroke-width="1" />`);
    out.push(`<text class="marker-label" x="6" y="-152" >ends</text>`);
  }
  // Strand guides.
  out.push(`<circle r="450" fill="none" stroke="${token('--hair')}" />`);
  // Features.
  for (const f of features) {
    const col = featureColour(f.key);
    const [r0, r1] = f.strand === '-' ? [428, 448] : [452, 472];
    for (const [b, e] of f.parts) {
      const a0 = angleOf(b, L), a1 = Math.max(angleOf(e, L), a0 + 0.0025);
      out.push(`<path class="feature${selected === f.i ? ' selected' : ''}" data-i="${f.i}" d="${sector(r0, r1, a0, a1)}" fill="${col}" fill-opacity="${f.pseudo ? 0.35 : 0.85}" />`);
    }
  }
  // GC against its mean, and the skew.
  const nb = tracks.gc.length;
  const finite = (xs) => xs.filter((v) => Number.isFinite(v));
  const dev = finite(tracks.gc.map((g) => g - tracks.mean_gc));
  const maxDev = Math.max(1e-6, ...dev.map(Math.abs));
  const maxSkew = Math.max(1e-6, ...finite(tracks.skew).map(Math.abs));
  const hi = token('--high'), lo = token('--low'), pos = token('--pos'), neg = token('--neg');
  const binAngle = (k) => angleOf(tracks.begin + ((tracks.end - tracks.begin) * k) / nb, L);
  out.push(`<circle r="385" fill="none" stroke="${token('--hair')}" />`);
  out.push(`<circle r="318" fill="none" stroke="${token('--hair')}" />`);
  for (let k = 0; k < nb; k++) {
    const a0 = binAngle(k), a1 = binAngle(k + 1) + 0.0004;
    const g = tracks.gc[k];
    if (Number.isFinite(g)) {
      const d = ((g - tracks.mean_gc) / maxDev) * 35;
      if (Math.abs(d) > 0.2) out.push(`<path d="${sector(d > 0 ? 385 : 385 + d, d > 0 ? 385 + d : 385, a0, a1)}" fill="${d > 0 ? hi : lo}" fill-opacity="0.8" />`);
    }
    const s = tracks.skew[k];
    if (Number.isFinite(s)) {
      const d = (s / maxSkew) * 28;
      if (Math.abs(d) > 0.2) out.push(`<path d="${sector(d > 0 ? 318 : 318 + d, d > 0 ? 318 + d : 318, a0, a1)}" fill="${d > 0 ? pos : neg}" fill-opacity="0.75" />`);
    }
  }
  // The cumulative skew.
  const cum = tracks.cumulative;
  const cmin = Math.min(...cum), cmax = Math.max(...cum), span = Math.max(1, cmax - cmin);
  const pts = cum.map((v, k) => pt(160 + (110 * (v - cmin)) / span, binAngle(k)));
  out.push(`<polygon points="${pts.join(' ')}" fill="${token('--accent')}" fill-opacity="0.06" stroke="${token('--accent')}" stroke-width="2" stroke-linejoin="round" />`);
  out.push(`<circle r="160" fill="none" stroke="${token('--hair')}" stroke-dasharray="2 4" />`);
  for (const [p, label, up] of [[tracks.origin, markers[0], false], [tracks.terminus, markers[1], true]]) {
    const a = angleOf(p, L);
    const r = up ? 282 : 148;
    out.push(`<path d="M${pt(r, a)}L${pt(r + (up ? -12 : 12), a - 0.03)}L${pt(r + (up ? -12 : 12), a + 0.03)}Z" fill="${ink}" />`);
    const lx = (up ? 296 : 128) * Math.cos(a), ly = (up ? 296 : 128) * Math.sin(a);
    out.push(`<text class="marker-label" x="${lx.toFixed(1)}" y="${(ly + 4).toFixed(1)}" text-anchor="middle">${label}</text>`);
  }
  // Repeats: chords joining the copies of each, inside the cumulative skew.
  for (const [p, q, colour, title] of links) {
    const a = angleOf(p, L), b = angleOf(q, L);
    out.push(`<path d="M${pt(150, a)}Q0,0 ${pt(150, b)}" fill="none" stroke="${colour}" stroke-width="1.2" stroke-opacity="0.55"><title>${esc(title)}</title></path>`);
  }
  // The window shown in the linear browser.
  if (window && window.end - window.begin < L) {
    const a0 = angleOf(window.begin, L), a1 = Math.max(angleOf(window.end, L), a0 + 0.006);
    out.push(`<path class="window-wedge" d="${sector(150, 478, a0, a1)}" />`);
  }
  // The middle.
  out.push(`<text class="center-title" text-anchor="middle" y="-6">${esc(title)}</text>`);
  out.push(`<text class="center-sub" text-anchor="middle" y="18">${esc(subtitle)}</text>`);
  return `<svg class="circle" viewBox="-615 -548 1230 1096" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="${esc(title)}: circular map">${out.join('')}</svg>`;
}

export const subtitleOf = (summary) => `${lengthText(summary.length)} · GC ${percent(summary.contigs[0]?.gc ?? NaN)}`;
