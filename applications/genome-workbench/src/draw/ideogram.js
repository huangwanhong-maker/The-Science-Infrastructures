// Chromosomes drawn as cytogeneticists draw them: rods banded by Giemsa's stain, narrowed at the
// centromere, capped at the telomeres; and painted, as a FISH paint colours each chromosome, so
// that a rearranged chromosome shows where its pieces came from.
//
//   paint(chrom)                    the paint of a chromosome (1–22, X, Y), as CSS
//   rodSvg(segments, ideogram, o)   one chromosome or derivative, its segments of source
//                                   chromosomes in order from pter, each with its bands
//   karyogramSvg({ideogram, copies, ...})
//                                   the karyotype: every chromosome's copies, in ISCN's order
//   imbalanceSvg({ideogram, segments})
//                                   the reference ideogram with losses to the left and gains to
//                                   the right, as copy number departs from two
//
// The ideogram is { chromosomes: [{ name, length, bands: [{ name, start, end, stain }] }] } in
// bases (0-based, half-open); stains are gneg, gpos25…gpos100, acen, gvar and stalk. A segment
// is { chrom, start, end, reverse }: the bases [start, end) of chrom, read from end to start
// when reversed.

import { esc } from '../dom.js';
import { token } from '../palette.js';

// ---- paints --------------------------------------------------------------------------------------
// OKLCH to sRGB (Ottosson 2020), for paints of equal chroma and spaced hue.
function oklch(L, C, h) {
  const a = C * Math.cos((h * Math.PI) / 180), b = C * Math.sin((h * Math.PI) / 180);
  const l_ = L + 0.3963377774 * a + 0.2158037573 * b, m_ = L - 0.1055613458 * a - 0.0638541728 * b, s_ = L - 0.0894841775 * a - 1.291485548 * b;
  const l = l_ ** 3, m = m_ ** 3, s = s_ ** 3;
  const lin = [4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s, -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s, -0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s];
  const g = (x) => { x = Math.min(1, Math.max(0, x)); return Math.round(255 * (x <= 0.0031308 ? 12.92 * x : 1.055 * x ** (1 / 2.4) - 0.055)); };
  return `rgb(${lin.map(g).join(',')})`;
}
export const CHROMOSOMES = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y'];
// Hues by the golden angle, so that neighbours differ most; lightness in three steps.
const PAINTS = CHROMOSOMES.map((_, i) => oklch([0.70, 0.62, 0.78][i % 3], 0.105, (18 + i * 137.508) % 360));
export function paint(chrom) {
  const k = CHROMOSOMES.indexOf(String(chrom).replace(/^chr/, ''));
  return k >= 0 ? PAINTS[k] : token('--other');
}

// The grey of a Giemsa stain, on the page's ink.
function stainFill(stain, night) {
  const ink = night ? [236, 222, 232] : [58, 36, 48];
  const paper = night ? [62, 46, 60] : [251, 246, 247];
  const mixed = (t) => `rgb(${ink.map((v, i) => Math.round(paper[i] + (v - paper[i]) * t)).join(',')})`;
  const m = /^gpos(\d+)$/.exec(stain || '');
  if (m) return mixed(0.18 + 0.72 * (Number(m[1]) / 100));
  if (stain === 'gpos') return mixed(0.6);
  if (stain === 'gvar') return mixed(0.32);
  if (stain === 'stalk') return mixed(0.12);
  if (stain === 'acen') return mixed(0.42);
  return mixed(0.03);  // gneg
}

const isNight = () => document.documentElement.dataset.theme === 'night';
const chromOf = (ideogram, name) => ideogram.chromosomes.find((c) => c.name === String(name).replace(/^chr/, ''));

// The bands of one segment, in drawing order (from the segment's first base read), clipped.
function segmentBands(ideogram, s) {
  const c = chromOf(ideogram, s.chrom);
  if (!c) return [];
  const out = [];
  for (const b of c.bands) {
    const lo = Math.max(b.start, s.start), hi = Math.min(b.end, s.end);
    if (hi > lo) out.push({ ...b, lo, hi });
  }
  if (s.reverse) out.reverse();
  return out;
}

// ---- one rod --------------------------------------------------------------------------------------
// A rod's geometry: x, y of its top, width w, scale (pixels a base), and whether to paint.
// Returns { svg, height, junctions: [y], centromeres: [y] }.
export function rodSvg(segments, ideogram, { x = 0, y = 0, w = 12, scale = 1e-6, painted = true, strength = 1, label = '', id = '' } = {}) {
  const night = isNight();
  const parts = [], paints = [], marks = [];
  let at = y;
  const total = segments.reduce((a, s) => a + Math.max(0, s.end - s.start), 0);
  const height = total * scale;
  // The narrowing: centromeric bands (acen) pinch the rod to a waist, stalks to a thread.
  const centromeres = [];
  segments.forEach((s, k) => {
    const bands = segmentBands(ideogram, s);
    const top = at;
    for (const b of bands) {
      const h = (b.hi - b.lo) * scale;
      let bw = w, bx = x;
      if (b.stain === 'acen') { bw = w * 0.55; bx = x + (w - bw) / 2; centromeres.push(at + h / 2); }
      else if (b.stain === 'stalk') { bw = w * 0.3; bx = x + (w - bw) / 2; }
      parts.push(`<rect x="${bx.toFixed(2)}" y="${at.toFixed(2)}" width="${bw.toFixed(2)}" height="${Math.max(0.3, h).toFixed(2)}" fill="${stainFill(b.stain, night)}"><title>${esc(`${String(s.chrom).replace(/^chr/, '')}${b.name}`)}</title></rect>`);
      if (b.stain === 'gvar') parts.push(`<rect x="${bx.toFixed(2)}" y="${at.toFixed(2)}" width="${bw.toFixed(2)}" height="${Math.max(0.3, h).toFixed(2)}" fill="url(#gvar-hatch)" />`);
      at += h;
    }
    if (painted) paints.push(`<rect x="${(x - 0.5).toFixed(2)}" y="${top.toFixed(2)}" width="${(w + 1).toFixed(2)}" height="${(at - top).toFixed(2)}" fill="${paint(s.chrom)}" />`);
    if (k > 0) marks.push(top);
  });
  const r = Math.min(w / 2, height / 2);
  const clip = `rod-${id || Math.random().toString(36).slice(2, 9)}`;
  const outline = `M${x + r},${y} h${w - 2 * r} a${r},${r} 0 0 1 ${r},${r} v${Math.max(0, height - 2 * r)} a${r},${r} 0 0 1 ${-r},${r} h${-(w - 2 * r)} a${r},${r} 0 0 1 ${-r},${-r} v${-Math.max(0, height - 2 * r)} a${r},${r} 0 0 1 ${r},${-r} z`;
  const out = [];
  out.push(`<clipPath id="${clip}"><path d="${outline}" /></clipPath>`);
  out.push(`<g clip-path="url(#${clip})">${parts.join('')}`);
  if (painted) out.push(`<g style="mix-blend-mode:${night ? 'screen' : 'multiply'}" opacity="${(night ? 0.5 : 0.62) * strength}">${paints.join('')}</g>`);
  out.push('</g>');
  // The waist: two notches at each centromere.
  for (const cy of centromeres) {
    out.push(`<path d="M${x - 1.5},${cy - 2.6} l${w * 0.24 + 1.5},2.6 l${-(w * 0.24 + 1.5)},2.6 z M${x + w + 1.5},${cy - 2.6} l${-(w * 0.24 + 1.5)},2.6 l${w * 0.24 + 1.5},2.6 z" fill="${token('--surface')}" />`);
  }
  out.push(`<path d="${outline}" fill="none" stroke="${token('--hair-3')}" stroke-width="0.8" />`);
  // Junctions between segments: a rose tick on each side, the "::" of ISCN's long form.
  for (const jy of marks) {
    out.push(`<path d="M${x - 4},${jy} h3 M${x + w + 1},${jy} h3" stroke="${token('--T')}" stroke-width="1.6" stroke-linecap="round" />`);
  }
  if (label) out.push(`<text class="tick-label" x="${x + w / 2}" y="${y + height + 13}" text-anchor="middle">${esc(label)}</text>`);
  return { svg: out.join(''), height, junctions: marks, centromeres };
}

const HATCH = () => `<defs><pattern id="gvar-hatch" width="3" height="3" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="3" stroke="${token('--ink-faint')}" stroke-width="0.8" /></pattern></defs>`;

// The whole chromosome as one segment.
export const wholeChromosome = (ideogram, name) => {
  const c = chromOf(ideogram, name);
  return c ? [{ chrom: c.name, start: 0, end: c.length, reverse: false }] : [];
};

// ---- the karyogram -------------------------------------------------------------------------------
// copies: { '1': [ { label, segments }, ... ], ... } — each chromosome's copies, the derivatives
// under the chromosome whose centromere they carry. Rows as the karyogram is laid out: 1–5, 6–12,
// 13–18, 19–22 with X and Y.
export function karyogramSvg({ ideogram, copies, width = 980, painted = true, highlight = null }) {
  const rows = [['1', '2', '3', '4', '5'], ['6', '7', '8', '9', '10', '11', '12'], ['13', '14', '15', '16', '17', '18'], ['19', '20', '21', '22', 'X', 'Y']];
  const longest = Math.max(...ideogram.chromosomes.map((c) => c.length));
  const rowHeight = (row) => Math.max(...row.map((n) => {
    const list = copies[n] || [];
    return Math.max(0, ...list.map((c) => c.segments.reduce((a, s) => a + s.end - s.start, 0)));
  })) ;
  const maxRod = 200, scale = maxRod / longest;
  const w = 14, gap = 6, groupGap = 30;
  const out = [HATCH()];
  let y = 18;
  rows.forEach((row) => {
    const groups = row.map((n) => ({ n, list: copies[n] || [] }));
    const groupWidth = (g) => Math.max(1, g.list.length) * (w + gap) - gap;
    const total = groups.reduce((a, g) => a + groupWidth(g), 0) + groupGap * (groups.length - 1);
    let x = Math.max(10, (width - total) / 2);
    const tallest = rowHeight(row) * scale;
    for (const g of groups) {
      const gw = groupWidth(g);
      const abnormal = g.list.some((c) => c.abnormal);
      if (g.list.length === 0) {
        out.push(`<text class="tick-label" x="${x + gw / 2}" y="${y + tallest + 16}" text-anchor="middle" style="opacity:.55">${esc(g.n)}</text>`);
        out.push(`<text class="tick-label" x="${x + gw / 2}" y="${y + tallest / 2}" text-anchor="middle" style="opacity:.55">—</text>`);
      }
      g.list.forEach((c, k) => {
        const h = c.segments.reduce((a, s) => a + s.end - s.start, 0) * scale;
        const top = y + tallest - h;  // aligned at the bottom, as karyograms are
        // the normal chromosomes painted faintly, so that the rearranged ones stand out
        const rod = rodSvg(c.segments, ideogram, { x: x + k * (w + gap), y: top, w, scale, painted, strength: c.abnormal ? 1 : 0.38, id: `k${g.n}-${k}` });
        const lit = highlight != null && c.key === highlight;
        out.push(`<g class="rod${c.abnormal ? ' abnormal' : ''}${lit ? ' lit' : ''}" data-key="${esc(c.key || '')}" style="cursor:${c.key ? 'pointer' : 'default'}"><title>${esc(c.label || g.n)}</title>${rod.svg}</g>`);
        if (c.abnormal) out.push(`<circle cx="${x + k * (w + gap) + w / 2}" cy="${top - 6}" r="2" fill="${token('--T')}" />`);
      });
      out.push(`<text class="tick-label" x="${x + gw / 2}" y="${y + tallest + 16}" text-anchor="middle" style="${abnormal ? `fill:${token('--accent-deep')};font-weight:600` : ''}">${esc(g.n)}</text>`);
      x += gw + groupGap;
    }
    y += tallest + 40;
  });
  return `<svg width="${width}" height="${y}" viewBox="0 0 ${width} ${y}" xmlns="http://www.w3.org/2000/svg">${out.join('')}</svg>`;
}

// ---- gains and losses on the reference ---------------------------------------------------------------
// segments: [{ chrom, start, end, copies, ploidy? }] where copies departs from the chromosome's
// normal number (2, or 1 for X and Y in a male when ploidy is given). Losses are drawn to the
// left of each rod in blue, gains to the right in rose, a bar's length its departure.
export function imbalanceSvg({ ideogram, segments, width = 980, sex = 'XX', painted = false }) {
  const names = CHROMOSOMES.filter((n) => chromOf(ideogram, n) && (n !== 'Y' || sex.includes('Y')));
  const longest = Math.max(...ideogram.chromosomes.map((c) => c.length));
  const rodH = 170, scale = rodH / longest, w = 9;
  const slot = Math.min(42, (width - 20) / names.length);
  const out = [HATCH()];
  const top = 16;
  const lossC = token('--low'), gainC = token('--T');
  names.forEach((n, i) => {
    const x = 10 + i * slot + slot / 2 - w / 2;
    const c = chromOf(ideogram, n);
    const rod = rodSvg(wholeChromosome(ideogram, n), ideogram, { x, y: top, w, scale, painted, id: `i${n}` });
    out.push(rod.svg);
    const normal = (n === 'X' && sex === 'XY') || n === 'Y' ? 1 : 2;
    for (const s of segments.filter((q) => String(q.chrom).replace(/^chr/, '') === n)) {
      const d = s.copies - (s.normal ?? normal);
      if (!d && !s.loh) continue;
      const y0 = top + s.start * scale, h = Math.max(1.2, (s.end - s.start) * scale);
      if (d < 0) out.push(`<rect x="${(x - 3 - 3 * Math.min(3, -d)).toFixed(1)}" y="${y0.toFixed(1)}" width="${(3 * Math.min(3, -d)).toFixed(1)}" height="${h.toFixed(1)}" fill="${lossC}"><title>${esc(`${n}: ${s.copies} copies ${s.label || ''}`)}</title></rect>`);
      if (d > 0) out.push(`<rect x="${(x + w + 3).toFixed(1)}" y="${y0.toFixed(1)}" width="${(3 * Math.min(3, d)).toFixed(1)}" height="${h.toFixed(1)}" fill="${gainC}"><title>${esc(`${n}: ${s.copies} copies ${s.label || ''}`)}</title></rect>`);
      if (s.loh && !d) out.push(`<rect x="${(x - 5).toFixed(1)}" y="${y0.toFixed(1)}" width="2" height="${h.toFixed(1)}" fill="${token('--high')}"><title>${esc(`${n}: copy-neutral loss of heterozygosity`)}</title></rect>`);
    }
    out.push(`<text class="tick-label" x="${x + w / 2}" y="${top + c.length * scale + 14}" text-anchor="middle">${esc(n)}</text>`);
  });
  const height = top + rodH + 26;
  return `<svg width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg">${out.join('')}</svg>`;
}
