// Tracks along a whole genome, its contigs side by side: points (read depth in bins, B-allele
// frequencies), step segments (copy number), spans (calls, regions of LOH, the truth), and arcs
// joining two places (a deletion's ends, an inversion's junctions, a translocation between
// chromosomes). The page's hand: hairlines, few ticks, colour only where it means something.
//
//   genomeTracksSvg({ contigs, tracks, arcs, width, window })
//
// contigs: [{ name, length }]; a place is { contig, pos } in bases. A track is
//   { kind: 'points',   label, height, y: [lo, hi], points: [{ contig, pos, value, colour }], lines: [v] }
//   { kind: 'segments', label, height, y: [lo, hi], segments: [{ contig, start, end, value, colour, title }], lines: [v] }
//   { kind: 'spans',    label, height, spans: [{ contig, start, end, colour, title, row }] }
// window, when given, is { contig, start, end }: the tracks are drawn over that stretch alone.

import { esc } from '../dom.js';
import { token } from '../palette.js';
import { tickStep, number } from '../format.js';

export function genomeTracksSvg({ contigs, tracks = [], arcs = [], width = 980, window = null, arcHeight = 70 }) {
  const m = { l: 128, r: 12, t: 8 };
  const W = width - m.l - m.r;
  const shown = window ? contigs.filter((c) => c.name === window.contig) : contigs;
  const gapPx = shown.length > 1 ? 6 : 0;
  const span = (c) => (window ? window.end - window.start : c.length);
  const total = shown.reduce((a, c) => a + span(c), 0);
  const px = (W - gapPx * (shown.length - 1)) / Math.max(1, total);
  const offset = new Map();
  let acc = 0;
  for (const c of shown) { offset.set(c.name, acc); acc += span(c) * px + gapPx; }
  const X = (contig, pos) => {
    const o = offset.get(contig);
    if (o === undefined) return null;
    const p = window ? pos - window.start : pos;
    if (window && (p < 0 || p > span(shown[0]))) return null;
    return m.l + o + p * px;
  };
  const faint = token('--ink-faint'), hair = token('--hair'), hair2 = token('--hair-2');
  const out = [];
  let y = m.t;
  // Arcs above the tracks.
  if (arcs.length) {
    const base = y + arcHeight;
    for (const a of arcs) {
      const x0 = X(a.a.contig, a.a.pos), x1 = X(a.b.contig, a.b.pos);
      if (x0 === null || x1 === null) continue;
      const h = Math.min(arcHeight - 4, 12 + Math.abs(x1 - x0) * 0.35);
      out.push(`<path d="M${x0.toFixed(1)},${base} C${x0.toFixed(1)},${(base - h * 1.3).toFixed(1)} ${x1.toFixed(1)},${(base - h * 1.3).toFixed(1)} ${x1.toFixed(1)},${base}" fill="none" stroke="${a.colour || token('--accent')}" stroke-width="${a.width || 1.4}" stroke-opacity="${a.opacity ?? 0.85}" ${a.dash ? `stroke-dasharray="${a.dash}"` : ''}>${a.title ? `<title>${esc(a.title)}</title>` : ''}</path>`);
    }
    y = base + 6;
  }
  for (const t of tracks) {
    const h = t.height || 90;
    out.push(`<text class="tick-label" x="${m.l - 30}" y="${(y + h / 2 + 4).toFixed(1)}" text-anchor="end" style="font-style:italic">${esc(t.label || '')}</text>`);
    // A band behind each contig, alternating, so that their borders read.
    shown.forEach((c, k) => {
      const x0 = X(c.name, window ? window.start : 0);
      out.push(`<rect x="${x0.toFixed(1)}" y="${y}" width="${(span(c) * px).toFixed(1)}" height="${h}" fill="${k % 2 ? token('--band') : 'transparent'}" />`);
    });
    if (t.kind === 'points' || t.kind === 'segments') {
      const [lo, hi] = t.y;
      const sy = (v) => y + h - ((Math.min(hi, Math.max(lo, v)) - lo) / (hi - lo || 1)) * h;
      for (const v of t.lines || []) out.push(`<line x1="${m.l}" y1="${sy(v).toFixed(1)}" x2="${m.l + W}" y2="${sy(v).toFixed(1)}" stroke="${hair}" stroke-dasharray="${v === (t.lines[0]) ? '' : '2 3'}" />`);
      for (const v of t.ticks || [lo, hi]) out.push(`<text class="tick-label" x="${m.l - 3}" y="${(sy(v) + 3.5).toFixed(1)}" text-anchor="end" style="font-size:10px">${esc(number(v, 2))}</text>`);
      if (t.kind === 'points') {
        for (const p of t.points) {
          const x = X(p.contig, p.pos);
          if (x === null || !Number.isFinite(p.value)) continue;
          out.push(`<rect x="${(x - 0.9).toFixed(1)}" y="${(sy(p.value) - 0.9).toFixed(1)}" width="1.8" height="1.8" fill="${p.colour || faint}" fill-opacity="${t.opacity ?? 0.6}" />`);
        }
      }
      // Segments: a track's own, or drawn over its points (a segmentation over the data it segments).
      {
        for (const s of t.segments || []) {
          const x0 = X(s.contig, s.start), x1 = X(s.contig, s.end);
          if (x0 === null || x1 === null) continue;
          out.push(`<line x1="${x0.toFixed(1)}" y1="${sy(s.value).toFixed(1)}" x2="${Math.max(x0 + 1, x1).toFixed(1)}" y2="${sy(s.value).toFixed(1)}" stroke="${s.colour || token('--accent')}" stroke-width="${s.width || 2.4}" stroke-linecap="round">${s.title ? `<title>${esc(s.title)}</title>` : ''}</line>`);
        }
      }
    } else if (t.kind === 'spans') {
      const rows = Math.max(1, ...t.spans.map((s) => (s.row ?? 0) + 1));
      const rh = Math.min(16, (h - 4) / rows);
      for (const s of t.spans) {
        const x0 = X(s.contig, s.start), x1 = X(s.contig, s.end);
        if (x0 === null || x1 === null) continue;
        const yy = y + 2 + (s.row ?? 0) * rh;
        out.push(`<rect class="${s.cls || ''}" data-key="${esc(s.key ?? '')}" x="${x0.toFixed(1)}" y="${yy.toFixed(1)}" width="${Math.max(2.2, x1 - x0).toFixed(1)}" height="${(rh - 3).toFixed(1)}" rx="2" fill="${s.colour || token('--accent')}" fill-opacity="${s.opacity ?? 0.8}" style="${s.key != null ? 'cursor:pointer' : ''}">${s.title ? `<title>${esc(s.title)}</title>` : ''}</rect>`);
      }
    }
    out.push(`<line x1="${m.l}" y1="${y + h}" x2="${m.l + W}" y2="${y + h}" stroke="${hair2}" />`);
    y += h + 10;
  }
  // The contigs' names, or a scale over the window.
  if (window) {
    const step = tickStep(window.end - window.start, 8);
    for (let p = Math.ceil(window.start / step) * step; p <= window.end; p += step) {
      const x = X(window.contig, p);
      out.push(`<line x1="${x.toFixed(1)}" y1="${y - 10}" x2="${x.toFixed(1)}" y2="${y - 6}" stroke="${faint}" />`);
      const q = p + (window.offset || 0);  // positions drawn may be a region's; labels its chromosome's
      out.push(`<text class="tick-label" x="${x.toFixed(1)}" y="${y + 4}" text-anchor="middle">${esc(q >= 1e6 ? `${number(q / 1e6, 4)} Mb` : q >= 1e3 ? `${number(q / 1e3, 4)} kb` : String(q))}</text>`);
    }
  } else {
    for (const c of shown) {
      const x = X(c.name, c.length / 2);
      if (c.length * px > 14) out.push(`<text class="tick-label" x="${x.toFixed(1)}" y="${y + 4}" text-anchor="middle">${esc(c.label || c.name)}</text>`);
    }
  }
  const height = y + 14;
  return `<svg width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg">${out.join('')}</svg>`;
}
