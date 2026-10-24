// The linear browser of a window [begin, end) of a contig, as an SVG string of stacked tracks:
// the scale, the features packed into rows (arrows by strand), the ORFs of the six frames, GC
// content about the genome's mean, GC skew, CpG observed/expected, the entropy of 3-mers, the
// letters themselves when there is room, and the hits of a search.

import { esc } from '../dom.js';
import { featureColour, token, BASE_COLOURS, CATEGORIES } from '../palette.js';
import { tickStep, tickLabel } from '../format.js';

export const LEFT = 74, RIGHT = 14;

// Features in rows so that none overlaps another in pixels (greedy, by start).
function pack(features, x, gap = 4) {
  const rows = [];
  const placed = [];
  const sorted = [...features].sort((a, b) => a.lo - b.lo || b.hi - a.hi);
  for (const f of sorted) {
    const x0 = x(f.lo), x1 = x(f.hi);
    let r = rows.findIndex((end) => end + gap <= x0);
    if (r < 0) { r = rows.length; rows.push(-Infinity); }
    rows[r] = Math.max(x1, x0 + 2);
    placed.push({ ...f, row: r });
  }
  return { rows: rows.length, placed };
}

export function linearSvg({ width, tracks, selected, hits, contig, show, extra = {} }) {
  const { begin, end } = tracks;
  const span = Math.max(1, end - begin);
  const w = Math.max(200, width);
  const inner = w - LEFT - RIGHT;
  const x = (p) => LEFT + ((p - begin) / span) * inner;
  const out = [];
  let y = 8;
  const ink = token('--ink'), faint = token('--ink-faint'), hair = token('--hair'), hair2 = token('--hair-2');
  const name = (text, yy, sub = '') => out.push(`<text class="track-name" x="8" y="${yy}">${esc(text)}${sub ? `<tspan class="sub" x="8" dy="11">${esc(sub)}</tspan>` : ''}</text>`);

  // The scale.
  const step = tickStep(span, Math.max(3, Math.floor(inner / 110)));
  out.push(`<line x1="${LEFT}" y1="${y + 18}" x2="${w - RIGHT}" y2="${y + 18}" stroke="${hair2}" />`);
  for (let p = Math.ceil(begin / step) * step; p <= end; p += step) {
    const xx = x(p);
    out.push(`<line x1="${xx.toFixed(1)}" y1="${y + 13}" x2="${xx.toFixed(1)}" y2="${y + 18}" stroke="${faint}" />`);
    out.push(`<text class="tick-label" x="${xx.toFixed(1)}" y="${y + 9}" text-anchor="middle">${esc(tickLabel(p, step))}</text>`);
  }
  y += 30;

  // Features.
  const items = [];
  for (const f of tracks.features) {
    const lo = Math.min(...f.parts.map((p) => p[0])), hi = Math.max(...f.parts.map((p) => p[1]));
    // A join across the origin is two pieces that are far apart: drawn as its parts' extents each.
    const wraps = f.parts.length > 1 && hi - lo > tracks.length / 2;
    if (wraps) for (const p of f.parts) items.push({ f, lo: p[0], hi: p[1], parts: [p] });
    else items.push({ f, lo, hi, parts: f.parts });
  }
  const visible = items.filter((it) => it.hi > begin && it.lo < end);
  const { rows, placed } = pack(visible, x);
  const rowH = 17, maxRows = 14;
  name('features', y + 11);
  // Too many to read one by one: the two strands as rows of marks, the forward above.
  if (rows > maxRows || span / inner > 400) {
    const h = 12;
    out.push(`<text class="tick-label" x="${LEFT - 6}" y="${y + 9}" text-anchor="end" style="font-size:9px">+</text>`);
    out.push(`<text class="tick-label" x="${LEFT - 6}" y="${y + h + 13}" text-anchor="end" style="font-size:9px">−</text>`);
    for (const it of visible) {
      const a = Math.max(LEFT, x(it.lo)), z = Math.min(w - RIGHT, x(it.hi));
      const top = it.f.strand === '-' ? y + h + 4 : y;
      out.push(`<rect class="feature${selected === it.f.i ? ' selected' : ''}" data-i="${it.f.i}" x="${a.toFixed(1)}" y="${top}" width="${Math.max(0.8, z - a).toFixed(2)}" height="${h}" fill="${featureColour(it.f.key)}" fill-opacity="0.75" />`);
    }
    out.push(`<text class="tick-label" x="${w - RIGHT}" y="${y + 2 * h + 18}" text-anchor="end">${visible.length} features: zoom in to read them</text>`);
    y += 2 * h + 26;
  }
  const shownRows = rows > maxRows || span / inner > 400 ? 0 : rows;
  for (const it of shownRows ? placed : []) {
    if (it.row >= maxRows) continue;
    const f = it.f;
    const top = y + it.row * rowH, h = 11;
    const col = featureColour(f.key);
    const x0 = Math.max(LEFT, x(it.lo)), x1 = Math.min(w - RIGHT, x(it.hi));
    if (it.parts.length > 1) out.push(`<line x1="${x0.toFixed(1)}" y1="${top + h / 2}" x2="${x1.toFixed(1)}" y2="${top + h / 2}" stroke="${col}" stroke-width="1" />`);
    for (const [b, e] of it.parts) {
      const a = Math.max(LEFT, x(b)), z = Math.min(w - RIGHT, x(e));
      if (z < LEFT || a > w - RIGHT) continue;
      const len = Math.max(1.2, z - a);
      // An arrow: a pointed end on the side the feature is read towards, when there is room.
      const tip = Math.min(6, len / 2);
      const d = f.strand === '-'
        ? `M${a.toFixed(1)},${top + h / 2}L${(a + tip).toFixed(1)},${top}H${(a + len).toFixed(1)}V${top + h}H${(a + tip).toFixed(1)}Z`
        : `M${a.toFixed(1)},${top}H${(a + len - tip).toFixed(1)}L${(a + len).toFixed(1)},${top + h / 2}L${(a + len - tip).toFixed(1)},${top + h}H${a.toFixed(1)}Z`;
      out.push(`<path class="feature${selected === f.i ? ' selected' : ''}" data-i="${f.i}" d="${d}" fill="${col}" fill-opacity="${f.pseudo ? 0.35 : 0.88}" />`);
    }
    const label = f.gene || f.name;
    if (label && x1 - x0 > label.length * 6.4 + 10) out.push(`<text class="feature-label" x="${((x0 + x1) / 2).toFixed(1)}" y="${top + h + 5.5}" text-anchor="middle" style="font-size:10px">${esc(label)}</text>`);
  }
  if (shownRows) y += shownRows * rowH + 8;
  else if (!visible.length) y += rowH;

  // The six frames.
  if (show.orfs && tracks.orfs) {
    const fh = 6, fg = 3;
    name('ORFs', y + 8);
    const labels = ['+1', '+2', '+3', '−1', '−2', '−3'];
    for (let k = 0; k < 6; k++) out.push(`<text class="tick-label" x="${LEFT - 6}" y="${y + k * (fh + fg) + 6}" text-anchor="end" style="font-size:8.5px">${labels[k]}</text>`);
    for (const [b, e, frame] of tracks.orfs) {
      const a = Math.max(LEFT, x(b)), z = Math.min(w - RIGHT, x(e));
      if (z <= a) continue;
      out.push(`<rect class="orf" x="${a.toFixed(1)}" y="${y + frame * (fh + fg)}" width="${Math.max(1, z - a).toFixed(1)}" height="${fh}" />`);
    }
    y += 6 * (fh + fg) + 10;
  }

  // What the sequence alone says: predicted genes, CpG islands, repeats.
  const inWindow = (b, e) => e > begin && b < end;
  if (extra.genes) {
    const h = 9;
    name('genes', y + 9, 'predicted');
    out.push(`<text class="tick-label" x="${LEFT - 6}" y="${y + 8}" text-anchor="end" style="font-size:9px">+</text><text class="tick-label" x="${LEFT - 6}" y="${y + h + 10}" text-anchor="end" style="font-size:9px">−</text>`);
    for (const [b, e, rev, score] of extra.genes) {
      if (!inWindow(b, e)) continue;
      const a = Math.max(LEFT, x(b)), z = Math.min(w - RIGHT, x(e));
      out.push(`<rect x="${a.toFixed(1)}" y="${rev ? y + h + 3 : y}" width="${Math.max(0.8, z - a).toFixed(1)}" height="${h - 1}" rx="1.5" fill="${token('--cds')}" fill-opacity="0.32" stroke="${token('--cds')}" stroke-width="0.6"><title>a predicted gene, ${e - b} bp, ${rev ? 'reverse' : 'forward'} strand, ${score} bits</title></rect>`);
    }
    y += 2 * h + 12;
  }
  if (extra.islands) {
    name('CpG', y + 9, 'islands');
    for (const [b, e, lo, gc, oe] of extra.islands) {
      if (!inWindow(b, e)) continue;
      const a = Math.max(LEFT, x(b)), z = Math.min(w - RIGHT, x(e));
      out.push(`<rect x="${a.toFixed(1)}" y="${y + 1}" width="${Math.max(1, z - a).toFixed(1)}" height="9" fill="${token('--G')}" fill-opacity="0.85"><title>a CpG island, ${e - b} bp, GC ${(gc * 100).toFixed(1)}%, CpG o/e ${oe}, ${lo} bits</title></rect>`);
    }
    y += 27;
  }
  if (extra.repeats) {
    name('repeats', y + 9);
    extra.repeats.forEach((r, k) => {
      for (const [p, rev] of r.at) {
        if (!inWindow(p, p + r.length)) continue;
        const a = Math.max(LEFT, x(p)), z = Math.min(w - RIGHT, x(p + r.length));
        out.push(`<rect x="${a.toFixed(1)}" y="${y + 1 + (rev ? 5 : 0)}" width="${Math.max(1, z - a).toFixed(1)}" height="5" fill="${CATEGORIES[k % CATEGORIES.length]}" fill-opacity="0.9"><title>a repeat of ${r.length} bp, ${r.at.length} copies${rev ? ' (this copy reversed)' : ''}</title></rect>`);
      }
    });
    y += 22;
  }

  // Composition: GC about the mean, skew, CpG o/e, entropy.
  const nb = tracks.gc.length;
  const bx = (k) => LEFT + (k / nb) * inner;
  const line = (values, top, h, lo, hi, colour, fill = false) => {
    let d = '';
    values.forEach((v, k) => {
      if (!Number.isFinite(v)) return;
      const yy = top + h - ((Math.min(hi, Math.max(lo, v)) - lo) / (hi - lo)) * h;
      d += `${d ? 'L' : 'M'}${(bx(k) + inner / nb / 2).toFixed(1)},${yy.toFixed(1)}`;
    });
    if (d) out.push(`<path d="${d}" fill="none" stroke="${colour}" stroke-width="1.3" stroke-linejoin="round" />`);
    if (fill && d) out.push(`<path d="${d}V${top + h}H${(bx(0) + inner / nb / 2).toFixed(1)}Z" fill="${colour}" fill-opacity="0.12" />`);
  };
  if (show.gc) {
    const h = 48;
    name('GC', y + 11);
    const gcs = tracks.gc.filter(Number.isFinite);
    const lo = Math.min(tracks.mean_gc - 0.05, ...gcs), hi = Math.max(tracks.mean_gc + 0.05, ...gcs);
    const my = y + h - ((tracks.mean_gc - lo) / (hi - lo)) * h;
    out.push(`<line x1="${LEFT}" y1="${my.toFixed(1)}" x2="${w - RIGHT}" y2="${my.toFixed(1)}" stroke="${hair2}" stroke-dasharray="3 3" />`);
    out.push(`<text class="tick-label" x="${LEFT - 6}" y="${(my + 3).toFixed(1)}" text-anchor="end">${(100 * tracks.mean_gc).toFixed(1)}%</text>`);
    // Bars above and below the mean.
    const hiC = token('--high'), loC = token('--low');
    tracks.gc.forEach((g, k) => {
      if (!Number.isFinite(g)) return;
      const yy = y + h - ((g - lo) / (hi - lo)) * h;
      out.push(`<rect x="${bx(k).toFixed(1)}" y="${Math.min(yy, my).toFixed(1)}" width="${Math.max(0.6, inner / nb).toFixed(2)}" height="${Math.max(0.4, Math.abs(yy - my)).toFixed(1)}" fill="${g >= tracks.mean_gc ? hiC : loC}" fill-opacity="0.7" />`);
    });
    y += h + 12;
  }
  if (show.skew) {
    const h = 36;
    name('GC skew', y + 11);
    const m = Math.max(1e-6, ...tracks.skew.filter(Number.isFinite).map(Math.abs));
    const mid = y + h / 2;
    out.push(`<line x1="${LEFT}" y1="${mid}" x2="${w - RIGHT}" y2="${mid}" stroke="${hair2}" />`);
    const pos = token('--pos'), neg = token('--neg');
    tracks.skew.forEach((s, k) => {
      if (!Number.isFinite(s) || s === 0) return;
      const hh = (Math.abs(s) / m) * (h / 2);
      out.push(`<rect x="${bx(k).toFixed(1)}" y="${(s > 0 ? mid - hh : mid).toFixed(1)}" width="${Math.max(0.6, inner / nb).toFixed(2)}" height="${hh.toFixed(1)}" fill="${s > 0 ? pos : neg}" fill-opacity="0.7" />`);
    });
    y += h + 12;
  }
  if (show.cpg) {
    const h = 32;
    name('CpG o/e', y + 11);
    const ty = y + h - (0.6 / 1.6) * h;
    out.push(`<line x1="${LEFT}" y1="${ty.toFixed(1)}" x2="${w - RIGHT}" y2="${ty.toFixed(1)}" stroke="${hair2}" stroke-dasharray="2 3" />`);
    out.push(`<text class="tick-label" x="${LEFT - 6}" y="${(ty + 3).toFixed(1)}" text-anchor="end">0.6</text>`);
    line(tracks.cpg, y, h, 0, 1.6, token('--ncrna'), true);
    y += h + 12;
  }
  if (show.entropy) {
    const h = 28;
    name('entropy', y + 11);
    const es = tracks.entropy.filter(Number.isFinite);
    line(tracks.entropy, y, h, Math.min(4, ...es), 6, token('--C'));
    y += h + 12;
  }

  // The letters.
  if (tracks.letters) {
    const h = 20;
    name('sequence', y + 13);
    const per = inner / span;
    const colours = BASE_COLOURS();
    for (let i = 0; i < tracks.letters.length; i++) {
      const b = tracks.letters[i];
      const xx = x(begin + i);
      if (per >= 8) out.push(`<text class="base-letter" x="${(xx + per / 2).toFixed(1)}" y="${y + 14}" fill="${colours[b] || colours.N}">${b}</text>`);
      else out.push(`<rect x="${xx.toFixed(2)}" y="${y + 3}" width="${Math.max(0.5, per).toFixed(2)}" height="${h - 6}" fill="${colours[b] || colours.N}" />`);
    }
    y += h + 8;
  }

  // Hits of a search.
  if (hits && hits.length) {
    name('search', y + 9);
    for (const [c, p, rev] of hits) {
      if (c !== contig || p < begin - 1000 || p > end) continue;
      const xx = x(p);
      out.push(`<path class="hit" d="M${xx.toFixed(1)},${y + 2}l4,7h-8z" fill="${rev ? token('--neg') : token('--pos')}" />`);
    }
    y += 16;
  }
  y += 6;
  return { svg: `<svg width="${w}" height="${y}" viewBox="0 0 ${w} ${y}" xmlns="http://www.w3.org/2000/svg">${out.join('')}</svg>`, height: y, x, position: (px) => begin + ((px - LEFT) / inner) * span };
}
