// A sequence logo (Schneider and Stephens 1990): at each position of a motif, the bases stacked
// by their probability, the stack as tall as the position's information, 2 − H bits for DNA with
// H = −Σ p log₂ p; the likeliest base on top. Used for promoters' boxes, the motifs a model has
// learnt, and position weight matrices.
//
//   logoSvg({ rows, width?, height?, title?, bits = true })
//
// rows: [[pA, pC, pG, pT], ...]; with bits = false every stack is full height (probabilities).

import { esc } from '../dom.js';
import { BASE_COLOURS, token } from '../palette.js';

const BASES = ['A', 'C', 'G', 'T'];

export function information(row) {
  const s = row.reduce((a, b) => a + Math.max(0, b), 0) || 1;
  let h = 0;
  for (const v of row) { const p = Math.max(0, v) / s; if (p > 0) h -= p * Math.log2(p); }
  return 2 - h;
}

export function logoSvg({ rows, width = null, height = 92, title = '', bits = true, start = 1 }) {
  const colours = BASE_COLOURS();
  const colW = 17, m = { l: 30, r: 8, t: title ? 20 : 8, b: 20 };
  const W = width || m.l + m.r + rows.length * colW;
  const H = height - m.t - m.b;
  const out = [];
  if (title) out.push(`<text x="${m.l}" y="13" style="font: italic 12.5px var(--serif); fill:${token('--ink-2')}">${esc(title)}</text>`);
  // The axis: 0, 1 and 2 bits.
  for (const b of bits ? [0, 1, 2] : [0, 0.5, 1]) {
    const y = m.t + H - (b / (bits ? 2 : 1)) * H;
    out.push(`<line x1="${m.l - 4}" y1="${y.toFixed(1)}" x2="${m.l}" y2="${y.toFixed(1)}" stroke="${token('--ink-faint')}" />`);
    out.push(`<text class="tick-label" x="${m.l - 7}" y="${(y + 3.5).toFixed(1)}" text-anchor="end">${b}</text>`);
  }
  out.push(`<line x1="${m.l}" y1="${m.t}" x2="${m.l}" y2="${m.t + H}" stroke="${token('--hair-2')}" />`);
  out.push(`<text class="tick-label" transform="translate(9 ${m.t + H / 2}) rotate(-90)" text-anchor="middle" style="font-style:italic">${bits ? 'bits' : 'p'}</text>`);
  rows.forEach((row, i) => {
    const s = row.reduce((a, b) => a + Math.max(0, b), 0) || 1;
    const p = row.map((v) => Math.max(0, v) / s);
    const total = bits ? (information(row) / 2) * H : H;
    const order = BASES.map((b, k) => [b, p[k]]).sort((a, b) => a[1] - b[1]);
    let y = m.t + H;
    const x = m.l + i * colW;
    for (const [b, q] of order) {
      const h = q * total;
      if (h < 0.4) continue;
      // A glyph of cap height 0.72 em scaled to h: font-size 100 and a vertical scale.
      const sy = h / 72, sx = (colW - 1.5) / 62;
      out.push(`<text transform="translate(${(x + colW / 2).toFixed(2)} ${y.toFixed(2)}) scale(${sx.toFixed(4)} ${sy.toFixed(4)})" text-anchor="middle" style="font: 700 100px 'Source Serif 4', Georgia, serif; fill:${colours[b]}">${b}</text>`);
      y -= h;
    }
    if (rows.length <= 40 || i % 5 === 0) out.push(`<text class="tick-label" x="${(x + colW / 2).toFixed(1)}" y="${m.t + H + 13}" text-anchor="middle">${start + i}</text>`);
  });
  return `<svg width="${W}" height="${height}" viewBox="0 0 ${W} ${height}" xmlns="http://www.w3.org/2000/svg">${out.join('')}</svg>`;
}
