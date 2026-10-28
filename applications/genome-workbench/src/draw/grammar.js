// The figures and the notation of group grammars (the Grammar view and the toolbox's grammar
// tools; domains/genomics/theory/group-grammars.md): the Klein group's square on the bases, the chain
// of breakings, parse trees, hairpins as stem and loop and as arcs over their sequence, bars of the
// symmetries a genome keeps, the genetic code on (C₂)⁶ as a grid, and the Wedderburn blocks of a
// group algebra. Each returns an SVG string drawn in the tokens of styles.css (as CSS variables, so
// that a change of theme recolours it), and every piece of text in it passes through esc().
//
// And the notation: grammars written one rule per line (S -> a S b | ε), groups as generators (the
// images of the alphabet's letters, ρ when the generator also reverses words), cyclotomic numbers.

import { esc } from '../dom.js';

const SUB = '₀₁₂₃₄₅₆₇₈₉', SUP = '⁰¹²³⁴⁵⁶⁷⁸⁹';
export const sub = (n) => String(n).replace(/\d/g, (d) => SUB[d]);
const sup = (n) => String(n).replace(/\d/g, (d) => SUP[d]);
const BASE = { A: 'var(--A)', C: 'var(--C)', G: 'var(--G)', T: 'var(--T)', N: 'var(--N)' };
const isBase = (c) => 'ACGT'.includes(c);

// ---- notation ------------------------------------------------------------------------------------------
// morphe's cyclotomic numbers as written ("1 - 2*ζ12^3 + 1/2*ζ12"), set as mathematics: ζ₁₂³, ½, −,
// and ζ₄ as i.
export function cyclotomic(s) {
  return String(s)
    .replace(/\*/g, '')
    .replace(/ζ4(?!\d)(\^\d+)?/g, (m, p) => (p ? `i${sup(p.slice(1))}` : 'i'))
    .replace(/ζ(\d+)(\^(\d+))?/g, (m, n, p, e) => `ζ${sub(n)}${e ? sup(e) : ''}`)
    .replace(/\b([123])\/([234])(?!\d)/g, (m, a, b) => ({ '1/2': '½', '1/3': '⅓', '2/3': '⅔', '1/4': '¼', '3/4': '¾' }[`${a}/${b}`] || m))
    .replace(/-/g, '−');
}

// Rules written one to a line (or separated by ;), alternatives by |, the arrow -> or →, the empty
// word ε: as the engine's [left, right] pairs.
export function rulesFromText(text) {
  const out = [];
  for (const raw of String(text).split(/[\n;]/)) {
    const line = raw.replace(/#.*$/, '').trim();
    if (!line) continue;
    const m = line.match(/^(.*?)\s*(?:->|→|::=)\s*(.*)$/);
    if (!m) throw new Error(`a rule needs an arrow: ${line}`);
    if (!m[1].trim()) throw new Error(`a rule has an empty left side: ${line}`);
    for (const alt of m[2].split('|')) out.push([m[1].trim(), alt.trim() === '' ? 'ε' : alt.trim()]);
  }
  if (!out.length) throw new Error('write at least one rule, such as S -> a S b | ε');
  return out;
}

// Generators written as the images of the letters, one to a line (or separated by ;), followed by
// ρ (or rev) when the generator reverses words: [[letters, reverses], ...].
export function generatorsFromText(text) {
  const out = [];
  for (const raw of String(text).split(/[\n;,]/)) {
    const line = raw.trim();
    if (!line) continue;
    const reverses = /(\s|^)(ρ|rev|reverse|reversal)$/i.test(line);
    const letters = line.replace(/(\s|^)(ρ|rev|reverse|reversal)$/i, '').trim();
    out.push([letters, reverses]);
  }
  return out;
}

// ---- parse trees -----------------------------------------------------------------------------------------
// A parse tree (the engine's nested nodes and leaves), laid out with the word along its foot: each
// leaf (or ε) a slot, a node centred over its children, the leaves dropped to the last row.
export function treeSvg(tree, { slot = 24, level = 44, maxWidth = 1400 } = {}) {
  const nodes = [];
  let slots = 0, depth = 0;
  const walk = (t, d, parent) => {
    const me = { t, d, parent, x: 0, kids: [] };
    nodes.push(me);
    depth = Math.max(depth, d);
    if (t.leaf) { me.x = slots++; return me; }
    if (!t.children.length) { me.eps = true; me.x = slots++; return me; }
    for (const c of t.children) me.kids.push(walk(c, d + 1, me));
    me.x = (me.kids[0].x + me.kids[me.kids.length - 1].x) / 2;
    return me;
  };
  walk(tree, 0, null);
  const s = Math.max(12, Math.min(slot, (maxWidth - 40) / Math.max(1, slots)));
  const W = Math.max(160, slots * s + 40), H = (depth + 1) * level + 40;
  const X = (n) => 20 + s / 2 + n.x * s;
  const Y = (n) => (n.t.leaf || n.eps ? (depth + 0.6) * level : n.d * level) + 22;
  const out = [`<svg class="g-tree" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">`];
  for (const n of nodes) {
    if (!n.parent) continue;
    const x0 = X(n.parent), y0 = Y(n.parent) + 9, x1 = X(n), y1 = Y(n) - 10;
    const dropped = n.t.leaf && n.d < depth;
    out.push(`<path d="M${x0.toFixed(1)},${y0.toFixed(1)} C${x0.toFixed(1)},${(y0 + 18).toFixed(1)} ${x1.toFixed(1)},${(Math.min(y1, y0 + level) - 16).toFixed(1)} ${x1.toFixed(1)},${(Math.min(y1, y0 + level) - 4).toFixed(1)}${dropped ? ` L${x1.toFixed(1)},${y1.toFixed(1)}` : ''}" style="fill:none;stroke:var(--hair-3);stroke-width:1.1"${dropped ? ' stroke-dasharray="2 3"' : ''}/>`);
  }
  for (const n of nodes) {
    const x = X(n), y = Y(n);
    if (n.t.leaf) {
      const c = n.t.letter;
      out.push(`<g><title>${esc(`${c} at ${n.t.at + 1}`)}</title><circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="9.5" style="fill:${isBase(c) ? BASE[c] : 'var(--accent)'};fill-opacity:0.2;stroke:${isBase(c) ? BASE[c] : 'var(--accent)'};stroke-width:1"/>`
        + `<text x="${x.toFixed(1)}" y="${(y + 4).toFixed(1)}" text-anchor="middle" style="font:13px var(--mono);fill:var(--ink)">${esc(c)}</text></g>`);
    } else if (n.eps) {
      out.push(`<g><title>${esc(n.t.text)}</title><text x="${x.toFixed(1)}" y="${(y - 12).toFixed(1)}" text-anchor="middle" style="font-size:12.5px;font-style:italic;fill:var(--ink-2)">${esc(n.t.lhs)}</text>`
        + `<text x="${x.toFixed(1)}" y="${(y + 4).toFixed(1)}" text-anchor="middle" style="font-size:13px;font-style:italic;fill:var(--ink-faint)">ε</text></g>`);
    } else {
      const w = Math.max(18, 7.5 * n.t.lhs.length + 10);
      out.push(`<g class="g-node"><title>${esc(`${n.t.text}  over ${n.t.begin + 1}..${n.t.end}`)}</title><rect x="${(x - w / 2).toFixed(1)}" y="${(y - 10).toFixed(1)}" width="${w.toFixed(1)}" height="19" rx="9.5" style="fill:var(--surface);stroke:var(--accent-ring)"/>`
        + `<text x="${x.toFixed(1)}" y="${(y + 4).toFixed(1)}" text-anchor="middle" style="font-size:12.5px;font-style:italic;fill:var(--ink)">${esc(n.t.lhs)}</text></g>`);
    }
  }
  out.push('</svg>');
  return { svg: out.join(''), width: W, height: H, leaves: slots };
}

// ---- the bases -------------------------------------------------------------------------------------------
// The Klein group on the bases: A and G above, T and C below, so that transitions (τ) run across,
// complements (κ) up and down, and σ along the diagonals. A subgroup's cosets are drawn as soft
// regions named by their IUPAC letters; the base point is ringed and every base carries the element
// that carries the base point to it.
export function kleinSquareSvg({ basePoint = 'A', torsor = [], cosets = [], size = 300 } = {}) {
  const P = { A: [70, 70], G: [230, 70], T: [70, 230], C: [230, 230] };
  const sc = size / 300;
  const at = (b) => P[b].map((v) => v * sc);
  const edge = (a, b, kind) => {
    const [x0, y0] = at(a), [x1, y1] = at(b);
    const style = { κ: 'stroke:var(--accent);stroke-width:2.4', τ: 'stroke:var(--G);stroke-width:2.2;stroke-dasharray:7 5', σ: 'stroke:var(--ink-faint);stroke-width:1.6;stroke-dasharray:1.5 4.5;stroke-linecap:round' }[kind];
    return `<line x1="${x0}" y1="${y0}" x2="${x1}" y2="${y1}" style="${style}"/>`;
  };
  const out = [`<svg class="g-klein" width="${size}" height="${size}" viewBox="0 0 ${size} ${size}">`];
  // the cosets of the chosen subgroup
  const tint = ['var(--accent-soft)', 'var(--warn-soft)', 'var(--good-soft)', 'var(--bad-soft)'];
  cosets.forEach((c, k) => {
    const pts = [...c.bases].map(at);
    const cx = pts.reduce((s, p) => s + p[0], 0) / pts.length, cy = pts.reduce((s, p) => s + p[1], 0) / pts.length;
    if (pts.length === 1) {
      out.push(`<circle cx="${cx}" cy="${cy}" r="${34 * sc}" style="fill:${tint[k % 4]}"/>`);
    } else if (pts.length === 2) {
      const [a, b] = pts, len = Math.hypot(b[0] - a[0], b[1] - a[1]), ang = (Math.atan2(b[1] - a[1], b[0] - a[0]) * 180) / Math.PI;
      out.push(`<rect x="${cx - len / 2 - 34 * sc}" y="${cy - 30 * sc}" width="${len + 68 * sc}" height="${60 * sc}" rx="${30 * sc}" transform="rotate(${ang} ${cx} ${cy})" style="fill:${tint[k % 4]}"/>`);
    } else {
      out.push(`<rect x="${22 * sc}" y="${22 * sc}" width="${256 * sc}" height="${256 * sc}" rx="${46 * sc}" style="fill:${tint[0]}"/>`);
    }
    // the coset's IUPAC letter: above or below a row, beside a column, inside a diagonal near its
    // upper end, at the centre for N; a single base is its own letter
    if (pts.length > 1) {
      let lx = 150 * sc, ly = 156 * sc;
      if (pts.length === 2) {
        const [a, b] = pts;
        if (a[1] === b[1]) { lx = cx; ly = cy + (cy < 150 * sc ? -40 : 50) * sc; }
        else if (a[0] === b[0]) { lx = cx + (cx < 150 * sc ? -46 : 46) * sc; ly = cy + 6 * sc; }
        else { const up = a[1] < b[1] ? a : b; lx = cx + (up[0] - cx) * 0.48; ly = cy + (up[1] - cy) * 0.48 + 6 * sc; }
      }
      out.push(`<text x="${lx}" y="${ly}" text-anchor="middle" style="font-size:${17 * sc}px;font-weight:600;fill:var(--ink-2);letter-spacing:0.05em">${esc(c.letter)}</text>`);
    }
  });
  out.push(edge('A', 'G', 'τ'), edge('T', 'C', 'τ'), edge('A', 'T', 'κ'), edge('G', 'C', 'κ'), edge('A', 'C', 'σ'), edge('G', 'T', 'σ'));
  for (const b of 'AGTC') {
    const [x, y] = at(b);
    const g = torsor.find((t) => t.base === b);
    const ring = b === basePoint ? `<circle cx="${x}" cy="${y}" r="${25 * sc}" style="fill:none;stroke:var(--ink);stroke-width:1.2;stroke-dasharray:3 3"/>` : '';
    out.push(`${ring}<circle cx="${x}" cy="${y}" r="${19 * sc}" style="fill:var(--surface);stroke:${BASE[b]};stroke-width:2.2"/>`
      + `<text x="${x}" y="${y + 6 * sc}" text-anchor="middle" style="font:600 ${17 * sc}px var(--mono);fill:${BASE[b]}">${b}</text>`);
    if (g) out.push(`<text x="${x}" y="${y + (y < 150 * sc ? -31 : 41) * sc}" text-anchor="middle" style="font-size:${12.5 * sc}px;font-style:italic;fill:var(--ink-faint)">${esc(g.element === '1' ? basePoint : `${g.element}·${basePoint}`)}</text>`);
  }
  out.push('</svg>');
  return out.join('');
}

// A chain of breakings, N at the root, the cosets of one subgroup, the bases at the leaves.
export function breakingTreeSvg(chain, { width = 300, height = 190 } = {}) {
  const leaves = 'ACGT'.split('').map((b, i) => ({ b, parent: chain.parent[i] }));
  // the bases ordered by their coset, so that no edge crosses
  leaves.sort((x, y) => x.parent - y.parent || 'ACGT'.indexOf(x.b) - 'ACGT'.indexOf(y.b));
  const yRoot = 26, yMid = 92, yLeaf = 160;
  const xLeaf = (k) => 40 + (k * (width - 80)) / 3;
  const mids = chain.middle.map((letter, q) => {
    const ks = leaves.map((l, k) => (l.parent === q ? k : -1)).filter((k) => k >= 0);
    return { letter, x: (xLeaf(ks[0]) + xLeaf(ks[ks.length - 1])) / 2 };
  });
  const out = [`<svg class="g-chain" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`];
  const rootX = width / 2;
  for (const m of mids) out.push(`<path d="M${rootX},${yRoot + 12} C${rootX},${yRoot + 40} ${m.x},${yMid - 40} ${m.x},${yMid - 12}" style="fill:none;stroke:var(--hair-3);stroke-width:1.2"/>`);
  leaves.forEach((l, k) => { const m = mids[l.parent]; out.push(`<path d="M${m.x},${yMid + 12} C${m.x},${yMid + 38} ${xLeaf(k)},${yLeaf - 38} ${xLeaf(k)},${yLeaf - 13}" style="fill:none;stroke:var(--hair-3);stroke-width:1.2"/>`); });
  const disc = (x, y, label, colour, r = 13) => `<circle cx="${x}" cy="${y}" r="${r}" style="fill:var(--surface);stroke:${colour};stroke-width:1.6"/><text x="${x}" y="${y + 5}" text-anchor="middle" style="font:600 14px var(--mono);fill:${colour}">${esc(label)}</text>`;
  out.push(disc(rootX, yRoot, 'N', 'var(--ink-2)'));
  for (const m of mids) out.push(disc(m.x, yMid, m.letter, 'var(--accent)'));
  leaves.forEach((l, k) => out.push(disc(xLeaf(k), yLeaf, l.b, BASE[l.b])));
  out.push(`<text x="${width - 6}" y="${yRoot + 4}" text-anchor="end" class="tick-label" style="font-style:italic">V₄ acts on 1 class</text>`,
    `<text x="${width - 6}" y="${yMid - 22}" text-anchor="end" class="tick-label" style="font-style:italic">${esc(chain.residual)} on 2</text>`,
    `<text x="${width - 6}" y="${yLeaf + 30}" text-anchor="end" class="tick-label" style="font-style:italic">V₄ on the 4 bases</text>`);
  out.push('</svg>');
  return out.join('');
}

// ---- bars ---------------------------------------------------------------------------------------------------
// Horizontal bars with their names at the left: rows [{label, sub, value, colour, note, marks}],
// linear from 0 or logarithmic.
export function hbarsSvg(rows, { width = 560, row = 30, label = 170, log = false, min = null, max = null, format = (v) => v.toFixed(3), axis = '' } = {}) {
  const vals = rows.map((r) => r.value).filter((v) => Number.isFinite(v) && (!log || v > 0));
  const hi = max ?? Math.max(...vals, log ? 1e-12 : 1e-12);
  const lo = log ? (min ?? Math.max(Math.min(...vals) / 3, hi * 1e-9)) : 0;
  const L = Math.log10;
  const W = width - label - 70;
  const sx = (v) => (log ? ((L(Math.max(v, lo)) - L(lo)) / (L(hi) - L(lo) || 1)) * W : (v / (hi || 1)) * W);
  const H = rows.length * row + (log ? 34 : 18);
  const out = [`<svg class="g-bars" width="${width}" height="${H}" viewBox="0 0 ${width} ${H}">`];
  if (log) {
    for (let e = Math.ceil(L(lo)); e <= Math.floor(L(hi)); e++) {
      const x = label + sx(10 ** e);
      out.push(`<line x1="${x.toFixed(1)}" y1="2" x2="${x.toFixed(1)}" y2="${rows.length * row + 2}" style="stroke:var(--hair)"/><text class="tick-label" x="${x.toFixed(1)}" y="${rows.length * row + 16}" text-anchor="middle">10${esc(sup(e).replace('-', '⁻'))}</text>`);
    }
    if (axis) out.push(`<text class="tick-label" x="${label + W / 2}" y="${rows.length * row + 30}" text-anchor="middle" style="font-style:italic">${esc(axis)}</text>`);
  }
  rows.forEach((r, k) => {
    const y = k * row + 4;
    const w = Number.isFinite(r.value) && (!log || r.value > 0) ? Math.max(1.5, sx(r.value)) : 0;
    out.push(`<g><title>${esc(r.title || `${r.label}: ${format(r.value)}`)}</title>`
      + `<text x="${label - 10}" y="${y + row / 2 + 1}" text-anchor="end" style="font-size:13.5px;fill:var(--ink)">${esc(r.label)}</text>`
      + (r.sub ? `<text x="${label - 10}" y="${y + row / 2 + 13}" text-anchor="end" style="font-size:10.5px;font-style:italic;fill:var(--ink-faint)">${esc(r.sub)}</text>` : '')
      + `<rect x="${label}" y="${y + 5}" width="${w.toFixed(1)}" height="${row - 12}" rx="3" style="fill:${r.colour || 'var(--accent)'};fill-opacity:${r.opacity ?? 0.78}"/>`
      + `<text x="${label + w + 6}" y="${y + row / 2 + 2}" style="font-size:12px;fill:var(--ink-2);font-variant-numeric:tabular-nums">${esc(r.note ?? format(r.value))}</text></g>`);
  });
  out.push('</svg>');
  return out.join('');
}

// ---- the genetic code on (C₂)⁶ -------------------------------------------------------------------------
// The 64 elements (or characters) of V₄³ as an 8 × 8 grid: blocks of 2 × 2 by the first position
// (rows) and the second (columns), and within a block the third position, its high bit down and its
// low bit across. names: the four labels of a position's two bits (for elements 1, σ, τ, κ; for
// characters ·, RY, MK, WS).
export function codeGridSvg(cells, { names, cell = 34, title = (c) => '', colour = () => 'var(--paper-3)', text = () => '', mark = () => false, selected = -1 } = {}) {
  const gap = 5, pad = 46;
  const pos = (m) => {
    const b1 = (m >> 4) & 3, b2 = (m >> 2) & 3, b3 = m & 3;
    return [pad + b2 * (2 * cell + gap) + (b3 & 1) * cell, pad + b1 * (2 * cell + gap) + (b3 >> 1) * cell];
  };
  const side = pad + 4 * (2 * cell + gap);
  const out = [`<svg class="g-code" width="${side + 8}" height="${side + 30}" viewBox="0 0 ${side + 8} ${side + 30}">`];
  for (let q = 0; q < 4; q++) {
    out.push(`<text x="${pad + q * (2 * cell + gap) + cell}" y="${pad - 12}" text-anchor="middle" style="font-size:13px;font-style:italic;fill:var(--ink-2)">${esc(names[q])}${sub(2)}</text>`);
    out.push(`<text x="${pad - 10}" y="${pad + q * (2 * cell + gap) + cell + 4}" text-anchor="end" style="font-size:13px;font-style:italic;fill:var(--ink-2)">${esc(names[q])}${sub(1)}</text>`);
  }
  for (const c of cells) {
    const [x, y] = pos(c.bits);
    const on = mark(c), sel = c.bits === selected;
    out.push(`<g data-bits="${c.bits}" style="cursor:pointer"><title>${esc(title(c))}</title><rect x="${x + 0.5}" y="${y + 0.5}" width="${cell - 1}" height="${cell - 1}" rx="3" style="fill:${colour(c)};stroke:${sel ? 'var(--ink)' : on ? 'var(--accent-deep)' : 'var(--hair)'};stroke-width:${sel || on ? 2 : 1}"/>`
      + (text(c) !== '' ? `<text x="${x + cell / 2}" y="${y + cell / 2 + 4}" text-anchor="middle" style="font-size:11.5px;fill:var(--ink);font-variant-numeric:tabular-nums;pointer-events:none">${esc(text(c))}</text>` : '') + '</g>');
  }
  // the key of the third position, within a block
  const kx = pad, ky = side + 8;
  out.push(`<text x="${kx}" y="${ky + 12}" class="tick-label" style="font-style:italic">within a block, the third position:  ${esc(names[0])} ${esc(names[1])} above,  ${esc(names[2])} ${esc(names[3])} below</text>`);
  out.push('</svg>');
  return out.join('');
}

// ---- group algebras --------------------------------------------------------------------------------------
// ℂ[G] ≅ ⊕ M_d(ℂ): each block a d × d square of cells on a diagonal, so that their areas sum to |G|.
export function wedderburnSvg(degrees, { unit = 16, label = true } = {}) {
  const total = degrees.reduce((s, d) => s + d, 0);
  const side = total * unit;
  const out = [`<svg class="g-wedderburn" width="${side + 2}" height="${side + 2}" viewBox="-1 -1 ${side + 2} ${side + 2}">`,
    `<rect x="0" y="0" width="${side}" height="${side}" style="fill:none;stroke:var(--hair)"/>`];
  let at = 0;
  for (const d of degrees) {
    const x = at * unit;
    out.push(`<rect x="${x}" y="${x}" width="${d * unit}" height="${d * unit}" rx="2" style="fill:var(--accent-soft);stroke:var(--accent);stroke-width:1.2"/>`);
    for (let k = 1; k < d; k++) out.push(`<line x1="${x + k * unit}" y1="${x}" x2="${x + k * unit}" y2="${x + d * unit}" style="stroke:var(--accent-ring)"/><line x1="${x}" y1="${x + k * unit}" x2="${x + d * unit}" y2="${x + k * unit}" style="stroke:var(--accent-ring)"/>`);
    if (label && d > 1 && d * unit >= 22) out.push(`<text x="${x + (d * unit) / 2}" y="${x + (d * unit) / 2 + 4}" text-anchor="middle" style="font-size:11px;font-style:italic;fill:var(--accent-deep)">M${sub(d)}</text>`);
    at += d;
  }
  out.push('</svg>');
  return out.join('');
}

// ---- hairpins ------------------------------------------------------------------------------------------------
const PAIRS_WC = new Set(['AT', 'TA', 'CG', 'GC']);

// A hairpin as stem and loop: its pairs as rungs, Watson–Crick solid and wobble dotted, the loop
// on a circle closed by the innermost pair.
export function stemLoopSvg(h, { gap = 24, r = 10 } = {}) {
  const s = h.stem, L = h.loop, letters = h.letters;
  const k = L + 2;  // the loop's bases and the closing pair on one circle
  const R = gap / (2 * Math.sin(Math.PI / k));
  const W = Math.max(2 * R + 2 * r + 60, gap + 2 * r + 110), cx = W / 2;
  const loopTop = 18 + r, cy = loopTop + R;
  const stemTop = cy + R * Math.cos(Math.PI / k);
  const H = stemTop + (s - 1) * gap + r + 30;
  const pts = new Array(letters.length);
  // the stem: left strand climbing to the loop, right strand coming down from it
  for (let q = 0; q < s; q++) {
    const y = stemTop + (s - 1 - q) * gap;
    pts[q] = [cx - gap / 2, y];
    pts[letters.length - 1 - q] = [cx + gap / 2, y];
  }
  // the loop: the vertices of a regular k-gon after the closing pair, clockwise over the top
  const a0 = Math.PI / 2 + Math.PI / k;
  for (let q = 0; q < L; q++) {
    const a = a0 + ((q + 1) * 2 * Math.PI) / k;
    pts[s + q] = [cx + R * Math.cos(a), cy + R * Math.sin(a)];
  }
  const out = [`<svg class="g-hairpin" width="${W.toFixed(0)}" height="${H.toFixed(0)}" viewBox="0 0 ${W.toFixed(0)} ${H.toFixed(0)}">`];
  // the backbone
  out.push(`<path d="M${pts.map((p) => `${p[0].toFixed(1)},${p[1].toFixed(1)}`).join(' L')}" style="fill:none;stroke:var(--hair-3);stroke-width:1.4"/>`);
  for (let q = 0; q < s; q++) {
    const [a, b] = [pts[q], pts[letters.length - 1 - q]];
    const wc = PAIRS_WC.has(letters[q] + letters[letters.length - 1 - q]);
    out.push(`<line x1="${a[0] + r}" y1="${a[1]}" x2="${b[0] - r}" y2="${b[1]}" style="stroke:var(--accent);stroke-width:${wc ? 2 : 1.6}"${wc ? '' : ' stroke-dasharray="2 3"'}/>`);
  }
  letters.split('').forEach((c, q) => {
    const [x, y] = pts[q];
    out.push(`<circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="${r}" style="fill:var(--surface);stroke:${BASE[c] || BASE.N};stroke-width:1.6"/><text x="${x.toFixed(1)}" y="${(y + 4).toFixed(1)}" text-anchor="middle" style="font:600 12px var(--mono);fill:${BASE[c] || BASE.N}">${esc(c)}</text>`);
  });
  out.push(`<text x="${cx - gap / 2 - r - 4}" y="${H - 10}" text-anchor="end" class="tick-label">5′ ${h.begin + 1}</text><text x="${cx + gap / 2 + r + 4}" y="${H - 10}" class="tick-label">${h.end} 3′</text>`);
  out.push('</svg>');
  return out.join('');
}

// A sequence as a line, its hairpins as arcs over it: each pair an arc from one base to its partner,
// so that a stem is a nest of arcs. The chosen hairpin is drawn in the accent, the others quieter.
export function arcSvg(letters, hairpins, { begin = 0, selected = -1, width = 980 } = {}) {
  const n = letters.length;
  const step = (width - 40) / Math.max(1, n);
  const X = (p) => 20 + (p - begin + 0.5) * step;
  const longest = Math.max(1, ...hairpins.map((h) => h.length));
  const top = Math.min(150, (longest * step) / 2 + 20);
  const base = top + 16;
  const H = base + (step >= 9 ? 40 : 30);
  const out = [`<svg class="g-arcs" width="${width}" height="${H.toFixed(0)}" viewBox="0 0 ${width} ${H.toFixed(0)}">`];
  out.push(`<line x1="${X(begin) - step / 2}" y1="${base}" x2="${X(begin + n - 1) + step / 2}" y2="${base}" style="stroke:var(--hair-2)"/>`);
  hairpins.forEach((h, k) => {
    const on = k === selected;
    out.push(`<g data-hairpin="${k}" style="cursor:pointer"><title>${esc(`${h.begin + 1}..${h.end}: ${h.stem} pairs, a loop of ${h.loop}`)}</title>`);
    for (const [i, j] of h.pairs) {
      const x0 = X(i), x1 = X(j), rr = (x1 - x0) / 2;
      out.push(`<path d="M${x0.toFixed(1)},${base - 2} A${rr.toFixed(1)},${Math.min(rr, top - 4).toFixed(1)} 0 0 1 ${x1.toFixed(1)},${base - 2}" style="fill:none;stroke:${on ? 'var(--accent)' : 'var(--ink-faint)'};stroke-opacity:${on ? 0.95 : 0.45};stroke-width:${on ? 1.6 : 1.1}"/>`);
    }
    // a wide transparent band over the span, to choose it
    out.push(`<rect x="${(X(h.begin) - step / 2).toFixed(1)}" y="${(base - Math.min(top, ((h.end - h.begin) * step) / 2) - 4).toFixed(1)}" width="${((h.end - h.begin) * step).toFixed(1)}" height="${(Math.min(top, ((h.end - h.begin) * step) / 2) + 4).toFixed(1)}" style="fill:transparent"/></g>`);
  });
  if (step >= 9) {
    letters.split('').forEach((c, q) => out.push(`<text x="${X(begin + q).toFixed(1)}" y="${base + 15}" text-anchor="middle" style="font:${Math.min(12, step + 1).toFixed(0)}px var(--mono);fill:${BASE[c] || BASE.N}">${esc(c)}</text>`));
  } else {
    letters.split('').forEach((c, q) => out.push(`<rect x="${(X(begin + q) - step / 2).toFixed(2)}" y="${base + 4}" width="${Math.max(0.6, step).toFixed(2)}" height="7" style="fill:${BASE[c] || BASE.N};fill-opacity:0.8"/>`));
  }
  out.push(`<text class="tick-label" x="20" y="${H - 4}">${begin + 1}</text><text class="tick-label" x="${width - 20}" y="${H - 4}" text-anchor="end">${begin + n}</text>`);
  out.push('</svg>');
  return out.join('');
}
