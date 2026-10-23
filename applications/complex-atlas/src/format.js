// How the Atlas writes what the engine computes: groups as ℤ² ⊕ ℤ/2, presentations as
// ⟨ a, b | a b a⁻¹ b⁻¹ ⟩, words with powers, chains as sums of cells, surfaces by name.

const SUP = { 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹', '-': '⁻' };
export const superscript = (n) => String(n).split('').map((c) => SUP[c] ?? c).join('');

// A finitely generated abelian group from its rank and torsion coefficients.
export function groupText(rank, torsion = []) {
  const parts = [];
  if (rank === 1) parts.push('ℤ');
  else if (rank > 1) parts.push(`ℤ${superscript(rank)}`);
  const counts = new Map();
  for (const t of torsion) counts.set(t, (counts.get(t) || 0) + 1);
  for (const [t, k] of [...counts].sort((a, b) => a[0] - b[0])) parts.push(k === 1 ? `ℤ/${t}` : `(ℤ/${t})${superscript(k)}`);
  return parts.length ? parts.join(' ⊕ ') : '0';
}

// A word in named generators (signed 1-based letters), runs of one letter written as powers:
// [1, 1, 2, -1] -> "a² b a⁻¹". The empty word is 1.
export function wordText(word, names, { powers = true } = {}) {
  if (!word || !word.length) return '1';
  const out = [];
  for (let i = 0; i < word.length;) {
    let j = i;
    if (powers) while (j + 1 < word.length && word[j + 1] === word[i]) j++;
    const k = word[i], n = j - i + 1, name = names[Math.abs(k) - 1] ?? `x${Math.abs(k)}`;
    const e = (k < 0 ? -1 : 1) * n;
    out.push(e === 1 ? name : `${name}${superscript(e)}`);
    i = j + 1;
  }
  return out.join(' ');
}

// ⟨ a, b | a b a⁻¹ b⁻¹ ⟩; the trivial group is ⟨ | ⟩.
export function presentationText(generators, relators) {
  const rels = relators.map((r) => wordText(r, generators)).join(', ');
  return `⟨ ${generators.join(', ')}${generators.length ? ' ' : ''}|${rels ? ` ${rels}` : ''} ⟩`;
}

// The names of the letters of a surface's normal form: a_i, b_i for an orientable surface, a_i
// for the crosscaps of a nonorientable one; with subscripts only when there is more than one.
export function normalFormNames(surface) {
  const n = Math.max(0, ...(surface.normal_form || []).map((k) => Math.abs(k)));
  const sub = (i) => String(i).split('').map((c) => '₀₁₂₃₄₅₆₇₈₉'[c]).join('');
  if (surface.orientable) {
    const g = n / 2;
    return Array.from({ length: n }, (_, i) => `${i % 2 ? 'b' : 'a'}${g > 1 ? sub(Math.floor(i / 2) + 1) : ''}`);
  }
  return Array.from({ length: n }, (_, i) => `a${n > 1 ? sub(i + 1) : ''}`);
}

export function surfaceText(s) {
  if (!s) return '';
  if (s.name) return s.name[0].toUpperCase() + s.name.slice(1);
  const kind = s.orientable ? `orientable surface of genus ${s.genus}` : `nonorientable surface with ${s.crosscaps} crosscaps`;
  return (s.boundary_components ? `${kind}, with ${s.boundary_components} boundary components` : kind).replace(/^./, (c) => c.toUpperCase());
}

const NOUNS = [['vertex', 'vertices'], ['edge', 'edges'], ['face', 'faces']];
export function countsText(counts) {
  return counts.map((n, k) => {
    const [one, many] = NOUNS[k] || [`${k}-cell`, `${k}-cells`];
    return `${n} ${n === 1 ? one : many}`;
  }).join(' · ');
}

// The name of a cell: a vertex by its number, an edge by its name, a face as f1, a k-cell as c1.
export function cellName(dim, index, edgeNames) {
  if (dim === 0) return `v${index}`;
  if (dim === 1) return edgeNames[index] ?? `e${index + 1}`;
  if (dim === 2) return `f${index + 1}`;
  return `c${index + 1}`;
}

// A chain [[cell, coefficient], ...] of dimension dim as a sum: "a + b − 2 f3".
export function chainText(chain, dim, edgeNames, { limit = 8 } = {}) {
  if (!chain.length) return '0';
  const terms = chain.slice(0, limit).map(([c, k], i) => {
    const name = cellName(dim, c, edgeNames), m = Math.abs(k);
    const body = m === 1 ? name : `${m} ${name}`;
    return i === 0 ? (k < 0 ? `−${body}` : body) : `${k < 0 ? '−' : '+'} ${body}`;
  });
  return terms.join(' ') + (chain.length > limit ? ` + … (${chain.length} cells)` : '');
}

// A number for reading: a few significant digits, a minus sign, and ± an error when given.
export function number(x, digits = 4) {
  if (x === null || x === undefined || Number.isNaN(x)) return '—';
  if (Number.isInteger(x)) return String(x).replace('-', '−');
  const a = Math.abs(x);
  const text = a !== 0 && (a < 1e-3 || a >= 1e6) ? x.toExponential(digits - 1) : x.toPrecision(digits).replace(/(\.\d*?)0+$/, '$1').replace(/\.$/, '');
  return text.replace('-', '−');
}
export const withError = (x, e) => (e ? `${number(x)} ± ${number(e, 2)}` : number(x));
