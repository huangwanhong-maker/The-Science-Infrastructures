// Cell complexes, drawn: an abstract CW complex realised in three dimensions, its 2-cells
// unfolded into a net, and its cells as a Hasse diagram; and the conversions from surface
// words and lists of facets to the Complex JSON of The Science Infrastructures.
//
// A complex (the protocol of applications/complex-atlas/README.md) has `vertices` (a count),
// `edges` ([tail, head] pairs), `faces` ({base, boundary}: a closed edge path from the base
// vertex, as signed 1-based edge indices, k for edge k-1 forward and -k backward), and
// `cells` (cells[k-3] the k-cells, each with its boundary [[index of a (k-1)-cell, coefficient]]);
// optionally `positions`, `labels`, `simplices`.
//
// - realise(): a complex has no geometry of its own, so one is made for it. Its vertices, the
//   midpoints of its edges and the centres of its faces are nodes of a graph (a node of an
//   edge joined to its ends, a node of a face to its edges and vertices), and the nodes are
//   placed in space by stress majorisation (Gansner, Koren and North 2004, Graph drawing by
//   stress majorization, LNCS 3383, 239-250): their distances in space are fitted to their
//   distances in the graph, weighted by the inverse square, from the initial positions of
//   pivot MDS (Brandes and Pich 2006, Eigensolver methods for progressive multidimensional
//   scaling of large data, LNCS 4372, 42-53). Each edge is a curve through its midpoint, so
//   that loops and parallel edges stay apart; each face is a disc spanned on its boundary,
//   fanned from its centre and relaxed toward a harmonic (soap-film) surface by Jacobi
//   smoothing of its interior rings. A face attached to a point is a balloon at that point.
// - net(): each 2-cell is a regular polygon with unit sides, read anticlockwise along its
//   attaching word; polygons are glued side to side along a spanning tree of the dual graph
//   wherever the gluing does not overlap what is already placed (separating axes for convex
//   polygons), as surfaces are cut and unfolded to classify them (Massey 1967, Algebraic
//   Topology: An Introduction, ch. 1). A one-face surface is its fundamental polygon.
// - hasse(): the cells by dimension, each joined to the cells of its boundary, ordered in each
//   layer by the barycentre heuristic (Sugiyama, Tagawa and Toda 1981, IEEE Trans. Syst. Man
//   Cybern. 11, 109-125).
//
// Everything here is deterministic: a layout depends only on the complex and its seed.

// ---- the complex --------------------------------------------------------------------------------

const letterEdge = (k) => Math.abs(k) - 1;
const faceList = (c) => c.faces || [];

// The number of cells in each dimension: [vertices, edges, faces, 3-cells, ...].
export function cellCounts(c) {
  const counts = [c.vertices, c.edges.length];
  const f = faceList(c).length, higher = c.cells || [];
  if (f || higher.some((h) => h.length)) counts.push(f);
  for (const h of higher) counts.push(h.length);
  while (counts.length > 1 && counts[counts.length - 1] === 0 && counts.length > 2) counts.pop();
  return counts;
}

export const eulerCharacteristic = (c) => cellCounts(c).reduce((s, n, k) => s + (k % 2 ? -n : n), 0);

// The names of the edges: the complex's labels, else a, b, c, ... (up to 26 edges) or e1, e2, ...
export function edgeNames(c) {
  const given = c.labels?.edges;
  return c.edges.map((_, i) => (given && given[i] != null && given[i] !== '' ? String(given[i])
    : c.edges.length <= 26 ? String.fromCharCode(97 + i) : `e${i + 1}`));
}

// A word (signed 1-based letters) as text: "a b a⁻¹ b⁻¹"; the empty word as "1".
export function wordText(word, names) {
  if (!word.length) return '1';
  return word.map((k) => `${names[letterEdge(k)] ?? `e${Math.abs(k)}`}${k < 0 ? '⁻¹' : ''}`).join(' ');
}

// What is wrong with a complex, as sentences; empty when it is well formed: indices in range,
// every attaching word a closed path from its base, positions finite, and ∂∂ = 0 above.
export function checkComplex(c) {
  const out = [];
  if (!c || typeof c !== 'object') return ['a complex is a JSON object'];
  const V = c.vertices;
  if (!Number.isInteger(V) || V < 0) return ['"vertices" must be a count, a whole number at least 0'];
  if (!Array.isArray(c.edges)) return ['"edges" must be a list of [tail, head] pairs'];
  c.edges.forEach((e, i) => {
    if (!Array.isArray(e) || e.length !== 2 || !e.every((v) => Number.isInteger(v) && v >= 0 && v < V)) out.push(`edge ${i + 1} must be [tail, head] with vertices below ${V}`);
  });
  if (out.length) return out;
  const faces = faceList(c);
  if (!Array.isArray(faces)) return ['"faces" must be a list'];
  faces.forEach((f, i) => {
    if (!f || !Number.isInteger(f.base) || f.base < 0 || f.base >= V) { out.push(`face ${i + 1} needs a base vertex below ${V}`); return; }
    if (!Array.isArray(f.boundary)) { out.push(`face ${i + 1} needs a boundary word`); return; }
    let at = f.base;
    for (const k of f.boundary) {
      if (!Number.isInteger(k) || k === 0 || Math.abs(k) > c.edges.length) { out.push(`face ${i + 1}: ${k} is not a signed edge index between 1 and ${c.edges.length}`); return; }
      const [t, h] = c.edges[letterEdge(k)];
      const from = k > 0 ? t : h;
      if (from !== at) { out.push(`face ${i + 1}: its boundary is not a path (edge ${Math.abs(k)} does not start at vertex ${at})`); return; }
      at = k > 0 ? h : t;
    }
    if (at !== f.base) out.push(`face ${i + 1}: its boundary does not close at its base`);
  });
  if (out.length) return out;
  // ∂∂ = 0 for the higher cells, the faces' boundaries taken abelian.
  let below = faces.map((f) => { const m = new Map(); for (const k of f.boundary) m.set(letterEdge(k), (m.get(letterEdge(k)) || 0) + Math.sign(k)); return m; });
  (c.cells || []).forEach((layer, j) => {
    const k = j + 3;
    const next = [];
    (layer || []).forEach((cell, i) => {
      const sum = new Map();
      for (const [idx, coef] of cell.boundary || []) {
        if (!Number.isInteger(idx) || idx < 0 || idx >= below.length || !Number.isInteger(coef)) { out.push(`${k}-cell ${i + 1}: [${idx}, ${coef}] does not name a ${k - 1}-cell`); continue; }
        for (const [e, v] of below[idx]) sum.set(e, (sum.get(e) || 0) + coef * v);
      }
      if ([...sum.values()].some((v) => v !== 0)) out.push(`${k}-cell ${i + 1}: the boundary of its boundary is not zero`);
      const m = new Map();
      for (const [idx, coef] of cell.boundary || []) m.set(idx, (m.get(idx) || 0) + coef);
      next.push(m);
    });
    below = next;
  });
  if (c.positions) {
    if (!Array.isArray(c.positions) || c.positions.length !== V || !c.positions.every((p) => Array.isArray(p) && p.length >= 2 && p.every(Number.isFinite)))
      out.push(`"positions" must give [x, y, z] for each of the ${V} vertices`);
  }
  return out;
}

// ---- surface words ------------------------------------------------------------------------------

// Polygons with their sides identified, written as words: "abAB" (a torus), "abAB cdCD" (two),
// "aa" (the projective plane), "a1b1A1B1a2b2A2B2". A letter with its digits is a generator;
// lowercase reads it forward, uppercase backward, as does a trailing ' or ⁻¹ or ^-1. Polygons
// are separated by spaces, commas, semicolons or lines; "1" is a disc attached to a point.
// The corners of the polygons are identified as the sides force (by union and find).
export function parseWords(text) {
  const polygons = String(text).trim().split(/[\s,;]+/).filter(Boolean);
  if (!polygons.length) throw new Error('no words: write a polygon as letters, such as abAB');
  const names = [], index = new Map(), words = [];
  for (const p of polygons) {
    if (p === '1') { words.push([]); continue; }
    const word = [];
    let i = 0;
    while (i < p.length) {
      const ch = p[i];
      if (!/[A-Za-z]/.test(ch)) throw new Error(`"${p}": ${ch} is not a letter`);
      let j = i + 1;
      while (j < p.length && /[0-9]/.test(p[j])) j++;
      let name = p.slice(i, j);
      let sign = ch === ch.toLowerCase() ? 1 : -1;
      name = name.toLowerCase();
      for (const suffix of ["'", '⁻¹', '^-1']) if (p.startsWith(suffix, j)) { sign = -sign; j += suffix.length; break; }
      if (!index.has(name)) { index.set(name, names.length); names.push(name); }
      word.push(sign * (index.get(name) + 1));
      i = j;
    }
    words.push(word);
  }
  // Nodes: each generator's tail and head, then each polygon's corners.
  const parent = [];
  const make = () => { parent.push(parent.length); return parent.length - 1; };
  const find = (x) => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
  const unite = (a, b) => { a = find(a); b = find(b); if (a !== b) parent[Math.max(a, b)] = Math.min(a, b); };
  const tail = names.map(() => make()), head = names.map(() => make());
  const bases = words.map((w) => {
    const corners = (w.length ? w : [0]).map(() => make());
    w.forEach((k, s) => {
      const g = letterEdge(k), a = corners[s], b = corners[(s + 1) % w.length];
      if (k > 0) { unite(a, tail[g]); unite(b, head[g]); } else { unite(a, head[g]); unite(b, tail[g]); }
    });
    return corners[0];
  });
  const vertexOf = new Map();
  const vertex = (node) => { const r = find(node); if (!vertexOf.has(r)) vertexOf.set(r, vertexOf.size); return vertexOf.get(r); };
  const faces = words.map((w, i) => ({ base: vertex(bases[i]), boundary: w }));
  const edges = names.map((_, g) => [vertex(tail[g]), vertex(head[g])]);
  return { name: 'words', title: polygons.join(' '), vertices: vertexOf.size, edges, faces, labels: { edges: names } };
}

// ---- facets -------------------------------------------------------------------------------------

// A simplicial complex from its facets: "0 1 2; 0 2 3" (or one facet per line). Vertices named by
// whole numbers keep their order; other names are numbered as they first appear. Every face of
// every facet becomes a cell: edges [a, b] with a < b oriented from a to b; a triangle [a, b, c]
// attached along a b, b c, then a c backward, so that ∂[abc] = [bc] − [ac] + [ab]; and each
// k-simplex, k >= 3, with the boundary Σ (−1)^i (its face without the i-th vertex).
export function parseFacets(text) {
  const rows = String(text).split(/[;\n]+/).map((r) => r.trim().split(/[\s,]+/).filter(Boolean)).filter((r) => r.length);
  if (!rows.length) throw new Error('no facets: write each facet as its vertices, such as 0 1 2; 0 2 3');
  const tokens = [...new Set(rows.flat())];
  const numeric = tokens.every((t) => /^-?\d+$/.test(t));
  if (numeric) tokens.sort((a, b) => Number(a) - Number(b));
  const id = new Map(tokens.map((t, i) => [t, i]));
  const key = (s) => s.join(',');
  const distinct = new Map();
  for (const r of rows) { const f = [...new Set(r.map((t) => id.get(t)))].sort((a, b) => a - b); distinct.set(key(f), f); }
  let facets = [...distinct.values()];
  if (facets.some((f) => f.length > 12)) throw new Error('a facet of more than 12 vertices (dimension 11) is too large to expand');
  const byDim = [], proper = new Set();
  for (const f of facets) {
    const n = f.length;
    for (let mask = 1; mask < 1 << n; mask++) {
      const s = f.filter((_, i) => mask & (1 << i)), k = key(s);
      (byDim[s.length - 1] ??= new Map()).set(k, s);
      if (s.length < n) proper.add(k);
    }
  }
  // Only the maximal facets are kept as "simplices".
  facets = facets.filter((f) => !proper.has(key(f)));
  const lists = byDim.map((m) => [...m.values()].sort((a, b) => { for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return a[i] - b[i]; return 0; }));
  const index = lists.map((l) => new Map(l.map((s, i) => [key(s), i])));
  const edges = (lists[1] || []).map(([a, b]) => [a, b]);
  const e = (a, b) => index[1].get(key([a, b])) + 1;
  const faces = (lists[2] || []).map(([a, b, c]) => ({ base: a, boundary: [e(a, b), e(b, c), -e(a, c)] }));
  const cells = [];
  for (let d = 3; d < lists.length; d++)
    cells.push(lists[d].map((s) => ({ boundary: s.map((_, i) => [index[d - 1].get(key(s.filter((__, j) => j !== i))), i % 2 ? -1 : 1]) })));
  const out = { name: 'facets', title: `${facets.length} facets`, vertices: tokens.length, edges, faces, simplices: facets, labels: { vertices: tokens } };
  if (cells.length) out.cells = cells;
  return out;
}

// ---- small vector algebra -----------------------------------------------------------------------

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = (a) => Math.hypot(a[0], a[1], a[2]);
const unit = (a, fallback = [1, 0, 0]) => { const n = norm(a); return n > 1e-12 ? mul(a, 1 / n) : fallback; };
// A unit vector perpendicular to a.
function perpendicular(a) {
  const u = unit(a);
  const axis = Math.abs(u[0]) < 0.6 ? [1, 0, 0] : Math.abs(u[1]) < 0.6 ? [0, 1, 0] : [0, 0, 1];
  return unit(cross(u, axis));
}
// A reproducible sequence of uniform numbers in [0, 1), from a seed (mulberry32).
function seeded(seed) {
  let s = (seed >>> 0) || 1;
  return () => {
    s = (s + 0x6d2b79f5) >>> 0;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// ---- layout by stress majorisation ----------------------------------------------------------------

// Distances in the graph from every node of a component (BFS), as rows of a Uint16Array.
function distances(adj, nodes) {
  const n = nodes.length, local = new Map(nodes.map((v, i) => [v, i]));
  const D = new Uint16Array(n * n).fill(65535);
  const queue = new Int32Array(n);
  for (let s = 0; s < n; s++) {
    const row = s * n;
    D[row + s] = 0;
    let head = 0, tail = 0;
    queue[tail++] = s;
    while (head < tail) {
      const u = queue[head++];
      for (const w of adj[nodes[u]]) {
        const j = local.get(w);
        if (D[row + j] === 65535) { D[row + j] = D[row + u] + 1; queue[tail++] = j; }
      }
    }
  }
  return D;
}

// Initial coordinates by pivot MDS: double-centred squared distances to k pivots, and the
// three leading eigenvectors of CᵀC by power iteration with deflation.
function pivotMds(D, n, random) {
  const k = Math.min(n, 40);
  const pivots = [0];
  const minDist = new Float64Array(n).fill(Infinity);
  while (pivots.length < k) {
    const p = pivots[pivots.length - 1];
    let best = -1, far = -1;
    for (let i = 0; i < n; i++) {
      minDist[i] = Math.min(minDist[i], D[p * n + i]);
      if (minDist[i] > far) { far = minDist[i]; best = i; }
    }
    if (far <= 0) break;
    pivots.push(best);
  }
  const m = pivots.length;
  const C = new Float64Array(n * m);
  const colMean = new Float64Array(m), rowMean = new Float64Array(n);
  let total = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) {
    const d = D[pivots[j] * n + i];
    const v = d * d;
    C[i * m + j] = v; colMean[j] += v / n; rowMean[i] += v / m; total += v / (n * m);
  }
  for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) C[i * m + j] = -0.5 * (C[i * m + j] - rowMean[i] - colMean[j] + total);
  // M = CᵀC, m × m.
  const M = new Float64Array(m * m);
  for (let a = 0; a < m; a++) for (let b = a; b < m; b++) {
    let s = 0;
    for (let i = 0; i < n; i++) s += C[i * m + a] * C[i * m + b];
    M[a * m + b] = M[b * m + a] = s;
  }
  const vecs = [];
  for (let d = 0; d < 3; d++) {
    let v = Array.from({ length: m }, () => random() - 0.5);
    for (let it = 0; it < 120; it++) {
      const w = new Array(m).fill(0);
      for (let a = 0; a < m; a++) { let s = 0; for (let b = 0; b < m; b++) s += M[a * m + b] * v[b]; w[a] = s; }
      for (const u of vecs) { const p = w.reduce((s, x, i) => s + x * u[i], 0); for (let a = 0; a < m; a++) w[a] -= p * u[a]; }
      const len = Math.hypot(...w);
      if (len < 1e-12) break;
      v = w.map((x) => x / len);
    }
    vecs.push(v);
  }
  const X = new Float64Array(n * 3);
  for (let i = 0; i < n; i++) for (let d = 0; d < 3; d++) {
    let s = 0;
    for (let j = 0; j < m; j++) s += C[i * m + j] * vecs[d][j];
    X[i * 3 + d] = s;
  }
  // Scale so that adjacent nodes sit about one unit apart, and break exact symmetries.
  let sum = 0, count = 0;
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) if (D[i * n + j] === 1) {
    sum += Math.hypot(X[i * 3] - X[j * 3], X[i * 3 + 1] - X[j * 3 + 1], X[i * 3 + 2] - X[j * 3 + 2]); count++;
  }
  const scale = count && sum > 1e-9 ? count / sum : 1;
  for (let i = 0; i < 3 * n; i++) X[i] = X[i] * scale + (random() - 0.5) * 0.05;
  return X;
}

// Stress majorisation, localised and in place (each node moved to the weighted mean of where
// each other node would have it), until the largest move is small.
function majorise(D, n, X, iterations, random) {
  for (let it = 0; it < iterations; it++) {
    let moved = 0;
    for (let i = 0; i < n; i++) {
      let nx = 0, ny = 0, nz = 0, den = 0;
      const xi = X[i * 3], yi = X[i * 3 + 1], zi = X[i * 3 + 2];
      for (let j = 0; j < n; j++) {
        if (j === i) continue;
        const d = D[i * n + j];
        if (d === 65535) continue;
        const w = 1 / (d * d);
        let dx = xi - X[j * 3], dy = yi - X[j * 3 + 1], dz = zi - X[j * 3 + 2];
        let len = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (len < 1e-9) { dx = random() - 0.5; dy = random() - 0.5; dz = random() - 0.5; len = Math.sqrt(dx * dx + dy * dy + dz * dz); }
        const f = d / len;
        nx += w * (X[j * 3] + f * dx); ny += w * (X[j * 3 + 1] + f * dy); nz += w * (X[j * 3 + 2] + f * dz);
        den += w;
      }
      if (!den) continue;
      nx /= den; ny /= den; nz /= den;
      moved = Math.max(moved, Math.sqrt((nx - xi) ** 2 + (ny - yi) ** 2 + (nz - zi) ** 2));
      X[i * 3] = nx; X[i * 3 + 1] = ny; X[i * 3 + 2] = nz;
    }
    if (moved < 1e-3) break;
  }
}

// A large component: distances from k pivots only (k breadth-first searches), pivot MDS on
// them, then stress majorisation restricted to each node's neighbours within two steps (the
// sparse stress of Ortmann, Klimenta and Brandes 2016, J. Graph Algorithms Appl. 21, 251-276,
// without its pivot terms).
function sparseLayout(adj, nodes, iterations, random) {
  const n = nodes.length, local = new Map(nodes.map((v, i) => [v, i]));
  const nbr = nodes.map((v) => adj[v].map((w) => local.get(w)));
  const bfs = (s) => {
    const d = new Uint16Array(n).fill(65535), q = new Int32Array(n);
    let h = 0, t = 0;
    d[s] = 0; q[t++] = s;
    while (h < t) { const u = q[h++]; for (const w of nbr[u]) if (d[w] === 65535) { d[w] = d[u] + 1; q[t++] = w; } }
    return d;
  };
  const k = Math.min(n, 40), rows = [], pivots = [0];
  const minDist = new Float64Array(n).fill(Infinity);
  while (rows.length < k) {
    const d = bfs(pivots[pivots.length - 1]);
    rows.push(d);
    let best = 0, far = -1;
    for (let i = 0; i < n; i++) { minDist[i] = Math.min(minDist[i], d[i]); if (minDist[i] > far) { far = minDist[i]; best = i; } }
    if (rows.length < k) pivots.push(best);
  }
  const m = rows.length, C = new Float64Array(n * m), colMean = new Float64Array(m), rowMean = new Float64Array(n);
  let total = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) { const v = rows[j][i] ** 2; C[i * m + j] = v; colMean[j] += v / n; rowMean[i] += v / m; total += v / (n * m); }
  for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) C[i * m + j] = -0.5 * (C[i * m + j] - rowMean[i] - colMean[j] + total);
  const vecs = [];
  for (let d = 0; d < 3; d++) {
    let v = Array.from({ length: m }, () => random() - 0.5);
    for (let it = 0; it < 100; it++) {
      const Cv = new Float64Array(n);
      for (let i = 0; i < n; i++) { let s = 0; for (let j = 0; j < m; j++) s += C[i * m + j] * v[j]; Cv[i] = s; }
      const w = new Array(m).fill(0);
      for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) w[j] += C[i * m + j] * Cv[i];
      for (const u of vecs) { const p = w.reduce((s, x, i) => s + x * u[i], 0); for (let a = 0; a < m; a++) w[a] -= p * u[a]; }
      const len = Math.hypot(...w);
      if (len < 1e-12) break;
      v = w.map((x) => x / len);
    }
    vecs.push(v);
  }
  const X = new Float64Array(n * 3);
  for (let i = 0; i < n; i++) for (let d = 0; d < 3; d++) { let s = 0; for (let j = 0; j < m; j++) s += C[i * m + j] * vecs[d][j]; X[i * 3 + d] = s; }
  let sum = 0, count = 0;
  for (let i = 0; i < n; i++) for (const j of nbr[i]) { sum += Math.hypot(X[i * 3] - X[j * 3], X[i * 3 + 1] - X[j * 3 + 1], X[i * 3 + 2] - X[j * 3 + 2]); count++; }
  const scale = count && sum > 1e-9 ? count / sum : 1;
  for (let i = 0; i < 3 * n; i++) X[i] = X[i] * scale + (random() - 0.5) * 0.05;
  // Two-step neighbourhoods with their distances.
  const near = nbr.map((list, i) => {
    const m2 = new Map(list.map((j) => [j, 1]));
    for (const j of list) for (const l of nbr[j]) if (l !== i && !m2.has(l)) m2.set(l, 2);
    return [...m2];
  });
  for (let it = 0; it < iterations; it++) {
    let moved = 0;
    for (let i = 0; i < n; i++) {
      let nx = 0, ny = 0, nz = 0, den = 0;
      for (const [j, d] of near[i]) {
        const w = 1 / (d * d);
        let dx = X[i * 3] - X[j * 3], dy = X[i * 3 + 1] - X[j * 3 + 1], dz = X[i * 3 + 2] - X[j * 3 + 2];
        let len = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (len < 1e-9) { dx = random() - 0.5; dy = random() - 0.5; dz = random() - 0.5; len = Math.sqrt(dx * dx + dy * dy + dz * dz); }
        nx += w * (X[j * 3] + (d / len) * dx); ny += w * (X[j * 3 + 1] + (d / len) * dy); nz += w * (X[j * 3 + 2] + (d / len) * dz);
        den += w;
      }
      if (!den) continue;
      nx /= den; ny /= den; nz /= den;
      moved = Math.max(moved, Math.hypot(nx - X[i * 3], ny - X[i * 3 + 1], nz - X[i * 3 + 2]));
      X[i * 3] = nx; X[i * 3 + 1] = ny; X[i * 3 + 2] = nz;
    }
    if (moved < 1e-3) break;
  }
  return X;
}

// Positions of the nodes (in `adj`) of each connected component, laid out alone and then set
// side by side along x. Components of more than `dense` nodes are laid out sparsely.
function layoutGraph(adj, iterations, random, dense = 900) {
  const N = adj.length, seen = new Uint8Array(N), out = new Array(N);
  let offset = 0;
  for (let s = 0; s < N; s++) {
    if (seen[s]) continue;
    const nodes = [s];
    seen[s] = 1;
    for (let h = 0; h < nodes.length; h++) for (const w of adj[nodes[h]]) if (!seen[w]) { seen[w] = 1; nodes.push(w); }
    const n = nodes.length;
    let X;
    if (n === 1) X = new Float64Array([0, 0, 0]);
    else if (n > dense) X = sparseLayout(adj, nodes, Math.min(iterations, 80), random);
    else {
      const D = distances(adj, nodes);
      X = pivotMds(D, n, random);
      majorise(D, n, X, iterations, random);
    }
    let lo = Infinity, hi = -Infinity, cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < n; i++) { cx += X[i * 3] / n; cy += X[i * 3 + 1] / n; cz += X[i * 3 + 2] / n; }
    for (let i = 0; i < n; i++) { lo = Math.min(lo, X[i * 3] - cx); hi = Math.max(hi, X[i * 3] - cx); }
    const shift = offset - lo;
    nodes.forEach((v, i) => { out[v] = [X[i * 3] - cx + shift, X[i * 3 + 1] - cy, X[i * 3 + 2] - cz]; });
    offset += hi - lo + 2;
  }
  return out;
}

// ---- realisation in three dimensions ---------------------------------------------------------------

// The geometry that shows a complex in space. Options: seed, iterations (of the layout), samples
// (points per edge curve), rings and smoothing (of the face discs), full (the most nodes laid
// out with the edges and faces; larger complexes lay out their vertices and place the rest),
// normalise (into the unit ball, default true).
//
// Returns {vertices: [[x,y,z]], edges: [{points, tail, head, loop, arrow: {at, dir}}],
// faces: [{positions: Float32Array, indices: Uint32Array, centre, balloon, boundary}],
// transform: {centre, scale}, unit} with every array indexed as the complex's cells.
export function realise(complex, options = {}) {
  const { seed = 1, iterations = 200, rings = 3, smoothing = 24, full = 1500, normalise = true } = options;
  const random = seeded(seed);
  const V = complex.vertices, edges = complex.edges, faces = faceList(complex), E = edges.length, F = faces.length;
  const samples = options.samples ?? (E > 2000 ? 6 : E > 400 ? 10 : 18);
  const given = Array.isArray(complex.positions) && complex.positions.length === V && V > 0;
  let vpos, mid = null, centre = null;

  if (given) {
    vpos = complex.positions.map((p) => [Number(p[0]), Number(p[1]), Number(p[2] ?? 0)]);
  } else if (!isRegular(complex) && V + E + F <= full) {
    // Vertices, edge midpoints, face centres: one graph, laid out together.
    const adj = Array.from({ length: V + E + F }, () => new Set());
    const link = (a, b) => { if (a !== b) { adj[a].add(b); adj[b].add(a); } };
    edges.forEach(([t, h], e) => { link(V + e, t); link(V + e, h); });
    faces.forEach((f, i) => {
      const node = V + E + i;
      if (!f.boundary.length) link(node, f.base);
      for (const k of f.boundary) { link(node, V + letterEdge(k)); }
    });
    const P = layoutGraph(adj.map((s) => [...s]), iterations, random);
    vpos = P.slice(0, V);
    mid = P.slice(V, V + E);
    centre = P.slice(V + E);
  } else {
    // Regular or large: the 1-skeleton laid out, the rest placed by it.
    const adj = Array.from({ length: V }, () => new Set());
    edges.forEach(([t, h]) => { if (t !== h) { adj[t].add(h); adj[h].add(t); } });
    vpos = layoutGraph(adj.map((s) => [...s]), Math.min(iterations, 60), random);
  }

  // A typical length: the median distance between the ends of the edges (or 1).
  const lengths = edges.filter(([t, h]) => t !== h).map(([t, h]) => norm(sub(vpos[t], vpos[h]))).filter((l) => l > 1e-9).sort((a, b) => a - b);
  const L = lengths.length ? lengths[lengths.length >> 1] : 1;
  // Where each vertex's neighbours lie, to point loops and balloons away from them.
  const away = vpos.map(() => [0, 0, 0]);
  const cog = vpos.reduce((s, p) => add(s, mul(p, 1 / Math.max(1, V))), [0, 0, 0]);
  edges.forEach(([t, h]) => { if (t !== h) { away[t] = add(away[t], sub(vpos[t], vpos[h])); away[h] = add(away[h], sub(vpos[h], vpos[t])); } });
  const outward = (v, k = 0) => {
    let d = away[v];
    if (norm(d) < 1e-9) d = sub(vpos[v], cog);
    if (norm(d) < 1e-9) d = [0, 1, 0];
    d = unit(d);
    if (!k) return d;
    // The k-th of several: turned about a perpendicular by multiples of 50 degrees.
    const p = perpendicular(d), q = cross(d, p), a = 0.9 * Math.ceil(k / 2) * (k % 2 ? 1 : -1);
    return unit(add(mul(d, Math.cos(a)), add(mul(p, Math.sin(a) * Math.cos(k)), mul(q, Math.sin(a) * Math.sin(k)))));
  };

  if (!mid) {
    // Midpoints: parallel edges and loops fanned apart.
    const groups = new Map();
    edges.forEach(([t, h], e) => { const g = `${Math.min(t, h)},${Math.max(t, h)}`; (groups.get(g) || groups.set(g, []).get(g)).push(e); });
    mid = new Array(E);
    for (const list of groups.values()) list.forEach((e, k) => {
      const [t, h] = edges[e];
      if (t === h) { mid[e] = add(vpos[t], mul(outward(t, k), 0.6 * L)); return; }
      const m = mul(add(vpos[t], vpos[h]), 0.5);
      const off = (k - (list.length - 1) / 2) * 0.35 * norm(sub(vpos[h], vpos[t]));
      mid[e] = off ? add(m, mul(perpendicular(sub(vpos[h], vpos[t])), off)) : m;
    });
  }

  // The size of an edge as drawn (a loop's, the diameter of its circle), for arrows and vertices.
  const drawn = edges.map(([t, h], e) => (t === h ? norm(sub(mid[e], vpos[t])) : norm(sub(vpos[h], vpos[t])))).filter((l) => l > 1e-9).sort((a, b) => a - b);
  const typical = drawn.length ? drawn[drawn.length >> 1] : L;

  // The edge curves: quadratic through the midpoint, or a circle through vertex and midpoint.
  const edgeGeometry = edges.map(([t, h], e) => {
    const points = [];
    if (t === h) {
      const a = vpos[t], m = mid[e];
      let axis = sub(m, a);
      if (norm(axis) < 1e-9) axis = mul(outward(t), 0.6 * L);
      const c = mul(add(a, add(a, axis)), 0.5), r = norm(axis) / 2, u = unit(axis), w = perpendicular(axis);
      for (let s = 0; s <= samples; s++) {
        const th = (2 * Math.PI * s) / samples;
        points.push(add(c, add(mul(u, -r * Math.cos(th)), mul(w, r * Math.sin(th)))));
      }
      return { points, tail: t, head: h, loop: true, arrow: { at: add(c, mul(u, r)), dir: w } };
    }
    const a = vpos[t], b = vpos[h], ctrl = sub(mul(mid[e], 2), mul(add(a, b), 0.5));
    for (let s = 0; s <= samples; s++) {
      const x = s / samples;
      points.push(add(add(mul(a, (1 - x) * (1 - x)), mul(ctrl, 2 * x * (1 - x))), mul(b, x * x)));
    }
    const dir = unit(sub(b, a));
    return { points, tail: t, head: h, loop: false, arrow: { at: add(add(mul(a, 0.25), mul(ctrl, 0.5)), mul(b, 0.25)), dir: unit(sub(b, a), dir) } };
  });

  // The face discs.
  const faceGeometry = faces.map((f, i) => {
    if (!f.boundary.length) return balloon(vpos[f.base], outward(f.base, i + 1), 0.4 * L);
    const boundary = [];
    for (const k of f.boundary) {
      const pts = edgeGeometry[letterEdge(k)].points;
      const seq = k > 0 ? pts : [...pts].reverse();
      for (let s = 0; s < seq.length - 1; s++) boundary.push(seq[s]);
    }
    const c = centre ? centre[i] : boundary.reduce((s, p) => add(s, mul(p, 1 / boundary.length)), [0, 0, 0]);
    // A triangle's fan is already flat; longer boundaries are relaxed into a film.
    return disc(boundary, c, f.boundary.length <= 3 ? 1 : rings, smoothing);
  });

  // Into the unit ball.
  const all = [...vpos, ...edgeGeometry.flatMap((g) => g.points)];
  for (const g of faceGeometry) for (let j = 0; j < g.positions.length; j += 3) all.push([g.positions[j], g.positions[j + 1], g.positions[j + 2]]);
  let cen = [0, 0, 0], scale = 1;
  if (normalise && all.length) {
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (const p of all) for (let d = 0; d < 3; d++) { lo[d] = Math.min(lo[d], p[d]); hi[d] = Math.max(hi[d], p[d]); }
    cen = mul(add(lo, hi), 0.5);
    let R = 1e-9;
    for (const p of all) R = Math.max(R, norm(sub(p, cen)));
    scale = 1 / R;
  }
  const T = (p) => mul(sub(p, cen), scale);
  for (const g of faceGeometry) {
    for (let j = 0; j < g.positions.length; j += 3) {
      const q = T([g.positions[j], g.positions[j + 1], g.positions[j + 2]]);
      g.positions[j] = q[0]; g.positions[j + 1] = q[1]; g.positions[j + 2] = q[2];
    }
    g.centre = T(g.centre);
    g.boundary = g.boundary.map(T);
  }
  return {
    vertices: vpos.map(T),
    edges: edgeGeometry.map((g) => ({ ...g, points: g.points.map(T), arrow: { at: T(g.arrow.at), dir: g.arrow.dir } })),
    faces: faceGeometry,
    transform: { centre: cen, scale },
    unit: typical * scale,
  };
}

// Whether every cell is embedded on its boundary as a polygon is: no loops, no two edges with the
// same ends, no face that passes a vertex twice. Such a complex is drawn from its vertices alone.
export function isRegular(c) {
  const pairs = new Set();
  for (const [t, h] of c.edges) {
    if (t === h) return false;
    const k = t < h ? `${t},${h}` : `${h},${t}`;
    if (pairs.has(k)) return false;
    pairs.add(k);
  }
  for (const f of faceList(c)) {
    if (!f.boundary.length) return false;
    const seen = new Set();
    for (const k of f.boundary) {
      const [t, h] = c.edges[letterEdge(k)];
      const v = k > 0 ? t : h;
      if (seen.has(v)) return false;
      seen.add(v);
    }
  }
  return true;
}

// A disc spanned on a closed boundary: `rings` rings from the centre out, the interior relaxed
// by Jacobi smoothing toward the mean of its neighbours (a discrete harmonic map).
function disc(boundary, c, rings, smoothing) {
  const m = boundary.length, K = Math.max(1, rings);
  const P = [c];
  for (let k = 1; k <= K; k++) for (let i = 0; i < m; i++) P.push(add(c, mul(sub(boundary[i], c), k / K)));
  const at = (k, i) => (k === 0 ? 0 : 1 + (k - 1) * m + ((i % m) + m) % m);
  for (let it = 0; it < smoothing && K > 1; it++) {
    const Q = P.map((p) => p);
    let s = [0, 0, 0];
    for (let i = 0; i < m; i++) s = add(s, P[at(1, i)]);
    Q[0] = mul(s, 1 / m);
    for (let k = 1; k < K; k++) for (let i = 0; i < m; i++) {
      const inner = k === 1 ? P[0] : P[at(k - 1, i)];
      Q[at(k, i)] = mul(add(add(P[at(k, i - 1)], P[at(k, i + 1)]), add(inner, P[at(k + 1, i)])), 0.25);
    }
    for (let j = 0; j < P.length; j++) P[j] = Q[j];
  }
  const indices = [];
  for (let i = 0; i < m; i++) indices.push(0, at(1, i), at(1, i + 1));
  for (let k = 1; k < K; k++) for (let i = 0; i < m; i++) {
    const a = at(k, i), b = at(k, i + 1), d = at(k + 1, i), e = at(k + 1, i + 1);
    indices.push(a, d, e, a, e, b);
  }
  const positions = new Float32Array(P.length * 3);
  P.forEach((p, j) => { positions[j * 3] = p[0]; positions[j * 3 + 1] = p[1]; positions[j * 3 + 2] = p[2]; });
  return { positions, indices: new Uint32Array(indices), centre: P[0], balloon: false, boundary };
}

// A sphere through the point p, bulging toward d: the 2-cell attached to a point.
function balloon(p, d, r) {
  const c = add(p, mul(d, r)), u = perpendicular(d), w = cross(d, u);
  const nLat = 10, nLon = 16, P = [];
  for (let a = 0; a <= nLat; a++) {
    const th = (Math.PI * a) / nLat;
    for (let b = 0; b < nLon; b++) {
      const ph = (2 * Math.PI * b) / nLon;
      // th = π lands on p (the pole opposite d).
      P.push(add(c, add(mul(d, r * Math.cos(th)), add(mul(u, r * Math.sin(th) * Math.cos(ph)), mul(w, r * Math.sin(th) * Math.sin(ph))))));
    }
  }
  const indices = [];
  for (let a = 0; a < nLat; a++) for (let b = 0; b < nLon; b++) {
    const i = a * nLon + b, j = a * nLon + ((b + 1) % nLon), k = i + nLon, l = j + nLon;
    indices.push(i, k, l, i, l, j);
  }
  const positions = new Float32Array(P.length * 3);
  P.forEach((q, j) => { positions[j * 3] = q[0]; positions[j * 3 + 1] = q[1]; positions[j * 3 + 2] = q[2]; });
  return { positions, indices: new Uint32Array(indices), centre: c, balloon: true, boundary: [p] };
}

// ---- the net of the 2-cells --------------------------------------------------------------------

const sub2 = (a, b) => [a[0] - b[0], a[1] - b[1]];
const add2 = (a, b) => [a[0] + b[0], a[1] + b[1]];
const mul2 = (a, s) => [a[0] * s, a[1] * s];
const cross2 = (a, b) => a[0] * b[1] - a[1] * b[0];
const close2 = (a, b, eps = 1e-6) => Math.abs(a[0] - b[0]) < eps && Math.abs(a[1] - b[1]) < eps;

// The corners of a regular polygon of m unit sides, anticlockwise, side 0 along the bottom.
function regular(m) {
  const R = 1 / (2 * Math.sin(Math.PI / m));
  return Array.from({ length: m }, (_, i) => {
    const a = -Math.PI / 2 - Math.PI / m + (2 * Math.PI * i) / m;
    return [R * Math.cos(a), R * Math.sin(a)];
  });
}

// The convex outline a polygon occupies, for the test of overlap.
function outline(poly) {
  if (poly.kind === 'polygon') return poly.points;
  const c = poly.centre, r = poly.radius;
  return Array.from({ length: 8 }, (_, i) => add2(c, [r * Math.cos((Math.PI * i) / 4), r * Math.sin((Math.PI * i) / 4)]));
}

// A convex outline shrunk slightly toward its centroid, so that shared sides do not overlap.
function shrunk(P) {
  const c = P.reduce((s, p) => add2(s, mul2(p, 1 / P.length)), [0, 0]);
  return P.map((p) => add2(c, mul2(sub2(p, c), 0.995)));
}

// Whether two convex polygons, already shrunk, overlap (separating axes).
function overlap(a, b) {
  for (const P of [a, b]) for (let i = 0; i < P.length; i++) {
    const e = sub2(P[(i + 1) % P.length], P[i]), n = [-e[1], e[0]];
    let amin = Infinity, amax = -Infinity, bmin = Infinity, bmax = -Infinity;
    for (const p of a) { const v = p[0] * n[0] + p[1] * n[1]; amin = Math.min(amin, v); amax = Math.max(amax, v); }
    for (const p of b) { const v = p[0] * n[0] + p[1] * n[1]; bmin = Math.min(bmin, v); bmax = Math.max(bmax, v); }
    if (amax <= bmin + 1e-9 || bmax <= amin + 1e-9) return false;
  }
  return true;
}

// Each 2-cell as a polygon, unfolded into nets. Returns {polygons: [{face, piece, kind:
// 'polygon' | 'digon' | 'monogon' | 'disc', points, centre, radius, mirrored, sides: [{index, edge,
// sign, from, to, tail, head, glued, curve}]}], pieces, single, bbox}. A side's `curve` is a
// control point (for the sides of a digon) or a circle {centre, r} (the side of a monogon).
export function net(complex, options = {}) {
  const faces = faceList(complex), gap = options.gap ?? 0.8;
  const make = (f, i) => {
    const m = f.boundary.length;
    if (m === 0) return { face: i, kind: 'disc', points: [], centre: [0, 0], radius: 0.35, sides: [], mirrored: false };
    if (m === 1) {
      const r = 0.45, c = [0, r];
      return { face: i, kind: 'monogon', points: [[0, 0]], centre: c, radius: r, mirrored: false,
        sides: [{ index: 0, edge: letterEdge(f.boundary[0]), sign: Math.sign(f.boundary[0]), from: [0, 0], to: [0, 0], curve: { centre: c, r } }] };
    }
    const pts = m === 2 ? [[-0.5, 0], [0.5, 0]] : regular(m);
    const sides = f.boundary.map((k, s) => ({ index: s, edge: letterEdge(k), sign: Math.sign(k), from: pts[s], to: pts[(s + 1) % m] }));
    if (m === 2) {
      sides[0].curve = [0, -0.4];
      sides[1].curve = [0, 0.4];
    }
    const centre = pts.reduce((s, p) => add2(s, mul2(p, 1 / m)), [0, 0]);
    return { face: i, kind: m === 2 ? 'digon' : 'polygon', points: pts, centre, radius: m === 2 ? 0.5 : 1 / (2 * Math.sin(Math.PI / m)), mirrored: false, sides };
  };
  const polys = faces.map(make);
  // Where each edge occurs among the sides of straight polygons.
  const occurs = new Map();
  polys.forEach((p) => { if (p.kind === 'polygon') for (const s of p.sides) (occurs.get(s.edge) || occurs.set(s.edge, []).get(s.edge)).push([p.face, s.index]); });
  const placed = new Array(polys.length).fill(null);
  const pieces = [];
  // The placed outlines of the present piece, in buckets of a grid two units wide.
  let buckets = new Map();
  const shapes = new Map();
  const cellsOf = (P) => {
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (const p of P) { x0 = Math.min(x0, p[0]); y0 = Math.min(y0, p[1]); x1 = Math.max(x1, p[0]); y1 = Math.max(y1, p[1]); }
    const out = [];
    for (let i = Math.floor(x0 / 2); i <= Math.floor(x1 / 2); i++) for (let j = Math.floor(y0 / 2); j <= Math.floor(y1 / 2); j++) out.push(`${i},${j}`);
    return out;
  };
  const remember = (face, P) => { const sh = shrunk(P); shapes.set(face, sh); for (const k of cellsOf(sh)) (buckets.get(k) || buckets.set(k, []).get(k)).push(face); };
  const collides = (P) => {
    const sh = shrunk(P), seen = new Set();
    for (const k of cellsOf(sh)) for (const f of buckets.get(k) || []) {
      if (seen.has(f)) continue;
      seen.add(f);
      if (overlap(sh, shapes.get(f))) return true;
    }
    return false;
  };
  const transform = (poly, map) => {
    const T = (p) => map(p);
    const out = { ...poly, points: poly.points.map(T), centre: T(poly.centre), sides: poly.sides.map((s) => ({ ...s, from: T(s.from), to: T(s.to) })) };
    return out;
  };
  // The rigid motion taking a to A and b to B (|ab| = |AB|), then mirrored in AB if `mirror`.
  const motion = (a, b, A, B, mirror) => {
    const u = sub2(b, a), U = sub2(B, A);
    const ang = Math.atan2(U[1], U[0]) - Math.atan2(u[1], u[0]), c = Math.cos(ang), s = Math.sin(ang);
    const dir = mul2(U, 1 / Math.hypot(U[0], U[1]));
    return (p) => {
      const q = sub2(p, a);
      let r = add2(A, [c * q[0] - s * q[1], s * q[0] + c * q[1]]);
      if (mirror) { const v = sub2(r, A), t = v[0] * dir[0] + v[1] * dir[1]; r = sub2(add2(A, mul2(mul2(dir, t), 2)), v); }
      return r;
    };
  };
  const ends = (side) => (side.sign > 0 ? [side.from, side.to] : [side.to, side.from]);
  // The placed faces whose outlines share a bucket with q's.
  const near = (q) => { const out = new Set(); for (const k of cellsOf(outline(q))) for (const f of buckets.get(k) || []) out.add(f); return [...out]; };
  // Glue a side of q to a side of p where they lie on one another, oriented alike.
  const markCoincident = (q, piece) => {
    for (const s of q.sides) {
      if (s.glued) continue;
      for (const pi of piece) {
        const p = placed[pi];
        if (p === q || p.kind !== 'polygon') continue;
        for (const t of p.sides) {
          if (t.glued || t.edge !== s.edge || (p === q && t === s)) continue;
          const [a, b] = ends(s), [A, B] = ends(t);
          if (close2(a, A) && close2(b, B)) { s.glued = t.glued = true; s.mate = [p.face, t.index]; t.mate = [q.face, s.index]; }
        }
      }
    }
  };
  for (let root = 0; root < polys.length; root++) {
    if (placed[root]) continue;
    const piece = [root];
    placed[root] = { ...polys[root], piece: pieces.length, sides: polys[root].sides.map((s) => ({ ...s, glued: false })) };
    buckets = new Map();
    remember(root, outline(placed[root]));
    if (polys[root].kind === 'polygon') {
      for (let h = 0; h < piece.length; h++) {
        const p = placed[piece[h]];
        for (const side of p.sides) {
          if (side.glued) continue;
          for (const [qi, j] of occurs.get(side.edge) || []) {
            if (placed[qi] || side.glued) continue;
            const q0 = polys[qi], t = q0.sides[j];
            const [a, b] = ends(t), [A, B] = ends(side);
            let q = transform(q0, motion(a, b, A, B, false));
            const line = sub2(B, A), sideOf = (pt) => cross2(line, sub2(pt, A));
            if (Math.sign(sideOf(q.centre)) === Math.sign(sideOf(p.centre))) { q = transform(q0, motion(a, b, A, B, true)); q.mirrored = true; }
            const out = outline(q);
            if (collides(out)) continue;
            q.piece = pieces.length;
            q.sides = q.sides.map((s) => ({ ...s, glued: false }));
            q.sides[j].glued = true; q.sides[j].mate = [p.face, side.index];
            side.glued = true; side.mate = [qi, j];
            placed[qi] = q;
            piece.push(qi);
            remember(qi, out);
            markCoincident(q, near(q));
          }
        }
      }
    }
    pieces.push(piece);
  }
  // The sides' tails and heads, and the pieces set side by side in rows.
  const boxes = pieces.map((piece) => {
    const lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
    for (const i of piece) for (const p of outline(placed[i])) { lo[0] = Math.min(lo[0], p[0]); lo[1] = Math.min(lo[1], p[1]); hi[0] = Math.max(hi[0], p[0]); hi[1] = Math.max(hi[1], p[1]); }
    return { lo, hi, w: hi[0] - lo[0], h: hi[1] - lo[1] };
  });
  const area = boxes.reduce((s, b) => s + (b.w + gap) * (b.h + gap), 0);
  const rowWidth = Math.max(Math.sqrt(area) * 1.6, ...boxes.map((b) => b.w));
  let x = 0, y = 0, rowH = 0;
  const bbox = { x0: Infinity, y0: Infinity, x1: -Infinity, y1: -Infinity };
  pieces.forEach((piece, k) => {
    const b = boxes[k];
    if (x > 0 && x + b.w > rowWidth) { x = 0; y += rowH + gap; rowH = 0; }
    const shift = [x - b.lo[0], y - b.lo[1]];
    for (const i of piece) {
      const p = placed[i];
      const T = (q) => add2(q, shift);
      p.points = p.points.map(T); p.centre = T(p.centre);
      for (const s of p.sides) {
        s.from = T(s.from); s.to = T(s.to);
        if (Array.isArray(s.curve)) s.curve = T(s.curve);
        else if (s.curve) s.curve = { centre: T(s.curve.centre), r: s.curve.r };
        [s.tail, s.head] = ends(s);
      }
      for (const q of outline(p)) { bbox.x0 = Math.min(bbox.x0, q[0]); bbox.y0 = Math.min(bbox.y0, q[1]); bbox.x1 = Math.max(bbox.x1, q[0]); bbox.y1 = Math.max(bbox.y1, q[1]); }
    }
    x += b.w + gap;
    rowH = Math.max(rowH, b.h);
  });
  // Digon sides bulge away from their centre: their control points were set before placement.
  for (const p of placed) if (p && p.kind === 'digon') {
    const [a, b] = [p.sides[0].from, p.sides[0].to];
    const mid = mul2(add2(a, b), 0.5), n = [-(b[1] - a[1]), b[0] - a[0]];
    p.sides[0].curve = add2(mid, mul2(n, -0.4));
    p.sides[1].curve = add2(mid, mul2(n, 0.4));
  }
  if (!polys.length) Object.assign(bbox, { x0: 0, y0: 0, x1: 1, y1: 1 });
  return { polygons: placed, pieces: pieces.length, single: polys.length > 0 && pieces.length === 1, bbox };
}

// The boundary of a piece of a net, read around: its unglued sides in order, each cycle of
// them a list of {face, index, edge, sign, from, to}, where sign is the letter as the cycle
// reads it (a side drawn in a mirrored polygon is read against its face's word). For a
// one-face surface, or a surface unfolded into one disc, the single cycle is its polygon schema.
export function outerSides(netResult, piece = 0) {
  const sides = [];
  for (const p of netResult.polygons) if (p && p.piece === piece) for (const s of p.sides) if (!s.glued) {
    const forward = !p.mirrored;
    sides.push({ face: p.face, index: s.index, edge: s.edge, sign: forward ? s.sign : -s.sign, from: forward ? s.from : s.to, to: forward ? s.to : s.from });
  }
  const key = (q) => `${q[0].toFixed(5)},${q[1].toFixed(5)}`;
  const starting = new Map();
  sides.forEach((s, i) => { (starting.get(key(s.from)) || starting.set(key(s.from), []).get(key(s.from))).push(i); });
  const used = new Uint8Array(sides.length), cycles = [];
  for (let i = 0; i < sides.length; i++) {
    if (used[i]) continue;
    const cycle = [];
    let j = i;
    while (j !== undefined && !used[j]) {
      used[j] = 1;
      cycle.push(sides[j]);
      j = (starting.get(key(sides[j].to)) || []).find((x) => !used[x]);
    }
    cycles.push(cycle);
  }
  return cycles;
}

// ---- the Hasse diagram ------------------------------------------------------------------------------

// The cells as nodes in layers by dimension (0 at the bottom), each joined to the cells of its
// boundary, with the multiplicity of the incidence; ordered by barycentres. Returns {layers:
// [[cell index, ...] by dimension, in order], nodes: [[{x, y}] by dimension], links: [{dim, cell,
// face, coefficient}] (dim the higher cell's dimension)}; x in [0, 1], y the dimension.
export function hasse(complex, options = {}) {
  const sweeps = options.sweeps ?? 12;
  const counts = cellCounts(complex);
  const links = [];
  complex.edges.forEach(([t, h], e) => {
    if (t === h) links.push({ dim: 1, cell: e, face: t, coefficient: 0 });
    else { links.push({ dim: 1, cell: e, face: t, coefficient: -1 }); links.push({ dim: 1, cell: e, face: h, coefficient: 1 }); }
  });
  faceList(complex).forEach((f, i) => {
    const m = new Map();
    for (const k of f.boundary) m.set(letterEdge(k), (m.get(letterEdge(k)) || 0) + Math.sign(k));
    if (!f.boundary.length) links.push({ dim: 2, cell: i, face: f.base, coefficient: 0, point: true });
    for (const [e, coef] of m) links.push({ dim: 2, cell: i, face: e, coefficient: coef });
  });
  (complex.cells || []).forEach((layer, j) => (layer || []).forEach((cell, i) => {
    const m = new Map();
    for (const [idx, coef] of cell.boundary || []) m.set(idx, (m.get(idx) || 0) + coef);
    for (const [idx, coef] of m) links.push({ dim: j + 3, cell: i, face: idx, coefficient: coef });
  }));
  const down = counts.map((n) => Array.from({ length: n }, () => []));
  const up = counts.map((n) => Array.from({ length: n }, () => []));
  for (const l of links) if (l.dim < counts.length && !l.point) { down[l.dim][l.cell].push(l.face); up[l.dim - 1][l.face].push(l.cell); }
  const order = counts.map((n) => Array.from({ length: n }, (_, i) => i));
  const pos = counts.map((n) => new Float64Array(n));
  const place = (d) => order[d].forEach((c, r) => { pos[d][c] = (r + 0.5) / order[d].length; });
  counts.forEach((_, d) => place(d));
  const sortBy = (d, nbr, other) => {
    const key = new Float64Array(counts[d]);
    for (let c = 0; c < counts[d]; c++) key[c] = nbr[d][c].length ? nbr[d][c].reduce((s, x) => s + pos[other][x], 0) / nbr[d][c].length : pos[d][c];
    order[d].sort((a, b) => key[a] - key[b] || a - b);
    place(d);
  };
  for (let s = 0; s < sweeps; s++) {
    for (let d = 1; d < counts.length; d++) sortBy(d, down, d - 1);
    for (let d = counts.length - 2; d >= 0; d--) sortBy(d, up, d + 1);
  }
  return { layers: order, nodes: pos.map((xs, d) => Array.from(xs, (x) => ({ x, y: d }))), links };
}

// ---- colour maps in the warm family ---------------------------------------------------------------

// sRGB <-> OKLab (Ottosson 2020), so that maps are interpolated where distances are perceptual.
const toLinear = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
const toGamma = (c) => (c <= 0.0031308 ? 12.92 * c : 1.055 * c ** (1 / 2.4) - 0.055);
function oklab([r, g, b]) {
  const R = toLinear(r), G = toLinear(g), B = toLinear(b);
  const l = Math.cbrt(0.4122214708 * R + 0.5363325363 * G + 0.0514459929 * B);
  const m = Math.cbrt(0.2119034982 * R + 0.6806995451 * G + 0.1073969566 * B);
  const s = Math.cbrt(0.0883024619 * R + 0.2817188376 * G + 0.6299787005 * B);
  return [0.2104542553 * l + 0.793617785 * m - 0.0040720468 * s, 1.9779984951 * l - 2.428592205 * m + 0.4505937099 * s, 0.0259040371 * l + 0.7827717662 * m - 0.808675766 * s];
}
function rgbOf([L, a, b]) {
  const l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3, m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3, s = (L - 0.0894841775 * a - 1.291485548 * b) ** 3;
  const clamp = (x) => Math.min(1, Math.max(0, toGamma(x)));
  return [clamp(4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s), clamp(-1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s), clamp(-0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s)];
}
const hexRgb = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16) / 255);
export const cssColour = (rgb) => `#${rgb.map((c) => Math.round(Math.min(1, Math.max(0, c)) * 255).toString(16).padStart(2, '0')).join('')}`;

function ramp(stops, cyclic = false) {
  const labs = stops.map((h) => oklab(hexRgb(h)));
  const n = cyclic ? labs.length : labs.length - 1;
  return (t) => {
    t = cyclic ? ((t % 1) + 1) % 1 : Math.min(1, Math.max(0, t));
    const x = t * n, i = Math.min(Math.floor(x), n - 1), f = x - i;
    const a = labs[i], b = labs[(i + 1) % labs.length];
    return rgbOf([a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f]);
  };
}

// For angles: plum, Old Rose, gold, blush lilac, and round again; θ in radians, any period of 2π.
const cyclicRamp = ramp(['#6b4a86', '#c08081', '#e9b872', '#d9b8d6'], true);
export const cyclicColour = (theta) => cyclicRamp(theta / (2 * Math.PI));
// For magnitudes in [0, 1]: plum night to Old Rose to cream.
export const sequentialColour = ramp(['#2a1f3a', '#6b4a86', '#c08081', '#f0d9b5']);
// For signed values in [-1, 1]: lavender through a pale rose to Old Rose and deep rust.
const divergingRamp = ramp(['#5d5aa0', '#a99ad0', '#f2e6ee', '#d38e8a', '#9c4a45']);
export const divergingColour = (t) => divergingRamp((t + 1) / 2);
// For telling cells apart: ten warm, distinct hues.
const CATEGORICAL = ['#e0a3b6', '#e9b872', '#9fc2b8', '#b9a2e0', '#e48f78', '#c9d38a', '#8fb3dd', '#d987c6', '#f0cf9a', '#a6d6c4'];
export const categoricalColour = (i) => hexRgb(CATEGORICAL[((i % CATEGORICAL.length) + CATEGORICAL.length) % CATEGORICAL.length]);

// ---- coloured graphs and quantum simplices ----------------------------------------------------------
//
// A coloured graph (the Feynman graphs and bubbles of tensor models and group field theories,
// morphe's ColouredGraph) is {colours, p, sigma}: p black vertices (the fields φ) and p white
// ones (φ̄), and for each colour c a permutation sigma[c] joining black b to white sigma[c][b].
// Vertices are numbered black 0..p-1, white p..2p-1; the edge of colour c at black b is the
// edge c p + b.

// The colour of index c of a coloured graph: colour 0, the propagators of a Feynman graph, a pale
// hairline cream; the others the warm categorical hues.
export const indexColour = (c, { propagator = true } = {}) => (propagator && c === 0 ? hexRgb('#efe4d8') : categoricalColour(propagator ? c - 1 : c));

// The 2-complex of a coloured graph: its vertices and edges, and a 2-cell bounded by each
// bicoloured cycle, the faces of colours i < j as the cycles of sigma_j⁻¹ sigma_i, each from
// its first black vertex along colour i and back along colour j. The numbering is morphe's
// ColouredGraph::two_complex(), so that the engine's answers and this agree cell by cell.
export function colouredGraphComplex(cg) {
  const { colours, p, sigma } = cg;
  const edges = [], names = [];
  for (let c = 0; c < colours; c++) for (let b = 0; b < p; b++) { edges.push([b, p + sigma[c][b]]); names.push(`c${c}.${b + 1}`); }
  const faces = [];
  for (let i = 0; i < colours; i++) for (let j = i + 1; j < colours; j++) {
    const inv = new Array(p);
    for (let b = 0; b < p; b++) inv[sigma[j][b]] = b;
    const seen = new Uint8Array(p);
    for (let b0 = 0; b0 < p; b0++) {
      if (seen[b0]) continue;
      const cyc = [];
      for (let b = b0; !seen[b]; b = inv[sigma[i][b]]) { seen[b] = 1; cyc.push(b); }
      const boundary = [];
      cyc.forEach((b, k) => { boundary.push(i * p + b + 1, -(j * p + cyc[(k + 1) % cyc.length] + 1)); });
      faces.push({ base: cyc[0], boundary });
    }
  }
  return { name: 'coloured graph', title: `2-complex of a ${colours}-coloured graph`, vertices: 2 * p, edges, faces, cells: [], labels: { edges: names }, coloured_graph: { colours, p, sigma } };
}

// A drawing of a small coloured graph in the plane: the vertices laid out by stress
// majorisation of their graph distances (as realise() does, then seen along the two principal
// axes), the parallel edges of a pair bent apart. Returns {nodes: [{x, y, black}], edges:
// [{u, v, colour, bend}]} with x and y in [0, 1].
export function colouredGraphLayout(cg, { seed = 1 } = {}) {
  const { colours, p, sigma } = cg, n = 2 * p;
  const adj = Array.from({ length: n }, () => []);
  const edges = [];
  for (let c = 0; c < colours; c++) for (let b = 0; b < p; b++) {
    const w = p + sigma[c][b];
    edges.push({ u: b, v: w, colour: c, bend: 0 });
    if (!adj[b].includes(w)) { adj[b].push(w); adj[w].push(b); }
  }
  let nodes;
  if (n === 2) nodes = [[0, 0], [1, 0]];
  else {
    const X = layoutGraph(adj, 300, seeded(seed));
    // The two principal axes of the points.
    const m = [0, 0, 0];
    for (const x of X) for (let k = 0; k < 3; k++) m[k] += x[k] / n;
    const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (const x of X) for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) C[a][b] += (x[a] - m[a]) * (x[b] - m[b]);
    const axis = (skip) => {
      let v = skip ? [skip[1], -skip[0], 0.3] : [1, 0.4, 0.2];
      for (let it = 0; it < 60; it++) {
        let w = [0, 1, 2].map((a) => C[a][0] * v[0] + C[a][1] * v[1] + C[a][2] * v[2]);
        if (skip) { const d = w[0] * skip[0] + w[1] * skip[1] + w[2] * skip[2]; w = w.map((x, k) => x - d * skip[k]); }
        const len = Math.hypot(...w) || 1;
        v = w.map((x) => x / len);
      }
      return v;
    };
    const a1 = axis(null), a2 = axis(a1);
    nodes = X.map((x) => [0, 1].map((k) => { const a = k ? a2 : a1; return (x[0] - m[0]) * a[0] + (x[1] - m[1]) * a[1] + (x[2] - m[2]) * a[2]; }));
    // Stress majorisation in the plane itself, from the projection, so that no two vertices
    // that the projection brought together stay together (pairs in different components are
    // left to the projection).
    const all = Array.from({ length: n }, (_, i) => i), D = distances(adj, all);
    for (let it = 0; it < 300; it++) {
      for (let i = 0; i < n; i++) {
        let sx = 0, sy = 0, den = 0;
        for (let j = 0; j < n; j++) {
          const d = D[i * n + j];
          if (j === i || d === 65535) continue;
          const w = 1 / (d * d);
          let dx = nodes[i][0] - nodes[j][0], dy = nodes[i][1] - nodes[j][1];
          let len = Math.hypot(dx, dy);
          if (len < 1e-9) { dx = Math.cos(i + j); dy = Math.sin(i + j); len = 1; }
          sx += w * (nodes[j][0] + (d / len) * dx); sy += w * (nodes[j][1] + (d / len) * dy); den += w;
        }
        if (den) nodes[i] = [sx / den, sy / den];
      }
    }
  }
  const xs = nodes.map((q) => q[0]), ys = nodes.map((q) => q[1]);
  const x0 = Math.min(...xs), y0 = Math.min(...ys), span = Math.max(Math.max(...xs) - x0, Math.max(...ys) - y0) || 1;
  const ox = (1 - (Math.max(...xs) - x0) / span) / 2, oy = (1 - (Math.max(...ys) - y0) / span) / 2;
  // Bends: the k-th of m edges between one pair spread symmetrically about the straight line.
  const groups = new Map();
  for (const e of edges) { const key = `${e.u}:${e.v}`; (groups.get(key) || groups.set(key, []).get(key)).push(e); }
  for (const list of groups.values()) list.forEach((e, k) => { e.bend = list.length > 1 ? (k - (list.length - 1) / 2) / Math.max(1, list.length - 1) : 0; });
  return { nodes: nodes.map((q, i) => ({ x: ox + (q[0] - x0) / span, y: oy + (q[1] - y0) / span, black: i < p })), edges };
}

// A tetrahedron with edge lengths [l01, l02, l03, l12, l13, l23], set in space: vertex 0 at the
// origin, 1 on the x axis, 2 in the xy plane, 3 above it, then centred. When the lengths admit
// no Euclidean tetrahedron (the Cayley-Menger determinant is not positive) the last vertex is
// left in the plane and `euclidean` is false.
export function embedTetrahedron(lengths) {
  const [a, b, c, d, e, f] = lengths;
  const p0 = [0, 0, 0], p1 = [a, 0, 0];
  const x2 = (a * a + b * b - d * d) / (2 * a), y2 = Math.sqrt(Math.max(0, b * b - x2 * x2));
  const x3 = (a * a + c * c - e * e) / (2 * a);
  const y3 = y2 > 1e-12 ? (b * b + c * c - f * f - 2 * x2 * x3) / (2 * y2) : 0;
  const z2 = c * c - x3 * x3 - y3 * y3;
  const euclidean = y2 > 1e-12 && z2 > 1e-12;
  const pts = [p0, p1, [x2, y2, 0], [x3, y3, Math.sqrt(Math.max(0, z2))]];
  const m = [0, 1, 2].map((k) => pts.reduce((s, q) => s + q[k], 0) / 4);
  return { positions: pts.map((q) => q.map((x, k) => x - m[k])), euclidean, volume: euclidean ? (a * y2 * Math.sqrt(z2)) / 6 : 0 };
}
