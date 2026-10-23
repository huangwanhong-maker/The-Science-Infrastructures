// What a selection marks: a chain of one dimension (a homology generator), a loop (a generator
// of the fundamental group, a Wilson loop), or one cell; read by every view.

// The cells a selection marks, by dimension: maps from index to coefficient (a loop marks its
// edges with the sign of their first traversal, and its vertices).
export function selectionSets(selection, complex) {
  const out = { vertices: new Map(), edges: new Map(), faces: new Map() };
  if (!selection) return out;
  if (selection.kind === 'chain') {
    const target = selection.dim === 0 ? out.vertices : selection.dim === 1 ? out.edges : selection.dim === 2 ? out.faces : null;
    if (target) for (const [c, k] of selection.chain) target.set(c, (target.get(c) || 0) + k);
    if (selection.dim === 1) for (const [c] of selection.chain) for (const v of complex.edges[c] || []) out.vertices.set(v, 1);
  } else if (selection.kind === 'loop' || selection.kind === 'word') {
    for (const k of selection.word) {
      const e = Math.abs(k) - 1;
      if (!out.edges.has(e)) out.edges.set(e, Math.sign(k));
      for (const v of complex.edges[e] || []) out.vertices.set(v, 1);
    }
  } else if (selection.kind === 'cell') {
    const target = [out.vertices, out.edges, out.faces][selection.dim];
    target?.set(selection.index, 1);
  }
  return out;
}
