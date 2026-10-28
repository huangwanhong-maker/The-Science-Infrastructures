// Formalism I: a grammar's rules (L ⇐ K ⇒ R, relabelling allowed) and its derivation space up to
// isomorphism, by depth, with the number of matches on each arc and of derivations reaching each
// class; the symmetry and accounting facts of each example.

const EXAMPLES = [
  ['stereocentre', 'Stereocentre (Case Study I)'],
  ['ethanol', 'Ethanol (Fig. grammar-syntax)'],
  ['glycine', 'Glycine backbones (Figs. grammar-growth, -tree)'],
  ['proton', 'Acid–base edit (Fig. proton-graph-edit)'],
];

export async function render(root, ctx) {
  const { h, busy, failure } = ctx.dom;
  const { draw } = ctx.graph;
  const stored = () => { try { return sessionStorage.getItem('grammar-example'); } catch { return null; } };
  let example = (ctx.args && ctx.args[0]) || stored() || 'stereocentre';
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Formalism I: grammars of graphs'),
    h('span', { class: 'lede' }, 'Rules are double pushouts on labelled graphs (L ⇐ K ⇒ R) with relabelling, so that an atom may change its charge and keep its identity. The derivation space is computed up to isomorphism, with every match counted.')));
  const keys = h('div', { class: 'keys' });
  const body = h('div', {});
  root.append(keys, body);
  const ruleView = (q) => {
    const kept = new Set(q.keptVertices.map(([l]) => l)), keptR = new Set(q.keptVertices.map(([, r]) => r));
    return h('div', { class: 'panel', style: 'margin: 8px 0' },
      h('h3', {}, q.name),
      h('div', { class: 'rule' },
        h('div', { class: 'graph-card' }, draw(q.left, { width: 300, height: 190, highlight: kept, ids: true }), h('div', { class: 'cap' }, 'L (kept in violet)')),
        h('span', { class: 'arrow' }, '⇒'),
        h('div', { class: 'graph-card' }, draw(q.right, { width: 300, height: 190, highlight: keptR, ids: true }), h('div', { class: 'cap' }, 'R'))),
      h('p', { class: 'caption' }, `|Aut L| = ${q.automorphismsOfL}${q.conditional ? '; with an application condition' : ''}. K: ${q.keptVertices.map(([l, r]) => `${l}↦${r}`).join(', ')}.`));
  };
  const spaceView = (space, title) => {
    const byDepth = new Map();
    space.classes.forEach((c, i) => { if (!byDepth.has(c.depth)) byDepth.set(c.depth, []); byDepth.get(c.depth).push([c, i]); });
    const rows = [...byDepth.entries()].sort((a, b) => a[0] - b[0]).map(([d, cs]) => h('div', {},
      h('h3', { class: 'muted', style: 'margin-top:10px' }, `after ${d} step${d === 1 ? '' : 's'}`),
      h('div', { class: 'derivation' }, ...cs.slice(0, 24).map(([c, i]) => {
        const incoming = space.arcs.filter((a) => a.to === i);
        const n = c.graph.vertices.length;
        return h('div', { class: 'graph-card' }, draw(c.graph, { width: Math.min(440, Math.max(200, 110 + 30 * n)), height: n > 8 ? 170 : 140, seed: 3 }),
          h('div', { class: 'cap' }, `class ${i}${c.terminal ? ' · terminal' : ''} · ${c.derivations} derivation${c.derivations === 1 ? '' : 's'}`),
          incoming.length ? h('div', { class: 'cap faint' }, incoming.map((a) => `from ${a.from} × ${a.multiplicity}`).join(', ')) : null);
      }), cs.length > 24 ? h('p', { class: 'caption' }, `and ${cs.length - 24} more`) : null)));
    return h('div', { class: 'section' }, h('h2', {}, title), ...rows);
  };
  const show = async () => {
    for (const k of keys.children) k.classList.toggle('on', k.dataset.id === example);
    busy(body);
    let g;
    try {
      g = await ctx.ask('grammar', { example });
    } catch (error) {
      return failure(body, error);
    }
    body.replaceChildren(h('h2', {}, g.title), ...g.rules.map(ruleView));
    if (example === 'stereocentre') {
      body.append(h('dl', { class: 'facts panel' },
        h('dt', {}, 'terminal classes'), h('dd', {}, `${g.terminalClasses}: the carbon labelled τ = +1 or τ = −1`),
        h('dt', {}, 'orbits of the reflection p'), h('dd', {}, `${g.orbits} (p exchanges the two)`),
        h('dt', {}, 'constitutional classes'), h('dd', {}, `${g.constitutionalClasses} (forgetting τ, they coincide)`),
        h('dt', {}, 'rules closed under p'), h('dd', {}, g.closed ? 'yes; derivations commute with p' : 'no'),
        h('dt', {}, 'with p₊ alone'), h('dd', {}, `closed: ${g.halfClosed ? 'yes' : 'no'}; equivariant: ${g.halfEquivariant ? 'yes' : 'no'} (Prop. 4.1 needs closure)`)));
      body.append(spaceView(g.space, 'Derivation space'));
    } else if (example === 'glycine') {
      body.append(h('div', { class: 'grid cols-2' },
        h('dl', { class: 'facts panel' }, h('dt', {}, 'labelled units'), h('dd', {}, `${g.labelled.terminalClasses} backbones, ${g.labelled.terminalDerivations} derivations, ${g.labelled.orbits} orbit of S₃`),
          h('dt', {}, 'equivariant under S₃'), h('dd', {}, g.labelled.equivariant ? 'yes, at every derivable graph' : 'no'),
          h('dt', {}, 'waters created'), h('dd', {}, String(g.labelled.water))),
        h('dl', { class: 'facts panel' }, h('dt', {}, 'identical units'), h('dd', {}, `${g.identical.terminalClasses} backbone, ${g.identical.terminalDerivations} derivations`),
          h('dt', {}, '3 Gly → Gly₃ + 2 H₂O'), h('dd', {}, g.balanced ? 'balanced in C, H, N, O' : 'unbalanced'))));
      body.append(h('p', { class: 'note' }, 'A peptide bond is directed (carboxyl to amine), so the six labelled backbones are the six orders of the units; each is reached by two derivations, as the manuscript\'s tree shows for one of them.'));
      body.append(spaceView(g.labelled.space, 'Derivation space, labelled units'));
    } else if (example === 'ethanol') {
      body.append(spaceView(g.space, 'Derivation'));
    } else if (example === 'proton') {
      const rows = Object.keys(g.before).map((k) => h('tr', {}, h('td', {}, k), h('td', {}, String(g.before[k])), h('td', {}, String(g.after[k]))));
      body.append(h('div', { class: 'grid cols-2 section' },
        h('div', {}, h('h3', {}, 'host'), draw(g.host, { width: 360, height: 220, ids: true, seed: 5 })),
        h('div', {}, h('h3', {}, 'derived graph (atoms traced)'), g.result ? draw(g.result, { width: 360, height: 220, ids: true, seed: 5 }) : 'no derivation')),
        h('div', { class: 'grid cols-2' },
          h('div', { class: 'panel' }, h('h3', {}, 'ledger'), h('table', { class: 'data' }, h('tr', {}, h('th', {}, ''), h('th', {}, 'before'), h('th', {}, 'after')), ...rows),
            h('p', { class: 'caption' }, `Species balance AN = 0, zᵀN = 0: ${g.balanced ? 'holds' : 'fails'}.`)),
          h('div', { class: 'panel' }, h('h3', {}, 'counting'),
            h('p', {}, `${g.matches} matches, |Aut L| = ${g.automorphismsOfL}, ${g.productClasses} product class: the two water hydrogens are exchanged by an automorphism of L, so one occurrence gives two matches.`),
            h('p', { class: 'caption' }, 'The two oxygens keep their identity and change their charges: a relabelling of kept vertices, which a span of label-preserving monomorphisms (Def. 4.2 as written) cannot express.'))));
    }
  };
  for (const [id, name] of EXAMPLES) {
    const b = h('button', { class: 'key', 'data-id': id }, name);
    b.addEventListener('click', () => { example = id; try { sessionStorage.setItem('grammar-example', id); } catch { /* private window */ } show(); });
    keys.append(b);
  }
  show();
}
