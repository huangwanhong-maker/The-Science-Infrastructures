// The other complexes of a record (the theory note, §7): the refined complex K̂ and its foams, the
// history complex with its closed shells filled (K³), the Dowker complex of who met whom with
// Atkin's Q-analysis, and the history complex filtered by time, its persistent homology drawn as
// bars: when each circle of encounters closes, and when a later meeting fills it.

import { el, svg } from './dom.js';
import { listText, plural, roundYears, yearName } from './format.js';

const SUP = { 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const groupText = (text) => String(text).replace(/\^(\d+)/g, (_, n) => [...n].map((d) => SUP[d]).join(''));
const homologyText = (h) => h.map((g) => `H${String(g.dimension).replace(/\d/g, (d) => '₀₁₂₃₄₅₆₇₈₉'[d])} = ${groupText(g.text)}`).join(', ');

export function renderModels(data, act) {
  const m = data.models;
  if (!m) return null;
  const inv = data.invariants;
  // A combining circumflex drifts in serif bold, so the hat is placed over the K by the style.
  const refined = el('span', {}, el('span', { class: 'khat', text: 'K' }), ', the refined complex');
  const rows = [
    ['K, the history complex', 'occasions · lives between them · lenses of ties', homologyText(inv.homology), 'H₁: circles of encounters no tie fills; H₂: closed shells of ties'],
    [refined, `presences · lives and rungs · lenses and ${plural(m.refined.triangles, 'triangle')}`, homologyText(m.refined.homology),
      `the carrier of the signed record; its ${m.refined.foams} independent ℤ/2 2-cycles are the spin foams`],
    ['K³, the shells filled', `K and ${plural(m.filled.counts[3] || 0, '3-cell')}, one on each minimal closed shell`, homologyText(m.filled.homology),
      'H₂: shells that are not filled episodes; H₃: quartets of mutually tied beings whose four shells close into a 3-sphere'],
    ['The Dowker complex', `a simplex for those present together: ${m.dowker.simplices.map((n, d) => `${n} of dimension ${d}`).filter((_, d) => m.dowker.simplices[d]).slice(0, 4).join(', ')}…`, homologyText(m.dowker.homology.filter((h) => h.rank)),
      `who met whom, filled by any co-presence; Atkin's Q-vector, the q-connected components for q = 0, 1, …: (${m.dowker.q.join(', ')})`],
  ];
  return el('section', { class: 'models' },
    el('h3', { text: 'The complexes of a record' }),
    el('p', { class: 'note', text: 'One record makes several complexes, each answering its own question (the theory note, §7). The history complex reads the topology of encounters; its refinement, which opens each occasion into its presences, carries the signs and their spin foams; with its closed shells filled, collective episodes become 3-cells; Dowker\'s complex of the presence relation, equivalent to its dual on the occasions, is the topology of co-presence alone.' }),
    el('div', { class: 'table-wrap' }, el('table', { class: 'record models-table' },
      el('thead', {}, el('tr', {}, el('th', { text: 'Complex' }), el('th', { text: 'Cells' }), el('th', { text: 'Homology' }), el('th', { text: 'What it says' }))),
      el('tbody', {}, ...rows.map(([name, cells, h, about]) => el('tr', {}, el('td', { class: 'strong' }, name), el('td', { class: 'faint', text: cells }), el('td', { class: 'group', text: h }), el('td', { text: about })))))),
    barcode(data, act));
}

// Persistence in time: one bar per class, from the occasion that creates it to the one that kills it.
function barcode(data, act) {
  const iv = (data.models.time?.intervals || []).filter((i) => i.dimension >= 1);
  if (!iv.length) return el('p', { class: 'note', text: 'The history in time has no circle of encounters and no shell.' });
  const occ = data.occasions, order = data.complex.vertices;
  const y0 = occ[order[0]].year, y1 = occ[order[order.length - 1]].year;
  const H1 = iv.filter((i) => i.dimension === 1), H2 = iv.filter((i) => i.dimension === 2);
  const W = 900, L = 150, R = 30, rowH = 9, gap = 22;
  const X = (y) => L + ((y - y0) / Math.max(1e-9, y1 - y0)) * (W - L - R);
  const height = 40 + (H1.length + H2.length) * rowH + 2 * gap;
  const s = svg('svg', { viewBox: `0 0 ${W} ${height}`, class: 'barcode-svg', role: 'img', 'aria-label': 'Persistent homology of the history in time' });
  const span = y1 - y0, step = span > 400 ? 100 : span > 150 ? 25 : span > 60 ? 10 : 5;
  for (const y of roundYears(y0, y1, step)) {
    s.append(svg('line', { x1: X(y), x2: X(y), y1: 16, y2: height - 8, class: 'gridline' }), svg('text', { x: X(y), y: 11, 'text-anchor': 'middle', class: 'tick' }, yearName(y)));
  }
  let row = 26;
  const group = (name, list, cls) => {
    s.append(svg('text', { x: L - 10, y: row + 7, 'text-anchor': 'end', class: 'bar-group' }, name));
    for (const i of list) {
      const b = occ[i.birth].year, d = i.death === null ? y1 : occ[i.death].year;
      const bar = svg('rect', { x: X(b), y: row, width: Math.max(2, X(d) - X(b)), height: rowH - 3, rx: 1.5, class: `bar ${cls}${i.death === null ? ' essential' : ''}` });
      const who = i.creator.kind === 'segment' ? `${data.beings[i.creator.being].name} reaches “${occ[i.creator.to].label}”` : i.creator.kind === 'lens' ? `the lens of ${data.ties[i.creator.tie].id} closes at “${occ[i.creator.to].label}”` : '';
      const end = i.killer ? (i.killer.kind === 'lens' ? `filled when ${data.ties[i.killer.tie].id} meets again at “${occ[i.killer.to].label}”` : 'killed') : 'never filled';
      bar.append(svg('title', {}, `${i.dimension === 1 ? 'A circle of encounters' : 'A closed shell'}, born ${occ[i.birth].date}: ${who}; ${end}.`));
      if (i.creator.kind === 'segment') bar.addEventListener('click', () => act.pick({ kind: 'being', id: data.beings[i.creator.being].id }));
      s.append(bar);
      row += rowH;
    }
    row += gap;
  };
  group(`H₁ · ${H1.length}`, H1, 'h1');
  group(`H₂ · ${H2.length}`, H2, 'h2');
  const filled = H1.filter((i) => i.death !== null);
  return el('div', { class: 'barcode' },
    el('h4', { text: 'The history in time: persistent homology' }),
    el('p', { class: 'note', text: `Each cell enters at its last occasion, and the homology of the history up to each date is followed over ℤ/2 (morphe's persistence). A circle of encounters is born when a being's arrival closes it, and dies when a later meeting's lens fills it; a bar reaching the right edge was never filled. Here ${plural(H1.length, 'circle')} closed and ${plural(filled.length, 'was', 'were')} filled${filled.length ? ` (${listText(filled.map((i) => { const t = data.ties[i.killer.tie]; return `${data.beings[t.a].short} and ${data.beings[t.b].short}`; }))})` : ''}; ${plural(H2.length, 'shell')} formed and none is filled, for the complex has no 3-cells.` }),
    el('div', { class: 'barcode-figure' }, s));
}
