// The spin foam of the history at G = ℤ₂ (the theory note, §6): the couplings fitted to the record;
// the labels of the character expansion and what each means; the cells most often in the foam; the
// memory of each tie across its life (a Wilson loop) against an isolated tie's; the power of each
// being over the reconstruction, computed on request; and SU(2) kinematics.

import { el, svg } from './dom.js';
import { colourOf } from './palette.js';
import { number, plural, listText } from './format.js';
import * as api from './api.js';

const swatch = (colour) => el('span', { class: 'swatch', style: `background:${colourOf(colour).hex}` });

export function renderFoamReadings(data, act, state) {
  const f = data.foam;
  if (!f) return null;
  const B = data.beings, T = data.ties, O = data.occasions;
  const lensCells = f.cells.filter((c) => c.kind === 'lens');
  const support = f.cells.filter((c) => c.support);
  // A cell whose holonomy is all but certain sits at the bound t/(1 + t), so one crowded occasion can fill the
  // list with its triangles: keep the likeliest cell of each occasion and of each tie.
  const seen = new Set();
  const top = [...support].sort((a, b) => b.spin - a.spin).filter((c) => {
    const key = c.kind === 'lens' ? `tie ${data.complex.faces[c.index].tie}` : `occasion ${f.triangles[c.index].occasion}`;
    return !seen.has(key) && seen.add(key);
  }).slice(0, 10);
  const describe = (c) => {
    if (c.kind === 'lens') {
      const face = data.complex.faces[c.index], t = T[face.tie];
      return el('span', {}, swatch(t.colour), ' ', el('button', { class: 'link', onclick: () => act.pick({ kind: 'tie', id: t.id }) }, `${B[t.a].short} — ${B[t.b].short}`),
        el('span', { class: 'faint', text: ` across ${O[face.from].date.slice(0, 7)} → ${O[face.to].date.slice(0, 7)}` }));
    }
    const tr = f.triangles[c.index];
    return el('span', {}, el('span', { class: 'tri', text: '△ ' }), listText(tr.beings.map((b) => B[b].short)), el('span', { class: 'faint', text: ' at ' }),
      el('button', { class: 'link', onclick: () => act.pick({ kind: 'occasion', id: O[tr.occasion].id }) }, O[tr.occasion].label));
  };
  const powerHost = el('div', { class: 'power' }, el('button', { class: 'button', onclick: () => loadPower(data, act, state, powerHost) }, 'Compute the power of each being'),
    el('span', { class: 'faint', text: ' the reconstruction without each being\'s ties, compared with the whole: a few seconds' }));
  if (state.power?.canon === data.canon.id) drawPower(data, act, state.power.answer, powerHost);
  const su2 = f.su2;
  const counts = new Map();
  su2.segments.forEach((s) => { if (s.live) counts.set(s.live, (counts.get(s.live) || 0) + 1); });

  return el('section', { class: 'reading foam-reading' },
    el('h3', { text: 'The spin foam of the history (ℤ₂)' }),
    el('p', { class: 'note', text: 'The signs are a configuration of a Gibbs weight on the refined complex: persistence β on each lens, balance β_Δ on each triangle of ties at an occasion. Its character expansion is a sum over spin foams, the closed ℤ/2 surfaces of the complex, and every label of the social spin foam is derived, not stipulated (the theory note, §6).' }),
    el('div', { class: 'reading two' },
      el('div', {},
        el('h4', { text: 'The couplings fitted to the record' }),
        el('div', { class: 'table-wrap' }, el('table', { class: 'record' },
          el('thead', {}, el('tr', {}, el('th', { text: 'Of' }), el('th', { class: 'num', text: 'Cells' }), el('th', { class: 'num', text: 'Flipped or frustrated' }), el('th', { class: 'num', text: 'tanh β' }), el('th', { class: 'num', text: 'β' }))),
          el('tbody', {}, ...f.couplings.map((c) => el('tr', { class: c.of === 'all' || c.of === 'balance' ? 'total' : '' },
            el('td', {}, c.of === 'balance' ? 'Balance, triangles at occasions' : c.of === 'all' ? 'All ties' : el('span', {}, swatch(c.of), ` ${colourOf(c.of).name}`)),
            el('td', { class: 'num', text: String(c.n) }), el('td', { class: 'num', text: String(c.k) }),
            el('td', { class: 'num', text: number(c.t, 3) }), el('td', { class: 'num', text: number(c.beta, 3) })))))),
        el('p', { class: 'note', text: 'Each from the cells whose signs are all attested, the share flipped or frustrated taken as (k + ½)/(n + 1), so that a coupling stays finite where nothing has yet changed; tanh β = 1 − 2p.' })),
      el('div', {},
        el('h4', { text: 'The labels, and what they mean here' }),
        el('dl', { class: 'meanings' },
          el('dt', { text: 'spin j_f ∈ {0, 1}' }), el('dd', { text: 'whether a cell is in the foam: a tie\'s change across a lens, or a triangle\'s balance, taking part in a closed configuration of changes' }),
          el('dt', { text: 'face amplitude tanh β' }), el('dd', { text: 'how strongly a tie carries its sign across one interval, or a triangle its balance' }),
          el('dt', { text: 'intertwiner, the Gauss law' }), el('dd', { text: 'every rung and life-segment lies on an even number of the foam\'s cells: changes close up through the beings that carry them' }),
          el('dt', { text: 'vertex amplitude' }), el('dd', { text: '1, the group being abelian: the agreement at an occasion is in its triangles' }),
          el('dt', { text: `the foams: ${f.foams}` }), el('dd', { text: `independent closed configurations, a basis of H₂ of the refined complex over ℤ/2; ${support.length} of its ${f.cells.length} cells lie in some foam` })))),
    el('div', { class: 'reading two' },
      el('div', {},
        el('h4', { text: 'The cells most often in the foam' }),
        el('p', { class: 'note', text: 'P(f ∈ S) = t(⟨h_f⟩ − t)/(1 − t²), from the theory\'s own ensemble; a cell in no closed surface is never in a foam, and none exceeds t/(1 + t), the bound a cell whose holonomy is all but certain reaches (Corollary 6.2). The likeliest cell of each occasion and of each tie.' }),
        el('ol', { class: 'spin-list' }, ...top.map((c) => el('li', {}, el('span', { class: 'spin', text: number(c.spin, 2) }), describe(c))))),
      el('div', {},
        el('h4', { text: 'The memory of a tie across its life' }),
        el('p', { class: 'note', text: 'A Wilson loop round the ribbon of each tie, from its first meeting to its last: ⟨u_first u_last⟩ in the theory, beside tanh^k β, an isolated tie\'s over the same k lenses. Where the bar outruns its mark, the triangles the tie belongs to carry its memory: balance as coherence.' }),
        livesChart(data, f.lives))),
    el('div', { class: 'reading two' },
      el('div', {},
        el('h4', { text: 'Power over the reconstruction' }),
        el('p', { class: 'note', text: 'The social spin foam\'s power (its Def. 10.3): the divergence, in nats, of the reconstruction without a being\'s ties from the reconstruction with them, over every other unattested sign; and the power of balance itself, the reconstruction with β_Δ = 0. A property of the fitted model, not a claim about what history would have done.' }),
        powerHost),
      el('div', {},
        el('h4', { text: 'SU(2) kinematics' }),
        el('p', { class: 'note', text: `If every live tie carried spin ½, a being holding n ties across a stretch of its life would have dim Inv(V_½^⊗n) ways to hold them together: the Catalan number C_{n/2}, none for n odd. Of the ${su2.carrying} stretches that carry ties, ${su2.admissible} are admissible at spin ½; the others need a tie of integer spin, a relation of another kind. The record does not yet measure SU(2), so this is kinematics, not a fit (the theory note, §6.5).` }),
        el('div', { class: 'bars' }, ...[...counts.entries()].sort((a, b) => a[0] - b[0]).map(([n, k]) => el('div', { class: 'bar-row' },
          el('span', { class: 'bar-name', text: `${n} ${n === 1 ? 'tie' : 'ties'}${n % 2 ? '' : `, ${catalan(n / 2)} channel${catalan(n / 2) === 1 ? '' : 's'}`}` }),
          el('span', { class: 'bar-track' }, el('span', { class: `bar ${n % 2 ? 'odd' : 'even'}`, style: `width:${(100 * k) / su2.carrying}%` })),
          el('span', { class: 'bar-value', text: String(k) })))))));
}

function catalan(n) { let c = 1; for (let k = 0; k < n; k++) c = (c * 2 * (2 * k + 1)) / (k + 2); return Math.round(c); }

function livesChart(data, lives) {
  const rows = [...lives].sort((a, b) => b.lenses - a.lenses || b.wilson - a.wilson).slice(0, 14);
  const W = 560, rowH = 20, L = 250, R = 30;
  const s = svg('svg', { viewBox: `0 0 ${W} ${rows.length * rowH + 24}`, class: 'lives-svg', role: 'img', 'aria-label': 'Wilson loops of the ties\' lives' });
  const X = (v) => L + Math.max(0, Math.min(1, v)) * (W - L - R);
  s.append(svg('line', { x1: X(0), x2: X(0), y1: 4, y2: rows.length * rowH + 4, class: 'gridline' }), svg('line', { x1: X(1), x2: X(1), y1: 4, y2: rows.length * rowH + 4, class: 'gridline' }),
    svg('text', { x: X(0), y: rows.length * rowH + 18, 'text-anchor': 'middle', class: 'tick' }, '0'), svg('text', { x: X(1), y: rows.length * rowH + 18, 'text-anchor': 'middle', class: 'tick' }, '1'));
  rows.forEach((l, i) => {
    const t = data.ties[l.tie], y = 6 + i * rowH;
    s.append(svg('text', { x: L - 8, y: y + 11, 'text-anchor': 'end', class: 'row-name' }, `${data.beings[t.a].short} — ${data.beings[t.b].short} · ${l.lenses}`),
      svg('rect', { x: X(0), y: y + 2, width: Math.max(1, X(l.wilson) - X(0)), height: 11, rx: 1.5, fill: colourOf(t.colour).hex, opacity: 0.75 }),
      svg('line', { x1: X(l.chain), x2: X(l.chain), y1: y, y2: y + 15, class: 'chain-mark' }));
  });
  return el('div', { class: 'lives' }, s, el('p', { class: 'faint small', text: 'Bars: the tie\'s memory in the theory; marks: an isolated tie\'s. The number after the names counts the lenses of its life.' }));
}

async function loadPower(data, act, state, host) {
  host.replaceChildren(el('p', { class: 'note', text: 'Reconstructing the history without each being in turn…' }));
  try {
    const answer = await api.power(data.canon.id);
    state.power = { canon: data.canon.id, answer };
    drawPower(data, act, answer, host);
  } catch (error) {
    host.replaceChildren(el('p', { class: 'error', text: error.message }));
  }
}

function drawPower(data, act, answer, host) {
  const rows = [...answer.beings].sort((a, b) => b.divergence - a.divergence).slice(0, 12);
  const max = Math.max(answer.balance.divergence, ...rows.map((r) => r.divergence), 1e-9);
  host.replaceChildren(
    el('div', { class: 'bars' },
      el('div', { class: 'bar-row total' }, el('span', { class: 'bar-name', text: 'Balance itself' }),
        el('span', { class: 'bar-track' }, el('span', { class: 'bar balance', style: `width:${(100 * answer.balance.divergence) / max}%` })),
        el('span', { class: 'bar-value', text: number(answer.balance.divergence, 2) })),
      ...rows.map((r) => el('div', { class: 'bar-row' },
        el('span', { class: 'bar-name' }, el('button', { class: 'link', onclick: () => act.pick({ kind: 'being', id: data.beings[r.being].id }) }, data.beings[r.being].name)),
        el('span', { class: 'bar-track' }, el('span', { class: 'bar', style: `width:${(100 * r.divergence) / max}%;background:${colourOf(data.beings[r.being].colour).hex}` })),
        el('span', { class: 'bar-value', text: number(r.divergence, 3) })))),
    el('p', { class: 'faint small', text: `Nats over ${plural(answer.balance.rungs, 'unattested sign')}; ${answer.sweeps} sweeps of the heat bath for each reconstruction${answer.computed < answer.considered ? `; for the ${answer.computed} beings with the most ties, of ${answer.considered}` : ''}.` }));
}
