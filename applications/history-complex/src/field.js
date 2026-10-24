// The field theory of a history (the theory note, §8): the history read as one Feynman graph of a
// group field theory on ℤ₂, the action whose expansion holds it, its amplitude, and the history in
// time as one term of a Dyson series, with the rates of its occasions and strands and the
// martingale checks of their fit. The engine's field request answers it, with the couplings the
// build fitted.

import { el, svg, clear } from './dom.js';
import { colourOf } from './palette.js';
import { number, plural, yearName } from './format.js';
import * as api from './api.js';

const swatch = (colour) => el('span', { class: 'swatch', style: `background:${colourOf(colour).hex}` });
// Nodes as a list in words: "A", "A and B", "A, B and C".
const joined = (nodes) => nodes.flatMap((n, i) => (i === 0 ? [n] : [i === nodes.length - 1 ? ' and ' : ', ', n]));

// The probability that a χ² variable on d degrees of freedom exceeds x: Q(d/2, x/2), from
// Q(½, y) = erfc √y and Q(1, y) = e^{−y} by Q(a + 1, y) = Q(a, y) + y^a e^{−y}/Γ(a + 1).
function chi2Tail(x, d) {
  if (!(d > 0) || !(x >= 0)) return NaN;
  const y = x / 2;
  const erfc = (z) => { // Abramowitz and Stegun 7.1.26, |error| < 1.5e-7
    const t = 1 / (1 + 0.3275911 * z);
    return t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * Math.exp(-z * z);
  };
  let a = d % 2 === 0 ? 1 : 0.5, q = d % 2 === 0 ? Math.exp(-y) : erfc(Math.sqrt(y));
  let gamma = d % 2 === 0 ? 1 : Math.sqrt(Math.PI);  // Γ(a)
  while (a < d / 2) {
    q += Math.pow(y, a) * Math.exp(-y) / (gamma * a);
    gamma *= a;
    a += 1;
  }
  return Math.min(1, Math.max(0, q));
}

// One interaction vertex: each being present a column, the segment arriving below and the one
// leaving above; the corners its strands turn at, coloured by the tie's discipline: closing ones
// below (a lens ends), opening ones above (the next begins), passing ones along a column.
function bubbleDrawing(data, b) {
  const beings = [...new Set(b.half_edges.map((h) => h.being))];
  const column = new Map(beings.map((x, i) => [x, i]));
  const W = Math.max(120, 36 + beings.length * 52), yin = 74, yout = 26;
  const X = (being) => 26 + column.get(being) * 52;
  const s = svg('svg', { viewBox: `0 0 ${W} 112`, class: 'bubble-svg', role: 'img', 'aria-label': 'An interaction vertex' });
  const passing = new Map();
  for (const c of b.corners) {
    const tie = data.ties[c.tie], hex = colourOf(tie.colour).hex;
    const [h1, h2] = c.ends.map((k) => b.half_edges[k]);
    if (c.kind === 'passing') {
      const n = passing.get(h1.being) || 0;
      passing.set(h1.being, n + 1);
      const x = X(h1.being) + (n % 2 ? 1 : -1) * Math.ceil(n / 2) * 3;
      s.append(svg('line', { x1: x, x2: x, y1: yin, y2: yout, stroke: hex, 'stroke-width': 1.4, 'stroke-dasharray': '3 2.5', opacity: 0.8 }));
      continue;
    }
    const x1 = X(h1.being), x2 = X(h2.being), y = c.kind === 'closing' ? yin : yout;
    const bow = (c.kind === 'closing' ? 1 : -1) * (12 + 5 * Math.abs(column.get(h1.being) - column.get(h2.being)));
    s.append(svg('path', { d: `M${x1} ${y} Q${(x1 + x2) / 2} ${y + bow} ${x2} ${y}`, fill: 'none', stroke: hex, 'stroke-width': 1.7 }));
  }
  for (const h of b.half_edges)
    s.append(svg('circle', { cx: X(h.being), cy: h.out ? yout : yin, r: 3.4, class: h.out ? 'half out' : 'half in' }));
  for (const x of beings) {
    const name = data.beings[x].short;
    s.append(svg('text', { x: X(x), y: 108, 'text-anchor': 'middle', class: 'bubble-name' }, name.length > 10 ? `${name.slice(0, 9)}…` : name));
  }
  return s;
}

export function renderField(host, data, act, state) {
  clear(host);
  const key = data.canon.id;
  if (state.field?.canon !== key) {
    host.append(el('h2', { class: 'view-title', text: 'The field theory of the history' }),
      el('p', { class: 'note', text: 'Reading the history as a Feynman graph and as a sample path…' }));
    if (state.fieldLoading !== key) {
      state.fieldLoading = key;
      api.field(key, api.couplingsOf(data))
        .then((answer) => {
          state.field = { canon: key, answer };
          state.fieldLoading = null;
          if (state.view === 'field' && state.data?.canon.id === key) renderField(host, state.data, act, state);
        })
        .catch((error) => {
          state.fieldLoading = null;
          host.replaceChildren(el('h2', { class: 'view-title', text: 'The field theory of the history' }), el('p', { class: 'error', text: error.message }));
        });
    }
    return;
  }
  const F = state.field.answer, foam = data.foam, P = F.process;
  const couplings = (foam?.couplings || []).filter((c) => c.of !== 'all');
  const lens = couplings.filter((c) => c.of !== 'balance'), balance = couplings.find((c) => c.of === 'balance');
  const occasions = F.occasions, types = F.types;
  const once = types.filter((t) => t.count === 1).length;
  const maxRank = Math.max(1, ...F.ranks.map((r) => r.segments));
  const tExact = (c) => (c.low === null || c.low === undefined ? number(c.t, 2) : `${number(c.t, 2)} (${number(c.low, 2)} to ${number(c.high, 2)})`);

  host.append(
    el('h2', { class: 'view-title', text: 'The field theory of the history' }),
    el('p', { class: 'lede narrow', text: 'A history is read here as one Feynman graph of a group field theory, and in time as one term of its Dyson series. The beings between occasions are the field\'s quanta, each with an argument for every tie it holds; the occasions are the interaction vertices; the lenses of the ties are the faces the strands close. The theory note\'s §8 states it; every identity below is checked by the domain on every canon.' }),

    // 1. The Feynman graph.
    el('section', { class: 'field-part' },
      el('h3', { text: '1. The history as a Feynman graph' }),
      el('div', { class: 'reading two' },
        el('div', {},
          el('dl', { class: 'figures' },
            el('dt', { text: 'Vertices: occasions' }), el('dd', { text: String(occasions) }),
            el('dt', { text: 'Propagators: segments of lives' }), el('dd', { text: String(F.segments) }),
            el('dt', { text: 'Strands: ties carried along them' }), el('dd', { text: String(F.strands) }),
            el('dt', { text: 'Faces: lenses, closed strands' }), el('dd', { text: String(F.lenses) })),
          el('p', { class: F.identity.holds ? 'note' : 'error', text: F.identity.holds
            ? `Proposition 8.1 holds here. Eliminating the δ's of the Feynman rules over ℤ₂ leaves one free strand variable for each of the ${F.lenses} lenses, and each lens's kernel takes exactly that lens's holonomy on the refined complex. So the amplitude of this graph is 2^(F + E) times the weight of the spin foam.`
            : 'Proposition 8.1 fails on this record: the δ\'s do not reduce to the lenses\' holonomies. That is an error to report.' }),
          el('h4', { text: 'The field\'s ranks' }),
          el('p', { class: 'note', text: 'A segment of a life carries one argument for each tie alive across it: a field φ(g₁, …, gₙ) of rank n, the n of the SU(2) kinematics in the readings.' }),
          el('div', { class: 'bars' }, ...F.ranks.map((r) => el('div', { class: 'bar-row' },
            el('span', { class: 'bar-name', text: `rank ${r.rank}` }),
            el('span', { class: 'bar-track' }, el('span', { class: 'bar', style: `width:${(100 * r.segments) / maxRank}%` })),
            el('span', { class: 'bar-value', text: String(r.segments) }))))),
        el('div', {},
          el('h4', { text: 'The interaction vertices' }),
          el('p', { class: 'note', text: `Each occasion's vertex is the link of the occasion in the complex: its arriving and leaving segments (dots below and above), and the corners its strands turn at, where a lens closes (below), the next opens (above) or a lens runs past (dashed). Up to isomorphism there are ${types.length} kinds of vertex among the ${occasions} occasions, and ${once} occur only once. Read literally, the action has almost as many couplings as the history has vertices.` }),
          el('div', { class: 'bubbles' }, ...types.slice(0, 8).map((t) => {
            const b = F.bubbles[t.example], o = data.occasions[b.occasion];
            return el('figure', { class: 'bubble' }, bubbleDrawing(data, b),
              el('figcaption', {},
                el('span', { class: 'num', text: `${t.count}×` }), ` ${plural(t.beings, 'being')}, |Aut| ${t.canonical ? t.automorphisms : '?'}; `,
                el('button', { class: 'link', onclick: () => act.pick({ kind: 'occasion', id: o.id }) }, o.label.length > 60 ? `${o.label.slice(0, 58)}…` : o.label)));
          }))))),

    // 2. The action.
    el('section', { class: 'field-part' },
      el('h3', { text: '2. The action whose expansion holds it' }),
      el('div', { class: 'formula', text: 'S[φ] = Σₙ ½ ⟨φₙ, Pₙ⁻¹ φₙ⟩ − Σ_B (λ_B / |Aut B|) V_B[φ]' }),
      el('dl', { class: 'meanings' },
        el('dt', { text: 'φₙ' }), el('dd', { text: 'a real field of rank n on ℤ₂ⁿ, invariant under the diagonal action (the closure of a quantum); one argument for each tie, coloured by its discipline' }),
        el('dt', { text: 'Pₙ' }), el('dd', {}, 'the propagator: the projector on invariant fields, with the persistence of each strand\'s discipline, ',
          ...(lens.length ? joined(lens.map((c) => el('span', {}, swatch(c.of), ` ${colourOf(c.of).name.toLowerCase()}, tanh β = ${tExact(c)}`))) : ['none fitted']), '.'),
        el('dt', { text: 'V_B' }), el('dd', { text: `the vertex of a bubble B: a δ at each passing and opening corner, the Gibbs kernel of its lens at each closing corner, and the balance of each triangle of ties met there, tanh β_Δ = ${balance ? tExact(balance) : '—'}` }),
        el('dt', { text: 'λ_B' }), el('dd', { text: 'the coupling of each kind of vertex: as good as unidentifiable from one history, where nearly every vertex is its own kind; §4 estimates the coarse ones, by the number of beings an occasion brings from before' })),
      el('p', { class: 'note', text: `In the expansion Z = ∫ Dφ e^(−S) = Σ_Γ Π λ / |Aut Γ| · A_Γ, this history is one term Γ. Its Feynman amplitude A_Γ is the spin foam of the readings. ${foam?.exact ? 'The couplings above are the maximum a posteriori of the exact likelihood of the record\'s signs, with their 95% intervals.' : 'The couplings above are the KT estimates: the exact likelihood is beyond reach here.'}` })),

    // 3. The amplitude.
    el('section', { class: 'field-part' },
      el('h3', { text: '3. Its amplitude' }),
      F.amplitude ? el('dl', { class: 'figures narrow' },
        el('dt', { text: 'log Z, the spin foam\'s weight (temporal gauge)' }), el('dd', { text: number(F.amplitude.log_z, 3) }),
        el('dt', { text: 'log A_Γ = log Z + (F + E) log 2' }), el('dd', { text: number(F.amplitude.log_amplitude, 3) }),
        el('dt', { text: 'log P(the attested signs)' }), el('dd', { text: number(F.amplitude.log_likelihood, 3) }),
        el('dt', { text: 'largest state of the contraction' }), el('dd', { text: `${F.amplitude.largest} characters` }))
        : el('p', { class: 'note', text: 'Beyond exact reach on this record: the boundary state of some cut is too large to contract. The readings sample it instead.' })),

    // 4. The Dyson series.
    el('section', { class: 'field-part' },
      el('h3', { text: '4. The history in time: one term of a Dyson series' }),
      el('p', { class: 'note narrow', text: 'Time-ordered, the history is a sample path. Its quanta are the beings between their first and last occasions; an occasion annihilates the segments arriving and creates those leaving, and the beings appearing for the first time are created from the vacuum. Its weight factorises into the law of the occasions, the law of the strands (which ties begin and end) and the law of the signs (§3). The first two are fitted here by maximum likelihood, as a marked point process.' }),
      el('div', { class: 'reading two' },
        el('div', {},
          el('h4', { text: 'The rates of the occasions' }),
          el('div', { class: 'table-wrap' }, el('table', { class: 'record' },
            el('thead', {}, el('tr', {}, el('th', { class: 'num', text: 'k arriving' }), el('th', { class: 'num', text: 'Occasions' }), el('th', { class: 'num', text: 'Exposure' }),
              el('th', { class: 'num', text: 'λ̂ per k-subset a year' }), el('th', { class: 'num', text: 'New, on average' }))),
            el('tbody', {}, ...P.rates.map((r) => el('tr', {},
              el('td', { class: 'num', text: String(r.k) }), el('td', { class: 'num', text: String(r.count) }), el('td', { class: 'num', text: number(r.exposure, 3) }),
              el('td', { class: 'num', text: r.count ? `${number(r.rate, 4)} ± ${number(r.error, 4)}` : `< ${number(r.upper, 4)}` }),
              el('td', { class: 'num', text: number(r.mean_new, 2) })))))),
          el('p', { class: 'note', text: 'An occasion bringing k beings from before happens at rate λ_k for each k of the beings then in the record: intensity λ_k C(n(τ), k). The exposure is ∫ C(n(τ), k) dτ, in k-subset-years; λ̂ = count / exposure, with its Poisson standard error.' })),
        el('div', {},
          el('h4', { text: 'The strands, and how ties bind' }),
          el('dl', { class: 'figures' },
            el('dt', { text: 'A tie begins, two continuing beings meeting untied' }), el('dd', { text: P.open.trials ? `${number(P.open.p, 3)} ± ${number(P.open.error, 3)} (${P.open.successes} of ${P.open.trials})` : '—' }),
            el('dt', { text: 'A tie met again ends there' }), el('dd', { text: P.close.trials ? `${number(P.close.p, 3)} ± ${number(P.close.error, 3)} (${P.close.successes} of ${P.close.trials})` : '—' }),
            el('dt', { text: 'A tied pair meets, a year' }), el('dd', { text: number(P.binding.tied, 3) }),
            el('dt', { text: 'An untied pair meets, a year' }), el('dd', { text: number(P.binding.untied, 4) }),
            el('dt', { text: 'Binding, ρ = μ₁/μ₀' }), el('dd', { text: P.binding.ratio ? `${number(P.binding.ratio, 1)} (${number(P.binding.low, 1)} to ${number(P.binding.high, 1)})` : '—' })),
          el('p', { class: 'note', text: 'Two beings in the record meet at rate μ₁ while a tie between them is live and μ₀ while none is: the ratio says how strongly a tie holds its beings together. The pairs at one occasion are not independent, so this is a composite likelihood and its interval is indicative.' }))),
      el('h4', { text: 'The martingale checks, era by era' }),
      el('p', { class: 'note narrow', text: 'The count of occasions less its compensator, ∫ intensity, has mean zero over any window and variance the compensator (Dynkin\'s formula, the stochastic counterpart of the Schwinger–Dyson equations). Fitted on the whole history, the rates are tested here on each third of it, where they were not fitted one by one.' }),
      el('div', { class: 'table-wrap' }, el('table', { class: 'record' },
        el('thead', {}, el('tr', {}, el('th', { text: 'Era' }), el('th', { class: 'num', text: 'Occasions' }), el('th', { class: 'num', text: 'Expected' }), el('th', { class: 'num', text: 'z' }),
          el('th', { class: 'num', text: 'Ties begun' }), el('th', { class: 'num', text: 'Expected' }), el('th', { class: 'num', text: 'z' }))),
        el('tbody', {}, ...P.eras.map((e) => el('tr', {},
          el('td', { text: `${yearName(Math.floor(e.from))} to ${yearName(Math.floor(e.to))}` }),
          el('td', { class: 'num', text: String(e.occasions) }), el('td', { class: 'num', text: number(e.expected, 1) }),
          el('td', { class: `num${Math.abs(e.z) > 2 ? ' strong' : ''}`, text: number(e.z, 2) }),
          el('td', { class: 'num', text: String(e.opened) }), el('td', { class: 'num', text: number(e.opened_expected, 1) }),
          el('td', { class: `num${Math.abs(e.opened_z) > 2 ? ' strong' : ''}`, text: number(e.opened_z, 2) })))))),
      el('p', { class: 'note', text: `χ² = ${number(P.chi2, 2)} on ${P.dof} degrees of freedom, p = ${number(chi2Tail(P.chi2, P.dof), 3)}. ${chi2Tail(P.chi2, P.dof) < 0.05 ? 'Constant rates are rejected: the canon\'s occasions are denser in some eras than its beings explain, as a canon that follows its sources\' attention would be.' : 'Constant rates are not rejected on this canon.'}` }),
      el('p', { class: 'note', text: `The log-likelihood of the history, part by part: occasions ${number(P.log_likelihood.occasions, 1)}, strands ${number(P.log_likelihood.strands, 1)}${F.amplitude ? `, signs ${number(F.amplitude.log_likelihood, 1)}` : ''}.` })),

    el('section', { class: 'field-part' },
      el('h3', { text: '5. What it licenses' }),
      el('p', { class: 'note narrow', text: 'The reading of a history as a Feynman graph is an identity of amplitudes, exact and checked. Its couplings are fitted to one history, recorded by one canon whose selection the canon states: they describe the record, not the past, and predict nothing until another record tests them.' })));
}
