// A result card: what was found, whether it was measured from records, proved for any street, or comes from a
// model with declared parameters, and how it was obtained, with its limits. Every analysis states its results
// through this card.

import { h } from './dom.js';

const KIND = {
  measured: 'measured from records',
  proved: 'holds of any street',
  model: 'a model with declared parameters, to explore a mechanism, not to predict',
};

export function finding({ title, kind = 'measured', statement, method, limits, against, extra }) {
  const ob = h('dl', { class: 'obligations' });
  const row = (k, v) => { if (v) ob.append(h('dt', {}, k), h('dd', {}, v)); };
  row('how', method);
  row('tested against', against);
  row('limits', limits);
  return h('div', { class: 'finding' },
    h('div', { class: 'finding-head' }, h('span', { class: `chip ${kind}`, title: KIND[kind] || '' }, kind), h('b', {}, title)),
    statement ? h('p', {}, statement) : null,
    extra || null,
    ob.children.length ? ob : null);
}
