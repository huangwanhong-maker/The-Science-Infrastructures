// The cards of a being, a tie and an occasion: what the record says of it, with its sources,
// confidence and provenance, and what the construction made of it.

import { el } from './dom.js';
import { colourOf } from './palette.js';
import { dateText, yearText, shortName, signText, number, roleText, plural, listText } from './format.js';

const swatch = (colour) => el('span', { class: 'swatch', style: `background:${colourOf(colour).hex}` });
const kv = (k, ...v) => el('div', { class: 'kv' }, el('span', { class: 'k', text: k }), el('span', { class: 'v' }, ...v));
const link = (text, onclick, title) => el('button', { class: 'link', onclick, title }, text);

function provenance(x) {
  const bits = [];
  if (x.drafted_by) bits.push(`drafted by ${x.drafted_by.replace(/^model:/, 'a language model, ')}`);
  bits.push(x.checked_by ? `checked by ${x.checked_by}` : 'not yet checked against its source');
  return bits.join('; ');
}
function sources(data, ids) {
  if (!ids?.length) return el('span', { class: 'faint', text: 'none named' });
  return el('ul', { class: 'cites' }, ...ids.map((id) => el('li', { text: data.sourceById.get(id)?.cite || id })));
}

export function occasionCard(data, o, act) {
  const oc = data.occasions[o];
  const ties = data.ties.map((t, i) => [t, i]).filter(([t]) => t.attested.some((a) => a.occasion === o));
  const triangles = data.readings.triangles.filter((t) => t.occasion === o);
  return el('div', { class: 'card' },
    el('p', { class: 'eyebrow', text: `Occasion · ${oc.type || 'occasion'}` }),
    el('h2', { text: oc.label }),
    el('p', { class: 'when', text: `${dateText(oc.date, oc.calendar)}${oc.place ? ` · ${oc.place}` : ''}` }),
    el('h4', { text: 'Present' }),
    el('ul', { class: 'people' }, ...oc.present.map((p) => el('li', {},
      swatch(data.beings[p.being].colour),
      link(data.beings[p.being].name, () => act.pick({ kind: 'being', id: data.beings[p.being].id })),
      el('span', { class: p.derived ? 'role derived' : 'role', text: roleText(p.role, p.derived) })))),
    ties.length ? el('h4', { text: 'Ties attested here' }) : null,
    ties.length ? el('ul', { class: 'people' }, ...ties.map(([t]) => {
      const s = t.attested.find((a) => a.occasion === o).sign;
      return el('li', {}, el('span', { class: `sign ${s < 0 ? 'neg' : 'pos'}`, text: signText(s) }),
        link(`${data.beings[t.a].short} — ${data.beings[t.b].short}`, () => act.pick({ kind: 'tie', id: t.id })),
        el('span', { class: 'role', text: t.type }));
    })) : null,
    triangles.length ? el('h4', { text: 'Triangles of ties here' }) : null,
    triangles.length ? el('ul', { class: 'people' }, ...triangles.map((t) => el('li', {},
      el('span', { class: t.product < 0 ? 'tri frustrated' : 'tri', text: t.product < 0 ? '▲' : '△' }),
      el('span', { text: listText(t.beings.map((b) => data.beings[b].short)) }),
      el('span', { class: 'role', text: `${t.product < 0 ? 'frustrated' : 'balanced'}${t.attested ? '' : ', with signs carried'}; expected ${number(t.expected, 2)}` })))) : null,
    el('h4', { text: 'Record' }),
    kv('Sources', sources(data, oc.sources)),
    kv('Confidence', number(oc.confidence, 2)),
    kv('Provenance', provenance(oc)),
    oc.note ? kv('Note', oc.note) : null,
    el('div', { class: 'row' }, el('button', { class: 'button', onclick: () => act.moment(oc.date) }, 'The moment at this date')));
}

export function tieCard(data, ti, act) {
  const t = data.ties[ti];
  const A = data.beings[t.a], B = data.beings[t.b];
  const faces = data.complex.faces.filter((f) => f.tie === ti);
  const occ = data.occasions;
  return el('div', { class: 'card' },
    el('p', { class: 'eyebrow' }, swatch(t.colour), ` Tie · ${colourOf(t.colour).name.toLowerCase()} · ${t.type}`),
    el('h2', {}, link(A.name, () => act.pick({ kind: 'being', id: A.id })), ' — ', link(B.name, () => act.pick({ kind: 'being', id: B.id }))),
    el('p', { class: 'when', text: `lives from ${dateText(occ[t.life[0]].date)} to ${dateText(occ[t.life[1]].date)}${t.ends !== null && t.ends !== undefined ? `, ending at “${occ[t.ends].label}”` : ''}` }),
    el('p', { class: 'note', text: `${plural(faces.length, 'lens', 'lenses')}, ${plural(faces.filter((f) => f.flip).length, 'flip')}: a lens is the 2-cell between two meetings of the two, bounded by their two lives.` }),
    el('h4', { text: 'At each meeting' }),
    el('table', { class: 'rungs-table' },
      el('tbody', {}, ...t.rungs.map((r) => el('tr', { class: r.attested ? '' : 'carried' },
        el('td', { class: 'num', text: yearText(occ[r.occasion].date) }),
        el('td', {}, link(occ[r.occasion].label, () => act.pick({ kind: 'occasion', id: occ[r.occasion].id }))),
        el('td', { class: `sign ${r.sign < 0 ? 'neg' : 'pos'}`, text: signText(r.sign) }),
        el('td', { class: 'faint', text: r.attested ? 'attested' : `carried, expected ${number(r.expected, 2)}` }))))),
    el('h4', { text: 'Record' }),
    kv('Sources', sources(data, t.sources)),
    kv('Confidence', number(t.confidence, 2)),
    kv('Provenance', provenance(t)),
    t.note ? kv('Note', t.note) : null);
}

export function beingCard(data, b, act, check) {
  const being = data.beings[b];
  const ties = data.ties.map((t, i) => [t, i]).filter(([t]) => t.a === b || t.b === b);
  const os = being.occasions.map((o) => data.occasions[o]);
  const c = check?.beings?.[being.id];
  return el('div', { class: 'card' },
    el('p', { class: 'eyebrow' }, swatch(being.colour), ` ${being.type || 'being'} · ${colourOf(being.colour).name.toLowerCase()}`),
    el('h2', { text: being.name }),
    el('p', { class: 'when', text: [being.lived ? `lived ${dateText(being.lived)}` : null, os.length ? `in the record from ${yearText(os[0].date)} to ${yearText(os[os.length - 1].date)}` : null].filter(Boolean).join(' · ') }),
    being.wikidata ? kv('Wikidata', el('a', { href: `https://www.wikidata.org/wiki/${being.wikidata}`, target: '_blank', rel: 'noopener', text: being.wikidata })) : kv('Wikidata', el('span', { class: 'faint', text: 'no identifier yet' })),
    c ? kv('Check', el('span', { class: `badge ${c.status}`, text: c.status }), c.detail ? el('span', { class: 'faint', text: ` ${c.detail}` }) : null) : null,
    el('h4', { text: plural(ties.length, 'tie') }),
    el('ul', { class: 'people' }, ...ties.map(([t]) => {
      const other = data.beings[t.a === b ? t.b : t.a];
      const last = t.attested[t.attested.length - 1];
      return el('li', {}, swatch(t.colour), link(other.name, () => act.pick({ kind: 'tie', id: t.id })),
        el('span', { class: 'role', text: `${t.type}, last ${signText(last?.sign)}` }));
    })),
    el('h4', { text: plural(os.length, 'occasion') }),
    el('ul', { class: 'occasions' }, ...being.occasions.map((o) => el('li', {},
      el('span', { class: 'num', text: yearText(data.occasions[o].date) }),
      link(data.occasions[o].label, () => act.pick({ kind: 'occasion', id: data.occasions[o].id }))))));
}

// The card of the canon, when nothing is followed.
export function canonCard(data) {
  const c = data.canon, [V, E, F] = data.invariants.counts;
  return el('div', { class: 'card' },
    el('p', { class: 'eyebrow', text: `Canon · ${c.status || 'seed'}` }),
    el('h2', { text: c.title }),
    c.about ? el('p', { class: 'lede', text: c.about }) : null,
    el('p', { class: 'note', text: `${plural(data.beings.length, 'being')}, ${plural(data.occasions.length, 'occasion')}, ${plural(data.ties.length, 'tie')}: a complex of ${V} vertices, ${E} edges and ${F} lenses.` }),
    data.problems.length ? el('p', { class: 'note', text: `The construction noted ${data.problems.filter((p) => p.severity === 'error').length ? `${data.problems.filter((p) => p.severity === 'error').length} errors and ` : ''}${data.problems.filter((p) => p.severity === 'warning').length} warnings (single meetings, overlapping dates): see the record's problems.` }) : null,
    c.drafted_by ? el('p', { class: 'note', text: `Drafted by ${c.drafted_by.replace(/^model:/, 'a language model, ')}${c.drafted ? ` on ${dateText(c.drafted)}` : ''}; coded by the codebook, version ${c.codebook ?? '?'}. Its assertions are proposals until checked against their sources.` }) : null);
}
