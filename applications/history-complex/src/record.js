// The record: every assertion of the canon beside its sources, confidence and provenance; the
// problems the construction found; and the check of the beings against Wikidata.

import { el, clear } from './dom.js';
import { colourOf } from './palette.js';
import { dateText, yearText, shortName, signText, number, roleText, listText } from './format.js';

const TABS = [['occasions', 'Occasions'], ['beings', 'Beings'], ['ties', 'Ties'], ['sources', 'Sources'], ['problems', 'Problems']];
const swatch = (colour) => el('span', { class: 'swatch', style: `background:${colourOf(colour).hex}` });

export function renderRecord(host, data, view, act, check) {
  clear(host);
  const tab = TABS.some(([id]) => id === view.tab) ? view.tab : 'occasions';
  const filter = el('input', { type: 'search', placeholder: 'Filter…', class: 'filter', value: act.filterText || '', 'aria-label': 'Filter the rows' });
  const body = el('div', { class: 'record-body' });
  const draw = () => {
    act.filterText = filter.value;
    const q = filter.value.trim().toLowerCase();
    body.replaceChildren(TABLES[tab](data, q, act, check));
  };
  filter.addEventListener('input', draw);
  const counts = { occasions: data.occasions.length, beings: data.beings.length, ties: data.ties.length, sources: data.sources.length, problems: data.problems.length };
  host.append(
    el('div', { class: 'record-head' },
      el('h2', { class: 'view-title', text: 'The record' }),
      el('a', { class: 'link', href: '/domains/history/codebook.md', target: '_blank', rel: 'noopener', text: 'The codebook' })),
    el('p', { class: 'lede narrow', text: 'Assertions, not facts: each says what its sources say, with a confidence, and who drafted and who checked it. The complex is built from these alone.' }),
    el('div', { class: 'record-tools' },
      el('nav', { class: 'subtabs', role: 'tablist' }, ...TABS.map(([id, name]) => el('button', { role: 'tab', 'aria-selected': String(id === tab), onclick: () => act.tab(id) }, `${name} `, el('span', { class: 'count', text: String(counts[id]) })))),
      filter),
    ...(tab === 'beings' ? [checkBar(data, act, check)] : []),
    body);
  draw();
}

const match = (q, ...texts) => !q || texts.some((t) => String(t ?? '').toLowerCase().includes(q));
const cites = (data, ids) => (ids || []).map((id) => data.sourceById.get(id)?.cite || id).join('; ');
const by = (x) => (x ? x.replace(/^model:/, 'model: ') : '—');
const pick = (act, kind, id, text) => el('button', { class: 'link', onclick: () => act.pick({ kind, id }) }, text);

const TABLES = {
  occasions(data, q, act) {
    const rows = data.complex.vertices.map((o) => data.occasions[o]).filter((o) =>
      match(q, o.label, o.place, o.date, o.type, cites(data, o.sources), ...o.present.map((p) => data.beings[p.being].name)));
    return table(['When', 'Occasion', 'Place', 'Present', 'Sources', 'Conf.', 'Drafted · checked'], rows.map((o) => [
      el('span', { title: o.date, text: dateText(o.date, o.calendar) }),
      el('span', {}, pick(act, 'occasion', o.id, o.label), el('span', { class: 'faint', text: ` · ${o.type}` })),
      o.place || '—',
      el('span', {}, ...o.present.flatMap((p, i) => [i ? ', ' : '', el('span', { class: p.derived ? 'derived' : '' }, data.beings[p.being].short, el('span', { class: 'faint', text: ` (${roleText(p.role, p.derived)})` }))])),
      el('span', { class: 'small', text: cites(data, o.sources) || '—' }),
      number(o.confidence, 2),
      el('span', { class: 'small', text: `${by(o.drafted_by)} · ${o.checked_by || 'unchecked'}` }),
    ]), [0, 5]);
  },
  beings(data, q, act, check) {
    const rows = data.beings.map((b, i) => [b, i]).filter(([b]) => match(q, b.name, b.id, b.type, b.wikidata, b.colour));
    return table(['Being', 'Kind', 'Lived', 'Discipline', 'Wikidata', 'Check', 'Occasions', 'Ties'], rows.map(([b, i]) => {
      const c = check?.beings?.[b.id];
      return [
        pick(act, 'being', b.id, b.name), b.type || '—', b.lived ? dateText(b.lived) : '—',
        el('span', {}, swatch(b.colour), ` ${colourOf(b.colour).name}`),
        b.wikidata ? el('a', { href: `https://www.wikidata.org/wiki/${b.wikidata}`, target: '_blank', rel: 'noopener', text: b.wikidata }) : el('span', { class: 'faint', text: '—' }),
        c ? el('span', { title: [c.label, c.description, c.detail].filter(Boolean).join(' — ') }, el('span', { class: `badge ${c.status}`, text: c.status }), c.status === 'proposed' && c.qid ? el('span', { class: 'faint', text: ` ${c.qid}` }) : null) : el('span', { class: 'faint', text: '—' }),
        String(b.occasions.length), String(data.ties.filter((t) => t.a === i || t.b === i).length),
      ];
    }), [6, 7]);
  },
  ties(data, q, act) {
    const rows = data.ties.filter((t) => match(q, t.id, t.type, t.colour, data.beings[t.a].name, data.beings[t.b].name));
    return table(['Tie', 'Type', 'Discipline', 'Signs at its attestations', 'Ends', 'Sources', 'Conf.'], rows.map((t) => [
      pick(act, 'tie', t.id, `${data.beings[t.a].short} — ${data.beings[t.b].short}`),
      t.type, el('span', {}, swatch(t.colour), ` ${colourOf(t.colour).name}`),
      el('span', { class: 'signs' }, ...t.attested.map((a) => el('span', { class: `sign ${a.sign < 0 ? 'neg' : 'pos'}`, title: data.occasions[a.occasion].label, text: `${yearText(data.occasions[a.occasion].date)} ${signText(a.sign)}` }))),
      t.ends !== null && t.ends !== undefined ? data.occasions[t.ends].label : el('span', { class: 'faint', text: 'at its last attestation' }),
      el('span', { class: 'small', text: cites(data, t.sources) || '—' }), number(t.confidence, 2),
    ]), [6]);
  },
  sources(data, q) {
    const uses = new Map();
    for (const x of [...data.occasions, ...data.ties]) for (const s of x.sources || []) uses.set(s, (uses.get(s) || 0) + 1);
    const rows = data.sources.filter((s) => match(q, s.id, s.cite));
    return table(['Source', 'Citation', 'Cited by'], rows.map((s) => [el('code', { text: s.id }), s.cite, String(uses.get(s.id) || 0)]), [2]);
  },
  problems(data, q) {
    const rows = data.problems.filter((p) => match(q, p.where, p.message, p.severity));
    if (!data.problems.length) return el('p', { class: 'note', text: 'The construction found nothing wrong with the record.' });
    return table(['', 'Where', 'What'], rows.map((p) => [el('span', { class: `badge ${p.severity}`, text: p.severity }), el('code', { text: p.where }), p.message]));
  },
};

function table(heads, rows, numeric = []) {
  if (!rows.length) return el('p', { class: 'note', text: 'Nothing matches.' });
  return el('div', { class: 'table-wrap' }, el('table', { class: 'record' },
    el('thead', {}, el('tr', {}, ...heads.map((h, i) => el('th', { class: numeric.includes(i) ? 'num' : '', text: h })))),
    el('tbody', {}, ...rows.map((r) => el('tr', {}, ...r.map((c, i) => el('td', { class: numeric.includes(i) ? 'num' : '' }, c)))))));
}

function checkBar(data, act, check) {
  const s = check?.summary;
  const status = check?.status;
  const text = status === 'running' ? 'Checking against Wikidata… (politely, about two seconds a being)'
    : status === 'done' && s ? `Checked against Wikidata on ${dateText(String(check.checked_at || '').slice(0, 10))}: ${listText(Object.entries(s).filter(([k]) => k.startsWith('beings_')).map(([k, v]) => `${v} ${k.slice(7)}`))}${s.deaths ? `; deaths: ${listText(Object.entries(s).filter(([k]) => k.startsWith('deaths_')).map(([k, v]) => `${v} ${k.slice(7)}`))}` : ''}.`
      : status === 'failed' ? `The last check failed: ${check.message || 'no reason given'}.`
        : check?.fixtures ? 'The fixtures hold no check; the server runs one.'
          : 'Not yet checked against Wikidata.';
  return el('div', { class: 'check-bar' },
    el('span', { text }),
    status === 'running' ? null : el('button', { class: 'button', disabled: check?.fixtures || undefined, onclick: () => act.runCheck() }, status === 'done' ? 'Check again' : 'Check against Wikidata'),
    el('span', { class: 'faint small', text: 'Wikidata (CC0) checks identities and years; the record is corrected by a person, never by the check.' }));
}
