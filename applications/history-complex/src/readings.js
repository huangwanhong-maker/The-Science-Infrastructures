// The readings: the signed ties read as a ℤ₂ field in temporal gauge. The persistence of ties by
// discipline (flips across lenses whose two ends are attested: p, the memory 1 − 2p = tanh β, and
// β), the triangles of ties at occasions, and how the disciplines' lenses meet.

import { el, clear } from './dom.js';
import { COLOURS, colourOf } from './palette.js';
import { dateText, shortName, number, listText, signText } from './format.js';
import { renderFoamReadings } from './foamreadings.js';

// β = artanh(memory): infinite when no lens flips (memory 1), none when the memory is gone (≤ 0).
const betaText = (e) => (!e.lenses ? '—' : e.memory >= 1 ? '∞' : e.beta === null || e.beta === undefined ? 'none' : number(e.beta, 3));
const swatch = (colour) => el('span', { class: 'swatch', style: `background:${colourOf(colour).hex}` });

export function renderReadings(host, data, act, state = {}) {
  clear(host);
  const r = data.readings;
  const per = r.persistence;
  const onlyFrustrated = el('input', { type: 'checkbox', checked: act.onlyFrustrated ?? true });
  const triangles = el('div');
  const drawTriangles = () => {
    act.onlyFrustrated = onlyFrustrated.checked;
    const rows = r.triangles.filter((t) => !onlyFrustrated.checked || t.product < 0);
    triangles.replaceChildren(rows.length ? el('div', { class: 'table-wrap' }, el('table', { class: 'record' },
      el('thead', {}, el('tr', {}, el('th', { text: 'When' }), el('th', { text: 'Occasion' }), el('th', { text: 'Beings' }), el('th', { class: 'num', text: 'Product' }), el('th', { class: 'num', text: 'Expected' }), el('th', { text: 'Signs' }))),
      el('tbody', {}, ...rows.map((t) => {
        const o = data.occasions[t.occasion];
        return el('tr', {},
          el('td', { text: dateText(o.date, o.calendar) }),
          el('td', {}, el('button', { class: 'link', onclick: () => act.pick({ kind: 'occasion', id: o.id }) }, o.label)),
          el('td', { text: listText(t.beings.map((b) => data.beings[b].short)) }),
          el('td', { class: `num sign ${t.product < 0 ? 'neg' : 'pos'}` }, t.product < 0 ? el('span', { class: 'tri frustrated', text: '▲ ' }) : null, signText(t.product)),
          el('td', { class: 'num', text: number(t.expected, 2) }),
          el('td', { class: 'faint', text: t.attested ? 'all attested' : 'some carried' }));
      })))) : el('p', { class: 'note', text: 'No triangle of ties at an occasion is frustrated.' }));
  };
  onlyFrustrated.addEventListener('change', drawTriangles);
  drawTriangles();

  const maxLenses = Math.max(1, ...Object.values(r.colours.lenses));
  host.append(
    el('h2', { class: 'view-title', text: 'The readings' }),
    el('p', { class: 'lede narrow', text: 'The signs of the ties read as a ℤ₂ gauge field in temporal gauge: a being keeps its frame between occasions, a lens carries the change of its tie between two meetings, and a triangle of ties at an occasion carries their product. Balance is flatness; frustration is curvature.' }),
    renderFoamReadings(data, act, state),
    el('section', { class: 'reading' },
      el('h3', { text: 'Persistence of ties, by discipline' }),
      el('p', { class: 'note', text: 'Counted over the lenses whose two ends are both attested: p is the share that flip; the memory 1 − 2p is tanh β, the correlation of a tie\'s sign across one lens in a ℤ₂ theory of stiffness β; a sign carried over k unattested meetings is trusted as memory^k.' }),
      el('div', { class: 'table-wrap' }, el('table', { class: 'record' },
        el('thead', {}, el('tr', {}, el('th', { text: 'Discipline' }), el('th', { class: 'num', text: 'Lenses' }), el('th', { class: 'num', text: 'Flips' }), el('th', { class: 'num', text: 'p' }), el('th', { class: 'num', text: 'Memory, tanh β' }), el('th', { class: 'num', text: 'β' }))),
        el('tbody', {}, ...per.map((e) => el('tr', { class: e.colour === 'all' ? 'total' : '' },
          el('td', {}, e.colour === 'all' ? 'All disciplines' : el('span', {}, swatch(e.colour), ` ${colourOf(e.colour).name}`)),
          el('td', { class: 'num', text: String(e.lenses) }), el('td', { class: 'num', text: String(e.flips) }),
          el('td', { class: 'num', text: number(e.p, 3) }), el('td', { class: 'num', text: number(e.memory, 3) }),
          el('td', { class: 'num', text: betaText(e) }))))))),
    el('section', { class: 'reading' },
      el('div', { class: 'reading-head' }, el('h3', { text: 'Triangles of ties at occasions' }), el('label', { class: 'check' }, onlyFrustrated, ' frustrated only')),
      el('p', { class: 'note', text: 'Three beings present at one occasion with three live ties among them: the spatial plaquettes of the history. Their product is −1 where the triangle is frustrated; expected, where some sign is carried rather than attested there.' }),
      triangles),
    el('section', { class: 'reading two' },
      el('div', {},
        el('h3', { text: 'Lenses by discipline' }),
        el('div', { class: 'bars' }, ...COLOURS.filter((c) => r.colours.lenses[c.id]).map((c) => el('div', { class: 'bar-row' },
          el('span', { class: 'bar-name', text: c.name }),
          el('span', { class: 'bar-track' }, el('span', { class: 'bar', style: `width:${(100 * r.colours.lenses[c.id]) / maxLenses}%;background:${c.hex}` })),
          el('span', { class: 'bar-value', text: String(r.colours.lenses[c.id]) }))))),
      el('div', {},
        el('h3', { text: 'Beings who carry several disciplines' }),
        r.colours.carriers.length ? el('ul', { class: 'people' }, ...r.colours.carriers.map((c) => el('li', {},
          el('button', { class: 'link', onclick: () => act.pick({ kind: 'being', id: data.beings[c.being].id }) }, data.beings[c.being].name),
          el('span', { class: 'chips' }, ...Object.entries(c.colours).map(([col, n]) => el('span', { class: 'chip-static' }, swatch(col), ` ${colourOf(col).name} ${n}`)))))) : el('p', { class: 'note', text: 'None.' }),
        el('h3', { text: 'Occasions where disciplines meet' }),
        r.colours.meetings.length ? el('ul', { class: 'people' }, ...r.colours.meetings.map((m) => {
          const o = data.occasions[m.occasion];
          return el('li', {}, el('button', { class: 'link', onclick: () => act.pick({ kind: 'occasion', id: o.id }) }, o.label),
            el('span', { class: 'chips' }, ...m.colours.map((col) => el('span', { class: 'chip-static' }, swatch(col), ` ${colourOf(col).name}`))));
        })) : el('p', { class: 'note', text: 'None.' }))));
}
