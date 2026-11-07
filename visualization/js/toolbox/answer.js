// A tool's answer, read: the blocks the engine sends (api/cpp/science/toolbox.hpp) as figures, tables, charts and notes. Values for
// places (the kinds "map", "rows" and "each") are for the application to carry to its map: it says how, by `places(block)`, and what
// it answers stands in the block's stead. Styled by each application's stylesheet (.tb-block, .figures, table.results, figure.chart).

import { el, csv, download } from './dom.js';
import { barChart, formatNumber, lineChart, matrixChart, scatterChart } from './charts.js';

const slug = (s) => s.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');

export function blockView(b, { name = 'answer', places = () => null } = {}) {
  switch (b.kind) {
    case 'figures':
      return el('div', { class: 'tb-block' }, el('div', { class: 'heading' }, el('span', { text: b.title })),
        el('div', { class: 'figures' }, b.items.map(([k, v, n]) => el('div', { class: 'figure-cell', title: n || '' }, el('div', { class: 'value', text: typeof v === 'number' ? formatNumber(v) : v ?? '—' }),
          el('div', { class: 'label', text: k })))));
    case 'table': {
      const copy = el('button', { class: 'note', text: 'copy', onclick: () => { navigator.clipboard?.writeText(csv(b.columns, b.rows)); copy.textContent = 'copied'; } });
      const save = el('button', { class: 'note', text: 'save CSV', onclick: () => download(`${name}-${slug(b.title)}.csv`, csv(b.columns, b.rows)) });
      return el('div', { class: 'tb-block' }, el('div', { class: 'heading' }, el('span', { text: b.title }), copy, save),
        el('div', { class: 'table-wrap' }, el('table', { class: 'results' },
          el('thead', {}, el('tr', {}, b.columns.map((c, i) => el('th', { class: i ? 'num' : '', text: c })))),
          el('tbody', {}, b.rows.map((r) => el('tr', {}, r.map((v, i) => el('td', { class: i && typeof v !== 'string' ? 'num' : '', text: typeof v === 'number' ? formatNumber(v) : v ?? '—' }))))))),
        b.note ? el('p', { class: 'note', text: b.note }) : null);
    }
    case 'lines': return el('div', { class: 'tb-block' }, lineChart(b));
    case 'bars': return el('div', { class: 'tb-block' }, barChart(b));
    case 'scatter': return el('div', { class: 'tb-block' }, scatterChart(b));
    case 'matrix': return el('div', { class: 'tb-block' }, matrixChart(b));
    case 'note': return el('p', { class: 'note tb-note', text: b.text });
    case 'map': case 'rows': case 'each': return places(b);
    default: return null;
  }
}
