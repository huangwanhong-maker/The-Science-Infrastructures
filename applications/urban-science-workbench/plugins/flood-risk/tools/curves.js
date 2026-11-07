// What a flood costs the city: the assessment's annual loss and exceedance curves, its own figures beside what is recomputed here from its
// event loss table, and which of its wards and buildings' uses carry it.

import { assessment, money } from '../shared.js';

export async function render(root, ctx) {
  const { chart, finding, dom: { h } } = ctx;
  const [S, P] = await Promise.all([ctx.ask('summary'), assessment()]);
  const pc = (a, b) => `${((100 * a) / b).toFixed(1)} %`;
  const exact = (v) => `$${v.toLocaleString('en', { minimumFractionDigits: 2, maximumFractionDigits: 2 })}`;
  root.append(h('h1', {}, 'What a flood costs the city'),
    h('p', { class: 'muted' }, `On average ${money(S.reported)} a year in Freetown, over ${S.events.toLocaleString('en')} modelled flood events in ${S.catalogue_years.toLocaleString('en')} simulated years. Most years pass with almost nothing, and a rare one costs a hundred times the average: the average hides what the curve shows.`));
  const x = S.return_period.map((t) => Math.log10(t));
  const series = [
    { name: 'the largest event of a year', x, y: S.oep.map((v) => v / 1e6), color: chart.category(0) },
    { name: 'all events of a year', x, y: S.aep.map((v) => v / 1e6), color: chart.category(1) },
    { name: 'read here from the events', x, y: S.read_from_events.map((v) => v / 1e6), color: chart.category(2), dash: '3 3' },
  ];
  root.append(chart.lines(series, { xlabel: 'once in so many years', ylabel: 'loss, millions of USD', xfmt: (t) => String(Math.round(10 ** t)), width: 460, height: 230 }));
  root.append(chart.legend(series));
  root.append(finding({
    title: 'The annual loss, recomputed', kind: 'measured',
    statement: `The sum of the event losses over the catalogue's ${S.catalogue_years.toLocaleString('en')} years gives ${exact(S.from_events)}; the assessment reports ${exact(S.reported)}.`,
    method: 'tyche::descriptive::average_annual_loss on the assessment\'s event loss table. The catalogue\'s length is inferred from the table itself: the sum over the reported average is 10,000.0.',
    against: 'The assessment\'s own aalcalc table.',
  }));
  const A = S.area_aep, O = S.area_oep;
  root.append(finding({
    title: 'And from its curves alone', kind: 'measured',
    statement: `Integrating the published curves over annual probability gives ${money(O.total)} from the largest event of a year and ${money(A.total)} from all events of a year, against ${money(S.reported)}. Of the latter, ${money(A.within)} lies within the curve's points, ${money(A.below)} below its first and ${money(A.beyond)} beyond its last, where a loss held at the last point is an assumption.`,
    method: 'tyche::descriptive::curve_area, by trapezoids in annual exceedance probability, loss falling to nothing at a probability of one.',
    limits: 'An area from a dozen points is only as good as they are: the tail beyond 1 in 1,000 years is what the curve cannot say. The event table is what to trust; this shows how far a published curve alone goes, as one often has only that.',
  }));
  const rows = Object.entries(S.by_occupancy).sort((a, b) => b[1] - a[1]);
  root.append(h('h2', {}, 'Who carries it'), h('table', { class: 'sources' }, h('tbody', {}, rows.map(([k, v]) => h('tr', {}, h('th', {}, k), h('td', {}, money(v)), h('td', {}, pc(v, S.reported)))))));
  const wards = [...P.wards].sort((a, b) => (b.aal || 0) - (a.aal || 0)).slice(0, 8);
  root.append(h('h2', {}, 'Which wards'), h('table', { class: 'sources' }, h('tbody', {}, wards.map((w) => h('tr', {}, h('th', {}, w.name), h('td', {}, money(w.aal || 0)))))));
  root.append(h('p', { class: 'note' }, S.credit, ' ', h('a', { href: '/data/products/rdls-freetown.json', target: '_blank', rel: 'noopener' }, 'Its record in the Risk Data Library Standard.')));
}
