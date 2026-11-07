// Where the walkers go, against their habit: on H3 cells holding Melbourne's counters, each cell's share of the day's walkers (the
// column's height) and how far that departs from the same weekday over the twelve weeks before (its colour), following the clock's date.
// The Urban Space Usage Index of eca-resilience, with counters where it has GPS pings: these are a sample of the centre's footpaths,
// not of the city, and the shares are of the counters reporting.

import { colourMap } from '/visualization/styles/index.js';

export async function layer(ctx) {
  const { gis, frame, clock, dom: { h } } = ctx;
  const id = ctx.layerId();
  const diverge = colourMap('sea-brick').colour;
  let last = null, dead = false, timer = null, data = null;
  const line = h('div', { class: 'ld-line' }, 'reading every counter…');
  const note = h('div', { class: 'ld-quiet' }, 'Height: the cell\'s share of the day\'s walkers at the counters. Colour: its departure from the same weekday over the twelve weeks before, in standard deviations.');
  const ramp = h('span', { class: 'ramp', style: `background:linear-gradient(90deg, ${Array.from({ length: 11 }, (_, i) => { const c = diverge(i / 10); return `rgb(${c.map((v) => Math.round(v * 255)).join(',')})`; }).join(',')})` });
  const ends = h('div', { class: 'ramp-ends' }, h('span', {}, '−3 σ'), h('span', {}, 'as usual'), h('span', {}, '+3 σ'));
  const source = h('div', { class: 'ld-quiet' }, 'City of Melbourne, pedestrian counting system; the index after eca-resilience');

  async function draw(date) {
    const a = await ctx.ask('usage', { date, weeks: 12, res: 9 });
    if (dead) return;
    data = a;
    const top = Math.max(...a.share);
    const cells = a.cells.map((c, k) => ({ cell: c, ring: a.rings[k].map(([lon, lat]) => frame.toLocal(lat, lon).slice(0, 2)), share: a.share[k], z: a.z[k], walkers: a.walkers[k], mean: a.mean[k], names: a.names[k] }));
    gis.addPrisms(id, cells, {
      ring: (f) => f.ring, height: (f) => 12 + (f.share / top) * 260, color: (f) => (f.z === null ? [0.5, 0.45, 0.55] : diverge(Math.min(1, Math.max(0, 0.5 + f.z / 6)))), z: 0.5, opacity: 0.92, shrink: 0.92, order: 3,
      pick: (f) => `${f.names.slice(0, 2).join(', ')}${f.names.length > 2 ? ` and ${f.names.length - 2} more` : ''}\n${(f.share * 100).toFixed(1)} % of the day's walkers (usually ${(f.mean * 100).toFixed(1)} %)${f.z === null ? '' : `\n${f.z >= 0 ? '+' : '−'}${Math.abs(f.z).toFixed(1)} σ`}`,
    });
    const unusual = cells.filter((c) => c.z !== null && Math.abs(c.z) >= 2).length;
    line.textContent = `${a.date}${a.moved ? ' (the nearest date with a record)' : ''}: ${a.panel} counters on ${a.cells.length} cells; ${unusual} more than 2 σ from their habit`;
  }
  const setTime = (c) => {
    if (c.date === last) return;
    last = c.date;
    clearTimeout(timer);
    timer = setTimeout(() => draw(c.date).catch((e) => { line.textContent = String(e.message || e); }), 180);
  };
  setTime(clock);
  return {
    legend: h('div', { class: 'card' }, h('b', {}, 'Where the walkers go'), line, ramp, ends, note, source),
    setTime,
    describe: (hit) => ({ kind: 'cell', title: hit.feature.names[0], lines: [`${(hit.feature.share * 100).toFixed(1)} % of the day's walkers at the counters, usually ${(hit.feature.mean * 100).toFixed(1)} %`, `${Math.round(hit.feature.walkers).toLocaleString('en')} walkers on ${data.date}`, hit.feature.z === null ? 'no spread in its baseline' : `${hit.feature.z >= 0 ? '+' : '−'}${Math.abs(hit.feature.z).toFixed(1)} standard deviations from its habit`, hit.feature.cell] }),
    dispose() { dead = true; clearTimeout(timer); gis.remove(id); },
  };
}
