// Shops and services along the named streets of the shown city (OpenStreetMap), coloured by use. A holding
// open at the clock's hour stands tall and full; a closed one low and pale; one that declares no hours in
// between. Opening is read from the declared intervals of the civil week, as the domain's street.hpp reads them.

export function openAt(hours, theta) {
  const t = ((theta % 168) + 168) % 168;
  return hours.some(([a, b]) => (t >= a && t < b) || (b > 168 && t + 168 >= a && t + 168 < b));
}

export async function layer(ctx) {
  const { gis, frame, ask, dom, G, area } = ctx;
  const { h } = dom;
  const streets = (await ask('streets')).filter((s) => s.area === area.id || (!s.area && s.city && area.name.includes(s.city)));
  const details = await Promise.all(streets.map((s) => ask('street', { street: s.id })));
  const uses = details[0] ? details[0].uses : [];
  const features = [];
  for (const d of details)
    for (const x of d.holdings)
      if (Number.isFinite(x.lat) && Number.isFinite(x.lon)) features.push({ ...x, street: d.id, streetName: d.name, xy: frame.toLocal(x.lat, x.lon).slice(0, 2) });
  const id = ctx.layerId();
  const cat = G.categorical();
  const pale = (c) => c.map((v) => v + (1 - v) * 0.62);
  let theta = ctx.clock.weekHour;
  const state = (f) => (!f.declared ? 'unknown' : openAt(f.hours, theta) ? 'open' : 'closed');
  gis.addPoints(id, features, {
    position: (f) => f.xy,
    radius: () => 2.4,
    height: (f) => ({ open: 11, closed: 2.2, unknown: 5 })[state(f)],
    color: (f) => (state(f) === 'closed' ? pale(cat(f.use)) : cat(f.use)),
    shadow: true,
    pick: (f) => `${f.name || '(unnamed)'} — ${uses[f.use]}\n${state(f) === 'unknown' ? 'no hours declared' : state(f)}${f.raw ? `: ${f.raw}` : ''}`,
  });
  const present = uses.map((_, i) => features.some((f) => f.use === i));
  const css = (rgb) => `rgb(${rgb.map((c) => Math.round(c * 255)).join(',')})`;
  const count = h('div', { class: 'faint' });
  const legend = h('div', { class: 'card' }, h('b', {}, 'Shops and services'),
    h('div', { class: 'legend' }, uses.map((u, i) => (present[i] ? h('span', {}, h('i', { class: 'dot', style: `background:${css(cat(i))}` }), u) : null))),
    count);
  const tally = () => {
    const open = features.filter((f) => state(f) === 'open').length, known = features.filter((f) => f.declared).length;
    count.textContent = `${open} of ${known} that declare hours are open now (tall); ${features.length - known} declare none; OpenStreetMap`;
  };
  tally();
  return {
    legend,
    setTime(clock) { theta = clock.weekHour; gis.restyle(id, {}); tally(); },
    describe(hit) {
      const f = hit.feature;
      const sel = { kind: 'holding', id: f.id, street: f.street, s: f.s, use: uses[f.use], hours: f.hours, name: f.name || '(unnamed)', title: f.name || '(unnamed)',
        lines: [`${uses[f.use]} · ${f.kind}`, f.raw ? `hours: ${f.raw}` : 'no hours declared', `on ${f.streetName}, ${Math.round(f.s)} m along it`] };
      sel.actions = [
        { title: 'the street\'s hours', run: () => ctx.openTool('hours', sel) },
        { title: 'pair it with another store', run: () => ctx.openTool('reaction', sel) },
        { title: 'the street\'s mix of uses', run: () => ctx.openTool('terrain', sel) },
      ];
      return sel;
    },
  };
}
