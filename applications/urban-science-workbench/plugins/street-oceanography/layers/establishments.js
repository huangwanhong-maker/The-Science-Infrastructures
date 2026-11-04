// Every business about Melbourne's sensors in the City of Melbourne's census of land use (CLUE, its latest
// year), coloured by use, placed by its offset from its sensor.

export async function layer(ctx) {
  const { gis, frame, ask, dom, G } = ctx;
  const { h } = dom;
  const d = await ask('establishments', {});
  const origin = Object.fromEntries(Object.entries(d.sensors).map(([k, [lat, lon]]) => [k, frame.toLocal(lat, lon)]));
  const features = d.points.map(([sid, x, y, use, name]) => ({ sid, use, name, xy: [origin[String(sid)][0] + x, origin[String(sid)][1] + y] }));
  const id = ctx.layerId();
  const cat = G.categorical();
  gis.addPoints(id, features, {
    position: (f) => f.xy, radius: () => 1.6, height: () => 3.5, base: () => 0.2,
    color: (f) => cat(f.use),
    pick: (f) => `${f.name || '(unnamed)'}\n${d.uses[f.use]}`,
  });
  const css = (rgb) => `rgb(${rgb.map((c) => Math.round(c * 255)).join(',')})`;
  const legend = h('div', { class: 'card' }, h('b', {}, `Businesses by use (${d.year})`),
    h('div', { class: 'legend' }, d.uses.map((u, i) => h('span', {}, h('i', { class: 'dot', style: `background:${css(cat(i))}` }), u))),
    h('div', { class: 'faint' }, `${features.length.toLocaleString('en')} establishments within 200 m of the sensors; City of Melbourne, CLUE`));
  return {
    legend,
    describe(hit) {
      const f = hit.feature;
      return { kind: 'establishment', title: f.name || '(unnamed)', lines: [d.uses[f.use], `about sensor ${f.sid}`] };
    },
  };
}
