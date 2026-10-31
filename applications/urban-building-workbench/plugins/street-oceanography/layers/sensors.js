// Pedestrians at Melbourne's sensors: a column at each, its height and colour the mean count at the clock's
// hour of the civil week over 2024–2026 (the engine's profiles of the City of Melbourne's hourly counts).

const PER_PERSON = 0.05;   // metres of column per person an hour

export async function layer(ctx) {
  const { gis, frame, ask, dom, G } = ctx;
  const { h } = dom;
  const list = await ask('profiles', { era: 'recent' });
  const id = ctx.layerId();
  const features = list.map((s) => ({ ...s, xy: frame.toLocal(s.lat, s.lon).slice(0, 2) }));
  const top = Math.max(1, ...features.flatMap((f) => f.profile.filter(Number.isFinite)));
  const colour = G.sequential(0, Math.log1p(top), 'gold', { from: 0.15, to: 1 });
  let theta = ctx.clock.weekHour;
  const at = (f) => {
    const k = Math.floor(theta) % 168, u = theta - Math.floor(theta);
    const a = f.profile[k], b = f.profile[(k + 1) % 168];
    return Number.isFinite(a) && Number.isFinite(b) ? a + u * (b - a) : (Number.isFinite(a) ? a : NaN);
  };
  gis.addPoints(id, features, {
    shape: 'box', shadow: true,
    position: (f) => f.xy,
    radius: () => 6,
    height: (f) => 2 + PER_PERSON * (Number.isFinite(at(f)) ? at(f) : 0),
    color: (f) => colour(Math.log1p(Math.max(0, at(f)))) || [0.7, 0.7, 0.7],
    pick: (f) => `${f.name}\n${Number.isFinite(at(f)) ? Math.round(at(f)).toLocaleString('en') : '—'} people an hour on average at this hour`,
  });
  const legend = h('div', { class: 'card' }, h('b', {}, 'Pedestrians at the sensors'),
    h('div', {}, `columns ${PER_PERSON * 1000} m tall per thousand people an hour; mean of ${Math.round(Math.min(...features.map((f) => f.weeks)))}+ weeks at each hour of the week, 2024–2026 (City of Melbourne)`));
  return {
    legend,
    setTime(clock) { theta = clock.weekHour; gis.restyle(id, {}); },
    describe(hit) {
      const f = hit.feature;
      const sel = { kind: 'sensor', id: f.id, name: f.name, title: f.name,
        lines: [`${Number.isFinite(at(f)) ? Math.round(at(f)).toLocaleString('en') : '—'} people an hour on average at this hour of the week`, `sensor ${f.id}`] };
      sel.actions = [
        { title: 'its rhythm', run: () => ctx.openTool('rhythm', sel) },
        { title: 'weather and holidays', run: () => ctx.openTool('response', sel) },
        { title: 'its mix of uses', run: () => ctx.openTool('terrain', sel) },
      ];
      return sel;
    },
  };
}
