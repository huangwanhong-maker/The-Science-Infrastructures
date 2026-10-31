// Kinds of streets: Melbourne's segments (its pedestrian sensors) placed by how finely their uses change (from
// the land-use census about them) and how selectively they answer the calendars (from their counts), and whether
// they fall into groups. A point picked here flies the map to its sensor.

export async function render(root, ctx) {
  const { h, busy, failure, fmt, css } = ctx.dom;
  const { scatter, category, legend } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Kinds of streets'),
    h('span', { class: 'lede' }, 'If streets differ in kind, they should fall into groups when placed by what they are measured to do.')));
  const out = h('div', {});
  root.append(out);
  busy(out, 'placing every outdoor sensor (once, about a minute)…');
  let d, rows, sensors;
  try { [d, rows, sensors] = await Promise.all([ctx.ask('plane'), ctx.ask('summary'), ctx.ask('sensors')]); } catch (e) { return failure(out, e); }
  const named = (n) => /Macaulay|Errol|Lygon|Bourke Street Mall|Southbank|Southern Cross|Flinders|Chinatown|Swanston|Docklands|QV/.test(n);
  const pts = d.points.map((p, i) => ({ x: p.compression, y: p.selectivity, color: named(p.name) ? css('--accent') : category(i % 10), r: named(p.name) ? 6 : 4.5,
    label: `${p.name}\nuses change ${fmt(p.compression, 3)} times in a 400 m visit; selectivity ${fmt(p.selectivity, 3)}\n(click to fly there)`, name: named(p.name) ? p.name.slice(0, 26) : '', id: p.id }));
  const plot = scatter(pts, { xlabel: 'changes of use a 400 m visit crosses', ylabel: 'selectivity of the rhythm', height: 420, width: 600, names: true });
  plot.querySelectorAll('circle').forEach((c, i) => {
    c.style.cursor = 'pointer';
    c.addEventListener('click', () => {
      const s = sensors.find((x) => String(x.id) === String(pts[i].id));
      if (!s || ctx.area.id !== 'melbourne') return;
      const [x, y] = ctx.frame.toLocal(s.lat, s.lon);
      ctx.gis.flyTo([x, y], 900);
      const sel = { kind: 'sensor', id: s.id, name: s.name, title: s.name, lines: [`sensor ${s.id}`] };
      sel.actions = [{ title: 'its rhythm', run: () => ctx.openTool('rhythm', sel) }, { title: 'weather and holidays', run: () => ctx.openTool('response', sel) }];
      ctx.select(sel);
    });
  });
  const missing = rows.filter((r) => !Number.isFinite(r.compression) || !Number.isFinite(r.selectivity));
  const title = { clustered: 'The streets fall into groups', unstructured: 'No groups: the streets fill the plane', regular: 'The streets spread evenly' }[d.reading] || 'Too few streets placed';
  out.replaceChildren(plot, legend([{ name: 'streets named on the plot', color: css('--accent'), dot: true }]),
    h('div', { class: 'stat-row' }, h('div', { class: 'count' }, h('b', {}, fmt(d.points.length)), h('span', {}, 'streets placed')),
      h('div', { class: 'count' }, h('b', {}, fmt(d.ratio, 3)), h('span', {}, 'nearest-neighbour ratio R')), h('div', { class: 'count' }, h('b', {}, fmt(d.p, 2)), h('span', {}, 'p'))),
    ctx.finding({ title, kind: 'measured',
      statement: `R compares the mean distance from each street to its nearest neighbour in the plane with what scattering at random would give (R < 1: groups; R ≈ 1: no structure; R > 1: an even spread). Here R = ${fmt(d.ratio, 3)} (z = ${fmt(d.z, 3)}).${missing.length ? ` ${missing.length} sensors are not placed (under 300 days of counts in 2024–2026, or under 20 businesses about them).` : ''}`,
      method: 'Clark and Evans\'s test with an edge correction, both axes scaled to their ranges.',
      limits: 'The sensors were placed for counting pedestrians, not chosen to span kinds of streets; the length of a visit (400 m) is declared, and does not change the grouping while it is the same for every street.' }));
}
