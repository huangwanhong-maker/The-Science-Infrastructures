// The district's events in spacetime (ADR 0067; the design, §12): the record's event complex drawn with time as
// height above the district, `scale` metres a year from the window's start.
//   The couplings (volumes): each building, street space and junction space a worldline rising from its anchor over
//   the years it stands, its vertices (the events and censuses it attends) marked along it.
//   The primitives (interfaces): each party wall, frontage and end of a street a sheet between the worldlines of
//   its two sides, cut into lenses at the events that meet it: the faces of the 2-complex, whose boundary runs up
//   one side and down the other.
//   The events: in the refined complex K̂ a bead at each coupling attending, joined by rungs along the interfaces
//   it meets; in the event complex K one bead at the event, its couplings' worldlines bending into it.
//   The moment: a plane at the clock's date, carrying the network of what stood then, nodes and links.
// The city beneath can be seen through: lower its opacity in the rail.

import { cellColour, css, ensureRecord, eventColour, EVENT_TYPES, kindOfCell, on, setView, state, z, yearOf } from '../shared.js';

export async function layer(ctx) {
  const { dom, gis } = ctx;
  const { h } = dom;
  const night = ctx.theme === 'night';
  const ids = {
    lines: ctx.layerId(':worldlines'), ticks: ctx.layerId(':vertices'), faces: ctx.layerId(':faces'), events: ctx.layerId(':events'),
    rungs: ctx.layerId(':rungs'), plane: ctx.layerId(':plane'), edge: ctx.layerId(':plane-edge'), nodes: ctx.layerId(':now-nodes'), links: ctx.layerId(':now-links'),
  };
  const status = h('div', { class: 'ld-quiet' }, 'building the record…');
  const nowLine = h('div', { class: 'ld-line' });
  const sw = (color, label, kind) => h('span', {}, h('i', { class: `ld-sw ${kind}`, style: `background:${color}` }), label);
  const legend = h('div', { class: 'card' }, h('b', {}, 'The district\'s events in spacetime'),
    h('div', { class: 'ld-quiet' }, 'volumes rise as worldlines, the interfaces between them as sheets, events as beads'),
    h('div', { class: 'ld-swatches' }, ...['building', 'street', 'junction'].map((t) => sw(css(cellColour(t, night)), t, 'line')),
      ...['party-wall', 'frontage'].map((t) => sw(css(cellColour(t, night)), t.replace('-', ' '), 'block'))),
    h('div', { class: 'ld-swatches' }, ...EVENT_TYPES.slice(0, 4).map((e) => sw(e.hue, e.name, 'dot'))),
    nowLine, status);
  let t = null;  // the clock's MJD

  const drawn = () => state.view;
  function build() {
    const R = state.record;
    if (!R) return;
    const V = drawn();
    const shown = (type) => V.show[kindOfCell(type)] !== false;
    const C = R.couplings;
    // Where each coupling stands at each event it attends in K: at the event's bead; elsewhere at its anchor.
    const bends = C.map(() => []);
    if (V.mode === 'events')
      for (const e of R.events) for (const c of e.attend) bends[c].push([e.from, e.at]);
    const lines = C.map((c, i) => ({ c, i })).filter(({ c }) => shown(c.type));
    gis.addTubes(ids.lines, lines, {
      path: ({ c, i }) => {
        const [x, y] = c.at, pts = [[x, y, z(c.life[0])]];
        for (const [when, at] of bends[i].sort((a, b) => a[0] - b[0])) pts.push([at[0], at[1], z(when)]);
        pts.push([x, y, z(c.life[1])]);
        return pts;
      },
      radius: ({ c }) => (c.type === 'building' ? 1.2 : c.type === 'junction' ? 0.35 : 0.5),
      color: ({ c }) => cellColour(c.type, night), sides: 6, opacity: 0.9,
    });
    // The censuses along the buildings' worldlines, each a small disc: the local ticks that observe them.
    const ticks = [];
    if (V.mode === 'refined' && V.censuses)
      for (const { c } of lines) for (const y of c.censuses || []) ticks.push({ c, when: Date.UTC(y, 0, 1) / 86400000 + 40587 });
    gis.addPoints(ids.ticks, ticks, { position: (f) => f.c.at, base: (f) => z(f.when) - 0.15, height: () => 0.3, radius: () => 1.9,
      color: (f) => cellColour(f.c.type, night) });
    // The faces: each interface's lenses, a thin sheet between its two worldlines.
    const faces = [];
    for (const [k, p] of R.primitives.entries()) {
      if (!shown(p.type) || !shown(C[p.a].type) || !shown(C[p.b].type)) continue;
      for (let m = 0; m + 1 < p.meets.length; m++) faces.push({ p, k, from: p.meets[m], to: p.meets[m + 1] });
    }
    gis.addSlabs(ids.faces, faces, {
      line: (f) => [C[f.p.a].at, C[f.p.b].at], width: () => 0.35, base: (f) => z(f.from), height: (f) => Math.max(0.05, z(f.to) - z(f.from)),
      color: (f) => cellColour(f.p.type, night), opacity: 0.2, order: 3,
    });
    // The events.
    const beads = [], rungs = [];
    for (const e of R.events) {
      if (V.mode === 'events') beads.push({ e, at: e.at, r: 5 });
      else {
        for (const c of e.attend) if (shown(C[c].type)) beads.push({ e, at: C[c].at, r: e.type === 'regime' ? 1.5 : 3.4, flat: e.type === 'regime' });
        for (const k of e.meets) {
          const p = R.primitives[k];
          if (shown(p.type)) rungs.push({ e, p });
        }
      }
    }
    gis.addPoints(ids.events, beads, { position: (b) => b.at, base: (b) => z(b.e.from) - (b.flat ? 0.2 : b.r / 2), height: (b) => (b.flat ? 0.4 : b.r), radius: (b) => b.r,
      color: (b) => eventColour(b.e.type) });
    gis.addTubes(ids.rungs, rungs, { path: ({ e, p }) => [[...C[p.a].at, z(e.from)], [...C[p.b].at, z(e.from)]], radius: () => 0.45,
      color: ({ e }) => eventColour(e.type), opacity: 0.9 });
    const B = state.box;
    gis.addFlats(ids.plane, [{ ring: [[B.x0, B.y0], [B.x1, B.y0], [B.x1, B.y1], [B.x0, B.y1]] }], {
      ring: (f) => f.ring, z: 0, color: () => (night ? [0.36, 0.52, 0.66] : [0.36, 0.56, 0.72]), opacity: 0.16, order: 6,
    });
    gis.addTubes(ids.edge, [{ ring: [[B.x0, B.y0], [B.x1, B.y0], [B.x1, B.y1], [B.x0, B.y1], [B.x0, B.y0]] }], {
      path: (f) => f.ring.map(([x, y]) => [x, y, 0]), radius: () => 0.9, color: () => (night ? [0.55, 0.72, 0.86] : [0.25, 0.48, 0.66]),
    });
    lift();
  }

  // The plane and the network of the moment: the couplings and interfaces whose lives hold the clock's date.
  function lift() {
    const R = state.record;
    if (!R || t === null) return;
    const V = drawn();
    const T = Math.min(R.window.to - 1e-3, Math.max(R.window.from, t));
    for (const id of [ids.plane, ids.edge]) { const l = gis.layers.get(id); if (l) l.object.position.y = z(T); }
    const alive = (life) => life[0] <= T && T < life[1];
    const C = R.couplings;
    const shown = (type) => V.show[kindOfCell(type)] !== false;
    const nodes = C.filter((c) => alive(c.life) && shown(c.type));
    const links = R.primitives.filter((p) => alive(p.life) && shown(p.type) && shown(C[p.a].type) && shown(C[p.b].type));
    gis.addPoints(ids.nodes, nodes, { position: (c) => c.at, base: () => z(T) - 1, height: () => 2, radius: (c) => (c.type === 'building' ? 2.4 : 1.6),
      color: (c) => cellColour(c.type, night) });
    gis.addTubes(ids.links, links, { path: (p) => [[...C[p.a].at, z(T)], [...C[p.b].at, z(T)]], radius: () => 0.35, color: (p) => cellColour(p.type, night) });
    const done = R.events.filter((e) => e.from <= T).length;
    nowLine.textContent = `${Math.floor(yearOf(T))}: ${nodes.length} volumes and ${links.length} interfaces standing; ${done} of ${R.events.length} events so far`;
  }

  // The view: from the side and a little above, at the middle of the record, from where the map was seen.
  function look() {
    const R = state.record, B = state.box;
    const mid = z((R.window.from + R.window.to) / 2), T = gis.toWorld((B.x0 + B.x1) / 2, (B.y0 + B.y1) / 2, mid);
    const o = gis.camera.position.clone().sub(gis.controls.target);
    const heading = Math.atan2(o.x, o.z), tilt = (66 * Math.PI) / 180, d = 1.25 * Math.max(B.x1 - B.x0, B.y1 - B.y0, 2 * mid);
    gis.controls.target.copy(T);
    gis.camera.position.set(T.x + d * Math.sin(tilt) * Math.sin(heading), T.y + d * Math.cos(tilt), T.z + d * Math.sin(tilt) * Math.cos(heading));
    gis.controls.update();
  }

  const summary = () => {
    const R = state.record, K = R.topology.K;
    const years = `${R.window.fromDate.slice(0, 4)}–${R.window.toDate.slice(0, 4)}`;
    return `${years} rising ${state.view.scale} m a year: ${R.network.buildings} buildings, ${R.network.streets} street spaces, ${R.network.junctions} junctions; ` +
      `${R.events.length} events and ${R.censuses.reduce((a, c) => a + c.buildings, 0)} census ticks; K has ${K.vertices} vertices, ${K.edges} edges, ${K.faces} faces`;
  };
  const off = on((what) => {
    if (what === 'view') build();
    if (what === 'error') status.textContent = state.error;
  });
  try {
    await ensureRecord(ctx);
    build();
    status.textContent = summary();
    look();
    // A foam is easier to read through a city half seen.
    if (ctx.cityOpacity && ctx.cityOpacity() > 0.5) ctx.setCityOpacity(0.35);
  } catch (e) {
    status.textContent = String(e.message || e);
  }
  return {
    legend,
    setTime(c) {
      t = c.instant / 86400000 + 40587;
      lift();
    },
    describe(hit) {
      const R = state.record;
      if (!R) return null;
      const C = R.couplings;
      const name = (c) => `${c.type}${c.name ? ` ${c.name}` : ''}`;
      if (hit.layer === ids.events || hit.layer === ids.rungs) {
        const e = hit.feature.e;
        const meets = {};
        for (const k of e.meets) meets[R.primitives[k].type] = (meets[R.primitives[k].type] || 0) + 1;
        return {
          kind: 'event', title: `${e.date} — ${e.kind || e.type}`,
          lines: [
            `${(EVENT_TYPES.find((x) => x.type === e.type) || { name: e.type }).name}; ${e.provenance}, from ${e.source || 'the record'}`,
            `attended by ${e.attend.length}: ${e.attend.slice(0, 4).map((c) => name(C[c])).join(', ')}${e.attend.length > 4 ? ', …' : ''}`,
            Object.keys(meets).length ? `meets ${Object.entries(meets).map(([k, n]) => `${n} ${k.replace('-', ' ')}${n === 1 ? '' : 's'}`).join(', ')}` : 'meets no interface: it cuts its worldlines alone',
          ],
        };
      }
      if (hit.layer === ids.lines) {
        const c = hit.feature.c;
        return { kind: 'coupling', title: name(c), lines: [`a coupling (a volume) from ${c.life[0] <= R.window.from + 1 ? 'before the record' : new Date((c.life[0] - 40587) * 86400000).getUTCFullYear()} to ${c.life[1] >= R.window.to - 1 ? 'its end' : new Date((c.life[1] - 40587) * 86400000).getUTCFullYear()}`,
          `${c.vertices - 2} events and censuses along its worldline`] };
      }
      if (hit.layer === ids.faces) {
        const f = hit.feature;
        return { kind: 'primitive', title: `${f.p.type.replace('-', ' ')} between ${name(C[f.p.a])} and ${name(C[f.p.b])}`,
          lines: [`a lens of the interface, ${Math.floor(yearOf(f.from))}–${Math.floor(yearOf(f.to))}, between two events that meet it`, f.p.area ? `${f.p.area.toLocaleString()} m²` : ''] };
      }
      return null;
    },
    dispose() { off(); },
  };
}

export { setView };
