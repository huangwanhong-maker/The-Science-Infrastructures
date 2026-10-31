// A building's rooms and their warmth (a model): the building picked in the analysis "A building's rooms", as a box
// of rooms on its floors fitted to its footprint, each room a 3-cell drawn as a block, coloured by its temperature
// at the clock's hour as the Sun on its facades and roof and the outside air warm it; stairs, doors and merged
// rooms appear from the hour they are opened. The heat is the engine's (urban-science, landscape/building.hpp).

import { on, state, unitRing } from '../shared.js';

export async function layer(ctx) {
  const { dom, gis, G } = ctx;
  const { h, fmt } = dom;
  const ids = { rooms: ctx.layerId(':rooms'), openings: ctx.layerId(':openings') };
  const status = h('div', { class: 'faint' }, 'pick a building and open "A building\'s rooms" to fill it with rooms.');
  const legend = h('div', { class: 'card' }, h('b', {}, 'A building\'s rooms (a model)'), ctx.chart.ramp('gold', 12, 32, 'room temperature, °C'), status);
  state.active.add('rooms');
  const warmth = G.sequential(12, 32, 'gold', { from: 0.05, to: 1 });
  let clock = ctx.clock;
  const show = () => {
    const B = state.building;
    if (!B) return;
    const hour = Math.max(0, Math.min(23, Math.floor(clock.hour)));
    const H = B.hours[hour];
    const spec = B.spec;
    // Each room's unit cubes, one solid each, at its floor.
    const features = [];
    for (const room of H.rooms) for (const u of room.units) features.push({ room, u, ring: unitRing(spec, u) });
    gis.addSolids(ids.rooms, features, { ring: (f) => f.ring, base: (f) => f.u[2] * spec.storey + 0.15, height: (f) => (f.u[2] + 1) * spec.storey - 0.15, color: (f) => warmth(f.room.temperature) || [0.8, 0.8, 0.8],
      pick: (f) => `a room (#${f.room.id}${f.room.units.length > 1 ? `, ${f.room.units.length} rooms merged` : ''}): ${fmt(f.room.temperature, 3)} °C` });
    const opened = new Map(B.events.map((e) => [e.kind === 'stairs' ? 'stair' : e.kind === 'doors' ? 'door' : e.kind, e.hour]));
    const shown = B.openings.filter((o) => o.type === 'entrance' || (opened.has(o.type) && hour >= opened.get(o.type)));
    gis.addSolids(ids.openings, shown, {
      ring: (o) => openingRing(spec, o),
      base: (o) => o.corner[2] * spec.storey - (o.span[2] ? 0 : 0.3), height: (o) => o.corner[2] * spec.storey + (o.span[2] ? spec.storey * 0.8 : 0.3),
      color: (o) => (o.type === 'stair' ? [0.79, 0.55, 0.42] : o.type === 'entrance' ? [0.36, 0.56, 0.73] : [0.6, 0.56, 0.47]), pick: (o) => `${o.type === 'stair' ? 'a stair' : o.type === 'door' ? 'a door' : 'the entrance'}` });
    status.textContent = `${B.name}: ${H.rooms.length} rooms on ${spec.floors} floors at ${String(hour).padStart(2, '0')}:00, the Sun bringing ${fmt(H.sun / 1000, 3)} kW; its slowest mode decays in ${fmt(B.timeConstants[0], 3)} h`;
  };
  const off = on((what) => { if (what === 'building') show(); });
  show();
  return {
    legend,
    setTime(c) { clock = c; show(); },
    dispose() { off(); state.active.delete('rooms'); },
  };
}

// A face of the grid as a thin footprint: a horizontal one (a stair) a square in the middle of its cell; a vertical
// one (a door, the entrance) a slab 0.6 m thick along the middle 40% of its wall.
function openingRing(spec, o) {
  const c = Math.cos(spec.heading), s = Math.sin(spec.heading), dx = spec.width / spec.nx, dy = spec.depth / spec.ny;
  const at = (u, v) => [spec.x + u * c - v * s, spec.y + u * s + v * c];
  const [i, j] = o.corner, [sx, sy, sz] = o.span;
  let u0, u1, v0, v1;
  if (!sz) u0 = (i + 0.3) * dx, u1 = (i + 0.7) * dx, v0 = (j + 0.3) * dy, v1 = (j + 0.7) * dy;
  else if (sx) u0 = (i + 0.3) * dx, u1 = (i + 0.7) * dx, v0 = j * dy - 0.3, v1 = j * dy + 0.3;
  else u0 = i * dx - 0.3, u1 = i * dx + 0.3, v0 = (j + 0.3) * dy, v1 = (j + 0.7) * dy;
  void sy;
  return [at(u0, v0), at(u1, v0), at(u1, v1), at(u0, v1)];
}
