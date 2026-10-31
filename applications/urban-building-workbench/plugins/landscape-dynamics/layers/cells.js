// The district as cells: the walking network of the Hoddle Grid's middle as a regular cell complex — junctions,
// segments of streets, lanes and malls, and the blocks they enclose — drawn on the map, each cell picked tells what
// it is; the two stores and where a lane between them would run. Once the lane is built (the analysis "The
// landscape and its events"), the cells it rewrote are drawn warm.

import { ensureDistrict, on, paintCells, state } from '../shared.js';

export async function layer(ctx) {
  const { dom, gis } = ctx;
  const { h } = dom;
  const ids = { blocks: ctx.layerId(':blocks'), edges: ctx.layerId(':edges'), vertices: ctx.layerId(':vertices'), stores: ctx.layerId(':stores'), lane: ctx.layerId(':lane') };
  const status = h('div', { class: 'faint' }, 'cutting the district from the city model…');
  const legend = h('div', { class: 'card' }, h('b', {}, 'The district as cells'), status);
  state.active.add('cells');
  const draw = () => {
    const D = state.district;
    if (!D) return;
    for (const id of Object.values(ids)) gis.remove(id);
    const cells = state.lane ? state.lane.cells : D.cells;
    const created = new Set(state.lane ? state.lane.steps.flatMap((s) => s.created) : []);
    paintCells(ctx, ids, cells, { highlight: created, lane: state.lane ? null : D.lane, stores: state.lane ? state.lane.stores : D.stores });
    const t = state.lane ? state.lane.after : D.topology;
    const [V, E, F] = t.counts;
    status.replaceChildren(
      h('div', {}, `${V} junctions, ${E} segments, ${F} blocks: Euler characteristic ${t.euler}`),
      h('div', {}, `${t.betti[0] === 1 ? 'connected' : `${t.betti[0]} pieces within the box`} (b₀ = ${t.betti[0]}); ${t.holes} block${t.holes === 1 ? '' : 's'} with a dead end left open, so b₁ = ${t.betti[1]}`),
      D.dropped && D.dropped.length ? h('div', { class: 'faint' }, `left out, crossing another way without a junction: ${[...new Set(D.dropped)].slice(0, 4).join(', ')}${D.dropped.length > 4 ? '…' : ''}`) : null,
      h('div', { class: 'faint' }, state.lane ? `the lane built: ${state.lane.steps.length} rewrites, the new cells warm` : D.lane && !D.lane.refused ? `a lane would run ${Math.round(D.lane.length)} m through the block between ${D.stores[0].name} and ${D.stores[1].name}` : D.lane ? `no lane: ${D.lane.refused}` : ''),
      h('div', { class: 'faint' }, 'OpenStreetMap\'s ways; a regular cell complex (morphe)'));
  };
  const off = on((what) => { if (what === 'district' || what === 'lane') draw(); });
  try {
    await ensureDistrict(ctx);
    draw();
  } catch (e) {
    status.textContent = String(e.message || e);
  }
  return {
    legend,
    describe(hit) {
      const f = hit.feature;
      if (hit.layer === ids.blocks) return { title: 'a block (a 2-cell)', lines: [`#${f.id}: ${f.edges} streets round it`, `${Math.round(f.area).toLocaleString('en')} m²`] };
      if (hit.layer === ids.edges) return { title: f.name || 'a segment (a 1-cell)', lines: [`#${f.id}, ${f.type}: ${Math.round(f.length)} m from junction #${f.tail} to #${f.head}`, f.blocks === 2 ? 'between two blocks' : f.blocks === 1 ? 'on the edge of one block' : 'bounding no block'] };
      if (hit.layer === ids.vertices) return { title: `a ${f.type} (a 0-cell)`, lines: [`#${f.id}, on ${f.degree} segments`] };
      if (hit.layer === ids.stores) return { title: f.name, lines: [f.use || 'a store', f.street ? `facing ${f.street}` : ''] };
      return null;
    },
    dispose() { off(); state.active.delete('cells'); },
  };
}
