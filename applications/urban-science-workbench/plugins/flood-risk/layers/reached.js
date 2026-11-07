// The buildings a flood reaches: the engine holds each of the assessment's modelled buildings against the map chosen; those in over 15 cm
// are summed on H3 cells and raised by their value.

import { chosen, money, onChosen } from '../shared.js';
import { colourMap } from '/visualization/styles/index.js';

export async function layer(ctx) {
  const { gis, frame, dom: { h } } = ctx;
  const id = ctx.layerId();
  const map = colourMap('plum').colour;
  const colour = (u) => map(ctx.theme === 'night' ? 0.2 + 0.8 * u : 1 - (0.15 + 0.8 * u));
  const line = h('div', { class: 'ld-line' }, 'holding every building against the map…');
  let token = 0, dead = false;
  async function draw() {
    const me = ++token;
    const a = await ctx.ask('exposed', { peril: chosen.peril, rp: chosen.rp, threshold: 0.15, res: 10 });
    if (dead || me !== token) return;
    const top = Math.max(1, ...a.value);
    const cells = a.cells.map((c, k) => ({ cell: c, ring: a.rings[k].map(([lon, lat]) => frame.toLocal(lat, lon).slice(0, 2)), count: a.count[k], value: a.value[k] }));
    gis.addPrisms(id, cells, {
      ring: (f) => f.ring, height: (f) => 8 + Math.sqrt(f.value / top) * 240, color: (f) => colour(Math.sqrt(f.value / top)), z: 0.5, opacity: 0.9, shrink: 0.9, order: 3,
      pick: (f) => `${f.count.toLocaleString('en')} buildings in the water here\n${money(f.value)} of value`,
    });
    const share = a.all.count ? (100 * a.reached.count) / a.all.count : 0;
    line.textContent = `${a.reached.count.toLocaleString('en')} buildings of ${a.all.count.toLocaleString('en')} (${share.toFixed(1)} %), ${money(a.reached.value)} of value in over 15 cm`;
  }
  const stop = onChosen(() => draw().catch((e) => { line.textContent = String(e.message || e); }));
  draw().catch((e) => { line.textContent = String(e.message || e); });
  return {
    legend: h('div', { class: 'card' }, h('b', {}, 'The buildings it reaches'), line, h('div', { class: 'ld-quiet' }, 'Modelled buildings of the assessment, on H3 cells of level 10. Value at replacement, in its USD. It is what stands in the water, not the loss.')),
    describe: (hit) => ({ kind: 'cell', title: `${hit.feature.count.toLocaleString('en')} buildings in the water`, lines: [`${money(hit.feature.value)} of value, at a flood of ${chosen.rp} years`, hit.feature.cell] }),
    dispose() { dead = true; stop(); gis.remove(id); },
  };
}
