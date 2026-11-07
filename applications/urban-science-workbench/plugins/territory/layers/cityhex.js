// A city's buildings on the same grid as the world's: what they enclose, summed on H3 cells of a chosen resolution by the engine, and
// raised as columns at a fifth of the buildings' scale, so that a block can be weighed against its neighbours. It takes the
// buildings in view, whichever way they came: a model's, or the world's tiles.

import { card, compact, tone } from '../shared.js';

const RINGS = [[8, '0.7 km²'], [9, '0.1 km²'], [10, '1.5 ha'], [11, '0.2 ha']];

// Every building in the scene as [lat, lon, volume].
function buildings(ctx) {
  const out = [];
  const { gis, frame } = ctx;
  const m = gis.layers.get('buildings');
  if (m) for (const b of m.features) {
    const ring = b.footprint || [];
    if (ring.length < 3) continue;
    let a = 0, cx = 0, cy = 0;
    for (let i = 0; i < ring.length; i++) { const [x0, y0] = ring[i], [x1, y1] = ring[(i + 1) % ring.length]; a += x0 * y1 - x1 * y0; cx += x0, cy += y0; }
    const g = frame.toGeodetic(cx / ring.length, cy / ring.length);
    out.push([g.lat, g.lon, (Math.abs(a) / 2) * (b.height_m || 9)]);
  }
  for (const [, l] of gis.layers) {
    if (l.kind !== 'built-solids') continue;
    const { boxes, rings, starts } = l.built, n = boxes.length / 5;
    for (let k = 0; k < n; k++) {
      const s = starts[2 * k], c = starts[2 * k + 1];
      let a = 0;
      for (let i = 0; i < c; i++) { const j = (i + 1) % c; a += rings[2 * (s + i)] * rings[2 * (s + j) + 1] - rings[2 * (s + j)] * rings[2 * (s + i) + 1]; }
      const g = frame.toGeodetic(boxes[5 * k], boxes[5 * k + 1]);
      out.push([g.lat, g.lon, (Math.abs(a) / 2) * (boxes[5 * k + 4] - boxes[5 * k + 3])]);
    }
  }
  return out;
}

export async function layer(ctx) {
  const { gis, frame, dom: { h } } = ctx;
  const id = ctx.layerId();
  const colour = tone('plum', ctx.theme, { from: 0.18 });
  let res = 10, dead = false, cells = [], seen = -1;
  const status = h('div', { class: 'ld-quiet' }, 'binning the buildings in view…');
  const ends = h('div', { class: 'ramp-ends' }, h('span', {}, '0'), h('span', {}, ''));

  async function draw(force) {
    if (dead) return;
    const points = buildings(ctx);
    if (!points.length) { status.textContent = 'no buildings in view yet'; return; }
    if (!force && points.length === seen) return;
    seen = points.length;
    const answer = await ctx.ask('bin_points', { res, points: points.map((p) => [+p[0].toFixed(6), +p[1].toFixed(6), Math.round(p[2])]) });
    if (dead) return;
    cells = answer.cells.map((c, k) => ({ cell: c, ring: answer.rings[k].map(([lon, lat]) => frame.toLocal(lat, lon).slice(0, 2)), volume: answer.weight[k], count: answer.count[k], area: answer.area_m2[k] }));
    const sorted = cells.map((c) => c.volume / c.area).sort((a, b) => a - b);
    const top = sorted[Math.floor(sorted.length * 0.99)] || 1;
    const size = Math.sqrt(cells[0].area);
    gis.addPrisms(id, cells, {
      ring: (f) => f.ring, height: (f) => Math.max(1, Math.min(1, f.volume / f.area / top) * size * 0.9), color: (f) => colour(Math.min(1, f.volume / f.area / top)), z: 0.4, opacity: 0.86, shrink: 0.9, order: 3,
      pick: (f) => `${compact(f.volume)} m³ in ${f.count} buildings\nas if ${(f.volume / f.area).toFixed(1)} m deep over the cell`,
    });
    ends.children[1].textContent = `${top.toFixed(1)} m³ per m²`;
    status.textContent = `${cells.length.toLocaleString('en')} cells of H3 level ${res}, from ${points.length.toLocaleString('en')} buildings`;
  }
  const keys = h('div', { class: 'keys-tight', title: 'the size of a cell' }, RINGS.map(([r, name]) => h('button', { class: r === res ? 'on' : '', onclick: (e) => {
    res = r;
    for (const b of e.currentTarget.parentElement.children) b.classList.toggle('on', b === e.currentTarget);
    draw(true);
  } }, name)));
  const timer = setInterval(() => draw(false).catch((e) => { status.textContent = String(e.message || e); }), 4000);
  draw(true).catch((e) => { status.textContent = String(e.message || e); });
  const legend = card(h, 'Built volume by hexagon', { colour, extra: [keys], note: 'What the buildings enclose, as a depth over the cell; a guess where a height is estimated.' });
  legend.insertBefore(ends, legend.querySelector('.ld-quiet'));
  legend.append(status);
  return {
    legend,
    describe: (hit) => ({ kind: 'cell', title: `A cell of level ${res}`, lines: [`${compact(hit.feature.volume)} m³ in ${hit.feature.count} buildings`, `as if ${(hit.feature.volume / hit.feature.area).toFixed(1)} m deep over the cell`, hit.feature.cell] }),
    dispose() { dead = true; clearInterval(timer); },
  };
}
