// How deep each flood reaches: a map of the assessment on the Earth's own surface, under the buildings, changing with the flood chosen.

import { assessment, chosen, depthColour, keys, onChosen } from '../shared.js';

const kept = new Map();
async function bytesOf(file) {
  if (!kept.has(file)) {
    if (kept.size > 3) kept.delete(kept.keys().next().value);
    kept.set(file, fetch(`/data/products/${file}`).then((r) => r.arrayBuffer()).then((b) => new Uint8Array(b)));
  }
  return kept.get(file);
}

export async function layer(ctx) {
  const { earth, dom: { h } } = ctx;
  const P = await assessment();
  const colour = depthColour(ctx.theme);
  const unit = P.depth_code.metres_per_unit;
  const ramp = (u) => { const d = u * 255 * unit; if (d <= 0.05) return [0, 0, 0, 0]; const c = colour(d); return [c[0], c[1], c[2], Math.min(0.94, 0.45 + d)]; };
  let dead = false;
  const line = h('div', { class: 'ld-line' }, '');
  async function draw() {
    const m = P.maps.find((q) => q.peril === chosen.peril && q.return_period === chosen.rp);
    line.textContent = `${P.perils[m.peril]}, once in ${m.return_period} years`;
    const data = await bytesOf(m.file);
    if (dead) return;
    earth.setDrape({ data, width: m.width, height: m.height, box: m.box, ramp, opacity: 0.94, floor: 0.002 });
  }
  const stop = onChosen(draw);
  await draw();
  const K = keys(h, P);
  const stops = Array.from({ length: 12 }, (_, i) => { const c = colour((i / 11) ** 2 * 4); return `rgb(${c.map((v) => Math.round(v * 255)).join(',')})`; }).join(',');
  return {
    legend: h('div', { class: 'card' }, h('b', {}, 'How deep each flood reaches'), K.el, line,
      h('span', { class: 'ramp', style: `background:linear-gradient(90deg, ${stops})` }), h('div', { class: 'ramp-ends' }, h('span', {}, '0'), h('span', {}, '1 m'), h('span', {}, '4 m')),
      h('div', { class: 'ld-quiet' }, P.credit)),
    dispose() { dead = true; stop(); K.stop(); earth.setDrape(null); },
  };
}
