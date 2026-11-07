// What the flood's layers share: the assessment's product, and which flood is chosen, so that the depth on the ground and the buildings it
// reaches are of the same one.

import { colourMap } from '/visualization/styles/index.js';

let product = null;
export function assessment() {
  if (!product) {
    product = fetch('/data/products/freetown-flood.json').then(async (r) => {
      if (!r.ok) throw new Error('The assessment is not built: python data/products/freetown_flood.py');
      return r.json();
    }).catch((e) => { product = null; throw e; });
  }
  return product;
}

const listeners = new Set();
export const chosen = { peril: 'FLRF', rp: 100 };
export const onChosen = (f) => { listeners.add(f); return () => listeners.delete(f); };
export function choose(changes) {
  Object.assign(chosen, changes);
  for (const f of listeners) f(chosen);
}

// The keys that choose it: the peril, and how rare a flood.
export function keys(h, { perils, return_periods }) {
  const row = (items, key, label) => h('div', { class: 'keys-tight', title: label }, items.map(([value, name]) => {
    const b = h('button', { class: chosen[key] === value ? 'on' : '', onclick: () => choose({ [key]: value }) }, name);
    b.dataset.value = String(value);
    return b;
  }));
  const a = row(Object.entries(perils).map(([k, n]) => [k, n]), 'peril', 'what floods'), c = row(return_periods.map((p) => [p, String(p)]), 'rp', 'once in so many years');
  const sync = () => {
    for (const b of a.children) b.classList.toggle('on', b.dataset.value === chosen.peril);
    for (const b of c.children) b.classList.toggle('on', +b.dataset.value === chosen.rp);
  };
  const stop = onChosen(sync);
  return { el: h('div', { class: 'flood-keys' }, a, h('div', { class: 'ld-quiet' }, 'once in so many years'), c), stop };
}

const seaMap = colourMap('sea').colour;
export const depthColour = (theme) => (d) => {
  const u = Math.min(1, Math.sqrt(d / 4));
  const c = seaMap(theme === 'night' ? 1 - (0.2 + 0.75 * u) : 0.2 + 0.75 * u);
  return c;
};
export const money = (v) => (v >= 1e6 ? `$${(v / 1e6).toFixed(v >= 1e7 ? 0 : 1)} m` : v >= 1e3 ? `$${Math.round(v / 1e3)}k` : `$${Math.round(v)}`);
