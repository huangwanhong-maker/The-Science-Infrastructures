// The card of a part: what is known of it (the anatomy domain's knowledge base, inherited along is-a where the
// part has no entry of its own), where it lies, what the domain measures of it, what it touches, its place
// in FMA's two trees, and for a vessel its calibre along its course.

import { systemColour } from './palette.js';
import { capital } from './tree.js';

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmt = (x, d = 0) => (x === null || x === undefined || !isFinite(x) ? '—' : x.toLocaleString('en-GB', { minimumFractionDigits: d, maximumFractionDigits: d }));
const SYSTEM_WORDS = {
  arteries: 'an artery', veins: 'a vein', nerves: 'a nerve', muscles: 'a muscle', skeleton: 'of the skeleton', heart: 'of the heart',
  brain: 'of the brain and spinal cord', skin: 'of the skin', fascia: 'fascia', eye: 'of the eye', ear: 'of the ear',
  respiratory: 'of the respiratory system', digestive: 'of the digestive system', urinary: 'of the urinary system',
  genital: 'of the genital system', endocrine: 'an endocrine gland', lymphoid: 'a lymphoid organ',
};

// The entry of a concept or of its nearest is-a ancestor with one.
export function knowledgeOf(atlas, knowledge, id) {
  for (let x = id, steps = 0; x && steps < 40; x = atlas.concepts.get(x)?.isa, steps++) {
    const e = knowledge.entries?.[x];
    if (e) return { entry: e, from: x === id ? null : x };
  }
  return null;
}

function topologyWords(t) {
  if (!t) return '';
  if (t.consistent) {
    const pieces = t.components === 1 ? 'one closed surface' : `${t.components} closed surfaces`;
    return t.genus ? `${pieces}, genus ${t.genus}` : `${pieces}, genus 0`;
  }
  if (t.surfaces === t.components && t.boundaryCircles) return `${t.components === 1 ? 'an open surface' : `${t.components} surfaces`} with ${t.boundaryCircles} boundary circle${t.boundaryCircles > 1 ? 's' : ''}`;
  return `not a closed surface${t.reason ? `: ${t.reason}` : ''}`;
}

function chain(atlas, id, rel) {
  const out = [];
  for (let x = id, steps = 0; x && steps < 40; x = atlas.concepts.get(x)?.[rel], steps++) out.unshift(x);
  return out;
}

export function renderCard(el, ctx) {
  const { atlas, knowledge, catalogue, selection, systemsByKey } = ctx;
  const c = atlas.concepts.get(selection.id);
  if (!c) { el.innerHTML = ''; return; }
  const elements = selection.elements;
  const single = elements.length === 1 ? catalogue.elements[elements[0]] : null;
  const k = knowledgeOf(atlas, knowledge, c.id);
  const sys = c.system || (single && single.system);
  const side = single?.side || '';
  const html = [];
  html.push(`<h1>${esc(capital(c.name))}</h1>`);
  if (k?.entry.latin) html.push(`<p class="latin">${esc(k.entry.latin)}</p>`);
  html.push('<div class="meta">');
  html.push(`<span>${esc(c.id)}</span>`);
  if (sys) html.push(`<span class="sys"><span class="sw" style="background:${systemColour(sys)}"></span>${esc(systemsByKey[sys]?.name || sys)}</span>`);
  if (elements.length > 1) html.push(`<span>${elements.length} parts</span>`);
  if (side) html.push(`<span>${side} side</span>`);
  html.push('</div>');
  html.push('<div class="actions">'
    + '<button class="key small" data-act="fly">Fly to</button>'
    + `<button class="key small" data-act="isolate">${ctx.isolated ? 'Show all' : 'Isolate'}</button>`
    + '<button class="key small" data-act="hide">Hide</button>'
    + '<button class="key small" data-act="cut">Cut through it</button>'
    + (single && single.measures.extents[0] > 3 * single.measures.extents[1] ? '<button class="key small" data-act="calibre">Calibre</button>' : '')
    + '</div>');

  if (k) {
    const e = k.entry;
    html.push('<section>');
    if (k.from) html.push(`<p class="seed">What follows is said of the ${esc(atlas.concepts.get(k.from)?.name || k.from)}, of which this is one.</p>`);
    html.push(`<p class="summary">${esc(e.summary)}</p>`);
    for (const [key, label] of [['function', 'What it does'], ['supply', 'Blood supply'], ['drainage', 'Venous drainage'], ['innervation', 'Nerves'], ['clinical', 'In the clinic']]) {
      if (e[key]) html.push(`<details${key === 'function' ? ' open' : ''}><summary><svg class="tw"><use href="#i-chev"/></svg>${label}</summary><p>${esc(e[key])}</p></details>`);
    }
    if (e.refs?.length) html.push(`<ol class="refs">${e.refs.map((r) => `<li>${esc(knowledge.references?.[r] || r)}</li>`).join('')}</ol>`);
    if (knowledge.status === 'seed') html.push('<p class="seed">A seed of the knowledge base, drafted from these references and not yet checked against them.</p>');
    html.push('</section>');
  } else {
    const words = sys ? SYSTEM_WORDS[sys] : null;
    const parent = atlas.concepts.get(c.partof) || atlas.concepts.get(c.isa);
    if (words || parent) html.push(`<section><p class="summary">${esc(capital(c.name))}${words ? `, ${words}` : ''}${parent ? `, ${c.partof ? 'part of' : 'a kind of'} the ${esc(parent.name)}` : ''}.</p></section>`);
  }

  const m = single ? single.measures : selection.part?.measures;
  const level = single ? single.level : selection.part?.level;
  html.push('<section><h3>Where it lies</h3><dl>');
  if (m) {
    const x = m.centre[0];
    const midline = Math.abs(x) < 8 ? 'on the midline' : `${fmt(Math.abs(x))} mm ${x > 0 ? 'left' : 'right'} of the midline`;
    html.push(`<dt>Level</dt><dd>${level ? `vertebra ${esc(level)}` : '—'}</dd>`);
    html.push(`<dt>Height</dt><dd>${fmt(m.centre[2] - catalogue.height.bottom)} mm above the soles</dd>`);
    html.push(`<dt>Side</dt><dd>${midline}</dd>`);
  } else html.push('<dt>…</dt><dd>measuring</dd>');
  html.push('</dl></section>');

  html.push('<section><h3>Measures</h3><dl>');
  if (m) {
    const published = single ? single.published : selection.part?.published;
    const part = selection.part;
    const partial = !single && part && m.volume === null && part.solids
      ? `at least ${fmt(part.solidVolume / 1000, 1)} cm³ <span class="small">(the ${part.solids} of ${part.elements.length} parts that are closed)</span>` : null;
    html.push(`<dt>Volume</dt><dd>${m.volume !== null ? `${fmt(m.volume / 1000, m.volume < 10000 ? 2 : 1)} cm³` : partial || '— <span class="small">(not a closed surface)</span>'}${published ? ` <span class="small">database ${fmt(published, published < 10 ? 2 : 1)}</span>` : ''}</dd>`);
    html.push(`<dt>Surface</dt><dd>${fmt(m.area / 100, 1)} cm²</dd>`);
    html.push(`<dt>Extent</dt><dd>${fmt(m.extents[0])} × ${fmt(m.extents[1])} × ${fmt(m.extents[2])} mm <span class="small">along its principal axes</span></dd>`);
    if (single) html.push(`<dt>Form</dt><dd>${esc(topologyWords(m.topology))}</dd>`);
  } else html.push('<dt>…</dt><dd>measuring</dd>');
  html.push('</dl><div class="profile" id="profile"></div></section>');

  const touches = new Map();
  for (const i of elements) for (const t of ctx.touching.get(i) || []) if (!elements.includes(t.j)) touches.set(t.j, Math.min(touches.get(t.j) ?? Infinity, t.d));
  if (touches.size) {
    const byConcept = new Map();
    for (const [j, d] of touches) {
      const cid = catalogue.elements[j].concept;
      if (!byConcept.has(cid) || byConcept.get(cid).d > d) byConcept.set(cid, { j, d });
    }
    const list = [...byConcept.entries()].sort((a, b) => a[1].d - b[1].d);
    html.push(`<section><h3>It touches</h3><div class="chips">${list.slice(0, 24).map(([cid, { j }]) => {
      const e = catalogue.elements[j];
      return `<span class="chip" data-concept="${esc(cid)}"><span class="sw" style="background:${systemColour(e.system)}"></span>${esc(capital(atlas.concepts.get(cid)?.name || e.name))}</span>`;
    }).join('')}</div>${list.length > 24 ? `<p class="seed">and ${list.length - 24} more</p>` : ''}<p class="seed">Parts within ${catalogue.touchDistance} mm of it.</p></section>`);
  }

  html.push('<section><h3>In the trees</h3>');
  for (const [rel, label] of [['partof', 'Part of'], ['isa', 'A kind of']]) {
    const ids = chain(atlas, c.id, rel);
    if (ids.length > 1) html.push(`<p class="crumbs"><span class="seed">${label}</span><br>${ids.map((x) => `<a data-concept="${esc(x)}" data-relation="${rel}">${esc(atlas.concepts.get(x)?.name || x)}</a>`).join('<span class="sep">›</span>')}</p>`);
  }
  html.push('</section>');
  el.innerHTML = html.join('');
}

// The calibre of a part along its long axis: the diameter of the circle of its cross-section's area.
export function renderProfile(el, profile) {
  const pts = profile.samples.filter((s) => s.area > 0);
  if (!pts.length) { el.innerHTML = '<p class="seed">No closed sections along its axis.</p>'; return; }
  const W = 300, H = 90, L = profile.length || 1;
  const dmax = Math.max(...pts.map((p) => p.diameter)) || 1;
  const x = (s) => 4 + (s / L) * (W - 8), y = (d) => H - 14 - (d / dmax) * (H - 26);
  const line = pts.map((p, k) => `${k ? 'L' : 'M'}${x(p.s).toFixed(1)},${y(p.diameter).toFixed(1)}`).join('');
  const area = `${line}L${x(pts[pts.length - 1].s).toFixed(1)},${H - 14}L${x(pts[0].s).toFixed(1)},${H - 14}Z`;
  const mean = pts.reduce((s, p) => s + p.diameter, 0) / pts.length;
  el.innerHTML = `<h3 style="margin-top:8px">Calibre along its length</h3><svg viewBox="0 0 ${W} ${H}" preserveAspectRatio="none">`
    + `<path class="area" d="${area}"/><path class="line" d="${line}"/><line class="axis" x1="4" x2="${W - 4}" y1="${H - 14}" y2="${H - 14}"/>`
    + `<text class="lbl" x="4" y="${H - 2}">0</text><text class="lbl" x="${W - 4}" y="${H - 2}" text-anchor="end">${fmt(L)} mm</text>`
    + `<text class="lbl" x="6" y="10">${fmt(dmax, 1)} mm</text></svg>`
    + `<p class="seed">The diameter of a circle of each cross-section's area, normal to its long axis; mean ${fmt(mean, 1)} mm. Sections integrate to ${fmt(profile.integral / 1000, 2)} cm³${profile.volume ? ` against its volume of ${fmt(profile.volume / 1000, 2)} cm³` : ''}.</p>`;
}
