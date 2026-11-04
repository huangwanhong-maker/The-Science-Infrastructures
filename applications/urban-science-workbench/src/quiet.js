// The workbench's small pieces of page, for every plug-in's analyses: tabs, segmented keys, switches, a row of
// figures, swatches, folds, findings in brief, sliders and hints, so that each analysis shows one thing at a time, its
// numbers before its prose, and keeps its explanations folded away. They are drawn with the workbench's h() and
// styled by its stylesheet (§7 of styles.css, the ld- classes, after Landscape Dynamics, which first used them).

// Tabs: a row of names, one panel shown. Returns { el, show(name) }; `remember` keeps the choice across openings,
// and show() without a name shows the one remembered. `lazy` draws no panel until the first show(), for panels that
// wait on a result.
const chosen = new Map();
export function tabs(h, items, { remember, onShow, lazy } = {}) {
  const bar = h('div', { class: 'ld-tabs', role: 'tablist' });
  const body = h('div', { class: 'ld-tab-body' });
  const show = (name) => {
    const want = name || (remember && chosen.get(remember));
    const it = items.find((x) => x.name === want) || items[0];
    if (remember) chosen.set(remember, it.name);
    for (const b of bar.children) b.classList.toggle('on', b.dataset.name === it.name);
    body.replaceChildren(typeof it.panel === 'function' ? it.panel() : it.panel);
    if (onShow) onShow(it.name);
  };
  for (const it of items) bar.append(h('button', { class: 'ld-tab', role: 'tab', 'data-name': it.name, onclick: () => show(it.name) }, it.name));
  if (!lazy) show();
  return { el: h('div', {}, bar, body), show };
}

// Segmented keys: a few choices, one on.
export function seg(h, options, value, onChange) {
  const el = h('div', { class: 'ld-seg' });
  const draw = (v) => el.replaceChildren(...options.map(([k, name]) => h('button', { class: k === v ? 'on' : '', onclick: () => { draw(k); onChange(k); } }, name)));
  draw(value);
  return el;
}

// A switch: a labelled on/off.
export function toggle(h, label, on, onChange) {
  const input = h('input', { type: 'checkbox', checked: on || null });
  input.addEventListener('change', () => onChange(input.checked));
  return h('label', { class: 'ld-switch' }, input, h('span', { class: 'ld-track' }), label);
}

// A row of figures: [{ value, label }].
export function figures(h, items) {
  return h('div', { class: 'ld-figures' }, items.filter(Boolean).map((x) => h('div', { class: 'ld-figure' }, h('b', {}, x.value), h('span', {}, x.label))));
}

// Swatches: [{ color, label, kind: 'block' | 'line' | 'dot' }].
export function swatches(h, items) {
  return h('div', { class: 'ld-swatches' }, items.map((x) => h('span', {}, h('i', { class: `ld-sw ${x.kind || 'block'}`, style: `background:${x.color}` }), x.label)));
}

// A fold: a quiet summary line that opens onto what it holds.
export function fold(h, summary, ...children) {
  return h('details', { class: 'ld-fold' }, h('summary', {}, summary), h('div', { class: 'ld-fold-body' }, ...children));
}

// A finding, briefly: its kind, title and statement shown; how it was obtained and its limits folded.
export function brief(ctx, { title, kind, statement, method, limits, against }) {
  const { h } = ctx.dom;
  const card = ctx.finding({ title, kind, statement });
  const more = [method && h('p', {}, h('i', {}, 'How. '), method), against && h('p', {}, h('i', {}, 'Tested against. '), against), limits && h('p', {}, h('i', {}, 'Limits. '), limits)].filter(Boolean);
  if (more.length) card.append(fold(h, 'how and its limits', ...more));
  return card;
}

// A slider with its value, compact.
export function slider(h, fmt, { label, min, max, step, value, unit = '', onChange }) {
  const input = h('input', { type: 'range', min, max, step, value });
  const out = h('output', {}, `${fmt(value)}${unit}`);
  input.addEventListener('input', () => { out.textContent = `${fmt(+input.value)}${unit}`; });
  input.addEventListener('change', () => onChange(+input.value));
  return label ? h('label', { class: 'ld-slider' }, h('span', {}, label), input, out) : h('label', { class: 'ld-slider bare' }, input, out);
}

// A hint: one quiet line.
export const hint = (h, text) => h('p', { class: 'ld-hint' }, text);
