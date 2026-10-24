// Circuits: a transfer network edited here, read by the engine (economics::transfer and
// circulation): every loop's holonomy and deficit in the tree gauge, the potential (log value
// units, a price system where every deficit vanishes), retention round directed cycles with
// flows, the decomposition of the flows, and their pairing with the transports.

import * as api from './api.js';
import { el, svg, clear, fixed, percent, signed, shareBar } from './dom.js';
import { PART, seriesColour } from './palette.js';
import { NETWORKS } from './presets.js';

const state = { preset: NETWORKS[0].id, network: structuredClone(NETWORKS[0].network), answer: null, fill: false };

// Faces: with `fill`, every triangle of positions whose three sides have legs.
function withFaces(net) {
  if (!state.fill) return net;
  const leg = (u, v) => {
    const i = net.legs.findIndex((l) => (l.from === u && l.to === v) || (l.from === v && l.to === u));
    return i < 0 ? null : [i, net.legs[i].from === u ? 1 : -1];
  };
  const n = net.positions.length, faces = [];
  for (let a = 0; a < n; a++) for (let b = a + 1; b < n; b++) for (let c = b + 1; c < n; c++) {
    const s = [leg(a, b), leg(b, c), leg(c, a)];
    if (s.every(Boolean)) faces.push(s);
  }
  return { ...net, faces };
}

function cycleText(path, net) {
  const names = [];
  for (const [e, sign] of path) {
    const l = net.legs[e];
    names.push(net.positions[sign > 0 ? l.from : l.to]);
  }
  return `${names.join(' → ')} → ${names[0]}`;
}

function drawing() {
  const net = state.network, A = state.answer;
  const n = net.positions.length, R = 120, cx = 220, cy = 160;
  const pos = net.positions.map((_, i) => [cx + R * Math.cos((2 * Math.PI * i) / n - Math.PI / 2), cy + R * Math.sin((2 * Math.PI * i) / n - Math.PI / 2)]);
  const g = svg('svg', { class: 'figure', width: 440, height: 330, viewBox: '0 0 440 330' },
    svg('defs', {}, svg('marker', { id: 'ci-head', viewBox: '0 0 10 10', refX: 9, refY: 5, markerWidth: 6, markerHeight: 6, orient: 'auto-start-reverse' }, svg('path', { d: 'M0,1L9,5L0,9z', fill: 'context-stroke' }))));
  net.legs.forEach((l, e) => {
    const [a, b] = [pos[l.from], pos[l.to]];
    const twin = net.legs.some((k, j) => j !== e && k.from === l.to && k.to === l.from);
    const bend = twin ? 0.22 : 0.08;
    const mx = (a[0] + b[0]) / 2 + (b[1] - a[1]) * bend, my = (a[1] + b[1]) / 2 - (b[0] - a[0]) * bend;
    const flowParts = A ? [A.flows.parts.driven[e], A.flows.parts.local[e] + A.flows.parts.global[e]] : [0, 0];
    const colour = (l.flow || 0) > 0 && Math.abs(flowParts[1]) > Math.abs(flowParts[0]) ? PART.local : PART.driven;
    const shrink = (p, q, d) => { const dx = q[0] - p[0], dy = q[1] - p[1], len = Math.hypot(dx, dy); return [p[0] + (dx * d) / len, p[1] + (dy * d) / len]; };
    const a2 = shrink(a, [mx, my], 12), b2 = shrink(b, [mx, my], 12);
    g.append(svg('path', { d: `M${a2[0]},${a2[1]}Q${mx},${my} ${b2[0]},${b2[1]}`, fill: 'none', stroke: colour, 'stroke-width': 1.6 + ((l.flow || 0) > 0 ? 1.4 : 0), 'marker-end': 'url(#ci-head)' }));
    g.append(svg('text', { x: mx, y: my, 'text-anchor': 'middle', text: `${fixed(l.ratio ?? 1, 3)}` }));
  });
  net.positions.forEach((p, i) => {
    g.append(svg('circle', { cx: pos[i][0], cy: pos[i][1], r: 8, fill: seriesColour(i) }));
    g.append(svg('text', { x: pos[i][0], y: pos[i][1] + (pos[i][1] < cy ? -14 : 24), 'text-anchor': 'middle', text: p }));
  });
  return g;
}

function editor() {
  const net = state.network;
  const options = (sel) => net.positions.map((p, i) => el('option', { value: i, selected: i === sel, text: p }));
  return el('div', {},
    el('h3', { text: 'Positions' }),
    el('div', { class: 'row-actions' }, net.positions.map((p, i) => el('input', { class: 'name', value: p, style: { width: '110px' }, onchange: (e) => { net.positions[i] = e.target.value; run(); } })),
      el('button', { onclick: () => { net.positions.push(`p${net.positions.length}`); render(); run(); } }, '+ a position')),
    el('h3', { text: 'Legs' }),
    el('table', { class: 'grid' }, el('thead', {}, el('tr', {}, ['from', 'to', 'transport U', 'flow', 'kind', ''].map((h) => el('th', { text: h })))),
      el('tbody', {}, net.legs.map((l, e) => el('tr', {},
        el('td', {}, el('select', { onchange: (ev) => { l.from = Number(ev.target.value); run(); } }, options(l.from))),
        el('td', {}, el('select', { onchange: (ev) => { l.to = Number(ev.target.value); run(); } }, options(l.to))),
        el('td', {}, el('input', { type: 'number', step: 'any', min: '0', value: l.ratio ?? 1, onchange: (ev) => { l.ratio = Number(ev.target.value); run(); } })),
        el('td', {}, el('input', { type: 'number', step: 'any', min: '0', value: l.flow ?? 0, onchange: (ev) => { l.flow = Number(ev.target.value); run(); } })),
        el('td', {}, el('input', { class: 'name', value: l.kind || '', style: { width: '90px' }, onchange: (ev) => { l.kind = ev.target.value; run(); } })),
        el('td', {}, el('button', { title: 'remove the leg', onclick: () => { net.legs.splice(e, 1); render(); run(); } }, '×')))))),
    el('div', { class: 'row-actions' },
      el('button', { onclick: () => { net.legs.push({ from: 0, to: Math.min(1, net.positions.length - 1), ratio: 1, flow: 0, kind: '' }); render(); run(); } }, '+ a leg'),
      el('label', {}, el('input', { type: 'checkbox', checked: state.fill, onchange: (e) => { state.fill = e.target.checked; run(); } }), ' fill every triangle of legs with a face')));
}

function results() {
  const box = document.getElementById('ci-panel');
  clear(box);
  const A = state.answer, net = state.network;
  if (!A) return;
  const flat = A.defect < 1e-12;
  box.append(
    el('h2', { text: 'Loops and their deficits' }),
    el('p', { class: 'lede', text: flat ? 'Every loop returns exactly what went in: there is a consistent system of prices, the potential below, and no arbitrage round any loop.' : 'Some loop does not return what went in. Its deficit is the share of value lost going round, and no re-description of value at the positions removes it.' }),
    el('table', {}, el('thead', {}, el('tr', {}, el('th', { text: 'cycle (tree gauge)' }), el('th', { class: 'num', text: 'λ' }), el('th', { class: 'num', text: 'δ' }))),
      el('tbody', {}, A.cycles.map((c) => el('tr', {}, el('td', { text: cycleText(c.path, net) }), el('td', { class: 'num', text: fixed(c.holonomy, 4) }), el('td', { class: 'num', text: fixed(c.deficit, 4) }))))),
    el('p', { class: 'note', text: 'One cycle for each leg outside a spanning tree: its holonomy λ = Π U^±1 and deficit δ = −ln λ are complete invariants of the network up to re-description of value at the positions.' }),
  );
  const retained = A.cycles.filter((c) => c.retention);
  if (retained.length) box.append(el('div', { class: 'section' }, el('h3', { text: 'Retention round directed cycles' }),
    retained.map((c) => el('p', {}, `${cycleText(c.path, net)}: ${c.retention.map((r) => fixed(r, 3)).join(' × ')} = ${fixed(c.retention.reduce((s, r) => s * r, 1), 4)} = λ`)),
    el('p', { class: 'note', text: 'Each position\'s retention is what it receives over what it passes on; their product is the holonomy whatever the flows. Which position "absorbs" the deficit depends on the flows, not on the network: it is not an invariant.' })));
  box.append(el('div', { class: 'section' }, el('h3', { text: flat ? 'Prices' : 'The potential of the tree' }),
    el('table', {}, el('tbody', {}, net.positions.map((p, i) => el('tr', {}, el('td', { text: p }), el('td', { class: 'num', text: signed(A.log_price[i], 4) }))))),
    el('p', { class: 'note', text: flat ? 'The potential φ with ln U = φ(head) − φ(tail): minus the log price, up to a constant, since a transport counts the units received for each unit given.' : 'Integrated along the spanning tree from its root; the chords then carry the deficits.' })));
  if ((A.curvature || []).length) box.append(el('p', {}, `Curvature on the faces: ${A.curvature.map((k) => signed(k, 4)).join(', ')} (the deficit round each face, negated; by Stokes they add up over any region).`));
  const F = A.flows, hasFlow = net.legs.some((l) => (l.flow || 0) > 0);
  if (hasFlow) box.append(el('div', { class: 'section' }, el('h3', { text: 'The flows' }),
    shareBar([{ label: 'driven', value: F.driven, colour: PART.driven }, { label: 'round faces', value: F.local, colour: PART.local }, { label: 'round holes', value: F.global, colour: PART.global }]),
    el('p', {}, `Pairing with the transports, ⟨F, ln U⟩ = ${signed(A.pairing.total, 4)}: ${signed(A.pairing.prices, 4)} by driven flow against prices, ${signed(A.pairing.local + A.pairing.global, 4)} by circulation against the deficits. Only what circulates meets a deficit.`)));
}

function render() {
  const main = document.getElementById('ci-main');
  clear(main);
  const preset = NETWORKS.find((p) => p.id === state.preset);
  main.append(
    el('h2', { text: 'Circuits of value' }),
    el('p', { class: 'lede', text: 'A leg carries value from one position to another with a transport U, the value registered at its head for each unit surrendered at its tail; re-describing value at the positions changes every U but no loop\'s product. These are the invariants of the user\'s paper on loop deficits, computed by the economics domain.' }),
    el('div', { class: 'presets' }, NETWORKS.map((p) => el('button', { class: p.id === state.preset ? 'active' : '', onclick: () => { state.preset = p.id; state.network = structuredClone(p.network); render(); run(); } }, p.title))),
    preset ? el('p', { class: 'note', text: preset.about }) : null,
    el('div', { id: 'ci-figure' }),
    editor(),
  );
}

async function run() {
  try {
    document.getElementById('status').textContent = 'reading the loops…';
    state.answer = await api.circuits(withFaces(state.network));
    document.getElementById('status').textContent = '';
  } catch (error) {
    document.getElementById('status').textContent = error.message;
    state.answer = null;
  }
  const fig = document.getElementById('ci-figure');
  if (fig) clear(fig).append(drawing());
  results();
}

export async function showCircuits() {
  if (!document.getElementById('ci-main').firstChild) render();
  await run();
}
