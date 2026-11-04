// The plate: a section of the body as the engine cut it, drawn as an engraved plate: each piece filled in
// its system's colour, the larger pieces named in the margins on fine leader lines, the subject's directions
// lettered at the edges as an imaging viewer letters them.

import { shade, systemColour } from './palette.js';
import { capital } from './tree.js';

const NS = 'http://www.w3.org/2000/svg';
const LETTER = { left: 'L', right: 'R', anterior: 'A', posterior: 'P', superior: 'S', inferior: 'I' };
const OPPOSITE = { left: 'right', right: 'left', anterior: 'posterior', posterior: 'anterior', superior: 'inferior', inferior: 'superior' };
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

function ringArea(flat) {
  let s = 0;
  for (let k = 0; k < flat.length; k += 2) {
    const j = (k + 2) % flat.length;
    s += flat[k] * flat[j + 1] - flat[j] * flat[k + 1];
  }
  return s / 2;
}

// A point well inside the largest loop: its centroid when inside, else the midpoint of the widest horizontal
// chord through its centroid's height.
function anchor(flat) {
  let a = 0, cx = 0, cy = 0;
  for (let k = 0; k < flat.length; k += 2) {
    const j = (k + 2) % flat.length;
    const f = flat[k] * flat[j + 1] - flat[j] * flat[k + 1];
    a += f; cx += (flat[k] + flat[j]) * f; cy += (flat[k + 1] + flat[j + 1]) * f;
  }
  if (a === 0) return [flat[0], flat[1]];
  cx /= 3 * a; cy /= 3 * a;
  const xs = [];
  for (let k = 0; k < flat.length; k += 2) {
    const j = (k + 2) % flat.length;
    const y0 = flat[k + 1], y1 = flat[j + 1];
    if ((y0 > cy) !== (y1 > cy)) xs.push(flat[k] + ((cy - y0) / (y1 - y0)) * (flat[j] - flat[k]));
  }
  xs.sort((p, q) => p - q);
  let best = null;
  for (let k = 0; k + 1 < xs.length; k += 2) if (!best || xs[k + 1] - xs[k] > best[1] - best[0]) best = [xs[k], xs[k + 1]];
  if (best && !(cx > best[0] && cx < best[1])) cx = (best[0] + best[1]) / 2;
  return [cx, cy];
}

export function renderPlate(svg, section, ctx) {
  const { catalogue, atlas, selected, onChoose, onHover } = ctx;
  const W = svg.clientWidth || 600, H = svg.clientHeight || 600;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.innerHTML = '';
  if (!section || !section.pieces.length) {
    const t = document.createElementNS(NS, 'text');
    t.setAttribute('x', W / 2); t.setAttribute('y', H / 2); t.setAttribute('text-anchor', 'middle'); t.setAttribute('class', 'orient');
    t.textContent = section ? 'Nothing shown is cut here.' : '';
    svg.appendChild(t);
    return;
  }
  // Fitted to every piece of the section (the trunk's outline when the skin is shown; the organs when it is
  // hidden, where the largest loop alone would be one organ), then zoomed and moved as the reader asks (ctx.view).
  let lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
  for (const l of section.pieces.flatMap((p) => [...p.loops, ...p.chains])) for (let k = 0; k < l.length; k += 2) {
    lo = [Math.min(lo[0], l[k]), Math.min(lo[1], l[k + 1])];
    hi = [Math.max(hi[0], l[k]), Math.max(hi[1], l[k + 1])];
  }
  // A transverse section's outline may join the arms to the trunk: its width is fitted to the trunk's.
  if (section.plane === 'transverse') { lo[0] = Math.max(lo[0], -200); hi[0] = Math.min(hi[0], 200); }
  const pad = 0.06 * Math.max(hi[0] - lo[0], hi[1] - lo[1]);
  lo = [lo[0] - pad, lo[1] - pad]; hi = [hi[0] + pad, hi[1] + pad];
  const view = ctx.view || { zoom: 1, du: 0, dv: 0 };
  const margin = { x: Math.min(158, W * 0.24), y: 30 };
  const scale = view.zoom * Math.min((W - 2 * margin.x) / (hi[0] - lo[0] || 1), (H - 2 * margin.y) / (hi[1] - lo[1] || 1));
  const cu = (lo[0] + hi[0]) / 2 + view.du, cv = (lo[1] + hi[1]) / 2 + view.dv;
  const ox = W / 2 - cu * scale, oy = H / 2 + cv * scale;
  ctx.transform = { scale, cu, cv, W, H };
  const clipId = 'plate-clip';
  const defs = document.createElementNS(NS, 'defs');
  defs.innerHTML = `<clipPath id="${clipId}"><rect x="${margin.x - 18}" y="0" width="${W - 2 * margin.x + 36}" height="${H}"/></clipPath>`;
  svg.appendChild(defs);
  const X = (u) => ox + u * scale, Y = (v) => oy - v * scale;
  const pathOf = (l, close) => {
    let d = '';
    for (let k = 0; k < l.length; k += 2) d += `${k ? 'L' : 'M'}${X(l[k]).toFixed(1)},${Y(l[k + 1]).toFixed(1)}`;
    return close ? `${d}Z` : d;
  };
  // Larger pieces first, so that small ones lie on top and can be reached.
  const pieces = [...section.pieces].sort((a, b) => b.area - a.area);
  const gPieces = document.createElementNS(NS, 'g');
  gPieces.setAttribute('clip-path', `url(#${clipId})`);
  for (const p of pieces) {
    const e = catalogue.elements[p.e];
    const colour = systemColour(e.system);
    if (p.loops.length) {
      const path = document.createElementNS(NS, 'path');
      path.setAttribute('d', p.loops.map((l) => pathOf(l, true)).join(''));
      path.setAttribute('class', `piece${selected.has(p.e) ? ' on' : ''}`);
      path.setAttribute('fill', colour);
      path.setAttribute('stroke', shade(colour));
      path.setAttribute('fill-rule', 'nonzero');
      path.addEventListener('click', () => onChoose(p.e));
      path.addEventListener('mousemove', (ev) => onHover(p.e, ev));
      path.addEventListener('mouseleave', () => onHover(null));
      gPieces.appendChild(path);
    }
    for (const c of p.chains) {
      const path = document.createElementNS(NS, 'path');
      path.setAttribute('d', pathOf(c, false));
      path.setAttribute('class', 'chain');
      path.setAttribute('stroke', shade(colour));
      gPieces.appendChild(path);
    }
  }
  svg.appendChild(gPieces);

  // Names: the largest pieces of distinct concepts within the frame, and the chosen ones, set in the margins.
  const seen = new Set();
  const named = [];
  for (const p of pieces) {
    if (!p.loops.length) continue;
    const e = catalogue.elements[p.e];
    const chosen = selected.has(p.e);
    if (seen.has(e.concept) && !chosen) continue;
    if (named.length >= 22 && !chosen) continue;
    if (p.area * scale * scale < 60 && !chosen) continue;
    const big = p.loops.reduce((a, l) => (ringArea(l) > ringArea(a) ? l : a), p.loops[0]);
    const [u, v] = anchor(big);
    const x = X(u), y = Y(v);
    if (x < margin.x - 10 || x > W - margin.x + 10 || y < 4 || y > H - 4) continue;
    seen.add(e.concept);
    named.push({ p, e, x, y });
  }
  const columns = { left: named.filter((n) => n.x < W / 2).sort((a, b) => a.y - b.y), right: named.filter((n) => n.x >= W / 2).sort((a, b) => a.y - b.y) };
  const gTags = document.createElementNS(NS, 'g');
  svg.appendChild(gTags);
  for (const [side, list] of Object.entries(columns)) {
    const gap = 17;
    const ys = list.map((n) => n.y);
    for (let k = 1; k < ys.length; k++) ys[k] = Math.max(ys[k], ys[k - 1] + gap);
    const over = ys.length ? ys[ys.length - 1] - (H - 14) : 0;
    if (over > 0) for (let k = ys.length - 1; k >= 0; k--) ys[k] = Math.min(ys[k] - over, k + 1 < ys.length ? ys[k + 1] - gap : Infinity);
    list.forEach((n, k) => {
      const text = document.createElementNS(NS, 'text');
      const tx = side === 'left' ? 12 : W - 12;
      text.setAttribute('x', tx); text.setAttribute('y', ys[k] + 4);
      text.setAttribute('text-anchor', side === 'left' ? 'start' : 'end');
      text.setAttribute('class', `tag${selected.has(n.p.e) ? ' on' : ''}`);
      let name = capital(atlas.concepts.get(n.e.concept)?.name || n.e.name);
      text.textContent = name;
      gTags.appendChild(text);
      const room = margin.x - 34;
      while (text.getComputedTextLength() > room && name.length > 6) { name = name.slice(0, -2); text.textContent = `${name.trimEnd()}…`; }
      const end = side === 'left' ? tx + text.getComputedTextLength() + 5 : tx - text.getComputedTextLength() - 5;
      const elbow = side === 'left' ? margin.x - 16 : W - margin.x + 16;
      const leader = document.createElementNS(NS, 'path');
      leader.setAttribute('class', 'leader');
      leader.setAttribute('d', `M${n.x.toFixed(1)},${n.y.toFixed(1)}L${elbow.toFixed(1)},${ys[k].toFixed(1)}L${end.toFixed(1)},${ys[k].toFixed(1)}`);
      const dot = document.createElementNS(NS, 'circle');
      dot.setAttribute('cx', n.x); dot.setAttribute('cy', n.y); dot.setAttribute('r', 1.7); dot.setAttribute('class', 'tag');
      text.addEventListener('click', () => onChoose(n.p.e));
      gTags.append(leader, dot);
    });
  }

  // The subject's directions at the edges.
  const letters = [[LETTER[section.right], W - margin.x + 8, 22], [LETTER[OPPOSITE[section.right]], margin.x - 8, 22],
    [LETTER[section.up], W / 2, 18], [LETTER[OPPOSITE[section.up]], W / 2, H - 10]];
  for (const [l, x, y] of letters) {
    const t = document.createElementNS(NS, 'text');
    t.setAttribute('x', x); t.setAttribute('y', y); t.setAttribute('text-anchor', 'middle'); t.setAttribute('class', 'orient');
    t.textContent = l;
    svg.appendChild(t);
  }
  // A scale bar of 50 mm.
  const bar = document.createElementNS(NS, 'g');
  const len = 50 * scale, bx = W - margin.x - len, by = H - 14;
  bar.innerHTML = `<line x1="${bx}" x2="${bx + len}" y1="${by}" y2="${by}" class="leader" style="stroke-width:1.2"/><text x="${bx + len / 2}" y="${by - 4}" text-anchor="middle" class="orient">50 mm</text>`;
  svg.appendChild(bar);
}

export function describeSection(section) {
  if (!section) return '';
  const plane = { transverse: 'Transverse', coronal: 'Coronal', sagittal: 'Sagittal' }[section.plane];
  return `${plane} section, seen from ${esc(section.seenFrom === 'inferior' ? 'below' : section.seenFrom === 'superior' ? 'above' : `the ${section.seenFrom === 'left' ? 'left' : section.seenFrom === 'right' ? 'right' : section.seenFrom === 'anterior' ? 'front' : 'back'}`)}`;
}

// Zoom with the ctrl key (or a pinch) about the pointer, and move by dragging; the plain wheel steps the level
// (shift: five steps); a double click on the ground fits the plate again.
export function attachPlateControls(svg, getCtx, redraw, step) {
  let drag = null;
  svg.addEventListener('wheel', (ev) => {
    ev.preventDefault();
    const ctx = getCtx();
    const t = ctx.transform;
    if ((ev.ctrlKey || ev.metaKey) && t) {
      const r = svg.getBoundingClientRect();
      const px = ev.clientX - r.left - t.W / 2, py = ev.clientY - r.top - t.H / 2;
      const u = t.cu + px / t.scale, v = t.cv - py / t.scale;
      const f = Math.exp(-ev.deltaY * 0.0015);
      const view = ctx.view;
      const zoom = Math.min(12, Math.max(0.6, view.zoom * f));
      const s2 = t.scale * (zoom / view.zoom);
      view.du += (u - px / s2) - t.cu;
      view.dv += (v + py / s2) - t.cv;
      view.zoom = zoom;
      redraw();
    } else step(ev.deltaY > 0 ? -1 : 1, ev.shiftKey);
  }, { passive: false });
  svg.addEventListener('pointerdown', (ev) => {
    if (ev.target.closest('.piece, .tag')) return;
    drag = [ev.clientX, ev.clientY];
    svg.setPointerCapture(ev.pointerId);
  });
  svg.addEventListener('pointermove', (ev) => {
    if (!drag) return;
    const ctx = getCtx();
    const t = ctx.transform;
    if (!t) return;
    ctx.view.du -= (ev.clientX - drag[0]) / t.scale;
    ctx.view.dv += (ev.clientY - drag[1]) / t.scale;
    drag = [ev.clientX, ev.clientY];
    redraw();
  });
  svg.addEventListener('pointerup', () => { drag = null; });
  svg.addEventListener('dblclick', (ev) => {
    if (ev.target.closest('.piece, .tag')) return;
    Object.assign(getCtx().view, { zoom: 1, du: 0, dv: 0 });
    redraw();
  });
}
