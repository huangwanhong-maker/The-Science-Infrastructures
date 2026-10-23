// The net: every 2-cell a regular polygon with unit sides, read anticlockwise along its attaching
// word, glued to its neighbours where that does not overlap (net() in visualization/js/complex.js).
// A side is drawn in its edge's colour with an arrow from tail to head and the edge's name; a
// side glued inside the net is dashed; a corner is dotted in its vertex's colour, so that the
// identifications of corners can be read. A field colours sides by link angle and polygons by
// plaquette. The drawing pans with the pointer and zooms with the wheel.

import { svg, clear } from '../dom.js';
import { categoricalColour, cyclicColour, divergingColour, cssColour } from '../../../../visualization/js/index.js';
import { selectionSets } from '../selection.js';

const MARK = '#f6cf7e', MARK_NEG = '#b9a2e0';

// The vertex at the start of each side of each face, walking its word from its base.
function corners(complex) {
  return (complex.faces || []).map((f) => {
    const out = [];
    let at = f.base;
    for (const k of f.boundary) {
      out.push(at);
      const [t, h] = complex.edges[Math.abs(k) - 1];
      at = k > 0 ? h : t;
    }
    return out.length ? out : [f.base];
  });
}

export function drawNet(host, complex, result, { names, style, fields, selection, onHover, fit }) {
  clear(host);
  const usable = fields && fields.link_angle?.length === complex.edges.length ? fields : null;
  const marks = selectionSets(selection, complex);
  const vertexOf = corners(complex);
  const many = result.polygons.length > 400;
  const root = svg('g', {});
  host.append(root);
  const b = result.bbox, pad = 0.7;
  if (fit || !host.__view) host.__view = { x: b.x0 - pad, y: b.y0 - pad, w: Math.max(1, b.x1 - b.x0 + 2 * pad), h: Math.max(1, b.y1 - b.y0 + 2 * pad) };
  setView(host);
  attachPanZoom(host);
  // The y axis of the net points up; SVG's points down.
  root.setAttribute('transform', `translate(0 ${b.y0 + b.y1}) scale(1 -1)`);
  const flipY = (p) => [p[0], b.y0 + b.y1 - p[1]];
  const text = (p, s, attrs) => svg('text', { x: flipY(p)[0], y: flipY(p)[1], ...attrs }, s);
  const labels = svg('g', {});
  const lw = 0.035 * style.edgeWidth / 2;

  for (const p of result.polygons) {
    if (!p) continue;
    const face = p.face;
    const plaq = usable?.plaquettes?.[face];
    const marked = marks.faces.has(face);
    const fill = marked ? MARK : plaq !== undefined ? cssColour(divergingColour(plaq)) : '#c08081';
    const opacity = marked ? 0.42 : plaq !== undefined ? 0.6 : style.faceOpacity * 0.55;
    const hover = () => onHover?.({ dim: 2, index: face });
    if (p.kind === 'disc') {
      root.append(svg('circle', { cx: p.centre[0], cy: p.centre[1], r: p.radius, fill, 'fill-opacity': opacity, onpointerenter: hover }));
    } else if (p.kind === 'monogon') {
      const c = p.sides[0].curve;
      root.append(svg('circle', { cx: c.centre[0], cy: c.centre[1], r: c.r, fill, 'fill-opacity': opacity, onpointerenter: hover }));
    } else if (p.kind === 'digon') {
      const [s0, s1] = p.sides;
      root.append(svg('path', { d: `M${s0.from} Q${s0.curve} ${s0.to} Q${s1.curve} ${s1.to}Z`, fill, 'fill-opacity': opacity, onpointerenter: hover }));
    } else {
      root.append(svg('path', { d: `M${p.points.map((q) => q.join(' ')).join('L')}Z`, fill, 'fill-opacity': opacity, onpointerenter: hover }));
    }
    if (!many) labels.append(text(p.centre, `f${face + 1}`, { class: 'net-face', 'font-size': 0.17, fill: '#9a8496', 'text-anchor': 'middle', 'dominant-baseline': 'middle', 'font-style': 'italic' }));

    for (const s of p.sides) {
      const k = marks.edges.get(s.edge);
      const colour = k !== undefined ? (k < 0 ? MARK_NEG : MARK) : usable ? cssColour(cyclicColour(usable.link_angle[s.edge])) : cssColour(categoricalColour(s.edge));
      const width = k !== undefined ? lw * 2.2 : lw;
      const common = { fill: 'none', stroke: colour, 'stroke-width': width, 'stroke-linecap': 'round', 'stroke-opacity': s.glued && k === undefined ? 0.45 : 1, 'stroke-dasharray': s.glued ? `${lw * 2} ${lw * 2.5}` : null, onpointerenter: () => onHover?.({ dim: 1, index: s.edge }) };
      let mid, dir;
      if (s.curve && !Array.isArray(s.curve)) {
        root.append(svg('circle', { cx: s.curve.centre[0], cy: s.curve.centre[1], r: s.curve.r, ...common }));
        mid = [s.curve.centre[0], s.curve.centre[1] + s.curve.r];
        dir = [s.sign > 0 ? -1 : 1, 0];
      } else if (s.curve) {
        root.append(svg('path', { d: `M${s.from} Q${s.curve} ${s.to}`, ...common }));
        mid = [0.25 * s.from[0] + 0.5 * s.curve[0] + 0.25 * s.to[0], 0.25 * s.from[1] + 0.5 * s.curve[1] + 0.25 * s.to[1]];
        dir = [s.head[0] - s.tail[0], s.head[1] - s.tail[1]];
      } else {
        root.append(svg('path', { d: `M${s.from}L${s.to}`, ...common }));
        mid = [(s.from[0] + s.to[0]) / 2, (s.from[1] + s.to[1]) / 2];
        dir = [s.head[0] - s.tail[0], s.head[1] - s.tail[1]];
      }
      const len = Math.hypot(dir[0], dir[1]) || 1;
      const u = [dir[0] / len, dir[1] / len], n = [-u[1], u[0]], a = 0.075;
      if (style.arrows) root.append(svg('path', { d: `M${mid[0] + u[0] * a} ${mid[1] + u[1] * a}L${mid[0] - u[0] * a + n[0] * a * 0.75} ${mid[1] - u[1] * a + n[1] * a * 0.75}L${mid[0] - u[0] * a - n[0] * a * 0.75} ${mid[1] - u[1] * a - n[1] * a * 0.75}Z`, fill: colour, 'fill-opacity': s.glued && k === undefined ? 0.5 : 1 }));
      if (style.labels && !many && !(s.glued && s.mate && (s.mate[0] < face || (s.mate[0] === face && s.mate[1] < s.index)))) {
        // An outer side is named outside its polygon; a glued pair once, small, beside its line.
        const out = [mid[0] - p.centre[0], mid[1] - p.centre[1]];
        const ol = Math.hypot(out[0], out[1]) || 1;
        const off = s.glued ? 0.08 : 0.15;
        const at = [mid[0] + (out[0] / ol) * off + n[0] * (s.glued ? 0.09 : 0), mid[1] + (out[1] / ol) * off + n[1] * (s.glued ? 0.09 : 0)];
        labels.append(text(at, names[s.edge] ?? `e${s.edge + 1}`, { 'font-size': s.glued ? 0.13 : 0.19, fill: colour, 'fill-opacity': s.glued ? 0.75 : 1, 'text-anchor': 'middle', 'dominant-baseline': 'middle', 'font-style': 'italic' }));
      }
    }
    // Corners, coloured by the vertex they are.
    const pts = p.kind === 'monogon' || p.kind === 'disc' ? [p.kind === 'disc' ? [p.centre[0], p.centre[1] - p.radius] : p.sides[0].from] : p.points;
    pts.forEach((q, i) => {
      const v = vertexOf[face][i];
      const marked = marks.vertices.has(v);
      root.append(svg('circle', { cx: q[0], cy: q[1], r: (marked ? 0.06 : 0.04) * style.vertexSize, fill: marked ? MARK : cssColour(categoricalColour(v + 5)), stroke: '#17111a', 'stroke-width': 0.012, onpointerenter: () => onHover?.({ dim: 0, index: v }) }));
    });
  }
  host.append(labels);
  host.onpointerleave = () => onHover?.(null);
}

function setView(host) {
  const v = host.__view;
  host.setAttribute('viewBox', `${v.x} ${v.y} ${v.w} ${v.h}`);
  host.setAttribute('preserveAspectRatio', 'xMidYMid meet');
}

// Pan by dragging, zoom about the pointer with the wheel; attached once per element.
export function attachPanZoom(host) {
  if (host.__panzoom) return;
  host.__panzoom = true;
  const toWorld = (e) => {
    const r = host.getBoundingClientRect(), v = host.__view;
    const s = Math.max(v.w / r.width, v.h / r.height);
    const ox = v.x + (v.w - r.width * s) / 2, oy = v.y + (v.h - r.height * s) / 2;
    return { x: ox + (e.clientX - r.left) * s, y: oy + (e.clientY - r.top) * s, s };
  };
  host.addEventListener('wheel', (e) => {
    e.preventDefault();
    const p = toWorld(e), v = host.__view, f = Math.exp(e.deltaY * 0.0015);
    host.__view = { x: p.x - (p.x - v.x) * f, y: p.y - (p.y - v.y) * f, w: v.w * f, h: v.h * f };
    setView(host);
  }, { passive: false });
  let drag = null;
  host.addEventListener('pointerdown', (e) => { drag = { x: e.clientX, y: e.clientY, v: { ...host.__view }, s: toWorld(e).s }; host.classList.add('dragging'); host.setPointerCapture(e.pointerId); });
  host.addEventListener('pointermove', (e) => {
    if (!drag) return;
    host.__view = { ...drag.v, x: drag.v.x - (e.clientX - drag.x) * drag.s, y: drag.v.y - (e.clientY - drag.y) * drag.s };
    setView(host);
  });
  const end = () => { drag = null; host.classList.remove('dragging'); };
  host.addEventListener('pointerup', end);
  host.addEventListener('pointercancel', end);
}
