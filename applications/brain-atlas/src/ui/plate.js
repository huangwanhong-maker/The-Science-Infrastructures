// Labels set out as on an engraved plate of an anatomical atlas: names in two
// columns beside the brain, each joined to its structure by a fine leader line
// ending in a dot on the part of the surface the viewer sees.

import { esc } from './card.js';

const GAP = 21;
const FONT = 'italic 13px "Source Serif 4", Georgia, serif';
const BOLD = '600 13px "Source Serif 4", Georgia, serif';
let measurer = null;
// The width of a label as drawn, in pixels.
function textWidth(text, bold = false) {
  measurer = measurer || document.createElement('canvas').getContext('2d');
  measurer.font = bold ? BOLD : FONT;
  return measurer.measureText(text).width + 20;
}

// candidates: [{id, side, name, colour, inst, samples: Float32Array (x, y, z, nx, ny, nz per sample), selected, weight}]
// A name is set once: of two candidates with one name (the two sides of a nucleus), the
// one better seen is labelled.
export function layoutPlate({ viewer, candidates, max = 12, svg, labelsEl }) {
  const cam = viewer.camera.position;
  const W = viewer.width, H = viewer.height;
  const scored = [];
  for (const c of candidates) {
    const s = c.samples;
    const off = c.inst.mesh.position;
    const pts = [];
    let sx = 0, sy = 0;
    for (let i = 0; i < s.length; i += 6) {
      const px = s[i] + off.x, py = s[i + 1] + off.y, pz = s[i + 2] + off.z;
      const facing = s[i + 3] * (cam.x - px) + s[i + 4] * (cam.y - py) + s[i + 5] * (cam.z - pz);
      if (facing <= 0) continue;
      if (viewer.cut && c.inst.mesh.material.clippingPlanes?.length && viewer.clipPlane.distanceToPoint({ x: px, y: py, z: pz }) < 0) continue;
      const q = viewer.project([s[i], s[i + 1], s[i + 2]], c.inst);
      if (!q || q.x < 0 || q.y < 0 || q.x > W || q.y > H) continue;
      pts.push([q.x, q.y, i]);
      sx += q.x; sy += q.y;
    }
    if (pts.length < 2 && !c.selected) continue;
    if (!pts.length) continue;
    scored.push({ c, pts, cx: sx / pts.length, cy: sy / pts.length, score: pts.length * (c.weight || 1) + (c.selected ? 1e9 : 0) });
  }
  scored.sort((a, b) => b.score - a.score);
  const fit = Math.max(4, Math.floor((H - 130) / GAP));
  const chosen = [];
  const named = new Set();
  for (const s of scored) {
    if (chosen.length >= Math.min(max, 2 * fit)) break;
    if (named.has(s.c.name)) continue;
    // The sample nearest the centre of the seen part, that is truly in view.
    const order = s.pts.map((p) => [(p[0] - s.cx) ** 2 + (p[1] - s.cy) ** 2, p]).sort((a, b) => a[0] - b[0]);
    let anchor = null;
    for (const [, p] of order.slice(0, 8)) {
      const i = p[2], smp = s.c.samples;
      if (viewer.visibleAt([smp[i], smp[i + 1], smp[i + 2]], s.c.inst)) { anchor = p; break; }
    }
    if (!anchor) continue;
    named.add(s.c.name);
    chosen.push({ ...s, ax: anchor[0], ay: anchor[1], w: textWidth(s.c.name, s.c.selected) });
  }
  if (!chosen.length) { svg.innerHTML = ''; labelsEl.innerHTML = ''; return []; }
  // Two columns beside the brain's extent on screen; a label goes to the side of its
  // anchor, and a crowded column gives its innermost labels to the other.
  const box = viewer.screenBox ? viewer.screenBox() : { x0: W * 0.25, x1: W * 0.75 };
  const mid = (box.x0 + box.x1) / 2;
  const cols = { left: [], right: [] };
  for (const s of chosen) (s.ax < mid ? cols.left : cols.right).push(s);
  const balance = (from, to) => {
    from.sort((a, b) => Math.abs(a.ax - mid) - Math.abs(b.ax - mid));
    while (from.length > fit || from.length > to.length + 4) to.push(from.shift());
  };
  if (cols.left.length > cols.right.length) balance(cols.left, cols.right); else balance(cols.right, cols.left);
  const wL = Math.max(0, ...cols.left.map((s) => s.w)), wR = Math.max(0, ...cols.right.map((s) => s.w));
  const leftX = Math.max(wL + 12, Math.min(box.x0 - 24, W / 2 - 40));
  const rightX = Math.min(W - wR - 12, Math.max(box.x1 + 24, W / 2 + 40));
  const spread = (list) => {
    list.sort((a, b) => a.ay - b.ay);
    let y = 34;
    for (const s of list) { s.ly = Math.max(s.ay, y); y = s.ly + GAP; }
    const over = y - GAP - (H - 104);
    if (over > 0) for (const s of list) s.ly -= over;
    for (let i = 1; i < list.length; i++) if (list[i].ly < list[i - 1].ly + GAP) list[i].ly = list[i - 1].ly + GAP;
  };
  spread(cols.left); spread(cols.right);
  const paths = [], labels = [];
  for (const [side, list] of Object.entries(cols)) {
    const x = side === 'left' ? leftX : rightX;
    for (const s of list) {
      const elbow = side === 'left' ? x + 14 : x - 14;
      paths.push(`<path d="M${s.ax.toFixed(1)},${s.ay.toFixed(1)} L${elbow.toFixed(1)},${s.ly.toFixed(1)} L${x.toFixed(1)},${s.ly.toFixed(1)}"${s.c.selected ? ' class="selected"' : ''}/><circle cx="${s.ax.toFixed(1)}" cy="${s.ay.toFixed(1)}" r="2.2"${s.c.selected ? ' class="selected"' : ''}/>`);
      labels.push(`<div class="plate-label ${side}${s.c.selected ? ' selected' : ''}" data-select="${esc(s.c.id)}" data-side="${esc(s.c.side || '')}" style="left:${(side === 'left' ? x - 4 : x + 4).toFixed(1)}px;top:${s.ly.toFixed(1)}px"><i style="background:${esc(s.c.colour || 'transparent')}"></i>${esc(s.c.name)}</div>`);
    }
  }
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.innerHTML = paths.join('');
  labelsEl.innerHTML = labels.join('');
  return chosen;
}

export function clearPlate(svg, labelsEl) { svg.innerHTML = ''; labelsEl.innerHTML = ''; }

// Up to n samples (position and normal) of a mesh's vertices, or of a subset of them.
export function samplesOf(inst, indices = null, n = 260) {
  const pos = inst.mesh.geometry.attributes.position.array, nrm = inst.mesh.geometry.attributes.normal.array;
  const count = indices ? indices.length : pos.length / 3;
  const step = Math.max(1, Math.floor(count / n));
  const out = new Float32Array(Math.ceil(count / step) * 6);
  let k = 0;
  for (let j = 0; j < count; j += step) {
    const i = indices ? indices[j] : j;
    out[k++] = pos[3 * i]; out[k++] = pos[3 * i + 1]; out[k++] = pos[3 * i + 2];
    out[k++] = nrm[3 * i]; out[k++] = nrm[3 * i + 1]; out[k++] = nrm[3 * i + 2];
  }
  return out.subarray(0, k);
}
