// Label placement without overlap.
//
// declutter: for HTML labels in a 3D scene, hides any label whose box on screen
//   would overlap a label of higher priority (entries are given in priority order).
// LabelPlacer: for canvas charts, tries a few positions around a point and keeps
//   the first one that is free.

import * as THREE from 'three';

const v = new THREE.Vector3();
const overlaps = (a, b) => a[0] < b[2] && a[2] > b[0] && a[1] < b[3] && a[3] > b[1];

// entries: [{ object, text, wanted, always, left }]
export function declutter(entries, camera, width, height) {
  camera.updateMatrixWorld();
  const placed = [];
  for (const e of entries) {
    if (!e.wanted) { e.object.visible = false; continue; }
    v.setFromMatrixPosition(e.object.matrixWorld).project(camera);
    if (v.z < -1 || v.z > 1) { e.object.visible = false; continue; }
    const x = ((v.x + 1) / 2) * width, y = ((1 - v.y) / 2) * height;
    const w = e.text.length * 7 + 16;
    const rect = e.left ? [x - w, y - 9, x, y + 9] : [x - 4, y - 9, x + w, y + 9];
    const clash = placed.some((r) => overlaps(rect, r));
    e.object.visible = e.always || !clash;
    if (e.object.visible) placed.push(rect);
  }
}

export class LabelPlacer {
  constructor(ctx) {
    this.ctx = ctx;
    this.rects = [];
  }

  reserve(x, y, r) {
    this.rects.push([x - r, y - r, x + r, y + r]);
  }

  // Draws text near (x, y), outside a symbol of radius r. Returns false if no
  // position was free and the label was not forced.
  place(text, x, y, r, { force = false, centred = false } = {}) {
    const ctx = this.ctx;
    const w = ctx.measureText(text).width;
    const h = 12;
    const candidates = centred
      ? [[x - w / 2, y + 4]]
      : [[x + r + 4, y - r + 2], [x + r + 4, y + r + h - 2], [x - r - 4 - w, y - r + 2],
        [x - r - 4 - w, y + r + h - 2], [x - w / 2, y + r + h + 2]];
    for (const [cx, cy] of candidates) {
      const rect = [cx - 1, cy - h, cx + w + 1, cy + 3];
      if (!this.rects.some((q) => overlaps(rect, q))) {
        ctx.fillText(text, cx, cy);
        this.rects.push(rect);
        return true;
      }
    }
    if (force) {
      const [cx, cy] = candidates[0];
      ctx.fillText(text, cx, cy);
      this.rects.push([cx - 1, cy - h, cx + w + 1, cy + 3]);
      return true;
    }
    return false;
  }
}
