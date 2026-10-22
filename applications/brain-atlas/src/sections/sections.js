// Sections through the brain: the template's T1 image in a warm ramp, every
// structure of the joined labels washed in its colour with crisp boundaries,
// and nothing outside the brain, so that the stage shows through.
//
// A section is drawn at `scale` pixels per voxel: the image bilinearly, the labels
// by their nearest voxel, and a boundary wherever two neighbouring pixels belong
// to different structures. The same drawing serves the section views and the cap
// of a cut in three dimensions.

const PLANES = {
  // For each plane: the axis held fixed, and the voxel axes running across (u) and up (v).
  axial: { axis: 2, u: 0, v: 1, name: 'Axial', letter: 'z' },
  coronal: { axis: 1, u: 0, v: 2, name: 'Coronal', letter: 'y' },
  sagittal: { axis: 0, u: 1, v: 2, name: 'Sagittal', letter: 'x' },
};
export { PLANES };

// The warm ramp of the T1 image: plum shadows, rose half-tones, ivory lights.
const RAMP = [[0, [20, 13, 22]], [0.18, [46, 30, 46]], [0.45, [128, 92, 112]], [0.72, [214, 186, 190]], [1, [252, 242, 232]]];
const ramp = new Uint8Array(256 * 3);
for (let i = 0; i < 256; i++) {
  const t = i / 255;
  let k = 0;
  while (k < RAMP.length - 2 && t > RAMP[k + 1][0]) k++;
  const [t0, c0] = RAMP[k], [t1, c1] = RAMP[k + 1];
  const f = (t - t0) / (t1 - t0);
  for (let c = 0; c < 3; c++) ramp[3 * i + c] = Math.round(c0[c] + (c1[c] - c0[c]) * f);
}

export class Sections {
  // t1: {dims, data (uint8 0..255 display values), affine}; labels: {dims, data (uint16), affine}, on the same grid.
  constructor({ t1, labels, palette }) {
    this.t1 = t1;
    this.labels = labels;
    this.dims = t1.dims;
    // The grid's axes run along the world's (RAS), with any voxel size: the human template's
    // 1 mm, the mouse's 0.05 mm.
    this.origin = [t1.affine[3], t1.affine[7], t1.affine[11]];
    this.spacing = [t1.affine[0], t1.affine[5], t1.affine[10]];
    this.palette = palette;          // Uint8Array, 4 per label value: r, g, b, a (overlay strength)
    this.overlay = 0.5;
    this.selected = null;            // Set of label values lit
    this.dimOthers = false;
  }

  voxelOf(p) { return [0, 1, 2].map((a) => (p[a] - this.origin[a]) / this.spacing[a]); }
  worldOf(ijk) { return [0, 1, 2].map((a) => ijk[a] * this.spacing[a] + this.origin[a]); }

  labelAt(p) {
    const [i, j, k] = this.voxelOf(p).map(Math.round);
    const [nx, ny, nz] = this.dims;
    if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz) return 0;
    return this.labels.data[i + nx * (j + ny * k)];
  }

  // The size in voxels of a plane's image: [width, height].
  size(plane) { const p = PLANES[plane]; return [this.dims[p.u], this.dims[p.v]]; }

  // Draw a section into a canvas. position: the world coordinate (mm) of the plane on its
  // axis. flipU mirrors across (the sagittal view seen from the left has anterior on the left).
  draw(canvas, plane, position, { scale = 3, flipU = false, boundaries = true, background = false } = {}) {
    const p = PLANES[plane];
    const [nx, ny, nz] = this.dims;
    const W = this.dims[p.u], H = this.dims[p.v];
    const fixed = Math.round((position - this.origin[p.axis]) / this.spacing[p.axis]);
    const w = W * scale, h = H * scale;
    if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
    const ctx = canvas.getContext('2d');
    const img = ctx.createImageData(w, h);
    const out = img.data;
    const stride = [1, nx, nx * ny];
    const t1 = this.t1.data, lab = this.labels.data, pal = this.palette;
    const inside = fixed >= 0 && fixed < this.dims[p.axis];
    if (!inside) { ctx.putImageData(img, 0, 0); return { width: W, height: H, fixed }; }
    const base = fixed * stride[p.axis];
    const su = stride[p.u], sv = stride[p.v];
    // Label of each output pixel, for the boundaries.
    const lab2 = new Uint16Array(w * h);
    for (let py = 0; py < h; py++) {
      const vf = (H - 1) - ((py + 0.5) / scale - 0.5);      // image top is the plane's highest v
      const v0 = Math.max(0, Math.min(H - 2, Math.floor(vf))), fv = Math.min(1, Math.max(0, vf - v0));
      const vn = Math.max(0, Math.min(H - 1, Math.round(vf)));
      for (let px = 0; px < w; px++) {
        const uRaw = (px + 0.5) / scale - 0.5;
        const uf = flipU ? (W - 1) - uRaw : uRaw;
        const u0 = Math.max(0, Math.min(W - 2, Math.floor(uf))), fu = Math.min(1, Math.max(0, uf - u0));
        const un = Math.max(0, Math.min(W - 1, Math.round(uf)));
        const a = base + u0 * su + v0 * sv;
        const g = (t1[a] * (1 - fu) + t1[a + su] * fu) * (1 - fv) + (t1[a + sv] * (1 - fu) + t1[a + su + sv] * fu) * fv;
        const L = lab[base + un * su + vn * sv];
        lab2[py * w + px] = L;
        const o = 4 * (py * w + px);
        const gi = Math.round(g);
        let r = ramp[3 * gi], gg = ramp[3 * gi + 1], b = ramp[3 * gi + 2];
        let alpha = 255;
        if (L) {
          let s = (pal[4 * L + 3] / 255) * this.overlay;
          const lit = this.selected && this.selected.has(L);
          if (lit) s = Math.min(0.9, s + 0.38);
          else if (this.selected && this.dimOthers) s *= 0.35;
          const lum = 0.35 + 0.65 * (gi / 255);
          r = r * (1 - s) + pal[4 * L] * lum * s;
          gg = gg * (1 - s) + pal[4 * L + 1] * lum * s;
          b = b * (1 - s) + pal[4 * L + 2] * lum * s;
        } else if (gi < 8 && !background) {
          alpha = Math.round(255 * (gi / 8));
        }
        out[o] = r; out[o + 1] = gg; out[o + 2] = b; out[o + 3] = alpha;
      }
    }
    if (boundaries) {
      for (let py = 0; py < h - 1; py++) {
        for (let px = 0; px < w - 1; px++) {
          const i = py * w + px, L = lab2[i];
          const R = lab2[i + 1], D = lab2[i + w];
          if (L === R && L === D) continue;
          const o = 4 * i;
          const lit = this.selected && (this.selected.has(L) || this.selected.has(R) || this.selected.has(D));
          if (lit) { out[o] = 246; out[o + 1] = 202; out[o + 2] = 132; out[o + 3] = 255; continue; }
          const k = (L && R && D) ? 0.62 : 0.5;
          out[o] = out[o] * k + 255 * (1 - k) * 0.35; out[o + 1] = out[o + 1] * k + 240 * (1 - k) * 0.35; out[o + 2] = out[o + 2] * k + 246 * (1 - k) * 0.35;
          out[o + 3] = Math.max(out[o + 3], 150);
        }
      }
    }
    ctx.putImageData(img, 0, 0);
    return { width: W, height: H, fixed };
  }

  // The world extent of a plane's image, for placing the cap of a cut: {u: [min, max], v: [min, max]} in mm.
  extent(plane) {
    const p = PLANES[plane];
    const s = this.spacing;
    return { u: [this.origin[p.u] - 0.5 * s[p.u], this.origin[p.u] + (this.dims[p.u] - 0.5) * s[p.u]], v: [this.origin[p.v] - 0.5 * s[p.v], this.origin[p.v] + (this.dims[p.v] - 0.5) * s[p.v]] };
  }
}

// A section view: a canvas fitted into its element, the section image drawn and the
// point marked; pointer and wheel events report world points.
export class SectionView {
  constructor(el, plane, sections, { flipU = false, onPoint, onHover, onSwap }) {
    this.el = el;
    this.plane = plane;
    this.sections = sections;
    this.flipU = flipU;
    this.canvas = el.querySelector('canvas');
    this.image = document.createElement('canvas');
    this.onPoint = onPoint; this.onHover = onHover; this.onSwap = onSwap;
    this.point = [0, 0, 0];
    this.bind();
    new ResizeObserver(() => this.render()).observe(el);
  }

  setPlane(plane, flipU) { this.plane = plane; this.flipU = flipU; this.el.dataset.plane = plane; this.redraw(); }

  // The mapping between the canvas (CSS pixels) and the plane's voxel coordinates.
  layout() {
    const r = this.el.getBoundingClientRect();
    const [W, H] = this.sections.size(this.plane);
    const s = Math.min((r.width - 24) / W, (r.height - 24) / H);
    const ox = (r.width - W * s) / 2, oy = (r.height - H * s) / 2;
    return { r, s, ox, oy, W, H };
  }

  toWorld(clientX, clientY) {
    const { r, s, ox, oy, W, H } = this.layout();
    const x = (clientX - r.left - ox) / s - 0.5, y = (clientY - r.top - oy) / s - 0.5;
    if (x < -0.5 || y < -0.5 || x > W - 0.5 || y > H - 0.5) return null;
    const p = PLANES[this.plane];
    const u = this.flipU ? (W - 1) - x : x, v = (H - 1) - y;
    const ijk = [0, 0, 0];
    ijk[p.u] = u; ijk[p.v] = v;
    ijk[p.axis] = (this.point[p.axis] - this.sections.origin[p.axis]) / this.sections.spacing[p.axis];
    const w = this.sections.worldOf(ijk);
    w[p.axis] = this.point[p.axis];
    return w;
  }

  toScreen(world) {
    const { s, ox, oy, W, H } = this.layout();
    const p = PLANES[this.plane];
    const ijk = this.sections.voxelOf(world);
    const u = this.flipU ? (W - 1) - ijk[p.u] : ijk[p.u];
    return { x: ox + (u + 0.5) * s, y: oy + ((H - 1) - ijk[p.v] + 0.5) * s };
  }

  bind() {
    let dragging = false;
    this.canvas.addEventListener('pointerdown', (e) => {
      const w = this.toWorld(e.clientX, e.clientY);
      if (!w) return;
      dragging = true;
      this.canvas.setPointerCapture(e.pointerId);
      this.onPoint(w, { click: true, plane: this.plane });
    });
    this.canvas.addEventListener('pointermove', (e) => {
      const w = this.toWorld(e.clientX, e.clientY);
      if (dragging && w) this.onPoint(w, { drag: true, plane: this.plane });
      else this.onHover(w, e);
    });
    this.canvas.addEventListener('pointerup', () => { dragging = false; });
    this.canvas.addEventListener('pointerleave', () => this.onHover(null));
    this.canvas.addEventListener('wheel', (e) => {
      e.preventDefault();
      const p = PLANES[this.plane];
      const q = this.point.slice();
      q[p.axis] += (e.deltaY > 0 ? -1 : 1) * (e.shiftKey ? 5 : 1) * this.sections.spacing[p.axis];
      this.onPoint(q, { wheel: true, plane: this.plane });
    }, { passive: false });
    const name = this.el.querySelector('.section-name');
    if (name) name.addEventListener('click', () => this.onSwap && this.onSwap(this.plane));
  }

  // Redraw the image (the slice changed), then the canvas.
  redraw() {
    const p = PLANES[this.plane];
    this.sections.draw(this.image, this.plane, this.point[p.axis], { scale: 3, flipU: this.flipU });
    this.render();
  }

  // Draw the image and the point into the canvas.
  render() {
    const { r, s, ox, oy, W, H } = this.layout();
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const c = this.canvas;
    if (c.width !== Math.round(r.width * dpr) || c.height !== Math.round(r.height * dpr)) { c.width = Math.round(r.width * dpr); c.height = Math.round(r.height * dpr); }
    const ctx = c.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, r.width, r.height);
    if (!this.image.width) return;
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(this.image, ox, oy, W * s, H * s);
    // The point: a gold cross with a gap around it.
    const q = this.toScreen(this.point);
    ctx.strokeStyle = 'rgba(242, 194, 124, 0.55)';
    ctx.lineWidth = 1;
    const gap = 9;
    ctx.beginPath();
    ctx.moveTo(ox, q.y); ctx.lineTo(q.x - gap, q.y); ctx.moveTo(q.x + gap, q.y); ctx.lineTo(ox + W * s, q.y);
    ctx.moveTo(q.x, oy); ctx.lineTo(q.x, q.y - gap); ctx.moveTo(q.x, q.y + gap); ctx.lineTo(q.x, oy + H * s);
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(q.x, q.y, 3.2, 0, 2 * Math.PI);
    ctx.fillStyle = '#f6c47e';
    ctx.fill();
    // Orientation letters.
    ctx.font = '600 11px "Source Serif 4", Georgia, serif';
    ctx.fillStyle = 'rgba(242, 194, 124, 0.75)';
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    const letters = { axial: ['L', 'R', 'A', 'P'], coronal: ['L', 'R', 'S', 'I'], sagittal: this.flipU ? ['A', 'P', 'S', 'I'] : ['P', 'A', 'S', 'I'] }[this.plane];
    ctx.fillText(letters[0], ox + 8, oy + H * s / 2);
    ctx.fillText(letters[1], ox + W * s - 8, oy + H * s / 2);
    ctx.fillText(letters[2], ox + W * s / 2, oy + 8);
    ctx.fillText(letters[3], ox + W * s / 2, oy + H * s - 8);
    const name = this.el.querySelector('.section-name'), pos = this.el.querySelector('.section-pos');
    if (name) name.textContent = PLANES[this.plane].name;
    if (pos) {
      const a = PLANES[this.plane].axis, fine = this.sections.spacing[a] < 0.5;
      pos.textContent = `${PLANES[this.plane].letter} = ${fine ? this.point[a].toFixed(2) : Math.round(this.point[a])} mm`;
    }
  }
}
