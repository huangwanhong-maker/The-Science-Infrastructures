// A scalp map as electroencephalographers draw it: the head seen from above, the nose at
// the top, flattened by the azimuthal equidistant projection about the vertex (the
// distance from the centre of the disc is the angle from the vertex), the values between
// the electrodes interpolated over the sphere, with contour lines, the electrodes as
// dots and the outline of the head, its nose and its ears.

const THETA_MAX = (125 * Math.PI) / 180;      // the disc reaches 125 degrees from the vertex
const RING = (Math.PI / 2) / THETA_MAX;       // the head's outline: the level of the ears

// The disc position (x right, y down; radius 1 at THETA_MAX) of a direction from the centre.
export function project(d) {
  const n = Math.hypot(d[0], d[1], d[2]) || 1;
  const z = d[2] / n, theta = Math.acos(Math.max(-1, Math.min(1, z)));
  const s = Math.hypot(d[0], d[1]) || 1e-12;
  const rho = theta / THETA_MAX;
  return [(rho * d[0]) / s, (-rho * d[1]) / s];
}

// field(direction) -> value; colour(value) -> [r, g, b] 0..255 or null; electrodes: [{name, direction, value}].
export function drawTopo(canvas, { field, colour, electrodes = [], contours = null, size = 240, title = '', unit = '' }) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const S = Math.round(size * dpr);
  if (canvas.width !== S) { canvas.width = S; canvas.height = S; }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.clearRect(0, 0, S, S);
  const R = S * 0.4, cx = S / 2, cy = S / 2 + S * 0.02;
  const maxRho = Math.min(1, Math.max(RING * 1.08, ...electrodes.map((e) => Math.hypot(...project(e.direction))) ) + 0.04);
  // The map, pixel by pixel, a quarter of the resolution, then scaled.
  const step = 1;
  const img = ctx.createImageData(S, S);
  const values = new Float32Array(Math.ceil(S / step) * Math.ceil(S / step)).fill(NaN);
  for (let py = 0; py < S; py += step) for (let px = 0; px < S; px += step) {
    const X = (px - cx) / R, Y = (py - cy) / R, rho = Math.hypot(X, Y);
    if (rho > maxRho) continue;
    const theta = rho * THETA_MAX, c = rho > 1e-9 ? X / rho : 1, s = rho > 1e-9 ? -Y / rho : 0;
    const d = [Math.sin(theta) * c, Math.sin(theta) * s, Math.cos(theta)];
    const v = field(d);
    values[(py / step) * Math.ceil(S / step) + px / step] = v;
  }
  const W = Math.ceil(S / step);
  const level = (v) => (contours ? Math.floor(v / contours) : 0);
  for (let py = 0; py < S; py++) for (let px = 0; px < S; px++) {
    const k = Math.floor(py / step) * W + Math.floor(px / step);
    const v = values[k];
    if (!Number.isFinite(v)) continue;
    const col = colour(v);
    if (!col) continue;
    let [r, g, b] = col;
    // A contour where the level changes against the neighbour to the right or below.
    if (contours && (px % step === 0 || py % step === 0)) {
      const right = values[k + 1], below = values[k + W];
      if ((Number.isFinite(right) && level(right) !== level(v)) || (Number.isFinite(below) && level(below) !== level(v))) { r *= 0.62; g *= 0.62; b *= 0.62; }
    }
    const o = 4 * (py * S + px);
    img.data[o] = r; img.data[o + 1] = g; img.data[o + 2] = b; img.data[o + 3] = 255;
  }
  ctx.putImageData(img, 0, 0);
  // The outline of the head, its nose and its ears.
  ctx.strokeStyle = 'rgba(246, 238, 245, 0.8)';
  ctx.lineWidth = 1.4 * dpr;
  ctx.beginPath(); ctx.arc(cx, cy, R * RING, 0, 2 * Math.PI); ctx.stroke();
  const hr = R * RING;
  ctx.beginPath(); ctx.moveTo(cx - 0.09 * hr, cy - 0.995 * hr); ctx.lineTo(cx, cy - 1.12 * hr); ctx.lineTo(cx + 0.09 * hr, cy - 0.995 * hr); ctx.stroke();
  for (const sgn of [-1, 1]) {
    ctx.beginPath();
    ctx.ellipse(cx + sgn * 1.04 * hr, cy, 0.06 * hr, 0.17 * hr, 0, 0, 2 * Math.PI);
    ctx.stroke();
  }
  // The electrodes.
  for (const e of electrodes) {
    const [X, Y] = project(e.direction);
    const x = cx + X * R, y = cy + Y * R;
    ctx.beginPath(); ctx.arc(x, y, 1.9 * dpr, 0, 2 * Math.PI);
    ctx.fillStyle = 'rgba(28, 20, 30, 0.85)'; ctx.fill();
    ctx.strokeStyle = 'rgba(255, 246, 232, 0.9)'; ctx.lineWidth = 0.8 * dpr; ctx.stroke();
  }
  if (title) {
    ctx.fillStyle = 'rgba(246, 238, 245, 0.75)';
    ctx.font = `italic ${11 * dpr}px "Source Serif 4", Georgia, serif`;
    ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    ctx.fillText(title, 4 * dpr, 2 * dpr);
    if (unit) { ctx.textAlign = 'right'; ctx.fillText(unit, S - 4 * dpr, 2 * dpr); }
  }
}
