// A stereonet of an earthquake's mechanism for the place card: the lower
// hemisphere in equal-area projection (the Schmidt net), north up, with the
// compressional quadrants coloured, the two nodal planes as great circles, and
// the P (pressure) and T (tension) axes.

import { fromEqualArea, radiation, planeTrace, toEqualArea } from '../core/moment.js';

// The beach ball as an image (a data URL), painted pixel by pixel.
function ballImage(m, colour, px) {
  const c = document.createElement('canvas');
  c.width = c.height = px;
  const ctx = c.getContext('2d');
  const img = ctx.createImageData(px, px);
  const [cr, cg, cb] = [1, 3, 5].map((i) => parseInt(colour.slice(i, i + 2), 16));
  const pale = [251, 241, 234];
  for (let j = 0; j < px; j++) {
    for (let i = 0; i < px; i++) {
      const x = ((i + 0.5) / px) * 2 - 1, y = 1 - ((j + 0.5) / px) * 2;
      const d = fromEqualArea(x, y);
      const o = 4 * (j * px + i);
      if (!d) continue;
      const push = radiation(m, d) > 0;
      const rgb = push ? [cr, cg, cb] : pale;
      img.data[o] = rgb[0]; img.data[o + 1] = rgb[1]; img.data[o + 2] = rgb[2]; img.data[o + 3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);
  return c.toDataURL('image/png');
}

// mech: from mechanism(); m: the tensor (NED); colour: the compressional colour.
export function stereonet(mech, m, colour) {
  const R = 64, C = 78, W = 156;
  const at = ([x, y]) => [C + x * R, C - y * R];
  const path = (pts) => pts.map((p, i) => `${i ? 'L' : 'M'}${at(p)[0].toFixed(1)},${at(p)[1].toFixed(1)}`).join('');
  const planes = mech.planes.map((pl) => `<path d="${path(planeTrace(pl))}" fill="none" stroke="#2b1d2d" stroke-width="1.3"/>`).join('');
  const axis = (name, a) => {
    const rad = (x) => (x * Math.PI) / 180;
    const dir = [Math.cos(rad(a.plunge)) * Math.cos(rad(a.azimuth)), Math.cos(rad(a.plunge)) * Math.sin(rad(a.azimuth)), Math.sin(rad(a.plunge))];
    const [x, y] = at(toEqualArea(dir));
    return `<circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="3.2" fill="#fff" stroke="#2b1d2d" stroke-width="1"/>`
      + `<text x="${(x + 5).toFixed(1)}" y="${(y - 4).toFixed(1)}" class="axis">${name}</text>`;
  };
  return `<svg class="stereonet" viewBox="0 0 ${W} ${W}" width="${W}" height="${W}" role="img" aria-label="The mechanism on a stereonet">
    <image href="${ballImage(m, colour, 2 * 2 * R)}" x="${C - R}" y="${C - R}" width="${2 * R}" height="${2 * R}"/>
    ${planes}
    <circle cx="${C}" cy="${C}" r="${R}" fill="none" stroke="#e9d6e4" stroke-width="1.2"/>
    <line x1="${C}" y1="${C - R - 6}" x2="${C}" y2="${C - R}" stroke="#e9d6e4"/>
    <text x="${C}" y="${C - R - 8}" class="north">N</text>
    ${axis('P', mech.axes.P)}${axis('T', mech.axes.T)}
  </svg>`;
}
