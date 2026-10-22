// A column of the Earth's crust for the place card: its layers to scale, from
// the surface down past the Moho, each with its P-wave speed, the Moho marked.

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const COLOURS = {
  water: '#4f78b8', ice: '#e9f0f7',
  'upper sediments': '#d9bf8c', 'middle sediments': '#c7a774', 'lower sediments': '#b38f60',
  'upper crystalline crust': '#c99aa6', 'middle crystalline crust': '#a7788b', 'lower crystalline crust': '#86596f',
  mantle: '#4d5b3f',
};

// c: from column(); returns an SVG string.
export function crustDiagram(c) {
  const W = 280, top = 10, left = 46, barW = 46, H = 210;
  const z0 = -c.surface;                       // depth of the surface (negative above sea level), km
  const zMax = c.moho + Math.max(8, c.moho * 0.15);
  const y = (z) => top + ((z - Math.min(0, z0)) / (zMax - Math.min(0, z0))) * (H - top - 14);
  let svg = '';
  for (const l of c.layers) {
    const y0 = y(-l.top), y1 = y(-l.bottom);
    svg += `<rect x="${left}" y="${y0.toFixed(1)}" width="${barW}" height="${Math.max(0.6, y1 - y0).toFixed(1)}" fill="${COLOURS[l.name]}"/>`;
    if (y1 - y0 >= 9) {
      svg += `<text x="${left + barW + 8}" y="${((y0 + y1) / 2 + 3.5).toFixed(1)}" class="layer">${esc(l.name)}<tspan class="v"> · ${l.vp.toFixed(1)} km/s</tspan></text>`;
    }
  }
  // The mantle below the Moho.
  svg += `<rect x="${left}" y="${y(c.moho).toFixed(1)}" width="${barW}" height="${(y(zMax) - y(c.moho)).toFixed(1)}" fill="${COLOURS.mantle}"/>`;
  svg += `<text x="${left + barW + 8}" y="${(y(c.moho) + 13).toFixed(1)}" class="layer">mantle<tspan class="v"> · Pn ${c.mantle.vp.toFixed(2)} km/s</tspan></text>`;
  svg += `<line x1="${left - 4}" x2="${left + barW + 4}" y1="${y(c.moho).toFixed(1)}" y2="${y(c.moho).toFixed(1)}" stroke="#f6c47e" stroke-width="1.6"/>`;
  svg += `<text x="${left - 7}" y="${(y(c.moho) + 3.5).toFixed(1)}" class="moho">Moho</text>`;
  // Depth scale, km below sea level.
  const step = zMax > 60 ? 20 : zMax > 25 ? 10 : 5;
  for (let z = 0; z <= zMax; z += step) {
    svg += `<line x1="${left - 3}" x2="${left}" y1="${y(z).toFixed(1)}" y2="${y(z).toFixed(1)}" class="tick"/><text x="${left - 6}" y="${(y(z) + 3.5).toFixed(1)}" class="depth">${z}</text>`;
  }
  svg += `<line x1="${left - 3}" x2="${left + barW + 3}" y1="${y(0).toFixed(1)}" y2="${y(0).toFixed(1)}" class="sea"/>`;
  return `<svg class="crust-diagram" viewBox="0 0 ${W} ${H}" role="img" aria-label="The crust beneath this place">${svg}<text x="${left - 6}" y="${H - 1}" class="depth">km</text></svg>`;
}
