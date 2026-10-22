// The climate diagram of Walter and Lieth (1960): a place's year in one figure.
//
// Months run along the axis, from January (or, south of the equator, from
// July, so that summer is in the middle). Temperature (left scale) and
// precipitation (right scale) are drawn so that 10 degrees C stand level with
// 20 mm: where the precipitation curve lies above the temperature curve the
// month is humid (hatched), where below it is dry (dotted), after Gaussen's
// rule P <= 2 T. Above 100 mm the precipitation scale is reduced tenfold, and
// the very wet season is filled. Months whose mean daily minimum is below 0 C,
// when frost is certain, are marked under the axis.
//
// Walter, H. and Lieth, H. (1960-1967), Klimadiagramm-Weltatlas, Jena: Fischer.

const MONTHS = 'JFMAMJJASOND';

// Precipitation (mm) in degrees of the temperature scale: 2 mm per degree up
// to 100 mm, then 20 mm per degree.
export const precipitationScale = (p) => (p <= 100 ? p / 2 : 50 + (p - 100) / 20);

// T, P, Tmin: twelve monthly values, January first. Returns an SVG string,
// whose height grows with the range drawn: a division of 10 C (20 mm, or 200 mm
// above 100 mm) is always `division` pixels, so that diagrams compare.
export function climateDiagram({ T, P, Tmin, lat, width = 288, division = 21 }) {
  const order = lat < 0 ? [6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5] : [...Array(12).keys()];
  const t = order.map((i) => T[i]), p = order.map((i) => precipitationScale(P[i]));
  const frost = order.map((i) => Tmin && Tmin[i] < 0);
  const lo = Math.min(0, Math.floor(Math.min(...t) / 10) * 10);
  const hi = Math.max(30, Math.ceil(Math.max(...t, ...p) / 10) * 10);
  // The temperature scale is labelled only as far as temperatures go.
  const tTop = Math.max(20, Math.ceil(Math.max(...t) / 10) * 10);
  const left = 26, right = 32, top = 18, bottom = 24;
  const w = width - left - right, h = ((hi - lo) / 10) * division, height = h + top + bottom;
  const X = (m) => left + (m / 12) * w;                 // m: months from the start of the first, 0..12
  const Y = (v) => top + ((hi - v) / (hi - lo)) * h;
  // The curves at month m (centres at m + 0.5), linear between centres and cyclic at the ends.
  const at = (a, m) => {
    const x = m - 0.5, i = Math.floor(x), f = x - i;
    return a[(i + 12) % 12] * (1 - f) + a[(i + 13) % 12] * f;
  };
  let svg = '';
  // Humid (hatched), very wet (filled) and dry (dotted) seasons.
  for (let px = 0.75; px < w; px += 2.5) {
    const m = (px / w) * 12, tv = at(t, m), pv = at(p, m), x = (left + px).toFixed(1);
    if (pv > tv) svg += `<line x1="${x}" x2="${x}" y1="${Y(tv).toFixed(1)}" y2="${Y(Math.min(pv, 50)).toFixed(1)}" class="wet"/>`;
    if (pv > 50) svg += `<line x1="${x}" x2="${x}" y1="${Y(50).toFixed(1)}" y2="${Y(pv).toFixed(1)}" class="very-wet"/>`;
  }
  for (let px = 1.5; px < w; px += 4) {
    const m = (px / w) * 12, tv = at(t, m), pv = at(p, m);
    for (let y = Y(pv) - 2; y > Y(tv) + 1; y -= 4) svg += `<circle cx="${(left + px).toFixed(1)}" cy="${y.toFixed(1)}" r="0.9" class="dry"/>`;
  }
  // Scales: 10 degrees to a division; 20 mm to a division, then 200 mm above 100 mm.
  const ticks = [];
  for (let v = lo; v <= hi; v += 10) ticks.push(v);
  svg += ticks.map((v) => `<line x1="${left}" x2="${left + w}" y1="${Y(v).toFixed(1)}" y2="${Y(v).toFixed(1)}" class="${v === 0 ? 'zero' : 'grid'}"/>`
    + (v <= tTop ? `<text x="${left - 4}" y="${(Y(v) + 3).toFixed(1)}" class="t-scale">${String(v).replace('-', '−')}</text>` : '')
    + (v >= 0 ? `<text x="${left + w + 4}" y="${(Y(v) + 3).toFixed(1)}" class="p-scale">${v <= 50 ? 2 * v : 100 + (v - 50) * 20}</text>` : '')).join('');
  // The curves, from the middle of the first month to the middle of the last, joined across the year's end.
  const curve = (a) => [0, ...a.map((_, i) => i + 0.5), 12].map((m, k) => `${k ? 'L' : 'M'}${X(m).toFixed(1)},${Y(at(a, m)).toFixed(1)}`).join('');
  svg += `<path d="${curve(p)}" class="p-line"/><path d="${curve(t)}" class="t-line"/>`;
  // Frost and the months.
  frost.forEach((f, i) => { if (f) svg += `<rect x="${(X(i) + 0.5).toFixed(1)}" y="${(Y(lo) + 2).toFixed(1)}" width="${(w / 12 - 1).toFixed(1)}" height="4" class="frost"/>`; });
  svg += order.map((i, k) => `<text x="${X(k + 0.5).toFixed(1)}" y="${(Y(lo) + 17).toFixed(1)}" class="month">${MONTHS[i]}</text>`).join('');
  svg += `<text x="${left - 4}" y="${(Y(tTop) - 9).toFixed(1)}" class="t-scale unit">°C</text><text x="${left + w + 4}" y="${top - 9}" class="p-scale unit">mm</text>`;
  return `<svg class="climate-diagram" viewBox="0 0 ${width} ${height}" role="img" aria-label="Climate diagram: monthly temperature and precipitation">${svg}</svg>`;
}
