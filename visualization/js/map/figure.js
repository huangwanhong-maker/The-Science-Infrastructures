// A map as a figure for publication: the drawing as shown, cloned with the page's computed styles
// written in, its moving arcs stilled and headed; composed with a title, a line of what is drawn, a
// legend (a ramp of colours with its ends, and keys of lines), and a line of source and link; the page's
// web fonts written in, so that it is set in the same type wherever it is opened. SVG, or PNG drawn from
// it at a multiple of its size. The figure is the map's own drawing, so the two cannot disagree.

const SVG = 'http://www.w3.org/2000/svg';
const KEEP = ['fill', 'fill-opacity', 'stroke', 'stroke-width', 'stroke-opacity', 'stroke-dasharray', 'stroke-linecap', 'stroke-linejoin',
  'paint-order', 'opacity', 'font-family', 'font-style', 'font-weight', 'font-variant', 'letter-spacing', 'vector-effect'];

// The SVG element cloned with its computed styles inlined and its classes dropped (in its definitions,
// those of patterns, whose fills are drawn); elements with the class `moving` (animated dashes) drawn
// whole, with the marker `head` at their end if there is one.
export function styledClone(svgEl, { moving = 'flowing', head = 'url(#head)' } = {}) {
  const copy = svgEl.cloneNode(true);
  const from = svgEl.querySelectorAll('*'), to = copy.querySelectorAll('*');
  from.forEach((node, n) => {
    if (!(node instanceof SVGGraphicsElement) || (node.closest('defs') && !node.closest('pattern'))) return;
    const s = getComputedStyle(node);
    const still = node.classList.contains(moving);
    to[n].setAttribute('style', KEEP.map((p) => `${p}:${p === 'stroke-dasharray' && still ? 'none' : s.getPropertyValue(p)}`).join(';'));
    if (still && head) to[n].setAttribute('marker-end', head);
    to[n].removeAttribute('class');
  });
  copy.removeAttribute('id');
  return copy;
}

// The page's Google Fonts, as @font-face rules with their files written in (the Latin subsets, at
// most four faces), or '' without the network.
let faces = null;
export function embeddedFonts(selector = 'link[href*="fonts.googleapis.com/css2"]') {
  faces ??= (async () => {
    const link = document.querySelector(selector);
    if (!link) return '';
    const css = await (await fetch(link.href)).text();
    const all = [...css.matchAll(/@font-face\s*{[^}]*}/g)].map((m) => m[0]);
    const latin = all.filter((f) => /unicode-range: U\+0000-00FF/.test(f) || !/unicode-range/.test(f));
    const out = [];
    for (const face of (latin.length ? latin : all).slice(0, 4)) {
      const url = face.match(/url\(([^)]+)\)/)?.[1];
      if (!url) continue;
      const bytes = new Uint8Array(await (await fetch(url)).arrayBuffer());
      let bin = '';
      for (let i = 0; i < bytes.length; i += 0x8000) bin += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
      out.push(face.replace(url, `data:font/woff2;base64,${btoa(bin)}`));
    }
    return out.join('\n');
  })().catch(() => '');
  return faces;
}

const el = (tag, attrs = {}, text) => {
  const n = document.createElementNS(SVG, tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== null && v !== undefined) n.setAttribute(k, v);
  if (text !== undefined) n.textContent = text;
  return n;
};

// The figure as SVG text. map: the SVG element shown; legend: {title, stops: [[offset, colour]], ends:
// [left, right], keys: [[colour, text]], hatched: [ground, line, text] (a swatch of the areas in
// dispute)}; footer: [left, right]; colours: {paper, ink, ink2, faint, border}.
export async function composeFigure({ map, title = '', subtitle = '', legend = null, footer = ['', ''], colours, width = 1800,
  font = "'Source Serif 4', Georgia, serif" }) {
  const vb = map.viewBox.baseVal;
  const pad = 48, top = subtitle ? 104 : 80, mapW = width - 2 * pad, mapH = Math.round((mapW * vb.height) / vb.width);
  const keyH = legend ? 150 : 60;
  const H = top + mapH + keyH;
  const text = (x, y, t, size, extra = {}) => el('text', { x, y, 'font-family': font, 'font-size': size, fill: colours.ink, ...extra }, t);
  const out = el('svg', { width, height: H, viewBox: `0 0 ${width} ${H}` });
  const fonts = await embeddedFonts();
  if (fonts) out.append(el('style', {}, fonts));
  const defs = el('defs');
  out.append(defs, el('rect', { width, height: H, fill: colours.paper }));
  if (title) out.append(text(pad, 56, title, 30, { 'font-weight': 600 }));
  if (subtitle) out.append(text(pad, 86, subtitle, 17, { fill: colours.ink2 }));
  const inner = styledClone(map);
  for (const [k, v] of Object.entries({ x: pad, y: top, width: mapW, height: mapH })) inner.setAttribute(k, v);
  out.append(inner);
  if (legend) {
    const ly = top + mapH + 30;
    const grad = el('linearGradient', { id: 'figure-ramp' });
    for (const [o, c] of legend.stops || []) grad.append(el('stop', { offset: o, 'stop-color': c }));
    defs.append(grad);
    if (legend.stops?.length) out.append(el('rect', { x: pad, y: ly, width: 360, height: 12, fill: 'url(#figure-ramp)' }));
    out.append(text(pad, ly - 10, legend.title || '', 15, { 'font-weight': 500 }));
    if (legend.ends) out.append(text(pad, ly + 32, legend.ends[0], 14, { fill: colours.ink2 }), text(pad + 360, ly + 32, legend.ends[1], 14, { fill: colours.ink2, 'text-anchor': 'end' }));
    (legend.keys || []).forEach(([c, t], k) => out.append(
      el('path', { d: `M${pad + 440} ${ly + 6 + 24 * k}h40`, stroke: c, 'stroke-width': 3, 'stroke-linecap': 'round' }),
      text(pad + 492, ly + 11 + 24 * k, t, 14, { fill: colours.ink2 })));
    if (legend.hatched) {
      const [ground, line, t] = legend.hatched;
      const hatch = el('pattern', { id: 'figure-hatch', patternUnits: 'userSpaceOnUse', width: 6, height: 6, patternTransform: 'rotate(45)' });
      hatch.append(el('rect', { width: 6, height: 6, fill: ground }), el('path', { d: 'M0 0V6', stroke: line, 'stroke-width': 1.6 }));
      defs.append(hatch);
      out.append(el('rect', { x: pad, y: ly + 50, width: 40, height: 14, fill: 'url(#figure-hatch)', stroke: colours.border, 'stroke-width': 1 }),
        text(pad + 52, ly + 62, t, 14, { fill: colours.ink2 }));
    }
  }
  out.append(el('path', { d: `M${pad} ${H - 40}H${width - pad}`, stroke: colours.border, 'stroke-width': 1 }),
    text(pad, H - 18, footer[0] || '', 13, { fill: colours.faint }),
    text(width - pad, H - 18, footer[1] || '', 13, { fill: colours.faint, 'text-anchor': 'end' }));
  return `<?xml version="1.0" encoding="utf-8"?>\n${new XMLSerializer().serializeToString(out)}\n`;
}

// The figure drawn into a PNG at `scale` times its size.
export function figurePng(svgText, scale = 2) {
  return new Promise((resolve, reject) => {
    const url = URL.createObjectURL(new Blob([svgText], { type: 'image/svg+xml' }));
    const img = new Image();
    img.onload = () => {
      const c = document.createElement('canvas');
      c.width = img.width * scale;
      c.height = img.height * scale;
      const g = c.getContext('2d');
      g.scale(scale, scale);
      g.drawImage(img, 0, 0);
      URL.revokeObjectURL(url);
      c.toBlob((b) => (b ? resolve(b) : reject(new Error('the figure could not be drawn'))), 'image/png');
    };
    img.onerror = () => reject(new Error('the figure could not be read'));
    img.src = url;
  });
}
