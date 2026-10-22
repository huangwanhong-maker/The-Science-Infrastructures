// Names on the map, as HTML text placed over the canvas, without overlaps.
//
// Items are given in order of priority. Each is placed at the screen position
// of its place, unless it would overlap an item already placed, or its place
// is out of view; then it is hidden. Elements are kept by key and reused, so
// moving the map only moves them.

const overlaps = (a, b) => a[0] < b[2] && a[2] > b[0] && a[1] < b[3] && a[3] > b[1];

export class Labels {
  constructor(container) {
    this.container = container;
    this.elements = new Map();
    this.canvas = document.createElement('canvas').getContext('2d');
    this.widths = new Map();
  }

  measure(text, font) {
    const key = `${font}|${text}`;
    let w = this.widths.get(key);
    if (w === undefined) {
      this.canvas.font = font;
      w = this.canvas.measureText(text).width;
      this.widths.set(key, w);
    }
    return w;
  }

  // items: [{ key, text, html, lon, lat, className, font, size, anchor: 'centre' | 'right', always, spacing, colour }]
  // (text is measured; html, if given, is shown instead)
  // screen(lon, lat) -> [x, y, visibility 0..1] in CSS pixels, or null.
  update(items, screen, width, height) {
    const placed = [];
    const seen = new Set();
    for (const item of items) {
      const s = screen(item.lon, item.lat);
      if (!s || s[2] <= 0) continue;
      const [x, y, vis] = s;
      if (x < -50 || y < -20 || x > width + 50 || y > height + 20) continue;
      const shown = item.upper ? item.text.toUpperCase() : item.text;
      const tw = this.measure(shown, item.font) + (item.tracking || 0) * shown.length;
      const h = item.size || 13;
      const rect = item.anchor === 'right'
        ? [x - 4, y - h * 0.6, x + tw + 12, y + h * 0.6]
        : [x - tw / 2 - 3, y - h * 0.62, x + tw / 2 + 3, y + h * 0.62];
      if (!item.always && placed.some((r) => overlaps(rect, r))) continue;
      placed.push(rect);
      seen.add(item.key);
      let el = this.elements.get(item.key);
      if (!el) {
        el = document.createElement('div');
        this.container.appendChild(el);
        this.elements.set(item.key, el);
      }
      const content = item.html || item.text;
      if (el._text !== content) { if (item.html) el.innerHTML = item.html; else el.textContent = item.text; el._text = content; }
      const cls = `label ${item.className}${item.anchor === 'right' ? ' right' : ''}${item.upper ? ' upper' : ''}`;
      if (el.className !== cls) el.className = cls;
      if ((el.style.color || undefined) !== item.colour) el.style.color = item.colour || '';
      el.style.transform = `translate(${Math.round(x)}px, ${Math.round(y)}px)`;
      el.style.opacity = vis < 1 ? vis.toFixed(2) : '';
      el.hidden = false;
    }
    for (const [key, el] of this.elements) {
      if (!seen.has(key)) el.hidden = true;
    }
  }

  clear() {
    for (const el of this.elements.values()) el.remove();
    this.elements.clear();
  }
}
