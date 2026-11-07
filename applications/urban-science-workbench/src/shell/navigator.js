// The navigator at the page's right: an altitude tape from the whole Earth to a doorstep, with a needle at the camera's
// range and the names of the scales, each a key that flies there; a compass whose ring turns the heading and whose
// drag tilts the view; a scale bar, read from the ground at the point of regard. The tape says where one is, and,
// beside each layer in the lenses, a small copy of it says where that layer means something.

import { h, svg } from '../dom.js';

export const SCALES = [
  { name: 'Earth', at: 2.0e7 }, { name: 'Continent', at: 5.0e6 }, { name: 'Region', at: 6.0e5 }, { name: 'Metropolis', at: 7.0e4 },
  { name: 'City', at: 1.2e4 }, { name: 'District', at: 2200 }, { name: 'Street', at: 380 }, { name: 'Doorstep', at: 70 },
];
const FAR = 6e7, NEAR = 30;
const DEG = Math.PI / 180;

// A range (m) to [0, 1] along the tape: 0 the farthest.
export const along = (range) => Math.min(1, Math.max(0, Math.log(FAR / range) / Math.log(FAR / NEAR)));
export const rangeAt = (u) => FAR * (NEAR / FAR) ** u;

export function scaleName(range) {
  let best = SCALES[0], d = Infinity;
  for (const s of SCALES) { const q = Math.abs(Math.log(s.at / range)); if (q < d) { d = q; best = s; } }
  return best.name;
}

export const distance = (m) => (m >= 1000 ? `${m >= 1e5 ? Math.round(m / 1000).toLocaleString('en') : (m / 1000).toFixed(m >= 1e4 ? 0 : 1)} km` : `${Math.round(m)} m`);

export function createNavigator({ ladderEl, compassEl, scaleEl, earth }) {
  let W = 118, H = 300;
  let last = {};
  const needle = svg('path', { class: 'needle' });
  const needleLine = svg('line', { class: 'needle-line' });
  const reading = svg('text', { class: 'reading' });
  const root = svg('svg', {});
  const labels = [];
  ladderEl.append(root);

  function draw() {
    const r = ladderEl.getBoundingClientRect();
    W = r.width || 118, H = r.height || 300;
    root.replaceChildren();
    labels.length = 0;
    const x = W - 3;
    const Y = (u) => 6 + (H - 12) * u;
    root.append(svg('line', { class: 'rule', x1: x, x2: x, y1: Y(0), y2: Y(1) }));
    for (let e = 2; e <= 8; e++) for (const m of [1, 2, 5]) {
      const q = m * 10 ** e;
      if (q > FAR || q < NEAR) continue;
      const y = Y(along(q));
      root.append(svg('line', { class: 't', x1: x - (m === 1 ? 7 : 3), x2: x, y1: y, y2: y }));
    }
    for (const s of SCALES) {
      const y = Y(along(s.at));
      root.append(svg('line', { class: 'rule', x1: x - 11, x2: x, y1: y, y2: y }));
      const t = svg('text', { x: x - 15, y: y + 4, onclick: (e) => { e.stopPropagation(); flyToScale(s); } }, s.name);
      t.dataset.scale = s.name;
      labels.push(t);
      root.append(t);
    }
    root.append(needleLine, needle, reading);
    last = {};
    update(earth.view.range);
  }

  function flyToScale(s) {
    const v = earth.view;
    const tilt = s.at > 3e5 ? 0 : s.at > 5e4 ? 20 * DEG : Math.max(v.tilt, 45 * DEG);
    earth.flyTo({ range: s.at, tilt }, 1.8);
  }

  function update(range) {
    if (range === last.range) return;
    last.range = range;
    const x = W - 3, u = along(range), y = 6 + (H - 12) * u;
    needle.setAttribute('d', `M${x - 1},${y} l-7,-3.6 v7.2 z`);
    needleLine.setAttribute('x1', x - 1), needleLine.setAttribute('x2', x + 1);
    needleLine.setAttribute('y1', y), needleLine.setAttribute('y2', y);
    const name = scaleName(range);
    for (const t of labels) t.classList.toggle('here', t.dataset.scale === name);
  }

  // Dragging along the tape sets the range.
  let dragging = false;
  const setFrom = (e) => {
    const r = ladderEl.getBoundingClientRect();
    earth.controls.flight = null;
    const v = rangeAt(Math.min(1, Math.max(0, (e.clientY - r.top - 6) / (r.height - 12))));
    earth.setView({ range: v });
  };
  ladderEl.addEventListener('pointerdown', (e) => { dragging = true; ladderEl.setPointerCapture(e.pointerId); setFrom(e); });
  ladderEl.addEventListener('pointermove', (e) => { if (dragging) setFrom(e); });
  ladderEl.addEventListener('pointerup', () => { dragging = false; });
  new ResizeObserver(draw).observe(ladderEl);

  // ---- the compass -------------------------------------------------------------------------------------------------
  const C = 42, R = 40;
  const rose = svg('g', {});
  for (let a = 0; a < 360; a += 30) {
    const major = a % 90 === 0, r0 = major ? R - 7 : R - 4;
    const s = Math.sin(a * DEG), c = Math.cos(a * DEG);
    rose.append(svg('line', { class: 'rose-t', x1: C + s * r0, y1: C - c * r0, x2: C + s * R, y2: C - c * R }));
  }
  rose.append(svg('path', { class: 'north', d: `M${C},${C - R + 1} l-3.4,8.5 h6.8 z` }));
  const tiltArc = svg('path', { class: 'tilt' });
  const horizon = svg('line', { class: 'horizon' });
  const compass = svg('svg', { viewBox: '0 0 84 84' },
    svg('circle', { class: 'ring', cx: C, cy: C, r: R }), rose, tiltArc, horizon,
    svg('circle', { class: 'hit', cx: C, cy: C, r: R }),
    svg('circle', { class: 'down', cx: C, cy: C, r: 6, onclick: (e) => { e.stopPropagation(); earth.flyTo({ tilt: 0, heading: 0 }, 0.9); } }, svg('title', {}, 'north up, looking straight down')));
  compassEl.append(compass);
  let turn = null;
  compass.addEventListener('pointerdown', (e) => {
    if (e.target.classList.contains('down')) return;
    compass.setPointerCapture(e.pointerId);
    earth.controls.flight = null;
    turn = { x: e.clientX, y: e.clientY, heading: earth.view.heading, tilt: earth.view.tilt };
  });
  compass.addEventListener('pointermove', (e) => {
    if (!turn) return;
    earth.setView({ heading: turn.heading - (e.clientX - turn.x) * 0.012, tilt: turn.tilt - (e.clientY - turn.y) * 0.012 });
  });
  compass.addEventListener('pointerup', () => { turn = null; });

  function updateCompass(v) {
    if (v.heading === last.heading && v.tilt === last.tilt) return;
    last.heading = v.heading, last.tilt = v.tilt;
    rose.setAttribute('transform', `rotate(${(-v.heading / DEG).toFixed(2)} ${C} ${C})`);
    // The tilt: a bar that falls from the centre's top as the view leans over from straight down.
    const t = v.tilt / (88 * DEG);
    const r = 24, a = 90 * t;
    const sx = Math.sin(a * DEG) * r, cy = Math.cos(a * DEG) * r;
    tiltArc.setAttribute('d', t > 0.01 ? `M${C - sx},${C + cy} A${r},${r} 0 0 0 ${C + sx},${C + cy}` : '');
    horizon.setAttribute('x1', C - 12), horizon.setAttribute('x2', C + 12);
    const y = C + 24 - 42 * t;
    horizon.setAttribute('y1', y), horizon.setAttribute('y2', y);
  }

  // ---- the scale bar -----------------------------------------------------------------------------------------------
  const barText = h('span', {}), bar = h('i', {});
  scaleEl.append(barText, bar);
  function updateScale() {
    const m = earth.groundMetresPerPixel();
    if (!Number.isFinite(m)) return;
    const target = m * 90, e = Math.floor(Math.log10(target)), f = target / 10 ** e;
    const nice = (f >= 5 ? 5 : f >= 2 ? 2 : 1) * 10 ** e;
    const px = Math.round(nice / m);
    if (px === last.px) return;
    last.px = px;
    barText.textContent = distance(nice);
    bar.style.width = `${px}px`;
  }

  return {
    update(v) {
      update(v.range);
      updateCompass(v);
      updateScale();
    },
    flyToScale,
  };
}
