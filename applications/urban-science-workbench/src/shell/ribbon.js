// The clock along the page's foot: the day itself. A band whose colour at each hour is the light of the real Sun at this
// place on this date (the plum of night, the rose of its rising and setting, the pale lilac of day), with the Sun's
// altitude drawn over it and a cursor at the clock's hour. Dragging along it sets the hour. Where a layer holds a record
// over years, a second scale of years stands above it (createYears).

import { h, svg } from '../dom.js';
import { instantOf, DAYS } from '../clock.js';
import { subsolarPoint, sunAltAz } from '/domains/astronomy/js/index.js';
import { ribbonColour } from './look.js';

const MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
const WIND = ['north', 'north-east', 'east', 'south-east', 'south', 'south-west', 'west', 'north-west'];

export function createRibbon(el, { theme }) {
  let clock = null, playing = null, profile = null, profileKey = '';
  const canvas = h('canvas', {});
  const slider = h('input', { type: 'range', min: 0, max: 23.75, step: 0.25, value: 12, 'aria-label': 'the hour' });
  const day = h('div', { class: 'day' }, canvas, slider);
  const play = h('button', { class: 'play', title: 'let the day run' }, '▶');
  const dateInput = h('input', { type: 'date', 'aria-label': 'the date' });
  const dateText = h('b', {});
  const zone = h('span', {});
  const when = h('div', { class: 'when' }, play, h('div', { class: 'date-box', style: 'position:relative' }, dateText, zone, dateInput));
  const time = h('div', { class: 'time' });
  const sunLine = h('div', { class: 'sun' });
  el.replaceChildren(when, day, h('div', { class: 'now' }, time, sunLine));

  const altitudeAt = (date, hour) => {
    const a = clock.area;
    return sunAltAz(a.lon, a.lat, subsolarPoint(instantOf(date, hour, a.timezone))).altitude;
  };

  function draw() {
    if (!clock) return;
    const r = day.getBoundingClientRect();
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    const W = Math.max(10, Math.round(r.width)), H = Math.max(10, Math.round(r.height));
    if (canvas.width !== W * dpr || canvas.height !== H * dpr) { canvas.width = W * dpr; canvas.height = H * dpr; }
    const g = canvas.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, W, H);
    const key = `${clock.area.id}|${clock.date}`;
    if (key !== profileKey) {
      profileKey = key;
      profile = Array.from({ length: 97 }, (_, i) => altitudeAt(clock.date, i / 4));
    }
    const band = H - 13, cs = getComputedStyle(document.documentElement);
    const ink = cs.getPropertyValue('--ink').trim(), faint = cs.getPropertyValue('--ink-faint').trim(), accent = cs.getPropertyValue('--accent').trim();
    const altAt = (u) => { const q = u * 96, i = Math.min(95, Math.floor(q)); return profile[i] + (profile[i + 1] - profile[i]) * (q - i); };
    for (let x = 0; x < W; x++) { g.fillStyle = ribbonColour(theme(), altAt(x / W)); g.fillRect(x, 0, 1.5, band); }
    // The Sun's path: its altitude over the band, the horizon at the middle of what it spans.
    const hi = Math.max(20, ...profile), lo = Math.min(-20, ...profile);
    const Y = (a) => band - 2 - ((a - lo) / (hi - lo)) * (band - 4);
    g.beginPath();
    for (let x = 0; x <= W; x += 2) { const y = Y(altAt(x / W)); if (x === 0) g.moveTo(x, y); else g.lineTo(x, y); }
    g.strokeStyle = 'rgba(30, 14, 40, 0.5)';
    g.lineWidth = 1;
    g.stroke();
    // The hours.
    g.font = "400 11px 'IBM Plex Mono', monospace";
    g.textBaseline = 'top';
    g.textAlign = 'center';
    for (let hr = 0; hr <= 24; hr += 3) {
      const x = (hr / 24) * W;
      g.fillStyle = faint;
      g.fillRect(Math.min(W - 1, x), band, 1, 3);
      if (hr > 0 && hr < 24) g.fillText(String(hr).padStart(2, '0'), x, band + 2);
    }
    // The cursor.
    const u = clock.hour / 24, cx = u * W;
    g.fillStyle = accent;
    g.fillRect(cx - 0.5, 0, 1, band + 3);
    g.save();
    g.translate(cx, Y(altAt(u)));
    g.rotate(Math.PI / 4);
    g.fillStyle = ink;
    g.strokeStyle = accent;
    g.fillRect(-2.5, -2.5, 5, 5);
    g.strokeRect(-2.5, -2.5, 5, 5);
    g.restore();
  }

  function show() {
    if (!clock) return;
    const [y, m, d] = clock.date.split('-').map(Number);
    dateText.textContent = `${DAYS[clock.weekday].slice(0, 3)} ${d} ${MONTHS[m - 1]} ${y}`;
    zone.textContent = clock.area.zoneNote || clock.area.timezone.replace(/_/g, ' ');
    if (dateInput.value !== clock.date) dateInput.value = clock.date;
    if (+slider.value !== clock.hour) slider.value = clock.hour;
    const hr = Math.floor(clock.hour), mi = Math.round((clock.hour - hr) * 60), off = clock.offset;
    time.replaceChildren(`${String(hr).padStart(2, '0')}:${String(mi).padStart(2, '0')}`, h('small', {}, `UTC${off >= 0 ? '+' : '−'}${Math.abs(off)}`));
    const s = clock.sun;
    sunLine.textContent = s.altitude > 0 ? `Sun ${Math.round(s.altitude)}° high, ${WIND[Math.round(((s.azimuth % 360) + 360) % 360 / 45) % 8]}` : `Sun ${Math.round(-s.altitude)}° below the horizon`;
    draw();
  }

  slider.addEventListener('input', () => clock && clock.set({ hour: +slider.value }));
  dateInput.addEventListener('change', () => clock && dateInput.value && clock.set({ date: dateInput.value }));
  dateText.parentElement.addEventListener('click', () => { try { dateInput.showPicker(); } catch { /* not offered */ } });
  play.addEventListener('click', () => {
    if (playing) { clearInterval(playing); playing = null; play.textContent = '▶'; return; }
    play.textContent = '❚❚';
    playing = setInterval(() => clock && clock.set({ hour: (clock.hour + 0.25) % 24 }), 400);
  });
  new ResizeObserver(draw).observe(day);

  return {
    bind(c) { clock = c; profileKey = ''; c.on(show); show(); },
    redraw() { show(); },
    stop() { if (playing) { clearInterval(playing); playing = null; play.textContent = '▶'; } },
  };
}

// A scale of years above the clock, for layers that hold a record over years. Only the last asked for is shown.
export function createYears() {
  let current = null, el = null;
  function close(owner) {
    if (!current || (owner && current !== owner)) return;
    if (current.timer) clearInterval(current.timer);
    current = null;
    if (el) { el.remove(); el = null; }
    document.body.classList.remove('has-years');
  }
  function open({ label, years, value, onChange }) {
    close();
    const first = years[0], last = years[years.length - 1];
    const me = { label };
    current = me;
    const big = h('b', {}, String(value));
    const input = h('input', { type: 'range', min: 0, max: years.length - 1, step: 1, value: Math.max(0, years.indexOf(value)), 'aria-label': label });
    const marks = svg('svg', {});
    const scrub = h('div', { class: 'scrub' }, marks, input);
    const play = h('button', { class: 'play', title: 'through the years', style: 'width:22px;height:22px;background:none;border:1px solid var(--hair-2);cursor:pointer;padding:0;font:400 11px var(--mono);color:var(--ink-2)' }, '▶');
    el = h('div', { class: 'years vellum' }, h('div', { style: 'display:flex;gap:8px;align-items:center' }, play, big, h('span', {}, label)), scrub);
    el.style.gridTemplateColumns = 'auto 1fr';
    document.body.append(el);
    document.body.classList.add('has-years');
    const draw = () => {
      const r = scrub.getBoundingClientRect();
      if (!r.width) return;
      marks.replaceChildren();
      const X = (k) => 4 + (r.width - 8) * (years.length > 1 ? k / (years.length - 1) : 0.5);
      const every = Math.ceil(years.length / Math.max(2, Math.floor(r.width / 58)));
      years.forEach((y, k) => {
        marks.append(svg('line', { x1: X(k), x2: X(k), y1: 3, y2: k % every === 0 ? 11 : 7 }));
        if (k % every === 0) marks.append(svg('text', { x: X(k), y: 23 }, String(y)));
      });
      marks.append(svg('line', { class: 'cursor', x1: X(+input.value), x2: X(+input.value), y1: 0, y2: 13 }));
    };
    const set = () => { const y = years[+input.value]; big.textContent = String(y); draw(); onChange(y); };
    input.addEventListener('input', set);
    play.addEventListener('click', () => {
      if (me.timer) { clearInterval(me.timer); me.timer = null; play.textContent = '▶'; return; }
      play.textContent = '❚❚';
      me.timer = setInterval(() => { input.value = (+input.value + 1) % years.length; set(); }, 650);
    });
    new ResizeObserver(draw).observe(scrub);
    draw();
    return { close: () => close(me), set(y) { input.value = Math.max(0, years.indexOf(y)); big.textContent = String(y); draw(); } };
  }
  return { open, close };
}
