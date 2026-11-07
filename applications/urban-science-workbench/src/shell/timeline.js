// The line of time of whatever is shown in space and time: its span with how much happens along it, a cursor to drag, how long a mark stays after it happens, and a
// key that lets it run. Moments are in UTC, and say so. The clock's day ribbon below keeps the hour of the place, and the Sun can be asked to follow this line.

import { h } from '../dom.js';
import { spacetime } from '/visualization/js/gis/index.js';

const { DAY, HOUR, YEAR, formatSpan, formatTime, timeTicks } = spacetime;
import { offsetHours } from '../clock.js';

const SPEEDS = [0.25, 1, 4];
const WINDOWS = (span) => [
  ['an instant', span / 240], ['a little', span / 60], ['a good while', span / 16], ['all so far', Infinity],
];
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

export function createTimeline({ onTime, onSun }) {
  let el = null, canvas = null, span = [0, 1], bars = null, now = 0, windowMs = 1, keep = 0, playing = 0, speed = 1, last = 0, sun = false, zone = null, weekly = false;
  // Moments read in UTC, or in a record's own zone where it names one.
  const shift = (t) => (zone ? t + offsetHours(t, zone) * HOUR : t);
  const WEEKDAYS = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'];
  const hm = (t) => new Date(t).toISOString().slice(11, 16);
  const parts = {};
  const X = (t, w) => 6 + ((t - span[0]) / Math.max(1, span[1] - span[0])) * (w - 12);

  function draw() {
    if (!canvas) return;
    const r = canvas.getBoundingClientRect();
    if (!r.width) return;
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    if (canvas.width !== Math.round(r.width * dpr)) { canvas.width = Math.round(r.width * dpr); canvas.height = Math.round(r.height * dpr); }
    const g = canvas.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, r.width, r.height);
    const base = r.height - 15, top = 3;
    const accent = css('--accent'), ink = css('--ink-faint'), ink2 = css('--ink-2');
    // How much happens: thin bars, those passed in the accent.
    if (bars) {
      const n = bars.length, most = Math.max(1, ...bars), w = (r.width - 12) / n;
      for (let k = 0; k < n; k++) {
        if (!bars[k]) continue;
        const v = Math.sqrt(bars[k] / most), x = 6 + k * w, t = span[0] + ((k + 0.5) / n) * (span[1] - span[0]);
        g.fillStyle = t <= now ? accent : ink;
        g.globalAlpha = t <= now && (now - t) <= windowMs ? 0.95 : 0.42;
        g.fillRect(x, base - v * (base - top), Math.max(1, w - 0.6), v * (base - top));
      }
      g.globalAlpha = 1;
    }
    g.strokeStyle = ink;
    g.lineWidth = 1;
    g.beginPath();
    g.moveTo(6, base + 0.5);
    g.lineTo(r.width - 6, base + 0.5);
    g.stroke();
    g.font = `400 11px ${css('--mono') || 'monospace'}`;
    g.fillStyle = ink2;
    g.textAlign = 'center';
    const back = shift(span[0]) - span[0];
    if (weekly) {
      for (let k = 0; k < 7; k++) {
        const x = X(span[0] + k * DAY, r.width);
        g.fillRect(x, base + 1, 1, 3);
        g.fillText(WEEKDAYS[new Date(span[0] + back + k * DAY).getUTCDay()], Math.min(r.width - 18, x + (r.width - 12) / 14), r.height - 1.5);
      }
    } else for (const { t, unit } of timeTicks(span[0] + back, span[1] + back, Math.max(2, Math.floor(r.width / 88)))) {
      const x = X(t - back, r.width);
      g.fillRect(x, base + 1, 1, 3);
      g.fillText(formatTime(t, unit), Math.min(r.width - 18, Math.max(18, x)), r.height - 1.5);
    }
    // What is alive now, and the moment.
    const x = X(now, r.width);
    if (Number.isFinite(windowMs)) {
      g.fillStyle = accent;
      g.globalAlpha = 0.13;
      const x0 = X(Math.max(span[0], now - windowMs), r.width);
      g.fillRect(x0, top, x - x0, base - top);
      g.globalAlpha = 1;
    }
    g.fillStyle = accent;
    g.fillRect(x - 0.5, top - 2, 1, base - top + 2);
    g.beginPath();
    g.moveTo(x, top - 3); g.lineTo(x + 3.5, top + 0.5); g.lineTo(x, top + 4); g.lineTo(x - 3.5, top + 0.5);
    g.closePath();
    g.fill();
  }

  function say() {
    if (!el) return;
    const at = shift(now);
    parts.now.textContent = weekly ? `${WEEKDAYS[new Date(at).getUTCDay()]} ${hm(at)}` : formatTime(at, 'full', span[1] - span[0]);
    parts.zone.textContent = weekly ? 'a mean week' : zone ? zone.split('/').pop().replace(/_/g, ' ') : 'UTC';
    parts.within.textContent = Number.isFinite(windowMs) ? `the last ${formatSpan(windowMs)}` : 'all so far';
    draw();
  }
  function tell() {
    onTime(now, windowMs, keep);
    say();
    if (sun && onSun) onSun(now);
  }
  function loop(at) {
    if (!playing) return;
    const dt = Math.min(0.1, (at - last) / 1000);
    last = at;
    // The whole line in about forty seconds, at the middle speed.
    now += (dt * speed * (span[1] - span[0])) / 40;
    if (now > span[1]) now = span[0];
    tell();
    playing = requestAnimationFrame(loop);
  }
  function setPlaying(on) {
    if (on && !playing) { last = performance.now(); playing = requestAnimationFrame(loop); }
    if (!on && playing) { cancelAnimationFrame(playing); playing = 0; }
    if (parts.play) parts.play.textContent = on ? '❚❚' : '▶';
  }
  function build() {
    canvas = h('canvas', { class: 'tl-line', 'aria-label': 'the line of time' });
    parts.now = h('b', { class: 'tl-now' });
    parts.within = h('span', { class: 'tl-within' });
    parts.zone = h('span', {});
    parts.play = h('button', { class: 'play', title: 'let it run', onclick: () => setPlaying(!playing) }, '▶');
    const speeds = h('div', { class: 'keys-tight tl-speed', title: 'how fast it runs' }, SPEEDS.map((s) => h('button', { class: s === speed ? 'on' : '', onclick: (e) => { speed = s; for (const b of e.currentTarget.parentElement.children) b.classList.toggle('on', b === e.currentTarget); } }, s === 0.25 ? 'slow' : s === 1 ? 'even' : 'fast')));
    parts.windows = h('div', { class: 'keys-tight tl-windows', title: 'how long a thing stays after it happens' });
    const sunKey = h('button', { class: 'tl-sun', title: 'let the Sun over the place follow this line', onclick: () => { sun = !sun; sunKey.classList.toggle('on', sun); if (sun && onSun) onSun(now); } }, 'the Sun follows');
    el = h('div', { class: 'timeline vellum' },
      h('div', { class: 'tl-left' }, parts.play, h('div', {}, parts.now, h('div', { class: 'tl-sub' }, parts.zone, ' · ', parts.within))),
      canvas,
      h('div', { class: 'tl-right' }, parts.windows, h('div', { class: 'tl-keys' }, speeds, sunKey)));
    document.body.append(el);
    document.body.classList.add('has-timeline');
    let dragging = false;
    const move = (e) => {
      const r = canvas.getBoundingClientRect();
      const u = Math.min(1, Math.max(0, (e.clientX - r.left - 6) / (r.width - 12)));
      now = span[0] + u * (span[1] - span[0]);
      tell();
    };
    canvas.addEventListener('pointerdown', (e) => { dragging = true; canvas.setPointerCapture(e.pointerId); setPlaying(false); move(e); });
    canvas.addEventListener('pointermove', (e) => { if (dragging) move(e); });
    canvas.addEventListener('pointerup', () => { dragging = false; });
    new ResizeObserver(draw).observe(canvas);
  }
  function windows(preferred) {
    const list = WINDOWS(span[1] - span[0]);
    const choose = preferred ?? list[2][1];
    windowMs = choose;
    parts.windows.replaceChildren(...list.map(([name, w]) => h('button', { class: w === choose || (preferred === undefined && name === 'a good while') ? 'on' : '', onclick: (e) => {
      windowMs = w;
      keep = w === Infinity ? 0.5 : keep && keep < 0.5 ? keep : 0.16;
      for (const b of e.currentTarget.parentElement.children) b.classList.toggle('on', b === e.currentTarget);
      tell();
    } }, name)));
  }

  return {
    // Open, or carry on over a new span: `bars` how much happens along it.
    set({ span: s, bars: b, window: w, keep: k = 0.16, now: at, zone: z = null, weekly: week = false }) {
      if (!el) build();
      zone = z, weekly = week;
      const fresh = span[0] !== s[0] || span[1] !== s[1];
      span = s;
      bars = b;
      keep = k;
      if (fresh || at !== undefined) {
        now = at ?? (span[0] + 0.02 * (span[1] - span[0]));
        windows(w);
      }
      tell();
    },
    close() {
      setPlaying(false);
      if (el) el.remove();
      el = canvas = null;
      document.body.classList.remove('has-timeline');
    },
    play: () => setPlaying(true),
    pause: () => setPlaying(false),
    seek(t) { now = Math.min(span[1], Math.max(span[0], t)); tell(); },
    get now() { return now; },
    get window() { return windowMs; },
    get keep() { return keep; },
    get playing() { return !!playing; },
    get isOpen() { return !!el; },
  };
}

export { DAY, HOUR, YEAR };
