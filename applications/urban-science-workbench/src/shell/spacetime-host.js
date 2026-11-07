// What is shown in space and time: the records a plug-in, a file of one's own or a script has given, each drawn on the map by the mark its form asks for, all on one
// line of time. A plug-in reaches it as ctx.st:
//
//   const shown = ctx.st.show(record, { mark, value, map, scale, opacity });    a record as spacetime.js sets out; refused with its reason if it cannot be read
//   ctx.st.records()     what is shown, in columns (spacetime.js, prepared)
//   ctx.st.now()         the moment (ms UTC), and ctx.st.window() how long a mark stays
//   shown.setStyle({...}); shown.remove();
//
// A record's marks stand in whichever frame serves it: on a city's own canvas among its buildings when it is about the open place at a city's scale, otherwise on the
// Earth's, relative to the eye. Which is decided again when the place changes or the camera leaves the city.

import { spacetime as ST, createBlooms, createThreads, localFrameAt, rampTexture } from '/visualization/js/gis/index.js';
import { cartesianDeg } from '/visualization/js/globe/figure.js';
import { colourMap, MAPS } from '/visualization/styles/index.js';
import { h } from '../dom.js';
import { createTimeline } from './timeline.js';

const DEG = Math.PI / 180;
const CITY_REACH = 300e3;        // a record across less than this can stand among a city's buildings
const MARKS = { events: ['blooms', 'columns'], tracks: ['threads'], flows: ['threads'], cells: ['columns'], grid: ['drape'] };
const MARK_NAMES = { blooms: 'blooms', columns: 'columns', threads: 'threads', drape: 'on the ground' };
const metres = (la0, lo0, la1, lo1) => Math.hypot((la0 - la1) * 111e3, (((lo0 - lo1 + 540) % 360) - 180) * 111e3 * Math.cos(la0 * DEG));

// Where a value falls: linearly, on a square root (counts, where a few cells hold most), on a logarithm where it spans decades, about zero where it changes sign.
function scaleFor([lo, hi], kind) {
  const signed = lo < 0 && hi > 0;
  const k = kind && kind !== 'auto' ? kind : signed ? 'diverging' : lo > 0 && hi / lo > 300 ? 'log' : 'linear';
  if (k === 'diverging') { const m = Math.max(-lo, hi) || 1; return { kind: k, u: (v) => 0.5 + v / (2 * m), low: -m, high: m }; }
  if (k === 'log') { const a = Math.log(Math.max(lo, 1e-12)), b = Math.log(hi); return { kind: k, u: (v) => (v > 0 ? (Math.log(v) - a) / (b - a) : 0), low: lo, high: hi }; }
  if (k === 'sqrt') { const a = Math.sqrt(Math.max(lo, 0)), b = Math.sqrt(Math.max(hi, 0)); return { kind: k, u: (v) => (b > a ? (Math.sqrt(Math.max(v, 0)) - a) / (b - a) : 0.5), low: lo, high: hi }; }
  return { kind: 'linear', u: (v) => (hi > lo ? (v - lo) / (hi - lo) : 0.5), low: lo, high: hi };
}
const fmt = (v) => { if (!Number.isFinite(v)) return '—'; const a = Math.abs(v); return a >= 1e6 ? `${(v / 1e6).toFixed(1)} m` : a >= 1e4 ? `${Math.round(v / 1e3)}k` : a >= 100 ? String(Math.round(v)) : a >= 1 ? v.toFixed(1) : a === 0 ? '0' : v.toPrecision(2); };

export function createSpaceTimeHost({ THREE, earth, gis, area, theme, cityFar, legendEl, openLegend }) {
  const shown = new Map();
  // A file dropped on the page, for whoever takes files; kept until someone does.
  let fileTaker = null, pendingFile = null;
  const acceptFiles = (take) => {
    fileTaker = take;
    if (pendingFile) { const f = pendingFile; pendingFile = null; take(f); }
    return () => { if (fileTaker === take) fileTaker = null; };
  };
  const dropped = (file) => { if (fileTaker) { fileTaker(file); return true; } pendingFile = file; return false; };
  let serial = 0;
  const timeline = createTimeline({
    onTime: (now, win, keep) => { for (const e of shown.values()) moveTo(e, now, win, keep); },
    onSun: (t) => {
      const c = area.clock?.();
      if (!c) return;
      const offset = c.offset ?? 0;
      const d = new Date(t + offset * 3600000);
      c.set({ date: d.toISOString().slice(0, 10), hour: d.getUTCHours() + d.getUTCMinutes() / 60 });
    },
  });

  // ---- colours -------------------------------------------------------------------------------------------------------------------------------------------

  // A map as drawn on this theme: from what stands off the ground up to its strongest, turned over where it would otherwise fade into the page.
  function toneOf(name) {
    const m = colourMap(name), night = theme() === 'night';
    const sequential = m.kind === 'sequential';
    const darkLow = name === 'plum';
    const flip = sequential && (night !== darkLow);
    return (u) => {
      const t = Math.min(1, Math.max(0, u));
      if (!sequential) return m.colour(t);
      const v = 0.14 + 0.86 * t;
      return m.colour(flip ? 1 - v : v);
    };
  }

  // ---- geometry -------------------------------------------------------------------------------------------------------------------------------------------

  function chooseScene(e) {
    const [s, w, n, ee] = e.rec.box;
    const across = metres(s, w, n, ee);
    const a = area.current();
    const c = [(s + n) / 2, (w + ee) / 2];
    return across < CITY_REACH && a && gis() && !cityFar() && metres(c[0], c[1], a.lat, a.lon) < 150e3 ? 'city' : 'earth';
  }
  function placer(where, centre) {
    if (where === 'city') {
      const c = earth.city, a = area.current();
      const frame = localFrameAt(a.lat, a.lon, c.height);
      const ground = (lat, lon) => (c.level ? c.height : earth.tree.heightAt(lon * DEG, lat * DEG) || c.height);
      return (lat, lon, up = 0) => { const p = frame.toLocal(lat, lon, ground(lat, lon) + up); return [p[0], p[2], -p[1]]; };
    }
    const fig = earth.fig;
    return (lat, lon, up = 0) => { const p = cartesianDeg(fig, lon, lat, up); return [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]]; };
  }
  const centreOf = (e) => { const [s, w, n, ee] = e.rec.box; return cartesianDeg(earth.fig, (w + ee) / 2, (s + n) / 2, 0); };
  function sphereOf(position) {
    let a = [Infinity, Infinity, Infinity], b = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < position.length; i += 3) for (let k = 0; k < 3; k++) { const v = position[i + k]; if (v < a[k]) a[k] = v; if (v > b[k]) b[k] = v; }
    return { centre: a.map((v, k) => (v + b[k]) / 2), radius: Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]) / 2 + 10 };
  }
  // The lift of marks off the ground, in whichever frame: above the terrain from orbit, just above the street in a city.
  const liftOf = (where, across) => (where === 'city' ? 2.5 : Math.max(1500, Math.min(14000, across * 0.004)));

  function valueColumn(e) {
    const v = e.style.value, r = v ? e.rec.ranges[v] : null;
    const scale = r ? scaleFor(r, e.style.scale) : null;
    return { name: v, scale, of: scale ? (x) => (Number.isFinite(x) ? Math.min(1, Math.max(0, scale.u(x))) : 0.5) : () => 0.72 };
  }

  function buildBlooms(e, where, centre, place) {
    const r = e.rec, n = r.count, position = new Float32Array(3 * n), time = new Float32Array(n), value = new Float32Array(n), size = new Float32Array(n);
    const [t0, t1] = r.span, T = Math.max(1, t1 - t0), V = valueColumn(e), col = V.name ? r.values[V.name] : null;
    const across = metres(r.box[0], r.box[1], r.box[2], r.box[3]), lift = liftOf(where, across);
    for (let k = 0; k < n; k++) {
      const p = place(r.lat[k], r.lon[k], lift);
      position[3 * k] = p[0], position[3 * k + 1] = p[1], position[3 * k + 2] = p[2];
      time[k] = (r.time[k] - t0) / T;
      const u = col ? V.of(col[k]) : 0.72;
      value[k] = u;
      size[k] = col ? u : 0.45;
    }
    return createBlooms(THREE, { position, time, value, size, sphere: sphereOf(position) }, { ramp: e.ramp, night: theme() === 'night', pixel: (e.style.pixel ?? 7) * Math.min(2, devicePixelRatio || 1) });
  }

  function buildThreads(e, where, centre, place) {
    const r = e.rec, [t0, t1] = r.span, T = Math.max(1, t1 - t0), V = valueColumn(e);
    const across = metres(r.box[0], r.box[1], r.box[2], r.box[3]), lift = liftOf(where, across);
    const pos = [], time = [], value = [];
    const push = (p, t, u) => { pos.push(p[0], p[1], p[2]); time.push((t - t0) / T); value.push(u); };
    if (r.kind === 'tracks') {
      for (const tr of r.tracks) {
        const col = V.name ? tr.values[V.name] : null;
        let prev = place(tr.lat[0], tr.lon[0], lift);
        for (let k = 1; k < tr.time.length; k++) {
          const p = place(tr.lat[k], tr.lon[k], lift);
          push(prev, tr.time[k - 1], col ? V.of(col[k - 1]) : 0.72);
          push(p, tr.time[k], col ? V.of(col[k]) : 0.72);
          prev = p;
        }
      }
    } else {
      const P = r.places, last = T / 60;
      for (const f of r.flows) {
        const a = P.get(f.from), b = P.get(f.to), dur = f.until ? f.until - f.time : last, d = metres(a.lat, a.lon, b.lat, b.lon), steps = 18;
        const u = V.of(f.value);
        let prev = null;
        for (let k = 0; k <= steps; k++) {
          const s = k / steps;
          // Along the great circle, near enough for a drawing, risen by a share of its length.
          const dlon = ((b.lon - a.lon + 540) % 360) - 180;
          const p = place(a.lat + (b.lat - a.lat) * s, a.lon + dlon * s, lift + d * 0.14 * 4 * s * (1 - s));
          if (prev) { push(prev.p, prev.t, u); push(p, f.time + s * dur, u); }
          prev = { p, t: f.time + s * dur };
        }
      }
    }
    const position = Float32Array.from(pos);
    return createThreads(THREE, { position, time: Float32Array.from(time), value: Float32Array.from(value), sphere: sphereOf(position) }, { ramp: e.ramp, night: theme() === 'night', through: where === 'city' ? 0.34 : 0, width: where === 'city' ? 2.2 : 1.1 });
  }

  // A small hexagon about a point, for what has a place and no shape of its own.
  const hexAbout = (lat, lon, radius) => Array.from({ length: 6 }, (_, k) => {
    const a = (k * Math.PI) / 3;
    return [lon + ((radius * Math.cos(a)) / 111e3) / Math.cos(lat * DEG), lat + (radius * Math.sin(a)) / 111e3];
  });

  function buildColumns(e, where) {
    const r = e.rec;
    const rings = r.kind === 'cells' ? r.cells.map((c) => c.ring) : Array.from({ length: r.count }, (_, k) => hexAbout(r.lat[k], r.lon[k], e.columnRadius ?? Math.max(30, metres(r.box[0], r.box[1], r.box[2], r.box[3]) / 240)));
    let size = 0;
    for (const q of rings.slice(0, 400)) size = Math.max(size, metres(q[0][1], q[0][0], q[Math.floor(q.length / 2)][1], q[Math.floor(q.length / 2)][0]));
    e.top = e.style.top ?? size * (where === 'city' ? 1.6 : 6);
    if (where === 'earth') {
      const set = earth.addColumns(e.id, rings.map((ring) => ({ ring })), { shrink: 0.88, opacity: e.style.opacity ?? 0.94 });
      return { prisms: set.prisms, earthSet: true };
    }
    const a = area.current(), frame = localFrameAt(a.lat, a.lon);
    const layer = gis().addPrisms(e.id, rings, { ring: (q) => q.map(([lon, lat]) => frame.toLocal(lat, lon).slice(0, 2)), height: () => -1, color: () => null, z: 0.6, opacity: e.style.opacity ?? 0.92, shrink: 0.88, order: 3, pick: false });
    return { prisms: layer.prisms };
  }

  // ---- what is drawn at a moment -----------------------------------------------------------------------------------------------------------------------

  function moveTo(e, now, win, keep) {
    const r = e.rec, T = Math.max(1, r.span[1] - r.span[0]);
    e.now = now, e.win = win, e.keep = keep;
    // Among a city's buildings, the way a thing has come stays plainer than it need across a continent.
    const kept = e.where === 'city' && e.draw !== 'columns' && r.kind !== 'events' ? Math.max(keep, 0.5) : keep;
    if (e.mark?.setTime) e.mark.setTime((now - r.span[0]) / T, Number.isFinite(win) ? win / T : 1e6, win === Infinity ? 0 : kept);
    if (e.draw === 'columns') writeColumns(e, now, win);
    if (e.draw === 'drape') writeDrape(e, now);
    if (e.counter && now - (e.counted || -Infinity) !== 0) count(e, now, win);
  }

  function count(e, now, win) {
    const r = e.rec;
    if (r.kind !== 'events') { e.counter.textContent = `${r.count.toLocaleString('en')} ${r.kind === 'tracks' ? 'tracks' : r.kind === 'cells' ? 'places' : r.kind === 'flows' ? 'flows' : 'cells'}`; return; }
    const b = ST.lowerBound(r.time, now + 1), a = Number.isFinite(win) ? ST.lowerBound(r.time, now - win) : 0;
    e.counter.textContent = `${(b - a).toLocaleString('en')} of ${r.count.toLocaleString('en')} now`;
  }

  // The step at a moment, and how far to the next.
  function stepAt(time, t) {
    const k = ST.lowerBound(time, t);
    if (k <= 0) return [0, 0];
    if (k >= time.length) return [time.length - 1, 0];
    return [k - 1, (t - time[k - 1]) / Math.max(1, time[k] - time[k - 1])];
  }

  function writeColumns(e, now, win) {
    const r = e.rec, n = r.count, V = valueColumn(e), lut = e.lut;
    const heights = e.h || (e.h = new Array(n)), colours = e.c || (e.c = new Array(n));
    if (r.kind === 'cells') {
      const col = V.name ? r.values[V.name] : null;
      if (!col) return;
      const [s, f] = stepAt(r.time, now), next = Math.min(r.time.length - 1, s + 1);
      for (let k = 0; k < n; k++) {
        const a = col[s * n + k], b = col[next * n + k];
        const v = Number.isFinite(a) && Number.isFinite(b) ? a + (b - a) * f : Number.isFinite(a) ? a : b;
        // A count of nothing raises nothing: a cell empty at this moment is not a low column.
        if (!Number.isFinite(v) || (v <= 0 && V.scale && V.scale.low >= 0)) { heights[k] = -1; colours[k] = null; continue; }
        const u = Math.min(1, Math.max(0, V.of(v)));
        heights[k] = e.top * (0.03 + u);
        colours[k] = lut[Math.round(u * 63)];
      }
    } else {
      // Events as columns: risen while they live, gone after.
      const col = V.name ? r.values[V.name] : null;
      for (let k = 0; k < n; k++) {
        const age = now - r.time[k];
        if (age < 0 || (Number.isFinite(win) && age > win)) { heights[k] = -1; colours[k] = null; continue; }
        const u = col ? V.of(col[k]) : 0.7, life = Number.isFinite(win) ? 1 - age / win : 1;
        heights[k] = e.top * (0.05 + u) * (0.25 + 0.75 * life);
        colours[k] = lut[Math.round(u * 63)];
      }
    }
    e.draws.prisms.set(heights, colours);
  }

  function writeDrape(e, now) {
    const r = e.rec, V = valueColumn(e);
    if (!V.name) return;
    const [s] = stepAt(r.time, now);
    if (e.step === s && e.drapeValue === V.name) return;
    e.step = s, e.drapeValue = V.name;
    const frames = e.frames || (e.frames = {});
    if (!frames[V.name]) frames[V.name] = r.values[V.name].map((f) => {
      const b = new Uint8Array(f.length);
      for (let k = 0; k < f.length; k++) b[k] = Number.isFinite(f[k]) ? 1 + Math.round(254 * Math.min(1, Math.max(0, V.of(f[k])))) : 0;
      return b;
    });
    if (!e.draped) {
      const tone = toneOf(e.style.map);
      earth.setDrape({ data: frames[V.name][s], width: r.width, height: r.height, box: r.box, ramp: (u) => { if (u <= 0.002) return [0, 0, 0, 0]; const c = tone((u * 255 - 1) / 254); return [c[0], c[1], c[2], Math.min(0.92, 1.7 * Math.sqrt(u))]; }, opacity: e.style.opacity ?? 0.92, floor: 0.002 });
      e.draped = true;
    } else earth.setDrapeData(frames[V.name][s]);
  }

  // ---- marks on and off -----------------------------------------------------------------------------------------------------------------------------------

  function clear(e) {
    if (e.mark) {
      if (e.where === 'earth') earth.removeObject(e.id); else gis()?.remove(e.id);
      e.mark = null;
    }
    if (e.draw === 'columns') {
      if (e.draws?.earthSet) earth.removeColumns(e.id); else gis()?.remove(e.id);
      e.h = e.c = null;
    }
    if (e.draw === 'drape' && e.draped) { earth.setDrape(null); e.draped = false; e.step = undefined; }
    e.draw = null;
    e.draws = null;
  }

  function build(e) {
    clear(e);
    const where = chooseScene(e);
    e.where = where;
    const tone = toneOf(e.style.map);
    e.ramp?.dispose?.();
    e.ramp = rampTexture(THREE, tone);
    e.lut = Array.from({ length: 64 }, (_, k) => tone(k / 63));
    const mark = e.style.mark;
    if (mark === 'blooms' || mark === 'threads') {
      const centre = where === 'earth' ? centreOf(e) : null;
      const place = placer(where, centre);
      e.mark = mark === 'blooms' ? buildBlooms(e, where, centre, place) : buildThreads(e, where, centre, place);
      if (e.style.opacity !== undefined) e.mark.setOpacity(e.style.opacity);
      if (where === 'earth') earth.addObject(e.id, e.mark.object, centre);
      else gis().addObject(e.id, e.mark.object);
    } else if (mark === 'columns') {
      e.draw = 'columns';
      e.draws = buildColumns(e, where);
    } else if (mark === 'drape') {
      e.draw = 'drape';
    }
    moveTo(e, timeline.now, timeline.window, timeline.keep);
  }

  function retime() {
    if (!shown.size) { timeline.close(); return; }
    let a = Infinity, b = -Infinity;
    for (const e of shown.values()) { a = Math.min(a, e.rec.span[0]); b = Math.max(b, e.rec.span[1]); }
    if (b <= a) b = a + ST.HOUR;
    const all = [...shown.values()];
    const bars = new Float32Array(180);
    for (const e of all) { const part = ST.histogram(e.rec.times, [a, b], 180); for (let k = 0; k < 180; k++) bars[k] += part[k]; }
    const asked = all[all.length - 1].style.window_ms;
    const newest = all[all.length - 1].style;
    timeline.set({ span: [a, b], bars, window: asked, keep: all.some((e) => e.rec.kind === 'events' || e.rec.kind === 'tracks') ? 0.16 : 0, zone: newest.zone || null, weekly: !!newest.weekly });
  }

  // ---- its card ---------------------------------------------------------------------------------------------------------------------------------------------

  function cardOf(e) {
    const r = e.rec;
    const ends = h('div', { class: 'ramp-ends' }, h('span', {}), h('span', {}));
    const ramp = h('span', { class: 'ramp' });
    e.counter = h('span', { class: 'st-count' });
    const paint = () => {
      const tone = toneOf(e.style.map), stops = Array.from({ length: 10 }, (_, i) => { const c = tone(i / 9); return `rgb(${c.map((v) => Math.round(v * 255)).join(',')})`; });
      ramp.style.background = `linear-gradient(90deg, ${stops.join(',')})`;
      const V = valueColumn(e);
      ramp.style.display = ends.style.display = V.name ? '' : 'none';
      if (V.scale) { ends.children[0].textContent = fmt(V.scale.low); ends.children[1].textContent = `${fmt(V.scale.high)}${r.units[V.name] ? ` ${r.units[V.name]}` : ''}`; }
    };
    const variable = r.variables.length ? h('select', { class: 'st-select', 'aria-label': 'the value shown', onchange: () => { e.style.value = variable.value || null; e.frames = null; e.draped = false; build(e); paint(); } },
      r.kind === 'events' || r.kind === 'tracks' ? h('option', { value: '' }, 'one colour') : null, r.variables.map((v) => h('option', { value: v, selected: v === e.style.value ? true : null }, v))) : null;
    const maps = Object.entries(MAPS).filter(([, m]) => m.kind === 'sequential' || m.kind === 'diverging').map(([k]) => k);
    const map = h('select', { class: 'st-select', 'aria-label': 'its colours', onchange: () => { e.style.map = map.value; e.frames = null; e.draped = false; build(e); paint(); } }, maps.map((m) => h('option', { value: m, selected: m === e.style.map ? true : null }, m)));
    const marks = MARKS[r.kind];
    const markKeys = marks.length > 1 && r.count <= 40000 ? h('div', { class: 'keys-tight' }, marks.map((m) => h('button', { class: m === e.style.mark ? 'on' : '', onclick: (ev) => { e.style.mark = m; for (const b of ev.currentTarget.parentElement.children) b.classList.toggle('on', b === ev.currentTarget); build(e); } }, MARK_NAMES[m]))) : null;
    const leftOut = Object.entries(r.left || {}).map(([why, n]) => `${n.toLocaleString('en')} left out: ${why}`).join('; ');
    const card = h('div', { class: 'card st-card' }, h('b', {}, r.title),
      h('button', { class: 'st-remove', title: 'take it away', onclick: () => remove(e.id) }, '×'),
      h('div', { class: 'st-line' }, e.counter),
      h('div', { class: 'st-keys' }, markKeys, variable, map), ramp, ends,
      r.note ? h('div', { class: 'ld-quiet' }, r.note) : null,
      leftOut ? h('div', { class: 'ld-quiet' }, leftOut) : null,
      r.report?.assumedUTC ? h('div', { class: 'ld-quiet' }, `${r.report.assumedUTC.toLocaleString('en')} times had no zone, and are read as UTC.`) : null,
      r.show?.credit || e.style.credit ? h('div', { class: 'ld-quiet' }, r.show?.credit || e.style.credit) : null);
    card.dataset.global = '1';
    paint();
    e.repaint = paint;
    return card;
  }

  // ---- reading what is under the pointer --------------------------------------------------------------------------------------------------------------

  function find(lat, lon) {
    const tolerance = Math.max(25, earth.groundMetresPerPixel() * 9);
    const when = (t) => ST.formatTime(t, 'full', 1);
    for (const e of [...shown.values()].reverse()) {
      const r = e.rec, now = timeline.now, win = timeline.window;
      if (r.kind === 'events') {
        const b = ST.lowerBound(r.time, now + 1), a = Number.isFinite(win) ? ST.lowerBound(r.time, now - win) : Math.max(0, b - 200000);
        let best = -1, d = tolerance;
        for (let k = a; k < b; k++) { if (Math.abs(r.lat[k] - lat) * 111e3 > d) continue; const q = metres(lat, lon, r.lat[k], r.lon[k]); if (q < d) { d = q; best = k; } }
        if (best >= 0) return { e, title: (r.labels && r.labels[best]) || r.title, lines: [new Date(r.time[best]).toISOString().slice(0, 16).replace('T', ' ') + ' UTC', ...r.variables.map((v) => `${v}: ${fmt(r.values[v][best])}${r.units[v] ? ` ${r.units[v]}` : ''}`)] };
      } else if (r.kind === 'tracks') {
        let best = null, d = tolerance;
        for (const t of r.tracks) { if (t.time[0] > now) continue; for (let k = 0; k < t.time.length; k++) { if (t.time[k] > now) break; const q = metres(lat, lon, t.lat[k], t.lon[k]); if (q < d) { d = q; best = { t, k }; } } }
        if (best) return { e, title: best.t.name || String(best.t.id), lines: [new Date(best.t.time[best.k]).toISOString().slice(0, 16).replace('T', ' ') + ' UTC', ...Object.keys(best.t.values).map((v) => `${v}: ${fmt(best.t.values[v][best.k])}${r.units[v] ? ` ${r.units[v]}` : ''}`)] };
      } else if (r.kind === 'cells') {
        const n = r.count, V = valueColumn(e), [s] = stepAt(r.time, now);
        for (let k = 0; k < n; k++) {
          const q = r.cells[k].ring;
          let inside = false;
          for (let i = 0, j = q.length - 1; i < q.length; j = i++) if ((q[i][1] > lat) !== (q[j][1] > lat) && lon < ((q[j][0] - q[i][0]) * (lat - q[i][1])) / (q[j][1] - q[i][1]) + q[i][0]) inside = !inside;
          if (inside) return { e, title: r.cells[k].name || String(r.cells[k].id), lines: [when(r.time[s]) , ...(V.name ? [`${V.name}: ${fmt(r.values[V.name][s * n + k])}${r.units[V.name] ? ` ${r.units[V.name]}` : ''}`] : [])] };
        }
      } else if (r.kind === 'grid') {
        const [south, west, north, east] = r.box;
        if (lat >= south && lat < north && lon >= west && lon < east && e.style.value) {
          const [s] = stepAt(r.time, now), i = Math.floor(((lon - west) / (east - west)) * r.width), j = Math.floor(((lat - south) / (north - south)) * r.height), v = r.values[e.style.value][s][j * r.width + i];
          if (Number.isFinite(v) && v !== 0) return { e, title: r.title, lines: [`${e.style.value}: ${fmt(v)}${r.units[e.style.value] ? ` ${r.units[e.style.value]}` : ''}`] };
        }
      }
    }
    return null;
  }

  // ---- the interface ---------------------------------------------------------------------------------------------------------------------------------------

  function show(record, style = {}) {
    const rec = ST.prepared(record);
    const id = `st:${++serial}`;
    const asked = { ...(rec.show || {}), ...style };
    const mark = MARKS[rec.kind].includes(asked.mark) ? asked.mark : MARKS[rec.kind][0];
    const value = asked.value && rec.variables.includes(asked.value) ? asked.value : (rec.kind === 'cells' || rec.kind === 'grid' ? rec.variables[0] : rec.variables.find((v) => /mag|mw|value|count|speed|wind|intensity|amount/i.test(v)) ?? rec.variables[0]) ?? null;
    const range = value ? rec.ranges[value] : null;
    const map = asked.map && MAPS[asked.map] ? asked.map : range && range[0] < 0 && range[1] > 0 ? 'lavender-rust' : 'plum';
    const e = { id, rec, raw: record, style: { ...asked, mark, value, map } };
    shown.set(id, e);
    e.card = cardOf(e);
    legendEl().append(e.card);
    openLegend(e.card);
    retime();
    build(e);
    return {
      id, record: rec,
      remove: () => remove(id),
      setStyle(changes) { Object.assign(e.style, changes); e.frames = null; e.draped = false; build(e); e.repaint?.(); },
    };
  }
  function remove(id) {
    const e = shown.get(id);
    if (!e) return;
    clear(e);
    e.ramp?.dispose?.();
    e.card?.remove();
    shown.delete(id);
    retime();
    for (const f of listeners) f();
  }
  const listeners = [];

  return {
    show, remove, acceptFiles, dropped,
    records: () => [...shown.values()].map((e) => ({ id: e.id, ...e.rec })),
    entries: () => [...shown.values()],
    raw: () => [...shown.values()].map((e) => e.raw),
    now: () => timeline.now,
    window: () => timeline.window,
    timeline,
    onChange: (f) => listeners.push(f),
    fit(id) {
      const e = id ? shown.get(id) : [...shown.values()].pop();
      if (!e) return;
      const [s, w, n, ee] = e.rec.box, across = metres(s, w, n, ee);
      earth.flyTo({ lat: ((s + n) / 2) * DEG, lon: ((w + ee) / 2) * DEG, height: 0, range: Math.min(2.2e7, Math.max(2500, across * 1.5)), tilt: across > 3e5 ? 0 : 50 * DEG }, 2.6);
    },
    // The place or the height changed: each stands where it should now.
    sceneChanged() {
      for (const e of shown.values()) {
        // What stood in the old city went with its scene.
        if (e.where === 'city') { e.mark = null; e.draw = null; e.draws = null; }
        build(e);
        e.repaint?.();
      }
    },
    needsRebuild() { for (const e of shown.values()) if (chooseScene(e) !== e.where) return true; return false; },
    hoverAt(lat, lon) { const f = find(lat, lon); return f ? [f.title, ...f.lines].join('\n') : null; },
    pickAt(lat, lon) { const f = find(lat, lon); return f ? { kind: 'spacetime', title: f.title, lines: f.lines, at: { lat, lon } } : null; },
  };
}

export const ST_DAY = ST.DAY;
