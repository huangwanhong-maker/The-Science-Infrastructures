// The loom: the history complex drawn with time across. Each being is a thread (a lane) from its
// first occasion to its last; each occasion a knot binding the threads of those present; each tie a
// ribbon of lenses between two threads, one lens (a 2-cell) between each two of their meetings,
// pinched at the knots: faint in its discipline's colour, hatched where the tie is negative, marked
// where it flips. Frustrated triangles of ties are marked above their knots in vermilion. Time is
// zoomed with the wheel and panned by dragging; a click on the scale of years sets the date.

import { svg, clear } from './dom.js';
import { COLOURS, colourOf } from './palette.js';
import { laneLayout, yearTicks, timeSpan } from './layout.js';
import { shortName, dateText, isoOfYear, yOfDate, yearName } from './format.js';

const LEFT = 176, TOP = 40, RIGHT = 18, BOTTOM = 26;

export class Loom {
  constructor(host, { onPick, onHover, onDate }) {
    this.host = host;
    this.onPick = onPick;
    this.onHover = onHover;
    this.onDate = onDate;
    this.root = svg('svg', { class: 'loom-svg' });
    host.append(this.root);
    this.domain = null;
    this.data = null;
    this.drag = null;
    this.root.addEventListener('wheel', (e) => this.wheel(e), { passive: false });
    this.root.addEventListener('pointerdown', (e) => this.down(e));
    this.root.addEventListener('pointermove', (e) => this.move(e));
    this.root.addEventListener('pointerup', (e) => this.up(e));
    this.root.addEventListener('pointerleave', () => { this.drag = null; this.onHover(null); });
    this.root.addEventListener('dblclick', () => { this.domain = null; this.render(); });
    new ResizeObserver(() => this.render()).observe(host);
  }

  set(data, view) {
    if (data !== this.data) this.domain = null;
    this.data = data;
    this.view = view;
    this.render();
  }

  // ---- scales ----
  x(year) {
    const [t0, t1] = this.domain;
    return LEFT + ((year - t0) / (t1 - t0)) * this.width;
  }
  yearAt(px) {
    const [t0, t1] = this.domain;
    return t0 + ((px - LEFT) / this.width) * (t1 - t0);
  }

  // ---- interaction ----
  wheel(e) {
    if (!this.data) return;
    e.preventDefault();
    const at = this.yearAt(e.offsetX);
    const k = Math.exp(Math.sign(e.deltaY) * 0.18);
    const [t0, t1] = this.domain, [f0, f1] = this.full;
    let a = at - (at - t0) * k, b = at + (t1 - at) * k;
    if (b - a > f1 - f0) [a, b] = [f0, f1];
    if (b - a < 0.25) return;
    this.domain = [a, b];
    this.render();
  }
  down(e) {
    if (!this.data || e.button !== 0) return;
    this.drag = { x: e.clientX, domain: [...this.domain], moved: false };
  }
  move(e) {
    if (this.drag) {
      const dx = e.clientX - this.drag.x;
      if (Math.abs(dx) > 3) this.drag.moved = true;
      if (this.drag.moved) {
        const [t0, t1] = this.drag.domain, shift = -(dx / this.width) * (t1 - t0);
        this.domain = [t0 + shift, t1 + shift];
        this.root.classList.add('dragging');
        this.render();
      }
    }
  }
  up(e) {
    const drag = this.drag;
    this.drag = null;
    this.root.classList.remove('dragging');
    if (drag?.moved) {
      // The click that ends a drag picks nothing.
      this.justDragged = true;
      setTimeout(() => { this.justDragged = false; }, 0);
      return;
    }
    if (e.offsetY < TOP && e.offsetX > LEFT && this.data) this.onDate(isoOfYear(this.yearAt(e.offsetX)));
  }

  // ---- drawing ----
  render() {
    const data = this.data, root = clear(this.root);
    if (!data) return;
    const view = this.view;
    const W = Math.max(480, this.host.clientWidth);
    const lanes = laneLayout(data, view.hide, { top: TOP + 6 });
    const H = lanes.height + BOTTOM;
    this.width = W - LEFT - RIGHT;
    this.full = timeSpan(data);
    if (!this.domain) this.domain = [...this.full];
    root.setAttribute('width', W);
    root.setAttribute('height', H);
    root.setAttribute('viewBox', `0 0 ${W} ${H}`);
    this.lanes = lanes;

    const focus = view.focus, sel = this.selection(focus);
    const defs = svg('defs', {},
      svg('clipPath', { id: 'loom-clip' }, svg('rect', { x: LEFT, y: 0, width: this.width, height: H })),
      ...COLOURS.map((c) => svg('pattern', { id: `hatch-${c.id}`, width: 6, height: 6, patternUnits: 'userSpaceOnUse', patternTransform: 'rotate(45)' },
        svg('rect', { width: 6, height: 6, fill: c.hex, 'fill-opacity': 0.07 }),
        svg('line', { x1: 0, y1: 0, x2: 0, y2: 6, stroke: c.hex, 'stroke-width': 1.3, 'stroke-opacity': 0.55 }))));
    root.append(defs);

    // Groups: a faint band per discipline behind its lanes (its name is with the lane labels).
    const bands = svg('g', { class: 'bands' });
    for (const g of lanes.groups) bands.append(svg('rect', { x: 0, y: g.top - 4, width: W, height: g.bottom - g.top + 8, class: 'band' }));
    root.append(bands);

    // The scale of years, and the faint lines of years beneath everything.
    const [t0, t1] = this.domain, { ticks } = yearTicks(t0, t1, this.width);
    const grid = svg('g', { class: 'grid', 'clip-path': 'url(#loom-clip)' });
    const axis = svg('g', { class: 'axis' }, svg('rect', { x: LEFT, y: 0, width: this.width, height: TOP - 6, class: 'axis-hit' }));
    for (const t of ticks) {
      const x = this.x(t);
      grid.append(svg('line', { x1: x, x2: x, y1: TOP - 4, y2: H - BOTTOM + 6 }));
      if (x >= LEFT - 1 && x <= LEFT + this.width + 1) axis.append(svg('text', { x, y: TOP - 14, 'text-anchor': 'middle' }, yearName(t)));
    }
    root.append(grid, axis);

    const plot = svg('g', { 'clip-path': 'url(#loom-clip)' });
    root.append(plot);
    const yOf = (b) => lanes.y.get(b);
    const occ = data.occasions;

    // Lenses.
    const lensLayer = svg('g', { class: 'lenses' });
    data.complex.faces.forEach((f, i) => {
      const t = data.ties[f.tie], ya = yOf(t.a), yb = yOf(t.b);
      if (ya === undefined || yb === undefined || view.hide.has(t.colour)) return;
      const x0 = this.x(occ[f.from].year), x1 = this.x(occ[f.to].year);
      if (x1 < LEFT - 2 || x0 > LEFT + this.width + 2) return;
      const neg = f.signs[0] < 0 && f.signs[1] < 0;
      const c = colourOf(t.colour);
      const state = sel.lens(i, f, t);
      const r = svg('rect', {
        x: x0, y: Math.min(ya, yb), width: Math.max(1, x1 - x0), height: Math.abs(yb - ya),
        class: `lens ${state}`, fill: neg ? `url(#hatch-${c.id})` : c.hex, 'data-face': i,
      });
      // A lens across many lanes is a large rectangle: it is lightened as the root of its span,
      // so that the near ties, which are most of them, keep the eye.
      const thin = 1 / Math.sqrt(Math.max(1, Math.abs(yb - ya) / 19));
      if (!neg) r.setAttribute('fill-opacity', state === 'on' ? 0.3 : state === 'off' ? 0.025 : 0.11 * thin);
      else r.setAttribute('opacity', state === 'on' ? 1 : state === 'off' ? 0.25 : 0.8 * thin);
      r.addEventListener('pointerenter', () => this.onHover({ kind: 'lens', face: i }));
      r.addEventListener('click', (e) => { if (!this.dragged(e)) this.onPick({ kind: 'tie', id: t.id }); });
      lensLayer.append(r);
      if (f.flip && x1 - x0 > 7) {
        const cx = (x0 + x1) / 2, cy = (ya + yb) / 2, s = 3.6;
        lensLayer.append(svg('path', { d: `M${cx - s} ${cy}L${cx} ${cy - s}L${cx + s} ${cy}L${cx} ${cy + s}Z`, class: `flip ${state}` }));
      }
    });
    plot.append(lensLayer);

    // Threads.
    const threads = svg('g', { class: 'threads' });
    for (const b of lanes.order) {
      const os = data.beings[b].occasions;
      if (!os.length) continue;
      const y = yOf(b), xa = this.x(occ[os[0]].year), xb = this.x(occ[os[os.length - 1]].year);
      const state = sel.being(b);
      const line = svg('line', { x1: xa, x2: Math.max(xb, xa + 1), y1: y, y2: y, stroke: colourOf(data.beings[b].colour).hex, class: `thread ${state}` });
      line.addEventListener('pointerenter', () => this.onHover({ kind: 'being', being: b }));
      line.addEventListener('click', (e) => { if (!this.dragged(e)) this.onPick({ kind: 'being', id: data.beings[b].id }); });
      threads.append(line);
    }
    plot.append(threads);

    // The cursor of the chosen date.
    if (view.date) {
      const y = yOfDate(view.date);
      if (y !== null) {
        const x = this.x(y);
        root.append(svg('g', { class: 'cursor' },
          svg('line', { x1: x, x2: x, y1: TOP - 8, y2: H - BOTTOM + 6 }),
          svg('text', { x: x + 4, y: TOP - 2 }, dateText(view.date))));
      }
    }

    // Knots.
    const knots = svg('g', { class: 'knots' });
    const frustrated = new Map();
    for (const t of data.readings.triangles) if (t.product < 0) {
      const list = frustrated.get(t.occasion) || [];
      list.push(t);
      frustrated.set(t.occasion, list);
    }
    for (const o of data.complex.vertices) {
      const oc = occ[o];
      const here = oc.present.filter((p) => yOf(p.being) !== undefined);
      if (!here.length) continue;
      const ys = here.map((p) => yOf(p.being)), y0 = Math.min(...ys), y1 = Math.max(...ys);
      const xa = this.x(oc.year), xb = this.x(oc.yearLast);
      if (xb < LEFT - 8 || xa > LEFT + this.width + 8) continue;
      const state = sel.occasion(o);
      const g = svg('g', { class: `knot ${state}` });
      if (xb - xa > 3) g.append(svg('rect', { x: xa, y: y0 - 4, width: xb - xa, height: y1 - y0 + 8, rx: 3, class: 'span' }));
      g.append(svg('line', { x1: xa, x2: xa, y1: y0, y2: y1, class: 'stem' }));
      for (const p of here) {
        const c = colourOf(data.beings[p.being].colour).hex;
        g.append(svg('circle', { cx: xa, cy: yOf(p.being), r: state === 'on' ? 4.2 : 3.3, class: p.derived ? 'dot derived' : 'dot', fill: p.derived ? 'var(--paper)' : c, stroke: p.derived ? c : 'var(--paper)' }));
      }
      const fr = frustrated.get(o);
      if (fr) {
        const attested = fr.some((t) => t.attested), s = 4.6, ty = y0 - 10;
        g.append(svg('path', { d: `M${xa - s} ${ty + s * 0.8}L${xa} ${ty - s}L${xa + s} ${ty + s * 0.8}Z`, class: attested ? 'frustration' : 'frustration carried' }));
      }
      g.append(svg('line', { x1: xa, x2: xa, y1: y0 - 14, y2: y1 + 6, class: 'hit' }));
      g.addEventListener('pointerenter', () => this.onHover({ kind: 'occasion', occasion: o }));
      g.addEventListener('click', (e) => { if (!this.dragged(e)) this.onPick({ kind: 'occasion', id: oc.id }); });
      knots.append(g);
    }
    plot.append(knots);

    // A tie followed: its signs at each meeting, beside the knots.
    if (focus?.kind === 'tie') {
      const ti = data.ties.findIndex((t) => t.id === focus.id), t = data.ties[ti];
      if (t && yOf(t.a) !== undefined && yOf(t.b) !== undefined) {
        const marks = svg('g', { class: 'rungs' });
        for (const r of t.rungs) {
          const x = this.x(occ[r.occasion].year), y = (yOf(t.a) + yOf(t.b)) / 2;
          marks.append(svg('text', { x: x + 5, y: y + 4, class: r.attested ? 'rung' : 'rung carried' }, r.sign > 0 ? '+' : '−'));
        }
        plot.append(marks);
      }
    }

    // A loop of H₁ drawn: along each being's lane between two occasions, and across at the turns.
    if (view.loop !== null && view.loop !== undefined) {
      const loop = data.invariants.loops[view.loop];
      if (loop) {
        let d = '', prev = null;
        for (const s of loop.steps) {
          const y = yOf(s.being);
          if (y === undefined) { prev = null; continue; }
          const xs = this.x(occ[s.from].year), xe = this.x(occ[s.to].year);
          // Joined to the step before at the knot where the loop turns from one being to the next.
          d += prev && Math.abs(prev - xs) < 0.5 ? `L${xs} ${y}` : `M${xs} ${y}`;
          d += `L${xe} ${y}`;
          prev = xe;
        }
        plot.append(svg('path', { d, class: 'loop' }));
      }
    }

    // Lane labels, over a paper margin, outside the clip.
    const labels = svg('g', { class: 'labels' }, svg('rect', { x: 0, y: TOP - 6, width: LEFT - 6, height: H, class: 'margin' }));
    for (const g of lanes.groups) {
      const c = colourOf(g.colour);
      labels.append(svg('text', { x: 24, y: g.top - 9, class: 'group-name' }, c.name),
        svg('rect', { x: 14, y: g.top - 4, width: 3, height: g.bottom - g.top + 8, fill: c.hex, rx: 1.5 }));
    }
    for (const b of lanes.order) {
      const being = data.beings[b];
      const t = svg('text', { x: LEFT - 12, y: yOf(b) + 4, 'text-anchor': 'end', class: `lane-name ${sel.being(b)}` }, being.short);
      t.addEventListener('pointerenter', () => this.onHover({ kind: 'being', being: b }));
      t.addEventListener('click', () => this.onPick({ kind: 'being', id: being.id }));
      labels.append(svg('title', {}, being.name), t);
    }
    root.append(labels);
  }

  dragged() { return Boolean(this.justDragged); }

  // What a followed being, tie or occasion brings forward: 'on', 'off' (pushed back) or ''.
  selection(focus) {
    const data = this.data;
    const none = { lens: () => '', being: () => '', occasion: () => '' };
    if (!focus) return none;
    if (focus.kind === 'being') {
      const b = data.beings.findIndex((x) => x.id === focus.id);
      if (b < 0) return none;
      const mine = new Set(data.beings[b].occasions);
      return {
        lens: (i, f, t) => (t.a === b || t.b === b ? 'on' : 'off'),
        being: (x) => (x === b ? 'on' : ''),
        occasion: (o) => (mine.has(o) ? 'on' : 'off'),
      };
    }
    if (focus.kind === 'tie') {
      const ti = data.ties.findIndex((t) => t.id === focus.id);
      if (ti < 0) return none;
      const t = data.ties[ti], joint = new Set(t.rungs.map((r) => r.occasion));
      return {
        lens: (i, f) => (f.tie === ti ? 'on' : 'off'),
        being: (x) => (x === t.a || x === t.b ? 'on' : ''),
        occasion: (o) => (joint.has(o) ? 'on' : 'off'),
      };
    }
    if (focus.kind === 'occasion') {
      const o = data.occasions.findIndex((x) => x.id === focus.id);
      if (o < 0) return none;
      const here = new Set(data.occasions[o].present.map((p) => p.being));
      return {
        lens: (i, f) => (f.from === o || f.to === o ? 'on' : 'off'),
        being: (x) => (here.has(x) ? 'on' : ''),
        occasion: (x) => (x === o ? 'on' : 'off'),
      };
    }
    return none;
  }
}
