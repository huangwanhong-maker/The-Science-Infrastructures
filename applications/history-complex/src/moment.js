// The moment: the network of beings and live ties at a date, the spin network of the boundary at
// that instant. Beings sit on a circle in the loom's order; a link is solid where the tie was last
// attested positive and dashed where negative, and fades with the meetings since it was last
// attested (|expected| = memory^k). Frustrated triangles are filled faintly in vermilion. Beneath,
// the frustration of the network at the end of each year: as last attested, and as expected with
// the persistence of each discipline.

import { el, svg, clear } from './dom.js';
import { colourOf, INK } from './palette.js';
import { laneLayout } from './layout.js';
import { shortName, dateText, number, plural, listText, isoOfYear, yOfDate, roundYears, yearName } from './format.js';

export class MomentView {
  constructor(host, { load, pick, setDate }) {
    this.host = host;
    this.load = load;
    this.pick = pick;
    this.setDate = setDate;
    this.cache = new Map();
    this.timer = null;
  }

  render(data, view) {
    this.data = data;
    this.view = view;
    const host = clear(this.host);
    const [y0, y1] = periodOf(data);
    const date = view.date || isoOfYear((y0 + y1) / 2).slice(0, 4) + '-12-31';
    const slider = el('input', { type: 'range', min: y0, max: y1, step: 1 / 12, value: yOfDate(date) ?? y0, 'aria-label': 'The date' });
    const input = el('input', { type: 'text', value: date, size: 11, 'aria-label': 'The date, as EDTF', class: 'date-input' });
    slider.addEventListener('input', () => { input.value = isoOfYear(Number(slider.value)); this.ask(input.value, false); });
    slider.addEventListener('change', () => this.setDate(input.value));
    input.addEventListener('change', () => { if (yOfDate(input.value) !== null) this.setDate(input.value); });
    this.circle = svg('svg', { class: 'moment-svg', viewBox: '-80 -90 800 780', role: 'img', 'aria-label': 'The network of the moment' });
    this.side = el('div', { class: 'moment-side' });
    this.series = el('div', { class: 'series' });
    host.append(
      el('div', { class: 'moment-head' },
        el('h2', { class: 'view-title', text: 'The moment' }),
        el('div', { class: 'date-control' }, el('span', { class: 'faint', text: yearName(Math.ceil(y0)) }), slider, el('span', { class: 'faint', text: yearName(Math.floor(y1)) }), input)),
      el('p', { class: 'lede narrow', text: 'The beings and the ties alive at one date: the spin network of the boundary of the history at that instant. A link is solid where its tie was last attested positive, dashed where negative, and fades with every meeting since it was last attested, by the persistence of its discipline. A frustrated triangle, whose three signs multiply to −1, is shaded.' }),
      el('div', { class: 'moment-body' }, el('div', { class: 'moment-figure' }, this.circle), this.side),
      this.series);
    this.drawSeries();
    this.ask(date, true);
  }

  ask(date, now) {
    clearTimeout(this.timer);
    const go = async () => {
      const key = date;
      let m = this.cache.get(key);
      if (!m) {
        try { m = await this.load(date); } catch (error) { this.side.replaceChildren(el('p', { class: 'error', text: error.message })); return; }
        this.cache.set(key, m);
      }
      this.draw(m, date);
    };
    if (now) go(); else this.timer = setTimeout(go, 160);
  }

  draw(m, date) {
    const data = this.data, g = clear(this.circle);
    const order = laneLayout(data, new Set()).order;
    const cx = 320, cy = 300, R = 210;
    const pos = new Map(order.map((b, i) => {
      const a = -Math.PI / 2 + (2 * Math.PI * i) / order.length;
      return [b, { x: cx + R * Math.cos(a), y: cy + R * Math.sin(a), a }];
    }));
    const live = new Set(m.nodes);
    // Frustrated triangles beneath the links.
    for (const t of m.triangles) if (t.product < 0) {
      const p = t.beings.map((b) => pos.get(b));
      g.append(svg('polygon', { points: p.map((q) => `${q.x},${q.y}`).join(' '), class: 'tri-fill' }));
    }
    // Links as chords bowed toward the centre.
    m.links.forEach((l, i) => {
      const a = pos.get(l.a), b = pos.get(l.b), tie = data.ties[l.tie];
      const post = m.posterior?.links?.[i];
      const certainty = post === undefined ? Math.abs(l.expected) : Math.abs(post);
      const mx = (a.x + b.x) / 2, my = (a.y + b.y) / 2, qx = mx + (cx - mx) * 0.35, qy = my + (cy - my) * 0.35;
      const path = svg('path', {
        d: `M${a.x} ${a.y}Q${qx} ${qy} ${b.x} ${b.y}`, class: l.sign < 0 ? 'chord neg' : 'chord',
        stroke: colourOf(tie.colour).hex, 'stroke-opacity': 0.2 + 0.8 * certainty,
      });
      path.append(svg('title', {}, `${data.beings[l.a].name} — ${data.beings[l.b].name}: ${tie.type}, last attested ${l.sign < 0 ? 'negative' : 'positive'} at “${data.occasions[l.since].label}” (${dateText(data.occasions[l.since].date)}); ${plural(l.lenses_since, 'meeting')} since${post === undefined ? `; expected ${number(l.expected, 2)}` : `; reconstructed ${number(post, 2)}`}`));
      path.addEventListener('click', () => this.pick({ kind: 'tie', id: tie.id }));
      g.append(path);
    });
    // Beings.
    for (const b of order) {
      const p = pos.get(b), on = live.has(b), being = data.beings[b];
      // Names run outward along their radius, so that neighbours never collide.
      const right = Math.cos(p.a) >= 0;
      const lx = cx + (R + 12) * Math.cos(p.a), ly = cy + (R + 12) * Math.sin(p.a);
      const deg = (p.a * 180) / Math.PI + (right ? 0 : 180);
      const node = svg('g', { class: on ? 'node' : 'node off' },
        svg('circle', { cx: p.x, cy: p.y, r: on ? 5.5 : 2.6, fill: on ? colourOf(being.colour).hex : 'var(--ink-faint)', stroke: 'var(--paper)', 'stroke-width': 2 }),
        svg('text', { x: lx, y: ly, 'text-anchor': right ? 'start' : 'end', 'dominant-baseline': 'middle', transform: `rotate(${deg.toFixed(1)} ${lx.toFixed(1)} ${ly.toFixed(1)})` }, being.short));
      node.append(svg('title', {}, being.name));
      node.addEventListener('click', () => this.pick({ kind: 'being', id: being.id }));
      g.append(node);
    }
    const frustrated = m.triangles.filter((t) => t.product < 0);
    this.side.replaceChildren(...[
      el('p', { class: 'eyebrow', text: 'At the end of' }),
      el('h3', { class: 'moment-date', text: dateText(m.date || date) }),
      m.substituted ? el('p', { class: 'note', text: `The fixtures hold the moments of ${listText(['1 October 1808', '15 November 1814'])} only: this is ${dateText(m.substituted)}.` }) : null,
      el('dl', { class: 'figures' },
        el('dt', { text: 'Beings with live ties' }), el('dd', { text: String(m.nodes.length) }),
        el('dt', { text: 'Live ties' }), el('dd', { text: String(m.links.length) }),
        el('dt', { text: 'Triangles' }), el('dd', { text: String(m.triangles.length) }),
        el('dt', { text: 'Frustrated, as last attested' }), el('dd', { text: `${m.frustrated} of ${m.triangles.length}` }),
        el('dt', { text: 'Frustrated, expected' }), el('dd', { text: number(m.expected, 2) }),
        m.posterior ? el('dt', { text: 'Frustrated, reconstructed' }) : null,
        m.posterior ? el('dd', { text: `${number(m.posterior.frustrated.mean, 2)} (${m.posterior.frustrated.lo}–${m.posterior.frustrated.hi})` }) : null),
      m.posterior ? el('p', { class: 'note', text: `Reconstructed: the signs the record does not attest at their last meeting are inferred from all the others at once, by the persistence of their ties and the balance of the triangles they close (the theory note, §6.4); ${m.posterior.exact ? 'computed exactly, the range holding the middle nine tenths of the law' : 'sampled by the heat bath, the range holding nine draws in ten'}.` }) : null,
      frustrated.length ? el('h4', { text: 'Frustrated triangles' })
        : m.triangles.length ? el('p', { class: 'note', text: `None of the ${plural(m.triangles.length, 'triangle')} is frustrated: here the network is balanced, its signs a flat ℤ₂ connection.` })
          : el('p', { class: 'note', text: 'There is no triangle of live ties at this date, so nothing can be frustrated: the record is silent, not balanced. (A tie lives only from its first attestation to its last.)' }),
      frustrated.length ? el('ul', { class: 'people' }, ...frustrated.map((t) => el('li', {},
        el('span', { class: 'tri frustrated', text: '▲' }),
        el('span', { text: listText(t.beings.map((b) => data.beings[b].short)) }),
        el('span', { class: 'role', text: m.posterior ? `reconstructed ${number(m.posterior.triangles[m.triangles.indexOf(t)], 2)}` : `expected ${number(t.expected, 2)}` })))) : null,
      m.networks?.length ? el('h4', { text: 'The spin networks of this cut' }) : null,
      m.networks?.length ? el('p', { class: 'note', text: 'In the theory\'s own ensemble, the foam crosses this cut along an even subgraph of the open ties: a ℤ₂ spin network. Each is shown with the weight of the foams that cross along it, computed exactly. The empty network is the weight of no foam crossing here (the theory note, Proposition 6.3).' }) : null,
      m.networks?.length ? el('ol', { class: 'networks' }, ...m.networks.map((nw) => el('li', {},
        el('span', { class: 'num', text: number(nw.weight, 3) }), ' ',
        nw.ties.length ? listText(nw.ties.map((t) => `${data.beings[data.ties[t].a].short}–${data.beings[data.ties[t].b].short}`)) : el('span', { class: 'faint', text: 'none: no foam crosses' })))) : null,
      el('div', { class: 'legend-block' },
        el('div', {}, svg('svg', { width: 34, height: 10 }, svg('line', { x1: 2, x2: 32, y1: 5, y2: 5, stroke: INK, 'stroke-width': 1.6 })), ' positive, last attested'),
        el('div', {}, svg('svg', { width: 34, height: 10 }, svg('line', { x1: 2, x2: 32, y1: 5, y2: 5, stroke: INK, 'stroke-width': 1.6, 'stroke-dasharray': '5 3' })), ' negative'),
        el('div', {}, svg('svg', { width: 34, height: 12 }, svg('polygon', { points: '4,11 17,1 30,11', class: 'tri-fill' })), ' a frustrated triangle')),
    ].filter(Boolean));
    this.markSeries(yOfDate(date));
  }

  // The frustration year by year: all triangles (a faint area), the frustrated as last attested
  // (a line) and as expected (a dashed line); one scale, triangles.
  drawSeries() {
    const data = this.data, ms = data.readings.moments;
    const host = clear(this.series);
    if (!ms.length) return;
    const W = 900, H = 190, L = 44, R = 150, T = 18, B = 26;
    const x0 = ms[0].year, x1 = ms[ms.length - 1].year + 1;
    const ymax = Math.max(1, ...ms.map((m) => m.triangles));
    const X = (y) => L + ((y - x0) / (x1 - x0)) * (W - L - R), Y = (v) => T + (1 - v / ymax) * (H - T - B);
    const at = (m) => X(m.year + 1);
    const s = svg('svg', { viewBox: `0 0 ${W} ${H}`, class: 'series-svg', role: 'img', 'aria-label': 'Frustration of the network at the end of each year' });
    for (let v = 0; v <= ymax; v += Math.max(1, Math.ceil(ymax / 4))) {
      s.append(svg('line', { x1: L, x2: W - R, y1: Y(v), y2: Y(v), class: 'gridline' }), svg('text', { x: L - 6, y: Y(v) + 4, 'text-anchor': 'end', class: 'tick' }, String(v)));
    }
    const span = x1 - x0, step = span > 1000 ? 200 : span > 400 ? 100 : span > 150 ? 25 : span > 60 ? 20 : span > 25 ? 5 : 2;
    for (const y of roundYears(x0, x1, step)) s.append(svg('text', { x: X(y), y: H - 8, 'text-anchor': 'middle', class: 'tick' }, yearName(y)));
    const line = (key) => ms.map((m, i) => `${i ? 'L' : 'M'}${at(m)} ${Y(m[key])}`).join('');
    s.append(svg('path', { d: `${line('triangles')}L${at(ms[ms.length - 1])} ${Y(0)}L${at(ms[0])} ${Y(0)}Z`, class: 'area' }));
    // The reconstruction: its mean and the band of nine draws in ten.
    const fm = (data.foam?.moments || []).filter((q) => q.year >= x0 && q.year < x1);
    if (fm.length) {
      const up = fm.map((q, i) => `${i ? 'L' : 'M'}${X(q.year + 1)} ${Y(q.hi)}`).join('');
      const down = [...fm].reverse().map((q) => `L${X(q.year + 1)} ${Y(q.lo)}`).join('');
      s.append(svg('path', { d: `${up}${down}Z`, class: 'band posterior' }),
        svg('path', { d: fm.map((q, i) => `${i ? 'L' : 'M'}${X(q.year + 1)} ${Y(q.mean)}`).join(''), class: 'line posterior' }));
    }
    s.append(svg('path', { d: line('expected'), class: 'line expected' }),
      svg('path', { d: line('frustrated'), class: 'line frustrated' }));
    const last = ms[ms.length - 1];
    s.append(svg('text', { x: at(last) + 8, y: Y(last.triangles) + 4, class: 'direct' }, 'triangles'),
      svg('text', { x: at(last) + 8, y: Y(last.frustrated) + 4, class: 'direct' }, 'frustrated'),
      svg('text', { x: at(last) + 8, y: Y(last.expected) + 16, class: 'direct faint' }, 'expected'));
    this.marker = svg('line', { y1: T - 6, y2: H - B, class: 'series-marker' });
    const hoverLine = svg('line', { y1: T, y2: H - B, class: 'series-hover', visibility: 'hidden' });
    const tip = el('div', { class: 'tip', hidden: true });
    const hit = svg('rect', { x: L, y: T, width: W - L - R, height: H - T - B, fill: 'transparent' });
    const nearest = (e) => {
      const r = s.getBoundingClientRect(), px = ((e.clientX - r.left) / r.width) * W;
      return ms.reduce((a, m) => (Math.abs(at(m) - px) < Math.abs(at(a) - px) ? m : a), ms[0]);
    };
    hit.addEventListener('pointermove', (e) => {
      const m = nearest(e), r = s.getBoundingClientRect();
      hoverLine.setAttribute('x1', at(m)); hoverLine.setAttribute('x2', at(m)); hoverLine.setAttribute('visibility', 'visible');
      tip.hidden = false;
      tip.style.left = `${(at(m) / W) * r.width + 12}px`;
      tip.replaceChildren(el('strong', { text: `End of ${m.year}` }), el('br'),
        `${plural(m.ties, 'live tie')}, ${plural(m.triangles, 'triangle')}`, el('br'),
        `frustrated ${m.frustrated} of ${m.triangles} (expected ${number(m.expected, 2)})`);
    });
    hit.addEventListener('pointerleave', () => { hoverLine.setAttribute('visibility', 'hidden'); tip.hidden = true; });
    hit.addEventListener('click', (e) => this.setDate(`${nearest(e).year}-12-31`));
    s.append(hoverLine, this.marker, hit);
    host.append(el('h4', { class: 'series-title', text: 'Frustration of the network, at the end of each year' }),
      el('div', { class: 'series-figure' }, s, tip),
      el('p', { class: 'note', text: 'Shaded grey: all triangles of live ties. Vermilion: those frustrated with every tie at its last attested sign. Mauve, with its band: the frustration of the reconstruction, the unattested signs inferred by persistence and balance together. Dashed: the older estimate, a carried sign trusted as far as its discipline’s memory. Click a year to see its moment.' }));
    this.seriesScale = { X, W, x0, x1 };
  }

  markSeries(y) {
    if (!this.marker || y === null || y === undefined) return;
    const x = this.seriesScale.X(y);
    this.marker.setAttribute('x1', x);
    this.marker.setAttribute('x2', x);
  }
}

function periodOf(data) {
  const ys = data.occasions.map((o) => o.year);
  return [Math.floor(Math.min(...ys)), Math.ceil(Math.max(...ys))];
}

