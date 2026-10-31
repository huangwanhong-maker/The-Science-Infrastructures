// Currents: the walkers each way past Melbourne's sensors that count each way, at the clock's hour. Each
// sensor's block of street carries two lanes of walkers, those walking forward on the left as left-hand
// traffic keeps (the walkers' handedness the engine tests). A lane holds as many walkers as walk there on
// average, N = λ L / v (Little's identity: λ the count toward that way an hour, L the lane's 120 m, v Weidmann's
// 1.34 m/s), moving at walking speed in real time, each walker drawn half again its size. A walker's colour
// says whether it walks toward the centre of the city's counts or away from it. An arrow on the street is the
// net current, its length growing with the flow. On a date the record holds, the counts are that day's;
// otherwise they are the mean at that hour of the week (the engine's current meters, the domain's
// flow/currents.hpp).

import { sensorBlocks } from './blocks.js';

const SPEED = 1.34;                       // m/s, a walker free (Weidmann 1993)
const INWARD = [0.8, 0.5, 0.36], OUTWARD = [0.25, 0.45, 0.63];   // coral toward the centre, sea blue away: apart in every vision
const hhmm = (x) => `${String(Math.floor(x)).padStart(2, '0')}:${String(Math.round((x % 1) * 60) % 60).padStart(2, '0')}`;
const css = (rgb) => `rgb(${rgb.map((c) => Math.round(255 * c)).join(',')})`;

export async function layer(ctx) {
  const { gis, ask, dom, G } = ctx;
  const { h, fmt } = dom;
  const swarmId = ctx.layerId(':walkers'), arrowId = ctx.layerId(':arrows');
  let clock = ctx.clock, features = [], sea = null, seaDate = null, centre = [0, 0], hand = null, stop = null;
  const status = h('div', { class: 'faint' }, 'reading every sensor\'s counts each way…');
  const finding = h('div', { class: 'faint', style: 'margin-top:3px' });
  const legend = h('div', { class: 'card' }, h('b', {}, 'Currents: walkers each way'),
    ctx.chart.legend([{ name: 'toward the centre of the counts', color: css(INWARD), dot: true }, { name: 'away from it', color: css(OUTWARD), dot: true }]),
    h('div', { class: 'faint' }, 'Each walker is one person, as many on each lane as walk there on average at this hour (N = λL/v), at walking speed in real time; the arrow is the net current.'),
    finding, status);

  // The rate toward each way (people an hour) at the clock: that day's count if the record holds it, else the
  // mean at that hour of the week.
  const rate = (f, w) => {
    const d = sea && sea.get(f.id);
    const day = d && (w ? d.way2 : d.way1);
    const x = Math.max(0, clock.hour - 0.5), i = Math.min(23, Math.floor(x)), u = x - i;
    if (day && Number.isFinite(day[i]) && Number.isFinite(day[Math.min(23, i + 1)])) return day[i] + u * (day[Math.min(23, i + 1)] - day[i]);
    const week = w ? f.way2 : f.way1, k = clock.weekday * 24 + clock.hour - 0.5, k0 = ((Math.floor(k) % 168) + 168) % 168, v = k - Math.floor(k);
    return week[k0] + v * (week[(k0 + 1) % 168] - week[k0]);
  };
  const net = (f) => rate(f, 0) - rate(f, 1);
  const inward = (f, w) => { const d = f.toward[w]; return d[0] * (centre[0] - f.foot[0]) + d[1] * (centre[1] - f.foot[1]) > 0; };

  // The lanes: forward along the block's line on its left edge, backward on its right.
  function lanes(f) {
    const e = G.edges(f.line, f.width / 2);
    const forward = (w) => f.toward[w][0] * f.axis[0] + f.toward[w][1] * f.axis[1] > 0;
    return [0, 1].map((w) => ({ line: forward(w) ? e.left : e.right, forward: forward(w), length: G.lineLength(f.line), phase: Math.random() * 50, w }));
  }

  const capacity = 9000;
  const xyz = new Float32Array(3 * capacity), rgb = new Float32Array(3 * capacity);
  const jitter = (j) => (((j * 2654435761) >>> 0) / 4294967296 - 0.5) * 1.1;
  let last = 0;
  function frame(now) {
    const dt = last ? Math.min(0.1, (now - last) / 1000) : 0;
    last = now;
    let n = 0;
    for (const f of features) {
      for (const lane of f.lanes) {
        const lam = Math.max(0, rate(f, lane.w)) / 3600;
        const N = Math.floor((lam * lane.length) / SPEED + 0.5);
        if (!N) continue;
        lane.phase = (lane.phase + SPEED * dt) % lane.length;
        const gap = lane.length / N, colour = inward(f, lane.w) ? INWARD : OUTWARD;
        for (let j = 0; j < N && n < capacity; j++) {
          let s = (lane.phase + j * gap) % lane.length;
          if (!lane.forward) s = lane.length - s;
          const p = G.pointAlong(lane.line, s), off = jitter(j + 7 * f.id + lane.w);
          xyz[3 * n] = p.x - p.dy * off, xyz[3 * n + 1] = p.y + p.dx * off, xyz[3 * n + 2] = 0.12;
          rgb[3 * n] = colour[0], rgb[3 * n + 1] = colour[1], rgb[3 * n + 2] = colour[2];
          n++;
        }
      }
    }
    gis.setSwarm(swarmId, n, xyz, rgb);
  }

  function arrows() {
    gis.restyle(arrowId, {});
  }

  async function loadDay(date) {
    if (seaDate === date) return;
    seaDate = date;
    try {
      const d = await ask('sea', { date });
      if (seaDate !== date) return;
      sea = new Map(d.sensors.map((x) => [x.id, x]));
      const any = d.sensors.some((x) => x.way1.some(Number.isFinite));
      status.textContent = any ? `${date}: that day's counts each way` : `${date} lies outside the record: the mean at each hour of the week (October 2024 to October 2026)`;
      arrows();
    } catch (e) {
      status.textContent = String(e.message || e);
    }
  }

  (async () => {
    try {
      const c = await ask('currents');
      hand = c.handedness;
      features = sensorBlocks(ctx, c.sensors).filter((f) => f.toward.length === 2 && f.toward.every(Boolean));
      let wsum = 0;
      for (const f of features) { centre[0] += f.transport * f.foot[0], centre[1] += f.transport * f.foot[1], wsum += f.transport; }
      centre = [centre[0] / wsum, centre[1] / wsum];
      for (const f of features) f.lanes = lanes(f);
      gis.addSwarm(swarmId, { capacity, radius: 0.6, height: 2.5 });
      gis.addArrows(arrowId, features, {
        position: (f) => [f.foot[0], f.foot[1]],
        direction: (f) => (net(f) >= 0 ? f.toward[0] : f.toward[1]),
        length: (f) => (Math.abs(net(f)) < 15 ? 0.01 : 6 + 0.03 * Math.abs(net(f))),
        width: (f) => Math.min(f.width * 0.9, 5 + 0.004 * Math.abs(net(f))),
        color: (f) => (inward(f, net(f) >= 0 ? 0 : 1) ? INWARD : OUTWARD),
        z: 0.5, opacity: 0.92,
        pick: (f) => `${f.name}\nnet ${Math.round(Math.abs(net(f)))} an hour toward ${f.directions[net(f) >= 0 ? 0 : 1] || ''}; ${Math.round(rate(f, 0))} toward ${f.directions[0]}, ${Math.round(rate(f, 1))} toward ${f.directions[1]}`,
      });
      if (hand && hand.footpaths) finding.textContent = `Walkers keep left: on ${hand.left} of ${hand.footpaths} footpaths the residual current runs as left-hand traffic would (sign test p = ${fmt(hand.pLeft, 2)}; footpaths of one street share their walkers).`;
      status.textContent = '';
      stop = gis.onFrame(frame);
      await loadDay(clock.date);
    } catch (e) {
      status.textContent = String(e.message || e);
    }
  })();

  return {
    legend,
    setTime(c) {
      clock = c;
      if (!features.length) return;
      if (c.date !== seaDate) loadDay(c.date);
      arrows();
    },
    describe(hit) {
      const f = hit.feature;
      const days = [0, 1].map((w) => Array.from({ length: 24 }, (_, hr) => { let s = 0; for (let d = 0; d < 5; d++) s += (w ? f.way2 : f.way1)[24 * d + hr]; return s / 5; }));
      const xs = Array.from({ length: 24 }, (_, k) => k + 0.5);
      const chart = ctx.chart.lines([
        { name: `toward ${f.directions[0]}`, x: xs, y: days[0], color: css(inward(f, 0) ? INWARD : OUTWARD), width: 2 },
        { name: `toward ${f.directions[1]}`, x: xs, y: days[1], color: css(inward(f, 1) ? INWARD : OUTWARD), width: 2, dash: '5 3' },
      ], { height: 150, width: 360, xlim: [0, 24], xlabel: 'hour of the mean weekday', ylabel: 'people an hour', marks: [{ x: clock.hour, label: hhmm(clock.hour) }] });
      const nowNet = net(f), irr = nowNet !== 0 && rate(f, 0) >= 1 && rate(f, 1) >= 1 ? nowNet * Math.log(rate(f, 0) / rate(f, 1)) : NaN;
      const sel = {
        kind: 'sensor', id: f.id, name: f.name, title: f.name,
        lines: [
          `now ${Math.round(rate(f, 0)).toLocaleString('en')} an hour toward ${f.directions[0]}, ${Math.round(rate(f, 1)).toLocaleString('en')} toward ${f.directions[1]}: a net ${Math.round(Math.abs(nowNet)).toLocaleString('en')} toward ${f.directions[nowNet >= 0 ? 0 : 1]}`,
          `on weekdays the flood runs toward ${f.directions[0]} at ${hhmm(f.floodAt)} (${Math.round(f.flood)} net an hour), the ebb toward ${f.directions[1]} at ${hhmm(f.ebbAt)} (${Math.round(-f.ebb)})`,
          `residual current ${Math.round(Math.abs(f.residual))} an hour toward ${f.directions[f.residual >= 0 ? 0 : 1]} (${fmt(Math.abs(100 * f.residual / f.transport), 2)}% of the transport); irreversibility now ${Number.isFinite(irr) ? fmt(irr, 3) : '—'} k_B an hour, ${fmt(f.sigma, 3)} on average`,
        ],
        extra: dom.h('div', {}, chart),
      };
      sel.actions = [
        { title: 'its rhythm', run: () => ctx.openTool('rhythm', sel) },
        { title: 'tides and currents of the city', run: () => ctx.openTool('sea', sel) },
      ];
      return sel;
    },
    dispose() { if (stop) stop(); },
  };
}
