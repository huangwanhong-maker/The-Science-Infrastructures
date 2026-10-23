// The foam: a history as a 2-complex in spacetime, drawn with three.js (the theory note, §2, §5–6).
//
// Time rises. Each being is a thread through the occasions it attends, in a place of its own on
// the ground plane: beings that meet and are tied often stand near one another, and the disciplines
// hold sectors of the ring. Between two of its occasions a thread drifts back toward its place, so
// that a life is seen to leave and return. An occasion is a knot where the threads of those present
// meet. Each tie's lenses, its 2-cells, are ribbons between its two threads from one meeting to the
// next, pinched at the meetings: coloured by discipline and darkening where the tie is negative, or
// by its spin in the foam, its flip, or the certainty of an inferred sign.
//
// Opened (the refined complex K̂), each occasion becomes the presences around it: the rungs of the
// ties there are short bars, the triangles of ties are small faces, frustrated ones in vermilion,
// and the lenses become squares between rungs. A cut at a date is the network of the moment, the
// boundary of the history at that instant: where the threads cross it, the live ties between them,
// and its frustrated triangles.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';
import { colourOf, FRUSTRATION, ACCENT } from './palette.js';
import { placeBeings, heights, HEIGHT, RING } from './foamlayout.js';
import { yOfDate, roundYears, yearName } from './format.js';

const NIGHT = '#241a26';
const IVORY = new THREE.Color('#f3e6d3');
const SMOKE = new THREE.Color('#1b141d');
const GOLD = new THREE.Color('#f0c674');
const MAUVE = new THREE.Color(ACCENT).lerp(new THREE.Color('#ffffff'), 0.25);
const VERMILION = new THREE.Color(FRUSTRATION).lerp(new THREE.Color('#ff7a5c'), 0.35);

// The height of a decimal year, and the year of a height, through the occasions' heights.
function yearToHeight(data, ys, year) {
  const order = data.complex.vertices;
  if (year <= data.occasions[order[0]].year) return ys[order[0]];
  for (let k = 1; k < order.length; k++) {
    const a = data.occasions[order[k - 1]], b = data.occasions[order[k]];
    if (year <= b.year) return b.year === a.year ? ys[order[k]] : ys[order[k - 1]] + ((year - a.year) / (b.year - a.year)) * (ys[order[k]] - ys[order[k - 1]]);
  }
  return ys[order[order.length - 1]];
}

export class FoamView {
  constructor(host, { onPick, onHover } = {}) {
    this.host = host;
    this.onPick = onPick;
    this.onHover = onHover;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    this.renderer.setClearColor(new THREE.Color(NIGHT), 1);
    host.append(this.renderer.domElement);
    this.labels = new CSS2DRenderer();
    this.labels.domElement.className = 'foam-labels';
    host.append(this.labels.domElement);
    this.scene = new THREE.Scene();
    this.scene.fog = new THREE.Fog(NIGHT, 60, 140);
    this.camera = new THREE.PerspectiveCamera(42, 1, 0.1, 400);
    this.camera.position.set(44, HEIGHT * 0.85, 50);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.target.set(0, HEIGHT * 0.48, 0);
    this.controls.enableDamping = true;
    this.controls.addEventListener('change', () => this.draw());
    this.scene.add(new THREE.AmbientLight(0xffffff, 0.9));
    const sun = new THREE.DirectionalLight(0xfff1e0, 1.1);
    sun.position.set(20, 60, 30);
    this.scene.add(sun);
    this.root = new THREE.Group();
    this.scene.add(this.root);
    this.raycaster = new THREE.Raycaster();
    this.pointer = new THREE.Vector2();
    this.pickables = [];
    this.key = null;
    const el = this.renderer.domElement;
    el.addEventListener('pointermove', (e) => this.hover(e));
    el.addEventListener('pointerleave', () => this.onHover?.(null));
    let down = null;
    el.addEventListener('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
    el.addEventListener('pointerup', (e) => { if (down && Math.hypot(e.clientX - down[0], e.clientY - down[1]) < 5) this.click(e); down = null; });
    this.resizer = new ResizeObserver(() => this.resize());
    this.resizer.observe(host);
    const tick = () => { this.frame = requestAnimationFrame(tick); if (this.controls.update()) this.draw(); this.step?.(); };
    tick();
    this.resize();
  }

  resize() {
    const w = this.host.clientWidth || 800, h = this.host.clientHeight || 600;
    this.renderer.setSize(w, h);
    this.labels.setSize(w, h);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.draw();
  }

  draw() {
    this.renderer.render(this.scene, this.camera);
    this.labels.render(this.scene, this.camera);
    this.declutter();
  }

  // Names that would overlap on the screen: the more present being keeps its name.
  declutter() {
    if (!this.threads) return;
    const named = this.threads.map((m, i) => (m && m.visible ? { el: m.userData.label, weight: (this.data?.beings[i].occasions.length || 0) + (m.userData.lit ? 1000 : 0) } : null)).filter(Boolean);
    named.sort((a, b) => b.weight - a.weight);
    const placed = [];
    for (const n of named) {
      n.el.style.visibility = 'visible';
      const r = n.el.getBoundingClientRect();
      if (!r.width) continue;
      const clash = placed.some((q) => r.left < q.right + 2 && r.right > q.left - 2 && r.top < q.bottom && r.bottom > q.top);
      if (clash) n.el.style.visibility = 'hidden';
      else placed.push(r);
    }
  }

  // The view of a canon under the page's state: { mode, opened, time, date, focus, loop, hide }.
  set(data, view) {
    const key = `${data.canon.id}|${view.time}|${view.opened}`;
    if (key !== this.key) { this.key = key; this.data = data; this.layout(data, view); this.build(data, view); }
    this.data = data;
    this.view = view;
    this.style(data, view);
    this.cut(data, view);
    this.aim(data, view);
    this.draw();
  }

  // When what is followed changes, the camera glides to it: an occasion, the middle of a tie's
  // lenses, a being's life, or the height of the cut.
  aim(data, view) {
    const f = view.focus;
    const key = `${f ? `${f.kind}:${f.id}` : ''}|${view.loop ?? ''}|${view.date || ''}|${this.key}`;
    if (key === this.aimed) return;
    const first = this.aimed === undefined;
    this.aimed = key;
    let target = null, distance = 30;
    if (f?.kind === 'occasion') {
      const o = data.occasions.findIndex((x) => x.id === f.id);
      if (o >= 0) { target = this.centre[o].clone().setY(this.y[o]); distance = 18; }
    } else if (f?.kind === 'tie') {
      const t = data.ties.findIndex((x) => x.id === f.id);
      const faces = data.complex.faces.filter((x) => x.tie === t);
      if (faces.length) {
        const ys = faces.flatMap((x) => [this.y[x.from], this.y[x.to]]);
        const mid = (Math.min(...ys) + Math.max(...ys)) / 2;
        const p = this.at(data.ties[t].a, mid) || this.centre[faces[0].from].clone().setY(mid);
        target = p; distance = 14 + (Math.max(...ys) - Math.min(...ys)) * 0.6;
      }
    } else if (f?.kind === 'being') {
      const b = data.beings.findIndex((x) => x.id === f.id);
      const segs = this.segments[b] || [];
      if (segs.length) { const mid = (segs[0].curve.v0.y + segs[segs.length - 1].curve.v2.y) / 2; target = this.at(b, mid); distance = 34; }
    } else if (view.date) {
      const m = String(view.date).match(/^(-?\d{1,4})/);
      if (m) { target = new THREE.Vector3(0, this.heightOfYear(Number(m[1]) + 0.5), 0); distance = 40; }
    }
    if (!target) { if (!first) this.glide(new THREE.Vector3(0, HEIGHT * 0.48, 0), 66); return; }
    if (first) return;
    this.glide(target, distance);
  }

  heightOfYear(year) { return yearToHeight(this.data, this.y, year); }

  glide(target, distance) {
    const from = this.controls.target.clone(), fromCam = this.camera.position.clone();
    const dir = fromCam.clone().sub(from).normalize();
    const toCam = target.clone().add(dir.multiplyScalar(distance));
    const t0 = performance.now(), T = 900;
    this.step = () => {
      const u = Math.min(1, (performance.now() - t0) / T), e = u * u * (3 - 2 * u);
      this.controls.target.lerpVectors(from, target, e);
      this.camera.position.lerpVectors(fromCam, toCam, e);
      this.controls.update();
      this.draw();
      if (u >= 1) this.step = null;
    };
  }

  // ---- the layout ---------------------------------------------------------------------------------
  layout(data, view) {
    this.home = placeBeings(data);
    this.y = heights(data, view.time);
    // Each presence: at its occasion's centre, or on a small ring round it when the occasion is opened.
    this.centre = data.occasions.map((o) => {
      const here = o.present.filter((p) => !p.derived);
      const list = here.length ? here : o.present;
      const x = list.reduce((s, p) => s + this.home[p.being].x, 0) / list.length;
      const z = list.reduce((s, p) => s + this.home[p.being].z, 0) / list.length;
      return new THREE.Vector3(x, 0, z);
    });
    this.presence = new Map();  // "being,occasion" -> point
    data.occasions.forEach((o, i) => {
      const c = this.centre[i].clone().setY(this.y[i]);
      const n = o.present.length;
      const r = view.opened && n > 1 ? 0.35 + 0.16 * n : 0;
      const ranked = o.present.map((p) => ({ b: p.being, a: Math.atan2(this.home[p.being].z - c.z, this.home[p.being].x - c.x) })).sort((p, q) => p.a - q.a);
      ranked.forEach(({ b }, k) => {
        const a = (2 * Math.PI * k) / n + 0.3;
        this.presence.set(`${b},${i}`, new THREE.Vector3(c.x + r * Math.cos(a), c.y, c.z + r * Math.sin(a)));
      });
    });
    // Each thread: quadratic Bézier segments between consecutive presences, the control point at the
    // middle height and drawn toward the being's place in proportion to the time apart. The height is
    // then linear in the parameter, so that a thread can be read at any height.
    this.segments = data.beings.map((b, i) => {
      const pts = b.occasions.map((o) => this.presence.get(`${i},${o}`));
      const out = [];
      for (let k = 0; k + 1 < pts.length; k++) {
        const p0 = pts[k], p1 = pts[k + 1];
        const pull = 0.7 * (1 - Math.exp(-Math.abs(p1.y - p0.y) / 6));
        const mid = p0.clone().add(p1).multiplyScalar(0.5);
        const c = new THREE.Vector3(mid.x + pull * (this.home[i].x - mid.x), mid.y, mid.z + pull * (this.home[i].z - mid.z));
        out.push({ from: b.occasions[k], to: b.occasions[k + 1], curve: new THREE.QuadraticBezierCurve3(p0, c, p1) });
      }
      return out;
    });
  }

  // A being's thread at a height (null outside its life).
  at(b, y) {
    for (const s of this.segments[b]) {
      const y0 = s.curve.v0.y, y1 = s.curve.v2.y;
      if (y >= Math.min(y0, y1) - 1e-9 && y <= Math.max(y0, y1) + 1e-9) return s.curve.getPoint(y1 === y0 ? 0 : (y - y0) / (y1 - y0));
    }
    return null;
  }

  // ---- the meshes ----------------------------------------------------------------------------------
  build(data, view) {
    for (const child of [...this.root.children]) { this.root.remove(child); child.traverse?.((o) => { o.geometry?.dispose(); if (o.element) o.element.remove(); }); }
    this.pickables = [];
    this.threads = [];
    this.lenses = [];
    this.knots = [];
    this.rungs = [];
    this.triangles = [];
    // Threads.
    data.beings.forEach((b, i) => {
      const segs = this.segments[i];
      if (!segs.length) { this.threads.push(null); return; }
      const path = new THREE.CurvePath();
      segs.forEach((s) => path.add(s.curve));
      const colour = new THREE.Color(colourOf(b.colour).hex).lerp(IVORY, 0.25);
      const mesh = new THREE.Mesh(new THREE.TubeGeometry(path, Math.max(12, segs.length * 14), 0.075, 6, false),
        new THREE.MeshLambertMaterial({ color: colour, transparent: true, opacity: 0.95, emissive: colour, emissiveIntensity: 0.25 }));
      mesh.userData = { kind: 'being', being: i };
      this.root.add(mesh);
      this.pickables.push(mesh);
      this.threads.push(mesh);
      // The name, where the life begins.
      const div = document.createElement('div');
      div.className = 'foam-name';
      div.textContent = b.short || b.name;
      const label = new CSS2DObject(div);
      const rest = segs.reduce((a, s) => (Math.abs(s.curve.v2.y - s.curve.v0.y) > Math.abs(a.curve.v2.y - a.curve.v0.y) ? s : a), segs[0]);
      const at = rest.curve.getPoint(0.5);
      const out = new THREE.Vector3(at.x, 0, at.z).normalize().multiplyScalar(0.9);
      label.position.copy(at).add(out);
      mesh.add(label);
      mesh.userData.label = div;
    });
    // Lenses: films between the two threads, sampled at the same heights and across, strong at
    // their boundary threads and faint in the middle, like soap films in a foam.
    const N = 18, M = 6;
    data.complex.faces.forEach((f, fi) => {
      const t = data.ties[f.tie];
      const y0 = this.y[f.from], y1 = this.y[f.to];
      const pos = [], col = [];
      const base = new THREE.Color(colourOf(t.colour).hex);
      const tone = (s) => (s > 0 ? base.clone().lerp(IVORY, 0.15) : base.clone().lerp(SMOKE, 0.62));
      const c0 = tone(f.signs[0]), c1 = tone(f.signs[1]);
      for (let k = 0; k <= N; k++) {
        const y = y0 + ((y1 - y0) * k) / N;
        const a = k === 0 ? this.presence.get(`${t.a},${f.from}`) : k === N ? this.presence.get(`${t.a},${f.to}`) : this.at(t.a, y);
        const b = k === 0 ? this.presence.get(`${t.b},${f.from}`) : k === N ? this.presence.get(`${t.b},${f.to}`) : this.at(t.b, y);
        if (!a || !b) return;
        const c = c0.clone().lerp(c1, k / N);
        for (let m = 0; m <= M; m++) {
          const w = m / M;
          pos.push(a.x + (b.x - a.x) * w, a.y + (b.y - a.y) * w, a.z + (b.z - a.z) * w);
          const edge = Math.pow(Math.abs(2 * w - 1), 2.2);
          col.push(c.r, c.g, c.b, 0.22 + 0.78 * edge);
        }
      }
      const index = [];
      for (let k = 0; k < N; k++)
        for (let m = 0; m < M; m++) {
          const i = k * (M + 1) + m, j = i + M + 1;
          index.push(i, i + 1, j, i + 1, j + 1, j);
        }
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
      g.setAttribute('color', new THREE.Float32BufferAttribute(col, 4));
      g.setIndex(index);
      const mesh = new THREE.Mesh(g, new THREE.MeshBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.3, side: THREE.DoubleSide, depthWrite: false }));
      mesh.userData = { kind: 'lens', face: fi, tie: f.tie, base: col.slice() };
      mesh.renderOrder = 1;
      this.root.add(mesh);
      this.pickables.push(mesh);
      this.lenses.push(mesh);
    });
    // Knots: the occasions.
    const knot = new THREE.SphereGeometry(1, 14, 10);
    data.occasions.forEach((o, i) => {
      const live = o.present.filter((p) => !p.derived).length;
      const mesh = new THREE.Mesh(knot, new THREE.MeshLambertMaterial({ color: IVORY, emissive: IVORY, emissiveIntensity: 0.35, transparent: true, opacity: view.opened ? 0.35 : 0.95 }));
      const r = 0.13 + 0.035 * live;
      mesh.scale.set(r, r, r);
      mesh.position.copy(this.centre[i]).setY(this.y[i]);
      mesh.userData = { kind: 'occasion', occasion: i };
      this.root.add(mesh);
      this.pickables.push(mesh);
      this.knots.push(mesh);
    });
    // Opened: rungs and the triangles of ties at each occasion.
    if (view.opened) {
      const foam = data.foam;
      data.ties.forEach((t, ti) => {
        t.rungs.forEach((g, k) => {
          const a = this.presence.get(`${t.a},${g.occasion}`), b = this.presence.get(`${t.b},${g.occasion}`);
          if (!a || !b) return;
          const p = foam?.ties[ti]?.posterior?.[k] ?? g.sign;
          const colour = (g.attested ? g.sign : Math.sign(p) || g.sign) > 0 ? GOLD : new THREE.Color('#a78bfa');
          const geom = new THREE.BufferGeometry().setFromPoints([a, b]);
          const mat = g.attested ? new THREE.LineBasicMaterial({ color: colour, transparent: true, opacity: 0.95 })
            : new THREE.LineDashedMaterial({ color: colour, dashSize: 0.12, gapSize: 0.1, transparent: true, opacity: 0.25 + 0.7 * Math.abs(p) });
          const line = new THREE.Line(geom, mat);
          if (!g.attested) line.computeLineDistances();
          line.userData = { kind: 'rung', tie: ti, rung: k };
          this.root.add(line);
          this.rungs.push(line);
        });
      });
      (foam?.triangles || []).forEach((T, i) => {
        const pts = T.beings.map((b) => this.presence.get(`${b},${T.occasion}`));
        if (pts.some((p) => !p)) return;
        const product = T.ties.reduce((s, t) => s * (data.ties[t].rungs.find((g) => g.occasion === T.occasion)?.sign ?? 1), 1);
        const g = new THREE.BufferGeometry().setFromPoints(pts);
        g.setIndex([0, 1, 2]);
        const mesh = new THREE.Mesh(g, new THREE.MeshBasicMaterial({ color: product < 0 ? VERMILION : IVORY, transparent: true, opacity: product < 0 ? 0.75 : 0.12, side: THREE.DoubleSide, depthWrite: false }));
        mesh.userData = { kind: 'triangle', triangle: i, occasion: T.occasion, product };
        mesh.renderOrder = 2;
        this.root.add(mesh);
        this.pickables.push(mesh);
        this.triangles.push(mesh);
      });
    } else {
      // Closed: a frustrated triangle at an occasion is a vermilion ring round its knot.
      const frustrated = new Set(data.readings.triangles.filter((t) => t.product < 0).map((t) => t.occasion));
      for (const o of frustrated) {
        const ring = new THREE.Mesh(new THREE.TorusGeometry(0.55, 0.05, 6, 32), new THREE.MeshBasicMaterial({ color: VERMILION }));
        ring.position.copy(this.centre[o]).setY(this.y[o]);
        ring.rotation.x = Math.PI / 2;
        ring.userData = { kind: 'occasion', occasion: o };
        this.root.add(ring);
      }
    }
    // The axis of time, outside the ring.
    const order = data.complex.vertices;
    const yr0 = data.occasions[order[0]].year, yr1 = data.occasions[order[order.length - 1]].year;
    const axis = new THREE.Group();
    const ax = RING + 6;
    axis.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(ax, 0, 0), new THREE.Vector3(ax, HEIGHT, 0)]), new THREE.LineBasicMaterial({ color: IVORY, transparent: true, opacity: 0.35 })));
    const span = yr1 - yr0;
    const step = span > 400 ? 100 : span > 150 ? 25 : span > 60 ? 10 : 5;
    let labelled = -Infinity;
    for (const yr of roundYears(yr0, yr1, step)) {
      const h = yearToHeight(data, this.y, yr);
      axis.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(ax - 0.4, h, 0), new THREE.Vector3(ax + 0.4, h, 0)]), new THREE.LineBasicMaterial({ color: IVORY, transparent: true, opacity: 0.5 })));
      // In the time of events an empty century shrinks to nothing: its years are marked, not named.
      if (h - labelled < 1.6) continue;
      labelled = h;
      const div = document.createElement('div');
      div.className = 'foam-year';
      div.textContent = yearName(yr);
      const label = new CSS2DObject(div);
      label.position.set(ax + 1.4, h, 0);
      axis.add(label);
    }
    this.root.add(axis);
    // The ground: the ring of the disciplines.
    const ground = new THREE.Mesh(new THREE.RingGeometry(RING - 0.05, RING + 0.05, 128), new THREE.MeshBasicMaterial({ color: IVORY, transparent: true, opacity: 0.18, side: THREE.DoubleSide }));
    ground.rotation.x = -Math.PI / 2;
    ground.position.y = -0.4;
    this.root.add(ground);
    // The cut, the moment's plane, built empty and filled by cut().
    this.cutGroup = new THREE.Group();
    this.root.add(this.cutGroup);
  }

  // ---- the look under the state ---------------------------------------------------------------------
  style(data, view) {
    const f = view.focus;
    const being = f?.kind === 'being' ? data.beings.findIndex((b) => b.id === f.id) : -1;
    const tie = f?.kind === 'tie' ? data.ties.findIndex((t) => t.id === f.id) : -1;
    const occasion = f?.kind === 'occasion' ? data.occasions.findIndex((o) => o.id === f.id) : -1;
    const loop = view.loop !== null && view.loop !== undefined ? data.invariants.loops[view.loop] : null;
    const loopBeings = new Set(loop ? loop.beings : []);
    const focused = being >= 0 || tie >= 0 || occasion >= 0 || loop;
    const hidden = (colour) => view.hide?.has(colour);
    const at = occasion >= 0 ? new Set(data.occasions[occasion].present.map((p) => p.being)) : null;
    const lit = (b) => (being >= 0 ? b === being : tie >= 0 ? b === data.ties[tie].a || b === data.ties[tie].b : at ? at.has(b) : loop ? loopBeings.has(b) : true);
    this.threads.forEach((m, i) => {
      if (!m) return;
      const on = lit(i) && !hidden(data.beings[i].colour);
      m.visible = !hidden(data.beings[i].colour);
      m.material.opacity = on ? 0.95 : 0.16;
      m.userData.label.style.opacity = on ? 1 : 0.25;
      m.userData.lit = focused && on;
    });
    const spins = new Map();
    (data.foam?.cells || []).forEach((c) => { if (c.kind === 'lens') spins.set(c.index, c); });
    this.lenses.forEach((m) => {
      const fi = m.userData.face, face = data.complex.faces[fi], t = data.ties[face.tie];
      m.visible = !hidden(t.colour);
      const own = tie >= 0 ? face.tie === tie : being >= 0 ? t.a === being || t.b === being : at ? at.has(t.a) && at.has(t.b) : loop ? false : true;
      let opacity = focused ? (own ? 0.8 : 0.05) : 0.34;
      const colours = m.geometry.getAttribute('color');
      const base = m.userData.base;
      if (view.mode === 'spin') {
        const c = spins.get(fi);
        const j = Math.max(0, Math.min(1, c?.spin ?? 0));
        const tint = GOLD.clone().lerp(IVORY, 0.2);
        for (let k = 0; k < colours.count; k++) colours.setXYZW(k, tint.r, tint.g, tint.b, base[4 * k + 3]);
        opacity = (focused && !own ? 0.04 : 0.03) + 0.9 * Math.sqrt(j);
      } else if (view.mode === 'flips') {
        const tint = face.flip ? GOLD : new THREE.Color('#6b5a6f');
        for (let k = 0; k < colours.count; k++) colours.setXYZW(k, tint.r, tint.g, tint.b, base[4 * k + 3]);
        opacity = face.flip ? (focused && !own ? 0.15 : 0.8) : focused && !own ? 0.03 : 0.08;
      } else if (view.mode === 'inference') {
        const post = data.foam?.ties[face.tie]?.posterior || [];
        const k0 = t.rungs.findIndex((g) => g.occasion === face.from), k1 = t.rungs.findIndex((g) => g.occasion === face.to);
        const certainty = (k) => (t.rungs[k]?.attested ? 1 : Math.abs(post[k] ?? 0));
        const latent = !t.rungs[k0]?.attested || !t.rungs[k1]?.attested;
        const tint = latent ? MAUVE : new THREE.Color('#8d7b91');
        for (let k = 0; k < colours.count; k++) colours.setXYZW(k, tint.r, tint.g, tint.b, base[4 * k + 3]);
        opacity = latent ? 0.2 + 0.7 * (1 - Math.min(certainty(k0), certainty(k1))) : 0.05;
      } else {
        for (let k = 0; k < colours.count; k++) colours.setXYZW(k, base[4 * k], base[4 * k + 1], base[4 * k + 2], base[4 * k + 3]);
      }
      colours.needsUpdate = true;
      m.material.opacity = opacity;
    });
    this.knots.forEach((m, i) => {
      const on = occasion >= 0 ? i === occasion : !focused || data.occasions[i].present.some((p) => lit(p.being));
      m.material.emissiveIntensity = i === occasion ? 1 : 0.35;
      m.material.opacity = (view.opened ? 0.35 : 0.95) * (on ? 1 : 0.25);
    });
    // The loop drawn over the threads.
    if (this.loopMesh) { this.root.remove(this.loopMesh); this.loopMesh.geometry.dispose(); this.loopMesh = null; }
    if (loop) {
      const path = new THREE.CurvePath();
      for (const s of loop.steps) {
        const seg = this.segments[s.being].find((x) => (x.from === s.from && x.to === s.to) || (x.from === s.to && x.to === s.from));
        if (!seg) continue;
        const c = seg.curve;
        path.add(seg.from === s.from ? c : new THREE.QuadraticBezierCurve3(c.v2, c.v1, c.v0));
      }
      if (path.curves.length) {
        this.loopMesh = new THREE.Mesh(new THREE.TubeGeometry(path, path.curves.length * 24, 0.17, 8, false), new THREE.MeshBasicMaterial({ color: MAUVE, transparent: true, opacity: 0.9 }));
        this.root.add(this.loopMesh);
      }
    }
  }

  // The cut at the chosen date: the network of the moment where the threads cross it.
  cut(data, view) {
    for (const child of [...this.cutGroup.children]) { this.cutGroup.remove(child); child.traverse?.((o) => { o.geometry?.dispose(); if (o.element) o.element.remove(); }); }
    if (!view.date) return;
    const year = yOfDate(view.date);
    if (year === null) return;
    const h = yearToHeight(data, this.y, year);
    const plane = new THREE.Mesh(new THREE.CircleGeometry(RING + 4, 96), new THREE.MeshBasicMaterial({ color: IVORY, transparent: true, opacity: 0.07, side: THREE.DoubleSide, depthWrite: false }));
    plane.rotation.x = -Math.PI / 2;
    plane.position.y = h;
    this.cutGroup.add(plane);
    const rim = new THREE.Mesh(new THREE.RingGeometry(RING + 3.9, RING + 4, 128), new THREE.MeshBasicMaterial({ color: MAUVE, transparent: true, opacity: 0.7, side: THREE.DoubleSide }));
    rim.rotation.x = -Math.PI / 2;
    rim.position.y = h;
    this.cutGroup.add(rim);
    // The live ties at the date, each at its last sign (the moment of the readings, theory §3.6).
    const jdn = Math.round(2451545 + (year - 2000) * 365.25);
    const live = [];
    data.ties.forEach((t, ti) => {
      if (!t.life || !t.rungs.length) return;
      if (data.occasions[t.life[0]].first > jdn || data.occasions[t.life[1]].last < jdn) return;
      let last = null;
      t.rungs.forEach((g, k) => { if (data.occasions[g.occasion].first <= jdn) last = k; });
      if (last === null) return;
      const a = this.at(t.a, h), b = this.at(t.b, h);
      if (!a || !b) return;
      const post = data.foam?.ties[ti]?.posterior?.[last];
      live.push({ t: ti, a: t.a, b: t.b, sign: t.rungs[last].sign, certainty: t.rungs[last].attested ? 1 : Math.abs(post ?? t.rungs[last].expected ?? 0), pa: a, pb: b });
    });
    const dots = new Set();
    for (const l of live) {
      if (view.hide?.has(data.ties[l.t].colour)) continue;
      const geom = new THREE.BufferGeometry().setFromPoints([l.pa, l.pb]);
      const line = new THREE.Line(geom, l.sign > 0
        ? new THREE.LineBasicMaterial({ color: GOLD, transparent: true, opacity: 0.35 + 0.6 * l.certainty })
        : new THREE.LineDashedMaterial({ color: '#b9a2e0', dashSize: 0.3, gapSize: 0.2, transparent: true, opacity: 0.35 + 0.6 * l.certainty }));
      if (l.sign < 0) line.computeLineDistances();
      this.cutGroup.add(line);
      dots.add(l.a); dots.add(l.b);
    }
    // Its frustrated triangles.
    const by = new Map(live.map((l) => [l.a < l.b ? `${l.a},${l.b}` : `${l.b},${l.a}`, l]));
    const nodes = [...dots];
    for (let i = 0; i < nodes.length; i++)
      for (let j = i + 1; j < nodes.length; j++)
        for (let k = j + 1; k < nodes.length; k++) {
          const key = (x, y) => (x < y ? `${x},${y}` : `${y},${x}`);
          const e = [by.get(key(nodes[i], nodes[j])), by.get(key(nodes[j], nodes[k])), by.get(key(nodes[i], nodes[k]))];
          if (e.some((x) => !x) || e[0].sign * e[1].sign * e[2].sign > 0) continue;
          const pts = [this.at(nodes[i], h), this.at(nodes[j], h), this.at(nodes[k], h)];
          const g = new THREE.BufferGeometry().setFromPoints(pts);
          g.setIndex([0, 1, 2]);
          this.cutGroup.add(new THREE.Mesh(g, new THREE.MeshBasicMaterial({ color: VERMILION, transparent: true, opacity: 0.4, side: THREE.DoubleSide, depthWrite: false })));
        }
    const dot = new THREE.SphereGeometry(0.16, 10, 8);
    for (const b of dots) {
      const p = this.at(b, h);
      if (!p) continue;
      const s = new THREE.Mesh(dot, new THREE.MeshBasicMaterial({ color: IVORY }));
      s.position.copy(p);
      this.cutGroup.add(s);
    }
    const div = document.createElement('div');
    div.className = 'foam-cut';
    div.textContent = view.dateText || view.date;
    const label = new CSS2DObject(div);
    label.position.set(-(RING + 4.5), h, 0);
    this.cutGroup.add(label);
  }

  // ---- the pointer ----------------------------------------------------------------------------------
  pick(e) {
    const r = this.renderer.domElement.getBoundingClientRect();
    this.pointer.set(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
    this.raycaster.setFromCamera(this.pointer, this.camera);
    const hits = this.raycaster.intersectObjects(this.pickables.filter((m) => m.visible), false);
    // Knots and threads before the translucent ribbons behind them.
    hits.sort((a, b) => (a.object.userData.kind === 'lens') - (b.object.userData.kind === 'lens') || a.distance - b.distance);
    return hits[0]?.object.userData || null;
  }
  hover(e) {
    const hit = this.pick(e);
    this.renderer.domElement.style.cursor = hit ? 'pointer' : 'grab';
    if (!hit) return this.onHover?.(null);
    if (hit.kind === 'lens') this.onHover?.({ kind: 'lens', face: hit.face });
    else if (hit.kind === 'triangle') this.onHover?.({ kind: 'triangle', triangle: hit.triangle });
    else this.onHover?.(hit);
  }
  click(e) {
    const hit = this.pick(e);
    const d = this.data;
    if (!hit || !d) return;
    if (hit.kind === 'being') this.onPick?.({ kind: 'being', id: d.beings[hit.being].id });
    else if (hit.kind === 'occasion' || hit.kind === 'triangle') this.onPick?.({ kind: 'occasion', id: d.occasions[hit.occasion].id });
    else if (hit.kind === 'lens') this.onPick?.({ kind: 'tie', id: d.ties[hit.tie].id });
  }

  // Turns the camera slowly round the history, or stops.
  spin(on) {
    this.controls.autoRotate = on;
    this.controls.autoRotateSpeed = 0.6;
  }
  reset() {
    this.camera.position.set(44, HEIGHT * 0.85, 50);
    this.controls.target.set(0, HEIGHT * 0.48, 0);
    this.controls.update();
    this.draw();
  }
}
