// A city in three dimensions with three.js: the scene of a city-scale GIS on the local frame (frame.js), its
// layers of features drawn as extruded solids (buildings), flat areas (parks, water), ribbons of a width
// (streets), points (holdings) and bars (a quantity at a place), lit by the Sun at a direction the caller
// gives (the real Sun of an instant: the astronomy domain's sunAltAz), with map controls, picking and hover.
// For what moves: slabs (water standing in a street's channel to a level, which can vary along it), swarms
// (walkers, placed afresh every frame), arrows lying on the ground (currents), arcs over the city (two
// places joined) and tubes along paths in space (worldlines, with time drawn as height), each restyled
// without being rebuilt.
//
// three.js and its MapControls are passed in (cityScene(THREE, { MapControls }, …)), as the globe's renderer
// takes three.js, so that the arithmetic (geometry.js, style.js) stays checkable in Node and another backend
// can be written beside this one. Local east, north, up (x, y, z) are three.js's (x, z, −y): a rotation, so
// orientation is kept, and a counter-clockwise ring seen from above stays counter-clockwise.

import { edges, extent, ribbon } from './geometry.js';

// sRGB to linear light, for colours written straight into buffers (the scene works in linear sRGB).
const linear = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);

export function cityScene(THREE, { MapControls }, container, options = {}) {
  const opts = { background: '#f3ece4', ground: '#ebe3d9', shadows: true, ...options };
  const renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: !!opts.preserveDrawingBuffer });
  renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
  renderer.shadowMap.enabled = !!opts.shadows;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 0.85;
  container.append(renderer.domElement);
  renderer.domElement.style.display = 'block';
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(opts.background);
  const camera = new THREE.PerspectiveCamera(42, 1, 1, 60000);
  const controls = new MapControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.12;
  controls.screenSpacePanning = false;
  controls.maxPolarAngle = Math.PI * 0.47;
  controls.minDistance = 30;
  controls.maxDistance = 20000;

  const hemi = new THREE.HemisphereLight(0xfff4ea, 0x8a7a70, 1.35);
  scene.add(hemi);
  const sun = new THREE.DirectionalLight(0xfff0dc, 2.2);
  sun.castShadow = !!opts.shadows;
  sun.shadow.mapSize.set(4096, 4096);
  sun.shadow.bias = -0.0004;
  sun.shadow.normalBias = 0.6;
  scene.add(sun, sun.target);

  const groundMat = new THREE.MeshStandardMaterial({ color: opts.ground, roughness: 1, metalness: 0 });
  const ground = new THREE.Mesh(new THREE.PlaneGeometry(1, 1), groundMat);
  ground.rotation.x = -Math.PI / 2;
  ground.receiveShadow = true;
  scene.add(ground);

  const toWorld = (x, y, z = 0) => new THREE.Vector3(x, z, -y);
  const layers = new Map();   // id → { kind, object, features, … }
  let bounds = { x0: -500, y0: -500, x1: 500, y1: 500 };
  let sunDir = { azimuth: (210 * Math.PI) / 180, altitude: (45 * Math.PI) / 180 };

  function fit(b, { tilt = 52, heading = 0 } = {}) {
    bounds = b;
    const cx = (b.x0 + b.x1) / 2, cy = (b.y0 + b.y1) / 2, size = Math.max(b.x1 - b.x0, b.y1 - b.y0, 200);
    ground.scale.set(size * 6, size * 6, 1);
    ground.position.set(cx, -0.05, -cy);
    const range = size * (opts.range ?? 0.62), t = (tilt * Math.PI) / 180, h = (heading * Math.PI) / 180;
    controls.target.copy(toWorld(cx, cy, 0));
    camera.position.copy(toWorld(cx - range * Math.sin(t) * Math.sin(h), cy - range * Math.sin(t) * Math.cos(h), range * Math.cos(t)));
    camera.far = size * 20;
    camera.updateProjectionMatrix();
    scene.fog = new THREE.Fog(opts.background, size * 2.2, size * 7);
    setSun(sunDir);
    controls.update();
  }

  // The Sun's light from an azimuth (from north through east) and an altitude, radians; below the horizon the
  // street is lit by the sky alone.
  function setSun({ azimuth, altitude }) {
    sunDir = { azimuth, altitude };
    const b = bounds, cx = (b.x0 + b.x1) / 2, cy = (b.y0 + b.y1) / 2, size = Math.max(b.x1 - b.x0, b.y1 - b.y0, 200);
    const up = Math.max(altitude, 0.02);
    const d = [Math.sin(azimuth) * Math.cos(up), Math.cos(azimuth) * Math.cos(up), Math.sin(up)];
    const dist = size * 1.5;
    sun.position.copy(toWorld(cx + d[0] * dist, cy + d[1] * dist, d[2] * dist));
    sun.target.position.copy(toWorld(cx, cy, 0));
    const cam = sun.shadow.camera;
    cam.left = -size * 0.75, cam.right = size * 0.75, cam.top = size * 0.75, cam.bottom = -size * 0.75;
    cam.near = 1, cam.far = dist * 2.5;
    cam.updateProjectionMatrix();
    const day = altitude > 0 ? Math.min(1, altitude / 0.35) : 0;
    sun.intensity = 1.9 * day;
    sun.castShadow = !!opts.shadows && day > 0;
    hemi.intensity = 0.35 + 0.75 * Math.max(day, altitude > -0.1 ? 0.35 : 0.15);
  }

  // ---- layers ---------------------------------------------------------------------------------------------

  function remove(id) {
    const l = layers.get(id);
    if (!l) return;
    scene.remove(l.object);
    l.object.traverse((o) => { if (o.geometry) o.geometry.dispose(); if (o.material) [].concat(o.material).forEach((m) => m.dispose()); });
    layers.delete(id);
  }

  function writeColour(attr, start, count, rgb) {
    for (let v = start; v < start + count; v++) attr.setXYZ(v, rgb[0], rgb[1], rgb[2]);
  }

  // Extruded solids: features with footprint (outer ring, counter-clockwise), holes, a base and a top (m).
  function addSolids(id, features, s) {
    remove(id);
    const pos = [], col = [], ranges = [], faceFeature = [];
    const tri = (a, b, c, k) => { pos.push(...a, ...b, ...c); faceFeature.push(k); };
    features.forEach((f, k) => {
      const outer = s.ring(f), holes = (s.holes ? s.holes(f) : []) || [];
      if (!outer || outer.length < 3) { ranges.push([pos.length / 3, 0]); return; }
      const z0 = s.base ? s.base(f) || 0 : 0, z1 = Math.max(z0 + 0.5, s.height(f));
      const start = pos.length / 3;
      const contour = outer.map(([x, y]) => new THREE.Vector2(x, y));
      const hv = holes.filter((h) => h.length >= 3).map((h) => h.map(([x, y]) => new THREE.Vector2(x, y)));
      const all = [...outer, ...holes.filter((h) => h.length >= 3).flat()];
      for (const [a, b, c] of THREE.ShapeUtils.triangulateShape(contour, hv)) {
        const w = (i) => toWorld(all[i][0], all[i][1], z1).toArray();
        tri(w(a), w(b), w(c), k);
      }
      for (const ring of [outer, ...holes]) {
        for (let i = 0; i < ring.length; i++) {
          const [ax, ay] = ring[i], [bx, by] = ring[(i + 1) % ring.length];
          const A0 = toWorld(ax, ay, z0).toArray(), B0 = toWorld(bx, by, z0).toArray(), A1 = toWorld(ax, ay, z1).toArray(), B1 = toWorld(bx, by, z1).toArray();
          tri(A0, B0, B1, k);
          tri(A0, B1, A1, k);
        }
      }
      const count = pos.length / 3 - start;
      const rgb = s.color(f, k);
      for (let v = 0; v < count; v++) col.push(rgb[0], rgb[1], rgb[2]);
      ranges.push([start, count]);
    });
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
    g.computeVertexNormals();
    const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.88, metalness: 0, side: THREE.DoubleSide, transparent: (s.opacity ?? 1) < 1, opacity: s.opacity ?? 1 }));
    m.castShadow = true;
    m.receiveShadow = true;
    scene.add(m);
    layers.set(id, { kind: 'solids', object: m, features, style: s, ranges, faceFeature: Int32Array.from(faceFeature) });
    return layers.get(id);
  }

  // Flat areas at a small height: features with a ring (and holes).
  function addFlats(id, features, s) {
    remove(id);
    const pos = [], col = [], ranges = [], faceFeature = [];
    const z = s.z ?? 0.02;
    features.forEach((f, k) => {
      const outer = s.ring(f), holes = (s.holes ? s.holes(f) : []) || [];
      const start = pos.length / 3;
      if (outer && outer.length >= 3) {
        const all = [...outer, ...holes.flat()];
        for (const [a, b, c] of THREE.ShapeUtils.triangulateShape(outer.map(([x, y]) => new THREE.Vector2(x, y)), holes.map((h) => h.map(([x, y]) => new THREE.Vector2(x, y))))) {
          for (const i of [a, b, c]) pos.push(...toWorld(all[i][0], all[i][1], z).toArray());
          faceFeature.push(k);
        }
      }
      const count = pos.length / 3 - start, rgb = s.color(f, k);
      for (let v = 0; v < count; v++) col.push(rgb[0], rgb[1], rgb[2]);
      ranges.push([start, count]);
    });
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
    g.computeVertexNormals();
    const clear = (s.opacity ?? 1) < 1;
    const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 1, metalness: 0, side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: -1,
      transparent: clear, opacity: s.opacity ?? 1, depthWrite: !clear }));
    m.receiveShadow = true;
    m.renderOrder = s.order ?? 0;
    scene.add(m);
    layers.set(id, { kind: 'flats', object: m, features, style: s, ranges, faceFeature: Int32Array.from(faceFeature) });
    return layers.get(id);
  }

  // Ribbons of a width along polylines, at a small height (streets). With `along`, each feature is split into
  // pieces of `step` metres coloured by along(f, s) — a quantity varying along a street.
  function addRibbons(id, features, s) {
    remove(id);
    const pos = [], col = [], ranges = [], faceFeature = [], pieces = [];
    const z = s.z ?? 0.06;
    features.forEach((f, k) => {
      const line = s.line(f), width = s.width(f);
      const start = pos.length / 3;
      if (line && line.length >= 2) {
        const t = ribbon(line, width);
        for (let i = 0; i < t.length; i += 2) pos.push(...toWorld(t[i], t[i + 1], z).toArray());
        for (let i = 0; i < t.length / 6; i++) faceFeature.push(k);
      }
      const count = pos.length / 3 - start, rgb = s.color(f, k);
      for (let v = 0; v < count; v++) col.push(rgb[0], rgb[1], rgb[2]);
      ranges.push([start, count]);
      pieces.push(line ? line.length - 1 : 0);
    });
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
    g.computeVertexNormals();
    const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.95, metalness: 0, side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: -2,
      transparent: (s.opacity ?? 1) < 1, opacity: s.opacity ?? 1 }));
    m.receiveShadow = true;
    m.renderOrder = s.order ?? 1;
    scene.add(m);
    layers.set(id, { kind: 'ribbons', object: m, features, style: s, ranges, faceFeature: Int32Array.from(faceFeature), pieces });
    return layers.get(id);
  }

  // Columns at points: features with a position; radius, height and colour per feature.
  function addPoints(id, features, s) {
    remove(id);
    const geo = s.shape === 'box' ? new THREE.BoxGeometry(1, 1, 1) : new THREE.CylinderGeometry(0.5, 0.5, 1, 14);
    const m = new THREE.InstancedMesh(geo, new THREE.MeshStandardMaterial({ roughness: 0.7, metalness: 0, transparent: (s.opacity ?? 1) < 1, opacity: s.opacity ?? 1 }), Math.max(1, features.length));
    m.count = features.length;
    m.castShadow = !!s.shadow;
    m.receiveShadow = true;
    const layer = { kind: 'points', object: m, features, style: s };
    layers.set(id, layer);
    scene.add(m);
    placePoints(layer);
    return layer;
  }

  function placePoints(layer) {
    const { object: m, features, style: s } = layer;
    const M = new THREE.Matrix4(), q = new THREE.Quaternion(), c = new THREE.Color();
    features.forEach((f, k) => {
      const [x, y] = s.position(f), r = s.radius ? s.radius(f, k) : 3, hgt = Math.max(0.01, s.height ? s.height(f, k) : 4), z0 = s.base ? s.base(f, k) : 0;
      M.compose(toWorld(x, y, z0 + hgt / 2), q, new THREE.Vector3(2 * r, hgt, 2 * r));
      m.setMatrixAt(k, M);
      const rgb = s.color(f, k) || [0.6, 0.6, 0.6];
      m.setColorAt(k, c.setRGB(rgb[0], rgb[1], rgb[2], THREE.SRGBColorSpace));
    });
    m.instanceMatrix.needsUpdate = true;
    if (m.instanceColor) m.instanceColor.needsUpdate = true;
    m.computeBoundingSphere();
  }

  // Restyle a layer without rebuilding its geometry: new colour (and, for points, height or radius) functions.
  function restyle(id, changes) {
    const l = layers.get(id);
    if (!l) return;
    Object.assign(l.style, changes);
    if (l.kind === 'points') return placePoints(l);
    if (l.kind === 'slabs') return writeSlabs(l);
    if (l.kind === 'arrows') return placeArrows(l);
    if (l.kind === 'arcs') return colourArcs(l);
    if (l.kind === 'swarm') return;
    const attr = l.object.geometry.getAttribute('color');
    l.features.forEach((f, k) => {
      const [start, count] = l.ranges[k];
      if (count) writeColour(attr, start, count, l.style.color(f, k));
    });
    attr.needsUpdate = true;
  }

  // ---- slabs, swarms, arrows and arcs ---------------------------------------------------------------------

  // Slabs: ribbons of a width along polylines, raised from base(f) (or z) to base + height, translucent, as
  // water stands in a street's channel. A feature's height is height(f), or heights(f)[i] at each vertex i
  // of its line; its colour color(f), or colors(f)[i]. restyle() moves the tops and recolours.
  function addSlabs(id, features, s) {
    remove(id);
    const pos = [], top = [], owner = [], station = [], faceFeature = [];
    features.forEach((f, k) => {
      const line = s.line(f);
      if (!line || line.length < 2) return;
      const { left, right } = edges(line, s.width(f));
      const n = line.length;
      const v = ([x, y], up, i) => { pos.push(x, 0, -y); top.push(up); owner.push(k); station.push(i); };
      const quad = (a, ua, ia, b, ub, ib, c, uc, ic, d, ud, id2) => {
        v(a, ua, ia), v(b, ub, ib), v(c, uc, ic), faceFeature.push(k);
        v(a, ua, ia), v(c, uc, ic), v(d, ud, id2), faceFeature.push(k);
      };
      for (let i = 1; i < n; i++) {
        quad(right[i - 1], 1, i - 1, right[i], 1, i, left[i], 1, i, left[i - 1], 1, i - 1);   // the surface
        quad(left[i - 1], 0, i - 1, left[i], 0, i, left[i], 1, i, left[i - 1], 1, i - 1);      // the left wall
        quad(right[i - 1], 0, i - 1, right[i], 0, i, right[i], 1, i, right[i - 1], 1, i - 1);  // the right wall
      }
      quad(left[0], 0, 0, right[0], 0, 0, right[0], 1, 0, left[0], 1, 0);                    // the ends
      quad(left[n - 1], 0, n - 1, right[n - 1], 0, n - 1, right[n - 1], 1, n - 1, left[n - 1], 1, n - 1);
    });
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.Float32BufferAttribute(new Float32Array(pos.length), 3));
    const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.32, metalness: 0.04, side: THREE.DoubleSide,
      transparent: true, opacity: s.opacity ?? 0.8, depthWrite: true }));
    m.receiveShadow = true;
    m.castShadow = !!s.shadow;
    m.renderOrder = s.order ?? 4;
    scene.add(m);
    const layer = { kind: 'slabs', object: m, features, style: s, top: Uint8Array.from(top), owner: Int32Array.from(owner), station: Int32Array.from(station),
      faceFeature: Int32Array.from(faceFeature) };
    layers.set(id, layer);
    writeSlabs(layer);
    return layer;
  }

  function writeSlabs(l) {
    const s = l.style, P = l.object.geometry.getAttribute('position'), C = l.object.geometry.getAttribute('color');
    const per = l.features.map((f, k) => ({
      base: s.base ? s.base(f, k) : (s.z ?? 0.15),
      h: s.heights ? s.heights(f, k) : Math.max(0.05, s.height ? s.height(f, k) : 1),
      c: s.colors ? s.colors(f, k) : s.color(f, k),
    }));
    const lin = new Map();
    const colour = (rgb) => { if (!lin.has(rgb)) lin.set(rgb, (rgb || [0.7, 0.7, 0.7]).map(linear)); return lin.get(rgb); };
    for (let v = 0; v < l.top.length; v++) {
      const q = per[l.owner[v]], i = l.station[v];
      const h = typeof q.h === 'number' ? q.h : Math.max(0.05, q.h[i] || 0);
      P.setY(v, q.base + (l.top[v] ? h : 0));
      const c = colour(s.colors ? q.c[i] : q.c);
      C.setXYZ(v, c[0], c[1], c[2]);
    }
    P.needsUpdate = true;
    C.needsUpdate = true;
    l.object.geometry.computeVertexNormals();
    l.object.geometry.computeBoundingSphere();
  }

  // A swarm: up to `capacity` walkers, small upright columns of a radius and height (m), placed by
  // setSwarm() every frame. Not picked.
  function addSwarm(id, s) {
    remove(id);
    const geo = new THREE.CylinderGeometry(0.5, 0.5, 1, 7);
    geo.translate(0, 0.5, 0);
    const capacity = Math.max(1, s.capacity || 1000);
    const m = new THREE.InstancedMesh(geo, new THREE.MeshStandardMaterial({ roughness: 0.6, metalness: 0 }), capacity);
    m.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(3 * capacity), 3);
    m.count = 0;
    m.frustumCulled = false;
    scene.add(m);
    const layer = { kind: 'swarm', object: m, features: [], style: { ...s, capacity, pick: false } };
    layers.set(id, layer);
    return layer;
  }

  // count walkers at xyz (east, north and the height of their feet, three numbers each), coloured rgb (sRGB in
  // [0, 1], three numbers each; omitted, the colours stay as they were).
  function setSwarm(id, count, xyz, rgb) {
    const l = layers.get(id);
    if (!l || l.kind !== 'swarm') return;
    const m = l.object, n = Math.min(count, l.style.capacity), a = m.instanceMatrix.array;
    const r = l.style.radius ?? 0.45, h = l.style.height ?? 1.7;
    for (let k = 0; k < n; k++) {
      const o = 16 * k;
      a[o] = 2 * r, a[o + 1] = 0, a[o + 2] = 0, a[o + 3] = 0;
      a[o + 4] = 0, a[o + 5] = h, a[o + 6] = 0, a[o + 7] = 0;
      a[o + 8] = 0, a[o + 9] = 0, a[o + 10] = 2 * r, a[o + 11] = 0;
      a[o + 12] = xyz[3 * k], a[o + 13] = xyz[3 * k + 2], a[o + 14] = -xyz[3 * k + 1], a[o + 15] = 1;
    }
    m.count = n;
    m.instanceMatrix.needsUpdate = true;
    if (rgb) {
      const c = m.instanceColor.array;
      for (let i = 0; i < 3 * n; i++) c[i] = linear(rgb[i]);
      m.instanceColor.needsUpdate = true;
    }
  }

  // Arrows lying on the ground: features with a position, a unit direction (east, north), a length and a
  // width (m), and a colour; at a height z (or base(f)).
  function addArrows(id, features, s) {
    remove(id);
    const shape = new THREE.Shape();
    shape.moveTo(-0.5, -0.16), shape.lineTo(0.1, -0.16), shape.lineTo(0.1, -0.5), shape.lineTo(0.5, 0);
    shape.lineTo(0.1, 0.5), shape.lineTo(0.1, 0.16), shape.lineTo(-0.5, 0.16), shape.closePath();
    const geo = new THREE.ShapeGeometry(shape);
    geo.rotateX(-Math.PI / 2);   // the shape's (x, y) becomes east and north
    const m = new THREE.InstancedMesh(geo, new THREE.MeshStandardMaterial({ roughness: 0.7, metalness: 0, side: THREE.DoubleSide, transparent: (s.opacity ?? 1) < 1,
      opacity: s.opacity ?? 1, polygonOffset: true, polygonOffsetFactor: -4 }), Math.max(1, features.length));
    m.count = features.length;
    m.renderOrder = s.order ?? 5;
    const layer = { kind: 'arrows', object: m, features, style: s };
    layers.set(id, layer);
    scene.add(m);
    placeArrows(layer);
    return layer;
  }

  function placeArrows(layer) {
    const { object: m, features, style: s } = layer;
    const M = new THREE.Matrix4(), q = new THREE.Quaternion(), c = new THREE.Color(), up = new THREE.Vector3(0, 1, 0);
    features.forEach((f, k) => {
      const [x, y] = s.position(f, k), [dx, dy] = s.direction(f, k);
      const len = Math.max(0.01, s.length(f, k)), wid = Math.max(0.01, s.width ? s.width(f, k) : len / 3);
      q.setFromAxisAngle(up, Math.atan2(dy, dx));
      M.compose(toWorld(x, y, s.base ? s.base(f, k) : (s.z ?? 0.4)), q, new THREE.Vector3(len, 1, wid));
      m.setMatrixAt(k, M);
      const rgb = s.color(f, k) || [0.6, 0.6, 0.6];
      m.setColorAt(k, c.setRGB(rgb[0], rgb[1], rgb[2], THREE.SRGBColorSpace));
    });
    m.instanceMatrix.needsUpdate = true;
    if (m.instanceColor) m.instanceColor.needsUpdate = true;
    m.computeBoundingSphere();
  }

  // Arcs over the city joining two places: from(f), to(f) (east, north), rising height(f) above base (m) at
  // their middle, tubes of a radius, coloured; restyle() recolours.
  function addArcs(id, features, s) {
    remove(id);
    const pos = [], nor = [], ranges = [], faceFeature = [];
    features.forEach((f, k) => {
      const [ax, ay] = s.from(f, k), [bx, by] = s.to(f, k);
      const z0 = s.base ?? 2, rise = s.height ? s.height(f, k) : Math.hypot(bx - ax, by - ay) / 3;
      const curve = new THREE.QuadraticBezierCurve3(toWorld(ax, ay, z0), toWorld((ax + bx) / 2, (ay + by) / 2, z0 + 2 * rise), toWorld(bx, by, z0));
      const g = new THREE.TubeGeometry(curve, 40, Math.max(0.05, s.radius ? s.radius(f, k) : 1), 8, false).toNonIndexed();
      const P = g.getAttribute('position'), N = g.getAttribute('normal');
      ranges.push([pos.length / 3, P.count]);
      for (let i = 0; i < P.count; i++) pos.push(P.getX(i), P.getY(i), P.getZ(i)), nor.push(N.getX(i), N.getY(i), N.getZ(i));
      for (let t = 0; t < P.count / 3; t++) faceFeature.push(k);
      g.dispose();
    });
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute('normal', new THREE.Float32BufferAttribute(nor, 3));
    g.setAttribute('color', new THREE.Float32BufferAttribute(new Float32Array(pos.length), 3));
    const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.45, metalness: 0.05, transparent: (s.opacity ?? 1) < 1, opacity: s.opacity ?? 1 }));
    m.castShadow = !!s.shadow;
    scene.add(m);
    const layer = { kind: 'arcs', object: m, features, style: s, ranges, faceFeature: Int32Array.from(faceFeature) };
    layers.set(id, layer);
    colourArcs(layer);
    return layer;
  }

  // Tubes along polylines in space: path(f) a list of [east, north, up] (m), of a radius, coloured; for
  // worldlines and links between events drawn with time as height. restyle() recolours.
  function addTubes(id, features, s) {
    remove(id);
    const pos = [], nor = [], ranges = [], faceFeature = [];
    features.forEach((f, k) => {
      const path = (s.path(f, k) || []).filter((p) => p.every(Number.isFinite));
      if (path.length < 2) { ranges.push([pos.length / 3, 0]); return; }
      const curve = new THREE.CurvePath();
      for (let i = 1; i < path.length; i++) curve.add(new THREE.LineCurve3(toWorld(...path[i - 1]), toWorld(...path[i])));
      const g = new THREE.TubeGeometry(curve, Math.max(1, path.length - 1) * (s.pieces ?? 1), Math.max(0.05, s.radius ? s.radius(f, k) : 0.5), s.sides ?? 6, false).toNonIndexed();
      const P = g.getAttribute('position'), N = g.getAttribute('normal');
      ranges.push([pos.length / 3, P.count]);
      for (let i = 0; i < P.count; i++) pos.push(P.getX(i), P.getY(i), P.getZ(i)), nor.push(N.getX(i), N.getY(i), N.getZ(i));
      for (let t = 0; t < P.count / 3; t++) faceFeature.push(k);
      g.dispose();
    });
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute('normal', new THREE.Float32BufferAttribute(nor, 3));
    g.setAttribute('color', new THREE.Float32BufferAttribute(new Float32Array(pos.length), 3));
    const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.55, metalness: 0, transparent: (s.opacity ?? 1) < 1, opacity: s.opacity ?? 1 }));
    m.castShadow = !!s.shadow;
    m.renderOrder = s.order ?? 2;
    scene.add(m);
    const layer = { kind: 'arcs', object: m, features, style: s, ranges, faceFeature: Int32Array.from(faceFeature) };
    layers.set(id, layer);
    colourArcs(layer);
    return layer;
  }

  function colourArcs(l) {
    const C = l.object.geometry.getAttribute('color');
    l.features.forEach((f, k) => {
      const [start, count] = l.ranges[k];
      writeColour(C, start, count, (l.style.color(f, k) || [0.6, 0.6, 0.6]).map(linear));
    });
    C.needsUpdate = true;
  }

  function setVisible(id, on) { const l = layers.get(id); if (l) l.object.visible = on; }

  // The page's colours for the sky and the ground (a day and a night theme).
  function setPalette({ background, ground: g }) {
    if (background) { scene.background.set(background); if (scene.fog) scene.fog.color.set(background); }
    if (g) groundMat.color.set(g);
  }

  // ---- picking and hover ------------------------------------------------------------------------------------

  const raycaster = new THREE.Raycaster();
  const pointer = new THREE.Vector2();
  function hit(event) {
    const r = renderer.domElement.getBoundingClientRect();
    pointer.set(((event.clientX - r.left) / r.width) * 2 - 1, -((event.clientY - r.top) / r.height) * 2 + 1);
    raycaster.setFromCamera(pointer, camera);
    const pickable = [...layers.entries()].filter(([, l]) => l.object.visible && l.style.pick !== false);
    const hits = raycaster.intersectObjects(pickable.map(([, l]) => l.object), false);
    for (const h of hits) {
      const entry = pickable.find(([, l]) => l.object === h.object);
      if (!entry) continue;
      const [id, l] = entry;
      const index = l.kind === 'points' || l.kind === 'arrows' ? h.instanceId : l.faceFeature[h.faceIndex];
      if (index === undefined || index < 0) continue;
      const p = h.point;
      return { layer: id, index, feature: l.features[index], point: [p.x, -p.z, p.y] };
    }
    const g = raycaster.ray.intersectPlane(new THREE.Plane(new THREE.Vector3(0, 1, 0), 0), new THREE.Vector3());
    return g ? { layer: null, point: [g.x, -g.z, 0] } : null;
  }
  const pickers = [], hoverers = [];
  let down = null;
  renderer.domElement.addEventListener('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
  renderer.domElement.addEventListener('pointerup', (e) => {
    if (!down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 5) return;
    const h = hit(e);
    for (const f of pickers) f(h, e);
  });
  let hoverQueued = null;
  renderer.domElement.addEventListener('pointermove', (e) => {
    if (!hoverers.length || e.buttons) return;
    hoverQueued = e;
  });
  renderer.domElement.addEventListener('pointerleave', () => { for (const f of hoverers) f(null, null); });

  // ---- the camera ----------------------------------------------------------------------------------------------

  let flight = null;
  function flyTo([x, y], range = null) {
    const target = toWorld(x, y, 0);
    const offset = camera.position.clone().sub(controls.target);
    if (range) offset.setLength(range);
    flight = { t0: performance.now(), from: controls.target.clone(), to: target, camFrom: camera.position.clone(), camTo: target.clone().add(offset) };
  }

  // ---- the loop ------------------------------------------------------------------------------------------------

  function resize() {
    const w = container.clientWidth || 800, hgt = container.clientHeight || 600;
    renderer.setSize(w, hgt, false);
    renderer.domElement.style.width = '100%';
    renderer.domElement.style.height = '100%';
    camera.aspect = w / hgt;
    camera.updateProjectionMatrix();
  }
  const observer = new ResizeObserver(resize);
  observer.observe(container);
  resize();
  let running = true;
  const before = [];
  function frame(now) {
    if (!running) return;
    if (flight) {
      const u = Math.min(1, (now - flight.t0) / 900), e = u < 0.5 ? 2 * u * u : 1 - 2 * (1 - u) * (1 - u);
      controls.target.lerpVectors(flight.from, flight.to, e);
      camera.position.lerpVectors(flight.camFrom, flight.camTo, e);
      if (u >= 1) flight = null;
    }
    for (const f of before) f(now);
    controls.update();
    if (hoverQueued) {
      const ev = hoverQueued;
      hoverQueued = null;
      const h = hit(ev);
      for (const f of hoverers) f(h, ev);
    }
    renderer.render(scene, camera);
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);

  return {
    THREE, scene, camera, renderer, controls, layers,
    toWorld, fit, setSun, setPalette, addSolids, addFlats, addRibbons, addPoints, addSlabs, addSwarm, setSwarm, addArrows, addArcs, addTubes, restyle, setVisible, remove, flyTo,
    extentOf: (points) => extent(points),
    onPick: (f) => pickers.push(f),
    onHover: (f) => hoverers.push(f),
    // A function called before each frame with the time (ms); returns the function that stops it.
    onFrame: (f) => { before.push(f); return () => { const i = before.indexOf(f); if (i >= 0) before.splice(i, 1); }; },
    resize,
    dispose() {
      running = false;
      observer.disconnect();
      for (const id of [...layers.keys()]) remove(id);
      controls.dispose();
      renderer.dispose();
      renderer.domElement.remove();
    },
  };
}
