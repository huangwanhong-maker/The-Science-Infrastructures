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
//
// A scene may stand alone (its own sky, ground and map controls) or ride above the Earth (earth.js): with `driver`
// its canvas is transparent, its ground only takes shadows, its camera follows the view the driver holds, and
// `eventElement` receives the gestures. Plug-ins are unaware of the difference: a move of `camera` and
// `controls.target` followed by `controls.update()` becomes a move of that view.

import { edges, extent, inRing, ribbon } from './geometry.js';
import { createPrisms } from './prisms.js';

// sRGB to linear light, for colours written straight into buffers (the scene works in linear sRGB).
const linear = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);

export function cityScene(THREE, { MapControls }, container, options = {}) {
  const opts = { background: '#f3ece4', ground: '#ebe3d9', shadows: true, driver: null, eventElement: null, shadowColour: '#1a0d24', ...options };
  const riding = !!opts.driver;
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: riding, preserveDrawingBuffer: !!opts.preserveDrawingBuffer });
  if (riding) renderer.setClearColor(0x000000, 0);
  renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
  renderer.shadowMap.enabled = !!opts.shadows;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 0.85;
  container.append(renderer.domElement);
  renderer.domElement.style.display = 'block';
  const scene = new THREE.Scene();
  scene.background = riding ? null : new THREE.Color(opts.background);
  const camera = new THREE.PerspectiveCamera(42, 1, 1, 60000);
  let controls;
  if (riding) {
    controls = { target: new THREE.Vector3(), update() { if (controls.placed) { opts.driver.place(camera.position, controls.target); controls.placed = false; } }, dispose() {}, placed: false };
    // A plug-in's move: it sets the camera and the target, then asks for an update.
    const update = controls.update;
    controls.update = () => { controls.placed = true; update(); };
  } else {
    controls = new MapControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.12;
    controls.screenSpacePanning = false;
    controls.maxPolarAngle = Math.PI * 0.47;
    controls.minDistance = 30;
    controls.maxDistance = 20000;
  }

  const hemi = new THREE.HemisphereLight(0xfff4ea, 0x8a7a70, 1.35);
  scene.add(hemi);
  const sun = new THREE.DirectionalLight(0xfff0dc, 2.2);
  sun.castShadow = !!opts.shadows;
  sun.shadow.mapSize.set(4096, 4096);
  sun.shadow.bias = -0.0004;
  sun.shadow.normalBias = 0.6;
  scene.add(sun, sun.target);

  const groundMat = riding ? new THREE.ShadowMaterial({ color: opts.shadowColour, opacity: opts.shadowOpacity ?? 0.34 }) : new THREE.MeshStandardMaterial({ color: opts.ground, roughness: 1, metalness: 0 });
  const ground = new THREE.Mesh(new THREE.PlaneGeometry(1, 1), groundMat);
  ground.rotation.x = -Math.PI / 2;
  ground.receiveShadow = true;
  scene.add(ground);

  const toWorld = (x, y, z = 0) => new THREE.Vector3(x, z, -y);
  const layers = new Map();   // id → { kind, object, features, … }
  let bounds = { x0: -500, y0: -500, x1: 500, y1: 500 };
  let groundBox = null;
  let sunDir = { azimuth: (210 * Math.PI) / 180, altitude: (45 * Math.PI) / 180 };

  // The camera's distance is the box's size times `range` (the view's, else the scene's, else 0.62): a city of tall
  // towers asks for a farther camera than one of low blocks.
  function fit(b, { tilt = 52, heading = 0, range: far, ...fitOptions } = {}) {
    bounds = b;
    const cx = (b.x0 + b.x1) / 2, cy = (b.y0 + b.y1) / 2, size = Math.max(b.x1 - b.x0, b.y1 - b.y0, 200);
    if (!groundBox) {
      ground.scale.set(size * 6, size * 6, 1);
      ground.position.set(cx, -0.05, -cy);
    }
    setSun(sunDir);
    if (riding) {
      if (!fitOptions.skipView) opts.driver.fit(b, { tilt, heading, range: far ?? opts.range ?? 0.62, ...fitOptions });
      return;
    }
    const range = size * (far ?? opts.range ?? 0.62), t = (tilt * Math.PI) / 180, h = (heading * Math.PI) / 180;
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
    sun.position.copy(toWorld(cx + d[0] * dist, cy + d[1] * dist, (b.z || 0) + d[2] * dist));
    sun.target.position.copy(toWorld(cx, cy, b.z || 0));
    const cam = sun.shadow.camera;
    cam.left = -size * 0.75, cam.right = size * 0.75, cam.top = size * 0.75, cam.bottom = -size * 0.75;
    cam.near = 1, cam.far = dist * 2.5;
    cam.updateProjectionMatrix();
    const day = altitude > 0 ? Math.min(1, altitude / 0.35) : 0;
    const light = opts.lightScale ?? 1;
    sun.intensity = 1.9 * day * light;
    sun.castShadow = !!opts.shadows && day > 0;
    hemi.intensity = (0.35 + 0.75 * Math.max(day, altitude > -0.1 ? 0.35 : 0.15)) * light;
    // The windows come on as the Sun goes down.
    facade.windows.value = 1 - Math.min(1, Math.max(0, (altitude + 0.06) / 0.2));
  }
  const facade = { windows: { value: 0 }, dissolve: { value: new THREE.Vector3(0, 0, 0) } };
  const DISSOLVE = `
  // The edge of what is held of the world dissolves in a stipple, as an engraving trails off, rather than ending in the tiles' straight lines.
  float bayer4(vec2 p) { vec2 q = mod(floor(p), 4.0); return (mod(q.x * 3.0 + q.y * 7.0 + q.x * q.y * 5.0, 16.0) + 0.5) / 16.0; }
  void dissolveAt(vec3 at) {
    if (uDissolve.z <= 0.0) return;
    float k = smoothstep(0.62 * uDissolve.z, uDissolve.z, length(at.xz - uDissolve.xy));
    if (k > bayer4(gl_FragCoord.xy)) discard;
  }
`;

  // ---- layers ---------------------------------------------------------------------------------------------

  function remove(id) {
    const l = layers.get(id);
    if (!l) return;
    scene.remove(l.object);
    l.object.traverse((o) => { if (o.geometry) o.geometry.dispose(); if (o.material) [].concat(o.material).forEach((m) => m.dispose()); });
    if (l.kind === 'prisms') l.prisms.dispose();
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
    const mat = new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.88, metalness: 0, side: THREE.DoubleSide, transparent: (s.opacity ?? 1) < 1, opacity: s.opacity ?? 1 });
    if (s.facade) withFacade(mat, s.facade);
    const m = new THREE.Mesh(g, mat);
    m.castShadow = true;
    m.receiveShadow = true;
    scene.add(m);
    layers.set(id, { kind: 'solids', object: m, features, style: s, ranges, faceFeature: Int32Array.from(faceFeature) });
    return layers.get(id);
  }

  // Walls that read as walls: darker at the street and lighter toward the roof, panes in rows of floors, and, as the Sun
  // goes down, some of those panes lit. It is a drawing's convention, in the page's colours (the floors are a nominal
  // 3.4 m and the bays 3.6 m), and says nothing of a real building's windows.
  // With `built` (geometry made elsewhere, addBuiltSolids) the shader finds a wall from its own derivatives and measures up from
  // the building's foot, which on real ground is at no one height; `form` then scales the building about its foot (rising as it comes).
  function withFacade(material, { pane = [1.0, 0.72, 0.86], warm = [1.0, 0.86, 0.62], glow = 0.6, foot = 0.62, built = false, form = null } = {}) {
    material.onBeforeCompile = (shader) => {
      shader.uniforms.uWindows = facade.windows;
      shader.uniforms.uDissolve = facade.dissolve;
      if (built) shader.uniforms.uForm = form;
      shader.vertexShader = shader.vertexShader
        .replace('#include <common>', `#include <common>\nvarying vec3 vWall;\nvarying float vUp;\n${built ? 'attribute float aRise;\nuniform float uForm;' : 'varying vec3 vWallN;'}`)
        .replace('#include <begin_vertex>', built
          ? '#include <begin_vertex>\ntransformed.y -= (1.0 - uForm) * aRise;\nvWall = transformed;\nvUp = aRise * uForm;'
          : '#include <begin_vertex>\nvWall = position;\nvUp = position.y;\nvWallN = normal;');
      shader.fragmentShader = shader.fragmentShader
        .replace('#include <common>', `#include <common>\nuniform float uWindows;\nuniform vec3 uDissolve;\nvarying vec3 vWall;\nvarying float vUp;\n${built ? DISSOLVE : 'varying vec3 vWallN;'}`)
        .replace('#include <emissivemap_fragment>', `#include <emissivemap_fragment>
  ${built ? 'dissolveAt(vWall);' : ''}
  {
    vec3 wn = ${built ? 'normalize(cross(dFdx(vWall), dFdy(vWall)))' : 'normalize(vWallN)'};
    if (abs(wn.y) < 0.5) {
      diffuseColor.rgb *= mix(${foot.toFixed(2)}, 1.08, smoothstep(0.0, 48.0, vUp));
      float along = dot(vWall, vec3(abs(wn.z), 0.0, abs(wn.x)));
      vec2 grid = vec2(along / 3.6, vUp / 3.4);
      // Panes too small to see give way to their average, so that a far wall does not shimmer.
      float near = 1.0 - smoothstep(0.3, 0.8, max(fwidth(grid.x), fwidth(grid.y)));
      vec2 cell = floor(grid), f = fract(grid);
      float pane = step(0.26, f.x) * step(f.x, 0.74) * step(0.30, f.y) * step(f.y, 0.76) * step(1.0, cell.y);
      float r = fract(sin(dot(cell, vec2(12.9898, 78.233)) + floor(vWall.x * 0.03 + vWall.z * 0.017) * 3.7) * 43758.5453);
      diffuseColor.rgb *= 1.0 - 0.13 * mix(0.2, pane, near);
      float lit = step(1.0 - 0.28 * uWindows, r);
      vec3 tone = mix(vec3(${pane.map((v) => v.toFixed(2)).join(', ')}), vec3(${warm.map((v) => v.toFixed(2)).join(', ')}), mix(0.5, fract(r * 7.0), near));
      totalEmissiveRadiance += tone * mix(0.056 * uWindows, pane * lit * uWindows, near) * ${glow.toFixed(2)};
    } else {
      diffuseColor.rgb *= 1.14;
    }
  }`);
    };
    material.customProgramCacheKey = () => `facade:${built}:${pane}:${warm}:${glow}:${foot}`;
  }

  // ---- geometry made elsewhere (a worker's, world-worker.js): typed arrays in the scene's axes (x east, y up, z south) ----

  const builtGeometry = (THREEattrs, index, sphere) => {
    const g = new THREE.BufferGeometry();
    for (const [name, attr] of THREEattrs) g.setAttribute(name, attr);
    g.setIndex(new THREE.BufferAttribute(index, 1));
    g.boundingSphere = new THREE.Sphere(new THREE.Vector3(...sphere.centre), sphere.radius);
    return g;
  };

  // Buildings: position, colour (bytes, four to a point), rise (each point's height above its building's foot), and for picking
  // each building's boxes [cx, cy, radius, z0, z1], rings and starts, and features. They rise from the ground as they come.
  function addBuiltSolids(id, b, s = {}) {
    remove(id);
    const g = builtGeometry([['position', new THREE.BufferAttribute(b.position, 3)], ['color', new THREE.BufferAttribute(b.colour, 4, true)], ['aRise', new THREE.BufferAttribute(b.rise, 1)]], b.index, b.sphere);
    const form = { value: s.grow === false ? (s.form ?? 1) : 0.02 };
    const mat = new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.88, metalness: 0, side: THREE.DoubleSide, flatShading: true });
    withFacade(mat, { ...(s.facade || {}), built: true, form });
    mat.userData.form = form;
    const m = new THREE.Mesh(g, mat);
    m.castShadow = s.shadows !== false;
    m.receiveShadow = true;
    m.renderOrder = 0;
    scene.add(m);
    layers.set(id, { kind: 'built-solids', object: m, features: b.features, style: { pick: s.pick }, built: b, form: s.form ?? 1, born: s.grow === false ? -1e9 : performance.now() });
    return layers.get(id);
  }

  // What is painted on the ground: drawn before everything that stands on it, in their order, writing no depth.
  function addBuiltPaint(id, b, s = {}) {
    remove(id);
    const g = builtGeometry([['position', new THREE.BufferAttribute(b.position, 3)], ['color', new THREE.BufferAttribute(b.colour, 4, true)]], b.index, b.sphere);
    const paint = new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 1, metalness: 0, side: THREE.DoubleSide, flatShading: true, depthWrite: false, depthTest: s.depthTest ?? false });
    paint.onBeforeCompile = (shader) => {
      shader.uniforms.uDissolve = facade.dissolve;
      shader.vertexShader = shader.vertexShader.replace('#include <common>', '#include <common>\nvarying vec3 vWall;').replace('#include <begin_vertex>', '#include <begin_vertex>\nvWall = position;');
      shader.fragmentShader = shader.fragmentShader
        .replace('#include <common>', `#include <common>\nuniform vec3 uDissolve;\nvarying vec3 vWall;\n${DISSOLVE}`)
        .replace('#include <emissivemap_fragment>', '#include <emissivemap_fragment>\ndissolveAt(vWall);');
    };
    paint.customProgramCacheKey = () => 'built-paint';
    const m = new THREE.Mesh(g, paint);
    m.renderOrder = s.order ?? -2;
    m.receiveShadow = true;
    scene.add(m);
    layers.set(id, { kind: 'built-paint', object: m, style: { pick: false } });
    return layers.get(id);
  }

  // A wave along ways: each way a ribbon whose points carry a time, drawn up to `setWave`'s now, with a crest at its front. For what spreads
  // over a network: the minutes from a place, the path of a flood.
  function addWave(id, lines, times, s) {
    remove(id);
    let vertices = 0, indices = 0;
    lines.forEach((l, k) => {
      const n = l.xyz.length / 3;
      if (n >= 2 && times[k].some((t) => t !== null)) { vertices += 2 * n; indices += 6 * (n - 1); }
    });
    if (!vertices) return null;
    const position = new Float32Array(3 * vertices), time = new Float32Array(vertices), index = new Uint32Array(indices);
    let v = 0, t = 0;
    const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
    lines.forEach((l, k) => {
      const n = l.xyz.length / 3;
      if (n < 2 || !times[k].some((q) => q !== null)) return;
      const line = [];
      for (let i = 0; i < n; i++) line.push([l.xyz[3 * i], l.xyz[3 * i + 1]]);
      const { left, right } = edges(line, s.width(l));
      const start = v;
      for (let i = 0; i < n; i++) {
        const z = l.xyz[3 * i + 2] + 0.25, tm = times[k][i] === null ? 1e9 : times[k][i];
        for (const p of [left[i], right[i]]) {
          position[3 * v] = p[0], position[3 * v + 1] = z, position[3 * v + 2] = -p[1];
          time[v] = tm;
          min[0] = Math.min(min[0], p[0]), max[0] = Math.max(max[0], p[0]), min[1] = Math.min(min[1], p[1]), max[1] = Math.max(max[1], p[1]), min[2] = Math.min(min[2], z), max[2] = Math.max(max[2], z);
          v++;
        }
      }
      for (let i = 0; i < n - 1; i++) {
        const a = start + 2 * i;
        index[t++] = a, index[t++] = a + 1, index[t++] = a + 3, index[t++] = a, index[t++] = a + 3, index[t++] = a + 2;
      }
    });
    const c = [(min[0] + max[0]) / 2, (min[2] + max[2]) / 2, -(min[1] + max[1]) / 2];
    const g = builtGeometry([['position', new THREE.BufferAttribute(position, 3)], ['aTime', new THREE.BufferAttribute(time, 1)]], index, { centre: c, radius: Math.hypot(max[0] - min[0], max[1] - min[1], max[2] - min[2]) / 2 + 1 });
    const linearOf = (rgb) => new THREE.Color().setRGB(rgb[0], rgb[1], rgb[2], THREE.SRGBColorSpace);
    const mat = new THREE.ShaderMaterial({
      uniforms: { uNow: { value: 0 }, uUpto: { value: s.upto }, uNear: { value: linearOf(s.colours.near) }, uFar: { value: linearOf(s.colours.far) }, uCrest: { value: linearOf(s.colours.crest) } },
      vertexShader: 'attribute float aTime; varying float vTime; void main() { vTime = aTime; gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0); }',
      fragmentShader: `uniform float uNow; uniform float uUpto; uniform vec3 uNear; uniform vec3 uFar; uniform vec3 uCrest; varying float vTime;
        void main() {
          if (vTime > uNow) discard;
          vec3 c = mix(uNear, uFar, smoothstep(0.0, 1.0, vTime / uUpto));
          float front = exp(-(uNow - vTime) / (uUpto * 0.018));
          gl_FragColor = vec4(mix(c, uCrest, front), 1.0);
          #include <colorspace_fragment>
        }`,
      depthWrite: false, side: THREE.DoubleSide,
    });
    const m = new THREE.Mesh(g, mat);
    m.renderOrder = -1;
    scene.add(m);
    layers.set(id, { kind: 'built-wave', object: m, style: { pick: false } });
    return layers.get(id);
  }
  function setWave(id, now) { const l = layers.get(id); if (l && l.kind === 'built-wave') l.object.material.uniforms.uNow.value = now; }

  // The ground under a tile, unseen: it takes shadows, and hides what stands behind a ridge (its depth a few metres low, so that
  // what lies on the ground is never cut by the coarseness of its grid).
  function addBuiltGround(id, b) {
    remove(id);
    const g = builtGeometry([['position', new THREE.BufferAttribute(b.position, 3)]], b.index, b.sphere);
    const group = new THREE.Group();
    const depth = new THREE.Mesh(g, new THREE.MeshBasicMaterial({ colorWrite: false, side: THREE.DoubleSide }));
    depth.position.y = -6;
    depth.renderOrder = -3;
    const shade = new THREE.Mesh(g, new THREE.ShadowMaterial({ color: opts.shadowColour, opacity: opts.shadowOpacity ?? 0.34, depthWrite: false }));
    shade.receiveShadow = true;
    shade.renderOrder = 2;
    shade.position.y = 0.15;
    group.add(depth, shade);
    scene.add(group);
    layers.set(id, { kind: 'built-ground', object: group, style: { pick: false } });
    return layers.get(id);
  }

  // An object of one's own, in the scene's axes: it is kept and let go with the layers, and picked by whoever made it.
  function addObject(id, object) {
    remove(id);
    scene.add(object);
    layers.set(id, { kind: 'object', object, style: { pick: false } });
    return layers.get(id);
  }

  // A building under a ray, by its box and then its footprint: tens of thousands are asked at the speed of a mouse.
  function pickBuilt(ray) {
    const o = [ray.origin.x, -ray.origin.z, ray.origin.y], d = [ray.direction.x, -ray.direction.z, ray.direction.y];
    const flat = Math.hypot(d[0], d[1]);
    let best = null, bestT = Infinity;
    for (const [id, l] of layers) {
      if (l.kind !== 'built-solids' || !l.object.visible) continue;
      const { boxes, rings, starts } = l.built, count = boxes.length / 5, form = l.object.material.userData.form.value;
      for (let k = 0; k < count; k++) {
        const cx = boxes[5 * k], cy = boxes[5 * k + 1], radius = boxes[5 * k + 2], z0 = boxes[5 * k + 3], z1 = z0 + (boxes[5 * k + 4] - z0) * form;
        let t0 = 0, t1 = Infinity;
        if (flat > 1e-6) {
          const ex = cx - o[0], ey = cy - o[1];
          const across = (ex * d[1] - ey * d[0]) / flat;
          if (Math.abs(across) > radius) continue;
          const mid = (ex * d[0] + ey * d[1]) / (flat * flat), half = Math.sqrt(radius * radius - across * across) / flat;
          t0 = mid - half, t1 = mid + half;
        } else if (Math.hypot(cx - o[0], cy - o[1]) > radius) continue;
        if (Math.abs(d[2]) > 1e-9) {
          const a = (z0 - o[2]) / d[2], b = (z1 - o[2]) / d[2];
          t0 = Math.max(t0, Math.min(a, b)), t1 = Math.min(t1, Math.max(a, b));
        }
        if (t0 >= t1 || t0 >= bestT || t1 <= 0) continue;
        const start = starts[2 * k], n = starts[2 * k + 1];
        const ring = [];
        for (let i = 0; i < n; i++) ring.push([rings[2 * (start + i)], rings[2 * (start + i) + 1]]);
        const steps = 14;
        for (let i = 0; i <= steps; i++) {
          const t = Math.max(0, t0) + ((Math.min(t1, bestT) - Math.max(0, t0)) * i) / steps;
          if (inRing([o[0] + d[0] * t, o[1] + d[1] * t], ring)) {
            if (t < bestT) { bestT = t; best = { layer: id, index: k, feature: l.features[k], point: [o[0] + d[0] * t, o[1] + d[1] * t, o[2] + d[2] * t], distance: t }; }
            break;
          }
        }
      }
    }
    return best;
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
    if (l.kind === 'prisms') return writePrisms(l);
    if (l.kind === 'raster') return writeRaster(l);
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

  // Columns on cells: features with a ring (local east, north), height(f, k) in metres and color(f, k), or null for a cell
  // not drawn; drawn by one mesh and a table (prisms.js), so that restyle() takes a new year at no cost. They rise as they arrive.
  function addPrisms(id, features, s) {
    remove(id);
    const prisms = createPrisms(THREE, features.map((f) => ({ ring: s.ring(f).map(([x, y]) => toWorld(x, y, s.z ?? 0.3).toArray()), up: [0, 1, 0] })), { opacity: s.opacity ?? 1, shrink: s.shrink ?? 0.9, order: s.order ?? 5 });
    scene.add(prisms.object);
    const layer = { kind: 'prisms', object: prisms.object, prisms, features, style: s, born: performance.now() };
    layers.set(id, layer);
    writePrisms(layer);
    return layer;
  }
  function writePrisms(l) {
    const s = l.style;
    l.prisms.set(l.features.map((f, k) => s.height(f, k)), l.features.map((f, k) => s.color(f, k)));
  }
  // A raster laid on the city: an image (its first row the south) over corners in the local frame, above the streets
  // and beneath the buildings. With `ramp` the image holds values in its red channel and is coloured here, and `cut`
  // hides values above it; restyle() moves `cut`, `floor` and `opacity`.
  function addRaster(id, s) {
    remove(id);
    const [x0, y0, x1, y1] = s.extent;
    const map = new THREE.Texture(s.image);
    map.flipY = false;
    map.colorSpace = s.ramp ? THREE.NoColorSpace : THREE.SRGBColorSpace;
    map.generateMipmaps = false;
    map.minFilter = map.magFilter = s.nearest ? THREE.NearestFilter : THREE.LinearFilter;
    map.needsUpdate = true;
    let ramp = null;
    if (s.ramp) {
      const data = new Uint8Array(4 * 256);
      for (let i = 0; i < 256; i++) {
        const c = s.ramp(i / 255) || [0, 0, 0, 0];
        data[4 * i] = Math.round(c[0] * 255), data[4 * i + 1] = Math.round(c[1] * 255), data[4 * i + 2] = Math.round(c[2] * 255), data[4 * i + 3] = Math.round((c[3] ?? 1) * 255);
      }
      ramp = new THREE.DataTexture(data, 256, 1, THREE.RGBAFormat);
      ramp.colorSpace = THREE.SRGBColorSpace;
      ramp.minFilter = ramp.magFilter = THREE.LinearFilter;
      ramp.needsUpdate = true;
    }
    const g = new THREE.PlaneGeometry(x1 - x0, y1 - y0);
    g.rotateX(-Math.PI / 2);   // the plane's (x, −z) becomes east and north
    const uv = g.getAttribute('uv');
    for (let i = 0; i < uv.count; i++) uv.setY(i, uv.getY(i));
    const m = new THREE.ShaderMaterial({
      transparent: true, depthWrite: false, side: THREE.DoubleSide,
      vertexShader: 'varying vec2 vUv; void main() { vUv = uv; gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0); }',
      fragmentShader: `uniform sampler2D map; uniform sampler2D ramp; uniform float hasRamp; uniform float opacity; uniform float cut; uniform float floorAt; varying vec2 vUv;
void main() {
  vec4 t = texture2D(map, vUv);
  vec4 c = t;
  if (hasRamp > 0.5) { if (t.r <= floorAt || t.r > cut) discard; c = texture2D(ramp, vec2(t.r, 0.5)); }
  gl_FragColor = vec4(c.rgb, c.a * opacity);
  #include <colorspace_fragment>
}`,
      uniforms: { map: { value: map }, ramp: { value: ramp || map }, hasRamp: { value: ramp ? 1 : 0 }, opacity: { value: s.opacity ?? 0.8 }, cut: { value: s.cut ?? 1 }, floorAt: { value: s.floor ?? 0 } },
    });
    const mesh = new THREE.Mesh(g, m);
    mesh.position.copy(toWorld((x0 + x1) / 2, (y0 + y1) / 2, s.z ?? 0.25));
    mesh.renderOrder = s.order ?? 3;
    scene.add(mesh);
    const layer = { kind: 'raster', object: mesh, features: [], style: { ...s, pick: false } };
    layers.set(id, layer);
    return layer;
  }
  function writeRaster(l) {
    const u = l.object.material.uniforms, s = l.style;
    if (s.cut !== undefined) u.cut.value = s.cut;
    if (s.floor !== undefined) u.floorAt.value = s.floor;
    if (s.opacity !== undefined) u.opacity.value = s.opacity;
  }

  function setVisible(id, on) { const l = layers.get(id); if (l) l.object.visible = on; }
  // A layer raised above the ground (m): the sheets of an exploded view.
  function setLift(id, z) {
    const l = layers.get(id);
    if (!l) return;
    if (l.kind === 'raster') { const zz = l.style.z ?? 0.25; l.object.position.y = zz + z; } else l.object.position.y = z;
  }

  // The page's colours for the sky and the ground (a day and a night theme).
  function setPalette({ background, ground: g }) {
    if (riding) return;
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
    const pickable = [...layers.entries()].filter(([, l]) => l.kind !== 'prisms' && !l.kind.startsWith('built') && l.object.visible && l.style.pick !== false);
    const hits = raycaster.intersectObjects(pickable.map(([, l]) => l.object), false);
    const built = pickBuilt(raycaster.ray);
    if (built && !(hits[0] && hits[0].distance < built.distance)) return built;
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
    if (!g) return null;
    // Columns are picked by the cell under the pointer on the ground: their mesh moves in its shader.
    for (const [id, l] of layers) {
      if (l.kind !== 'prisms' || !l.object.visible || l.style.pick === false) continue;
      const k = l.features.findIndex((f) => inRing([g.x, -g.z], l.style.ring(f)));
      if (k >= 0 && l.style.height(l.features[k], k) >= 0) { l.prisms.setPicked(k); return { layer: id, index: k, feature: l.features[k], point: [g.x, -g.z, 0] }; }
      l.prisms.setPicked(-1);
    }
    return { layer: null, point: [g.x, -g.z, 0] };
  }
  const pickers = [], hoverers = [];
  let down = null;
  const events = opts.eventElement || renderer.domElement;
  const listen = [];
  const on = (name, f) => { events.addEventListener(name, f); listen.push([name, f]); };
  on('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
  on('pointerup', (e) => {
    if (!down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 5) return;
    const h = state.paused ? null : hit(e);
    for (const f of pickers) f(h, e);
  });
  let hoverQueued = null;
  on('pointermove', (e) => {
    if (!hoverers.length || e.buttons) return;
    hoverQueued = e;
  });
  on('pointerleave', () => { for (const f of hoverers) f(null, null); });

  // ---- the camera ----------------------------------------------------------------------------------------------

  let flight = null;
  // Fly to a point, keeping the camera's heading; at `range` metres if given, and at `tilt` degrees from the
  // vertical if given (else the camera's own).
  function flyTo([x, y], range = null, tilt = null) {
    if (riding) return opts.driver.flyTo([x, y], range, tilt);
    const target = toWorld(x, y, 0);
    const offset = camera.position.clone().sub(controls.target);
    if (tilt !== null) {
      const t = (tilt * Math.PI) / 180, length = offset.length(), flat = Math.hypot(offset.x, offset.z) || 1;
      offset.set((offset.x / flat) * Math.sin(t) * length, Math.cos(t) * length, (offset.z / flat) * Math.sin(t) * length);
    }
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
  const state = { paused: false };   // while the city is out of sight (from far above), nothing of it is drawn
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
    if (riding) opts.driver.sync(camera, controls.target);
    else controls.update();
    for (const l of layers.values()) {
      if (l.kind === 'prisms') {
        const t = (now - l.born) / 1000;
        l.prisms.setGrow(t >= 1 ? 1 : 1 - (1 - Math.max(0, t)) ** 3);
      } else if (l.kind === 'built-solids' && l.born > 0) {
        const t = (now - l.born) / 1100;
        l.object.material.userData.form.value = l.form * (t >= 1 ? 1 : 1 - (1 - Math.max(0, t)) ** 3);
        if (t >= 1) l.born = -1;
      }
    }
    if (hoverQueued) {
      const ev = hoverQueued;
      hoverQueued = null;
      const h = state.paused ? null : hit(ev);
      for (const f of hoverers) f(h, ev);
    }
    if (!state.paused) renderer.render(scene, camera);
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);

  return {
    THREE, scene, camera, renderer, controls, layers,
    toWorld, fit, setSun, setPalette, addSolids, addFlats, addRibbons, addPoints, addPrisms, addRaster, addSlabs, addSwarm, setSwarm, addArrows, addArcs, addTubes, addBuiltSolids, addBuiltPaint, addBuiltGround, addWave, setWave, addObject, restyle, setVisible, setLift, remove, flyTo,
    // Where what is drawn of the world's tiles trails off: about a point of the local frame, within a radius (m); none at 0.
    setDissolve(x, y, radius) { facade.dissolve.value.set(x, -y, radius); },
    // The box the Sun's shadows are computed over, which may follow the view, apart from the box the camera was fitted to.
    setShadowBox(b) { bounds = b; setSun(sunDir); },
    // The flat ground that takes shadows: over a box of the local frame only, or nowhere (where each tile brings its own).
    setGroundBox(b) {
      ground.visible = !!b;
      groundBox = b;
      if (b) { ground.scale.set(b.x1 - b.x0, b.y1 - b.y0, 1); ground.position.set((b.x0 + b.x1) / 2, -0.05, -(b.y0 + b.y1) / 2); }
    },
    // The form of the buildings made elsewhere: 1 as built, lower to shrink them about their feet.
    setBuiltForm(id, f) { const l = layers.get(id); if (l && l.kind === 'built-solids') { l.form = f; if (l.born <= 0) l.object.material.userData.form.value = f; } },
    setPaused: (v) => { state.paused = v; },
    extentOf: (points) => extent(points),
    onPick: (f) => pickers.push(f),
    onHover: (f) => hoverers.push(f),
    // A function called before each frame with the time (ms); returns the function that stops it.
    onFrame: (f) => { before.push(f); return () => { const i = before.indexOf(f); if (i >= 0) before.splice(i, 1); }; },
    resize,
    dispose() {
      running = false;
      observer.disconnect();
      for (const [name, f] of listen) events.removeEventListener(name, f);
      for (const id of [...layers.keys()]) remove(id);
      controls.dispose();
      renderer.dispose();
      // Browsers keep few live contexts, and lose the oldest: a city left gives its own back at once.
      renderer.forceContextLoss();
      renderer.domElement.remove();
    },
  };
}
