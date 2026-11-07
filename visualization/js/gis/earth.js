// The Earth beneath a city: one view from the whole globe down to a street. The globe (terrain and imagery by
// level of detail, the shared layer's `globe`) is drawn on a canvas of its own, in the body-fixed frame with the
// camera at the origin and a logarithmic depth buffer; a city (renderer-three.js, in its local east, north, up
// frame) is drawn above it on a second, transparent canvas, with a camera derived each frame from the same view.
// One set of gestures, the globe's, moves both. The city keeps its own depth and its polygon offsets, and whatever
// a plug-in draws in the local frame needs nothing of the globe; the price is that a city is always in front of the
// ground, which is why the ground beneath it is made level (`flat`), feathered back into the real terrain.
//
// Beside the tiles this module holds what is drawn in the body-fixed frame: sets of columns (prisms.js) over a
// country or the whole Earth, placed relative to the camera in double precision, a raster draped on the ground, and
// a backdrop: the sky by the Sun's altitude near the ground, space from above. Colours are the application's.
//
// three.js is passed in, as the globe's renderer takes it.

import { createPrisms } from './prisms.js';
import { cartesian, figure, geodetic, normal, rayHit } from '../globe/figure.js';
import { cameraOf, GlobeControls, localFrame } from '../globe/controls.js';
import { IMAGERY, imageryLayer } from '../globe/imagery.js';
import { TileTree } from '../globe/lod.js';
import { createGlobeRenderer } from '../globe/renderer-three.js';
import { TerrainWorkers } from '../globe/workers.js';

const DEG = Math.PI / 180;
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const unit = (v) => { const n = Math.hypot(v[0], v[1], v[2]) || 1; return [v[0] / n, v[1] / n, v[2] / n]; };
const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const smooth = (a, b, x) => { const t = clamp((x - a) / (b - a), 0, 1); return t * t * (3 - 2 * t); };

const BACKDROP_VERTEX = /* glsl */ `
uniform vec3 camF;
uniform vec3 camR;
uniform vec3 camU;
uniform float tanY;
uniform float aspect;
varying vec3 vRay;
void main() {
  vRay = camF + position.x * tanY * aspect * camR + position.y * tanY * camU;
  gl_Position = vec4(position.xy, 1.0, 1.0);
}
`;
const BACKDROP_FRAGMENT = /* glsl */ `
uniform vec3 zenith;
uniform vec3 horizon;
uniform vec3 below;
uniform vec3 space;
uniform vec3 sunTint;
uniform vec3 sunDir;
uniform vec3 upDir;
uniform float air;
uniform float sunUp;
uniform float grain;
varying vec3 vRay;
float noise(vec3 p) { return fract(sin(dot(p, vec3(12.9898, 78.233, 37.719))) * 43758.5453); }
void main() {
  vec3 d = normalize(vRay);
  float e = dot(d, upDir);
  vec3 sky = mix(horizon, zenith, pow(clamp(e, 0.0, 1.0), 0.46));
  sky = mix(sky, below, smoothstep(0.0, -0.16, e));
  float c = max(dot(d, sunDir), 0.0);
  sky += sunTint * sunUp * (pow(c, 1400.0) * 1.6 + pow(c, 60.0) * 0.14 + pow(c, 7.0) * 0.10);
  // The sky's glow about the limb of the Sun's side, from far off a faint wash in the space behind the globe.
  vec3 colour = mix(space, sky, air);
  // Fine points of light in space: a suggestion of stars, far thinner than a catalogue's.
  vec3 q = floor(d * 1500.0);
  float star = step(0.99955, noise(q));
  float tone = noise(q + 3.1);
  colour += grain * (1.0 - air) * star * tone * tone * vec3(0.95, 0.86, 1.0) * 0.55;
  gl_FragColor = vec4(colour, 1.0);
  #include <colorspace_fragment>
}
`;

export function createEarth(THREE, container, options = {}) {
  const opts = { eventElement: container, base: 'eox-s2cloudless-2016', night: 'gibs-blackmarble', maxLevel: 16, fov: 40, preserveDrawingBuffer: false, ...options };
  const fig = figure('earth');
  const canvas = document.createElement('canvas');
  canvas.className = 'earth-canvas';
  container.prepend(canvas);
  const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, logarithmicDepthBuffer: true, powerPreference: 'high-performance', preserveDrawingBuffer: !!opts.preserveDrawingBuffer });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(opts.fov, 1, 1, 2e10);

  const backdrop = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), new THREE.ShaderMaterial({
    vertexShader: BACKDROP_VERTEX, fragmentShader: BACKDROP_FRAGMENT, depthTest: false, depthWrite: false,
    uniforms: Object.fromEntries(['camF', 'camR', 'camU', 'sunDir', 'upDir'].map((k) => [k, { value: new THREE.Vector3() }])
      .concat(['zenith', 'horizon', 'below', 'space', 'sunTint'].map((k) => [k, { value: new THREE.Color() }]))
      .concat([['tanY', { value: 1 }], ['aspect', { value: 1 }], ['air', { value: 0 }], ['sunUp', { value: 1 }], ['grain', { value: 1 }]])),
  }));
  backdrop.frustumCulled = false;
  backdrop.renderOrder = -100;
  scene.add(backdrop);

  const gpu = createGlobeRenderer(THREE, { fig, lighting: 'earth', atmosphere: true });
  scene.add(gpu.group);
  const workers = new TerrainWorkers();
  const layers = () => [opts.base, opts.night].filter((id) => id && IMAGERY[id]).map((id) => imageryLayer(id, IMAGERY[id]));
  const tree = new TileTree({
    fig, maxLevel: opts.maxLevel, heightRange: [-11000, 8900], layers: layers(),
    loadTerrain: (z, x, y, signal) => workers.build({ url: opts.terrainUrl ? opts.terrainUrl(z, x, y) : null, figure: 'earth', z, x, y, exaggeration: 1, surface: 'visible' }, signal),
    on: {
      ready: (t) => gpu.ready(t), released: (t) => gpu.released(t),
      imageReady: (t, layer, image) => gpu.imageReady(t, layer, image), imageReleased: (t, layer) => gpu.imageReleased(t, layer),
    },
    tolerance: 2.5, maxLoads: 12, maxImageLoads: 24, maxTiles: 900,
  });

  const size = { w: 1, h: 1 };
  const rayAt = (cam, clientX, clientY) => {
    const rect = canvas.getBoundingClientRect();
    const x = ((clientX - rect.left) / rect.width) * 2 - 1, y = -(((clientY - rect.top) / rect.height) * 2 - 1);
    const t = Math.tan((camera.fov * DEG) / 2);
    const f = cam.forward, u = cam.up;
    const r = [f[1] * u[2] - f[2] * u[1], f[2] * u[0] - f[0] * u[2], f[0] * u[1] - f[1] * u[0]];
    return { origin: cam.position, direction: unit([0, 1, 2].map((k) => f[k] + x * t * camera.aspect * r[k] + y * t * u[k])) };
  };
  // A city's frame, when one is bound: the level ground beneath it.
  let city = null;
  const heightAt = (lon, lat) => {
    if (city && city.level) {
      const p = sub(cartesian(fig, lon, lat, city.height), city.centre);
      if (Math.abs(dot(p, city.east)) < city.halfEast && Math.abs(dot(p, city.north)) < city.halfNorth) return city.height;
    }
    return tree.heightAt(lon, lat) || 0;
  };
  let changed = () => {};
  const controls = new GlobeControls(opts.eventElement, {
    fig, heightAt, rayAt: (x, y) => rayAt(controls.camera(), x, y),
    onChange: () => changed(controls.view),
  });
  controls.maxRange = 9 * fig.a;
  controls.view = { lon: 150 * DEG, lat: -18 * DEG, height: 0, range: 3.1 * fig.a, heading: 0, tilt: 0 };

  function resize() {
    const w = container.clientWidth || 800, h = container.clientHeight || 600;
    size.w = w, size.h = h;
    renderer.setSize(w, h, false);
    canvas.style.width = '100%';
    canvas.style.height = '100%';
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }
  const observer = new ResizeObserver(resize);
  observer.observe(container);
  resize();

  // ---- what is drawn in the body-fixed frame ---------------------------------------------------------------

  const sets = new Map();     // id → { prisms, centre, born, scale }
  let drape = null;
  let look = { grade: null, haze: null, air: null, sky: null, space: null, tint: null, tintAmount: 0, grain: 1 };
  let sunBody = [1, 0, 0];

  // Columns over rings of [lon, lat] (degrees): { id, cells: [{ ring }], heights, colours, fade: [near, far] (camera ranges, m, at which it is seen) }.
  function addColumns(id, cells, { opacity = 1, shrink = 0.88, fade = null, grow = true } = {}) {
    removeColumns(id);
    if (!cells.length) return null;
    // The set's centre: the mean of its cells, so that single precision holds their positions over a region.
    const c = [0, 0, 0];
    const rings = cells.map((cell) => cell.ring.map(([lon, lat]) => cartesian(fig, lon * DEG, lat * DEG, 0)));
    let n = 0;
    for (const r of rings) for (const p of r) { c[0] += p[0], c[1] += p[1], c[2] += p[2]; n++; }
    const centre = c.map((x) => x / n);
    const prisms = createPrisms(THREE, cells.map((cell, k) => {
      const m = cell.ring.reduce((a, q) => [a[0] + q[0], a[1] + q[1]], [0, 0]).map((x) => x / cell.ring.length);
      return { ring: rings[k].map((p) => sub(p, centre)), up: normal(m[0] * DEG, m[1] * DEG) };
    }), { opacity, shrink });
    prisms.object.matrixAutoUpdate = false;
    scene.add(prisms.object);
    const entry = { prisms, centre, born: grow ? performance.now() : -1e9, fade, opacity, cells };
    sets.set(id, entry);
    return entry;
  }
  function removeColumns(id) {
    const s = sets.get(id);
    if (!s) return;
    scene.remove(s.prisms.object);
    s.prisms.dispose();
    sets.delete(id);
  }

  const rampTexture = (ramp) => {
    const data = new Uint8Array(4 * 256);
    for (let i = 0; i < 256; i++) {
      const c = ramp(i / 255) || [0, 0, 0, 0];
      data[4 * i] = Math.round(c[0] * 255), data[4 * i + 1] = Math.round(c[1] * 255), data[4 * i + 2] = Math.round(c[2] * 255), data[4 * i + 3] = Math.round((c[3] ?? 1) * 255);
    }
    const t = new THREE.DataTexture(data, 256, 1, THREE.RGBAFormat);
    t.colorSpace = THREE.SRGBColorSpace;
    t.minFilter = t.magFilter = THREE.LinearFilter;
    t.needsUpdate = true;
    return t;
  };

  // A raster on the ground: image the map of values (its red channel, read as it is, north up), box [south, west, north, east] (degrees),
  // ramp(u) → [r, g, b, a?] for u in [0, 1]; nearest: classes, not magnitudes.
  // Or values as bytes { data (Uint8Array), width, height }, rows from the south.
  function setDrape(spec) {
    if (drape) { drape.map.dispose(); drape.ramp.dispose(); drape = null; }
    if (!spec) return;
    const map = spec.data ? new THREE.DataTexture(spec.data, spec.width, spec.height, THREE.RedFormat, THREE.UnsignedByteType) : new THREE.Texture(spec.image);
    if (spec.data) map.unpackAlignment = 1;
    map.flipY = false;   // the bitmap comes flipped (createImageBitmap's imageOrientation): its first row is the south
    map.colorSpace = THREE.NoColorSpace;
    map.generateMipmaps = false;
    map.minFilter = map.magFilter = spec.nearest ? THREE.NearestFilter : THREE.LinearFilter;
    map.wrapS = map.wrapT = THREE.ClampToEdgeWrapping;
    map.needsUpdate = true;
    const [s, w, n, e] = spec.box;
    drape = { map, ramp: rampTexture(spec.ramp), box: [w, s, e - w, n - s], opacity: spec.opacity ?? 0.85, cut: spec.cut ?? 1, floor: spec.floor ?? 0 };
  }

  // ---- a city's frame --------------------------------------------------------------------------------------

  // With `level` false the ground stays as it is, and the city is expected to follow it.
  function bindCity({ lat, lon, height = 0, box = null, feather = 1800, level = true }) {
    const [E, N, U] = localFrame(lon * DEG, lat * DEG);
    const origin = cartesian(fig, lon * DEG, lat * DEG, height);
    const [x0, y0, x1, y1] = box || [-3000, -3000, 3000, 3000];
    const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
    city = {
      lat, lon, height, origin, east: E, north: N, up: U, level, cx, cy,
      centre: [0, 1, 2].map((k) => origin[k] + cx * E[k] + cy * N[k]),
      halfEast: (x1 - x0) / 2 + 150, halfNorth: (y1 - y0) / 2 + 150, feather,
    };
    return city;
  }
  // The local east, north, up of a body-fixed point in the city's frame.
  const toLocal = (p) => { const d = sub(p, city.origin); return [dot(d, city.east), dot(d, city.north), dot(d, city.up)]; };
  const fromLocal = (x, y, z = 0) => [0, 1, 2].map((k) => city.origin[k] + x * city.east[k] + y * city.north[k] + z * city.up[k]);

  // A view of the city: its point of regard (local east, north), range, tilt and heading (radians).
  function localView(x, y, range, tilt, heading) {
    const g = geodetic(fig, ...fromLocal(x, y, 0));
    return { lon: g.lon, lat: g.lat, height: city.height, range, tilt, heading };
  }

  // The driver a city scene is given: its camera follows this view, and its fit, flights and the plug-ins' own
  // moves of the camera turn into moves of the view.
  function driver() {
    return {
      fit(b, { tilt = 52, heading = 0, range: far, instant = false, seconds = 3 } = {}) {
        const cx = (b.x0 + b.x1) / 2, cy = (b.y0 + b.y1) / 2, s = Math.max(b.x1 - b.x0, b.y1 - b.y0, 200);
        const v = localView(cx, cy, s * (far ?? 0.62) * 1.25, tilt * DEG, heading * DEG);
        if (instant) controls.setView(v); else controls.flyTo(v, seconds);
      },
      flyTo([x, y], range = null, tilt = null) {
        const g = geodetic(fig, ...fromLocal(x, y, 0));
        const v = { lon: g.lon, lat: g.lat, height: city.height };
        if (range) v.range = range;
        if (tilt !== null) v.tilt = tilt * DEG;
        controls.flyTo(v, 1.1);
      },
      // A plug-in moved the camera and its target (three.js's local axes: x east, y up, z south).
      place(position, target) {
        const d = [position.x - target.x, position.y - target.y, position.z - target.z];
        const range = Math.hypot(...d) || 1;
        controls.setView(localView(target.x, -target.z, range, Math.acos(clamp(d[1] / range, -1, 1)), Math.atan2(-d[0], d[2])));
      },
      sync(cam3, targetLocal) {
        const c = controls.camera();
        const [x, y, z] = toLocal(c.position);
        const f = [dot(c.forward, city.east), dot(c.forward, city.north), dot(c.forward, city.up)];
        const u = [dot(c.up, city.east), dot(c.up, city.north), dot(c.up, city.up)];
        cam3.position.set(x, z, -y);
        cam3.up.set(u[0], u[2], -u[1]);
        cam3.lookAt(x + f[0], z + f[2], -y - f[1]);
        const near = Math.max(0.5, Math.min(2000, Math.abs(z) * 0.08));
        const far = Math.max(60000, controls.view.range * 12);
        if (Math.abs(cam3.near - near) / near > 0.15 || cam3.far !== far || cam3.fov !== camera.fov || cam3.aspect !== camera.aspect) {
          cam3.near = near, cam3.far = far, cam3.fov = camera.fov, cam3.aspect = camera.aspect;
          cam3.updateProjectionMatrix();
        }
        const t = toLocal(cartesian(fig, controls.view.lon, controls.view.lat, controls.view.height));
        targetLocal.set(t[0], t[2], -t[1]);
      },
    };
  }

  // ---- a frame ---------------------------------------------------------------------------------------------

  const frustum = new THREE.Frustum(), projView = new THREE.Matrix4(), viewRot = new THREE.Matrix3();
  const I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  const sv = new THREE.Vector3();

  function frame(now) {
    controls.tick(now);
    const c = controls.camera();
    const P = c.position;
    camera.position.set(0, 0, 0);
    camera.up.set(...c.up);
    camera.lookAt(...c.forward);
    const g = geodetic(fig, ...P);
    const altitude = Math.max(1, g.h - (tree.heightAt(g.lon, g.lat) || 0));
    const near = Math.max(0.5, Math.min(5e4, altitude * 0.02));
    if (Math.abs(camera.near - near) / near > 0.2) { camera.near = near; camera.updateProjectionMatrix(); }
    camera.updateMatrixWorld();
    projView.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse);
    frustum.setFromProjectionMatrix(projView);
    const rel = [-P[0], -P[1], -P[2]];
    const planes = frustum.planes.map((p) => [p.normal.x, p.normal.y, p.normal.z, p.constant + p.normal.x * rel[0] + p.normal.y * rel[1] + p.normal.z * rel[2]]);
    const selection = tree.update({ camera: P, planes, height: size.h, fovY: camera.fov * DEG });
    viewRot.setFromMatrix4(camera.matrixWorldInverse);
    sv.set(...sunBody).applyMatrix3(viewRot);
    // Near the ground the air thickens: haze over distance, and the sky behind.
    const air = 1 - smooth(25000, 900000, altitude);
    gpu.draw(selection, {
      // With nightAll, the whole visible side is the night side: the Sun stands behind the Earth from wherever one looks.
      bodyRelative: rel, rotation: I, sunBody: look.nightAll ? unit(P).map((v) => -v) : sunBody, base: opts.base, night: opts.night, overlay: 0, exposure: 1,
      sunView: [sv.x, sv.y, sv.z], atmosphere: true, grade: look.grade, air: look.air,
      // Haze over distances long against the camera's height: the ground beneath stays clear, the horizon goes.
      haze: look.haze ? { colour: look.haze, density: (air * (look.hazeScale ?? 1)) / (30 * altitude + 15000) } : null,
      flat: city && city.level ? { centre: city.centre, east: city.east, north: city.north, halfEast: city.halfEast, halfNorth: city.halfNorth, feather: city.feather, height: city.height, amount: city.amount ?? 1, tint: look.tint, tintAmount: look.tintAmount } : null,
      drape,
    });
    for (const s of sets.values()) {
      const r = sub(s.centre, P);
      s.prisms.object.matrix.set(1, 0, 0, r[0], 0, 1, 0, r[1], 0, 0, 1, r[2], 0, 0, 0, 1);
      s.prisms.object.matrixWorldNeedsUpdate = true;
      s.prisms.setSun(sv.x, sv.y, sv.z);
      const range = controls.view.range;
      const seen = s.fade ? smooth(s.fade[0] * 0.6, s.fade[0], range) * (1 - smooth(s.fade[1], s.fade[1] * 1.7, range)) : 1;
      s.prisms.setOpacity(seen * s.opacity);
      const t = (now - s.born) / 1100;
      s.prisms.setGrow(t >= 1 ? 1 : 1 - (1 - Math.max(0, t)) ** 3);
    }
    // The backdrop.
    const b = backdrop.material.uniforms;
    const forward = c.forward, up = c.up;
    const right = [forward[1] * up[2] - forward[2] * up[1], forward[2] * up[0] - forward[0] * up[2], forward[0] * up[1] - forward[1] * up[0]];
    b.camF.value.set(...forward), b.camR.value.set(...right), b.camU.value.set(...up);
    b.tanY.value = Math.tan((camera.fov * DEG) / 2), b.aspect.value = camera.aspect;
    const n = normal(g.lon, g.lat);
    b.upDir.value.set(...n), b.sunDir.value.set(...sunBody);
    const sunUp = dot(n, sunBody);
    b.air.value = air, b.sunUp.value = smooth(-0.12, 0.08, sunUp);
    if (look.sky) {
      const sky = look.sky(sunUp);
      b.zenith.value.setRGB(...sky.zenith), b.horizon.value.setRGB(...sky.horizon), b.below.value.setRGB(...sky.below), b.sunTint.value.setRGB(...sky.sun);
      if (look.haze) look.haze = sky.horizon;
    }
    if (look.space) b.space.value.setRGB(...look.space);
    b.grain.value = look.grain ?? 1;
    renderer.render(scene, camera);
  }

  // Where a point of the ground stands on the screen, or null behind the globe. For place marks.
  function project(lonDeg, latDeg, h = 0) {
    const c = controls.camera();
    const p = cartesian(fig, lonDeg * DEG, latDeg * DEG, h);
    const d = sub(p, c.position);
    const n = normal(lonDeg * DEG, latDeg * DEG);
    if (dot(d, n) > 0) return null;   // beyond the horizon
    const f = c.forward, u = c.up;
    const r = [f[1] * u[2] - f[2] * u[1], f[2] * u[0] - f[0] * u[2], f[0] * u[1] - f[1] * u[0]];
    const z = dot(d, f);
    if (z <= 0) return null;
    const t = Math.tan((camera.fov * DEG) / 2);
    return { x: (0.5 + (0.5 * dot(d, r)) / (z * t * camera.aspect)) * size.w, y: (0.5 - (0.5 * dot(d, u)) / (z * t)) * size.h, distance: Math.hypot(...d) };
  }

  // The place of the ground under a pointer, or null off the globe: { lon, lat } in degrees.
  function groundAt(clientX, clientY) {
    const { origin, direction } = rayAt(controls.camera(), clientX, clientY);
    // Over the real ground: the height at the point of regard is near enough for a pointer over a tile.
    const hit = rayHit(fig, origin, direction, city ? (city.level ? city.height : controls.view.height) : 0);
    if (!hit) return null;
    const g = geodetic(fig, ...hit);
    return { lon: g.lon / DEG, lat: g.lat / DEG };
  }

  return {
    THREE, fig, canvas, renderer, scene, camera, controls, tree, frame, resize, project, groundAt, bindCity, driver,
    addColumns, removeColumns, columns: sets, setDrape, setDrapeCut: (cut) => { if (drape) drape.cut = cut; }, hasDrape: () => !!drape, setDrapeOpacity: (a) => { if (drape) drape.opacity = a; },
    onChange: (f) => { changed = f; },
    setLook: (changes) => { look = { ...look, ...changes }; },
    // The Sun, from its subsolar point (degrees).
    setSubsolar: (lat, lon) => { sunBody = normal(lon * DEG, lat * DEG); },
    get city() { return city; },
    get view() { return controls.view; },
    groundMetresPerPixel() {
      // At the point of regard.
      const d = Math.hypot(...sub(controls.camera().position, cartesian(fig, controls.view.lon, controls.view.lat, controls.view.height)));
      return (2 * d * Math.tan((camera.fov * DEG) / 2)) / size.h;
    },
    setView: (v) => controls.setView(v),
    flyTo: (v, seconds) => controls.flyTo(v, seconds),
    dispose() {
      observer.disconnect();
      workers.terminate();
      gpu.dispose();
      renderer.dispose();
      canvas.remove();
    },
  };
}
