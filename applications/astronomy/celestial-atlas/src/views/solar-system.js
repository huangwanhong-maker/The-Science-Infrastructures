// The Solar System in three dimensions, heliocentric, in the ecliptic frame of J2000.
//
// Distances are true. Body sizes may be exaggerated so that the planets stay
// visible (the factor is shown to the user); the Moon's distance from the Earth
// is then exaggerated by a smaller factor so that it does not sit inside the
// enlarged Earth. At a size factor of 1 everything is to scale.
//
// Bodies carry maps of their surfaces (or cloud tops), turned by the IAU
// rotation models (src/core/rotation.js) to their orientation at the instant
// shown; the Earth shows its lights on the night side; the sky behind the stars
// is the Milky Way.
//
// The Earth is turned instead by Greenwich mean sidereal time and precession,
// the frame in which the Sky Chart places an observer (src/core/sky.js), so that
// a place chosen by clicking the Earth is the place whose sky the chart shows.
// (The IAU model of the Earth differs from it by 0.05 to 0.26 degrees in
// longitude between 2000 and 2100.)

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';
import { BODIES, BY_KEY } from '../core/bodies.js';
import { AU_KM, DEG, ICRS_TO_ECLIPTIC, ICRS_TO_GALACTIC, apply, multiply, transpose, sub, rotZ } from '../core/frames.js';
import { geocentricLatitude, geodeticLatitude } from '../core/geodesy.js';
import { bodyFixedMatrix } from '../core/rotation.js';
import { BODY_FROM_LOCAL, quaternionOf, skyRotation, skySphere, nightLights } from '../render/maps.js';
import { BodyDetail } from '../render/detail.js';
import { osculatingOrbit } from '../core/sky.js';
import { starPoints, magnitudeStyle } from '../render/points.js';
import { starRgb } from '../render/colour.js';
import { declutter } from '../render/declutter.js';

const SKY_RADIUS = 5e4;          // AU; radius of the star sphere around the camera
const SUN_SIZE_CAP = 20;         // the Sun is never enlarged more than this
const STAR_LIMIT = 7.0;          // faintest magnitude on the star sphere
const ORBIT_REFRESH_S = 6 * 3600;
const GRID_RADII = [1, 2, 5, 10, 20, 30, 50];
const GRID_LABELS = new Set([1, 5, 10, 30, 50]);
const LABEL_PRIORITY = ['sun', 'jupiter', 'saturn', 'earth', 'venus', 'mars', 'mercury',
  'uranus', 'neptune', 'moon', 'pluto'];

// Ecliptic (x to the equinox, z to the north ecliptic pole) to scene (y up).
const toScene = (e, out = new THREE.Vector3()) => out.set(e[0], e[2], -e[1]);
const SCENE_FROM_ECLIPTIC = [1, 0, 0, 0, 0, 1, 0, -1, 0];
const SCENE_FROM_ICRS = multiply(SCENE_FROM_ECLIPTIC, ICRS_TO_ECLIPTIC);
const icrsToEclipticAu = (v) => {
  const e = apply(ICRS_TO_ECLIPTIC, v);
  return [e[0] / AU_KM, e[1] / AU_KM, e[2] / AU_KM];
};

export function moonDistanceFactor(sizeScale) {
  return Math.max(1, sizeScale * 0.025);
}

function radialTexture(stops, size = 256) {
  const canvas = document.createElement('canvas');
  canvas.width = canvas.height = size;
  const ctx = canvas.getContext('2d');
  const g = ctx.createRadialGradient(size / 2, size / 2, 0, size / 2, size / 2, size / 2);
  for (const [at, colour] of stops) g.addColorStop(at, colour);
  ctx.fillStyle = g;
  ctx.fillRect(0, 0, size, size);
  const texture = new THREE.CanvasTexture(canvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

function sceneLabel(className, html, onClick) {
  const el = document.createElement(onClick ? 'button' : 'div');
  el.className = className;
  el.innerHTML = html;
  if (onClick) {
    el.type = 'button';
    el.addEventListener('click', (event) => { event.stopPropagation(); onClick(); });
  }
  return { el, object: new CSS2DObject(el) };
}

export class SolarSystemView {
  constructor({ renderer, container, stars, constellations, onSelect }) {
    this.renderer = renderer;
    this.onSelect = onSelect;
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(45, 1, 1e-7, 1e6);
    this.camera.position.set(0, 2.6, 3.4);
    this.controls = new OrbitControls(this.camera, renderer.domElement);
    Object.assign(this.controls, { enableDamping: true, dampingFactor: 0.08, zoomSpeed: 1.3,
      minDistance: 2e-5, maxDistance: 400 });
    this.labels = new CSS2DRenderer();
    this.labels.domElement.className = 'label-layer';
    container.appendChild(this.labels.domElement);

    this.options = { sizeScale: 600, orbits: true, labels: true, stars: true, constellations: false, grid: true, milkyWay: true };
    this.items = {};
    this.focusKey = null;
    this.selectedKey = null;
    this.lastOrbitEt = -Infinity;
    this.tween = null;
    this.lastFocus = new THREE.Vector3();
    this.raycaster = new THREE.Raycaster();

    this.buildLights();
    this.buildBodies();
    this.buildGrid();
    this.buildSky(stars, constellations);
  }

  buildLights() {
    this.scene.add(new THREE.AmbientLight(0xf6eaf4, 0.07));
    const sunlight = new THREE.PointLight(0xfff1dc, 2.8, 0, 0);
    this.scene.add(sunlight);
  }

  buildBodies() {
    // Fine enough that the limb stays round when a body fills the view.
    const sphere = new THREE.SphereGeometry(1, 192, 96);
    for (const b of BODIES) {
      const group = new THREE.Group();
      const material = b.key === 'sun'
        ? new THREE.MeshBasicMaterial({ color: b.colour })
        : new THREE.MeshStandardMaterial({ color: b.colour, roughness: 0.92, metalness: 0 });
      const mesh = new THREE.Mesh(sphere, material);
      mesh.userData.key = b.key;
      group.add(mesh);
      if (b.rings) mesh.add(this.makeRings(b));
      if (b.key === 'sun') {
        const glow = new THREE.Sprite(new THREE.SpriteMaterial({
          map: radialTexture([[0, 'rgba(255,236,196,0.95)'], [0.18, 'rgba(255,200,130,0.45)'],
            [0.45, 'rgba(240,150,110,0.12)'], [1, 'rgba(240,150,110,0)']]),
          blending: THREE.AdditiveBlending, depthWrite: false, transparent: true,
        }));
        group.add(glow);
        this.glow = glow;
      }
      // The Moon's label sits to the left of its dot, away from the Earth's.
      const left = b.key === 'moon';
      const label = sceneLabel(`body-label${left ? ' left' : ''}`,
        `<span class="dot" style="background:${b.colour}"></span><span class="name">${b.name}</span>`,
        () => this.onSelect(b.key, { focus: true }));
      label.object.center.set(left ? 1 : 0, 0.5);
      group.add(label.object);
      this.scene.add(group);

      let orbit = null;
      if (b.key !== 'sun') {
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(256 * 3), 3));
        orbit = new THREE.LineLoop(geometry, new THREE.LineBasicMaterial({
          color: b.colour, transparent: true, opacity: b.parent ? 0.5 : 0.34 }));
        orbit.frustumCulled = false;
        this.scene.add(orbit);   // the Moon's orbit is moved to the Earth each frame
      }
      this.items[b.key] = { body: b, group, mesh, label, orbit };
    }
  }

  makeRings(b) {
    const inner = b.rings.inner / b.radius, outer = b.rings.outer / b.radius;
    this.ringRadii = [inner, outer];
    // Radial bands of Saturn's C, B and A rings and the Cassini division (km).
    const f = (km) => km / b.rings.outer;
    const texture = radialTexture([
      [0, 'rgba(0,0,0,0)'], [f(74500), 'rgba(210,190,160,0.10)'], [f(92000), 'rgba(222,203,168,0.35)'],
      [f(95000), 'rgba(236,218,180,0.85)'], [f(117000), 'rgba(230,212,176,0.80)'],
      [f(117580), 'rgba(40,30,30,0.05)'], [f(122170), 'rgba(40,30,30,0.05)'],
      [f(122300), 'rgba(222,203,168,0.6)'], [f(136000), 'rgba(210,192,160,0.5)'],
      [1, 'rgba(210,192,160,0)'],
    ], 512);
    const geometry = new THREE.RingGeometry(inner, outer, 192, 1);
    // Texture coordinates by radius, for a radial map of the rings.
    const pos = geometry.attributes.position, uv = geometry.attributes.uv;
    for (let i = 0; i < pos.count; i++) uv.setXY(i, (Math.hypot(pos.getX(i), pos.getY(i)) - inner) / (outer - inner), 0.5);
    const ring = new THREE.Mesh(geometry,
      new THREE.MeshBasicMaterial({ map: texture, transparent: true, side: THREE.DoubleSide, depthWrite: false }));
    // In the planet's equatorial plane (the body's local x-z plane).
    ring.rotation.x = -Math.PI / 2;
    this.ring = ring;
    return ring;
  }

  buildGrid() {
    this.grid = new THREE.Group();
    const material = new THREE.LineBasicMaterial({ color: 0xf1e2ee, transparent: true, opacity: 0.075 });
    for (const r of GRID_RADII) {
      const points = [];
      for (let k = 0; k <= 192; k++) {
        const a = (2 * Math.PI * k) / 192;
        points.push(new THREE.Vector3(r * Math.cos(a), 0, -r * Math.sin(a)));
      }
      this.grid.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material));
      if (GRID_LABELS.has(r)) {
        const label = sceneLabel('grid-label', `${r} AU`);
        label.object.position.set(r * Math.cos(-0.35), 0, -r * Math.sin(-0.35));
        this.grid.add(label.object);
      }
    }
    const equinox = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(0, 0, 0), new THREE.Vector3(60, 0, 0)]),
      new THREE.LineBasicMaterial({ color: 0xca8dbe, transparent: true, opacity: 0.28 }));
    this.grid.add(equinox);
    const mark = sceneLabel('grid-label accent', '♈︎ equinox');
    mark.object.position.set(1.15, 0, 0);
    this.grid.add(mark.object);
    this.scene.add(this.grid);
  }

  buildSky(stars, constellations) {
    this.sky = new THREE.Group();
    const chosen = [];
    for (const i of stars.byBrightness) {
      if (stars.mag[i] > STAR_LIMIT) break;
      chosen.push(i);
    }
    const positions = new Float32Array(chosen.length * 3);
    const colours = new Float32Array(chosen.length * 3);
    const sizes = new Float32Array(chosen.length);
    const v = new THREE.Vector3();
    chosen.forEach((i, k) => {
      const e = apply(ICRS_TO_ECLIPTIC, [stars.unit[3 * i], stars.unit[3 * i + 1], stars.unit[3 * i + 2]]);
      toScene(e, v).multiplyScalar(SKY_RADIUS).toArray(positions, 3 * k);
      const { size, intensity } = magnitudeStyle(stars.mag[i], STAR_LIMIT);
      const [r, g, b] = starRgb(stars.ci[i]);
      colours.set([r * intensity, g * intensity, b * intensity], 3 * k);
      sizes[k] = size;
    });
    this.starPoints = starPoints(positions, colours, sizes, this.renderer.getPixelRatio());
    this.sky.add(this.starPoints);

    const segs = constellations.segments;
    const linePositions = new Float32Array(segs.length);
    for (let k = 0; k < segs.length; k += 3) {
      const e = apply(ICRS_TO_ECLIPTIC, [segs[k], segs[k + 1], segs[k + 2]]);
      toScene(e, v).multiplyScalar(SKY_RADIUS * 0.995).toArray(linePositions, k);
    }
    const lineGeometry = new THREE.BufferGeometry();
    lineGeometry.setAttribute('position', new THREE.BufferAttribute(linePositions, 3));
    this.constellationLines = new THREE.LineSegments(lineGeometry,
      new THREE.LineBasicMaterial({ color: 0xd9b3d0, transparent: true, opacity: 0.16, depthWrite: false }));
    this.constellationLines.frustumCulled = false;
    this.sky.add(this.constellationLines);
    this.scene.add(this.sky);
  }

  // Maps from src/render/maps.js: surfaces, Saturn's rings, the Earth by day and night, the sky.
  setMaps(maps) {
    if (!maps) return;
    for (const [key, item] of Object.entries(this.items)) {
      const texture = key === 'earth' ? maps.earthDay || maps.bodies.earth : maps.bodies[key];
      if (!texture) continue;
      const m = item.mesh.material;
      m.map = texture;
      m.color.set(key === 'sun' ? '#fff3e0' : '#ffffff');
      if (key === 'earth' && maps.earthNight) {
        m.emissiveMap = maps.earthNight;
        m.emissive = new THREE.Color(1.0, 0.82, 0.58);
        m.emissiveIntensity = 1.25;
        nightLights(m);
        this.earthMaterial = m;
      }
      m.needsUpdate = true;
    }
    // Finer maps of the Moon, Mars and Mercury, fetched when one fills the view.
    this.detail = new BodyDetail(this.renderer);
    for (const [key, live] of Object.entries(maps.live || {})) {
      if (this.items[key] && this.items[key].mesh.material.map) this.detail.attach(key, this.items[key].mesh, live);
    }
    if (maps.rings && this.ring) {
      this.ring.material.map = maps.rings;
      this.ring.material.needsUpdate = true;
    }
    if (maps.sky) {
      this.milkyWay = skySphere(maps.sky, SKY_RADIUS * 1.02, 0.8);
      this.milkyWay.quaternion.copy(skyRotation(multiply(SCENE_FROM_ICRS, transpose(ICRS_TO_GALACTIC))));
      this.sky.add(this.milkyWay);
    }
    this.mapped = true;
  }

  setOptions(options) {
    const sizeChanged = options.sizeScale !== undefined && options.sizeScale !== this.options.sizeScale;
    Object.assign(this.options, options);
    if (sizeChanged) this.lastOrbitEt = -Infinity;
  }

  displayRadius(key) {
    const b = BY_KEY[key];
    const factor = key === 'sun' ? Math.min(this.options.sizeScale, SUN_SIZE_CAP) : this.options.sizeScale;
    return (b.radius / AU_KM) * factor;
  }

  bodyPosition(snap, key, out = new THREE.Vector3()) {
    if (key === 'moon') {
      const earth = this.bodyPosition(snap, 'earth', new THREE.Vector3());
      const rel = sub(snap.bodies.moon.ssb.p, snap.bodies.earth.ssb.p);
      const k = moonDistanceFactor(this.options.sizeScale);
      return toScene(icrsToEclipticAu(rel), out).multiplyScalar(k).add(earth);
    }
    return toScene(icrsToEclipticAu(snap.bodies[key].helio.p), out);
  }

  updateOrbits(snap) {
    const k = moonDistanceFactor(this.options.sizeScale);
    const v = new THREE.Vector3();
    for (const b of BODIES) {
      const item = this.items[b.key];
      if (!item.orbit) continue;
      const points = osculatingOrbit(snap, b.key, 256);
      if (!points) continue;
      const attr = item.orbit.geometry.attributes.position;
      const factor = b.parent ? k : 1;
      for (let j = 0; j < 256; j++) {
        const e = icrsToEclipticAu([points[3 * j], points[3 * j + 1], points[3 * j + 2]]);
        toScene(e, v).multiplyScalar(factor).toArray(attr.array, 3 * j);
      }
      attr.needsUpdate = true;
    }
    this.lastOrbitEt = snap.et;
  }

  focus(key) {
    this.focusKey = key;
    const target = this.controls.target.clone();
    const distance = this.camera.position.distanceTo(target);
    const desired = key === 'sun' ? 4.5
      : Math.max(this.displayRadius(key) * 8, key === 'moon' ? 0.004 : 0.02);
    this.tween = { from: target, distFrom: distance, distTo: desired, t: 0 };
  }

  select(key) {
    this.selectedKey = key;
    for (const [k, item] of Object.entries(this.items)) {
      item.label.el.classList.toggle('selected', k === key);
    }
  }

  pick(ndc) {
    this.raycaster.setFromCamera(ndc, this.camera);
    const hits = this.raycaster.intersectObjects(Object.values(this.items).map((i) => i.mesh), false);
    return hits.length ? hits[0].object.userData.key : null;
  }

  // The place on the Earth under the pointer, { lat, lon } in radians (geodetic,
  // WGS 84), or null when the pointer is not on the Earth.
  pickPlace(ndc) {
    this.raycaster.setFromCamera(ndc, this.camera);
    const mesh = this.items.earth.mesh;
    const hit = this.raycaster.intersectObject(mesh, false)[0];
    if (!hit) return null;
    const p = mesh.worldToLocal(hit.point.clone()).normalize();
    // Body-fixed axes from the sphere's local ones: (x, y, z)_body = (x, -z, y)_local.
    const lon = Math.atan2(-p.z, p.x), psi = Math.asin(Math.max(-1, Math.min(1, p.y)));
    return { lat: geodeticLatitude(psi), lon };
  }

  // A mark on the Earth at the observer's place (radians).
  setObserver(place) {
    if (!this.observerMark) {
      const canvas = document.createElement('canvas');
      canvas.width = canvas.height = 64;
      const g = canvas.getContext('2d');
      g.strokeStyle = '#CA8DBE'; g.lineWidth = 7;
      g.beginPath(); g.arc(32, 32, 22, 0, 2 * Math.PI); g.stroke();
      g.fillStyle = '#fdf4fa';
      g.beginPath(); g.arc(32, 32, 8, 0, 2 * Math.PI); g.fill();
      const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(canvas), sizeAttenuation: false, depthWrite: false }));
      sprite.scale.set(0.022, 0.022, 1);
      sprite.renderOrder = 5;
      this.items.earth.group.add(sprite);
      this.observerMark = sprite;
    }
    const psi = geocentricLatitude(place.lat);
    // Local direction of the place: (x, y, z)_local = (x, z, -y)_body.
    this.observerLocal = new THREE.Vector3(Math.cos(psi) * Math.cos(place.lon), Math.sin(psi), -Math.cos(psi) * Math.sin(place.lon));
  }

  update(snap, dt) {
    const pos = new THREE.Vector3();
    const jdTdb = 2451545.0 + snap.et / 86400;
    for (const b of BODIES) {
      const item = this.items[b.key];
      this.bodyPosition(snap, b.key, pos);
      item.group.position.copy(pos);
      item.mesh.scale.setScalar(this.displayRadius(b.key));
      // The body's orientation: local -> body-fixed -> ICRF -> scene.
      const B = b.key === 'earth' ? multiply(rotZ(snap.gmst), snap.precession) : bodyFixedMatrix(b.key, jdTdb);
      if (B) quaternionOf(multiply(SCENE_FROM_ICRS, multiply(transpose(B), BODY_FROM_LOCAL)), item.mesh.quaternion);
      if (item.orbit) {
        item.orbit.visible = this.options.orbits;
        if (b.parent) item.orbit.position.copy(this.items[b.parent].group.position);
      }
    }
    this.glow.scale.setScalar(this.displayRadius('sun') * 9);
    if (this.observerMark) {
      const earth = this.items.earth.mesh;
      this.observerMark.position.copy(this.observerLocal).applyQuaternion(earth.quaternion).multiplyScalar(this.displayRadius('earth') * 1.004);
    }
    if (Math.abs(snap.et - this.lastOrbitEt) > ORBIT_REFRESH_S) this.updateOrbits(snap);

    this.grid.visible = this.options.grid;
    this.starPoints.visible = this.options.stars;
    if (this.milkyWay) this.milkyWay.visible = this.options.stars && this.options.milkyWay;
    this.constellationLines.visible = this.options.stars && this.options.constellations;

    if (this.focusKey) {
      const target = this.items[this.focusKey].group.position;
      if (this.tween) {
        this.tween.t = Math.min(1, this.tween.t + dt / 0.9);
        const k = this.tween.t < 0.5 ? 4 * this.tween.t ** 3 : 1 - (-2 * this.tween.t + 2) ** 3 / 2;
        const direction = this.camera.position.clone().sub(this.controls.target).normalize();
        this.controls.target.lerpVectors(this.tween.from, target, k);
        const distance = this.tween.distFrom + (this.tween.distTo - this.tween.distFrom) * k;
        this.camera.position.copy(this.controls.target).addScaledVector(direction, distance);
        if (this.tween.t >= 1) this.tween = null;
      } else {
        this.camera.position.add(pos.copy(target).sub(this.lastFocus));
        this.controls.target.copy(target);
      }
      this.lastFocus.copy(target);
    }
    // Close to a focused body, down to a few hundredths of its radius above it.
    this.controls.minDistance = this.focusKey ? this.displayRadius(this.focusKey) * 1.03 : 2e-5;
    this.controls.update();
    this.sky.position.copy(this.camera.position);

    this.scene.updateMatrixWorld();
    const selected = this.selectedKey;
    const order = selected ? [selected, ...LABEL_PRIORITY.filter((k) => k !== selected)] : LABEL_PRIORITY;
    declutter(order.map((k) => ({
      object: this.items[k].label.object,
      text: BY_KEY[k].name,
      wanted: this.options.labels || k === selected,
      always: k === selected,
      left: k === 'moon',
    })), this.camera, this.width, this.height);
  }

  resize(width, height) {
    this.width = width;
    this.height = height;
    this.camera.aspect = width / height;
    this.camera.updateProjectionMatrix();
    this.labels.setSize(width, height);
  }

  setActive(active) {
    this.controls.enabled = active;
    this.labels.domElement.style.display = active ? '' : 'none';
  }

  render() {
    if (this.earthMaterial) {
      // The direction of the Sun from the Earth, in view space, for the lights at night.
      const earth = this.items.earth.group.position, sun = this.items.sun.group.position;
      this.earthMaterial.userData.sunView.value.copy(sun).sub(earth).normalize().transformDirection(this.camera.matrixWorldInverse);
    }
    if (this.detail) this.detail.update(this.camera, this.width, this.height, this.focusKey);
    this.renderer.render(this.scene, this.camera);
    this.labels.render(this.scene, this.camera);
  }

  legend() {
    const s = this.options.sizeScale;
    if (s <= 1) return 'All sizes and distances to scale';
    const moon = moonDistanceFactor(s);
    return `Bodies enlarged ×${Math.round(s)} (Sun ×${Math.min(s, SUN_SIZE_CAP)})`
      + (moon > 1 ? `, Moon's distance ×${moon.toFixed(moon < 10 ? 1 : 0)}` : '');
  }
}
