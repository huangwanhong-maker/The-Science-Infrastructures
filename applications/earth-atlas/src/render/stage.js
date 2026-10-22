// The stage: the WebGL canvas, the camera for the globe or the flat map, and
// the gestures that move them.
//
// World coordinates are Earth-fixed for the globe (unit sphere; x toward 0 E,
// y toward 90 E, z north) and the plane of the projection for a flat map
// (units of the Earth's radius). The stage turns places into world positions
// for the layers (a projector), and screen positions into places for the probe.

import * as THREE from 'three';
import { DEG, toVector, wrapLongitude } from '../core/sphere.js';
import { PROJECTIONS, extent, forward, inverse, outline } from '../core/projections.js';
import { createSurfaceMaterial } from './surface.js';

const FOV = 28;

function globeGeometry(nlon = 256, nlat = 128) {
  const pos = new Float32Array((nlon + 1) * (nlat + 1) * 3);
  let k = 0;
  for (let j = 0; j <= nlat; j++) {
    const lat = 90 - (180 * j) / nlat;
    for (let i = 0; i <= nlon; i++) {
      const v = toVector(-180 + (360 * i) / nlon, lat);
      pos[k++] = v[0]; pos[k++] = v[1]; pos[k++] = v[2];
    }
  }
  const index = [];
  for (let j = 0; j < nlat; j++) {
    for (let i = 0; i < nlon; i++) {
      const a = j * (nlon + 1) + i, b = a + nlon + 1;
      index.push(a, b, a + 1, a + 1, b, b + 1);
    }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setIndex(index);
  return g;
}

function haloMaterial() {
  return new THREE.ShaderMaterial({
    glslVersion: THREE.GLSL3,
    transparent: true, depthWrite: false, side: THREE.BackSide, blending: THREE.AdditiveBlending,
    uniforms: { uColor: { value: new THREE.Color('#d9a6cc') } },
    vertexShader: /* glsl */ `
      out vec3 vNormal;
      out vec3 vView;
      void main() {
        vec4 mv = modelViewMatrix * vec4(position, 1.0);
        vNormal = normalize(normalMatrix * normal);
        vView = normalize(-mv.xyz);
        gl_Position = projectionMatrix * mv;
      }`,
    fragmentShader: /* glsl */ `
      in vec3 vNormal;
      in vec3 vView;
      uniform vec3 uColor;
      out vec4 fragColor;
      void main() {
        float rim = clamp(1.0 + dot(vNormal, vView) * 1.35, 0.0, 1.0);
        fragColor = vec4(uColor, pow(rim, 3.0) * 0.3);
      }`,
  });
}

export class Stage {
  constructor(canvas) {
    this.canvas = canvas;
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.setClearColor(0x000000, 0);
    this.scene = new THREE.Scene();
    this.material = createSurfaceMaterial();
    this.globe = new THREE.Mesh(globeGeometry(), this.material);
    this.halo = new THREE.Mesh(new THREE.SphereGeometry(1.06, 96, 48), haloMaterial());
    this.halo.renderOrder = -1;
    this.plane = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), this.material);
    this.scene.add(this.globe, this.halo, this.plane);
    this.overlay = new THREE.Group();
    this.scene.add(this.overlay);

    this.perspective = new THREE.PerspectiveCamera(FOV, 1, 0.01, 100);
    this.perspective.up.set(0, 0, 1);
    this.orthographic = new THREE.OrthographicCamera(-1, 1, 1, -1, -10, 10);
    this.resolution = new THREE.Vector2(1, 1);

    this.mode = 'globe';
    this.lambda0 = 0;
    this.globeView = { lon: 135, lat: 22, dist: 4.6 };
    this.mapView = { x: 0, y: 0, zoom: 1 };
    this.velocity = { lon: 0, lat: 0 };
    this.listeners = { change: [], hover: [], click: [], projection: [] };
    this.revision = 0;
    this.resize();
    this.setProjection('globe');
    this.bindGestures();
    new ResizeObserver(() => this.resize()).observe(canvas.parentElement);
  }

  on(event, fn) { this.listeners[event].push(fn); }
  emit(event, ...args) { for (const fn of this.listeners[event]) fn(...args); }

  get camera() { return this.mode === 'globe' ? this.perspective : this.orthographic; }
  get projection() { return PROJECTIONS[this.mode] || null; }

  resize() {
    const el = this.canvas.parentElement;
    const w = Math.max(1, el.clientWidth), h = Math.max(1, el.clientHeight);
    this.width = w; this.height = h;
    this.renderer.setSize(w, h, false);
    this.renderer.getDrawingBufferSize(this.resolution);
    this.perspective.aspect = w / h;
    this.perspective.updateProjectionMatrix();
    this.updateCamera();
  }

  // ---- projection ---------------------------------------------------------

  setProjection(mode) {
    const wasGlobe = this.mode === 'globe';
    this.mode = mode;
    const u = this.material.uniforms;
    this.globe.visible = this.halo.visible = mode === 'globe';
    this.plane.visible = mode !== 'globe';
    if (mode === 'globe') {
      u.uProjection.value = 0;
    } else {
      const p = this.projection;
      const [xm, ym] = extent(p);
      // The flat map is centred on the meridian the globe was facing.
      if (wasGlobe) this.lambda0 = Math.round(this.globeView.lon / 5) * 5;
      u.uProjection.value = p.id === 'equirectangular' ? 1 : 2;
      u.uLambda0.value = this.lambda0 * DEG;
      u.uExtent.value.set(xm, ym);
      this.plane.scale.set(xm * 1.001, ym * 1.001, 1);
      this.mapView = { x: 0, y: 0, zoom: 1 };
    }
    this.updateCamera();
    this.emit('projection', mode);
  }

  // A projector for the layers: world position of a place, and whether a
  // segment between two world positions crosses the map's edge.
  projector() {
    if (this.mode === 'globe') {
      return {
        globe: true,
        project: (lon, lat, out, r = 1.0015) => { toVector(lon, lat, out); out[0] *= r; out[1] *= r; out[2] *= r; return true; },
        breaks: () => false,
      };
    }
    const p = this.projection, l0 = this.lambda0, [xm] = extent(p);
    return {
      globe: false,
      project: (lon, lat, out) => {
        const xy = forward(p, lon, lat, l0);
        out[0] = xy[0]; out[1] = xy[1]; out[2] = 0.002;
        return true;
      },
      breaks: (a, b) => Math.abs(a[0] - b[0]) > xm * 0.5,
    };
  }

  mapOutline() {
    return this.mode === 'globe' ? null : outline(this.projection);
  }

  // ---- camera -------------------------------------------------------------

  updateCamera() {
    if (this.mode === 'globe') {
      const v = this.globeView;
      v.lat = Math.max(-85, Math.min(85, v.lat));
      v.lon = wrapLongitude(v.lon);
      v.dist = Math.max(1.04, Math.min(9, v.dist));
      const p = toVector(v.lon, v.lat);
      this.perspective.position.set(p[0] * v.dist, p[1] * v.dist, p[2] * v.dist);
      this.perspective.lookAt(0, 0, 0);
      this.perspective.near = Math.max(0.001, (v.dist - 1) * 0.5);
      this.perspective.far = v.dist + 2;
      this.perspective.updateProjectionMatrix();
    } else {
      const [xm, ym] = extent(this.projection);
      const m = this.mapView;
      m.zoom = Math.max(0.8, Math.min(40, m.zoom));
      // Zoom 1 shows the whole map with a margin.
      const fit = Math.max(ym * 1.06, (xm * 1.03 * this.height) / this.width);
      const halfH = fit / m.zoom, halfW = (halfH * this.width) / this.height;
      m.x = Math.max(-xm, Math.min(xm, m.x));
      m.y = Math.max(-ym, Math.min(ym, m.y));
      Object.assign(this.orthographic, { left: m.x - halfW, right: m.x + halfW, top: m.y + halfH, bottom: m.y - halfH });
      this.orthographic.position.set(0, 0, 5);
      this.orthographic.updateProjectionMatrix();
    }
    // At once, not at the next render: names and arrows are placed on the screen
    // (screen()) earlier in the frame than the scene is drawn.
    this.camera.updateMatrixWorld();
    this.revision++;
    this.emit('change');
  }

  // Pixels per radian of the Earth's surface at the centre of the view.
  pixelsPerRadian() {
    if (this.mode === 'globe') {
      return this.height / 2 / Math.tan((FOV / 2) * DEG) / (this.globeView.dist - 1);
    }
    return this.height / (this.orthographic.top - this.orthographic.bottom);
  }

  pixelsPerMetre() { return this.pixelsPerRadian() / 6371008.8; }

  // The zoom level of a web map (256-pixel tiles) with the same scale, for
  // choosing which names to show.
  zoomLevel() { return Math.log2((2 * Math.PI * this.pixelsPerRadian()) / 256); }

  // ---- places and the screen ------------------------------------------------

  // Screen position (CSS pixels) of a place, with its visibility: fades out
  // toward the limb of the globe. Returns [x, y, visibility] or null.
  screen(lon, lat) {
    const v = new THREE.Vector3();
    let vis = 1;
    if (this.mode === 'globe') {
      const p = toVector(lon, lat);
      v.set(p[0], p[1], p[2]);
      const c = this.perspective.position;
      const facing = (v.x * c.x + v.y * c.y + v.z * c.z) / c.length();       // cos of angle from the view axis
      const limb = 1 / this.globeView.dist;
      vis = Math.max(0, Math.min(1, (facing - limb) / 0.12));
      if (vis <= 0) return null;
    } else {
      const xy = forward(this.projection, lon, lat, this.lambda0);
      v.set(xy[0], xy[1], 0);
    }
    v.project(this.camera);
    return [((v.x + 1) / 2) * this.width, ((1 - v.y) / 2) * this.height, vis];
  }

  // The place under a screen position (CSS pixels), or null.
  place(x, y) {
    const ndc = new THREE.Vector2((x / this.width) * 2 - 1, 1 - (y / this.height) * 2);
    if (this.mode === 'globe') {
      const ray = new THREE.Raycaster();
      ray.setFromCamera(ndc, this.perspective);
      const o = ray.ray.origin, d = ray.ray.direction;
      const b = o.dot(d), c = o.lengthSq() - 1, disc = b * b - c;
      if (disc < 0) return null;
      const t = -b - Math.sqrt(disc);
      const p = o.clone().addScaledVector(d, t);
      return [Math.atan2(p.y, p.x) / DEG, Math.asin(Math.max(-1, Math.min(1, p.z))) / DEG];
    }
    const w = new THREE.Vector3(ndc.x, ndc.y, 0).unproject(this.orthographic);
    return inverse(this.projection, w.x, w.y, this.lambda0);
  }

  // A sampler of random places in view, uniform in area, for particles.
  spawner() {
    if (this.mode === 'globe') {
      const v = this.globeView;
      const cap = Math.acos(1 / v.dist) * 0.98;
      const centre = toVector(v.lon, v.lat);
      const east = [-Math.sin(v.lon * DEG), Math.cos(v.lon * DEG), 0];
      const north = [centre[1] * east[2] - centre[2] * east[1], centre[2] * east[0] - centre[0] * east[2], centre[0] * east[1] - centre[1] * east[0]];
      const zmin = Math.cos(cap);
      return (random) => {
        const z = zmin + (1 - zmin) * random(), r = Math.sqrt(1 - z * z), a = 2 * Math.PI * random();
        const x = r * Math.cos(a), y = r * Math.sin(a);
        const p = [0, 1, 2].map((k) => z * centre[k] + x * east[k] + y * north[k]);
        return [Math.atan2(p[1], p[0]) / DEG, Math.asin(Math.max(-1, Math.min(1, p[2]))) / DEG];
      };
    }
    const o = this.orthographic, p = this.projection, l0 = this.lambda0;
    return (random) => {
      for (let i = 0; i < 12; i++) {
        const x = o.left + (o.right - o.left) * random(), y = o.bottom + (o.top - o.bottom) * random();
        const q = inverse(p, x, y, l0);
        if (q) return q;
      }
      return [l0, 0];
    };
  }

  // Fraction of the screen covered by the Earth, for scaling particle numbers.
  coverage() {
    if (this.mode === 'globe') {
      const r = this.pixelsPerRadian() * Math.tan(Math.asin(1 / this.globeView.dist)) * (this.globeView.dist - 1) / 1;
      return Math.min(1, (Math.PI * r * r) / (this.width * this.height));
    }
    return 0.75;
  }

  // ---- gestures -------------------------------------------------------------

  bindGestures() {
    const el = this.canvas;
    const pointers = new Map();
    let last = null, moved = 0, pinch = 0, lastTime = 0;
    el.addEventListener('pointerdown', (e) => {
      el.setPointerCapture(e.pointerId);
      pointers.set(e.pointerId, [e.clientX, e.clientY]);
      last = [e.clientX, e.clientY]; moved = 0; lastTime = performance.now();
      this.velocity.lon = this.velocity.lat = 0;
      this.dragging = true;
      if (pointers.size === 2) pinch = this.pinchDistance(pointers);
    });
    el.addEventListener('pointermove', (e) => {
      const r = el.getBoundingClientRect();
      if (!pointers.has(e.pointerId)) {
        this.emit('hover', this.place(e.clientX - r.left, e.clientY - r.top), e);
        return;
      }
      pointers.set(e.pointerId, [e.clientX, e.clientY]);
      if (pointers.size === 2) {
        const d = this.pinchDistance(pointers);
        if (pinch > 0) this.zoomBy(pinch / d, null);
        pinch = d;
        return;
      }
      const dx = e.clientX - last[0], dy = e.clientY - last[1];
      moved += Math.abs(dx) + Math.abs(dy);
      last = [e.clientX, e.clientY];
      const now = performance.now(), dt = Math.max(1, now - lastTime);
      lastTime = now;
      this.pan(dx, dy, dt);
    });
    const end = (e) => {
      const r = el.getBoundingClientRect();
      if (pointers.has(e.pointerId) && pointers.size === 1 && moved < 4) {
        this.emit('click', this.place(e.clientX - r.left, e.clientY - r.top), e);
      }
      pointers.delete(e.pointerId);
      this.dragging = pointers.size > 0;
      if (pointers.size < 2) pinch = 0;
      if (performance.now() - lastTime > 80) this.velocity.lon = this.velocity.lat = 0;
    };
    el.addEventListener('pointerup', end);
    el.addEventListener('pointercancel', end);
    el.addEventListener('pointerleave', () => this.emit('hover', null));
    el.addEventListener('wheel', (e) => {
      e.preventDefault();
      const r = el.getBoundingClientRect();
      this.zoomBy(Math.exp(e.deltaY * (e.deltaMode === 1 ? 0.05 : 0.0015)), [e.clientX - r.left, e.clientY - r.top]);
    }, { passive: false });
    el.addEventListener('dblclick', (e) => {
      const r = el.getBoundingClientRect();
      const place = this.place(e.clientX - r.left, e.clientY - r.top);
      if (place && this.mode === 'globe') { this.globeView.lon = place[0]; this.globeView.lat = place[1]; }
      this.zoomBy(0.55, [e.clientX - r.left, e.clientY - r.top]);
    });
  }

  pinchDistance(pointers) {
    const [a, b] = [...pointers.values()];
    return Math.hypot(a[0] - b[0], a[1] - b[1]);
  }

  pan(dx, dy, dt = 16) {
    if (this.mode === 'globe') {
      const k = 1 / this.pixelsPerRadian() / DEG;
      const v = this.globeView;
      const dlon = (-dx * k) / Math.max(0.2, Math.cos(v.lat * DEG)), dlat = dy * k;
      v.lon += dlon; v.lat += dlat;
      this.velocity.lon = (dlon / dt) * 16; this.velocity.lat = (dlat / dt) * 16;
    } else {
      const k = (this.orthographic.top - this.orthographic.bottom) / this.height;
      this.mapView.x -= dx * k; this.mapView.y += dy * k;
    }
    this.updateCamera();
  }

  // Zooms by a factor (< 1 zooms in) about a screen point.
  zoomBy(factor, at) {
    if (this.mode === 'globe') {
      this.globeView.dist = 1 + (this.globeView.dist - 1) * factor;
    } else {
      const o = this.orthographic;
      const before = at && new THREE.Vector3((at[0] / this.width) * 2 - 1, 1 - (at[1] / this.height) * 2, 0).unproject(o);
      this.mapView.zoom /= factor;
      this.updateCamera();
      if (before) {
        const after = new THREE.Vector3((at[0] / this.width) * 2 - 1, 1 - (at[1] / this.height) * 2, 0).unproject(o);
        this.mapView.x += before.x - after.x; this.mapView.y += before.y - after.y;
      }
    }
    this.updateCamera();
  }

  // Inertia after a drag of the globe; returns true while moving.
  coast() {
    if (this.mode !== 'globe' || this.dragging) return false;
    const v = this.velocity;
    if (Math.abs(v.lon) + Math.abs(v.lat) < 0.002) { v.lon = v.lat = 0; return false; }
    this.globeView.lon += v.lon; this.globeView.lat += v.lat;
    v.lon *= 0.92; v.lat *= 0.92;
    this.updateCamera();
    return true;
  }

  render() {
    this.renderer.render(this.scene, this.camera);
  }
}
