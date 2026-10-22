// The stage: one WebGL canvas drawing up to four views of the atlas: the three
// anatomical slices through the current point, and the brain in the round.
//
// The template and the label index are 3D textures; a slice is a quad through
// the point, whose fragments sample them by their atlas coordinates (RAS mm),
// so that the same material draws a plane of any orientation and the planes
// of the three-dimensional view. The views follow the neurological convention:
// the subject's left is on the viewer's left, the sagittal view is from the left.
//
// Surfaces (the brain's outline, a region, the scalp), data points and
// electrodes are meshes in the same scene; the planes live on layers 1-3, the
// rest on layer 4, and each view's camera enables what it shows.

import * as THREE from 'three';
import { invert, worldBounds } from '../../../../domains/neuroscience/js/index.js';

export const VIEWS = ['axial', 'coronal', 'sagittal', 'three'];
const AXIS = { axial: 2, coronal: 1, sagittal: 0 };
const LAYER = { axial: 1, coronal: 2, sagittal: 3, three: 4 };

function sliceMaterial(textures) {
  return new THREE.ShaderMaterial({
    glslVersion: THREE.GLSL3,
    side: THREE.DoubleSide,
    uniforms: {
      uVolume: { value: textures.volume }, uLabels: { value: textures.labels }, uPalette: { value: textures.palette },
      uAtlasToVoxel: { value: new THREE.Matrix4() }, uDims: { value: new THREE.Vector3(1, 1, 1) },
      uLow: { value: 0 }, uHigh: { value: 1 }, uOverlay: { value: 0.55 }, uSelected: { value: -1 },
      uDiscardEmpty: { value: false }, uHasLabels: { value: false }, uTint: { value: new THREE.Color('#fff2ea') },
    },
    vertexShader: /* glsl */ `
      out vec3 vAtlas;
      void main() {
        vAtlas = position;
        gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
      }`,
    fragmentShader: /* glsl */ `
      precision highp float;
      precision highp usampler3D;
      precision highp sampler3D;
      in vec3 vAtlas;
      uniform sampler3D uVolume;
      uniform usampler3D uLabels;
      uniform sampler2D uPalette;
      uniform mat4 uAtlasToVoxel;
      uniform vec3 uDims;
      uniform float uLow, uHigh, uOverlay;
      uniform int uSelected;
      uniform bool uDiscardEmpty, uHasLabels;
      uniform vec3 uTint;
      out vec4 fragColor;
      void main() {
        vec3 v = (uAtlasToVoxel * vec4(vAtlas, 1.0)).xyz;
        vec3 tc = (v + 0.5) / uDims;
        if (any(lessThan(tc, vec3(0.0))) || any(greaterThan(tc, vec3(1.0)))) discard;
        float g = texture(uVolume, tc).r;
        float grey = clamp((g - uLow) / max(uHigh - uLow, 1e-4), 0.0, 1.0);
        grey = pow(grey, 0.9);
        uint label = 0u;
        if (uHasLabels) {
          ivec3 iv = clamp(ivec3(floor(v + 0.5)), ivec3(0), ivec3(uDims) - 1);
          label = texelFetch(uLabels, iv, 0).r;
        }
        if (uDiscardEmpty && grey <= 0.01 && label == 0u) discard;
        vec3 colour = grey * uTint;
        if (label != 0u) {
          int l = int(label);
          vec4 pc = texelFetch(uPalette, ivec2(l & 255, l >> 8), 0);
          colour = mix(colour, pc.rgb * (0.55 + 0.6 * grey), pc.a * uOverlay);
          if (l == uSelected) colour = mix(colour, vec3(1.0, 0.95, 0.86), 0.45);
        }
        fragColor = vec4(colour, 1.0);
      }`,
  });
}

const SURFACE_COLOURS = { brain: '#d9a6cc', scalp: '#f4dbc8', region: '#f6c47e' };

export class Stage {
  constructor(canvas, viewElements) {
    this.canvas = canvas;
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.setClearColor(0x000000, 0);
    this.renderer.autoClear = false;
    this.scene = new THREE.Scene();
    this.group = new THREE.Group();            // the atlas, scaled to a head when asked
    this.scene.add(this.group);
    this.views = viewElements;                 // {axial: element, ...}
    this.layout = 'four';
    this.point = [0, 0, 0];
    this.bounds = { min: [-100, -100, -100], max: [100, 100, 100] };
    this.listeners = { point: [], hover: [], render: [] };
    this.dirty = true;
    this.pickables = [];
    this.textures = this.makeTextures();
    this.planes = {};
    for (const view of ['axial', 'coronal', 'sagittal']) {
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(12), 3));
      geometry.setIndex([0, 1, 2, 0, 2, 3]);
      const mesh = new THREE.Mesh(geometry, sliceMaterial(this.textures));
      mesh.layers.set(LAYER[view]);
      mesh.userData.view = view;
      mesh.frustumCulled = false;
      this.group.add(mesh);
      this.planes[view] = mesh;
      this.pickables.push(mesh);
    }
    this.surfaces = {};
    this.markers = new THREE.Group();
    this.markers.layers.set(LAYER.three);
    this.group.add(this.markers);
    this.lights();
    this.cameras = {
      axial: new THREE.OrthographicCamera(-1, 1, 1, -1, -2000, 2000),
      coronal: new THREE.OrthographicCamera(-1, 1, 1, -1, -2000, 2000),
      sagittal: new THREE.OrthographicCamera(-1, 1, 1, -1, -2000, 2000),
      three: new THREE.PerspectiveCamera(30, 1, 1, 5000),
    };
    for (const v of ['axial', 'coronal', 'sagittal']) this.cameras[v].layers.set(LAYER[v]);
    this.cameras.three.layers.enable(LAYER.axial); this.cameras.three.layers.enable(LAYER.coronal); this.cameras.three.layers.enable(LAYER.sagittal); this.cameras.three.layers.enable(LAYER.three);
    this.cameras.three.layers.enable(0);
    this.orbit = { azimuth: 35, elevation: 22, distance: 520, target: [0, -15, 10] };
    this.zoom2d = { axial: 1, coronal: 1, sagittal: 1 };   // mm per pixel is set by fit
    this.pan2d = { axial: [0, 0], coronal: [0, 0], sagittal: [0, 0] };
    this.bindGestures();
    new ResizeObserver(() => this.resize()).observe(canvas.parentElement);
    this.resize();
    this.animate = this.animate.bind(this);
    requestAnimationFrame(this.animate);
  }

  on(event, fn) { this.listeners[event].push(fn); }
  emit(event, ...args) { for (const fn of this.listeners[event]) fn(...args); }
  invalidate() { this.dirty = true; }

  makeTextures() {
    const volume = new THREE.Data3DTexture(new Uint8Array(8), 2, 2, 2);
    volume.format = THREE.RedFormat; volume.type = THREE.UnsignedByteType;
    volume.minFilter = THREE.LinearFilter; volume.magFilter = THREE.LinearFilter; volume.unpackAlignment = 1; volume.needsUpdate = true;
    const labels = new THREE.Data3DTexture(new Uint16Array(8), 2, 2, 2);
    labels.format = THREE.RedIntegerFormat; labels.type = THREE.UnsignedShortType; labels.internalFormat = 'R16UI';
    labels.minFilter = THREE.NearestFilter; labels.magFilter = THREE.NearestFilter; labels.unpackAlignment = 1; labels.needsUpdate = true;
    const palette = new THREE.DataTexture(new Uint8Array(256 * 256 * 4), 256, 256, THREE.RGBAFormat, THREE.UnsignedByteType);
    palette.minFilter = THREE.NearestFilter; palette.magFilter = THREE.NearestFilter; palette.needsUpdate = true;
    return { volume, labels, palette };
  }

  lights() {
    const key = new THREE.DirectionalLight(0xfff1e6, 1.6);
    key.position.set(0.6, -0.5, 1.2);
    const fill = new THREE.DirectionalLight(0xd9a6cc, 0.5);
    fill.position.set(-0.8, 0.3, -0.4);
    const ambient = new THREE.AmbientLight(0xfff6f0, 0.55);
    for (const l of [key, fill, ambient]) { l.layers.enable(LAYER.three); this.scene.add(l); }
  }

  // ---- the data ----------------------------------------------------------------------------

  // A 3D texture is made anew when its size changes (a texture once uploaded
  // keeps its dimensions), and every plane's material takes the new one.
  replaceTexture(name, data, dims) {
    const old = this.textures[name];
    const t = new THREE.Data3DTexture(data, dims[0], dims[1], dims[2]);
    t.format = old.format; t.type = old.type; if (old.internalFormat) t.internalFormat = old.internalFormat;
    t.minFilter = old.minFilter; t.magFilter = old.magFilter; t.unpackAlignment = 1; t.needsUpdate = true;
    this.textures[name] = t;
    const uniform = { volume: 'uVolume', labels: 'uLabels' }[name];
    for (const p of Object.values(this.planes)) p.material.uniforms[uniform].value = t;
    old.dispose();
    return t;
  }

  // The template: 8-bit intensities on its grid, with the affine; the bounds set the views.
  setVolume(display, dims, affine) {
    this.replaceTexture('volume', display, dims);
    this.dims = dims.slice();
    this.affine = affine.slice();
    const inv = invert(affine);
    const m = new THREE.Matrix4().set(...inv);
    for (const plane of Object.values(this.planes)) {
      plane.material.uniforms.uAtlasToVoxel.value.copy(m);
      plane.material.uniforms.uDims.value.set(dims[0], dims[1], dims[2]);
    }
    this.bounds = worldBounds({ dims, affine, data: display });
    const c = [0, 1, 2].map((a) => (this.bounds.min[a] + this.bounds.max[a]) / 2);
    const size = Math.max(...[0, 1, 2].map((a) => this.bounds.max[a] - this.bounds.min[a]));
    this.orbit.target = c;
    this.orbit.distance = size * 2.6;
    this.cameras.three.near = size * 0.05; this.cameras.three.far = size * 20;
    this.cameras.three.updateProjectionMatrix();
    for (const v of ['axial', 'coronal', 'sagittal']) { this.pan2d[v] = [0, 0]; this.zoom2d[v] = 1; }
    this.fit();
    this.setPoint(c);
  }

  setWindow(low, high) { for (const p of Object.values(this.planes)) { p.material.uniforms.uLow.value = low; p.material.uniforms.uHigh.value = high; } this.invalidate(); }
  setOverlay(alpha) { for (const p of Object.values(this.planes)) p.material.uniforms.uOverlay.value = alpha; this.invalidate(); }
  setSelected(index) { for (const p of Object.values(this.planes)) p.material.uniforms.uSelected.value = index ?? -1; this.invalidate(); }

  // The label index on the template's grid (null for none).
  setLabels(index) {
    if (index) this.replaceTexture('labels', index, this.dims);
    for (const p of Object.values(this.planes)) p.material.uniforms.uHasLabels.value = Boolean(index);
    this.invalidate();
  }

  // The palette: rgba bytes for 65536 indices (index 0 unused).
  setPalette(rgba) {
    this.textures.palette.image.data.set(rgba.subarray(0, 256 * 256 * 4));
    this.textures.palette.needsUpdate = true;
    this.invalidate();
  }

  // A surface mesh (positions, normals, indices in atlas mm), by kind; null removes it.
  setSurface(kind, mesh, { colour = SURFACE_COLOURS[kind], opacity = 1, colours = null } = {}) {
    if (this.surfaces[kind]) { this.group.remove(this.surfaces[kind]); this.surfaces[kind].geometry.dispose(); this.surfaces[kind].material.dispose(); delete this.surfaces[kind]; }
    this.pickables = this.pickables.filter((p) => !p.userData.surface);
    if (!mesh) { this.invalidate(); return; }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(mesh.positions, 3));
    geometry.setAttribute('normal', new THREE.BufferAttribute(mesh.normals, 3));
    if (colours) geometry.setAttribute('color', new THREE.BufferAttribute(colours, 3));
    geometry.setIndex(new THREE.BufferAttribute(mesh.indices, 1));
    const material = new THREE.MeshPhongMaterial({ color: new THREE.Color(colour), transparent: opacity < 1, opacity, depthWrite: opacity >= 1, shininess: 18, specular: new THREE.Color('#2a1d2c'), vertexColors: Boolean(colours), side: THREE.FrontSide });
    const object = new THREE.Mesh(geometry, material);
    object.layers.set(LAYER.three);
    object.renderOrder = opacity < 1 ? 2 : 1;
    object.userData.surface = kind;
    this.group.add(object);
    this.surfaces[kind] = object;
    if (opacity >= 1 || kind === 'brain') this.pickables.push(object);
    this.invalidate();
  }

  surfaceVisible(kind, visible) { if (this.surfaces[kind]) { this.surfaces[kind].visible = visible; this.invalidate(); } }

  // Colours per vertex (Float32Array, 3 per vertex, 0..1) on a surface, or none; a
  // coloured surface is drawn opaque.
  setSurfaceColours(kind, colours, { opacity = null } = {}) {
    const mesh = this.surfaces[kind];
    if (!mesh) return;
    const g = mesh.geometry, m = mesh.material;
    if (colours) { g.setAttribute('color', new THREE.BufferAttribute(colours, 3)); g.attributes.color.needsUpdate = true; }
    else if (g.attributes.color) g.deleteAttribute('color');
    m.vertexColors = Boolean(colours);
    const alpha = opacity ?? (colours ? 1 : m.userData.baseOpacity ?? m.opacity);
    if (m.userData.baseOpacity === undefined) m.userData.baseOpacity = m.opacity;
    m.opacity = alpha; m.transparent = alpha < 1; m.depthWrite = alpha >= 1;
    m.needsUpdate = true;
    this.invalidate();
  }

  // Spheres at atlas positions: [{position, colour, radius, name, kind}].
  setMarkers(markers) {
    while (this.markers.children.length) { const m = this.markers.children.pop(); m.geometry.dispose(); m.material.dispose(); }
    this.markerList = markers;
    for (const m of markers) {
      const sphere = new THREE.Mesh(new THREE.SphereGeometry(m.radius || 2.5, 18, 12), new THREE.MeshPhongMaterial({ color: new THREE.Color(m.colour || '#e8a0a0'), shininess: 30, emissive: new THREE.Color(m.colour || '#e8a0a0'), emissiveIntensity: 0.25 }));
      sphere.position.set(...m.position);
      sphere.layers.set(LAYER.three);
      sphere.userData.marker = m;
      this.markers.add(sphere);
    }
    this.invalidate();
  }

  // The atlas scaled by (sx, sy, sz) about a centre, in atlas mm.
  setScale(scales, centre) {
    const [sx, sy, sz] = scales;
    const c = centre || [0, 0, 0];
    this.group.scale.set(sx, sy, sz);
    this.group.position.set(c[0] * (1 - sx), c[1] * (1 - sy), c[2] * (1 - sz));
    this.group.updateMatrixWorld(true);
    this.invalidate();
  }

  // ---- the point ---------------------------------------------------------------------------

  setPoint(p, { from = null } = {}) {
    this.point = p.map((x) => (Number.isFinite(x) ? x : 0));
    this.updatePlanes();
    this.invalidate();
    this.emit('point', this.point, from);
  }

  updatePlanes() {
    const b = this.bounds, p = this.point, pad = 2;
    const lo = [b.min[0] - pad, b.min[1] - pad, b.min[2] - pad], hi = [b.max[0] + pad, b.max[1] + pad, b.max[2] + pad];
    const set = (view, corners) => {
      const a = this.planes[view].geometry.attributes.position;
      corners.forEach((c, i) => a.setXYZ(i, c[0], c[1], c[2]));
      a.needsUpdate = true;
      this.planes[view].geometry.computeBoundingSphere();
    };
    set('axial', [[lo[0], lo[1], p[2]], [hi[0], lo[1], p[2]], [hi[0], hi[1], p[2]], [lo[0], hi[1], p[2]]]);
    set('coronal', [[lo[0], p[1], lo[2]], [hi[0], p[1], lo[2]], [hi[0], p[1], hi[2]], [lo[0], p[1], hi[2]]]);
    set('sagittal', [[p[0], lo[1], lo[2]], [p[0], hi[1], lo[2]], [p[0], hi[1], hi[2]], [p[0], lo[1], hi[2]]]);
  }

  // ---- views and cameras -------------------------------------------------------------------

  setLayout(layout) { this.layout = layout; this.resize(); }

  resize() {
    const el = this.canvas.parentElement;
    const w = Math.max(1, el.clientWidth), h = Math.max(1, el.clientHeight);
    this.width = w; this.height = h;
    this.renderer.setSize(w, h, false);
    this.fit();
    this.invalidate();
  }

  // The rectangle of a view in canvas pixels (y up, as WebGL has it), or null when hidden.
  rect(view) {
    const el = this.views[view];
    if (!el || el.offsetParent === null) return null;
    const r = el.getBoundingClientRect(), c = this.canvas.getBoundingClientRect();
    if (r.width < 2 || r.height < 2) return null;
    return { x: r.left - c.left, y: c.bottom - r.bottom, w: r.width, h: r.height, top: r.top - c.top };
  }

  // Each orthographic camera framed on the volume, keeping its zoom and pan.
  fit() {
    const b = this.bounds;
    const extent = [0, 1, 2].map((a) => b.max[a] - b.min[a]);
    for (const view of ['axial', 'coronal', 'sagittal']) {
      const r = this.rect(view);
      if (!r) continue;
      const axes = { axial: [0, 1], coronal: [0, 2], sagittal: [1, 2] }[view];
      const mmPerPx = Math.max(extent[axes[0]] / r.w, extent[axes[1]] / r.h) * 1.08 / this.zoom2d[view];
      const cam = this.cameras[view];
      const halfW = (r.w * mmPerPx) / 2, halfH = (r.h * mmPerPx) / 2;
      cam.left = -halfW; cam.right = halfW; cam.top = halfH; cam.bottom = -halfH;
      cam.updateProjectionMatrix();
      this.mmPerPx = this.mmPerPx || {};
      this.mmPerPx[view] = mmPerPx;
    }
  }

  placeCameras() {
    const c = [0, 1, 2].map((a) => (this.bounds.min[a] + this.bounds.max[a]) / 2);
    const far = 1000;
    const place = (view, position, up) => {
      const cam = this.cameras[view];
      const pan = this.pan2d[view];
      const look = c.slice();
      // Pan in the view's screen axes.
      const axes = { axial: [0, 1], coronal: [0, 2], sagittal: [1, 2] }[view];
      look[axes[0]] += pan[0]; look[axes[1]] += pan[1];
      cam.position.set(look[0] + position[0] * far, look[1] + position[1] * far, look[2] + position[2] * far);
      cam.up.set(...up);
      cam.lookAt(look[0], look[1], look[2]);
      cam.updateMatrixWorld();
    };
    place('axial', [0, 0, 1], [0, 1, 0]);        // from above: x right, y up
    place('coronal', [0, -1, 0], [0, 0, 1]);     // from behind: x right, z up (the subject's left on the left)
    place('sagittal', [-1, 0, 0], [0, 0, 1]);    // from the left: anterior to the left
    const o = this.orbit, DEG = Math.PI / 180;
    const cam = this.cameras.three;
    const t = o.target;
    cam.position.set(
      t[0] + o.distance * Math.cos(o.elevation * DEG) * Math.sin(o.azimuth * DEG),
      t[1] - o.distance * Math.cos(o.elevation * DEG) * Math.cos(o.azimuth * DEG),
      t[2] + o.distance * Math.sin(o.elevation * DEG),
    );
    cam.up.set(0, 0, 1);
    cam.lookAt(t[0], t[1], t[2]);
    cam.updateMatrixWorld();
  }

  // Screen position (pixels within the view) of an atlas point, in a view; null when hidden.
  project(view, p) {
    const r = this.rect(view);
    if (!r) return null;
    const v = new THREE.Vector3(p[0], p[1], p[2]).applyMatrix4(this.group.matrixWorld).project(this.cameras[view]);
    return { x: ((v.x + 1) / 2) * r.w, y: ((1 - v.y) / 2) * r.h, depth: v.z, inside: v.z < 1 };
  }

  // The atlas point under screen coordinates (pixels within a 2D view), on the view's slice plane.
  unproject(view, x, y) {
    const r = this.rect(view);
    if (!r || view === 'three') return null;
    const cam = this.cameras[view];
    const v = new THREE.Vector3((x / r.w) * 2 - 1, 1 - (y / r.h) * 2, 0).unproject(cam);
    const inv = this.group.matrixWorld.clone().invert();
    v.applyMatrix4(inv);
    const p = [v.x, v.y, v.z];
    p[AXIS[view]] = this.point[AXIS[view]];
    return p;
  }

  // The atlas point hit by a ray through screen coordinates of the 3D view: a plane, a surface or a marker.
  pick(x, y) {
    const r = this.rect('three');
    if (!r) return null;
    const raycaster = new THREE.Raycaster();
    raycaster.layers.enableAll();
    raycaster.setFromCamera(new THREE.Vector2((x / r.w) * 2 - 1, 1 - (y / r.h) * 2), this.cameras.three);
    const hits = raycaster.intersectObjects([...this.pickables, ...this.markers.children], false);
    const hit = hits.find((h) => h.object.visible && (h.object.userData.view || h.object.userData.surface === 'region' || h.object.userData.marker || h.object.userData.surface === 'brain'));
    if (!hit) return null;
    const local = hit.point.clone().applyMatrix4(this.group.matrixWorld.clone().invert());
    return { point: [local.x, local.y, local.z], marker: hit.object.userData.marker || null, surface: hit.object.userData.surface || null };
  }

  // ---- gestures ----------------------------------------------------------------------------

  bindGestures() {
    for (const view of VIEWS) {
      const el = this.views[view];
      if (!el) continue;
      let drag = null;
      el.addEventListener('pointerdown', (e) => {
        if (e.button !== 0) return;
        el.setPointerCapture(e.pointerId);
        const r = el.getBoundingClientRect();
        drag = { x: e.clientX, y: e.clientY, moved: false, shift: e.shiftKey, start: { ...this.orbit, pan: this.pan2d[view] ? this.pan2d[view].slice() : null } };
        if (view !== 'three') {
          if (e.shiftKey) return;
          const p = this.unproject(view, e.clientX - r.left, e.clientY - r.top);
          if (p) this.setPoint(p, { from: view });
        }
      });
      el.addEventListener('pointermove', (e) => {
        const r = el.getBoundingClientRect();
        const x = e.clientX - r.left, y = e.clientY - r.top;
        if (drag) {
          const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
          if (Math.abs(dx) + Math.abs(dy) > 2) drag.moved = true;
          if (view === 'three') {
            if (drag.shift) {
              const s = this.orbit.distance / 600;
              const cam = this.cameras.three;
              const right = new THREE.Vector3().setFromMatrixColumn(cam.matrixWorld, 0), up = new THREE.Vector3().setFromMatrixColumn(cam.matrixWorld, 1);
              const t = drag.start.target;
              this.orbit.target = [t[0] - dx * s * right.x + dy * s * up.x, t[1] - dx * s * right.y + dy * s * up.y, t[2] - dx * s * right.z + dy * s * up.z];
            } else {
              this.orbit.azimuth = drag.start.azimuth - dx * 0.4;
              this.orbit.elevation = Math.max(-89, Math.min(89, drag.start.elevation + dy * 0.4));
            }
            this.invalidate();
          } else if (drag.shift) {
            const mm = this.mmPerPx ? this.mmPerPx[view] : 1;
            this.pan2d[view] = [drag.start.pan[0] - dx * mm, drag.start.pan[1] + dy * mm];
            this.invalidate();
          } else {
            const p = this.unproject(view, x, y);
            if (p) this.setPoint(p, { from: view });
          }
        } else {
          this.emit('hover', view, x, y);
        }
      });
      const end = (e) => {
        if (!drag) return;
        if (view === 'three' && !drag.moved) {
          const r = el.getBoundingClientRect();
          const hit = this.pick(e.clientX - r.left, e.clientY - r.top);
          if (hit) this.setPoint(hit.point, { from: 'three', marker: hit.marker });
        }
        drag = null;
      };
      el.addEventListener('pointerup', end);
      el.addEventListener('pointercancel', () => { drag = null; });
      el.addEventListener('pointerleave', () => this.emit('hover', view, null, null));
      el.addEventListener('wheel', (e) => {
        e.preventDefault();
        if (view === 'three') {
          this.orbit.distance = Math.max(20, Math.min(5000, this.orbit.distance * Math.exp(e.deltaY * 0.0012)));
        } else if (e.ctrlKey || e.metaKey) {
          this.zoom2d[view] = Math.max(0.5, Math.min(12, this.zoom2d[view] * Math.exp(-e.deltaY * 0.0015)));
          this.fit();
        } else {
          const step = (e.shiftKey ? 5 : 1) * (e.deltaY > 0 ? -1 : 1);
          const p = this.point.slice();
          p[AXIS[view]] = Math.round((p[AXIS[view]] + step) * 100) / 100;
          this.setPoint(p, { from: view });
        }
        this.invalidate();
      }, { passive: false });
      el.addEventListener('dblclick', () => { if (view !== 'three') { this.zoom2d[view] = 1; this.pan2d[view] = [0, 0]; this.fit(); this.invalidate(); } });
    }
  }

  // ---- rendering ---------------------------------------------------------------------------

  animate() {
    requestAnimationFrame(this.animate);
    if (!this.dirty) return;
    this.dirty = false;
    this.render();
  }

  render() {
    const gl = this.renderer;
    gl.setScissorTest(true);
    this.placeCameras();
    this.group.updateMatrixWorld(true);
    const dpr = gl.getPixelRatio();
    gl.setScissor(0, 0, this.width * dpr, this.height * dpr);
    gl.setViewport(0, 0, this.width * dpr, this.height * dpr);
    gl.clear(true, true, true);
    for (const view of VIEWS) {
      const r = this.rect(view);
      if (!r) continue;
      const cam = this.cameras[view];
      if (view === 'three') { cam.aspect = r.w / r.h; cam.updateProjectionMatrix(); }
      gl.setViewport(r.x * dpr, r.y * dpr, r.w * dpr, r.h * dpr);
      gl.setScissor(r.x * dpr, r.y * dpr, r.w * dpr, r.h * dpr);
      for (const plane of Object.values(this.planes)) plane.material.uniforms.uDiscardEmpty.value = view === 'three';
      gl.render(this.scene, cam);
    }
    this.emit('render');
  }
}
