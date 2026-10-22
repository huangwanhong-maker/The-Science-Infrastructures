// The anatomical stage: the brain's surfaces drawn in three dimensions.
//
// World coordinates are RAS millimetres of the MNI 152 2009c space, the z axis up.
// Every surface of the anatomy product becomes an instance: a mesh with its own
// material, its base colour per vertex shaded by the baked openness (sulci and
// creases darker), lit by a studio environment and by key and rim lights that
// travel with the camera, so that the brain is always lit from the viewer's side.
//
// The viewer draws; the application decides what is shown and how coloured.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';
import { computeBoundsTree, disposeBoundsTree, acceleratedRaycast } from 'three-mesh-bvh';

THREE.BufferGeometry.prototype.computeBoundsTree = computeBoundsTree;
THREE.BufferGeometry.prototype.disposeBoundsTree = disposeBoundsTree;
THREE.Mesh.prototype.raycast = acceleratedRaycast;

const AXES = { sagittal: 0, coronal: 1, axial: 2 };
const ease = (t) => (t < 0.5 ? 4 * t * t * t : 1 - (-2 * t + 2) ** 3 / 2);

// A fresnel material: the cortex as glass, or any surface as a ghost.
function glassMaterial({ colour = '#f0cfe4', base = 0.035, rim = 0.5, power = 2.4 } = {}) {
  return new THREE.ShaderMaterial({
    clipping: true, transparent: true, depthWrite: false, side: THREE.FrontSide,
    uniforms: { uColor: { value: new THREE.Color(colour) }, uBase: { value: base }, uRim: { value: rim }, uPower: { value: power } },
    vertexShader: /* glsl */ `
      #include <common>
      #include <clipping_planes_pars_vertex>
      varying vec3 vN;
      varying vec3 vV;
      void main() {
        vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
        vN = normalize(normalMatrix * normal);
        vV = -mvPosition.xyz;
        gl_Position = projectionMatrix * mvPosition;
        #include <clipping_planes_vertex>
      }`,
    fragmentShader: /* glsl */ `
      #include <common>
      #include <clipping_planes_pars_fragment>
      uniform vec3 uColor;
      uniform float uBase, uRim, uPower;
      varying vec3 vN;
      varying vec3 vV;
      void main() {
        #include <clipping_planes_fragment>
        float f = pow(1.0 - abs(dot(normalize(vN), normalize(vV))), uPower);
        gl_FragColor = vec4(uColor * (0.75 + 0.5 * f), uBase + uRim * f);
        #include <colorspace_fragment>
      }`,
  });
}

function tissueMaterial(kind) {
  if (kind === 'ventricle') {
    return new THREE.MeshPhysicalMaterial({
      vertexColors: true, roughness: 0.18, metalness: 0, clearcoat: 1, clearcoatRoughness: 0.15,
      transparent: true, opacity: 0.62, depthWrite: false, envMapIntensity: 0.9, side: THREE.DoubleSide,
    });
  }
  if (kind === 'tract') {
    return new THREE.MeshPhysicalMaterial({
      vertexColors: true, roughness: 0.36, metalness: 0, clearcoat: 0.55, clearcoatRoughness: 0.3,
      sheen: 0.35, sheenRoughness: 0.4, sheenColor: new THREE.Color('#fff4e6'), envMapIntensity: 0.6,
    });
  }
  return new THREE.MeshPhysicalMaterial({
    vertexColors: true, roughness: 0.58, metalness: 0, clearcoat: 0.1, clearcoatRoughness: 0.55,
    sheen: 0.45, sheenRoughness: 0.55, sheenColor: new THREE.Color('#ffe2ec'), envMapIntensity: 0.5,
  });
}

export class Viewer {
  constructor(canvas, { onHover, onPick, onView } = {}) {
    this.canvas = canvas;
    this.onHover = onHover || (() => {});
    this.onPick = onPick || (() => {});
    this.onView = onView || (() => {});
    const renderer = this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true, powerPreference: 'high-performance' });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    renderer.setClearColor(0x000000, 0);
    renderer.toneMapping = THREE.NeutralToneMapping;
    renderer.toneMappingExposure = 1.02;
    renderer.outputColorSpace = THREE.SRGBColorSpace;
    renderer.localClippingEnabled = true;
    this.scene = new THREE.Scene();
    const pmrem = new THREE.PMREMGenerator(renderer);
    this.scene.environment = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
    this.scene.environmentIntensity = 0.55;
    this.camera = new THREE.PerspectiveCamera(28, 1, 5, 5000);
    this.camera.up.set(0, 0, 1);
    this.scene.add(this.camera);
    // Lights that travel with the camera: a warm key from the upper left, a lilac fill, a rim from behind.
    const key = new THREE.DirectionalLight('#fff0e6', 1.9);
    key.position.set(-0.6, 0.8, 0.9);
    const fill = new THREE.DirectionalLight('#d6c2ee', 0.32);
    fill.position.set(0.9, -0.3, 0.5);
    const rim = new THREE.DirectionalLight('#ffd6e8', 0.95);
    rim.position.set(0.15, 0.7, -1);
    this.camera.add(key, fill, rim);
    this.scene.add(new THREE.HemisphereLight('#fff4ef', '#2e1e2b', 0.22));
    this.root = new THREE.Group();
    this.scene.add(this.root);
    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.09;
    this.controls.rotateSpeed = 0.75;
    this.controls.zoomSpeed = 0.9;
    this.controls.minDistance = 40;
    this.controls.maxDistance = 1600;
    this.controls.addEventListener('change', () => { this.dirty = true; this.moving = true; this.onView('change'); });
    this.controls.addEventListener('start', () => { this.tween = null; this.onView('start'); });
    this.controls.addEventListener('end', () => { this.moving = false; this.onView('end'); });
    this.instances = new Map();          // key -> instance
    this.cut = null;
    this.clipPlane = new THREE.Plane(new THREE.Vector3(1, 0, 0), 0);
    this.explode = 0;
    this.explodeTarget = 0;
    this.glass = glassMaterial();
    this.ghost = glassMaterial({ colour: '#e9d6e6', base: 0.012, rim: 0.2, power: 2.0 });
    this.raycaster = new THREE.Raycaster();
    this.raycaster.firstHitOnly = true;
    this.pointer = null;
    this.dirty = true;
    this.bindPointer();
    new ResizeObserver(() => this.resize()).observe(canvas.parentElement);
    this.resize();
    this.loop = this.loop.bind(this);
    requestAnimationFrame(this.loop);
  }

  // ---- instances ----------------------------------------------------------------------------

  // A mesh of the anatomy product. entry: the catalogue's record; bytes: its pack.
  addMesh(entry, bytes, scale, { kind = 'tissue', id = null, side = null, colour = '#dddddd', clipped = false, render = 0 } = {}) {
    const n = entry.vertices;
    const base = bytes.byteOffset;
    const pos16 = new Int16Array(bytes.buffer, base + entry.position, n * 3);
    const pos = new Float32Array(n * 3);
    for (let i = 0; i < pos.length; i++) pos[i] = pos16[i] * scale;
    const n8 = new Int8Array(bytes.buffer, base + entry.normal, n * 4);
    const nrm = new Float32Array(n * 3);
    for (let i = 0; i < n; i++) { nrm[3 * i] = n8[4 * i] / 127; nrm[3 * i + 1] = n8[4 * i + 1] / 127; nrm[3 * i + 2] = n8[4 * i + 2] / 127; }
    const open = new Uint8Array(bytes.buffer, base + entry.openness, n).slice();
    const index = entry.indexType === 'uint16'
      ? new Uint16Array(bytes.buffer, base + entry.index, entry.triangles * 3).slice()
      : new Uint32Array(bytes.buffer, base + entry.index, entry.triangles * 3).slice();
    const labels = entry.label !== undefined ? new Uint8Array(bytes.buffer, base + entry.label, n).slice() : null;
    const depthRaw = entry.depth !== undefined ? new Uint8Array(bytes.buffer, base + entry.depth, n) : null;
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    g.setAttribute('normal', new THREE.BufferAttribute(nrm, 3));
    g.setAttribute('color', new THREE.BufferAttribute(new Float32Array(n * 3), 3));
    g.setIndex(new THREE.BufferAttribute(index, 1));
    g.computeBoundingBox();
    g.computeBoundingSphere();
    const material = tissueMaterial(kind);
    const mesh = new THREE.Mesh(g, material);
    mesh.renderOrder = render;
    mesh.userData.key = entry.key;
    // Shading: the sulcal depth (large surfaces) and the openness, each normalised over the
    // mesh between its 3rd and 97th percentiles, darken the vertex colours.
    const range = (arr) => {
      const sorted = Array.from(arr).sort((a, b) => a - b);
      const lo = sorted[Math.floor(0.03 * (n - 1))];
      return [lo, Math.max(lo + 1, sorted[Math.floor(0.97 * (n - 1))])];
    };
    const [olo, ohi] = range(open);
    const [dlo, dhi] = depthRaw ? range(depthRaw) : [0, 1];
    const shade = new Float32Array(n);
    const aoStrength = kind === 'cortex' ? 0.38 : kind === 'ventricle' ? 0.15 : 0.4;
    const depthStrength = kind === 'cortex' ? 0.66 : 0.5;
    for (let i = 0; i < n; i++) {
      const t = Math.min(1, Math.max(0, (open[i] - olo) / (ohi - olo)));
      let s = 1 - aoStrength * (1 - t) ** 1.3;
      if (depthRaw) {
        const d = Math.min(1, Math.max(0, (depthRaw[i] - dlo) / (dhi - dlo)));
        s *= 1 - depthStrength * d * d * (3 - 2 * d);
      }
      shade[i] = s;
    }
    const centroid = new THREE.Vector3();
    for (let i = 0; i < n; i++) centroid.x += pos[3 * i], centroid.y += pos[3 * i + 1], centroid.z += pos[3 * i + 2];
    centroid.divideScalar(n);
    const inst = { key: entry.key, entry, mesh, material, kind, id, side, colour, labels, shade, centroid, clipped, visible: true, ghost: false, offset: new THREE.Vector3(), labelIndex: null };
    this.instances.set(entry.key, inst);
    this.root.add(mesh);
    this.paint(inst);
    this.applyClip(inst);
    this.dirty = true;
    return inst;
  }

  // Colour an instance: one colour, or one per vertex from fn(vertexIndex) -> [r, g, b] (sRGB 0..1).
  paint(inst, fn = null, { emphasis = null } = {}) {
    const c = inst.mesh.geometry.attributes.color;
    const arr = c.array, shade = inst.shade, n = shade.length;
    const tmp = new THREE.Color();
    if (!fn) {
      tmp.set(inst.colour);
      for (let i = 0; i < n; i++) { const s = shade[i]; arr[3 * i] = tmp.r * s; arr[3 * i + 1] = tmp.g * s; arr[3 * i + 2] = tmp.b * s; }
    } else {
      const cache = new Map();
      for (let i = 0; i < n; i++) {
        const rgb = fn(i);
        let lin = cache.get(rgb);
        if (!lin) { tmp.setRGB(rgb[0], rgb[1], rgb[2], THREE.SRGBColorSpace); lin = [tmp.r, tmp.g, tmp.b]; cache.set(rgb, lin); }
        let s = shade[i];
        if (emphasis) s *= emphasis(i);
        arr[3 * i] = lin[0] * s; arr[3 * i + 1] = lin[1] * s; arr[3 * i + 2] = lin[2] * s;
      }
    }
    c.needsUpdate = true;
    this.dirty = true;
  }

  // Vertices on a border between labels (or between the keys keyOf gives), for engraving.
  borders(inst, keyOf = (i) => inst.labels[i]) {
    const idx = inst.mesh.geometry.index.array;
    const out = new Uint8Array(inst.shade.length);
    for (let t = 0; t < idx.length; t += 3) {
      const a = idx[t], b = idx[t + 1], c = idx[t + 2];
      const ka = keyOf(a), kb = keyOf(b), kc = keyOf(c);
      if (ka !== kb || kb !== kc) { out[a] = 1; out[b] = 1; out[c] = 1; }
    }
    return out;
  }

  // The vertex indices of each label of a labelled surface (the cortex's gyri).
  labelIndex(inst) {
    if (inst.labelIndex || !inst.labels) return inst.labelIndex;
    const groups = new Map();
    inst.labels.forEach((l, i) => { if (!groups.has(l)) groups.set(l, []); groups.get(l).push(i); });
    inst.labelIndex = new Map([...groups].map(([l, list]) => [l, Uint32Array.from(list)]));
    return inst.labelIndex;
  }

  // A mesh that is not of the anatomy product (the scalp, a region of an atlas): positions,
  // normals and indices in atlas millimetres. It is an instance like the others, without
  // labels or baked shading; colours may be given per vertex later with paint().
  addFreeMesh(key, { positions, normals, indices }, { kind = 'tissue', id = null, side = null, colour = '#e8d2c6', opacity = 1, pickable = true, clipped = false, render = 0 } = {}) {
    const old = this.instances.get(key);
    if (old) this.removeInstance(old);
    const n = positions.length / 3;
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    if (normals) g.setAttribute('normal', new THREE.BufferAttribute(normals, 3)); else g.computeVertexNormals();
    g.setAttribute('color', new THREE.BufferAttribute(new Float32Array(n * 3), 3));
    g.setIndex(new THREE.BufferAttribute(indices, 1));
    g.computeBoundingBox();
    g.computeBoundingSphere();
    const material = kind === 'skin'
      ? new THREE.MeshPhysicalMaterial({ vertexColors: true, roughness: 0.62, clearcoat: 0.25, clearcoatRoughness: 0.5, sheen: 0.4, sheenColor: new THREE.Color('#ffd9d0'), transparent: opacity < 1, opacity, depthWrite: opacity >= 1, side: THREE.DoubleSide })
      : tissueMaterial(kind);
    if (kind !== 'skin' && opacity < 1) { material.transparent = true; material.opacity = opacity; material.depthWrite = false; }
    const mesh = new THREE.Mesh(g, material);
    mesh.renderOrder = render;
    mesh.userData.key = key;
    const centroid = new THREE.Vector3();
    for (let i = 0; i < n; i++) centroid.x += positions[3 * i], centroid.y += positions[3 * i + 1], centroid.z += positions[3 * i + 2];
    centroid.divideScalar(n);
    const inst = { key, entry: null, mesh, material, kind, id, side, colour, labels: null, shade: new Float32Array(n).fill(1), centroid, clipped, visible: true, ghost: false, offset: new THREE.Vector3(), labelIndex: null, pickable, free: true, style: 'solid' };
    this.instances.set(key, inst);
    this.root.add(mesh);
    this.paint(inst);
    this.applyClip(inst);
    this.dirty = true;
    return inst;
  }

  removeInstance(inst) {
    this.root.remove(inst.mesh);
    inst.mesh.geometry.dispose();
    if (inst.free) inst.material.dispose();
    this.instances.delete(inst.key);
    this.dirty = true;
  }

  // Markers: small spheres at points (electrodes, landmarks, data points), each
  // {position, colour, radius}; drawn above the surfaces they sit on.
  setMarkers(list) {
    if (this.markers) { this.root.remove(this.markers); this.markers.geometry.dispose(); this.markers.material.dispose(); this.markers = null; }
    this.markerList = list || [];
    if (!this.markerList.length) { this.dirty = true; return; }
    const geo = new THREE.SphereGeometry(1, 20, 14);
    const mat = new THREE.MeshPhysicalMaterial({ roughness: 0.35, clearcoat: 0.6, clearcoatRoughness: 0.25, metalness: 0 });
    const m = new THREE.InstancedMesh(geo, mat, this.markerList.length);
    const M = new THREE.Matrix4(), c = new THREE.Color();
    this.markerList.forEach((k, i) => {
      M.makeScale(k.radius || 2.5, k.radius || 2.5, k.radius || 2.5).setPosition(...k.position);
      m.setMatrixAt(i, M);
      m.setColorAt(i, c.set(k.colour || '#f6c47e'));
    });
    m.renderOrder = 5;
    this.markers = m;
    this.root.add(m);
    this.dirty = true;
  }

  // The marker under a screen point: its index in the list, or -1.
  pickMarker(clientX, clientY) {
    if (!this.markers) return -1;
    const r = this.canvas.getBoundingClientRect();
    const ndc = new THREE.Vector2(((clientX - r.left) / r.width) * 2 - 1, -((clientY - r.top) / r.height) * 2 + 1);
    this.raycaster.setFromCamera(ndc, this.camera);
    const hits = this.raycaster.intersectObject(this.markers, false);
    return hits.length ? hits[0].instanceId : -1;
  }

  // ---- the scale of a subject's head --------------------------------------------------------

  // Scale everything drawn about a centre (atlas mm), by [sx, sy, sz]: a template brought
  // to a subject's head. Picking and projection convert between the two.
  setScale(scales, centre) {
    this.root.scale.set(...scales);
    this.root.position.set(centre[0] * (1 - scales[0]), centre[1] * (1 - scales[1]), centre[2] * (1 - scales[2]));
    this.root.updateMatrixWorld(true);
    this.setCut(this.cut);
    this.dirty = true;
  }

  setVisible(inst, visible) {
    inst.visible = visible;
    inst.mesh.visible = visible;
    this.dirty = true;
  }

  // Draw an instance solid, as glass (the cortex), or as a ghost (dimmed behind a focus).
  setStyle(inst, style) {
    inst.style = style;
    inst.mesh.material = style === 'glass' ? this.glass : style === 'ghost' ? this.ghost : inst.material;
    inst.mesh.renderOrder = style === 'solid' ? (inst.kind === 'ventricle' ? 2 : 0) : 3;
    this.applyClip(inst);
    this.dirty = true;
  }

  setEmphasis(inst, amount, colour = null) {
    const m = inst.material;
    if (!m.emissive) return;
    m.emissive.set(colour || inst.colour);
    m.emissiveIntensity = amount;
    this.dirty = true;
  }

  // ---- the cut --------------------------------------------------------------------------------

  // cut: null, or {axis: 'sagittal' | 'coronal' | 'axial', position (mm), keep: 1 | -1, all: bool}
  setCut(cut) {
    this.cut = cut;
    if (cut) {
      const normal = new THREE.Vector3();
      normal.setComponent(AXES[cut.axis], cut.keep);
      // Plane: normal . x + constant >= 0 is kept, in atlas coordinates, then carried to the world.
      this.clipPlane.set(normal, -cut.keep * cut.position);
      this.root.updateMatrixWorld(true);
      this.clipPlane.applyMatrix4(this.root.matrixWorld);
    }
    for (const inst of this.instances.values()) this.applyClip(inst);
    this.updateCap();
    this.dirty = true;
  }

  applyClip(inst) {
    const on = this.cut && (inst.clipped || this.cut.all);
    const planes = on ? [this.clipPlane] : [];
    inst.material.clippingPlanes = planes;
    if (inst.mesh.material !== inst.material) {
      // Shared glass and ghost materials are clipped whenever there is a cut.
      inst.mesh.material.clippingPlanes = this.cut ? [this.clipPlane] : [];
    }
    inst.material.needsUpdate = true;
  }

  // The section at the cut, drawn as its cap: canvas is the section's image, bounds its extent.
  setCapImage(canvas, rect) {
    if (!this.cap) {
      const tex = new THREE.CanvasTexture(canvas);
      tex.colorSpace = THREE.SRGBColorSpace;
      tex.anisotropy = 4;
      const m = new THREE.MeshBasicMaterial({ map: tex, transparent: true, side: THREE.DoubleSide, toneMapped: false, depthWrite: true, alphaTest: 0.05 });
      this.cap = new THREE.Mesh(new THREE.PlaneGeometry(1, 1), m);
      this.cap.renderOrder = 1;
      this.root.add(this.cap);
    } else if (this.cap.material.map.image !== canvas) {
      this.cap.material.map.image = canvas;
    }
    this.cap.material.map.needsUpdate = true;
    this.capRect = rect;
    this.updateCap();
  }

  updateCap() {
    if (!this.cap) return;
    const c = this.cut;
    this.cap.visible = Boolean(c && c.cap && this.capRect);
    if (!this.cap.visible) { this.dirty = true; return; }
    const r = this.capRect;   // {u: [min, max], v: [min, max]} in the plane's two axes, mm
    const w = r.u[1] - r.u[0], h = r.v[1] - r.v[0];
    const cu = (r.u[0] + r.u[1]) / 2, cv = (r.v[0] + r.v[1]) / 2;
    this.cap.scale.set(w, h, 1);
    this.cap.rotation.set(0, 0, 0);
    const eps = 0.15 * -c.keep;  // a hair to the removed side, so the cap is not clipped by its own plane's surfaces
    if (c.axis === 'axial') { this.cap.position.set(cu, cv, c.position + eps); }
    if (c.axis === 'coronal') { this.cap.rotation.x = Math.PI / 2; this.cap.position.set(cu, c.position + eps, cv); }
    if (c.axis === 'sagittal') { this.cap.rotation.set(Math.PI / 2, Math.PI / 2, 0); this.cap.position.set(c.position + eps, cu, cv); }
    this.dirty = true;
  }

  // ---- explosion --------------------------------------------------------------------------

  setExplode(amount) { this.explodeTarget = amount; this.dirty = true; }

  // ---- the camera ---------------------------------------------------------------------------

  // Fly to look at target from a direction (unit vector from target to the camera) at a distance.
  flyTo({ target, direction, distance, up = [0, 0, 1], duration = 900 }) {
    const from = { pos: this.camera.position.clone(), target: this.controls.target.clone(), up: this.camera.up.clone() };
    const t = target ? new THREE.Vector3(...target) : this.controls.target.clone();
    const dir = direction ? new THREE.Vector3(...direction).normalize() : this.camera.position.clone().sub(this.controls.target).normalize();
    const dist = distance || this.camera.position.distanceTo(this.controls.target);
    const to = { pos: t.clone().add(dir.multiplyScalar(dist)), target: t, up: new THREE.Vector3(...up) };
    this.tween = { from, to, start: performance.now(), duration };
    this.dirty = true;
  }

  // The distance at which a sphere of radius r fills the view.
  fitDistance(radius, margin = 1.25) {
    const fov = THREE.MathUtils.degToRad(this.camera.fov);
    const aspect = Math.max(0.2, this.camera.aspect);
    const fit = Math.min(fov / 2, Math.atan(Math.tan(fov / 2) * aspect));
    return (radius * margin) / Math.sin(fit);
  }

  frame(box, { direction = null, margin = 1.25, duration = 900 } = {}) {
    const centre = new THREE.Vector3(), size = new THREE.Vector3();
    box.getCenter(centre); box.getSize(size);
    const radius = Math.max(6, size.length() / 2);
    this.flyTo({ target: centre.toArray(), direction, distance: this.fitDistance(radius, margin), duration });
  }

  view(direction, up = [0, 0, 1]) {
    this.flyTo({ direction, up, duration: 750 });
  }

  // ---- picking ------------------------------------------------------------------------------

  bindPointer() {
    const el = this.canvas;
    let down = null;
    el.addEventListener('pointerdown', (e) => { down = { x: e.clientX, y: e.clientY, t: performance.now() }; });
    el.addEventListener('pointermove', (e) => {
      this.pointer = { x: e.clientX, y: e.clientY };
      this.hoverWanted = true;
    });
    el.addEventListener('pointerleave', () => { this.pointer = null; this.hoverWanted = true; });
    el.addEventListener('pointerup', (e) => {
      if (!down) return;
      const moved = Math.hypot(e.clientX - down.x, e.clientY - down.y);
      const quick = performance.now() - down.t < 450;
      down = null;
      if (moved < 5 && quick) this.onPick(this.pick(e.clientX, e.clientY), e);
    });
    el.addEventListener('dblclick', (e) => this.onPick(this.pick(e.clientX, e.clientY), e, true));
  }

  // The first visible surface under a screen point: {inst, point, label} or null.
  pick(clientX, clientY) {
    const r = this.canvas.getBoundingClientRect();
    const ndc = new THREE.Vector2(((clientX - r.left) / r.width) * 2 - 1, -((clientY - r.top) / r.height) * 2 + 1);
    this.raycaster.setFromCamera(ndc, this.camera);
    const targets = [];
    for (const inst of this.instances.values()) if (inst.visible && inst.style !== 'ghost' && inst.pickable !== false) targets.push(inst.mesh);
    for (const m of targets) if (!m.geometry.boundsTree) m.geometry.computeBoundsTree();
    let hits = this.raycaster.intersectObjects(targets, false);
    if (this.cut) hits = hits.filter((h) => !(h.object.material.clippingPlanes?.length) || this.clipPlane.distanceToPoint(h.point) >= -0.01);
    // A glass cortex lets the pointer through to what lies inside, unless nothing does.
    const solid = hits.find((h) => this.instanceOf(h.object).style !== 'glass');
    const hit = solid || hits[0];
    if (!hit) return null;
    const inst = this.instanceOf(hit.object);
    let label = null;
    if (inst.labels && hit.face) {
      const { a, b, c } = hit.face;
      const ls = [inst.labels[a], inst.labels[b], inst.labels[c]];
      label = ls[0] === ls[1] || ls[0] === ls[2] ? ls[0] : ls[1] === ls[2] ? ls[1] : ls[0];
    }
    const local = this.root.worldToLocal(hit.point.clone()).sub(inst.mesh.position);
    return { inst, point: local.toArray(), label, distance: hit.distance, face: hit.face };
  }

  instanceOf(mesh) { return this.instances.get(mesh.userData.key); }

  // Whether a world point on an instance is seen from the camera (not hidden behind another surface).
  visibleAt(point, inst) {
    const p = this.root.localToWorld(new THREE.Vector3(...point).add(inst.mesh.position));
    const dir = p.clone().sub(this.camera.position);
    const dist = dir.length();
    this.raycaster.set(this.camera.position, dir.normalize());
    const targets = [];
    for (const i of this.instances.values()) if (i.visible && i.style === 'solid' && i.pickable !== false) targets.push(i.mesh);
    for (const m of targets) if (!m.geometry.boundsTree) m.geometry.computeBoundsTree();
    // The face of a cut hides what lies behind it.
    if (this.cap && this.cap.visible) targets.push(this.cap);
    const hits = this.raycaster.intersectObjects(targets, false).filter((h) => h.object === this.cap || !this.cut || !(h.object.material.clippingPlanes?.length) || this.clipPlane.distanceToPoint(h.point) >= -0.01);
    if (!hits.length) return true;
    return hits[0].distance >= dist - 1.2;
  }

  // Screen position (CSS pixels within the canvas) of a world point; null when behind the camera.
  project(point, inst = null) {
    const p = new THREE.Vector3(...point);
    if (inst) p.add(inst.mesh.position);
    this.root.localToWorld(p);
    p.project(this.camera);
    if (p.z > 1) return null;
    return { x: ((p.x + 1) / 2) * this.width, y: ((1 - p.y) / 2) * this.height };
  }

  // The camera's axes in world coordinates: right, up and toward the viewer.
  basis() {
    const m = this.camera.matrixWorld.elements;
    return { right: [m[0], m[1], m[2]], up: [m[4], m[5], m[6]], back: [m[8], m[9], m[10]] };
  }

  // ---- the loop -----------------------------------------------------------------------------

  resize() {
    const el = this.canvas.parentElement;
    const w = Math.max(1, el.clientWidth), h = Math.max(1, el.clientHeight);
    this.width = w; this.height = h;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.dirty = true;
  }

  loop(now) {
    requestAnimationFrame(this.loop);
    if (this.tween) {
      const { from, to, start, duration } = this.tween;
      const t = Math.min(1, (now - start) / duration), k = ease(t);
      this.camera.position.lerpVectors(from.pos, to.pos, k);
      this.controls.target.lerpVectors(from.target, to.target, k);
      this.camera.up.lerpVectors(from.up, to.up, k).normalize();
      this.camera.lookAt(this.controls.target);
      this.dirty = true;
      if (t >= 1) { this.tween = null; this.onView('end'); } else this.onView('change');
    }
    if (Math.abs(this.explode - this.explodeTarget) > 1e-3) {
      this.explode += (this.explodeTarget - this.explode) * 0.14;
      if (Math.abs(this.explode - this.explodeTarget) < 1e-3) this.explode = this.explodeTarget;
      for (const inst of this.instances.values()) inst.mesh.position.copy(inst.offset).multiplyScalar(this.explode);
      this.dirty = true;
      this.onView(this.explode === this.explodeTarget ? 'end' : 'change');
    }
    const changed = this.controls.update();
    if (this.hoverWanted) {
      this.hoverWanted = false;
      if (!this.moving && !this.tween) this.onHover(this.pointer ? this.pick(this.pointer.x, this.pointer.y) : null, this.pointer);
    }
    if (this.dirty || changed) {
      this.dirty = false;
      this.renderer.render(this.scene, this.camera);
      if (this.afterRender) this.afterRender();
    }
  }
}
