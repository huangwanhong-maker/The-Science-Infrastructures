// The round view: the body's meshes in three.js, coloured by system; glass, hidden and isolated parts; the
// chosen parts lit in rose gold; picking by ray against a bounding-volume hierarchy built for each mesh the
// first time a ray reaches it; a cut by a cardinal plane whose face is the engine's section, drawn exactly
// where the domain cut the meshes; a compass that names the subject's own directions on the screen.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';
import { acceleratedRaycast, computeBoundsTree } from 'three-mesh-bvh';

import { lpsToWorld } from './data.js';
import { systemColour } from './palette.js';

THREE.Mesh.prototype.raycast = acceleratedRaycast;
THREE.BufferGeometry.prototype.computeBoundsTree = computeBoundsTree;

const ROSE = new THREE.Color('#c99278');
const AXIS = { sagittal: 0, coronal: 1, transverse: 2 };

export class Viewer {
  constructor(container) {
    this.container = container;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    this.renderer.localClippingEnabled = true;
    this.renderer.toneMapping = THREE.NeutralToneMapping;
    this.renderer.toneMappingExposure = 1.0;
    container.appendChild(this.renderer.domElement);

    this.scene = new THREE.Scene();
    const pmrem = new THREE.PMREMGenerator(this.renderer);
    this.scene.environment = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
    this.scene.environmentIntensity = 0.55;
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0xc8ccd2, 1.1));
    const key = new THREE.DirectionalLight(0xffffff, 1.6);
    key.position.set(-900, 1400, 1600);
    const fill = new THREE.DirectionalLight(0xf3eee9, 0.55);
    fill.position.set(1200, 300, 900);
    const rim = new THREE.DirectionalLight(0xe8eef6, 0.7);
    rim.position.set(0, 900, -1800);
    this.scene.add(key, fill, rim);

    this.camera = new THREE.PerspectiveCamera(28, 1, 10, 30000);
    this.camera.position.set(0, 120, 3900);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.target.set(0, 60, 0);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.12;
    this.controls.zoomSpeed = 0.9;
    this.controls.minDistance = 60;
    this.controls.maxDistance = 9000;
    this.controls.addEventListener('change', () => this.request());

    this.root = new THREE.Group();
    this.scene.add(this.root);
    this.caps = new THREE.Group();
    this.scene.add(this.caps);

    this.meshes = [];            // per element
    this.system = [];            // per element
    this.styles = {};            // system → 'solid' | 'glass' | 'off'
    this.hidden = new Set();     // elements hidden by hand, or copies
    this.isolated = null;        // a Set of elements drawn alone, the rest as ghosts
    this.selected = new Set();
    this.materials = {};
    this.ghost = new THREE.MeshStandardMaterial({ color: 0xb8bec6, transparent: true, opacity: 0.045, depthWrite: false, roughness: 0.5 });
    this.clip = null;            // THREE.Plane when cut
    this.raycaster = new THREE.Raycaster();
    this.raycaster.firstHitOnly = true;
    this.onCompass = null;

    this.resize();
    new ResizeObserver(() => this.resize()).observe(container);
    this.animate = this.animate.bind(this);
    this.dirty = true;
    requestAnimationFrame(this.animate);
  }

  request() { this.dirty = true; }

  resize() {
    const w = this.container.clientWidth || 1, h = this.container.clientHeight || 1;
    this.renderer.setSize(w, h, false);
    this.renderer.domElement.style.width = `${w}px`;
    this.renderer.domElement.style.height = `${h}px`;
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.request();
  }

  animate(t) {
    requestAnimationFrame(this.animate);
    if (this.flight) this.step(t);
    const moved = this.controls.update();
    if (!this.dirty && !moved && !this.flight) return;
    this.dirty = false;
    this.renderer.render(this.scene, this.camera);
    this.compass();
  }

  material(system, style) {
    const key = `${system}:${style}`;
    if (!this.materials[key]) {
      const colour = new THREE.Color(systemColour(system));
      const m = style === 'glass'
        ? new THREE.MeshStandardMaterial({ color: colour, transparent: true, opacity: system === 'skin' ? 0.1 : 0.2, depthWrite: false, roughness: 0.35 })
        : new THREE.MeshStandardMaterial({ color: colour, roughness: system === 'skeleton' ? 0.78 : 0.58, metalness: 0 });
      m.clippingPlanes = this.clip ? [this.clip] : [];
      this.materials[key] = m;
    }
    return this.materials[key];
  }

  addElement(i, decoded, system, { inward = false, solid = true } = {}) {
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(decoded.position, 3));
    let normal = decoded.normal;
    let index = decoded.index;
    if (inward) {
      normal = normal.map((x) => -x);
      for (let t = 0; t < index.length; t += 3) { const a = index[t + 1]; index[t + 1] = index[t + 2]; index[t + 2] = a; }
    }
    g.setAttribute('normal', new THREE.BufferAttribute(normal, 3, true));
    g.setIndex(new THREE.BufferAttribute(index, 1));
    g.computeBoundingSphere();
    g.computeBoundingBox();
    const mesh = new THREE.Mesh(g, this.material(system, 'solid'));
    mesh.userData = { i, system, solid };
    mesh.matrixAutoUpdate = false;
    this.meshes[i] = mesh;
    this.system[i] = system;
    this.root.add(mesh);
    this.restyle(i);
    this.request();
  }

  setStyle(system, style) {
    this.styles[system] = style;
    for (let i = 0; i < this.meshes.length; i++) if (this.system[i] === system) this.restyle(i);
    this.request();
  }

  restyle(i) {
    const mesh = this.meshes[i];
    if (!mesh) return;
    const style = this.styles[this.system[i]] || 'solid';
    let visible = style !== 'off' && !this.hidden.has(i);
    let material;
    if (this.isolated) {
      if (this.isolated.has(i)) { visible = !this.hidden.has(i); material = this.material(this.system[i], 'solid'); }
      else material = this.ghost;
    } else {
      material = this.material(this.system[i], style === 'glass' ? 'glass' : 'solid');
    }
    if (this.selected.has(i) && visible) {
      if (!mesh.userData.lit) {
        mesh.userData.lit = material.clone();
        mesh.userData.lit.emissive = ROSE.clone();
        mesh.userData.lit.emissiveIntensity = 0.42;
      }
      mesh.userData.lit.clippingPlanes = this.clip ? [this.clip] : [];
      mesh.userData.lit.opacity = material.opacity;
      mesh.userData.lit.transparent = material.transparent;
      mesh.material = mesh.userData.lit;
    } else {
      mesh.material = material;
    }
    if (!mesh.userData.solid) mesh.material.side = THREE.DoubleSide;
    mesh.visible = visible;
  }

  restyleAll() {
    for (let i = 0; i < this.meshes.length; i++) this.restyle(i);
    this.request();
  }

  // The chosen parts are lit in rose gold where they are seen, and drawn faintly through whatever lies in
  // front of them, so that a part inside the body is found at once.
  select(elements) {
    for (const i of this.selected) if (this.meshes[i]) { this.meshes[i].userData.lit = null; }
    this.selected = new Set(elements);
    for (const o of this.overlays || []) this.scene.remove(o);
    this.overlays = [];
    if (!this.through) this.through = new THREE.MeshBasicMaterial({ color: 0xfff1ea, transparent: true, opacity: 0.16, depthTest: false, depthWrite: false });
    this.through.clippingPlanes = this.clip ? [this.clip] : [];
    for (const i of elements) {
      const mesh = this.meshes[i];
      if (!mesh) continue;
      const o = new THREE.Mesh(mesh.geometry, this.through);
      o.renderOrder = 1000;
      o.matrixAutoUpdate = false;
      o.raycast = () => {};
      this.overlays.push(o);
      this.scene.add(o);
    }
    this.restyleAll();
  }

  isolate(elements) {
    this.isolated = elements && elements.length ? new Set(elements) : null;
    this.restyleAll();
  }

  isVisible(i) { return !!(this.meshes[i] && this.meshes[i].visible); }

  // ---- picking ----------------------------------------------------------------------------------------------

  pick(clientX, clientY) {
    const r = this.renderer.domElement.getBoundingClientRect();
    const ndc = new THREE.Vector2(((clientX - r.left) / r.width) * 2 - 1, -((clientY - r.top) / r.height) * 2 + 1);
    this.raycaster.setFromCamera(ndc, this.camera);
    const ray = this.raycaster.ray;
    const sphere = new THREE.Sphere();
    const candidates = [[], []];  // opaque, then glass and ghosts
    for (const mesh of this.meshes) {
      if (!mesh || !mesh.visible) continue;
      sphere.copy(mesh.geometry.boundingSphere);
      if (!ray.intersectsSphere(sphere)) continue;
      if (!mesh.geometry.boundsTree) mesh.geometry.computeBoundsTree();
      candidates[mesh.material.transparent ? 1 : 0].push(mesh);
    }
    for (const list of candidates) {
      let best = null;
      for (const mesh of list) {
        const hits = this.raycaster.intersectObject(mesh, false);
        for (const h of hits) {
          if (this.clip && this.clip.distanceToPoint(h.point) < 0) continue;
          if (!best || h.distance < best.distance) best = h;
          break;
        }
      }
      if (this.clip && this.capMeshes) {
        const capHits = this.raycaster.intersectObjects(this.capMeshes, false);
        if (capHits.length && (!best || capHits[0].distance < best.distance)) best = { ...capHits[0], object: { userData: { i: capHits[0].object.userData.i } } };
      }
      if (best) return { i: best.object.userData.i, point: best.point };
    }
    return null;
  }

  // ---- the camera -------------------------------------------------------------------------------------------

  frame(box, { duration = 900, direction = null } = {}) {
    const centre = box.getCenter(new THREE.Vector3());
    const size = Math.max(420, box.getSize(new THREE.Vector3()).length());
    const dir = direction ? direction.clone().normalize() : this.camera.position.clone().sub(this.controls.target).normalize();
    const dist = Math.max(180, (size / 2) / Math.tan((this.camera.fov * Math.PI) / 360) * 1.08);
    this.flyTo(centre, centre.clone().add(dir.multiplyScalar(dist)), duration);
  }

  flyTo(target, position, duration = 900) {
    this.flight = { t0: null, duration, from: [this.controls.target.clone(), this.camera.position.clone()], to: [target.clone(), position.clone()] };
    this.request();
  }

  step(t) {
    const f = this.flight;
    if (f.t0 === null) f.t0 = t;
    const u = Math.min(1, (t - f.t0) / f.duration);
    const e = u < 0.5 ? 4 * u * u * u : 1 - Math.pow(-2 * u + 2, 3) / 2;
    this.controls.target.lerpVectors(f.from[0], f.to[0], e);
    this.camera.position.lerpVectors(f.from[1], f.to[1], e);
    if (u >= 1) this.flight = null;
    this.dirty = true;
  }

  boxOf(elements) {
    const box = new THREE.Box3();
    for (const i of elements) if (this.meshes[i]) box.union(this.meshes[i].geometry.boundingBox);
    return box;
  }

  view() {
    const r = (v) => v.toArray().map((x) => Math.round(x));
    return { p: r(this.camera.position), t: r(this.controls.target) };
  }

  setView(v) {
    if (!v) return;
    this.camera.position.fromArray(v.p);
    this.controls.target.fromArray(v.t);
    this.controls.update();
    this.request();
  }

  // ---- the cut ----------------------------------------------------------------------------------------------

  // Keeps the half-space on the side `keep` (+1 or −1 along the plane's LPS axis) of the plane at `level`.
  setCut(plane, level, keep = -1) {
    if (!plane) {
      this.clip = null;
    } else {
      const n = [0, 0, 0];
      n[AXIS[plane]] = keep;
      const point = [0, 0, 0];
      point[AXIS[plane]] = level;
      const nw = new THREE.Vector3(...lpsToWorld(n)).sub(new THREE.Vector3(...lpsToWorld([0, 0, 0]))).normalize();
      const pw = new THREE.Vector3(...lpsToWorld(point));
      this.clip = new THREE.Plane().setFromNormalAndCoplanarPoint(nw, pw);
    }
    for (const m of Object.values(this.materials)) m.clippingPlanes = this.clip ? [this.clip] : [];
    this.ghost.clippingPlanes = this.clip ? [this.clip] : [];
    for (const mesh of this.meshes) if (mesh && mesh.userData.lit) mesh.userData.lit.clippingPlanes = this.clip ? [this.clip] : [];
    if (this.through) this.through.clippingPlanes = this.clip ? [this.clip] : [];
    if (!plane) this.setCaps(null);
    this.request();
  }

  // The face of the cut: the engine's section, its pieces' loops (u, v in the display's LPS basis) laid on the plane.
  setCaps(section, colourOf = (i) => systemColour(this.system[i])) {
    for (const c of this.caps.children) c.geometry.dispose();
    this.caps.clear();
    this.capMeshes = [];
    if (!section) { this.request(); return; }
    const U = section.u, V = section.v, a = AXIS[section.plane];
    const at = (u, v) => {
      const p = [u * U[0] + v * V[0], u * U[1] + v * V[1], u * U[2] + v * V[2]];
      p[a] = section.level;
      return new THREE.Vector3(...lpsToWorld(p));
    };
    // Larger pieces first: the caps lie in one plane, and a later one is drawn over an earlier.
    const pieces = [...section.pieces].sort((a, b) => b.area - a.area);
    let order = 0;
    for (const piece of pieces) {
      if (!this.isVisible(piece.e) && !(this.isolated && this.meshes[piece.e])) continue;
      if (this.isolated && !this.isolated.has(piece.e)) continue;
      const loops = piece.loops.map((flat) => { const l = []; for (let k = 0; k < flat.length; k += 2) l.push(new THREE.Vector2(flat[k], flat[k + 1])); return l; });
      const area = (l) => { let s = 0; for (let k = 0; k < l.length; k++) { const p = l[k], q = l[(k + 1) % l.length]; s += p.x * q.y - q.x * p.y; } return s / 2; };
      const outers = loops.filter((l) => area(l) > 0), holes = loops.filter((l) => area(l) < 0);
      const groups = outers.map((o) => ({ outer: o, holes: [] }));
      for (const h of holes) {
        let best = null, bestArea = Infinity;
        for (const g of groups) if (inside(h[0], g.outer) && area(g.outer) < bestArea) { best = g; bestArea = area(g.outer); }
        if (best) best.holes.push(h);
      }
      const positions = [];
      for (const g of groups) {
        const tris = THREE.ShapeUtils.triangulateShape(g.outer, g.holes);
        const all = g.outer.concat(...g.holes);
        for (const t of tris) for (const k of t) { const p = at(all[k].x, all[k].y); positions.push(p.x, p.y, p.z); }
      }
      if (!positions.length) continue;
      const geom = new THREE.BufferGeometry();
      geom.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
      const colour = new THREE.Color(colourOf(piece.e)).lerp(new THREE.Color(0xffffff), 0.12);
      const lit = this.selected.has(piece.e);
      const mat = new THREE.MeshBasicMaterial({ color: lit ? colour.clone().lerp(ROSE, 0.45) : colour, side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: -1, polygonOffsetUnits: -1 });
      const m = new THREE.Mesh(geom, mat);
      m.userData.i = piece.e;
      m.renderOrder = order++;
      this.caps.add(m);
      this.capMeshes.push(m);
    }
    this.request();
  }

  // ---- the compass ------------------------------------------------------------------------------------------

  compass() {
    if (!this.onCompass) return;
    const m = this.camera.matrixWorld.elements;
    const right = [m[0], m[1], m[2]], up = [m[4], m[5], m[6]];
    // World x = left, y = superior, z = anterior.
    const name = (v) => {
      const k = [0, 1, 2].reduce((a, b) => (Math.abs(v[b]) > Math.abs(v[a]) ? b : a), 0);
      return [['L', 'R'], ['S', 'I'], ['A', 'P']][k][v[k] > 0 ? 0 : 1];
    };
    const opp = { L: 'R', R: 'L', S: 'I', I: 'S', A: 'P', P: 'A' };
    const r = name(right), u = name(up);
    this.onCompass({ right: r, left: opp[r], up: u, down: opp[u] });
  }

  snapshot() { this.renderer.render(this.scene, this.camera); return this.renderer.domElement.toDataURL('image/png'); }
}

function inside(p, ring) {
  let c = false;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    const a = ring[i], b = ring[j];
    if ((a.y > p.y) !== (b.y > p.y) && p.x < ((b.x - a.x) * (p.y - a.y)) / (b.y - a.y) + a.x) c = !c;
  }
  return c;
}
