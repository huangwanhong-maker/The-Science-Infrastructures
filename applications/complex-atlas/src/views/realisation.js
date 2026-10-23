// The realisation: a complex in space with three.js. Its faces are translucent discs (realise()
// in visualization/js/complex.js), its edges curves of a chosen width with an arrow for their
// direction, its vertices small spheres; a field colours the edges by the angle of their link
// variables and the faces by their plaquettes; a selection (a chain or a loop) is drawn over it.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { LineSegments2 } from 'three/addons/lines/LineSegments2.js';
import { LineSegmentsGeometry } from 'three/addons/lines/LineSegmentsGeometry.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';
import { cyclicColour, divergingColour } from '../../../../visualization/js/index.js';
import { selectionSets } from '../selection.js';

const FACE = new THREE.Color('#c08081');
const EDGE = new THREE.Color('#ecd4e2');
const VERTEX = new THREE.Color('#f2c27c');
const MARK = new THREE.Color('#f6cf7e');
const MARK_NEG = new THREE.Color('#b9a2e0');

export class Realisation {
  constructor(host, { onHover } = {}) {
    this.host = host;
    this.onHover = onHover;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    this.renderer.setClearColor(0x000000, 0);
    host.append(this.renderer.domElement);
    this.labels = new CSS2DRenderer();
    this.labels.domElement.className = 'labels';
    host.append(this.labels.domElement);
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(36, 1, 0.01, 100);
    this.camera.position.set(1.7, 1.9, 2.0);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.08;
    this.scene.add(new THREE.HemisphereLight(0xf6e7f0, 0x2a1f2f, 1.4));
    const key = new THREE.DirectionalLight(0xfff1e6, 1.3);
    key.position.set(2, 3, 4);
    this.camera.add(key);
    this.scene.add(this.camera);
    this.group = new THREE.Group();
    this.scene.add(this.group);
    this.raycaster = new THREE.Raycaster();
    this.pointer = null;
    this.renderer.domElement.addEventListener('pointermove', (e) => {
      const r = this.renderer.domElement.getBoundingClientRect();
      this.pointer = new THREE.Vector2(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
    });
    this.renderer.domElement.addEventListener('pointerleave', () => { this.pointer = null; this.onHover?.(null); });
    new ResizeObserver(() => this.resize()).observe(host);
    this.resize();
    const loop = () => {
      this.frame = requestAnimationFrame(loop);
      this.controls.update();
      this.pick();
      this.renderer.render(this.scene, this.camera);
      this.labels.render(this.scene, this.camera);
    };
    loop();
  }

  resize() {
    const w = Math.max(1, this.host.clientWidth), h = Math.max(1, this.host.clientHeight);
    this.renderer.setSize(w, h);
    this.labels.setSize(w, h);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    for (const m of [this.edgeMaterial, this.markMaterial]) m?.resolution.set(w, h);
  }

  // The complex, its geometry (from realise()) and the style: {vertexSize, edgeWidth, faceOpacity,
  // labels, arrows}. The camera is kept unless `fit` is set.
  show(complex, geometry, style, fit = true) {
    this.complex = complex;
    this.geometry = geometry;
    this.style = style;
    this.rebuild();
    if (fit) { this.camera.position.set(1.7, 1.9, 2.0); this.controls.target.set(0, 0, 0); }
  }
  setStyle(style) { this.style = style; this.rebuild(); }
  setFields(fields) { this.fields = fields; this.rebuild(); }
  setNames(names) { this.names = names; }
  highlight(selection) { this.selection = selection; this.rebuild(); }

  clearGroup() {
    for (const o of [...this.group.children]) {
      o.traverse((x) => {
        x.geometry?.dispose?.();
        if (x.material) (Array.isArray(x.material) ? x.material : [x.material]).forEach((m) => m.dispose());
        if (x.isCSS2DObject) x.element.remove();
      });
      this.group.remove(o);
    }
  }

  rebuild() {
    this.clearGroup();
    const g = this.geometry, c = this.complex, st = this.style;
    if (!g || !c) return;
    const unit = Math.max(0.02, Math.min(0.3, g.unit || 0.2));
    const fields = this.fields && this.fields.link_angle?.length === c.edges.length ? this.fields : null;
    const marks = selectionSets(this.selection, c);

    // Faces, one mesh; the triangles of each face remembered for picking.
    const pos = [], col = [], idx = [];
    this.faceOfTriangle = [];
    g.faces.forEach((f, i) => {
      const base = pos.length / 3;
      const colour = fields && fields.plaquettes?.[i] !== undefined ? new THREE.Color(...divergingColour(fields.plaquettes[i])) : FACE;
      for (let j = 0; j < f.positions.length; j++) pos.push(f.positions[j]);
      for (let j = 0; j < f.positions.length / 3; j++) col.push(colour.r, colour.g, colour.b);
      for (let j = 0; j < f.indices.length; j++) idx.push(base + f.indices[j]);
      for (let j = 0; j < f.indices.length / 3; j++) this.faceOfTriangle.push(i);
    });
    if (idx.length) {
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
      geo.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
      geo.setIndex(idx);
      geo.computeVertexNormals();
      const mat = new THREE.MeshStandardMaterial({ vertexColors: true, transparent: true, opacity: st.faceOpacity, side: THREE.DoubleSide, depthWrite: false, roughness: 0.75, metalness: 0 });
      this.faceMesh = new THREE.Mesh(geo, mat);
      this.faceMesh.renderOrder = 1;
      this.group.add(this.faceMesh);
    } else this.faceMesh = null;
    // Marked faces, over the others.
    if (marks.faces.size) {
      const mp = [], mi = [];
      for (const [i] of marks.faces) {
        const f = g.faces[i], base = mp.length / 3;
        mp.push(...f.positions);
        for (const j of f.indices) mi.push(base + j);
      }
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.Float32BufferAttribute(mp, 3));
      geo.setIndex(mi);
      geo.computeVertexNormals();
      const m = new THREE.Mesh(geo, new THREE.MeshStandardMaterial({ color: MARK, transparent: true, opacity: Math.min(0.85, st.faceOpacity + 0.3), side: THREE.DoubleSide, depthWrite: false, roughness: 0.6 }));
      m.renderOrder = 2;
      this.group.add(m);
    }

    // Edges, as wide lines.
    const segs = [], ecol = [];
    this.edgeOfSegment = [];
    g.edges.forEach((e, i) => {
      const colour = fields ? new THREE.Color(...cyclicColour(fields.link_angle[i])) : EDGE;
      for (let s = 0; s + 1 < e.points.length; s++) {
        segs.push(...e.points[s], ...e.points[s + 1]);
        ecol.push(colour.r, colour.g, colour.b, colour.r, colour.g, colour.b);
        this.edgeOfSegment.push(i);
      }
    });
    this.edgeMaterial = new LineMaterial({ linewidth: st.edgeWidth, vertexColors: true, transparent: true, opacity: 0.95, worldUnits: false });
    this.edgeMaterial.resolution.set(this.host.clientWidth, this.host.clientHeight);
    if (segs.length) {
      const geo = new LineSegmentsGeometry();
      geo.setPositions(segs);
      geo.setColors(ecol);
      this.edgeLines = new LineSegments2(geo, this.edgeMaterial);
      this.edgeLines.renderOrder = 3;
      this.group.add(this.edgeLines);
    } else this.edgeLines = null;
    // Marked edges, wider, gold for a positive coefficient and lilac for a negative one.
    if (marks.edges.size) {
      const ms = [], mc = [];
      for (const [i, k] of marks.edges) {
        const pts = g.edges[i].points, colour = k < 0 ? MARK_NEG : MARK;
        for (let s = 0; s + 1 < pts.length; s++) { ms.push(...pts[s], ...pts[s + 1]); mc.push(colour.r, colour.g, colour.b, colour.r, colour.g, colour.b); }
      }
      const geo = new LineSegmentsGeometry();
      geo.setPositions(ms);
      geo.setColors(mc);
      this.markMaterial = new LineMaterial({ linewidth: st.edgeWidth * 2.4 + 1, vertexColors: true, worldUnits: false });
      this.markMaterial.resolution.set(this.host.clientWidth, this.host.clientHeight);
      const m = new LineSegments2(geo, this.markMaterial);
      m.renderOrder = 4;
      this.group.add(m);
    } else this.markMaterial = null;

    // Arrows at the middle of each edge, pointing from tail to head (against it for a marked
    // edge whose coefficient is negative, so that a loop reads in its own direction).
    if (st.arrows && g.edges.length) {
      const cone = new THREE.ConeGeometry(0.32, 1, 12);
      cone.translate(0, 0, 0);
      const mesh = new THREE.InstancedMesh(cone, new THREE.MeshStandardMaterial({ roughness: 0.5 }), g.edges.length);
      const m4 = new THREE.Matrix4(), q = new THREE.Quaternion(), up = new THREE.Vector3(0, 1, 0), s = unit * 0.24 * st.arrowSize;
      g.edges.forEach((e, i) => {
        const k = marks.edges.get(i);
        const dir = new THREE.Vector3(...e.arrow.dir).multiplyScalar(k < 0 ? -1 : 1);
        q.setFromUnitVectors(up, dir.normalize());
        m4.compose(new THREE.Vector3(...e.arrow.at), q, new THREE.Vector3(s, s, s));
        mesh.setMatrixAt(i, m4);
        const colour = k !== undefined ? (k < 0 ? MARK_NEG : MARK) : fields ? new THREE.Color(...cyclicColour(fields.link_angle[i])) : EDGE;
        mesh.setColorAt(i, colour);
      });
      mesh.renderOrder = 5;
      this.group.add(mesh);
    }

    // Vertices.
    if (g.vertices.length) {
      const r = unit * 0.075 * st.vertexSize;
      const sphere = new THREE.SphereGeometry(1, 16, 12);
      this.vertexMesh = new THREE.InstancedMesh(sphere, new THREE.MeshStandardMaterial({ roughness: 0.45 }), g.vertices.length);
      const m4 = new THREE.Matrix4();
      g.vertices.forEach((p, i) => {
        const big = marks.vertices.has(i) ? 1.8 : 1;
        m4.compose(new THREE.Vector3(...p), new THREE.Quaternion(), new THREE.Vector3(r * big, r * big, r * big));
        this.vertexMesh.setMatrixAt(i, m4);
        this.vertexMesh.setColorAt(i, marks.vertices.has(i) ? MARK : VERTEX);
      });
      this.vertexMesh.renderOrder = 6;
      this.group.add(this.vertexMesh);
    } else this.vertexMesh = null;

    // Labels, where there are few enough to read.
    if (st.labels) {
      const names = this.names || [];
      const label = (text, at, cls) => {
        const div = document.createElement('div');
        div.className = `label3d ${cls}`;
        div.textContent = text;
        const o = new CSS2DObject(div);
        o.position.set(...at);
        this.group.add(o);
      };
      if (g.edges.length <= 80) g.edges.forEach((e, i) => {
        const off = new THREE.Vector3(...e.arrow.at).normalize().multiplyScalar(unit * 0.25);
        label(names[i] ?? `e${i + 1}`, [e.arrow.at[0] + off.x, e.arrow.at[1] + off.y, e.arrow.at[2] + off.z], 'edge');
      });
      if (g.vertices.length <= 60) g.vertices.forEach((p, i) => label(`v${i}`, [p[0], p[1] + unit * 0.22, p[2]], 'vertex'));
      if (g.faces.length <= 40) g.faces.forEach((f, i) => label(`f${i + 1}`, f.centre, 'face'));
    }
  }

  // What lies under the pointer, reported to onHover once per change.
  pick() {
    if (!this.pointer || !this.onHover) return;
    this.raycaster.setFromCamera(this.pointer, this.camera);
    this.raycaster.params.Line2 = { threshold: 6 };
    const targets = [this.vertexMesh, this.edgeLines, this.faceMesh].filter(Boolean);
    const hit = this.raycaster.intersectObjects(targets, false)[0];
    let what = null;
    if (hit) {
      if (hit.object === this.vertexMesh) what = { dim: 0, index: hit.instanceId };
      else if (hit.object === this.edgeLines) what = { dim: 1, index: this.edgeOfSegment[hit.faceIndex] };
      else if (hit.object === this.faceMesh) what = { dim: 2, index: this.faceOfTriangle[hit.faceIndex] };
    }
    const key = what ? `${what.dim}:${what.index}` : '';
    if (key !== this.hovered) { this.hovered = key; this.onHover(what); }
  }

  // A picture of the view, as a data URL.
  snapshot() { return this.renderer.domElement.toDataURL('image/png'); }
}
