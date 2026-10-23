// A quantum tetrahedron in space: the tetrahedron dual to a spin-foam vertex, its six edges of
// lengths j + ½ (Ponzano and Regge 1968), set in space from those lengths (embedTetrahedron in
// visualization/js/complex.js), its faces translucent, its edges coloured and labelled by spin.
// When the lengths admit no Euclidean tetrahedron it is drawn flattened, and says so.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';
import { embedTetrahedron } from '../../../../visualization/js/index.js';

const PAIRS = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];

export class TetrahedronView {
  constructor(host) {
    this.host = host;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    host.append(this.renderer.domElement);
    this.labels = new CSS2DRenderer();
    this.labels.domElement.className = 'labels';
    host.append(this.labels.domElement);
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(34, 1, 0.01, 100);
    this.camera.position.set(0, 0.6, 3.2);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.autoRotate = true;
    this.controls.autoRotateSpeed = 1.2;
    this.scene.add(new THREE.HemisphereLight(0xf6e7f0, 0x2a1f2f, 1.6));
    const key = new THREE.DirectionalLight(0xfff1e6, 1.2);
    key.position.set(2, 3, 4);
    this.camera.add(key);
    this.scene.add(this.camera);
    this.group = new THREE.Group();
    this.scene.add(this.group);
    new ResizeObserver(() => this.resize()).observe(host);
    this.resize();
    const loop = () => {
      this.frame = requestAnimationFrame(loop);
      this.controls.update();
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
  }

  // lengths [l01, l02, l03, l12, l13, l23]; colours [[r, g, b]] per edge; labels per edge.
  show(lengths, { colours, labels }) {
    for (const o of [...this.group.children]) {
      o.traverse((x) => { x.geometry?.dispose?.(); x.material?.dispose?.(); if (x.isCSS2DObject) x.element.remove(); });
      this.group.remove(o);
    }
    const T = embedTetrahedron(lengths);
    const scale = 1.6 / Math.max(...lengths);
    const P = T.positions.map((p) => new THREE.Vector3(...p).multiplyScalar(scale));
    const faces = new THREE.BufferGeometry().setFromPoints([0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3].map((i) => P[i]));
    faces.computeVertexNormals();
    this.group.add(new THREE.Mesh(faces, new THREE.MeshStandardMaterial({ color: '#c08081', transparent: true, opacity: 0.28, side: THREE.DoubleSide, depthWrite: false, roughness: 0.7 })));
    PAIRS.forEach(([a, b], k) => {
      const from = P[a], to = P[b], dir = new THREE.Vector3().subVectors(to, from), len = dir.length();
      const tube = new THREE.Mesh(new THREE.CylinderGeometry(0.022, 0.022, len, 12), new THREE.MeshStandardMaterial({ color: new THREE.Color(...colours[k]), roughness: 0.4 }));
      tube.position.copy(from).addScaledVector(dir, 0.5);
      tube.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.clone().normalize());
      this.group.add(tube);
      const div = document.createElement('div');
      div.className = 'label3d spin';
      div.textContent = labels[k];
      const tag = new CSS2DObject(div);
      tag.position.copy(from).addScaledVector(dir, 0.5).multiplyScalar(1.12);
      this.group.add(tag);
    });
    for (const p of P) {
      const s = new THREE.Mesh(new THREE.SphereGeometry(0.05, 16, 12), new THREE.MeshStandardMaterial({ color: '#f2c27c', roughness: 0.45 }));
      s.position.copy(p);
      this.group.add(s);
    }
    return T;
  }

  dispose() {
    cancelAnimationFrame(this.frame);
    this.renderer.dispose();
    this.host.replaceChildren();
  }
}
