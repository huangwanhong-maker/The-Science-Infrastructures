// Molecules in three dimensions (three.js 0.170.0): atoms as spheres sized by covalent radius and
// coloured by element, bonds as cylinders (double bonds as two), labels, and, optionally, symmetry
// elements: rotation axes as lines labelled Cn, mirror planes as faint discs, the centre of inversion
// as a ring; and an arrow for a vector response (a dipole).

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';
import { css, h } from '../dom.js';

const v3 = (a) => new THREE.Vector3(a[0], a[1], a[2]);

function label(text, cls = '') {
  const el = h('span', { style: `font: 12px var(--serif); color: var(--ink); background: color-mix(in srgb, var(--surface) 72%, transparent); padding: 0 3px; border-radius: 3px; ${cls}` }, text);
  return new CSS2DObject(el);
}

// molecules: [{molecule, offset: [x,y,z]}]; o: {elements: pointGroup JSON, dipole: bool, mirror: {normal, point}}
export function view(container, molecules, o = {}) {
  const width = container.clientWidth || 600, height = container.clientHeight || 420;
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(36, width / height, 0.1, 200);
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setSize(width, height);
  container.append(renderer.domElement);
  const labels = new CSS2DRenderer();
  labels.setSize(width, height);
  Object.assign(labels.domElement.style, { position: 'absolute', left: '0', top: '0', pointerEvents: 'none' });
  container.append(labels.domElement);
  scene.add(new THREE.HemisphereLight(0xffffff, 0x9988aa, 1.6));
  const sun = new THREE.DirectionalLight(0xffffff, 1.5);
  sun.position.set(4, 6, 8);
  scene.add(sun);

  const group = new THREE.Group();
  scene.add(group);
  const ink = new THREE.Color(css('--ink-2'));
  let extent = 1;
  for (const { molecule: m, offset = [0, 0, 0], title } of molecules) {
    const off = v3(offset);
    const g = new THREE.Group();
    g.position.copy(off);
    m.atoms.forEach((a) => {
      const color = new THREE.Color(css(`--atom-${a.element}`) || '#999');
      const mesh = new THREE.Mesh(new THREE.SphereGeometry(0.22 + a.radius * 0.32, 32, 20),
        new THREE.MeshStandardMaterial({ color, roughness: 0.55, metalness: 0.0 }));
      mesh.position.copy(v3(a.x));
      g.add(mesh);
      if (a.element !== 'C' || o.labelCarbon) {
        const l = label(a.element + (a.charge > 0 ? '+' : a.charge < 0 ? '−' : ''));
        l.position.copy(v3(a.x)).add(new THREE.Vector3(0, 0.42 + a.radius * 0.3, 0));
        g.add(l);
      }
      extent = Math.max(extent, v3(a.x).add(off).length());
    });
    for (const [i, j, order] of m.bonds) {
      const A = v3(m.atoms[i].x), B = v3(m.atoms[j].x);
      const d = B.clone().sub(A), len = d.length();
      const side = new THREE.Vector3().crossVectors(d, Math.abs(d.z) < 0.9 * len ? new THREE.Vector3(0, 0, 1) : new THREE.Vector3(1, 0, 0)).normalize().multiplyScalar(0.07);
      const shifts = order === 2 ? [side, side.clone().negate()] : [new THREE.Vector3()];
      for (const s of shifts) {
        const cyl = new THREE.Mesh(new THREE.CylinderGeometry(order === 2 ? 0.045 : 0.07, order === 2 ? 0.045 : 0.07, len, 14),
          new THREE.MeshStandardMaterial({ color: ink, roughness: 0.7 }));
        cyl.position.copy(A.clone().add(B).multiplyScalar(0.5).add(s));
        cyl.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), d.clone().normalize());
        g.add(cyl);
      }
    }
    if (title) {
      const t = label(title, 'font-style: italic; font-size: 13px;');
      t.position.set(0, -2.1, 0);
      g.add(t);
    }
    group.add(g);
  }
  // Symmetry elements of the first molecule, about its centroid.
  const pg = o.elements;
  if (pg && pg.elements) {
    const c = v3(pg.centroid || [0, 0, 0]).add(v3((molecules[0].offset) || [0, 0, 0]));
    const accent = new THREE.Color(css('--accent')), minor = new THREE.Color(css('--minor'));
    const L = extent * 1.25;
    for (const e of pg.elements) {
      const axis = v3(e.axis).normalize();
      if (e.kind === 'rotation' && o.show?.axes !== false) {
        const geo = new THREE.BufferGeometry().setFromPoints([c.clone().addScaledVector(axis, -L), c.clone().addScaledVector(axis, L)]);
        group.add(new THREE.Line(geo, new THREE.LineBasicMaterial({ color: accent })));
        const t = label(`C${e.order}`, `color: var(--accent-deep);`);
        t.position.copy(c.clone().addScaledVector(axis, L + 0.25));
        group.add(t);
      } else if (e.kind === 'mirror' && o.show?.mirrors !== false) {
        const disc = new THREE.Mesh(new THREE.CircleGeometry(extent * 1.1, 48),
          new THREE.MeshBasicMaterial({ color: minor, transparent: true, opacity: 0.12, side: THREE.DoubleSide, depthWrite: false }));
        disc.position.copy(c);
        disc.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), axis);
        group.add(disc);
      } else if (e.kind === 'inversion' && o.show?.inversion !== false) {
        const ring = new THREE.Mesh(new THREE.TorusGeometry(0.16, 0.035, 10, 32), new THREE.MeshStandardMaterial({ color: minor }));
        ring.position.copy(c);
        group.add(ring);
      } else if (e.kind === 'rotoreflection' && o.show?.axes !== false) {
        const t = label(`S${e.order}`, `color: var(--minor);`);
        t.position.copy(c.clone().addScaledVector(axis, -L - 0.25));
        group.add(t);
      }
    }
  }
  if (o.mirror) {
    const geo = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(0, -extent * 0.6, 0), new THREE.Vector3(0, extent * 0.6, 0)]);
    const line = new THREE.Line(geo, new THREE.LineDashedMaterial({ color: new THREE.Color(css('--minor')), dashSize: 0.12, gapSize: 0.08 }));
    line.computeLineDistances();
    group.add(line);
  }
  if (o.dipole) {
    const arrow = new THREE.ArrowHelper(v3(o.dipole.direction).normalize(), v3(o.dipole.origin || [0, 0, 0]), o.dipole.length || 1.4, new THREE.Color(css('--minor')), 0.28, 0.16);
    group.add(arrow);
  }

  // An oblique view, so that a planar molecule does not hide its axes and planes edge-on.
  camera.position.copy(v3(o.view || [0.42, 0.34, 1]).normalize().multiplyScalar(extent * 3.6));
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  let alive = true;
  const tick = () => {
    if (!alive) return;
    controls.update();
    renderer.render(scene, camera);
    labels.render(scene, camera);
    requestAnimationFrame(tick);
  };
  tick();
  const resize = () => {
    const w = container.clientWidth, hh = container.clientHeight;
    if (!w || !hh) return;
    camera.aspect = w / hh;
    camera.updateProjectionMatrix();
    renderer.setSize(w, hh);
    labels.setSize(w, hh);
  };
  window.addEventListener('resize', resize);
  return () => {
    alive = false;
    window.removeEventListener('resize', resize);
    controls.dispose();
    renderer.dispose();
  };
}
