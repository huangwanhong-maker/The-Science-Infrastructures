// The Sun's stellar neighbourhood in three dimensions, from catalogue distances
// (Hipparcos through HYG, or Gaia), in galactic coordinates (parsecs): x towards
// the galactic centre, y in the direction of galactic rotation, z towards the
// north galactic pole.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';
import { ICRS_TO_GALACTIC, apply } from '../core/frames.js';
import { starPoints, magnitudeStyle } from '../render/points.js';
import { starColour } from '../render/colour.js';
import { declutter } from '../render/declutter.js';
import { skyRotation, skySphere } from '../render/maps.js';

export const LY_PER_PC = 3.261563777;
const RING_RADII = [5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000];
const MAX_NAMED_LABELS = 32;
const NEAREST_LABELS = 6;

const toScene = (g, out = new THREE.Vector3()) => out.set(g[0], g[2], -g[1]);
const SCENE_FROM_GALACTIC = [1, 0, 0, 0, 0, 1, 0, -1, 0];

function label(className, html, onClick) {
  const el = document.createElement(onClick ? 'button' : 'div');
  el.className = className;
  el.innerHTML = html;
  if (onClick) {
    el.type = 'button';
    el.addEventListener('click', (event) => { event.stopPropagation(); onClick(); });
  }
  const object = new CSS2DObject(el);
  object.center.set(0, 0.5);
  return object;
}

export function absoluteMagnitude(mag, distPc) {
  return mag - 5 * Math.log10(distPc / 10);
}

export class NeighbourhoodView {
  constructor({ renderer, container, stars, onSelect }) {
    this.renderer = renderer;
    this.stars = stars;
    this.onSelect = onSelect;
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(50, 1, 0.001, 1e5);
    this.controls = new OrbitControls(this.camera, renderer.domElement);
    Object.assign(this.controls, { enableDamping: true, dampingFactor: 0.08, zoomSpeed: 1.2, minDistance: 0.05 });
    this.labels = new CSS2DRenderer();
    this.labels.domElement.className = 'label-layer';
    container.appendChild(this.labels.domElement);
    this.options = { labels: true, grid: true };
    this.radius = null;
    this.content = new THREE.Group();
    this.scene.add(this.content);
    this.raycaster = new THREE.Raycaster();

    const sun = label('body-label sun', '<span class="dot" style="background:#ffcf7a"></span><span class="name">Sun</span>');
    this.scene.add(sun);
    this.sunLabel = sun;
    this.marker = label('star-marker', '');
    this.marker.center.set(0.5, 0.5);
    this.marker.visible = false;
    this.scene.add(this.marker);
    this.setRadius(20);
  }

  galactic(i) {
    const s = this.stars;
    return apply(ICRS_TO_GALACTIC, [s.xyz[3 * i], s.xyz[3 * i + 1], s.xyz[3 * i + 2]]);
  }

  setRadius(radius) {
    const previous = this.radius;
    this.radius = radius;
    this.rebuild();
    if (previous) {
      this.camera.position.multiplyScalar(radius / previous);
    } else {
      this.camera.position.set(radius * 1.25, radius * 0.75, radius * 1.45);
    }
    this.controls.maxDistance = radius * 8;
    this.controls.target.set(0, 0, 0);
  }

  setStars(stars) {
    this.stars = stars;
    this.select(null);
    this.rebuild();
  }

  rebuild() {
    const radius = this.radius;
    // Release GPU buffers, and remove label elements: three.js removes a
    // CSS2DObject's element only when that object itself is removed, not when
    // an ancestor group is cleared.
    this.content.traverse((o) => {
      o.geometry?.dispose();
      o.material?.dispose();
      if (o.isCSS2DObject) o.element.remove();
    });
    this.content.clear();
    const s = this.stars;
    const indices = [];
    for (let i = 0; i < s.count; i++) {
      if (s.dist[i] <= radius) indices.push(i);
    }
    this.indices = indices;

    const positions = new Float32Array(indices.length * 3);
    const colours = new Float32Array(indices.length * 3);
    const sizes = new Float32Array(indices.length);
    const v = new THREE.Vector3();
    // With many stars, additive blending saturates; dim and shrink the points
    // as the count grows, so that density still reads as structure.
    const density = Math.min(1, Math.sqrt(3000 / Math.max(indices.length, 1)));
    const dim = Math.max(0.12, density), shrink = Math.max(0.55, density ** 0.3);
    indices.forEach((i, k) => {
      toScene(this.galactic(i), v).toArray(positions, 3 * k);
      const { size, intensity } = magnitudeStyle(absoluteMagnitude(s.mag[i], s.dist[i]), 15);
      const [r, g, b] = starColour(s.photometry, s.ci[i]);
      const k3 = intensity * dim;
      colours.set([r * k3, g * k3, b * k3], 3 * k);
      sizes[k] = size * shrink;
    });
    this.points = starPoints(positions, colours, sizes, this.renderer.getPixelRatio());
    this.content.add(this.points);

    // Labels: the brightest named stars, and the nearest stars whatever their names.
    this.labelGroup = new THREE.Group();
    const named = indices.filter((i) => s.starName(i)).sort((a, b) => s.mag[a] - s.mag[b])
      .slice(0, MAX_NAMED_LABELS);
    const nearest = [...indices].sort((a, b) => s.dist[a] - s.dist[b]).slice(0, NEAREST_LABELS);
    const nearestSet = new Set(nearest);
    this.labelEntries = [];
    for (const i of new Set([...named, ...nearest])) {
      const text = s.label(i) || (nearestSet.has(i) ? s.designation(i) : '');
      if (!text) continue;
      const object = label('star-label', `<span class="name">${text}</span>`, () => this.onSelect(i));
      toScene(this.galactic(i), object.position);
      this.labelGroup.add(object);
      this.labelEntries.push({ object, text, i });
    }
    this.content.add(this.labelGroup);

    this.grid = new THREE.Group();
    this.ringLabels = [];
    const material = new THREE.LineBasicMaterial({ color: 0xf1e2ee, transparent: true, opacity: 0.09 });
    for (const r of RING_RADII.filter((x) => x <= radius * 1.05)) {
      const pts = [];
      for (let k = 0; k <= 160; k++) {
        const a = (2 * Math.PI * k) / 160;
        pts.push(new THREE.Vector3(r * Math.cos(a), 0, -r * Math.sin(a)));
      }
      this.grid.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), material));
      const text = `${r.toLocaleString('en')} pc · ${Math.round(r * LY_PER_PC).toLocaleString('en')} ly`;
      const tag = label('grid-label', text);
      tag.position.set(r * Math.cos(-2.2), 0, -r * Math.sin(-2.2));
      this.grid.add(tag);
      this.ringLabels.push({ object: tag, text });
    }
    const reach = radius * 1.1;
    const axis = new THREE.Line(new THREE.BufferGeometry().setFromPoints(
      [new THREE.Vector3(0, 0, 0), new THREE.Vector3(reach, 0, 0)]),
      new THREE.LineBasicMaterial({ color: 0xca8dbe, transparent: true, opacity: 0.35 }));
    this.grid.add(axis);
    const centre = label('grid-label accent', '→ galactic centre');
    centre.position.set(radius * 0.8, 0, 0);
    this.grid.add(centre);
    const pole = new THREE.Line(new THREE.BufferGeometry().setFromPoints(
      [new THREE.Vector3(0, -reach * 0.5, 0), new THREE.Vector3(0, reach * 0.5, 0)]),
      new THREE.LineBasicMaterial({ color: 0xf1e2ee, transparent: true, opacity: 0.12 }));
    this.grid.add(pole);
    const north = label('grid-label', 'north galactic pole');
    north.position.set(0, reach * 0.5, 0);
    this.grid.add(north);
    this.content.add(this.grid);
  }

  setOptions(options) {
    Object.assign(this.options, options);
  }

  select(index) {
    this.selected = index;
    if (index == null) { this.marker.visible = false; return; }
    toScene(this.galactic(index), this.marker.position);
    this.marker.visible = true;
  }

  pick(ndc) {
    this.raycaster.params.Points.threshold = this.radius * 0.012;
    this.raycaster.setFromCamera(ndc, this.camera);
    const hits = this.raycaster.intersectObject(this.points, false);
    if (!hits.length) return null;
    hits.sort((a, b) => a.distanceToRay - b.distanceToRay);
    return this.indices[hits[0].index];
  }

  // The Milky Way behind the stars, as seen from the Sun.
  setSky(texture) {
    if (!texture || this.skyMesh) return;
    this.skyMesh = skySphere(texture, 8e4, 0.7);
    this.skyMesh.quaternion.copy(skyRotation(SCENE_FROM_GALACTIC));
    this.scene.add(this.skyMesh);
  }

  update() {
    this.grid.visible = this.options.grid;
    if (this.skyMesh) {
      this.skyMesh.position.copy(this.camera.position);
      this.skyMesh.visible = this.options.milkyWay !== false;
    }
    this.controls.update();
    this.scene.updateMatrixWorld();
    declutter([
      { object: this.sunLabel, text: 'Sun', wanted: true, always: true },
      ...this.labelEntries.map((e) => ({
        object: e.object, text: e.text,
        wanted: this.options.labels || e.i === this.selected, always: e.i === this.selected,
      })),
      ...this.ringLabels.map((e) => ({ object: e.object, text: e.text, wanted: this.options.grid })),
    ], this.camera, this.width, this.height);
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
    this.renderer.render(this.scene, this.camera);
    this.labels.render(this.scene, this.camera);
  }

  legend() {
    return `${this.indices.length.toLocaleString('en')} stars within ${this.radius.toLocaleString('en')} pc `
      + `(${Math.round(this.radius * LY_PER_PC).toLocaleString('en')} light-years) · ${this.stars.title}`;
  }
}
