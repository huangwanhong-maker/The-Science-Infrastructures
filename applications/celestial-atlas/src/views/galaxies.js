// Galaxies in three dimensions, in supergalactic coordinates, with the relations
// between them:
//   Local Volume (tidal catalogue): each galaxy is linked to its main disturber,
//     the neighbour whose tidal force on it is greatest; galaxies whose tidal
//     index is positive form groups.
//   Local Universe (group catalogue): each galaxy is linked to the centre of the
//     group it belongs to.
// The supergalactic plane (SGZ = 0) is horizontal; the Milky Way is at the origin.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';
import { starPoints } from '../render/points.js';
import { ramp, TIDAL_STOPS, TYPE_STOPS, GROUP_STOPS } from '../render/colour.js';
import { declutter } from '../render/declutter.js';

const toScene = (x, y, z, out = new THREE.Vector3()) => out.set(x, z, -y);
const SCALES = {
  tidal: { rings: [1, 2, 5, 10], camera: [7, 5.5, 9], pick: 0.08, unit: 'Mpc' },
  groups: { rings: [25, 50, 100, 150, 200], camera: [150, 120, 190], pick: 1.5, unit: 'h⁻¹ Mpc' },
};
const MAX_LABELS = 48;

function label(className, html, onClick, centre = [0, 0.5]) {
  const el = document.createElement(onClick ? 'button' : 'div');
  el.className = className;
  el.innerHTML = html;
  if (onClick) {
    el.type = 'button';
    el.addEventListener('click', (event) => { event.stopPropagation(); onClick(); });
  }
  const object = new CSS2DObject(el);
  object.center.set(...centre);
  return object;
}

const escapeHtml = (s) => String(s).replace(/[&<>"']/g, (c) => `&#${c.charCodeAt(0)};`);

export class GalaxiesView {
  constructor({ renderer, container, onSelect }) {
    this.renderer = renderer;
    this.onSelect = onSelect;
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(50, 1, 0.001, 1e5);
    this.controls = new OrbitControls(this.camera, renderer.domElement);
    Object.assign(this.controls, { enableDamping: true, dampingFactor: 0.08, zoomSpeed: 1.2, minDistance: 0.02 });
    this.labels = new CSS2DRenderer();
    this.labels.domElement.className = 'label-layer';
    container.appendChild(this.labels.domElement);
    this.content = new THREE.Group();
    this.scene.add(this.content);
    this.raycaster = new THREE.Raycaster();
    this.options = { relations: 'published', field: false, colourBy: 'relation', minMembers: 3, labels: true, grid: true };
    this.cat = null;
    this.selected = null;
    this.tween = null;
    this.labelEntries = [];
    this.marker = label('star-marker', '', null, [0.5, 0.5]);
    this.marker.visible = false;
    this.scene.add(this.marker);
    this.selectionLabel = label('body-label selected', '', null);
    this.selectionLabel.visible = false;
    this.scene.add(this.selectionLabel);
  }

  setCatalogue(cat) {
    this.cat = cat;
    this.select(null);
    this.rebuild(true);
  }

  setOptions(options) {
    const before = JSON.stringify(this.options);
    Object.assign(this.options, options);
    if (this.cat && JSON.stringify(this.options) !== before) this.rebuild(false);
  }

  relationSet() {
    return this.cat.kind === 'tidal' ? this.cat.relations[this.options.relations] : null;
  }

  // Colour of galaxy i under the current options, as [r, g, b] in 0..1.
  colourOf(i) {
    const c = this.cat;
    if (this.options.colourBy === 'type') return ramp(TYPE_STOPS, c.type[i]);
    if (c.kind === 'tidal') return ramp(TIDAL_STOPS, this.relationSet().theta1[i]);
    return ramp(GROUP_STOPS, Math.log10(Math.max(1, c.groups.nHere[c.groupOf[i]])));
  }

  scenePosition(i, out = new THREE.Vector3()) {
    const p = this.cat.positions;
    return toScene(p[3 * i], p[3 * i + 1], p[3 * i + 2], out);
  }

  rebuild(resetCamera) {
    this.content.traverse((o) => {
      o.geometry?.dispose();
      o.material?.dispose();
      if (o.isCSS2DObject) o.element.remove();
    });
    this.content.clear();
    this.labelEntries = [];
    const c = this.cat, n = c.count, tidal = c.kind === 'tidal';
    const scale = SCALES[c.kind];

    // Galaxies.
    const positions = new Float32Array(3 * n), colours = new Float32Array(3 * n), sizes = new Float32Array(n);
    const v = new THREE.Vector3();
    for (let i = 0; i < n; i++) {
      this.scenePosition(i, v).toArray(positions, 3 * i);
      const L = c.logL[i];
      const size = tidal
        ? (Number.isFinite(L) ? Math.min(11, Math.max(2.4, 2.4 + (L - 6.5) * 1.6)) : 2.6)
        : (Number.isFinite(L) ? Math.min(5, Math.max(1.0, 1.0 + (L - 9.2) * 1.4)) : 1.0);
      const lone = !tidal && c.groups.nHere[c.groupOf[i]] < 2;
      // Additive points saturate where they crowd; the group map is dimmer.
      const k = tidal ? 1 : lone ? 0.3 : 0.5;
      const [r, g, b] = this.colourOf(i);
      colours.set([r * k, g * k, b * k], 3 * i);
      sizes[i] = size;
    }
    this.points = starPoints(positions, colours, sizes, this.renderer.getPixelRatio());
    this.content.add(this.points);

    // Relations.
    const segments = [], segmentColours = [];
    const a = new THREE.Vector3(), b = new THREE.Vector3();
    if (tidal) {
      const rel = this.relationSet();
      for (let i = 0; i < n; i++) {
        const d = rel.disturber[i], theta = rel.theta1[i];
        if (d < 0 || !Number.isFinite(theta) || (!this.options.field && theta <= 0)) continue;
        this.scenePosition(i, a); this.scenePosition(d, b);
        segments.push(a.x, a.y, a.z, b.x, b.y, b.z);
        const [r, g, bl] = ramp(TIDAL_STOPS, theta);
        const k = theta > 0 ? 0.42 : 0.16;
        segmentColours.push(r * k, g * k, bl * k, r * k * 0.45, g * k * 0.45, bl * k * 0.45);
      }
    } else {
      const gr = c.groups;
      for (let i = 0; i < n; i++) {
        const j = c.groupOf[i];
        if (gr.nHere[j] < Math.max(2, this.options.minMembers)) continue;
        this.scenePosition(i, a);
        toScene(gr.centre[3 * j], gr.centre[3 * j + 1], gr.centre[3 * j + 2], b);
        segments.push(a.x, a.y, a.z, b.x, b.y, b.z);
        const [r, g, bl] = ramp(GROUP_STOPS, Math.log10(gr.nHere[j]));
        segmentColours.push(r * 0.2, g * 0.2, bl * 0.2, r * 0.08, g * 0.08, bl * 0.08);
      }
    }
    this.relationCount = segments.length / 6;
    if (segments.length) {
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.Float32BufferAttribute(segments, 3));
      geometry.setAttribute('color', new THREE.Float32BufferAttribute(segmentColours, 3));
      const lines = new THREE.LineSegments(geometry, new THREE.LineBasicMaterial({
        vertexColors: true, transparent: true, depthWrite: false, blending: THREE.AdditiveBlending }));
      lines.frustumCulled = false;
      this.content.add(lines);
    }

    // Labels: groups, largest first.
    const labelGroup = new THREE.Group();
    if (tidal) {
      const rel = this.relationSet();
      for (const g of rel.groups.groups.filter((x) => x.size >= 3).slice(0, MAX_LABELS)) {
        const text = `${c.names[g.leader]} group`;
        const object = label('star-label', `<span class="name">${escapeHtml(text)}</span>`, () => this.onSelect(g.leader));
        this.scenePosition(g.leader, object.position);
        labelGroup.add(object);
        this.labelEntries.push({ object, text });
      }
    } else {
      const gr = c.groups;
      const ranked = [...gr.names.keys()].sort((x, y) =>
        (gr.names.get(y).endsWith('Cluster') - gr.names.get(x).endsWith('Cluster')) || gr.nHere[y] - gr.nHere[x]);
      for (const j of ranked.slice(0, MAX_LABELS)) {
        const text = gr.names.get(j);
        const target = gr.brightest[j] >= 0 ? gr.brightest[j] : null;
        const object = label('star-label', `<span class="name">${escapeHtml(text)}</span>`,
          target === null ? null : () => this.onSelect(target));
        toScene(gr.centre[3 * j], gr.centre[3 * j + 1], gr.centre[3 * j + 2], object.position);
        labelGroup.add(object);
        this.labelEntries.push({ object, text });
      }
      const home = label('body-label sun', '<span class="dot" style="background:#ffcf7a"></span><span class="name">Milky Way</span>');
      labelGroup.add(home);
      this.labelEntries.unshift({ object: home, text: 'Milky Way', always: true });
    }
    this.content.add(labelGroup);

    // Supergalactic plane: rings and the SGX and SGY axes.
    this.grid = new THREE.Group();
    this.ringLabels = [];
    const ringMaterial = new THREE.LineBasicMaterial({ color: 0xf1e2ee, transparent: true, opacity: 0.08 });
    for (const r of scale.rings) {
      const pts = [];
      for (let k = 0; k <= 160; k++) {
        const t = (2 * Math.PI * k) / 160;
        pts.push(new THREE.Vector3(r * Math.cos(t), 0, -r * Math.sin(t)));
      }
      this.grid.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), ringMaterial));
      const text = `${r} ${scale.unit}`;
      const tag = label('grid-label', text);
      tag.position.set(r * Math.cos(-2.4), 0, -r * Math.sin(-2.4));
      this.grid.add(tag);
      this.ringLabels.push({ object: tag, text });
    }
    const reach = scale.rings[scale.rings.length - 1] * 1.08;
    for (const [dir, text] of [[[reach, 0, 0], 'SGX'], [[0, 0, -reach], 'SGY'], [[0, reach * 0.4, 0], 'SGZ']]) {
      this.grid.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(
        [new THREE.Vector3(0, 0, 0), new THREE.Vector3(...dir)]),
        new THREE.LineBasicMaterial({ color: 0xca8dbe, transparent: true, opacity: 0.3 })));
      const tag = label('grid-label accent', text);
      tag.position.set(...dir);
      this.grid.add(tag);
    }
    this.content.add(this.grid);

    if (resetCamera) {
      this.camera.position.set(...scale.camera);
      this.controls.target.set(0, 0, 0);
      this.controls.maxDistance = reach * 6;
      this.tween = null;
    }
    if (this.selected != null) this.select(this.selected);
  }

  select(index) {
    this.selected = index;
    if (index == null || !this.cat) {
      this.marker.visible = false;
      this.selectionLabel.visible = false;
      return;
    }
    const name = this.cat.kind === 'tidal' ? this.cat.names[index] : this.cat.galaxyName(index);
    this.scenePosition(index, this.marker.position);
    this.selectionLabel.element.innerHTML = `<span class="name">${escapeHtml(name)}</span>`;
    this.selectionLabel.position.copy(this.marker.position);
    this.selectionLabel.center.set(-0.15, 0.5);
    this.marker.visible = true;
    this.selectionLabel.visible = true;
  }

  // Flies the camera to look at a point, keeping its direction.
  flyTo(target, distance) {
    this.tween = {
      from: this.controls.target.clone(), to: target.clone(),
      distFrom: this.camera.position.distanceTo(this.controls.target), distTo: distance, t: 0,
    };
  }

  focusGalaxy(index) {
    const radius = this.cat.kind === 'tidal' ? 1.2 : 12;
    this.flyTo(this.scenePosition(index), radius);
  }

  focusGroupCentre(j) {
    const gr = this.cat.groups;
    const r = Number.isFinite(gr.r2t[j]) && gr.r2t[j] > 0 ? gr.r2t[j] : 1;
    this.flyTo(toScene(gr.centre[3 * j], gr.centre[3 * j + 1], gr.centre[3 * j + 2]), Math.max(8, r * 9));
  }

  pick(ndc) {
    if (!this.points) return null;
    const distance = this.camera.position.distanceTo(this.controls.target);
    this.raycaster.params.Points.threshold = Math.max(SCALES[this.cat.kind].pick * 0.2, distance * 0.012);
    this.raycaster.setFromCamera(ndc, this.camera);
    const hits = this.raycaster.intersectObject(this.points, false);
    if (!hits.length) return null;
    hits.sort((x, y) => x.distanceToRay - y.distanceToRay);
    return hits[0].index;
  }

  update(snap, dt) {
    if (!this.cat) return;
    if (this.tween) {
      this.tween.t = Math.min(1, this.tween.t + dt / 0.9);
      const t = this.tween.t, k = t < 0.5 ? 4 * t ** 3 : 1 - (-2 * t + 2) ** 3 / 2;
      const direction = this.camera.position.clone().sub(this.controls.target).normalize();
      this.controls.target.lerpVectors(this.tween.from, this.tween.to, k);
      this.camera.position.copy(this.controls.target)
        .addScaledVector(direction, this.tween.distFrom + (this.tween.distTo - this.tween.distFrom) * k);
      if (t >= 1) this.tween = null;
    }
    this.grid.visible = this.options.grid;
    this.controls.update();
    this.scene.updateMatrixWorld();
    declutter([
      ...(this.selected != null ? [{ object: this.selectionLabel, text: this.selectionLabel.element.textContent,
        wanted: true, always: true }] : []),
      ...this.labelEntries.map((e) => ({ object: e.object, text: e.text, wanted: e.always || this.options.labels,
        always: e.always })),
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
    const c = this.cat;
    if (!c) return '';
    if (c.kind === 'tidal') {
      const rel = this.relationSet();
      const groups = rel.groups.groups.filter((g) => g.size >= 3).length;
      return `${c.count} galaxies · ${groups} tidal groups of three or more · `
        + `${this.relationCount} relations drawn · ${rel.label}`;
    }
    return `${c.count.toLocaleString('en')} galaxies in ${c.groups.count.toLocaleString('en')} groups · `
      + `distances from velocities, ${c.units}`;
  }
}
