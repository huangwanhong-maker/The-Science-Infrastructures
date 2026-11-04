// Fabrica's Reconstruct mode (ADR 0070): vessels reconstructed from images by the medical-imaging domain. A case
// is chosen on the left; the engine reconstructs it (once; it keeps the result), and the page draws the vessels
// in the round, the reference beside them as glass, the centrelines by radius, and the slices of the image with
// what was found and what the reference holds. The measures and the agreement go to the card.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

import { ask } from './data.js';
import { SYSTEM_STYLE } from './palette.js';

const $ = (id) => document.getElementById(id);
const VEIN = '#3d5a9e';
const COLOURS = { found: '#9aa3ad', arteries: SYSTEM_STYLE.arteries?.colour || '#cc4a52', veins: VEIN, reference: '#b98b74' };

const state = { cases: [], current: null, result: null, k: 0, slices: 1, busy: false, show: { found: true, reference: true, centrelines: false } };
let renderer, scene, camera, controls, group, started = false;

function css(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }

function startViewer() {
  if (started) return;
  started = true;
  const host = $('recon-view');
  renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(2, window.devicePixelRatio));
  host.appendChild(renderer.domElement);
  scene = new THREE.Scene();
  camera = new THREE.PerspectiveCamera(35, 1, 1, 5000);
  camera.up.set(0, 0, 1);
  scene.add(new THREE.HemisphereLight(0xffffff, 0x8a8f96, 1.6));
  const sun = new THREE.DirectionalLight(0xffffff, 1.4);
  sun.position.set(-1, -2, 3);
  scene.add(sun);
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  group = new THREE.Group();
  scene.add(group);
  new ResizeObserver(resize).observe(host);
  resize();
  renderer.setAnimationLoop(() => { controls.update(); renderer.render(scene, camera); });
}

function resize() {
  if (!renderer) return;
  const host = $('recon-view');
  const w = host.clientWidth || 1, h = host.clientHeight || 1;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}

// A surface from the engine (LPS millimetres) as a mesh: the page draws LPS as it comes, z up.
function mesh(surface, colour, { glass = false } = {}) {
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(surface.points, 3));
  g.setIndex(surface.triangles);
  g.computeVertexNormals();
  const m = new THREE.MeshStandardMaterial({ color: colour, roughness: 0.55, metalness: 0.05, transparent: glass, opacity: glass ? 0.28 : 1,
    depthWrite: !glass, side: glass ? THREE.DoubleSide : THREE.FrontSide });
  // LPS is right-handed (RAS turned half a turn about z), so the engine's outward triangles stay outward here.
  return new THREE.Mesh(g, m);
}

function centrelines(tree) {
  const out = new THREE.Group();
  const rmax = Math.max(...tree.segments.map((s) => s.radius), 1);
  for (const s of tree.segments) {
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(s.points, 3));
    const t = s.radius / rmax;
    const colour = new THREE.Color(css('--silver-deep') || '#7b838d').lerp(new THREE.Color(css('--rose') || '#b98b74'), t);
    out.add(new THREE.Line(g, new THREE.LineBasicMaterial({ color: colour })));
  }
  return out;
}

function draw() {
  if (!group) return;
  group.clear();
  const r = state.result;
  if (!r) return;
  const S = r.surfaces;
  if (state.show.found) {
    if (S.found) group.add(mesh(S.found, COLOURS.found));
    if (S.arteries) group.add(mesh(S.arteries, COLOURS.arteries));
    if (S.veins) group.add(mesh(S.veins, COLOURS.veins));
  }
  if (state.show.reference && S.reference) group.add(mesh(S.reference, COLOURS.reference, { glass: true }));
  if (state.show.centrelines) group.add(centrelines(r.tree));
}

function frame() {
  const box = new THREE.Box3().setFromObject(group);
  if (box.isEmpty()) return;
  const c = box.getCenter(new THREE.Vector3()), size = box.getSize(new THREE.Vector3()).length();
  controls.target.copy(c);
  // seen from the front (anterior is −y in LPS), a little from the right and above
  camera.position.copy(c).add(new THREE.Vector3(-0.35, -1.0, 0.35).normalize().multiplyScalar(1.7 * size));
  camera.near = size / 100, camera.far = size * 20;
  camera.updateProjectionMatrix();
}

const fmt = (x, d = 3) => (x === null || x === undefined || Number.isNaN(x) ? '—' : Number(x).toFixed(d));

function card() {
  const r = state.result, c = state.cases.find((x) => x.id === state.current);
  const el = $('recon-card');
  if (!r) { el.innerHTML = `<p class="lead">${c ? c.title : 'Choose a case.'}</p><p>${c ? c.note : ''}</p>`; return; }
  const a = r.agreement, t = r.tree, g = r.geometry;
  const rows = (pairs) => pairs.map(([k, v]) => `<tr><th>${k}</th><td>${v}</td></tr>`).join('');
  let html = `<p class="lead">${c.title}</p><p>${c.note}</p>
    <h3>Image</h3><table class="facts">${rows([['grid', `${g.size.join(' × ')}`], ['spacing', `${g.spacing.map((x) => fmt(x, 2)).join(' × ')} mm`]])}</table>
    <h3>Agreement with the ${r.kind === 'simulated' ? 'truth' : 'experts’ labels'}</h3>
    <table class="facts">${rows([['Dice', fmt(a.dice)], ['clDice', `${fmt(a.cldice)} <span class="note">(centreline precision ${fmt(a.precision)}, sensitivity ${fmt(a.sensitivity)})</span>`],
      ['surface distance', `${fmt(a.mean_distance, 2)} mm mean, ${fmt(a.p95_distance, 2)} mm at the 95th percentile`]])}</table>`;
  if (r.kind === 'simulated') {
    html += `<h3>Integration</h3><table class="facts">${rows([['registration', `${fmt(r.registration.error, 2)} mm over the CT’s box (NMI ${fmt(r.registration.similarity)})`],
      ['arteries', `Dice ${fmt(r.arteries.dice)}`], ['veins', `Dice ${fmt(r.veins.dice)}`]])}</table>`;
  }
  html += `<h3>The vascular tree</h3><table class="facts">${rows([['branches', t.segments.length], ['junctions', t.junctions], ['ends', t.ends],
    ['length', `${fmt(t.length, 0)} mm`], ['components', t.components], ['loops', t.cycles],
    ['Murray’s exponent', t.murray === null ? '— <span class="note">(not determined by these bifurcations)</span>' : `${fmt(t.murray, 2)} <span class="note">(fitted; 3 is Murray’s)</span>`]])}</table>
    <p class="credit">${r.credit}</p>`;
  el.innerHTML = html;
}

async function slice() {
  if (!state.result) return;
  const s = await ask('slice', { case: state.current, k: state.k });
  state.slices = s.slices;
  $('recon-k').max = String(s.slices - 1);
  $('recon-k').value = String(s.k);
  $('recon-where').textContent = `slice ${s.k + 1} of ${s.slices}, z = ${s.z} mm`;
  const W = s.width, H = s.height;
  const px = Uint8Array.from(atob(s.pixels), (ch) => ch.charCodeAt(0)), codes = Uint8Array.from(atob(s.codes), (ch) => ch.charCodeAt(0));
  const canvas = $('recon-canvas');
  canvas.width = W, canvas.height = H;
  const ctx = canvas.getContext('2d');
  const img = ctx.createImageData(W, H);
  // radiological display: the patient's right on the left; anterior up. LPS: +x left, +y posterior.
  const d = s.geometry.direction;
  const flipX = d[0][0] < 0, flipY = d[1][1] < 0;
  const tint = (hex) => [parseInt(hex.slice(1, 3), 16), parseInt(hex.slice(3, 5), 16), parseInt(hex.slice(5, 7), 16)];
  const cf = tint(COLOURS.found), cr = tint(COLOURS.reference), ca = tint(COLOURS.arteries), cv = tint(COLOURS.veins);
  for (let j = 0; j < H; j++)
    for (let i = 0; i < W; i++) {
      const src = i + W * j, x = flipX ? W - 1 - i : i, y = flipY ? H - 1 - j : j, o = 4 * (x + W * y);
      let r = px[src], g = px[src], b = px[src];
      const c = codes[src];
      const mix = (col, w) => { r = r * (1 - w) + col[0] * w; g = g * (1 - w) + col[1] * w; b = b * (1 - w) + col[2] * w; };
      if (c & 16) { r *= 0.55; g *= 0.55; b *= 0.6; }
      if (c & 4) mix(ca, 0.6); else if (c & 8) mix(cv, 0.6); else if (c & 1) mix(cf, 0.5);
      if ((c & 2) && !(c & 1) && !(c & 4) && !(c & 8)) mix(cr, 0.7);
      img.data[o] = r, img.data[o + 1] = g, img.data[o + 2] = b, img.data[o + 3] = 255;
    }
  ctx.putImageData(img, 0, 0);
}

function renderCases() {
  $('recon-cases').innerHTML = state.cases.map((c) => `<li><button class="case ${c.id === state.current ? 'on' : ''}" data-case="${c.id}">
    <span class="t">${c.title}</span><span class="k">${c.kind === 'simulated' ? 'simulated · CT + MRA' : 'real · CT'}</span></button></li>`).join('');
  for (const b of $('recon-cases').querySelectorAll('button')) b.addEventListener('click', () => choose(b.dataset.case));
}

async function choose(id) {
  if (state.busy) return;
  state.current = id;
  state.result = null;
  renderCases();
  card();
  draw();
  state.busy = true;
  $('recon-status').textContent = 'Reconstructing: vesselness across scales, the lumen at its half-maximum, centrelines and surfaces. The first time takes a minute or two; the engine keeps the result.';
  try {
    state.result = await ask('reconstruct', { case: id });
    state.k = Math.floor(state.result.geometry.size[2] / 2);
    $('recon-status').textContent = '';
    draw();
    frame();
    card();
    await slice();
  } catch (e) {
    $('recon-status').textContent = `The engine could not reconstruct this case: ${e.message || e}`;
  } finally {
    state.busy = false;
  }
}

export async function openReconstruct() {
  startViewer();
  resize();
  if (!state.cases.length) {
    try {
      state.cases = await ask('cases', {});
    } catch (e) {
      $('recon-status').textContent = `The engine is not ready: ${e.message || e}`;
      return;
    }
    renderCases();
    card();
  }
}

export function attachReconstruct() {
  for (const [key, id] of [['found', 'recon-show-found'], ['reference', 'recon-show-reference'], ['centrelines', 'recon-show-lines']]) {
    $(id).addEventListener('change', (e) => { state.show[key] = e.target.checked; draw(); });
  }
  $('recon-k').addEventListener('input', (e) => { state.k = Number(e.target.value); clearTimeout(state.timer); state.timer = setTimeout(slice, 60); });
}
