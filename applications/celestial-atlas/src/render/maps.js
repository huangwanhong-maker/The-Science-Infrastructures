// Maps of the Sun, the Moon, the planets and the sky, from the imagery products
// of the data library (data/products/body_maps.py and imagery.py), as textures,
// and the orientations that put them in place.
//
// A body's map is in east longitude with 0 at its centre and north at the top.
// three.js wraps a texture round a sphere with the image's centre on the
// sphere's +x axis and its top at +y, the image's right-hand side toward -z.
// The body-fixed frame (x to the prime meridian, z to the north pole) is
// therefore (x, y, z)_body = (x, -z, y)_local.
//
// The sky map is in galactic coordinates with longitude 0 at the centre,
// increasing to the left, and the south galactic pole at the top; its local
// frame is (X, -Z, Y)_galactic, by the same reasoning.

import * as THREE from 'three';
import { multiply, transpose } from '../core/frames.js';

// Row-major 3 x 3 matrices.
export const BODY_FROM_LOCAL = [1, 0, 0, 0, 0, -1, 0, 1, 0];
export const LOCAL_FROM_GALACTIC = [1, 0, 0, 0, 0, -1, 0, 1, 0];

// A rotation as a three.js quaternion, from a row-major 3 x 3 matrix.
export function quaternionOf(m, out = new THREE.Quaternion()) {
  const m4 = new THREE.Matrix4().set(m[0], m[1], m[2], 0, m[3], m[4], m[5], 0, m[6], m[7], m[8], 0, 0, 0, 0, 1);
  return out.setFromRotationMatrix(m4);
}

// The rotation that puts the sky map in a scene whose axes are `sceneFromGalactic` of galactic.
export function skyRotation(sceneFromGalactic) {
  return quaternionOf(multiply(sceneFromGalactic, transpose(LOCAL_FROM_GALACTIC)));
}

async function json(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
  return r.json();
}

// Loads the maps that exist: { bodies: {key: texture}, rings, sky, earthDay,
// earthNight, live: {key: tile service for detail}, credits }.
export async function loadMaps(base, renderer) {
  let index;
  try { index = await json(`${base}/products/index.json`); } catch { return null; }
  const loader = new THREE.TextureLoader();
  const anisotropy = Math.min(8, renderer.capabilities.getMaxAnisotropy());
  const credits = new Set();
  const maps = { bodies: {}, rings: null, sky: null, earthDay: null, earthNight: null, live: {}, credits };
  const load = async (entry) => {
    try {
      const meta = await json(`${base}/products/${entry.file}`);
      if (meta.live && meta.body) maps.live[meta.body] = meta.live;
      const texture = await loader.loadAsync(`${base}/products/${meta.file}`);
      texture.colorSpace = THREE.SRGBColorSpace;
      texture.anisotropy = anisotropy;
      if (meta.source && meta.source.credit) credits.add(meta.source.credit);
      return texture;
    } catch {
      return null;
    }
  };
  await Promise.all((index.imagery || []).map(async (e) => {
    if (e.role === 'body') maps.bodies[e.body] = await load(e);
    else if (e.role === 'rings') maps.rings = await load(e);
    else if (e.role === 'sky') maps.sky = await load(e);
    else if (e.role === 'day' && e.name === 'blue-marble') maps.earthDay = await load(e);
    else if (e.role === 'night-lights') maps.earthNight = await load(e);
  }));
  return maps;
}

// A sphere seen from inside, carrying the sky map.
export function skySphere(texture, radius, brightness = 0.85) {
  const mesh = new THREE.Mesh(new THREE.SphereGeometry(radius, 96, 48),
    new THREE.MeshBasicMaterial({ map: texture, side: THREE.BackSide, depthWrite: false, color: new THREE.Color(brightness, brightness, brightness) }));
  mesh.renderOrder = -10;
  mesh.frustumCulled = false;
  return mesh;
}

// Lights on the night side of a body: the emissive map is shown only where the
// Sun is below the horizon (uSunView: unit vector toward the Sun, view space).
export function nightLights(material) {
  material.userData.sunView = { value: new THREE.Vector3(1, 0, 0) };
  material.onBeforeCompile = (shader) => {
    shader.uniforms.uSunView = material.userData.sunView;
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', '#include <common>\nuniform vec3 uSunView;')
      .replace('#include <emissivemap_fragment>',
        '#include <emissivemap_fragment>\n  totalEmissiveRadiance *= smoothstep(0.10, -0.12, dot(normal, uSunView));');
  };
  material.needsUpdate = true;
}
