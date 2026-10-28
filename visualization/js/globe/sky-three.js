// The sky behind the globes, drawn with three.js: the stars of a catalogue at their directions, and
// the Sun's disc with its glare.
//
// The world's axes are those of the frame the page draws in (the GCRS or the ICRS: the same axes),
// with the camera at the origin (renderer-three.js); a star is a direction, drawn far beyond
// everything else. Its brightness follows its magnitude, a star of magnitude m having the flux
// 10^(−0.4 m) relative to magnitude 0, shown as a point whose size and opacity grow with the flux
// (a perceptual mapping, not a photometric one); its colour is that of a blackbody at the
// temperature its colour index gives (stars.js). Aberration and parallax are not applied: at the
// resolution of a globe's view they are below a pixel.

import { starColour } from '../stars.js';

const STAR_VERTEX = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_vertex>
attribute float flux;
attribute vec3 colour;
uniform float scale;
uniform float limit;
varying vec3 vColour;
varying float vAlpha;
void main() {
  vColour = colour;
  float f = flux * scale;
  vAlpha = clamp(f, 0.0, 1.0);
  gl_PointSize = clamp(1.2 + 2.6 * sqrt(f), 1.0, 7.0);
  vec4 mv = modelViewMatrix * vec4(position, 1.0);
  gl_Position = projectionMatrix * mv;
  #include <logdepthbuf_vertex>
}
`;
const STAR_FRAGMENT = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_fragment>
varying vec3 vColour;
varying float vAlpha;
void main() {
  #include <logdepthbuf_fragment>
  vec2 d = gl_PointCoord - 0.5;
  float r = length(d) * 2.0;
  float a = vAlpha * smoothstep(1.0, 0.25, r);
  if (a < 0.01) discard;
  gl_FragColor = vec4(vColour, a);
}
`;

const SUN_FRAGMENT = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_fragment>
varying vec2 vUv;
uniform float disc; // the disc's radius as a fraction of the sprite's
void main() {
  #include <logdepthbuf_fragment>
  float r = length(vUv - 0.5) * 2.0;
  float core = smoothstep(disc * 1.02, disc * 0.96, r);
  float glare = pow(max(1.0 - r, 0.0), 4.0) * 0.55 + pow(max(1.0 - r, 0.0), 14.0) * 0.6;
  vec3 c = vec3(1.0, 0.97, 0.90) * core + vec3(1.0, 0.86, 0.66) * glare;
  float a = max(core, glare);
  if (a < 0.003) discard;
  gl_FragColor = vec4(c, a);
}
`;
const SUN_VERTEX = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_vertex>
varying vec2 vUv;
void main() {
  vUv = uv;
  vec4 mv = modelViewMatrix * vec4(position, 1.0);
  gl_Position = projectionMatrix * mv;
  #include <logdepthbuf_vertex>
}
`;

// The Milky Way: a map of the sky in galactic coordinates (longitude 0 at the map's centre and
// increasing to the left, latitude from the map's top), drawn on a sphere about the camera. Each
// fragment's direction in the world's axes (the ICRS's) is turned to galactic axes by the matrix
// `toGalactic` (rows: SOFA's iauIcrs2g) and the map is read there.
const MILKY_VERTEX = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_vertex>
varying vec3 vDir;
void main() {
  vDir = position;
  vec4 mv = modelViewMatrix * vec4(position, 1.0);
  gl_Position = projectionMatrix * mv;
  #include <logdepthbuf_vertex>
}
`;
const MILKY_FRAGMENT = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_fragment>
uniform sampler2D map;
uniform mat3 toGalactic;
uniform float brightness;
uniform float southUp;
varying vec3 vDir;
void main() {
  #include <logdepthbuf_fragment>
  vec3 g = normalize(toGalactic * normalize(vDir));
  float l = atan(g.y, g.x);
  float b = asin(clamp(g.z, -1.0, 1.0));
  float u = fract(0.5 - l / (2.0 * PI));
  float v = southUp > 0.5 ? 0.5 - b / PI : 0.5 + b / PI;
  vec3 c = texture2D(map, vec2(u, v)).rgb;
  gl_FragColor = vec4(c * c * brightness, 1.0);
}
`;

// stars: { unit: Float32Array of directions (3 per star), mag, ci, photometry, byBrightness? };
// limit: the faintest magnitude drawn; distance: how far the sky is drawn (m); milkyWay: an image (or
// ImageBitmap) of the galactic map, with `southUp` when the map's top is the south galactic pole (as
// Solar System Scope's is), and `toGalactic`, the ICRS → galactic matrix (rows).
export function createSky(THREE, { stars = null, limit = 6.5, distance = 1e14, milkyWay = null } = {}) {
  const group = new THREE.Group();
  group.matrixAutoUpdate = false;
  let milky = null;
  if (milkyWay?.image) {
    const texture = new THREE.Texture(milkyWay.image);
    texture.needsUpdate = true;
    texture.colorSpace = THREE.SRGBColorSpace;
    texture.wrapS = THREE.RepeatWrapping;
    const m = milkyWay.toGalactic;
    const toGalactic = new THREE.Matrix3().set(m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2]);
    milky = new THREE.Mesh(new THREE.SphereGeometry(distance * 0.9, 64, 32), new THREE.ShaderMaterial({
      vertexShader: MILKY_VERTEX, fragmentShader: MILKY_FRAGMENT, side: THREE.BackSide, depthWrite: false,
      uniforms: { map: { value: texture }, toGalactic: { value: toGalactic }, brightness: { value: 0.35 }, southUp: { value: milkyWay.southUp ? 1 : 0 } },
    }));
    milky.frustumCulled = false;
    milky.renderOrder = -3;
    group.add(milky);
  }
  let points = null;
  if (stars) {
    const chosen = [];
    for (let i = 0; i < stars.mag.length; i++) if (stars.mag[i] <= limit) chosen.push(i);
    const position = new Float32Array(chosen.length * 3);
    const colour = new Float32Array(chosen.length * 3);
    const flux = new Float32Array(chosen.length);
    chosen.forEach((i, k) => {
      for (let c = 0; c < 3; c++) position[3 * k + c] = stars.unit[3 * i + c] * distance;
      const [r, g, b] = starColour(stars.photometry || 'johnson', stars.ci[i]);
      colour.set([r, g, b], 3 * k);
      flux[k] = 10 ** (-0.4 * (stars.mag[i] - 1.0));
    });
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(position, 3));
    geometry.setAttribute('colour', new THREE.BufferAttribute(colour, 3));
    geometry.setAttribute('flux', new THREE.BufferAttribute(flux, 1));
    points = new THREE.Points(geometry, new THREE.ShaderMaterial({
      vertexShader: STAR_VERTEX, fragmentShader: STAR_FRAGMENT, transparent: true, depthWrite: false,
      blending: THREE.AdditiveBlending, uniforms: { scale: { value: 1 }, limit: { value: limit } },
    }));
    points.frustumCulled = false;
    points.renderOrder = -2;
    group.add(points);
  }
  const sun = new THREE.Mesh(new THREE.PlaneGeometry(1, 1), new THREE.ShaderMaterial({
    vertexShader: SUN_VERTEX, fragmentShader: SUN_FRAGMENT, transparent: true, depthWrite: false, side: THREE.DoubleSide,
    blending: THREE.AdditiveBlending, uniforms: { disc: { value: 0.1 } },
  }));
  sun.frustumCulled = false;
  sun.matrixAutoUpdate = false;
  sun.renderOrder = -1;
  group.add(sun);

  return {
    group,
    milky,
    count: points ? points.geometry.attributes.flux.count : 0,
    // sun: { direction: unit vector in world axes, radius: its angular radius (rad), camera };
    // brightness: 0..1, dimming the stars (by day, or near a bright globe).
    update({ sunDirection, sunRadius = 0.00465, camera, brightness = 1, sun: drawSun = true }) {
      if (points) points.material.uniforms.scale.value = brightness;
      if (milky) milky.material.uniforms.brightness.value = 0.35 * brightness;
      if (!drawSun) { sun.visible = false; return; }
      if (!sunDirection) { sun.visible = false; return; }
      sun.visible = true;
      // A square facing the camera, at a distance d, whose half side subtends 10 radii of the disc.
      const d = distance * 0.5;
      const half = Math.tan(sunRadius * 10) * d;
      const z = new THREE.Vector3(...sunDirection).normalize();
      const up = camera ? new THREE.Vector3(0, 1, 0).applyQuaternion(camera.quaternion) : new THREE.Vector3(0, 0, 1);
      const x = new THREE.Vector3().crossVectors(up, z).normalize();
      if (x.lengthSq() < 1e-12) x.set(1, 0, 0);
      const y = new THREE.Vector3().crossVectors(z, x);
      const c = z.clone().multiplyScalar(d);
      // Facing the camera at the origin: its normal is −z.
      sun.matrix.set(-x.x * 2 * half, y.x * 2 * half, -z.x, c.x, -x.y * 2 * half, y.y * 2 * half, -z.y, c.y, -x.z * 2 * half, y.z * 2 * half, -z.z, c.z, 0, 0, 0, 1);
      sun.matrixWorldNeedsUpdate = true;
      sun.material.uniforms.disc.value = 0.1;
    },
    dispose() {
      if (points) { points.geometry.dispose(); points.material.dispose(); }
      if (milky) { milky.geometry.dispose(); milky.material.uniforms.map.value.dispose(); milky.material.dispose(); }
      sun.geometry.dispose();
      sun.material.dispose();
    },
  };
}

// The stars of the data library's product (data/products/star_catalog.py: float32 records of
// ra, dec (rad), magnitude, colour index, ...), as directions for createSky.
export async function loadStarProduct(base, name = 'stars-v1', limit = 6.5) {
  const meta = await (await fetch(`${base}/products/${name}.json`)).json();
  const record = new Float32Array(await (await fetch(`${base}/products/${meta.file || `${name}.bin`}`)).arrayBuffer());
  const stride = meta.fields.length, n = meta.count;
  const unit = [], mag = [], ci = [];
  for (let i = 0; i < n; i++) {
    const o = i * stride;
    const m = record[o + 2];
    if (!(m <= limit)) continue;
    const ra = record[o], dec = record[o + 1], c = Math.cos(dec);
    unit.push(c * Math.cos(ra), c * Math.sin(ra), Math.sin(dec));
    mag.push(m);
    ci.push(record[o + 3]);
  }
  return { unit: Float32Array.from(unit), mag: Float32Array.from(mag), ci: Float32Array.from(ci), photometry: meta.photometry || 'johnson', title: meta.title, credit: meta.source };
}
