// Marks that live in time, on the graphics card: the moment moves by a uniform, so that scrubbing through a hundred thousand events costs nothing in the page. Each
// point or vertex carries its time on a normalised line (0 to 1 over a record's span, which single precision holds to a part in sixteen million), its value on
// [0, 1] read through a colour ramp, and a size. Geometry is relative to a centre the caller keeps, in whichever frame the marks are drawn: the Earth's relative
// to the eye, or a city's local frame.
//
//   blooms    events: each opens as a small point with a ring widening about it, and fades, like a stone in still water; what has passed can stay as dust
//   threads   tracks and flows: fine strips of a stated width on the screen, drawn up to the moment, brightest at the front, which runs between samples, trailing off behind
//
// three.js is passed in, as everywhere in this layer.

const LOG_V = `#include <common>
#include <logdepthbuf_pars_vertex>`;
const LOG_F = `#include <common>
#include <logdepthbuf_pars_fragment>`;

const BLOOM_V = /* glsl */ `${LOG_V}
attribute float aTime;
attribute float aValue;
attribute float aSize;
uniform float uNow;
uniform float uWindow;
uniform float uKeep;
uniform float uPixel;
varying float vAge;
varying float vValue;
void main() {
  float age = uNow - aTime;
  bool past = age > uWindow;
  if (age < 0.0 || (past && uKeep <= 0.0)) {
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
    gl_PointSize = 0.0;
    vAge = 2.0;
    vValue = 0.0;
    return;
  }
  vAge = past ? 2.0 : age / max(uWindow, 1e-9);
  vValue = aValue;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
  float base = uPixel * (0.45 + 1.55 * aSize);
  gl_PointSize = past ? max(1.5, base * 0.55) : base * (1.0 + 3.2 * sqrt(vAge));
  #include <logdepthbuf_vertex>
}
`;
const BLOOM_F = /* glsl */ `${LOG_F}
uniform sampler2D uRamp;
uniform float uOpacity;
uniform float uKeep;
varying float vAge;
varying float vValue;
void main() {
  vec2 p = gl_PointCoord * 2.0 - 1.0;
  float r = length(p);
  if (r > 1.0) discard;
  vec3 c = texture2D(uRamp, vec2(vValue, 0.5)).rgb;
  float a;
  if (vAge >= 1.5) {
    a = (1.0 - smoothstep(0.55, 1.0, r)) * uKeep;
  } else {
    float e = fwidth(r) * 1.3;
    float life = 1.0 - vAge;
    float core = (1.0 - smoothstep(0.0, 0.34 / (1.0 + 3.2 * sqrt(vAge)) + e, r)) * mix(0.55, 1.0, life);
    float ring = (1.0 - smoothstep(0.0, 0.05 + e, abs(r - 0.93))) * life * life * 0.8;
    a = max(core, ring);
    c = mix(c, vec3(1.0), 0.32 * life * life * core);
  }
  if (a <= 0.004) discard;
  gl_FragColor = vec4(c, a * uOpacity);
  #include <logdepthbuf_fragment>
  #include <colorspace_fragment>
}
`;

// A thread is a strip of a stated width on the screen, since the line the graphics card draws is one device pixel, which a fine display can hardly show. Each end of a segment
// is two vertices, either side of the line through it and the segment's other end. The width is in device pixels; a soft edge is drawn in its outer half.
const THREAD_V = /* glsl */ `${LOG_V}
attribute vec3 aOther;
attribute float aAcross;
attribute float aFlip;
attribute float aTime;
attribute float aValue;
uniform vec2 uViewport;
uniform float uWidth;
varying float vTime;
varying float vValue;
varying float vEdge;
void main() {
  vTime = aTime;
  vValue = aValue;
  vEdge = aAcross;
  vec4 here = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
  vec4 there = projectionMatrix * modelViewMatrix * vec4(aOther, 1.0);
  vec2 d = here.w > 0.0 && there.w > 0.0 ? (there.xy / there.w - here.xy / here.w) * uViewport : vec2(1.0, 0.0);
  float len = length(d);
  vec2 n = len > 1e-4 ? vec2(-d.y, d.x) / len : vec2(0.0, 1.0);
  here.xy += n * aFlip * aAcross * uWidth / uViewport * here.w;
  gl_Position = here;
  #include <logdepthbuf_vertex>
}
`;
const THREAD_F = /* glsl */ `${LOG_F}
uniform sampler2D uRamp;
uniform float uNow;
uniform float uWindow;
uniform float uKeep;
uniform float uOpacity;
varying float vTime;
varying float vValue;
varying float vEdge;
void main() {
  float age = uNow - vTime;
  if (age < 0.0) discard;
  float u = age / max(uWindow, 1e-9);
  float life = u < 1.0 ? (1.0 - u) : 0.0;
  float a = max(life * life, uKeep) * (1.0 - smoothstep(0.5, 1.0, abs(vEdge)));
  if (a <= 0.004) discard;
  vec3 c = texture2D(uRamp, vec2(vValue, 0.5)).rgb;
  c = mix(c, vec3(1.0), 0.7 * exp(-u * 26.0));
  gl_FragColor = vec4(c, a * uOpacity);
  #ifdef THROUGH
  gl_FragColor.a *= THROUGH;
  #endif
  #include <logdepthbuf_fragment>
  #include <colorspace_fragment>
}
`;

function uniforms(THREE, ramp, extra = {}) {
  return { uNow: { value: 0 }, uWindow: { value: 0.1 }, uKeep: { value: 0 }, uOpacity: { value: 1 }, uRamp: { value: ramp }, ...extra };
}
const handle = (object, material) => ({
  object,
  material,
  // The moment and how long a mark lives, on the normalised line; `keep` is how bright what has passed stays (0: gone).
  setTime(now, window, keep = material.uniforms.uKeep.value) {
    const u = material.uniforms;
    u.uNow.value = now, u.uWindow.value = window, u.uKeep.value = keep;
  },
  setOpacity(a) { material.uniforms.uOpacity.value = a; },
  setRamp(ramp) { material.uniforms.uRamp.value = ramp; },
  dispose() {
    object.geometry.dispose();
    material.dispose();
  },
});

// A ramp of 256 colours from a function u → [r, g, b] (sRGB, each 0 to 1).
export function rampTexture(THREE, colour) {
  const data = new Uint8Array(4 * 256);
  for (let i = 0; i < 256; i++) {
    const c = colour(i / 255);
    data[4 * i] = Math.round(c[0] * 255), data[4 * i + 1] = Math.round(c[1] * 255), data[4 * i + 2] = Math.round(c[2] * 255), data[4 * i + 3] = 255;
  }
  const t = new THREE.DataTexture(data, 256, 1, THREE.RGBAFormat);
  t.colorSpace = THREE.SRGBColorSpace;
  t.minFilter = t.magFilter = THREE.LinearFilter;
  t.needsUpdate = true;
  return t;
}

// Events: positions (x, y, z) and, for each, time, value and size.
export function createBlooms(THREE, { position, time, value, size, sphere }, { ramp, night = true, pixel = 7 } = {}) {
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(position, 3));
  g.setAttribute('aTime', new THREE.BufferAttribute(time, 1));
  g.setAttribute('aValue', new THREE.BufferAttribute(value, 1));
  g.setAttribute('aSize', new THREE.BufferAttribute(size, 1));
  if (sphere) g.boundingSphere = new THREE.Sphere(new THREE.Vector3(...sphere.centre), sphere.radius);
  const material = new THREE.ShaderMaterial({
    vertexShader: BLOOM_V, fragmentShader: BLOOM_F,
    uniforms: uniforms(THREE, ramp, { uPixel: { value: pixel } }),
    transparent: true, depthWrite: false, blending: night ? THREE.AdditiveBlending : THREE.NormalBlending,
  });
  const object = new THREE.Points(g, material);
  object.frustumCulled = !!sphere;
  object.renderOrder = 6;
  return { ...handle(object, material), setPixel(p) { material.uniforms.uPixel.value = p; } };
}

// Tracks and flows: pairs of vertices (x, y, z), each with its time and value.
// position, time, value: two vertices for each segment, as line segments are given. width: of the thread on the screen, in CSS pixels.
// through: how bright a thread stays where something stands in front of it (0: hidden there), for threads among a city's buildings.
export function createThreads(THREE, { position, time, value, sphere }, { ramp, night = true, through = 0, width = 1.3 } = {}) {
  const segments = position.length / 6, verts = 4 * segments;
  const P = new Float32Array(3 * verts), O = new Float32Array(3 * verts), T = new Float32Array(verts), V = new Float32Array(verts), A = new Float32Array(verts), F = new Float32Array(verts);
  const index = new Uint32Array(6 * segments);
  for (let k = 0; k < segments; k++) {
    const a = 2 * k, b = a + 1, v0 = 4 * k;
    [[a, b, 1, 1], [a, b, -1, 1], [b, a, 1, -1], [b, a, -1, -1]].forEach(([end, other, across, flip], i) => {
      const v = v0 + i;
      for (let c = 0; c < 3; c++) { P[3 * v + c] = position[3 * end + c]; O[3 * v + c] = position[3 * other + c]; }
      T[v] = time[end], V[v] = value[end], A[v] = across, F[v] = flip;
    });
    index.set([v0, v0 + 1, v0 + 2, v0 + 1, v0 + 3, v0 + 2], 6 * k);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(P, 3));
  g.setAttribute('aOther', new THREE.BufferAttribute(O, 3));
  g.setAttribute('aAcross', new THREE.BufferAttribute(A, 1));
  g.setAttribute('aFlip', new THREE.BufferAttribute(F, 1));
  g.setAttribute('aTime', new THREE.BufferAttribute(T, 1));
  g.setAttribute('aValue', new THREE.BufferAttribute(V, 1));
  g.setIndex(new THREE.BufferAttribute(index, 1));
  if (sphere) g.boundingSphere = new THREE.Sphere(new THREE.Vector3(...sphere.centre), sphere.radius);
  const u = uniforms(THREE, ramp, { uViewport: { value: new THREE.Vector2(1, 1) }, uWidth: { value: width } });
  const material = new THREE.ShaderMaterial({
    vertexShader: THREAD_V, fragmentShader: THREAD_F, uniforms: u,
    transparent: true, depthWrite: false, side: THREE.DoubleSide, blending: night ? THREE.AdditiveBlending : THREE.NormalBlending,
  });
  const object = new THREE.Mesh(g, material);
  object.frustumCulled = !!sphere;
  object.renderOrder = 5;
  const buffer = new THREE.Vector2();
  // A pixel more than asked, the soft edge's.
  object.onBeforeRender = (renderer) => { renderer.getDrawingBufferSize(buffer); u.uViewport.value.copy(buffer); u.uWidth.value = (width + 1) * renderer.getPixelRatio(); };
  if (!through) return handle(object, material);
  // The same strip again, drawn only where something stands in front of it, and sharing the first's moment.
  const behind = new THREE.ShaderMaterial({
    vertexShader: THREAD_V, fragmentShader: THREAD_F, uniforms: u, defines: { THROUGH: through.toFixed(3) },
    transparent: true, depthWrite: false, depthFunc: THREE.GreaterDepth, side: THREE.DoubleSide, blending: material.blending,
  });
  const ghost = new THREE.Mesh(g, behind);
  ghost.frustumCulled = object.frustumCulled;
  ghost.renderOrder = 5;
  object.add(ghost);
  const made = handle(object, material);
  return { ...made, dispose() { behind.dispose(); made.dispose(); } };
}
