// Line segments of a chosen width in pixels: coastlines, boundaries, rivers,
// the graticule, isolines, and the fading trails of particles.
//
// Each segment is an instance of a quad, widened on the screen in the vertex
// shader, so lines keep their width at every zoom and on every projection.
// Per segment: its ends (world coordinates), a birth frame and an opacity
// factor, and its position along its line (for dashes). Trails fade with the
// number of frames since their birth.

import * as THREE from 'three';

const vertexShader = /* glsl */ `
in vec2 corner;
in vec3 aA;
in vec3 aB;
in vec4 aInfo;            // birth frame, opacity factor, position along the line at A and at B
uniform vec2 uResolution; // drawing buffer, pixels
uniform float uWidth;     // pixels
uniform float uFrame;
uniform float uTrail;     // frames over which a trail fades; 0 for lines that do not fade
uniform float uOpacity;
out float vAlpha;
out float vAlong;
out float vSide;
void main() {
  vec4 ca = projectionMatrix * modelViewMatrix * vec4(aA, 1.0);
  vec4 cb = projectionMatrix * modelViewMatrix * vec4(aB, 1.0);
  vec2 hres = uResolution * 0.5;
  vec2 sa = ca.xy / ca.w * hres, sb = cb.xy / cb.w * hres;
  vec2 d = sb - sa;
  float len = length(d);
  vec2 dir = len > 1e-6 ? d / len : vec2(1.0, 0.0);
  vec2 n = vec2(-dir.y, dir.x);
  vec4 c = corner.x < 0.5 ? ca : cb;
  float w = uWidth * 0.5 + 0.5;
  c.xy += (n * corner.y * w + dir * (corner.x < 0.5 ? -0.5 : 0.5) * min(w, 1.0)) / hres * c.w;
  gl_Position = c;
  float fade = uTrail > 0.0 ? clamp(1.0 - (uFrame - aInfo.x) / uTrail, 0.0, 1.0) : 1.0;
  vAlpha = uOpacity * aInfo.y * fade;
  vAlong = mix(aInfo.z, aInfo.w, corner.x);
  vSide = corner.y * w;
}`;

const fragmentShader = /* glsl */ `
in float vAlpha;
in float vAlong;
in float vSide;
uniform vec3 uColor;
uniform float uWidth;
uniform float uDash;      // period along the line (radians of arc); 0 for solid
out vec4 fragColor;
void main() {
  if (vAlpha <= 0.002) discard;
  if (uDash > 0.0 && fract(vAlong / uDash) > 0.6) discard;
  // Coverage across the line: full within the width, falling off over one pixel.
  float edge = clamp(uWidth * 0.5 + 0.5 - abs(vSide), 0.0, 1.0);
  fragColor = vec4(uColor, vAlpha * edge * min(1.0, uWidth));
}`;

export class SegmentLayer {
  constructor({ capacity = 1024, color = '#ffffff', width = 1, opacity = 1, trail = 0, dash = 0, order = 2 } = {}) {
    const geometry = new THREE.InstancedBufferGeometry();
    geometry.setAttribute('corner', new THREE.BufferAttribute(new Float32Array([0, -1, 0, 1, 1, -1, 1, 1]), 2));
    geometry.setIndex([0, 2, 1, 1, 2, 3]);
    // The quad's own position attribute is unused but required by three.js.
    geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(12), 3));
    this.geometry = geometry;
    this.material = new THREE.ShaderMaterial({
      glslVersion: THREE.GLSL3,
      vertexShader, fragmentShader,
      transparent: true, depthWrite: false, depthTest: true,
      uniforms: {
        uResolution: { value: new THREE.Vector2(1, 1) },
        uWidth: { value: width },
        uFrame: { value: 0 },
        uTrail: { value: trail },
        uOpacity: { value: opacity },
        uColor: { value: new THREE.Color(color) },
        uDash: { value: dash },
      },
    });
    this.mesh = new THREE.Mesh(geometry, this.material);
    this.mesh.frustumCulled = false;
    this.mesh.renderOrder = order;
    this.allocate(capacity);
  }

  allocate(capacity) {
    this.capacity = capacity;
    this.a = new THREE.InstancedBufferAttribute(new Float32Array(3 * capacity), 3);
    this.b = new THREE.InstancedBufferAttribute(new Float32Array(3 * capacity), 3);
    this.info = new THREE.InstancedBufferAttribute(new Float32Array(4 * capacity), 4);
    for (const attr of [this.a, this.b, this.info]) attr.setUsage(THREE.DynamicDrawUsage);
    this.geometry.setAttribute('aA', this.a);
    this.geometry.setAttribute('aB', this.b);
    this.geometry.setAttribute('aInfo', this.info);
    this.geometry.instanceCount = 0;
  }

  // Replaces all segments: a and b are Float32Arrays of xyz triples; info of quadruples.
  set(a, b, info) {
    const n = a.length / 3;
    if (n > this.capacity) this.allocate(Math.ceil(n * 1.25));
    this.a.array.set(a); this.b.array.set(b); this.info.array.set(info);
    for (const attr of [this.a, this.b, this.info]) { attr.clearUpdateRanges(); attr.needsUpdate = true; }
    this.geometry.instanceCount = n;
  }

  // Writes n segments starting at instance `start` (for trails held in a ring of slots).
  write(start, n, a, b, info) {
    this.a.array.set(a.subarray(0, 3 * n), 3 * start);
    this.b.array.set(b.subarray(0, 3 * n), 3 * start);
    this.info.array.set(info.subarray(0, 4 * n), 4 * start);
    for (const [attr, k] of [[this.a, 3], [this.b, 3], [this.info, 4]]) {
      attr.clearUpdateRanges();
      attr.addUpdateRange(k * start, k * n);
      attr.needsUpdate = true;
    }
  }

  clear() { this.geometry.instanceCount = 0; }

  set visible(v) { this.mesh.visible = v; }
  get visible() { return this.mesh.visible; }

  style({ color, width, opacity, dash, trail } = {}) {
    const u = this.material.uniforms;
    if (color !== undefined) u.uColor.value.set(color);
    if (width !== undefined) u.uWidth.value = width;
    if (opacity !== undefined) u.uOpacity.value = opacity;
    if (dash !== undefined) u.uDash.value = dash;
    if (trail !== undefined) u.uTrail.value = trail;
  }

  frame(resolution, frameNumber) {
    this.material.uniforms.uResolution.value.copy(resolution);
    this.material.uniforms.uFrame.value = frameNumber;
  }
}

// Segments in world coordinates from polylines in longitude and latitude.
// project(lon, lat, out) writes xyz and returns false where the place is not
// drawn; breaks(p, q) is true where a segment would cross the map's edge.
// Returns { a, b, info } for SegmentLayer.set.
export function polylineSegments(coords, runs, project, breaks, { alpha = 1, subdivide = 0 } = {}) {
  const a = [], b = [], info = [];
  const p = [0, 0, 0], q = [0, 0, 0];
  const R = Math.PI / 180;
  for (const [first, count] of runs) {
    let along = 0;
    for (let k = 0; k + 1 < count; k++) {
      const i = first + k, j = i + 1;
      const lon0 = coords[2 * i], lat0 = coords[2 * i + 1], lon1 = coords[2 * j], lat1 = coords[2 * j + 1];
      const span = Math.max(Math.abs(lon1 - lon0), Math.abs(lat1 - lat0));
      const parts = subdivide > 0 ? Math.max(1, Math.ceil(span / subdivide)) : 1;
      for (let s = 0; s < parts; s++) {
        const t0 = s / parts, t1 = (s + 1) / parts;
        const la = lon0 + (lon1 - lon0) * t0, pa = lat0 + (lat1 - lat0) * t0;
        const lb = lon0 + (lon1 - lon0) * t1, pb = lat0 + (lat1 - lat0) * t1;
        const step = Math.hypot((lb - la) * Math.cos(((pa + pb) / 2) * R), pb - pa) * R;
        if (project(la, pa, p) && project(lb, pb, q) && !breaks(p, q)) {
          a.push(p[0], p[1], p[2]); b.push(q[0], q[1], q[2]);
          info.push(0, alpha, along, along + step);
        }
        along += step;
      }
    }
  }
  return { a: new Float32Array(a), b: new Float32Array(b), info: new Float32Array(info) };
}

// Segments from a flat list of (lon0, lat0, lon1, lat1), as from isolines().
export function pairSegments(pairs, project, breaks, alpha = 1) {
  const n = pairs.length / 4;
  const a = new Float32Array(3 * n), b = new Float32Array(3 * n), info = new Float32Array(4 * n);
  const p = [0, 0, 0], q = [0, 0, 0];
  let m = 0;
  for (let i = 0; i < n; i++) {
    if (!project(pairs[4 * i], pairs[4 * i + 1], p) || !project(pairs[4 * i + 2], pairs[4 * i + 3], q) || breaks(p, q)) continue;
    a.set(p, 3 * m); b.set(q, 3 * m); info.set([0, alpha, 0, 0], 4 * m);
    m++;
  }
  return { a: a.subarray(0, 3 * m), b: b.subarray(0, 3 * m), info: info.subarray(0, 4 * m) };
}
