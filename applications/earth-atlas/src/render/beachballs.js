// Beach balls: earthquake mechanisms drawn where they happened. Each is a disc
// facing the screen, the lower hemisphere of the source seen from above in
// equal-area projection (src/core/moment.js), turned so that its north is the
// map's north there. A fragment is in a compressional quadrant, and coloured,
// where x^T M x > 0 for the direction x it stands for; the dilatational
// quadrants are pale.

import * as THREE from 'three';

const vertexShader = /* glsl */ `
in vec2 corner;
in vec3 aPos;
in float aSize;            // radius, pixels
in vec4 aM1;               // M_nn, M_ne, M_nd, M_ee (north, east, down)
in vec2 aM2;               // M_ed, M_dd
in vec2 aNorth;            // the direction of north on the screen (unit, y up)
in vec4 aColor;            // the compressional colour and the opacity
uniform vec2 uResolution;
out vec2 vOffset;
out float vSize;
out vec4 vM1;
out vec2 vM2;
out vec2 vNorth;
out vec4 vColor;
void main() {
  vec4 c = projectionMatrix * modelViewMatrix * vec4(aPos, 1.0);
  float s = aSize + 1.5;
  c.xy += corner * s / (uResolution * 0.5) * c.w;
  gl_Position = c;
  vOffset = corner * s;
  vSize = aSize; vM1 = aM1; vM2 = aM2; vNorth = aNorth; vColor = aColor;
}`;

const fragmentShader = /* glsl */ `
in vec2 vOffset;
in float vSize;
in vec4 vM1;
in vec2 vM2;
in vec2 vNorth;
in vec4 vColor;
uniform vec3 uPale;
uniform vec3 uRing;
out vec4 fragColor;
void main() {
  float r = length(vOffset);
  if (r > vSize + 1.0) discard;
  vec2 east = vec2(vNorth.y, -vNorth.x);
  float x = dot(vOffset, east) / vSize, y = dot(vOffset, vNorth) / vSize;
  float rr = min(1.0, length(vec2(x, y)));
  float theta = 2.0 * asin(rr * 0.70710678);
  float az = atan(x, y);
  vec3 d = vec3(sin(theta) * cos(az), sin(theta) * sin(az), cos(theta));
  float a = vM1.x * d.x * d.x + vM1.w * d.y * d.y + vM2.y * d.z * d.z
          + 2.0 * (vM1.y * d.x * d.y + vM1.z * d.x * d.z + vM2.x * d.y * d.z);
  // Soften the edge between quadrants over about a pixel.
  float w = fwidth(a) + 1e-6;
  vec3 col = mix(uPale, vColor.rgb, smoothstep(-w, w, a));
  float ring = clamp(1.0 - abs(r - (vSize - 0.4)) / 0.9, 0.0, 1.0);
  col = mix(col, uRing, ring * 0.85);
  float alpha = clamp(vSize + 0.6 - r, 0.0, 1.0) * vColor.a;
  fragColor = vec4(col, alpha);
}`;

export class BeachBallLayer {
  constructor({ capacity = 512, order = 5, pale = '#fbf1ea', ring = '#2b1d2d' } = {}) {
    const g = new THREE.InstancedBufferGeometry();
    g.setAttribute('corner', new THREE.BufferAttribute(new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]), 2));
    g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(12), 3));
    g.setIndex([0, 1, 2, 2, 1, 3]);
    this.geometry = g;
    this.material = new THREE.ShaderMaterial({
      glslVersion: THREE.GLSL3, vertexShader, fragmentShader,
      transparent: true, depthWrite: false, depthTest: true,
      uniforms: { uResolution: { value: new THREE.Vector2(1, 1) }, uPale: { value: new THREE.Color(pale) }, uRing: { value: new THREE.Color(ring) } },
    });
    this.mesh = new THREE.Mesh(g, this.material);
    this.mesh.frustumCulled = false;
    this.mesh.renderOrder = order;
    this.allocate(capacity);
  }

  allocate(n) {
    this.capacity = n;
    const add = (name, size) => {
      const a = new THREE.InstancedBufferAttribute(new Float32Array(size * n), size);
      this.geometry.setAttribute(name, a);
      return a;
    };
    this.attrs = { aPos: add('aPos', 3), aSize: add('aSize', 1), aM1: add('aM1', 4), aM2: add('aM2', 2), aNorth: add('aNorth', 2), aColor: add('aColor', 4) };
    this.geometry.instanceCount = 0;
  }

  // arrays: { aPos, aSize, aM1, aM2, aNorth, aColor } as Float32Arrays of n instances.
  set(arrays, n) {
    if (n > this.capacity) this.allocate(Math.ceil(n * 1.25));
    for (const [k, a] of Object.entries(this.attrs)) { a.array.set(arrays[k].subarray(0, n * a.itemSize)); a.needsUpdate = true; }
    this.geometry.instanceCount = n;
  }

  clear() { this.geometry.instanceCount = 0; }

  frame(resolution) { this.material.uniforms.uResolution.value.copy(resolution); }
}
