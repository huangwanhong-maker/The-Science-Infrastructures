// Round markers of a size in pixels at places: events such as earthquakes, and
// the pinned place. Each marker is a quad facing the screen, drawn as a soft
// disc inside a ring.

import * as THREE from 'three';

const vertexShader = /* glsl */ `
in vec2 corner;
in vec3 aPos;
in vec4 aColor;
in float aSize;            // radius, pixels
uniform vec2 uResolution;
out vec2 vOffset;
out vec4 vColor;
out float vSize;
void main() {
  vec4 c = projectionMatrix * modelViewMatrix * vec4(aPos, 1.0);
  float s = aSize + 1.5;
  c.xy += corner * s / (uResolution * 0.5) * c.w;
  gl_Position = c;
  vOffset = corner * s;
  vColor = aColor;
  vSize = aSize;
}`;

const fragmentShader = /* glsl */ `
in vec2 vOffset;
in vec4 vColor;
in float vSize;
uniform float uFill;
out vec4 fragColor;
void main() {
  float r = length(vOffset);
  float inside = clamp(vSize + 0.5 - r, 0.0, 1.0);
  float ring = clamp(1.4 - abs(r - (vSize - 0.7)), 0.0, 1.0);
  float a = max(ring * 0.95, inside * uFill) * vColor.a;
  if (a <= 0.004) discard;
  fragColor = vec4(vColor.rgb, a);
}`;

export class MarkerLayer {
  constructor({ capacity = 256, fill = 0.38, order = 4 } = {}) {
    const g = new THREE.InstancedBufferGeometry();
    g.setAttribute('corner', new THREE.BufferAttribute(new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]), 2));
    g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(12), 3));
    g.setIndex([0, 1, 2, 2, 1, 3]);
    this.geometry = g;
    this.material = new THREE.ShaderMaterial({
      glslVersion: THREE.GLSL3, vertexShader, fragmentShader,
      transparent: true, depthWrite: false, depthTest: true,
      uniforms: { uResolution: { value: new THREE.Vector2(1, 1) }, uFill: { value: fill } },
    });
    this.mesh = new THREE.Mesh(g, this.material);
    this.mesh.frustumCulled = false;
    this.mesh.renderOrder = order;
    this.allocate(capacity);
  }

  allocate(n) {
    this.capacity = n;
    this.pos = new THREE.InstancedBufferAttribute(new Float32Array(3 * n), 3);
    this.color = new THREE.InstancedBufferAttribute(new Float32Array(4 * n), 4);
    this.size = new THREE.InstancedBufferAttribute(new Float32Array(n), 1);
    this.geometry.setAttribute('aPos', this.pos);
    this.geometry.setAttribute('aColor', this.color);
    this.geometry.setAttribute('aSize', this.size);
    this.geometry.instanceCount = 0;
  }

  set(positions, colors, sizes) {
    const n = sizes.length;
    if (n > this.capacity) this.allocate(Math.ceil(n * 1.25));
    this.pos.array.set(positions); this.color.array.set(colors); this.size.array.set(sizes);
    for (const a of [this.pos, this.color, this.size]) a.needsUpdate = true;
    this.geometry.instanceCount = n;
  }

  clear() { this.geometry.instanceCount = 0; }

  frame(resolution) { this.material.uniforms.uResolution.value.copy(resolution); }

  set visible(v) { this.mesh.visible = v; }
  get visible() { return this.mesh.visible; }
}
