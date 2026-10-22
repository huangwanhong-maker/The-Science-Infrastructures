// Star points for the 3D views: per-star colour and size, soft round discs,
// additive blending, compatible with the logarithmic depth buffer.

import * as THREE from 'three';

const vertexShader = /* glsl */ `
  #include <common>
  #include <logdepthbuf_pars_vertex>
  attribute vec3 starColour;
  attribute float size;
  uniform float pixelRatio;
  uniform float sizeScale;
  varying vec3 vColour;
  void main() {
    vColour = starColour;
    vec4 mv = modelViewMatrix * vec4(position, 1.0);
    gl_Position = projectionMatrix * mv;
    gl_PointSize = max(size * sizeScale, 1.0) * pixelRatio;
    #include <logdepthbuf_vertex>
  }
`;

const fragmentShader = /* glsl */ `
  #include <common>
  #include <logdepthbuf_pars_fragment>
  varying vec3 vColour;
  void main() {
    #include <logdepthbuf_fragment>
    float r = length(gl_PointCoord - 0.5) * 2.0;
    float a = smoothstep(1.0, 0.0, r);
    a *= a;
    if (a < 0.02) discard;
    gl_FragColor = vec4(vColour * a, a);
  }
`;

export function starPoints(positions, colours, sizes, pixelRatio) {
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute('starColour', new THREE.BufferAttribute(colours, 3));
  geometry.setAttribute('size', new THREE.BufferAttribute(sizes, 1));
  const material = new THREE.ShaderMaterial({
    uniforms: { pixelRatio: { value: pixelRatio }, sizeScale: { value: 1 } },
    vertexShader,
    fragmentShader,
    transparent: true,
    depthWrite: false,
    blending: THREE.AdditiveBlending,
  });
  const points = new THREE.Points(geometry, material);
  points.frustumCulled = false;
  return points;
}

// Point size (px) and intensity for an apparent or absolute magnitude, relative
// to a reference magnitude that maps to the smallest visible point.
export function magnitudeStyle(m, reference) {
  const flux = 10 ** (-0.4 * (m - reference));
  const size = Math.min(8, Math.max(1.2, 1.2 * flux ** 0.25));
  const intensity = Math.min(1, Math.max(0.28, 0.28 + 0.18 * Math.log10(Math.max(flux, 1))));
  return { size, intensity };
}
