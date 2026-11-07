// Columns standing on the ground: cells (the hexagons of a discrete global grid, or any rings) raised to
// heights and coloured, a value at a place read at a glance. They are drawn by one mesh whose vertices name
// their cell; each cell's height and colour are a texel of a table, so that a new year of a record, or a new
// weighting of an index, is a new table and not a new mesh. The columns' walls and tops are shaded flat, from
// the screen-space derivatives of their positions, which needs no normals in the buffers.
//
// Coordinates are the caller's: base rings and an `up` for each cell in the object's own axes, so that one
// module serves the globe (body axes, relative to the set's centre) and the city (three.js's local axes).
// three.js is passed in. A cell whose height is negative is not drawn (no record there). `grow` carries all
// heights from 0 to 1 as a set arrives; `opacity` fades a set outside the scales it means anything at.

const VERTEX = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_vertex>
attribute float cell;
attribute float top;
attribute vec3 riseAlong;
uniform highp sampler2D table;
uniform float tableWidth;
uniform float grow;
uniform float lift;
uniform float picked;
varying vec3 vColour;
varying vec3 vView;
varying float vTop;
varying float vPicked;
void main() {
  int i = int(cell + 0.5);
  int w = int(tableWidth);
  vec4 t = texelFetch(table, ivec2(i % w, i / w), 0);
  if (t.a < 0.0) {
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
    vColour = vec3(0.0); vView = vec3(0.0); vTop = 0.0; vPicked = 0.0;
    return;
  }
  vec3 p = position + riseAlong * (lift + top * t.a * grow);
  vec4 mv = modelViewMatrix * vec4(p, 1.0);
  vView = mv.xyz;
  vColour = t.rgb;
  vTop = top;
  vPicked = abs(cell - picked) < 0.5 ? 1.0 : 0.0;
  gl_Position = projectionMatrix * mv;
  #include <logdepthbuf_vertex>
}
`;

const FRAGMENT = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_fragment>
uniform vec3 sunView;
uniform float opacity;
uniform float shade;
varying vec3 vColour;
varying vec3 vView;
varying float vTop;
varying float vPicked;
void main() {
  #include <logdepthbuf_fragment>
  vec3 n = normalize(cross(dFdx(vView), dFdy(vView)));
  float sun = max(dot(n, normalize(sunView)), 0.0);
  // The foot of a wall sinks into shade and its head catches the light, which reads as height without any line.
  float wall = mix(0.52, 1.0, smoothstep(0.0, 1.0, vTop));
  vec3 c = vColour * wall * (1.0 - shade + shade * (0.45 + 0.75 * sun));
  c = mix(c, vec3(1.0, 0.94, 0.97), vPicked * 0.55);
  gl_FragColor = vec4(c, opacity);
  #include <colorspace_fragment>
}
`;

export const linearOf = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);

// cells: [{ ring: [[x, y, z], …] (the base, in the object's axes), up: [x, y, z] (unit) }]; shrink pulls each ring
// toward its centre, a little light between columns.
export function createPrisms(THREE, cells, { opacity = 1, shrink = 0.9, shade = 1, order = 3 } = {}) {
  let verts = 0;
  for (const c of cells) verts += 3 * (c.ring.length - 2) + 6 * c.ring.length;
  const position = new Float32Array(3 * verts), rise = new Float32Array(3 * verts), cell = new Float32Array(verts), top = new Float32Array(verts);
  let v = 0;
  const put = (p, u, k, t) => {
    position[3 * v] = p[0], position[3 * v + 1] = p[1], position[3 * v + 2] = p[2];
    rise[3 * v] = u[0], rise[3 * v + 1] = u[1], rise[3 * v + 2] = u[2];
    cell[v] = k, top[v] = t;
    v++;
  };
  cells.forEach((c, k) => {
    const n = c.ring.length, m = [0, 0, 0];
    for (const p of c.ring) for (let i = 0; i < 3; i++) m[i] += p[i] / n;
    const r = c.ring.map((p) => [0, 1, 2].map((i) => m[i] + (p[i] - m[i]) * shrink));
    for (let i = 1; i < n - 1; i++) { put(r[0], c.up, k, 1); put(r[i], c.up, k, 1); put(r[i + 1], c.up, k, 1); }
    for (let i = 0; i < n; i++) {
      const a = r[i], b = r[(i + 1) % n];
      put(a, c.up, k, 0); put(b, c.up, k, 0); put(b, c.up, k, 1);
      put(a, c.up, k, 0); put(b, c.up, k, 1); put(a, c.up, k, 1);
    }
  });
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(position, 3));
  g.setAttribute('riseAlong', new THREE.BufferAttribute(rise, 3));
  g.setAttribute('cell', new THREE.BufferAttribute(cell, 1));
  g.setAttribute('top', new THREE.BufferAttribute(top, 1));
  const width = Math.max(1, Math.ceil(Math.sqrt(cells.length))), hgt = Math.max(1, Math.ceil(cells.length / width));
  const data = new Float32Array(4 * width * hgt).fill(-1);
  const table = new THREE.DataTexture(data, width, hgt, THREE.RGBAFormat, THREE.FloatType);
  table.minFilter = table.magFilter = THREE.NearestFilter;
  table.needsUpdate = true;
  const material = new THREE.ShaderMaterial({
    vertexShader: VERTEX, fragmentShader: FRAGMENT, side: THREE.DoubleSide,
    transparent: opacity < 1,
    uniforms: {
      table: { value: table }, tableWidth: { value: width }, grow: { value: 1 }, lift: { value: 0 }, picked: { value: -1 },
      sunView: { value: new THREE.Vector3(0, 0, 1) }, opacity: { value: opacity }, shade: { value: shade },
    },
  });
  const mesh = new THREE.Mesh(g, material);
  mesh.frustumCulled = false;
  mesh.renderOrder = order;
  return {
    object: mesh, material, count: cells.length,
    // heights: m (negative: not drawn); colours: [r, g, b] in sRGB, or null for a cell not drawn.
    set(heights, colours) {
      for (let k = 0; k < cells.length; k++) {
        const c = colours[k], h = heights[k];
        if (!c || !Number.isFinite(h) || h < 0) { data[4 * k + 3] = -1; continue; }
        data[4 * k] = linearOf(c[0]), data[4 * k + 1] = linearOf(c[1]), data[4 * k + 2] = linearOf(c[2]), data[4 * k + 3] = h;
      }
      table.needsUpdate = true;
    },
    setOpacity(a) { material.uniforms.opacity.value = a; material.transparent = a < 1; material.depthWrite = a > 0.99; mesh.visible = a > 0.004; },
    setGrow(t) { material.uniforms.grow.value = t; },
    setLift(m) { material.uniforms.lift.value = m; },
    setPicked(k) { material.uniforms.picked.value = k ?? -1; },
    setSun(x, y, z) { material.uniforms.sunView.value.set(x, y, z); },
    dispose() { g.dispose(); table.dispose(); material.dispose(); },
  };
}

// Columns rising over about a second: t in [0, 1] to a height's share, with a soft landing.
export const rise = (t) => { const u = Math.min(1, Math.max(0, t)); return 1 - (1 - u) ** 3; };
