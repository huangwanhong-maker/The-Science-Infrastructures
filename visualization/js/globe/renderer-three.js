// A globe drawn with three.js: the tiles a TileTree (lod.js) selects, as meshes on the graphics
// card, lit by the Sun and textured with imagery, each placed relative to the camera.
//
// three.js is passed in (createGlobeRenderer(THREE, ...)) rather than imported, so that this
// module, like the rest of the layer, loads in Node for its checks, and so that another backend
// (WebGPU without three.js, a native renderer) can take the same tree and meshes.
//
// Precision. Positions in the body-fixed frame reach 4 × 10⁸ m (the Moon from the Earth), beyond
// single precision's millimetres; the camera is therefore the origin of the scene, and each tile's
// matrix carries its centre relative to the camera, computed in double precision on the CPU
// (the caller gives the body's position relative to the camera; the tile's centre is rotated and
// added here, in JavaScript's doubles). The depth buffer is logarithmic, so that the Moon and a
// pebble a metre away share it (three.js's logarithmicDepthBuffer, Kemen 2012).
//
// Lighting. In the body-fixed frame, by the Sun's direction there. The Earth: the imagery's colour
// lit by the terrain's normal and by the sky (a hemispherical term, more on ground facing up), the night side with the lights of the Black Marble
// when that layer is on; the sea (provenance 6) with the Sun's glint, a Blinn–Phong highlight on the
// figure's normal (the sea's mean surface, without waves). The Moon: the lunar-Lambert law, (1 − A) μ₀ + A · 2μ₀ / (μ₀ + μ), with
// A = 0.6 (McEwen 1991, Icarus 92, 298), so that the full Moon is evenly bright to its limb.
// Without imagery a tile is coloured by height. Two overlays may be laid on the base, each with its
// opacity: a field of data (the sea's temperature, rain) and a reference (coastlines, labels), the
// first lit as the ground is, the second not, so that it reads at night. The provenance of the
// terrain can be laid over: amber where the source filled a gap, blue where the tile is finer than
// its data, violet where detail is synthetic.

const VERTEX = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_vertex>
attribute vec3 up;
attribute float provenance;
attribute float height;
uniform vec3 relCentre;
varying vec2 vUv;
varying vec3 vNormal;
varying vec3 vUp;
varying vec3 vRel;
varying float vProvenance;
varying float vHeight;
void main() {
  vUv = uv;
  vNormal = normal;
  vUp = up;
  vRel = relCentre + position;
  vProvenance = provenance;
  vHeight = height;
  vec4 mv = modelViewMatrix * vec4(position, 1.0);
  gl_Position = projectionMatrix * mv;
  #include <logdepthbuf_vertex>
}
`;

const FRAGMENT = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_fragment>
uniform sampler2D baseMap;
uniform vec3 baseRegion;
uniform float hasBase;
uniform sampler2D nightMap;
uniform vec3 nightRegion;
uniform float hasNight;
uniform sampler2D dataMap;
uniform vec3 dataRegion;
uniform float dataOpacity;
uniform sampler2D refMap;
uniform vec3 refRegion;
uniform float refOpacity;
uniform vec3 sunDir;
uniform float lighting;   // 0 the Earth, 1 a rough body (the Moon, the planets), 2 the Sun (emissive)
uniform float overlay;    // 0 to 1: the provenance laid over
uniform float exposure;
varying vec2 vUv;
varying vec3 vNormal;
varying vec3 vUp;
varying vec3 vRel;
varying float vProvenance;
varying float vHeight;

vec3 byHeight(float h) {
  if (lighting > 0.5) return vec3(0.55, 0.53, 0.50) + 0.00002 * h;
  if (h < 0.0) return mix(vec3(0.05, 0.13, 0.26), vec3(0.20, 0.36, 0.52), clamp(1.0 + h / 6000.0, 0.0, 1.0));
  return mix(vec3(0.33, 0.40, 0.27), vec3(0.62, 0.55, 0.46), clamp(h / 3500.0, 0.0, 1.0));
}
vec3 provenanceColour(float p) {
  if (p > 5.5) return vec3(-1.0);               // the sea's surface: drawn as it is
  if (p > 4.5) return vec3(0.62, 0.36, 0.78);   // synthetic
  if (p > 3.5) return vec3(0.35, 0.58, 0.95);   // finer than the data
  if (p > 1.5 && p < 2.5) return vec3(0.98, 0.66, 0.22); // filled by the source
  return vec3(-1.0);
}
vec2 inRegion(vec3 region) { return vec2(region.x + vUv.x * region.z, 1.0 - region.y - region.z + vUv.y * region.z); }

void main() {
  #include <logdepthbuf_fragment>
  vec3 n = normalize(vNormal);
  vec3 up = normalize(vUp);
  vec3 s = normalize(sunDir);
  vec3 albedo = hasBase > 0.5 ? texture2D(baseMap, inRegion(baseRegion)).rgb : byHeight(vHeight);
  vec3 colour;
  if (lighting > 1.5) {
    // The Sun shines by itself, darkened toward its limb: I(μ) = 1 − u (1 − μ), u = 0.6, the order of
    // the visible's linear coefficient (a perceptual rendering, not a photometric one).
    vec3 v = normalize(-vRel);
    float mu = max(dot(n, v), 0.0);
    colour = albedo * (1.0 - 0.6 * (1.0 - mu)) * 1.6;
  } else if (lighting > 0.5) {
    vec3 v = normalize(-vRel);
    float mu0 = max(dot(n, s), 0.0), mu = max(dot(n, v), 0.0);
    float A = 0.6;
    float f = (1.0 - A) * mu0 + A * 2.0 * mu0 / (mu0 + mu + 1e-4);
    colour = albedo * (0.015 + 1.15 * f);
  } else {
    float lambert = max(dot(n, s), 0.0);
    float daylight = smoothstep(-0.10, 0.12, dot(up, s));
    // The Sun on the slope, and the sky's light, more on ground that faces up: the imagery carries
    // its own shadows already, so the added shading is kept gentle.
    float sky = 0.55 + 0.45 * dot(n, up);
    colour = albedo * (0.30 * sky + 0.80 * lambert);
    if (vProvenance > 5.5) {
      // The sea: the Sun's glint on its mean surface, and a little of the sky's light.
      vec3 v = normalize(-vRel);
      vec3 h = normalize(s + v);
      float glint = pow(max(dot(up, h), 0.0), 180.0) * 0.9 + pow(max(dot(up, h), 0.0), 18.0) * 0.06;
      colour += vec3(1.0, 0.94, 0.86) * glint * smoothstep(0.0, 0.08, dot(up, s));
    }
    if (hasNight > 0.5) {
      vec3 lights = texture2D(nightMap, inRegion(nightRegion)).rgb;
      colour += lights * lights * 1.6 * (1.0 - daylight);
    }
    colour = mix(albedo * 0.025, colour, max(daylight, 0.0) * 0.9 + 0.1);
  }
  if (dataOpacity > 0.0) {
    vec4 d = texture2D(dataMap, inRegion(dataRegion));
    float lit = lighting > 0.5 ? 1.0 : 0.35 + 0.65 * smoothstep(-0.10, 0.12, dot(up, s));
    colour = mix(colour, d.rgb * lit, d.a * dataOpacity);
  }
  if (refOpacity > 0.0) {
    vec4 r = texture2D(refMap, inRegion(refRegion));
    colour = mix(colour, r.rgb, r.a * refOpacity);
  }
  vec3 pc = provenanceColour(vProvenance);
  if (overlay > 0.0 && pc.r >= 0.0) colour = mix(colour, pc * (0.35 + 0.65 * max(dot(n, s), 0.25)), 0.55 * overlay);
  gl_FragColor = vec4(colour * exposure, 1.0);
  #include <colorspace_fragment>
}
`;

const ATMOSPHERE_VERTEX = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_vertex>
varying vec3 vNormal;
varying vec3 vView;
void main() {
  vec4 mv = modelViewMatrix * vec4(position, 1.0);
  vNormal = normalize(normalMatrix * normal);
  vView = normalize(-mv.xyz);
  gl_Position = projectionMatrix * mv;
  #include <logdepthbuf_vertex>
}
`;
const ATMOSPHERE_FRAGMENT = /* glsl */ `
#include <common>
#include <logdepthbuf_pars_fragment>
uniform vec3 sunView;
varying vec3 vNormal;
varying vec3 vView;
void main() {
  #include <logdepthbuf_fragment>
  float rim = 1.0 - abs(dot(normalize(vNormal), normalize(vView)));
  float glow = pow(rim, 3.2);
  float lit = smoothstep(-0.35, 0.45, dot(normalize(vNormal), normalize(sunView)));
  vec3 colour = mix(vec3(0.30, 0.50, 0.95), vec3(0.85, 0.62, 0.55), pow(1.0 - lit, 2.0)) * glow * (0.12 + 1.1 * lit);
  gl_FragColor = vec4(colour, glow * (0.15 + 0.85 * lit));
  #include <colorspace_fragment>
}
`;

// lighting: 'earth' (land and sea with the night's lights), 'moon' or 'body' (a rough surface, lunar-Lambert),
// 'sun' (emissive, darkened toward the limb).
export function createGlobeRenderer(THREE, { fig, lighting = 'earth', atmosphere = false, maxAnisotropy = 8 } = {}) {
  const group = new THREE.Group();
  group.matrixAutoUpdate = false;
  const meshes = new Set();
  const blank = new THREE.DataTexture(new Uint8Array([128, 128, 128, 255]), 1, 1);
  blank.needsUpdate = true;

  function material() {
    return new THREE.ShaderMaterial({
      vertexShader: VERTEX,
      fragmentShader: FRAGMENT,
      side: THREE.DoubleSide,
      uniforms: {
        baseMap: { value: blank }, baseRegion: { value: new THREE.Vector3(0, 0, 1) }, hasBase: { value: 0 },
        nightMap: { value: blank }, nightRegion: { value: new THREE.Vector3(0, 0, 1) }, hasNight: { value: 0 },
        dataMap: { value: blank }, dataRegion: { value: new THREE.Vector3(0, 0, 1) }, dataOpacity: { value: 0 },
        refMap: { value: blank }, refRegion: { value: new THREE.Vector3(0, 0, 1) }, refOpacity: { value: 0 },
        sunDir: { value: new THREE.Vector3(1, 0, 0) }, relCentre: { value: new THREE.Vector3() },
        lighting: { value: lighting === 'sun' ? 2 : lighting === 'moon' || lighting === 'body' ? 1 : 0 }, overlay: { value: 0 }, exposure: { value: 1 },
      },
    });
  }

  let shell = null;
  if (atmosphere) {
    const geometry = new THREE.SphereGeometry(1, 96, 48);
    shell = new THREE.Mesh(geometry, new THREE.ShaderMaterial({
      vertexShader: ATMOSPHERE_VERTEX, fragmentShader: ATMOSPHERE_FRAGMENT,
      uniforms: { sunView: { value: new THREE.Vector3(1, 0, 0) } },
      side: THREE.BackSide, transparent: true, depthWrite: false, blending: THREE.AdditiveBlending,
    }));
    shell.matrixAutoUpdate = false;
    shell.frustumCulled = false;
    group.add(shell);
  }

  const api = {
    group,
    // TileTree's events.
    ready(tile) {
      const m = tile.mesh;
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.BufferAttribute(m.position, 3));
      g.setAttribute('normal', new THREE.BufferAttribute(m.normal, 3));
      g.setAttribute('up', new THREE.BufferAttribute(m.up, 3));
      g.setAttribute('uv', new THREE.BufferAttribute(m.uv, 2));
      g.setAttribute('provenance', new THREE.BufferAttribute(m.provenance, 1));
      g.setAttribute('height', new THREE.BufferAttribute(m.height, 1));
      g.setIndex(new THREE.BufferAttribute(m.index, 1));
      // One material, in an array, so that the geometry's groups can choose the quarters drawn.
      const mesh = new THREE.Mesh(g, [material()]);
      mesh.userData.quarters = m.quarters;
      mesh.userData.total = m.index.length;
      mesh.matrixAutoUpdate = false;
      mesh.frustumCulled = false;
      mesh.visible = false;
      group.add(mesh);
      meshes.add(mesh);
      tile.gpu = mesh;
    },
    released(tile) {
      const mesh = tile.gpu;
      if (!mesh) return;
      group.remove(mesh);
      meshes.delete(mesh);
      mesh.geometry.dispose();
      mesh.material[0].dispose();
      tile.gpu = null;
    },
    imageReady(tile, layer, image) {
      const texture = new THREE.Texture(image);
      texture.colorSpace = THREE.SRGBColorSpace;
      texture.flipY = false;
      texture.wrapS = texture.wrapT = THREE.ClampToEdgeWrapping;
      texture.anisotropy = maxAnisotropy;
      texture.generateMipmaps = true;
      texture.minFilter = THREE.LinearMipmapLinearFilter;
      texture.needsUpdate = true;
      tile.images.get(layer).gpu = texture;
    },
    imageReleased(tile, layer) {
      const img = tile.images.get(layer);
      img?.gpu?.dispose();
      if (img?.image?.close) img.image.close();
    },

    // Draw the selected tiles. frame: { bodyRelative: [x, y, z] the body's centre relative to the
    // camera, world axes (m, double); rotation: 3 x 3 body-fixed → world (row-major arrays);
    // sunBody: the Sun's direction in the body-fixed frame; base, night: layer ids; data, reference:
    // { layer, opacity } or null; overlay: the provenance's opacity, 0..1; sunView: the Sun's
    // direction in world axes (for the atmosphere) }.
    draw(selection, frame) {
      for (const mesh of meshes) mesh.visible = false;
      const R = frame.rotation, c = frame.bodyRelative;
      const sun = new THREE.Vector3(...frame.sunBody);
      for (const { tile, images, quarters } of selection) {
        const mesh = tile.gpu;
        if (!mesh) continue;
        // The whole tile, or the quarters of it its children do not yet cover.
        const g = mesh.geometry;
        g.clearGroups();
        if (!quarters || !mesh.userData.quarters) g.addGroup(0, mesh.userData.total, 0);
        else for (const q of quarters) g.addGroup(...mesh.userData.quarters[q], 0);
        const ctr = tile.mesh.centre;
        // The tile's centre relative to the camera, world axes, in double precision.
        const dx = c[0] + R[0][0] * ctr[0] + R[0][1] * ctr[1] + R[0][2] * ctr[2];
        const dy = c[1] + R[1][0] * ctr[0] + R[1][1] * ctr[1] + R[1][2] * ctr[2];
        const dz = c[2] + R[2][0] * ctr[0] + R[2][1] * ctr[1] + R[2][2] * ctr[2];
        mesh.matrix.set(R[0][0], R[0][1], R[0][2], dx, R[1][0], R[1][1], R[1][2], dy, R[2][0], R[2][1], R[2][2], dz, 0, 0, 0, 1);
        mesh.matrixWorldNeedsUpdate = true;
        const u = mesh.material[0].uniforms;
        u.sunDir.value.copy(sun);
        // The tile's centre relative to the camera in body axes (for the direction of view).
        u.relCentre.value.set(R[0][0] * dx + R[1][0] * dy + R[2][0] * dz, R[0][1] * dx + R[1][1] * dy + R[2][1] * dz, R[0][2] * dx + R[1][2] * dy + R[2][2] * dz);
        u.overlay.value = frame.overlay || 0;
        u.exposure.value = frame.exposure ?? 1;
        const texture = (id) => {
          const i = id && images.find((im) => im && im.layer === id);
          const tex = i && i.owner.images.get(i.layer)?.gpu;
          return tex ? { tex, region: i.region } : null;
        };
        const base = texture(frame.base), night = texture(frame.night);
        const data = frame.data && texture(frame.data.layer), ref = frame.reference && texture(frame.reference.layer);
        u.hasBase.value = base ? 1 : 0;
        if (base) { u.baseMap.value = base.tex; u.baseRegion.value.set(...base.region); }
        u.hasNight.value = night ? 1 : 0;
        if (night) { u.nightMap.value = night.tex; u.nightRegion.value.set(...night.region); }
        u.dataOpacity.value = data ? frame.data.opacity ?? 1 : 0;
        if (data) { u.dataMap.value = data.tex; u.dataRegion.value.set(...data.region); }
        u.refOpacity.value = ref ? frame.reference.opacity ?? 1 : 0;
        if (ref) { u.refMap.value = ref.tex; u.refRegion.value.set(...ref.region); }
        mesh.visible = true;
      }
      if (shell) {
        const r = fig.a + 60000, rb = fig.b + 60000;
        shell.matrix.set(R[0][0] * r, R[0][1] * r, R[0][2] * rb, c[0], R[1][0] * r, R[1][1] * r, R[1][2] * rb, c[1], R[2][0] * r, R[2][1] * r, R[2][2] * rb, c[2], 0, 0, 0, 1);
        shell.matrixWorldNeedsUpdate = true;
        if (frame.sunView) shell.material.uniforms.sunView.value.set(...frame.sunView);
        shell.visible = frame.atmosphere !== false;
      }
    },
    dispose() {
      for (const mesh of meshes) { mesh.geometry.dispose(); mesh.material[0].dispose(); }
      meshes.clear();
      blank.dispose();
      if (shell) { shell.geometry.dispose(); shell.material.dispose(); }
    },
  };
  return api;
}
