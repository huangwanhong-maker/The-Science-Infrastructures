// The Earth's surface as drawn by the GPU: the base map (relief or imagery,
// with finer detail where the user looks), lights at night, forecast clouds,
// one data layer (a field, a difference of two fields, a choropleth of
// countries, or a map of classes and how they changed between two periods),
// and day and night, on the globe or on a flat map.
//
// One fragment shader serves every projection. For the globe it reads the
// place from the sphere's surface; for a flat map it inverts the projection at
// each pixel (transcribing src/core/projections.js). It then samples, at that
// place, the base, the data and the Sun's altitude. Fields are sampled from
// float textures holding the grid exactly as stored, with the bilinear rule of
// src/core/grid.js (weights renormalised beside missing values) and linear
// interpolation between two times, so the GPU draws what the probe reports.

import * as THREE from 'three';

export const GPU_MISSING = -1e30;

const vertexShader = /* glsl */ `
out vec3 vEarth;
out vec2 vPlane;
void main() {
  vEarth = position;
  vPlane = (modelMatrix * vec4(position, 1.0)).xy;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}`;

const fragmentShader = /* glsl */ `
precision highp float;
precision highp int;
precision highp usampler2D;

in vec3 vEarth;
in vec2 vPlane;
out vec4 fragColor;

uniform int uProjection;          // 0 globe, 1 equirectangular, 2 Equal Earth
uniform float uLambda0;           // central meridian, radians
uniform vec2 uExtent;             // half-width and half-height of the map

uniform sampler2D uBase;          // relief: RGB tinted and shaded; A: shading x 127.5
uniform sampler2D uImagery;       // imagery of the surface by day (RGB), black where there is none
uniform float uImageryOn;
uniform sampler2D uLights;        // lights at night
uniform float uLightsOn;
uniform sampler2D uDetail;        // finer imagery or relief of one region
uniform vec4 uDetailBox;          // west longitude, south latitude, width, height (degrees)
uniform float uDetailOn;
uniform float uDetailKind;        // 0 imagery (colour); 1 shaded relief (R: grey, G: 1 on land, 0 on water)
uniform float uShadeFlat;         // the grey of flat land in the shaded relief, 0 to 1

uniform int uMode;                // 0 none, 1 field, 2 speed of a vector, 3 countries, 4 difference A - B,
                                  // 5 classes, 6 classes where they differ from a second map
uniform sampler2D uA0, uA1, uB0, uB1;
uniform float uMix, uMixB;        // weights of the later time, for A and for B
uniform vec4 uGrid, uGridB;       // lon0, dlon, lat0, dlat (degrees)
uniform ivec3 uGridN, uGridNB;    // nlon, nlat, periodic
uniform sampler2D uLut;
uniform vec3 uRange;              // low, high (stored units), scale: 0 linear, 1 logarithmic, 2 square root
uniform float uOpacity;
uniform float uReliefOnData;
uniform usampler2D uCountries;    // 1 + index of the country at each pixel, 0 for none
uniform sampler2D uCountryColours;
uniform usampler2D uK0, uK1;      // class of each cell (0 for none): the map shown, and another
uniform vec4 uGridK, uGridKB;
uniform ivec3 uGridNK, uGridNKB;
uniform sampler2D uClassColours;  // RGBA for each class value

uniform sampler2D uC0, uC1;       // cloud cover (%)
uniform float uMixC, uCloudOn;
uniform vec4 uGridC;
uniform ivec3 uGridNC;

uniform vec3 uSun;                // unit vector toward the Sun, Earth-fixed
uniform float uNight;             // 0 none, 1 full
uniform vec3 uNightTint;

const float PI = 3.141592653589793;
const float MISSING = -1.0e29;
const float A1 = 1.340264, A2 = -0.081106, A3 = 0.000893, A4 = 0.003796, M = 0.8660254037844386;

vec2 equalEarthInverse(vec2 p) {
  float t = p.y;
  for (int i = 0; i < 6; i++) {
    float t2 = t * t, t6 = t2 * t2 * t2;
    float fy = t * (A1 + A2 * t2 + t6 * (A3 + A4 * t2)) - p.y;
    float fpy = A1 + 3.0 * A2 * t2 + t6 * (7.0 * A3 + 9.0 * A4 * t2);
    t -= fy / fpy;
  }
  float t2 = t * t, t6 = t2 * t2 * t2;
  float lambda = M * p.x * (A1 + 3.0 * A2 * t2 + t6 * (7.0 * A3 + 9.0 * A4 * t2)) / cos(t);
  return vec2(lambda, asin(clamp(sin(t) / M, -1.0, 1.0)));
}

float cell(sampler2D t, int c, int r) { return texelFetch(t, ivec2(c, r), 0).r; }

// Bilinear interpolation with renormalised weights (src/core/grid.js).
float fieldAt(sampler2D t, vec2 ll, vec4 grid, ivec3 gridN) {
  float nlon = float(gridN.x), nlat = float(gridN.y);
  float x = (ll.x - grid.x) / grid.y;
  if (gridN.z == 1) x = mod(x, nlon);
  else if (x < -0.5 || x > nlon - 0.5) return MISSING;
  float y = (grid.z - ll.y) / grid.w;
  if (y < -0.5 || y > nlat - 0.5) return MISSING;
  y = clamp(y, 0.0, nlat - 1.0);
  int c0 = int(floor(x)), r0 = int(floor(y));
  float fx = x - float(c0), fy = y - float(r0);
  int c1 = c0 + 1, r1 = min(r0 + 1, gridN.y - 1);
  if (c1 >= gridN.x) c1 = gridN.z == 1 ? c1 - gridN.x : gridN.x - 1;
  c0 = max(c0, 0);
  vec4 v = vec4(cell(t, c0, r0), cell(t, c1, r0), cell(t, c0, r1), cell(t, c1, r1));
  vec4 w = vec4((1.0 - fx) * (1.0 - fy), fx * (1.0 - fy), (1.0 - fx) * fy, fx * fy);
  vec4 ok = vec4(greaterThan(v, vec4(MISSING)));
  float weight = dot(w, ok);
  if (weight < 0.5) return MISSING;
  return dot(w * ok, v) / weight;
}

// The class of the cell containing a place: no interpolation between classes.
uint classAt(usampler2D t, vec2 ll, vec4 grid, ivec3 gridN) {
  float x = floor((ll.x - grid.x) / grid.y + 0.5);
  if (gridN.z == 1) x = mod(x, float(gridN.x));
  else if (x < 0.0 || x > float(gridN.x) - 1.0) return 0u;
  float y = floor((grid.z - ll.y) / grid.w + 0.5);
  if (y < 0.0 || y > float(gridN.y) - 1.0) return 0u;
  return texelFetch(t, ivec2(int(x), int(y)), 0).r;
}

float atTime(sampler2D t0, sampler2D t1, float mixT, vec2 ll, vec4 grid, ivec3 gridN) {
  float a = fieldAt(t0, ll, grid, gridN);
  if (mixT <= 0.0 || a < MISSING) return a;
  float b = fieldAt(t1, ll, grid, gridN);
  return b < MISSING ? MISSING : mix(a, b, mixT);
}

vec4 colourOf(float v) {
  float t;
  if (uRange.z > 0.5 && uRange.z < 1.5) t = (log(max(v, 1e-30)) - log(uRange.x)) / (log(uRange.y) - log(uRange.x));
  else if (uRange.z > 1.5) t = (sqrt(max(v, 0.0)) - sqrt(max(uRange.x, 0.0))) / (sqrt(uRange.y) - sqrt(max(uRange.x, 0.0)));
  else t = (v - uRange.x) / (uRange.y - uRange.x);
  return texture(uLut, vec2(clamp(t, 0.0, 1.0) * (255.0 / 256.0) + 0.5 / 256.0, 0.5));
}

// Darkening by the Sun's altitude, in steps at sunset and at the ends of civil,
// nautical and astronomical twilight (-0.833, -6, -12, -18 degrees).
float darkness(float alt) {
  float d = 0.30 * smoothstep(-0.533, -1.133, alt);
  d += 0.26 * smoothstep(-5.75, -6.25, alt);
  d += 0.24 * smoothstep(-11.75, -12.25, alt);
  d += 0.20 * smoothstep(-17.75, -18.25, alt);
  return d;
}

// A texture sampled with gradients taken across the seam at 180 degrees, so
// that the mipmap level does not jump there.
vec4 seamless(sampler2D t, vec2 uv) {
  vec2 uvs = vec2(fract(uv.x + 0.5), uv.y);
  vec2 dx = dFdx(uv), dy = dFdy(uv), dxs = dFdx(uvs), dys = dFdy(uvs);
  if (abs(dxs.x) < abs(dx.x)) dx.x = dxs.x;
  if (abs(dys.x) < abs(dy.x)) dy.x = dys.x;
  return textureGrad(t, uv, dx, dy);
}

void main() {
  vec3 earth;
  vec2 ll;
  if (uProjection == 0) {
    earth = normalize(vEarth);
    ll = vec2(atan(earth.y, earth.x), asin(clamp(earth.z, -1.0, 1.0))) * (180.0 / PI);
  } else {
    if (abs(vPlane.y) > uExtent.y) discard;
    vec2 lp = uProjection == 1 ? vPlane : equalEarthInverse(vPlane);
    if (abs(lp.x) > PI) discard;
    float lon = mod(lp.x + uLambda0 + PI, 2.0 * PI) - PI;
    earth = vec3(cos(lp.y) * cos(lon), cos(lp.y) * sin(lon), sin(lp.y));
    ll = vec2(lon, lp.y) * (180.0 / PI);
  }

  vec2 uv = vec2((ll.x + 180.0) / 360.0, (90.0 - ll.y) / 180.0);
  vec4 relief = seamless(uBase, uv);
  vec3 colour = relief.rgb;
  float shade = relief.a * 2.0;
  float shadeData = shade;          // the shading applied to data drawn over the relief
  if (uImageryOn > 0.5) {
    vec3 img = seamless(uImagery, uv).rgb;
    // Where the imagery has no data (black), the relief shows through.
    float has = smoothstep(0.012, 0.05, max(img.r, max(img.g, img.b)));
    colour = mix(colour, img, has);
  }
  if (uDetailOn > 0.0) {
    float dl = mod(ll.x - uDetailBox.x + 360.0, 360.0);
    vec2 duv = vec2(dl / uDetailBox.z, (uDetailBox.y + uDetailBox.w - ll.y) / uDetailBox.w);
    if (duv.x >= 0.0 && duv.x <= 1.0 && duv.y >= 0.0 && duv.y <= 1.0) {
      vec3 d = texture(uDetail, duv).rgb;
      float edge = smoothstep(0.0, 0.04, min(min(duv.x, 1.0 - duv.x), min(duv.y, 1.0 - duv.y)));
      if (uDetailKind > 0.5) {
        // Fine shading on land: the relief map's own, coarser shading divided out
        // and the fine one applied, in the same form 1 + (s - 1) k (src/render/relief.js).
        float land = d.g * edge * uDetailOn;
        float fine = d.r / uShadeFlat;
        float coarse = 1.0 + (shade - 1.0) * 0.85;
        colour = mix(colour, colour / max(coarse, 0.25) * (1.0 + (fine - 1.0) * 0.9), land);
        shadeData = mix(shade, fine, land);
      } else {
        float has = smoothstep(0.012, 0.05, max(d.r, max(d.g, d.b)));
        colour = mix(colour, d, edge * has * uDetailOn);
      }
    }
  }

  float alt = asin(clamp(dot(earth, uSun), -1.0, 1.0)) * (180.0 / PI);
  float dark = darkness(alt) * uNight;
  // Night darkens the map fully, and data drawn over it by half, so that
  // values stay readable on the night side.
  colour *= mix(vec3(1.0), uNightTint, dark);
  vec3 dataNight = mix(vec3(1.0), uNightTint, dark * 0.5);
  if (uLightsOn > 0.5 && dark > 0.05) {
    vec3 lights = seamless(uLights, uv).rgb;
    float lum = dot(lights, vec3(0.30, 0.55, 0.15));
    colour += vec3(1.0, 0.80, 0.52) * pow(lum, 1.35) * 1.5 * smoothstep(0.2, 0.75, dark / 1.0);
  }

  if (uCloudOn > 0.0) {
    float c = atTime(uC0, uC1, uMixC, ll, uGridC, uGridNC);
    if (c > MISSING) {
      float a = pow(smoothstep(15.0, 100.0, c), 1.4) * 0.66 * uCloudOn;
      vec3 cloud = vec3(0.97, 0.96, 0.99) * mix(vec3(1.0), uNightTint * 0.8, dark);
      colour = mix(colour, cloud, a);
    }
  }

  if (uMode == 1 || uMode == 2 || uMode == 4) {
    float v = atTime(uA0, uA1, uMix, ll, uGrid, uGridN);
    if (uMode == 2 && v > MISSING) {
      float w = atTime(uB0, uB1, uMix, ll, uGrid, uGridN);
      v = w > MISSING ? length(vec2(v, w)) : MISSING;
    } else if (uMode == 4 && v > MISSING) {
      float w = atTime(uB0, uB1, uMixB, ll, uGridB, uGridNB);
      v = w > MISSING ? v - w : MISSING;
    }
    if (v > MISSING) {
      vec4 c = colourOf(v);
      vec3 shaded = c.rgb * mix(1.0, shadeData, uReliefOnData) * dataNight;
      colour = mix(colour, shaded, c.a * uOpacity);
    }
  } else if (uMode == 5 || uMode == 6) {
    uint k = classAt(uK0, ll, uGridK, uGridNK);
    if (k > 0u) {
      vec4 c = texelFetch(uClassColours, ivec2(int(k), 0), 0);
      if (uMode == 6 && classAt(uK1, ll, uGridKB, uGridNKB) == k) {
        // Unchanged: a quiet wash of the class, so that the changes stand out.
        c.rgb = mix(c.rgb, vec3(dot(c.rgb, vec3(0.30, 0.55, 0.15))), 0.7);
        c.a *= 0.3;
      }
      vec3 shaded = c.rgb * mix(1.0, shadeData, uReliefOnData) * dataNight;
      colour = mix(colour, shaded, c.a * uOpacity);
    }
  } else if (uMode == 3) {
    ivec2 size = textureSize(uCountries, 0);
    ivec2 px = ivec2(clamp(uv * vec2(size), vec2(0.0), vec2(size) - 1.0));
    uint id = texelFetch(uCountries, px, 0).r;
    if (id > 0u) {
      vec4 c = texelFetch(uCountryColours, ivec2(int(id), 0), 0);
      vec3 shaded = c.rgb * mix(1.0, shadeData, uReliefOnData) * dataNight;
      colour = mix(colour, shaded, c.a * uOpacity);
    }
  }
  fragColor = vec4(colour, 1.0);
}`;

function floatTexture(width, height, data) {
  const t = new THREE.DataTexture(data, width, height, THREE.RedFormat, THREE.FloatType);
  t.minFilter = t.magFilter = THREE.NearestFilter;
  t.generateMipmaps = false;
  t.needsUpdate = true;
  return t;
}

const empty = () => floatTexture(1, 1, new Float32Array([GPU_MISSING]));
const blank = () => {
  const t = new THREE.DataTexture(new Uint8Array([0, 0, 0, 255]), 1, 1, THREE.RGBAFormat);
  t.needsUpdate = true;
  return t;
};

export function createSurfaceMaterial() {
  const lut = new THREE.DataTexture(new Uint8Array(256 * 4), 256, 1, THREE.RGBAFormat);
  lut.minFilter = lut.magFilter = THREE.LinearFilter;
  lut.needsUpdate = true;
  const countries = new THREE.DataTexture(new Uint16Array(1), 1, 1, THREE.RedIntegerFormat, THREE.UnsignedShortType);
  countries.internalFormat = 'R16UI';
  countries.needsUpdate = true;
  const colours = new THREE.DataTexture(new Uint8Array(4), 1, 1, THREE.RGBAFormat);
  colours.needsUpdate = true;
  const base = new THREE.DataTexture(new Uint8Array([60, 50, 70, 128]), 1, 1, THREE.RGBAFormat);
  base.needsUpdate = true;
  const noClasses = () => classTexture({ values: new Uint8Array(4), grid: { nlon: 4, nlat: 1 } });
  const classColours = new THREE.DataTexture(new Uint8Array(256 * 4), 256, 1, THREE.RGBAFormat);
  classColours.needsUpdate = true;
  return new THREE.ShaderMaterial({
    glslVersion: THREE.GLSL3,
    vertexShader,
    fragmentShader,
    uniforms: {
      uProjection: { value: 0 },
      uLambda0: { value: 0 },
      uExtent: { value: new THREE.Vector2(Math.PI, Math.PI / 2) },
      uBase: { value: base },
      uImagery: { value: blank() }, uImageryOn: { value: 0 },
      uLights: { value: blank() }, uLightsOn: { value: 0 },
      uDetail: { value: blank() }, uDetailBox: { value: new THREE.Vector4(0, 0, 1, 1) }, uDetailOn: { value: 0 },
      uDetailKind: { value: 0 }, uShadeFlat: { value: 144 / 255 },
      uMode: { value: 0 },
      uA0: { value: empty() }, uA1: { value: empty() }, uB0: { value: empty() }, uB1: { value: empty() },
      uMix: { value: 0 }, uMixB: { value: 0 },
      uGrid: { value: new THREE.Vector4(0, 1, 0, 1) }, uGridB: { value: new THREE.Vector4(0, 1, 0, 1) },
      uGridN: { value: new THREE.Vector3(1, 1, 0) }, uGridNB: { value: new THREE.Vector3(1, 1, 0) },
      uLut: { value: lut },
      uRange: { value: new THREE.Vector3(0, 1, 0) },
      uOpacity: { value: 0.88 },
      uReliefOnData: { value: 0.55 },
      uCountries: { value: countries },
      uCountryColours: { value: colours },
      uK0: { value: noClasses() }, uK1: { value: noClasses() },
      uGridK: { value: new THREE.Vector4(0, 1, 0, 1) }, uGridKB: { value: new THREE.Vector4(0, 1, 0, 1) },
      uGridNK: { value: new THREE.Vector3(1, 1, 0) }, uGridNKB: { value: new THREE.Vector3(1, 1, 0) },
      uClassColours: { value: classColours },
      uC0: { value: empty() }, uC1: { value: empty() }, uMixC: { value: 0 }, uCloudOn: { value: 0 },
      uGridC: { value: new THREE.Vector4(0, 1, 0, 1) }, uGridNC: { value: new THREE.Vector3(1, 1, 0) },
      uSun: { value: new THREE.Vector3(1, 0, 0) },
      uNight: { value: 1 },
      uNightTint: { value: new THREE.Color(0.30, 0.27, 0.47) },
    },
  });
}

export function setBase(material, rgba, width, height) {
  const t = new THREE.DataTexture(rgba, width, height, THREE.RGBAFormat);
  t.wrapS = THREE.RepeatWrapping;
  t.wrapT = THREE.ClampToEdgeWrapping;
  t.minFilter = THREE.LinearMipmapLinearFilter;
  t.magFilter = THREE.LinearFilter;
  t.generateMipmaps = true;
  t.anisotropy = 4;
  t.colorSpace = THREE.NoColorSpace;
  t.needsUpdate = true;
  const old = material.uniforms.uBase.value;
  material.uniforms.uBase.value = t;
  old.dispose();
}

// An image (HTMLImageElement or canvas) as a global equirectangular texture.
// Its first row is the north, as in the field textures; three.js would flip an
// image on upload unless told not to, mirroring it north to south.
export function imageTexture(image, renderer, global = true) {
  const t = new THREE.Texture(image);
  t.flipY = false;
  t.wrapS = global ? THREE.RepeatWrapping : THREE.ClampToEdgeWrapping;
  t.wrapT = THREE.ClampToEdgeWrapping;
  t.minFilter = THREE.LinearMipmapLinearFilter;
  t.magFilter = THREE.LinearFilter;
  t.generateMipmaps = true;
  t.anisotropy = renderer ? Math.min(8, renderer.capabilities.getMaxAnisotropy()) : 4;
  t.colorSpace = THREE.NoColorSpace;
  t.needsUpdate = true;
  return t;
}

// Sets a sampler uniform, disposing of the texture it held.
export function setTexture(material, name, texture) {
  const u = material.uniforms[name];
  if (u.value && u.value !== texture && u.value.dispose) u.value.dispose();
  u.value = texture || blank();
}

export function setLut(material, table) {
  const t = material.uniforms.uLut.value;
  t.image.data.set(table);
  t.needsUpdate = true;
}

export function setCountryRaster(material, ids, width, height) {
  const t = new THREE.DataTexture(ids, width, height, THREE.RedIntegerFormat, THREE.UnsignedShortType);
  t.internalFormat = 'R16UI';
  t.minFilter = t.magFilter = THREE.NearestFilter;
  t.generateMipmaps = false;
  t.needsUpdate = true;
  material.uniforms.uCountries.value.dispose();
  material.uniforms.uCountries.value = t;
}

// colours: Uint8Array of RGBA for country index 0..n (entry 0 unused).
export function setCountryColours(material, colours) {
  const n = colours.length / 4;
  const t = new THREE.DataTexture(colours, n, 1, THREE.RGBAFormat);
  t.minFilter = t.magFilter = THREE.NearestFilter;
  t.needsUpdate = true;
  material.uniforms.uCountryColours.value.dispose();
  material.uniforms.uCountryColours.value = t;
}

// A class map as an integer texture (one byte per cell).
export function classTexture(raster) {
  const g = raster.grid;
  const t = new THREE.DataTexture(raster.values, g.nlon, g.nlat, THREE.RedIntegerFormat, THREE.UnsignedByteType);
  t.internalFormat = 'R8UI';
  t.minFilter = t.magFilter = THREE.NearestFilter;
  t.generateMipmaps = false;
  t.unpackAlignment = 1;
  t.needsUpdate = true;
  return t;
}

// Points the material at a class map (mode 5), or at the cells where it differs
// from another (mode 6). colours: Uint8Array of 256 RGBA entries, by class value.
export function setClasses(material, texture, raster, colours, other = null) {
  const u = material.uniforms;
  u.uK0.value = texture;
  u.uGridK.value.set(raster.grid.lon0, raster.grid.dlon, raster.grid.lat0, raster.grid.dlat);
  u.uGridNK.value.set(raster.grid.nlon, raster.grid.nlat, raster.grid.periodic ? 1 : 0);
  if (other) {
    u.uK1.value = other.texture;
    const g = other.raster.grid;
    u.uGridKB.value.set(g.lon0, g.dlon, g.lat0, g.dlat);
    u.uGridNKB.value.set(g.nlon, g.nlat, g.periodic ? 1 : 0);
  }
  u.uClassColours.value.image.data.set(colours);
  u.uClassColours.value.needsUpdate = true;
  u.uMode.value = other ? 6 : 5;
}

// A cache of field textures, keyed by product, variable and time index.
export class FieldTextures {
  constructor(limit = 32) {
    this.limit = limit;
    this.map = new Map();
  }

  get(product, name, index) {
    const key = `${product.meta.name}/${name}/${index}`;
    let t = this.map.get(key);
    if (t) { this.map.delete(key); this.map.set(key, t); return t; }
    const values = product.slice(name, Math.max(0, index));
    const data = new Float32Array(values.length);
    for (let i = 0; i < values.length; i++) data[i] = values[i] === values[i] ? values[i] : GPU_MISSING;
    t = floatTexture(product.grid.nlon, product.grid.nlat, data);
    this.map.set(key, t);
    if (this.map.size > this.limit) {
      const [k, old] = this.map.entries().next().value;
      this.map.delete(k);
      old.dispose();
    }
    return t;
  }

  // Forgets the textures of a product whose values have changed.
  forget(productName) {
    for (const [k, t] of this.map) if (k.startsWith(`${productName}/`)) { t.dispose(); this.map.delete(k); }
  }
}

function gridUniforms(u, suffix, product) {
  const g = product.grid;
  u[`uGrid${suffix}`].value.set(g.lon0, g.dlon, g.lat0, g.dlat);
  u[`uGridN${suffix}`].value.set(g.nlon, g.nlat, g.periodic ? 1 : 0);
}

// Points the material at a field: one variable (mode 1), the two components of
// a vector (mode 2), or the difference of two variables of two products (mode 4).
export function setField(material, textures, product, names, bracket, second = null) {
  const u = material.uniforms;
  const i0 = Math.max(0, bracket.i0), i1 = Math.max(0, bracket.i1);
  u.uA0.value = textures.get(product, names[0], i0);
  u.uA1.value = textures.get(product, names[0], i1);
  u.uMix.value = bracket.w || 0;
  gridUniforms(u, '', product);
  if (second) {
    const j0 = Math.max(0, second.bracket.i0), j1 = Math.max(0, second.bracket.i1);
    u.uB0.value = textures.get(second.product, second.name, j0);
    u.uB1.value = textures.get(second.product, second.name, j1);
    u.uMixB.value = second.bracket.w || 0;
    gridUniforms(u, 'B', second.product);
    u.uMode.value = 4;
  } else if (names[1]) {
    u.uB0.value = textures.get(product, names[1], i0);
    u.uB1.value = textures.get(product, names[1], i1);
    u.uMode.value = 2;
  } else {
    u.uMode.value = 1;
  }
}

export function setClouds(material, textures, product, name, bracket, on) {
  const u = material.uniforms;
  u.uCloudOn.value = on ? 1 : 0;
  if (!on || !product) return;
  u.uC0.value = textures.get(product, name, Math.max(0, bracket.i0));
  u.uC1.value = textures.get(product, name, Math.max(0, bracket.i1));
  u.uMixC.value = bracket.w || 0;
  gridUniforms(u, 'C', product);
}
