// Finer maps of a body seen from close by.
//
// When the focused body fills much of the view, the part of its surface in
// view is fetched from the tile service that its map product declares (NASA's
// Trek services for the Moon, Mars and Mercury), at about the screen's
// resolution, and drawn over the global map by the body's own shader. The
// region is found by casting rays from the screen onto the body and turning the
// points into body-fixed longitude and latitude (src/render/maps.js), so the
// detail follows the body's rotation like its map.
//
// The tiles are equirectangular: level z has 2^(z+1) columns and 2^z rows of
// 256-pixel tiles, from 180 W and 90 N, in east longitude.

import * as THREE from 'three';
import { tilesFor } from '../core/tiles.js';

const GRID = 8;                        // rays across the screen, each way

function loadImage(url) {
  return new Promise((resolve, reject) => {
    const img = new Image();
    img.crossOrigin = 'anonymous';
    img.onload = () => resolve(img);
    img.onerror = () => reject(new Error(url));
    img.src = url;
  });
}

// Lets a body's material draw a detail texture over its map, in a box of
// body-fixed longitude and latitude: [west, south, width, height] in degrees.
function patchMaterial(material, uniforms) {
  const previous = material.onBeforeCompile;
  material.onBeforeCompile = (shader, renderer) => {
    if (previous) previous(shader, renderer);
    Object.assign(shader.uniforms, uniforms);
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', '#include <common>\nuniform sampler2D uDetail;\nuniform vec4 uDetailBox;\nuniform float uDetailOn;')
      .replace('#include <map_fragment>', `
#ifdef USE_MAP
  vec4 sampledDiffuseColor = texture2D( map, vMapUv );
  if ( uDetailOn > 0.0 ) {
    // The map's u runs east from 180 W, its v north from 90 S.
    vec2 ll = vec2( vMapUv.x * 360.0 - 180.0, vMapUv.y * 180.0 - 90.0 );
    float dl = mod( ll.x - uDetailBox.x + 360.0, 360.0 );
    vec2 duv = vec2( dl / uDetailBox.z, ( ll.y - uDetailBox.y ) / uDetailBox.w );
    if ( duv.x >= 0.0 && duv.x <= 1.0 && duv.y >= 0.0 && duv.y <= 1.0 ) {
      vec4 d = texture2D( uDetail, duv );
      float edge = smoothstep( 0.0, 0.03, min( min( duv.x, 1.0 - duv.x ), min( duv.y, 1.0 - duv.y ) ) );
      sampledDiffuseColor.rgb = mix( sampledDiffuseColor.rgb, d.rgb, d.a * edge * uDetailOn );
    }
  }
  diffuseColor *= sampledDiffuseColor;
#endif`);
  };
  material.needsUpdate = true;
}

export class BodyDetail {
  constructor(renderer) {
    this.renderer = renderer;
    this.bodies = new Map();          // key -> { mesh, live, uniforms }
    this.raycaster = new THREE.Raycaster();
    this.pending = null;
    this.timer = 0;
    this.lastKey = '';
    this.credit = '';
  }

  attach(key, mesh, live) {
    const uniforms = {
      uDetail: { value: new THREE.Texture() },
      uDetailBox: { value: new THREE.Vector4(0, 0, 1, 1) },
      uDetailOn: { value: 0 },
    };
    patchMaterial(mesh.material, uniforms);
    this.bodies.set(key, { mesh, live, uniforms, target: 0, box: null });
  }

  // The region of the body in view, [west, south, width, height] in degrees,
  // and the screen's pixels per degree across it; or null.
  viewRegion(body, camera, width, height) {
    const mesh = body.mesh;
    const centre = new THREE.Vector3().setFromMatrixPosition(mesh.matrixWorld);
    const radius = mesh.scale.x;
    const d = camera.position.distanceTo(centre);
    if (d <= radius) return null;
    const angular = Math.asin(Math.min(1, radius / d));
    const pixels = (angular / THREE.MathUtils.degToRad(camera.fov / 2)) * (height / 2);
    if (pixels < 0.3 * height) return null;
    const points = [];
    const ndc = new THREE.Vector2(), local = new THREE.Vector3();
    for (let i = 0; i <= GRID; i++) {
      for (let j = 0; j <= GRID; j++) {
        ndc.set((2 * i) / GRID - 1, (2 * j) / GRID - 1);
        this.raycaster.setFromCamera(ndc, camera);
        const hit = this.raycaster.intersectObject(mesh, false)[0];
        if (!hit) continue;
        mesh.worldToLocal(local.copy(hit.point)).normalize();
        // Body-fixed axes from the sphere's local ones: (x, y, z)_body = (x, -z, y)_local.
        const lon = THREE.MathUtils.radToDeg(Math.atan2(-local.z, local.x));
        const lat = THREE.MathUtils.radToDeg(Math.asin(Math.max(-1, Math.min(1, local.y))));
        points.push([lon, lat]);
      }
    }
    if (points.length < 6) return null;
    const lon0 = points[Math.floor(points.length / 2)][0];
    let w = Infinity, e = -Infinity, s = Infinity, n = -Infinity;
    for (const [lon, lat] of points) {
      const dl = ((lon - lon0 + 540) % 360) - 180;
      w = Math.min(w, dl); e = Math.max(e, dl);
      s = Math.min(s, lat); n = Math.max(n, lat);
    }
    // Rays between the samples may reach a little further.
    const mx = (e - w) / GRID, my = (n - s) / GRID;
    let west = lon0 + w - mx, span = e - w + 2 * mx;
    s = Math.max(-90, s - my); n = Math.min(90, n + my);
    if (n > 85 || s < -85) { west = -180; span = 360; }
    west = ((west + 540) % 360) - 180;
    if (span < 360 && west + span > 180) {
      // The service does not wrap round 180 degrees: keep the larger side.
      if (180 - west >= west + span - 180) span = 180 - west;
      else { span = west + span - 180; west = -180; }
    }
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const ppd = (Math.min(width, 2 * pixels) * dpr) / Math.max(1e-6, Math.min(span, 180));
    return { box: [west, s, Math.min(span, 360), n - s], ppd };
  }

  // Called every frame; asks for detail once the view has settled.
  update(camera, width, height, focusKey) {
    for (const [key, body] of this.bodies) {
      if (key !== focusKey) body.target = 0;
    }
    const body = focusKey && this.bodies.get(focusKey);
    if (body) {
      const region = this.viewRegion(body, camera, width, height);
      if (!region) body.target = 0;
      else {
        const key = `${focusKey}:${region.box.map((x) => x.toFixed(1)).join(',')}:${Math.round(Math.log2(region.ppd))}`;
        if (key !== this.lastKey) {
          this.lastKey = key;
          clearTimeout(this.timer);
          this.timer = setTimeout(() => this.fetch(focusKey, body, region), 500);
        }
      }
    }
    let changing = false;
    for (const b of this.bodies.values()) {
      const u = b.uniforms.uDetailOn, t = b.target;
      if (Math.abs(u.value - t) > 0.01) { u.value += (t - u.value) * 0.12; changing = true; } else u.value = t;
    }
    return changing;
  }

  async fetch(key, body, { box, ppd }) {
    const live = body.live;
    const { z, c0, c1, r0, r1, box: covered } = tilesFor(box, ppd, live);
    // Already showing this region at this level or finer.
    if (body.box && body.box.z >= z && body.box.c0 <= c0 && body.box.c1 >= c1 && body.box.r0 <= r0 && body.box.r1 >= r1) {
      body.target = 1;
      return;
    }
    const ticket = (this.pending = {});
    const canvas = document.createElement('canvas');
    canvas.width = (c1 - c0 + 1) * live.tile;
    canvas.height = (r1 - r0 + 1) * live.tile;
    const g = canvas.getContext('2d');
    const jobs = [];
    for (let r = r0; r <= r1; r++) {
      for (let c = c0; c <= c1; c++) {
        const url = live.template.replace('{z}', z).replace('{row}', r).replace('{col}', c);
        jobs.push(loadImage(url).then((img) => g.drawImage(img, (c - c0) * live.tile, (r - r0) * live.tile)).catch(() => null));
      }
    }
    await Promise.all(jobs);
    if (this.pending !== ticket) return;
    const texture = new THREE.CanvasTexture(canvas);
    texture.colorSpace = THREE.SRGBColorSpace;
    texture.anisotropy = Math.min(8, this.renderer.capabilities.getMaxAnisotropy());
    const u = body.uniforms;
    if (u.uDetail.value) u.uDetail.value.dispose();
    u.uDetail.value = texture;
    u.uDetailBox.value.set(...covered);
    body.box = { z, c0, c1, r0, r1 };
    body.target = 1;
    this.credit = live.credit;
  }
}
