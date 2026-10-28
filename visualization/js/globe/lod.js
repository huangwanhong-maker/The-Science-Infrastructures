// The level of detail of a globe: which tiles to draw each frame, and which to load.
//
// The tiles form a quadtree over the body (tiling.js). Each frame the tree is walked from its
// roots; a tile is refined into its four children when its geometric error, projected on the
// screen, exceeds a tolerance of a few pixels,
//   ε_screen = ε · h / (2 d tan(φ/2)),
// with ε the tile's geometric error (the size of its cells), h the viewport's height in pixels, d
// the distance from the camera to the tile and φ the vertical field of view (Cozzi and Ring 2011,
// 3D Engine Design for Virtual Globes, §12.2; Cesium's quadtree). A tile outside the view frustum
// or below the horizon (figure.js, `occluded`, of the tile's horizon culling point) is not drawn and
// not refined. Children are drawn
// only when all four are covered by ready tiles; until then their parent is drawn in their place,
// so that the globe never has holes and refines as its data arrive. A child whose terrain failed
// (none exists so deep, or its source did not answer) leaves its parent drawn until it is tried
// again. The walk does not stop at a tile not yet loaded: judged by its proxy bounds it goes on for a
// few levels, so that a view brought close at once asks for all the levels it needs together
// (Cesium's "skipping levels of detail", in spirit), the nearer to the tiles already drawn first.
//
// Loading is by priority: the tiles that most exceed their tolerance first, coarser before finer,
// at most a few at a time; a request whose tile has left the view is cancelled. The tiles kept in
// memory are bounded, the least recently drawn released first. Imagery is loaded per tile and
// layer; a tile without its own image draws the part of its nearest ancestor's (tiling.js,
// `subregion`), so that imagery sharpens in place.
//
// Nothing here draws: a renderer (renderer-three.js) is told when a tile becomes ready or is
// released, and is given each frame the tiles to draw with their images. The tree works in the
// body-fixed frame of its body.

import { cartesian, horizonCullingPoint, occluded } from './figure.js';
import { children, key, roots, subregion, tileBounds, cellDegrees, GRID } from './tiling.js';

const DEG = Math.PI / 180;

export class Tile {
  // heightRange: the heights the tile's terrain may span, for its proxy bounds until it is loaded;
  // inset: the depth of the body's lowest ground, for the horizon test.
  constructor(z, x, y, parent, fig, heightRange, inset = -heightRange[0] + 1) {
    this.z = z, this.x = x, this.y = y;
    this.key = key(z, x, y);
    this.parent = parent;
    this.children = null;
    this.state = 'new'; // new, loading, ready, failed
    this.mesh = null;   // the built mesh (terrain.js)
    this.gpu = null;    // what the renderer made of it
    this.images = new Map(); // layer id → { state, image, gpu }
    this.lastUsed = -1;
    this.abort = null;
    this.priority = 0;
    this.setProxy(fig, heightRange, inset);
    this.geometricError = cellDegrees(z, GRID) * DEG * fig.a;
  }

  // Proxy bounds before the terrain is known: the tile between the lowest and highest heights it may
  // have; for the horizon test, a 5 x 5 grid of its points at the highest, stood for by one point.
  setProxy(fig, heightRange, inset) {
    const { z, x, y } = this;
    const [west, south, east, north] = tileBounds(z, x, y);
    const pts = [];
    const grid = [0, 0.25, 0.5, 0.75, 1];
    for (const h of heightRange)
      for (const fx of grid) for (const fy of grid) pts.push(cartesian(fig, (west + (east - west) * fx) * DEG, (north - (north - south) * fy) * DEG, h));
    const c = [0, 0, 0];
    for (const p of pts) for (let k = 0; k < 3; k++) c[k] += p[k] / pts.length;
    let r = 0;
    for (const p of pts) r = Math.max(r, Math.hypot(p[0] - c[0], p[1] - c[1], p[2] - c[2]));
    const probes = pts.slice(25);
    // A tile spanning more than a hemisphere's width (the roots) bulges beyond its corners.
    this.proxy = {
      centre: c, radius: z === 0 ? fig.a * 1.05 : r * 1.02 + (east - west > 45 ? fig.a * 0.1 : 0), probes,
      occludee: z === 0 ? null : horizonCullingPoint(fig, probes, inset),
    };
    this.proxyRange = heightRange;
  }

  get bounds() {
    return this.mesh ? { centre: this.mesh.centre, radius: this.mesh.radius, probes: this.mesh.probes, occludee: this.mesh.occludee } : this.proxy;
  }
}

// The heights a tile not yet loaded may span: its nearest loaded ancestor's, widened by a quarter of
// their span and 100 m (finer data reach higher peaks and deeper valleys than a coarser grid's
// samples), within the body's; the body's when no ancestor is loaded.
function estimatedRange(parent, bodyRange) {
  for (let a = parent; a; a = a.parent) {
    if (a.mesh) {
      const margin = 0.25 * (a.mesh.max - a.mesh.min) + 100;
      return [Math.max(bodyRange[0], a.mesh.min - margin), Math.min(bodyRange[1], a.mesh.max + margin)];
    }
  }
  return bodyRange;
}

export class TileTree {
  // fig: the body's figure; maxLevel: the deepest level; heightRange: [lowest, highest] heights of
  // the body (m); loadTerrain(z, x, y, signal) → Promise of a built mesh; layers: [{ id, maxLevel,
  // load(z, x, y, signal) → Promise of an image }]; on: { ready(tile), released(tile),
  // imageReady(tile, layer), imageReleased(tile, layer) }.
  constructor({ fig, maxLevel, heightRange, loadTerrain, layers = [], on = {}, tolerance = 2.5, maxLoads = 6, maxImageLoads = 8, maxTiles = 700, falloff = 0.1 }) {
    Object.assign(this, { fig, maxLevel, heightRange, loadTerrain, layers, on, tolerance, maxLoads, maxImageLoads, maxTiles });
    this.retryFrames = 3600;
    this.skip = 6; // levels walked below the nearest ready tile
    this.walkBudget = 600; // tiles not yet loaded walked through in a frame
    this.stage = 3; // every third level on the way down is loaded as a stage
    this.falloff = falloff; // how fast detail falls off toward the horizon (0: not at all)
    this.tiles = new Map();
    this.frame = 0;
    this.loading = 0;
    this.imageLoading = 0;
    this.roots = roots().map(([z, x, y]) => this.tile(z, x, y, null));
    this.stats = { drawn: 0, ready: 0, loading: 0, culled: 0 };
  }

  tile(z, x, y, parent) {
    const k = key(z, x, y);
    let t = this.tiles.get(k);
    if (!t) {
      t = new Tile(z, x, y, parent, this.fig, estimatedRange(parent, this.heightRange), -this.heightRange[0] + 1);
      this.tiles.set(k, t);
    }
    return t;
  }

  // When a tile's terrain arrives, its descendants not yet loaded are bounded by it.
  narrow(t) {
    const stack = [...(t.children || [])];
    while (stack.length) {
      const c = stack.pop();
      if (c.mesh) continue;
      const range = estimatedRange(c.parent, this.heightRange);
      if (range[0] !== c.proxyRange[0] || range[1] !== c.proxyRange[1]) c.setProxy(this.fig, range, -this.heightRange[0] + 1);
      if (c.children) stack.push(...c.children);
    }
  }

  // The terrain's height as drawn (m, with the mesh's exaggeration) at a longitude and latitude
  // (radians), from the finest tile loaded there, bilinear; 0 where none is.
  heightAt(lon, lat) {
    const lonD = Math.atan2(Math.sin(lon), Math.cos(lon)) / DEG, latD = lat / DEG;
    let t = this.roots[lonD < 0 ? 0 : 1];
    let best = t.state === 'ready' ? t : null;
    while (t.children) {
      const [, , x0, y0] = [0, 0, 2 * t.x, 2 * t.y];
      const [w, s, e, n] = tileBounds(t.z, t.x, t.y);
      const cx = lonD >= (w + e) / 2 ? 1 : 0, cy = latD <= (s + n) / 2 ? 1 : 0;
      const c = t.children.find((k) => k.x === x0 + cx && k.y === y0 + cy);
      if (!c || c.state !== 'ready') break;
      best = t = c;
    }
    if (!best || !best.mesh) return 0;
    const n = best.mesh.header?.n || 65;
    const [w, s, e, nn] = tileBounds(best.z, best.x, best.y);
    const fx = Math.min(n - 1.000001, Math.max(0, ((lonD - w) / (e - w)) * (n - 1)));
    const fy = Math.min(n - 1.000001, Math.max(0, ((nn - latD) / (nn - s)) * (n - 1)));
    const i = Math.floor(fx), j = Math.floor(fy), a = fx - i, b = fy - j;
    const h = best.mesh.height;
    const at = (ii, jj) => h[jj * n + ii];
    return ((1 - b) * ((1 - a) * at(i, j) + a * at(i + 1, j)) + b * ((1 - a) * at(i, j + 1) + a * at(i + 1, j + 1))) * (best.mesh.exaggeration ?? 1);
  }

  setLayers(layers) {
    // Release images of layers no longer shown.
    const ids = new Set(layers.map((l) => l.id));
    for (const t of this.tiles.values()) {
      for (const [id, img] of t.images) {
        if (!ids.has(id)) {
          if (img.abort) img.abort.abort();
          if (img.state === 'ready') this.on.imageReleased?.(t, id);
          t.images.delete(id);
        }
      }
    }
    this.layers = layers;
  }

  // view: { camera: [x, y, z] body-fixed (m), planes: [[nx, ny, nz, d], ...] body-fixed (n·p + d ≥ 0
  // inside), height (px), fovY (rad) }. Returns the tiles to draw, each with its images and the
  // quarters of it to draw (null for all of it).
  update(view) {
    this.frame++;
    const wanted = [];
    const st = { drawn: 0, culled: 0 };
    const visible = (t) => {
      const b = t.bounds;
      for (const [nx, ny, nz, d] of view.planes) if (nx * b.centre[0] + ny * b.centre[1] + nz * b.centre[2] + d < -b.radius) return false;
      if (b.occludee && occluded(this.fig, view.camera, b.occludee, -this.heightRange[0] + 1)) return false;
      return true;
    };
    // The camera's height above the figure: |p| (1 − 1/s), s the point's scaled norm (s = 1 on the figure).
    const [cx, cy, cz] = view.camera;
    const sNorm = Math.sqrt((cx * cx + cy * cy) / (this.fig.a * this.fig.a) + (cz * cz) / (this.fig.b * this.fig.b));
    const altitude = Math.max(50, Math.hypot(cx, cy, cz) * (1 - 1 / sNorm));
    // The screen-space error of a tile, lessened toward the horizon: a tile d away counts as if it were
    // 1 + falloff (d / altitude − 1) times farther, so that an oblique view does not spend its tiles
    // where haze and foreshortening hide them (Cesium's dynamic screen-space error, in spirit).
    const error = (t) => {
      const b = t.bounds;
      const d = Math.max(1, Math.hypot(b.centre[0] - cx, b.centre[1] - cy, b.centre[2] - cz) - b.radius);
      const fade = 1 + this.falloff * Math.max(0, d / altitude - 1);
      return (t.geometricError * view.height) / (2 * d * fade * Math.tan(view.fovY / 2));
    };
    const want = (t, priority) => {
      t.lastUsed = this.frame;
      if (t.state === 'new' || t.state === 'queued') {
        t.state = 'queued';
        t.priority = priority;
        wanted.push(t);
      }
    };
    // What covers t's part of the body at the detail the view asks for: a list of [tile, quarter]
    // (quarter −1 for the whole tile, else 0..3, the part of it its child of that index covers), or
    // null when t is not loaded and its part cannot yet be covered by loaded tiles below it. A loaded
    // tile always covers its part: by its children's where they can, and by its own quarters where
    // they cannot yet, so that each quarter refines as soon as its own tiles have come.
    // `since` counts the levels below the nearest loaded ancestor: the walk goes on through tiles not
    // yet loaded, judged by their proxy bounds, for `skip` levels, so that a view that comes close at
    // once asks for its tiles at once rather than one level after another.
    let walked = 0; // tiles not yet loaded walked through this frame
    const resolve = (t, since) => {
      t.lastUsed = this.frame;
      if (!visible(t)) { st.culled++; return []; }
      const ready = t.state === 'ready';
      const e = error(t);
      const refine = t.z < this.maxLevel && e > this.tolerance;
      const descend = refine && (ready || (since < this.skip && walked++ < this.walkBudget));
      if (descend) {
        if (!t.children) t.children = children(t.z, t.x, t.y).map(([z, x, y]) => this.tile(z, x, y, t));
        const parts = [];
        let covered = true;
        t.children.forEach((c, q) => {
          // A failure is tried again after a while (a source may have been away); until then the
          // parent stands in for it.
          if (c.state === 'failed' && this.frame - c.failedAt > this.retryFrames) c.state = 'new';
          let part = null;
          if (c.state === 'failed') { c.lastUsed = this.frame; part = visible(c) ? null : []; }
          else part = resolve(c, ready ? 1 : since + 1);
          if (part === null) {
            covered = false;
            if (ready) parts.push([[t, q]]);
          } else parts.push(part);
        });
        if (covered || ready) {
          if (ready) t.screenError = e;
          return parts.flat();
        }
      }
      if (ready) {
        t.screenError = e;
        return [[t, -1]];
      }
      // Asked for: the tiles where the walk ends (the view's detail), and every third level on the
      // way as a stage, so that a view brought close sharpens in a few steps; nearer the tiles
      // already drawn first, and among those the larger error first.
      if (!descend || since === 1 || since % this.stage === 0) want(t, 1e6 / since + Math.min(e, 1e5) * 1e-3);
      return null;
    };
    let cover = [];
    for (const r of this.roots) cover = cover.concat(resolve(r, 1) || []);
    // One entry for each tile drawn, with the quarters drawn (null: the whole tile).
    const byTile = new Map();
    for (const [t, q] of cover) {
      if (!byTile.has(t)) byTile.set(t, q < 0 ? null : []);
      const quarters = byTile.get(t);
      if (quarters && q >= 0) quarters.push(q);
    }
    const draw = [...byTile.keys()];
    st.drawn = draw.length;
    // Their images are asked for, the largest error first.
    for (const t of [...draw].sort((a, b) => b.screenError - a.screenError)) this.wantImages(t, t.screenError);
    // Start the most wanted loads.
    wanted.sort((a, b) => b.priority - a.priority);
    for (const t of wanted) {
      if (this.loading >= this.maxLoads) break;
      this.startLoad(t);
    }
    // Cancel loads that are no longer wanted.
    for (const t of this.tiles.values()) {
      if (t.state === 'loading' && t.lastUsed < this.frame - 30 && t.abort) {
        t.abort.abort();
        this.aborted = (this.aborted || 0) + 1;
      }
      if (t.state === 'queued' && t.lastUsed < this.frame) t.state = 'new';
    }
    this.evict();
    if (this.frame % 120 === 0) this.prune();
    // The images each drawn tile uses: its own, or the part of an ancestor's.
    const result = draw.map((t) => ({ tile: t, quarters: byTile.get(t), images: this.layers.map((l) => this.imageFor(t, l)) }));
    this.stats = { drawn: draw.length, culled: st.culled, loading: this.loading + this.imageLoading, tiles: this.tiles.size };
    return result;
  }

  startLoad(t) {
    t.state = 'loading';
    t.abort = new AbortController();
    this.loading++;
    this.loadTerrain(t.z, t.x, t.y, t.abort.signal).then((mesh) => {
      this.loading--;
      t.abort = null;
      if (!mesh) { t.state = 'new'; return; }
      mesh.occludee = t.z === 0 ? null : horizonCullingPoint(this.fig, mesh.probes, -this.heightRange[0] + 1);
      t.mesh = mesh;
      t.state = 'ready';
      this.narrow(t);
      this.on.ready?.(t);
    }, (err) => {
      this.loading--;
      t.abort = null;
      t.state = err && err.name === 'AbortError' ? 'new' : 'failed';
      if (t.state === 'failed') { t.error = String(err && err.message || err); t.failedAt = this.frame; }
    });
  }

  wantImages(t, e) {
    for (const layer of this.layers) {
      if (t.z > layer.maxLevel) continue; // drawn from an ancestor's image
      let img = t.images.get(layer.id);
      if (!img) t.images.set(layer.id, (img = { state: 'new' }));
      if (img.state === 'new' && this.imageLoading < this.maxImageLoads) {
        img.state = 'loading';
        img.abort = new AbortController();
        this.imageLoading++;
        layer.load(t.z, t.x, t.y, img.abort.signal).then((image) => {
          this.imageLoading--;
          img.abort = null;
          if (!image) { img.state = 'failed'; return; }
          img.state = 'ready';
          img.image = image;
          this.on.imageReady?.(t, layer.id, image);
        }, (err) => {
          this.imageLoading--;
          img.abort = null;
          img.state = err && err.name === 'AbortError' ? 'new' : 'failed';
        });
      }
    }
  }

  // The image to draw a tile with for a layer: [tile whose image, u0, v0, size], or null.
  imageFor(t, layer) {
    for (let a = t; a; a = a.parent) {
      if (a.z > layer.maxLevel) continue;
      const img = a.images.get(layer.id);
      if (img && img.state === 'ready') return { owner: a, region: subregion(t.z, t.x, t.y, a.z), layer: layer.id };
    }
    return null;
  }

  evict() {
    let count = 0;
    for (const t of this.tiles.values()) if (t.state === 'ready') count++;
    if (count <= this.maxTiles) return;
    const ready = [...this.tiles.values()].filter((t) => t.state === 'ready' && t.z > 0);
    ready.sort((a, b) => a.lastUsed - b.lastUsed);
    const recent = (t) => t.lastUsed >= this.frame - 2 || (t.children || []).some(recent);
    for (const t of ready) {
      if (count <= this.maxTiles || t.lastUsed >= this.frame - 2) break;
      if (!this.tiles.has(t.key) || !t.parent) continue;
      const group = t.parent.children || [t];
      if (group.some(recent)) continue; // a sibling or a descendant is still drawn
      const before = this.tiles.size;
      this.release(t);
      count -= before - this.tiles.size;
    }
  }

  // Groups of siblings none of which has been drawn or walked for a while (ten seconds at sixty
  // frames a second) are released, ready or not, so that the tree does not keep what the view left.
  prune(stale = 600) {
    for (const t of [...this.tiles.values()]) {
      if (!t.children || !this.tiles.has(t.key)) continue;
      if (t.children.every((c) => c.lastUsed < this.frame - stale && c.state !== 'loading')) this.release(t.children[0]);
    }
  }

  // A tile is released with its siblings, so that a tile has either no children or all four.
  release(t) {
    if (!t.parent) return;
    const group = t.parent.children || [t];
    for (const s of group) this.releaseOne(s);
    t.parent.children = null;
  }
  releaseOne(t) {
    if (t.children) {
      for (const c of t.children) this.releaseOne(c);
      t.children = null;
    }
    if (t.abort) t.abort.abort();
    for (const [id, img] of t.images) {
      if (img.abort) img.abort.abort();
      if (img.state === 'ready') this.on.imageReleased?.(t, id);
    }
    t.images.clear();
    if (t.state === 'ready') this.on.released?.(t);
    t.mesh = null;
    t.state = 'new';
    this.tiles.delete(t.key);
  }
}
