// Checks of the shared visualization layer's provisional JavaScript: isolines and
// extrema, particles carried by a flow, the units in which quantities are shown,
// and the tiles of a map service.
// Run from the repository root:  node visualization/tests/run-tests.mjs

import { suite, seeded } from '../../development/tools/checks.mjs';
import {
  isolines, levels, extrema, ParticleSystem, uniformSpawn, displayUnit, toShown, toStored, format, tilesFor,
} from '../js/index.js';
import { DEG, EARTH_RADIUS_M } from '../../domains/geodesy/js/index.js';

const { check, assert, done } = suite();

// ---- contours --------------------------------------------------------------
{
  // The isoline of f = distance from (0, 0) is a circle (level chosen off the grid values).
  const grid = { lon0: -5, dlon: 0.1, nlon: 101, lat0: 5, dlat: 0.1, nlat: 101, periodic: false };
  const values = new Float32Array(grid.nlon * grid.nlat);
  for (let r = 0; r < grid.nlat; r++) for (let c = 0; c < grid.nlon; c++) values[r * grid.nlon + c] = Math.hypot(grid.lon0 + c * grid.dlon, grid.lat0 - r * grid.dlat);
  const R0 = 2.95;
  const seg = isolines(values, grid, R0);
  let worst = 0, length = 0;
  for (let i = 0; i < seg.length; i += 4) {
    worst = Math.max(worst, Math.abs(Math.hypot(seg[i], seg[i + 1]) - R0), Math.abs(Math.hypot(seg[i + 2], seg[i + 3]) - R0));
    length += Math.hypot(seg[i + 2] - seg[i], seg[i + 3] - seg[i + 1]);
  }
  check('contour: points of the isoline lie on the level', worst, 2e-3, 'deg');
  check('contour: length of the circle 2 pi r', length / (2 * Math.PI * R0) - 1, 1e-3);
  // Every segment end is shared by exactly one other segment: the line is closed.
  const ends = new Map();
  for (let i = 0; i < seg.length; i += 2) {
    const k = `${seg[i].toFixed(5)},${seg[i + 1].toFixed(5)}`;
    ends.set(k, (ends.get(k) || 0) + 1);
  }
  assert('contour: the isoline is closed', [...ends.values()].every((n) => n === 2), `${ends.size} vertices`);
  assert('contour: levels at multiples of a step', JSON.stringify(levels([3.2, 9.9], 2)) === '[4,6,8]');
  // A single high and a single low.
  const g = { lon0: -180, dlon: 2, nlon: 180, lat0: 89, dlat: 2, nlat: 90, periodic: true };
  const p = new Float32Array(g.nlon * g.nlat);
  for (let r = 0; r < g.nlat; r++) for (let c = 0; c < g.nlon; c++) {
    const lon = g.lon0 + c * g.dlon, lat = g.lat0 - r * g.dlat;
    p[r * g.nlon + c] = 101300 + 2000 * Math.exp(-((lon - 40) ** 2 + (lat - 45) ** 2) / 200) - 2500 * Math.exp(-((lon + 30) ** 2 + (lat + 51) ** 2) / 150);
  }
  const ex = extrema(p, g, { radius: 6, prominence: 300 });
  const high = ex.find((e) => e.kind === 'high'), low = ex.find((e) => e.kind === 'low');
  assert('contour: one high and one low found', ex.length === 2 && high && low, JSON.stringify(ex.map((e) => [e.kind, e.lon, e.lat])));
  if (high && low) check('contour: positions of the high and the low', Math.hypot(high.lon - 40, high.lat - 45) + Math.hypot(low.lon + 30, low.lat + 51), 1e-9, 'deg');
}

// ---- particles -------------------------------------------------------------
{
  // A uniform eastward flow of 10 m/s moves a particle on the equator by 10 dt metres.
  const ps = new ParticleSystem(1, { minLife: 1000, maxLife: 1001, random: seeded(1) });
  ps.reset(() => [0, 0]);
  ps.lon[0] = 0; ps.lat[0] = 0; ps.age[0] = 0;
  for (let i = 0; i < 100; i++) ps.step(() => [10, 0], 3600, () => [0, 0]);
  check('flow: uniform flow along the equator, 100 h', ps.lon[0] - (10 * 3600 * 100) / EARTH_RADIUS_M / DEG, 1e-9, 'deg');
  // Solid-body rotation about the pole: a particle keeps its latitude (the midpoint rule).
  const omega = 1e-5;   // rad/s
  const ring = new ParticleSystem(1, { minLife: 100000, maxLife: 100001, random: seeded(2) });
  ring.reset(() => [0, 60]);
  ring.lon[0] = 0; ring.lat[0] = 60; ring.age[0] = 0;
  for (let i = 0; i < 2000; i++) ring.step((lon, lat) => [omega * EARTH_RADIUS_M * Math.cos(lat * DEG), 0], 600, () => [0, 60]);
  check('flow: solid-body rotation keeps the latitude (2000 steps)', ring.lat[0] - 60, 0.05, 'deg');
  // A particle dies where the field has no value.
  const dead = new ParticleSystem(1, { random: seeded(3) });
  dead.reset(() => [5, 5]);
  dead.step(() => null, 60, () => [7, 7]);
  assert('flow: rebirth where there is no flow', dead.alive[0] === 0 && dead.lon[0] === 7);
  // Uniform spawning: the mean of sin(latitude) is near 0.
  const spawn = uniformSpawn();
  let s = 0;
  const rnd = seeded(4);
  for (let i = 0; i < 20000; i++) s += Math.sin(spawn(rnd)[1] * DEG);
  check('flow: spawning uniform in area', s / 20000, 0.01);
}

// ---- units ------------------------------------------------------------------
{
  const k = displayUnit('K');
  check('units: 300 K in degrees Celsius', toShown(300, k) - 26.85, 1e-9);
  check('units: back to kelvin', toStored(26.85, k) - 300, 1e-9);
  check('units: 1e-3 kg m-2 s-1 is 3.6 mm/h', toShown(1e-3, displayUnit('kg m-2 s-1')) - 3.6, 1e-12);
  check('units: 10 m/s in knots', toShown(10, displayUnit('m s-1', 'kn')) - 19.4384, 1e-4);
  assert('units: numbers for reading', format(124516650) === '125 M' && format(3.14159) === '3.14' && format(1013.2) === '1013');
}

// ---- tiles ------------------------------------------------------------------
// Tiles of the Trek map services, for detail at close range: the tiles chosen
// cover the region asked for, at a level no coarser than needed (unless the
// count of tiles forbids it), and never more than the limit.
{
  const live = { tile: 256, levels: 8 };
  let uncovered = 0, tooMany = 0, coarse = 0;
  let seed = 9;
  const rnd = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296);
  for (let k = 0; k < 2000; k++) {
    const width = 1 + 200 * rnd() ** 2, height = Math.min(170, 1 + 120 * rnd() ** 2);
    const west = -180 + (360 - Math.min(width, 359)) * rnd(), south = -90 + (180 - height) * rnd();
    const box = [west, south, Math.min(width, 180 - west), height];
    const ppd = 2 + 600 * rnd() ** 3;
    const t = tilesFor(box, ppd, live);
    const [w, s, cw, ch] = t.box;
    if (w > box[0] + 1e-9 || s > box[1] + 1e-9 || w + cw < box[0] + box[2] - 1e-9 || s + ch < box[1] + box[3] - 1e-9) uncovered++;
    if ((t.c1 - t.c0 + 1) * (t.r1 - t.r0 + 1) > 48 && t.z > 0) tooMany++;
    const needed = Math.min(live.levels, Math.max(0, Math.ceil(Math.log2((ppd * 360) / 512))));
    if (t.z < needed && tilesFor(box, ppd, live, Infinity).z === t.z) coarse++;
    if (t.z === needed && needed < live.levels && (2 ** (t.z + 1) * 256) / 360 < ppd * 0.999) coarse++;
  }
  check('tiles cover the region asked for (2000 regions)', uncovered, 0, 'regions');
  check('tiles: never more than 48', tooMany, 0, 'regions');
  check('tiles: the level is the one needed, or limited by the count', coarse, 0, 'regions');
  const t = tilesFor([-10, -5, 20, 10], 100, live);
  check('tiles: 20 x 10 degrees at 100 pixels per degree is level 6, 8 x 4 tiles',
    Math.abs(t.z - 6) + Math.abs(t.c1 - t.c0 + 1 - 8) + Math.abs(t.r1 - t.r0 + 1 - 4), 0, '');
}

// ---- isosurfaces by surface nets ----------------------------------------------------------
{
  const { isosurface, meshMeasures } = await import('../js/index.js');
  // A ball of radius 12 voxels of 2 mm: a closed surface, outward normals, the area and volume of a sphere of 24 mm.
  const n = 40, c = 19.5, r = 12, R = 2 * r;
  const data = new Uint8Array(n * n * n);
  for (let k = 0; k < n; k++) for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) if (Math.hypot(i - c, j - c, k - c) <= r) data[i + n * (j + n * k)] = 1;
  const affine = [2, 0, 0, -39, 0, 2, 0, -39, 0, 0, 2, -39, 0, 0, 0, 1];
  const m = isosurface({ dims: [n, n, n], data, affine }, 0.5, { smoothing: 2 });
  const mm = meshMeasures(m.positions, m.indices);
  check('isosurface: the area of a ball', Math.abs(mm.area / (4 * Math.PI * R * R) - 1), 0.02);
  check('isosurface: the volume enclosed, positive for outward faces', Math.abs(mm.volume / ((4 / 3) * Math.PI * R ** 3) - 1), 0.03);
  assert('isosurface: closed, with the Euler characteristic of a sphere', mm.euler === 2 && m.triangleCount === mm.faces, `${mm.euler}`);
  let inward = 0, far = 0;
  for (let v = 0; v < m.positions.length; v += 3) {
    const d = m.positions[v] * m.normals[v] + m.positions[v + 1] * m.normals[v + 1] + m.positions[v + 2] * m.normals[v + 2];
    if (d < 0) inward++;
    if (Math.abs(Math.hypot(m.positions[v], m.positions[v + 1], m.positions[v + 2]) - R) > 1.5) far++;
  }
  check('isosurface: every normal points outward', inward, 0, 'vertices');
  check('isosurface: every vertex within 1.5 mm of the sphere', far, 0, 'vertices');
  const empty = isosurface({ dims: [4, 4, 4], data: new Uint8Array(64), affine }, 0.5);
  assert('isosurface: nothing from an empty volume', empty.vertexCount === 0 && empty.triangleCount === 0);
}

// ---- cell complexes: words, facets, realisation, nets and Hasse diagrams ----------------------
{
  const { parseWords, parseFacets, checkComplex, eulerCharacteristic, cellCounts, realise, net, outerSides, hasse, cyclicColour, sequentialColour, divergingColour } = await import('../js/index.js');
  // Words: corners identified as the sides force.
  const torus = parseWords('abAB'), rp2 = parseWords('aa'), klein = parseWords('abAb'), disc = parseWords('abc'), sphere = parseWords('aA');
  assert('complex: abAB is one vertex, two loops and a square attached along a b a⁻¹ b⁻¹',
    torus.vertices === 1 && torus.edges.length === 2 && JSON.stringify(torus.faces[0].boundary) === '[1,2,-1,-2]');
  assert('complex: the Euler characteristics of words (torus 0, RP² 1, Klein 0, disc 1, sphere 2)',
    [torus, rp2, klein, disc, sphere].map(eulerCharacteristic).join() === '0,1,0,1,2');
  const genus = [0, 1, 2, 3, 4].map((g) => eulerCharacteristic(g ? parseWords(Array.from({ length: g }, (_, i) => `a${i}b${i}A${i}B${i}`).join('')) : sphere));
  assert('complex: χ = 2 − 2g for the standard words of genus 0 to 4', genus.join() === '2,0,-2,-4,-6', genus.join());
  assert('complex: two polygons with no letter in common have two vertices', parseWords('abAB cdCD').vertices === 2);
  assert('complex: "1" is a disc attached to a point', parseWords('1').faces[0].boundary.length === 0 && eulerCharacteristic(parseWords('1')) === 2);
  // Facets: every face, oriented; ∂∂ = 0 in every dimension.
  const sims = [1, 2, 3, 4, 5].map((n) => parseFacets(Array.from({ length: n + 1 }, (_, i) => i).join(' ')));
  const binom = (n, k) => { let r = 1; for (let i = 1; i <= k; i++) r = (r * (n - k + i)) / i; return r; };
  assert('complex: the n-simplex has C(n+1, k+1) cells of dimension k, n = 1 to 5',
    sims.every((c, n) => cellCounts(c).every((m, k) => m === binom(n + 2, k + 1))), sims.map((c) => cellCounts(c).join('/')).join('  '));
  assert('complex: every word of a simplex closes and ∂∂ = 0 in every dimension', sims.every((c) => checkComplex(c).length === 0), sims.map((c) => checkComplex(c).join('; ')).join(' | '));
  const tet = parseFacets('0 1 2; 0 1 3; 0 2 3; 1 2 3');
  assert('complex: the boundary of the tetrahedron, χ = 2, its triangles kept as simplices', eulerCharacteristic(tet) === 2 && tet.simplices.length === 4);
  assert('complex: faces of facets are not kept as simplices', parseFacets('0 1 2; 0 1; 2').simplices.length === 1);
  const t7 = parseFacets(Array.from({ length: 7 }, (_, i) => `${i} ${(i + 1) % 7} ${(i + 3) % 7}; ${i} ${(i + 2) % 7} ${(i + 3) % 7}`).join('; '));
  const uses = new Map();
  for (const f of t7.faces) for (const k of f.boundary) uses.set(Math.abs(k), (uses.get(Math.abs(k)) || 0) + 1);
  assert('complex: the 7-vertex torus has 21 edges, each in two triangles, χ = 0', t7.edges.length === 21 && [...uses.values()].every((u) => u === 2) && eulerCharacteristic(t7) === 0);
  assert('complex: a boundary that does not close is refused', checkComplex({ vertices: 2, edges: [[0, 1]], faces: [{ base: 0, boundary: [1] }] }).length === 1);
  assert('complex: a 3-cell whose boundary has a boundary is refused',
    checkComplex({ ...tet, cells: [[{ boundary: [[0, 1], [1, 1], [2, 1], [3, 1]] }]] }).some((p) => p.includes('boundary of its boundary')));
  // Realisation: deterministic, in the unit ball, loops kept apart, given positions respected.
  const g1 = realise(torus, { seed: 3 }), g2 = realise(torus, { seed: 3 });
  assert('complex: a realisation is the same from the same seed', JSON.stringify(g1.vertices) === JSON.stringify(g2.vertices) && JSON.stringify(g1.edges[1].points) === JSON.stringify(g2.edges[1].points));
  const inBall = (g) => [...g.vertices, ...g.edges.flatMap((e) => e.points)].every((p) => Math.hypot(...p) <= 1 + 1e-9) && g.faces.every((f) => { for (let i = 0; i < f.positions.length; i += 3) if (Math.hypot(f.positions[i], f.positions[i + 1], f.positions[i + 2]) > 1 + 1e-6) return false; return true; });
  assert('complex: every point of a realisation lies in the unit ball', inBall(g1) && inBall(realise(t7)));
  const a0 = g1.edges[0].arrow.at, a1 = g1.edges[1].arrow.at;
  check('complex: the two loops of the torus are drawn apart (shortfall of their midpoints\' distance below 0.2)', Math.max(0, 0.2 - Math.hypot(a0[0] - a1[0], a0[1] - a1[1], a0[2] - a1[2])), 0);
  assert('complex: a loop is a closed curve through its vertex', g1.edges.every((e) => Math.hypot(...e.points[0].map((x, i) => x - e.points.at(-1)[i])) < 1e-9 && Math.hypot(...e.points[0].map((x, i) => x - g1.vertices[0][i])) < 1e-9));
  const square = { vertices: 4, edges: [[0, 1], [1, 2], [2, 3], [3, 0]], faces: [{ base: 0, boundary: [1, 2, 3, 4] }], positions: [[0, 0, 0], [2, 0, 0], [2, 1, 0], [0, 1, 0]] };
  const gs = realise(square), d = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
  check('complex: given positions are kept up to scale (ratio of two sides, 2)', d(gs.vertices[0], gs.vertices[1]) / d(gs.vertices[1], gs.vertices[2]) - 2, 1e-9);
  const cap = realise(parseWords('1'));
  assert('complex: a disc attached to a point is a balloon through it', cap.faces[0].balloon && cap.faces[0].indices.length > 0);
  const counts = (g) => [g.vertices.length, g.edges.length, g.faces.length].join('/');
  const n30 = 30, fac = [];
  for (let i = 0; i < n30; i++) for (let j = 0; j < n30; j++) {
    const a = i * n30 + j, b = i * n30 + ((j + 1) % n30), c = ((i + 1) % n30) * n30 + j, e = ((i + 1) % n30) * n30 + ((j + 1) % n30);
    fac.push(`${a} ${b} ${e}`, `${a} ${e} ${c}`);
  }
  assert('complex: a large complex (a 30 × 30 triangulated torus) is laid out by its 1-skeleton, every cell drawn', counts(realise(parseFacets(fac.join(';')))) === `${n30 * n30}/${3 * n30 * n30}/${2 * n30 * n30}`);
  // Nets: unit sides; glued sides lie on one another, alike oriented; no two polygons overlap.
  const nets = [torus, klein, tet, t7, parseWords('abAB cdCD'), parseWords('abcdeABCDE')].map((c) => [c, net(c)]);
  let side = 0, glue = 0;
  for (const [, n] of nets) for (const p of n.polygons) for (const s of p.sides) {
    if (p.kind === 'polygon') side = Math.max(side, Math.abs(Math.hypot(s.to[0] - s.from[0], s.to[1] - s.from[1]) - 1));
    if (s.glued && s.mate) {
      const t = n.polygons[s.mate[0]].sides[s.mate[1]];
      glue = Math.max(glue, Math.hypot(s.tail[0] - t.tail[0], s.tail[1] - t.tail[1]), Math.hypot(s.head[0] - t.head[0], s.head[1] - t.head[1]));
    }
  }
  check('complex: every side of a net polygon has unit length', side, 1e-9);
  // Separating axes, written again here: two convex polygons of one piece share at most a boundary.
  const apart = (A, B) => [A, B].some((P) => P.some((p, i) => {
    const q = P[(i + 1) % P.length], nx = p[1] - q[1], ny = q[0] - p[0];
    const proj = (Q) => Q.map((r) => r[0] * nx + r[1] * ny);
    return Math.max(...proj(A)) <= Math.min(...proj(B)) + 1e-7 || Math.max(...proj(B)) <= Math.min(...proj(A)) + 1e-7;
  }));
  let overlaps = 0;
  for (const [, n] of nets) {
    const ps = n.polygons.filter((p) => p.kind === 'polygon');
    for (let i = 0; i < ps.length; i++) for (let j = i + 1; j < ps.length; j++) if (ps[i].piece === ps[j].piece && !apart(ps[i].points, ps[j].points)) overlaps++;
  }
  check('complex: no two polygons of a piece of a net overlap (pairs)', overlaps, 0);
  check('complex: glued sides coincide, tail on tail and head on head', glue, 1e-9);
  const outer = outerSides(nets[0][1]);
  assert('complex: the torus unfolds to its fundamental polygon a b a⁻¹ b⁻¹', nets[0][1].single && outer.length === 1 && outer[0].map((s) => s.sign * (s.edge + 1)).join() === '1,2,-1,-2');
  const tn = nets[2][1], glued = tn.polygons.flatMap((p) => p.sides).filter((s) => s.glued).length;
  assert('complex: the tetrahedron unfolds into one net, its outline one circle of six sides', tn.single && glued >= 6 && outerSides(tn).length === 1 && outerSides(tn)[0].length === 6);
  assert('complex: the 7-vertex torus unfolds into one piece', nets[3][1].single, `${nets[3][1].pieces} pieces`);
  assert('complex: two separate polygons are two pieces', nets[4][1].pieces === 2);
  // Hasse diagrams: layers by dimension, links between neighbouring layers.
  const h = hasse(tet);
  assert('complex: the Hasse diagram of the tetrahedron has layers 4, 6, 4 and 2·6 + 3·4 links', h.layers.map((l) => l.length).join() === '4,6,4' && h.links.length === 24);
  assert('complex: incidences of multiplicity two are kept (the projective plane meets its edge twice)', hasse(rp2).links.some((l) => l.dim === 2 && l.coefficient === 2));
  // Colour maps.
  const same = (a, b) => a.every((x, i) => Math.abs(x - b[i]) < 1e-9);
  assert('complex: the cyclic colour map closes on itself', same(cyclicColour(0), cyclicColour(2 * Math.PI)) && same(cyclicColour(-1), cyclicColour(2 * Math.PI - 1)));
  assert('complex: colours stay in the unit cube', [0, 0.3, 0.7, 1].every((t) => [sequentialColour(t), divergingColour(2 * t - 1), cyclicColour(6 * t)].every((c) => c.every((x) => x >= 0 && x <= 1))));
}

// ---- coloured graphs and quantum simplices ----------------------------------------------------------
{
  const { colouredGraphComplex, colouredGraphLayout, embedTetrahedron, indexColour, checkComplex, eulerCharacteristic, cellCounts } = await import('../js/index.js');
  // The elementary melon of four colours: 2 vertices, 4 edges, 6 bigons (one per pair of colours).
  const melon = colouredGraphComplex({ colours: 4, p: 1, sigma: [[0], [0], [0], [0]] });
  assert('complex: the 2-complex of the 4-coloured melon has 2 vertices, 4 edges, 6 faces, χ = 4, and is whole',
    cellCounts(melon).join() === '2,4,6' && eulerCharacteristic(melon) === 4 && checkComplex(melon).length === 0);
  // A random 4-coloured graph: one face per cycle of sigma_j⁻¹ sigma_i, each a closed path alternating colours.
  const random = seeded(7);
  const perm = (p) => { const s = [...Array(p).keys()]; for (let i = p - 1; i > 0; i--) { const j = Math.floor(random() * (i + 1)); [s[i], s[j]] = [s[j], s[i]]; } return s; };
  let bad = 0;
  for (let t = 0; t < 30; t++) {
    const p = 1 + (t % 6), cg = { colours: 3 + (t % 3), p, sigma: [] };
    for (let c = 0; c < cg.colours; c++) cg.sigma.push(c ? perm(p) : [...Array(p).keys()]);
    const X = colouredGraphComplex(cg);
    // Every edge of colour c lies on (colours − 1) faces, once each.
    const uses = new Array(X.edges.length).fill(0);
    for (const f of X.faces) for (const k of f.boundary) uses[Math.abs(k) - 1]++;
    if (checkComplex(X).length || uses.some((u) => u !== cg.colours - 1)) bad++;
  }
  assert('complex: the 2-complexes of 30 random coloured graphs are whole, each edge on (colours − 1) faces', bad === 0, `${bad} wrong`);
  const L = colouredGraphLayout({ colours: 5, p: 4, sigma: [[0, 1, 2, 3], [1, 0, 3, 2], [0, 1, 2, 3], [2, 3, 0, 1], [0, 1, 3, 2]] });
  let near = Infinity;
  for (let i = 0; i < L.nodes.length; i++) for (let j = 0; j < i; j++) near = Math.min(near, Math.hypot(L.nodes[i].x - L.nodes[j].x, L.nodes[i].y - L.nodes[j].y));
  assert('complex: a coloured graph is drawn in the unit square with its vertices apart and its parallel edges bent apart',
    L.nodes.every((n) => n.x >= 0 && n.x <= 1 && n.y >= 0 && n.y <= 1) && near > 0.1 && L.edges.filter((e) => e.u === 0 && e.v === 4).map((e) => e.bend).join() === '-0.5,0,0.5', `nearest ${near}`);
  const regular = embedTetrahedron([1, 1, 1, 1, 1, 1]);
  const d = (a, b) => Math.hypot(...a.map((x, k) => x - b[k]));
  const pairs = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
  const lengths = [1.5, 2, 2.5, 1.8, 2.2, 1.7], T = embedTetrahedron(lengths);
  check('complex: a tetrahedron set in space from its edge lengths keeps them', Math.max(...pairs.map(([a, b], k) => Math.abs(d(T.positions[a], T.positions[b]) - lengths[k]))), 1e-9);
  check('complex: the regular tetrahedron of unit edge has volume 1/(6√2)', regular.volume - 1 / (6 * Math.SQRT2), 1e-12);
  assert('complex: lengths that admit no tetrahedron are flagged', !embedTetrahedron([1, 1, 1, 1, 1, 3]).euclidean && regular.euclidean);
  assert('complex: colour 0 of a Feynman graph is the propagators\' pale cream, the others warm hues', indexColour(0).every((x) => x > 0.8) && indexColour(1).join() !== indexColour(2).join());
}

done();
