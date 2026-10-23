// The server's interface (serve.py): the engine's state, its catalogue of families, complexes
// built and analysed, gauge fields and state sums on them, the generation of spacetime by field
// theories (theories, expansions, large N, condensates, random geometries, dimensions, spin foams),
// and the recipe that builds the engine.
//
// With ?fixtures in the address the page answers from tests/fixtures/ instead: invented answers
// for a few complexes, so that the page can be seen and developed without the engine.

import { colouredGraphComplex } from '../../../visualization/js/index.js';

const FIXTURES = typeof location !== 'undefined' && new URLSearchParams(location.search).has('fixtures');
export const fixtures = FIXTURES;

async function json(response) {
  const text = await response.text();
  let body;
  try { body = JSON.parse(text); } catch { body = { message: text.slice(0, 300) || response.statusText }; }
  if (!response.ok) {
    const error = new Error(body.message || `${response.status} ${response.statusText}`);
    error.status = response.status;
    error.body = body;
    throw error;
  }
  return body;
}

const post = (path, body) => fetch(path, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) }).then(json);
const fixture = (name) => fetch(`tests/fixtures/${name}.json`).then(json);
const noFixture = (what) => () => { throw new Error(`the fixtures hold ${what} only for their own complexes (torus, surface, projective plane, lattice torus, boundary of a simplex)`); };

export const engineState = () => (FIXTURES
  ? Promise.resolve({ built: true, current: true, toolchain: 'fixtures', fixtures: true, limits: { cells: 250000, sweeps: 200000, r: 9 } })
  : fetch('/api/engine').then(json));
export const catalogue = () => (FIXTURES ? fixture('catalogue') : fetch('/api/catalogue').then(json));
export const setup = () => fetch('/setup').then(json);
export const buildEngine = () => fetch('/setup/run?again=1', { method: 'POST' }).then(json);
export const cancel = () => fetch('/api/cancel', { method: 'POST' }).then(json);

// A complex of a family, with its parameters.
export const build = (family, params) => (FIXTURES
  ? fixture(`build-${family}`).catch(() => fixture('build-torus'))
  : post('/api/build', { family, params }));
// Its invariants.
export const analyse = (complex) => (FIXTURES
  ? fixture(`analyse-${complex.name}`).catch(noFixture('analyses'))
  : post('/api/analyse', { complex }));
// A gauge field on it: {complex, group, beta, sweeps, thermalise, seed, update, loops}.
export const gauge = (request) => (FIXTURES
  ? fixture(`gauge-${request.complex.name}`).catch(noFixture('gauge runs'))
  : post('/api/gauge', request));
// The Turaev-Viro state sum of a triangulated 3-manifold at the root of unity e^{iπ/r}.
export const turaevViro = (complex, r) => (FIXTURES
  ? fixture(`turaev-viro-${complex.name}`).catch(noFixture('state sums'))
  : post('/api/turaev-viro', { complex, r }));

// ---- the generation of spacetime by field theories (physis decision 0006) ---------------------------

const noGenerationFixture = (what) => () => { throw new Error(`the fixtures hold ${what} only for a few theories (the matrix model, the quartic melonic tensors, the coloured simplicial model, Boulatov's GFT)`); };

// The theories, with their actions and bubbles.
export const theories = () => (FIXTURES ? fixture('theories') : fetch('/api/theories').then(json));
// The Feynman graphs of a theory at an order, with their N-powers and weights, and its series.
export const expand = (theory, params, order) => (FIXTURES
  ? fixture(`expand-${theory}-${order}`).catch(noGenerationFixture('expansions'))
  : post('/api/expand', { theory, params, order }));
// The 2-complex (or the pseudomanifold) of a coloured graph. Without the engine the page makes the
// 2-complex itself (the same cells, numbered as the engine numbers them); the pseudomanifold needs it.
export const graphComplex = (coloured_graph, view = '2-complex') => (FIXTURES
  ? (view === '2-complex' ? Promise.resolve({ complex: colouredGraphComplex(coloured_graph) }) : Promise.reject(new Error('the pseudomanifold of a graph needs the engine')))
  : post('/api/graph-complex', { coloured_graph, view }));
// The large-N phase of a theory: an eigenvalue density, a melonic resummation, a gap equation,
// or, for a group field theory, the offer of its condensate.
export const largeN = (theory, params) => (FIXTURES
  ? (theory.startsWith('gft-') ? Promise.resolve({ kind: 'condensate' }) : fixture(`large-n-${theory}`).catch(noGenerationFixture('large-N phases')))
  : post('/api/large-n', { theory, params }));
// A condensate's cosmology: {spins: [{two_j, Q, E}], G, phi: [min, max], points}.
export const condensate = (request) => (FIXTURES ? fixture('condensate') : post('/api/condensate', request));
// A random geometry of a large-N phase: {kind: 'quadrangulation' | 'melonic', size, D, seed}.
export const randomGeometry = (request) => (FIXTURES
  ? fixture(`random-geometry-${request.kind}`)
  : post('/api/random-geometry', request));
// The spectral and Hausdorff dimensions of a complex's 1-skeleton or dual graph.
export const dimensions = (complex, graph = 'skeleton', options = {}) => (FIXTURES
  ? fixture(`dimensions-${complex.name}`).catch(noFixture('dimensions'))
  : post('/api/dimensions', { complex, graph, steps: 200, starts: 64, seed: 1, ...options }));
// A spin foam on the dual of a closed triangulated 3-manifold (or a 4-coloured graph).
// (Without the engine, any 4-coloured graph is answered with the melonic graph's fixture.)
export const spinFoam = (request) => (FIXTURES
  ? fixture(`spin-foam-${request.complex.name}`).catch(() => (request.complex.coloured_graph?.colours === 4 ? fixture('spin-foam-melon') : noFixture('spin foams')()))
  : post('/api/spin-foam', request));
