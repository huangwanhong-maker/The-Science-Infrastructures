// The server's interface (serve.py): the engine's state, its catalogue of families, complexes
// built and analysed, gauge fields and state sums on them, and the recipe that builds the engine.
//
// With ?fixtures in the address the page answers from tests/fixtures/ instead: invented answers
// for a few complexes, so that the page can be seen and developed without the engine.

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
