// The server's interface (serve.py; the protocol is the README's §3): the engine's state, the
// canons, a canon built into its complex with its invariants and readings, the network of the
// moment at a date, the check of a canon against Wikidata, and the recipe that builds the engine.
//
// With ?fixtures in the address the page answers from tests/fixtures/ instead: the engine's own
// answers for the first canon, and its moments at two dates, so that the page can be seen and
// developed without the engine.

const FIXTURES = typeof location !== 'undefined' && new URLSearchParams(location.search).has('fixtures');
export const fixtures = FIXTURES;
// The dates of the moments the fixtures hold.
export const FIXTURE_DATES = ['1808-10-01', '1814-11-15'];

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

export const engineState = () => (FIXTURES
  ? Promise.resolve({ built: true, current: true, toolchain: 'fixtures', fixtures: true })
  : fetch('/api/engine').then(json));
export const canons = () => (FIXTURES ? fixture('canons') : fetch('/api/canons').then(json));
export const setup = () => fetch('/setup').then(json);
export const buildEngine = () => fetch('/setup/run?again=1', { method: 'POST' }).then(json);

// A canon built: its record, complex, invariants and readings.
export const build = (canon) => (FIXTURES ? fixture(`build-${canon}`) : post('/api/build', { canon }));

// The power of each being over the reconstruction (the engine's power request).
export const power = (canon) => (FIXTURES ? fixture(`power-${canon}`) : post('/api/power', { canon }));

// The network of the moment at a date. The fixtures hold two dates; any other is answered with the
// nearer of them, and says so.
export async function moment(canon, date) {
  if (!FIXTURES) return post('/api/moment', { canon, date });
  const t = Date.parse(date.length === 4 ? `${date}-12-31` : date.length === 7 ? `${date}-28` : date);
  const near = FIXTURE_DATES.reduce((a, b) => (Math.abs(Date.parse(a) - t) <= Math.abs(Date.parse(b) - t) ? a : b));
  const m = await fixture(`moment-${canon}-${near}`);
  return near === date ? m : { ...m, substituted: near };
}

// The check of a canon's beings against Wikidata, and a new one.
export const check = (canon) => (FIXTURES ? Promise.resolve({ status: 'none', fixtures: true }) : fetch(`/api/check?canon=${encodeURIComponent(canon)}`).then(json));
export const runCheck = (canon) => (FIXTURES ? Promise.reject(new Error('checking needs the server, not the fixtures')) : post('/api/check/run', { canon }));
