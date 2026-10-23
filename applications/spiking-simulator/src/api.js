// The server's interface (serve.py): the engine's state, its catalogue of models, runs, and the
// recipe that builds the engine.

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

export const engineState = () => fetch('/api/engine').then(json);
export const catalogue = () => fetch('/api/models').then(json).then((b) => b.models);
export const setup = () => fetch('/setup').then(json);
export const build = () => fetch('/setup/run?again=1', { method: 'POST' }).then(json);
export const cancel = () => fetch('/api/cancel', { method: 'POST' }).then(json);

// Runs a specification; its results as the engine writes them (times in seconds, SI units).
export function run(spec, threads) {
  const body = JSON.stringify(threads ? { ...spec, threads } : spec);
  return fetch('/api/run', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body }).then(json);
}
