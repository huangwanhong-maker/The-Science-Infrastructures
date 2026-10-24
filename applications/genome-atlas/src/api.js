// The server's interface (serve.py; the protocol is the README's §3): the engine and its
// workers, the genomes there are, one operation, a job cut into shards, and the recipe.

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

export const engineState = () => fetch('/api/engine').then(json);
export const genomes = () => fetch('/api/genomes').then(json);
export const setup = () => fetch('/setup').then(json);
export const runSetup = (again = false) => fetch(`/setup/run${again ? '?again=1' : ''}`, { method: 'POST' }).then(json);
export const restart = () => post('/api/restart', {});

// One operation on one worker: run('genome.tracks', {genome, begin, end, bins}).
export const run = (op, params = {}) => post('/api/run', { ...params, op });
// A job cut into shards across the workers, answered with its result and its trace (_job).
export const job = (op, params = {}) => post('/api/job', { ...params, op });

// Requests whose answers are kept, by their JSON, for views that ask the same thing often.
export function cached(limit = 48) {
  const store = new Map();
  return (op, params) => {
    const key = JSON.stringify([op, params]);
    if (store.has(key)) {
      const v = store.get(key);
      store.delete(key);
      store.set(key, v);
      return v;
    }
    const p = run(op, params).catch((e) => { store.delete(key); throw e; });
    store.set(key, p);
    while (store.size > limit) store.delete(store.keys().next().value);
    return p;
  };
}
