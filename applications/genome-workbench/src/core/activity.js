// The activity of the page: every request to the engine, while it runs and after, for the status
// bar and its tray. api.js reports each request here; listeners are told of every change.

const entries = [];
const listeners = new Set();
let next = 1;
const KEEP = 200;

export function begin(kind, op, params = {}) {
  const e = { id: next++, kind, op, about: describe(params), started: performance.now(), at: new Date(), ended: null, ok: null, error: null, worker: null, seconds: null };
  entries.unshift(e);
  if (entries.length > KEEP) entries.length = KEEP;
  emit();
  return e;
}

export function end(e, result, error = null) {
  e.ended = performance.now();
  e.ok = !error;
  e.error = error ? error.message : null;
  if (result && typeof result === 'object') {
    e.worker = result._worker || result._job?.trace?.[0]?.worker || null;
    e.seconds = result._seconds ?? result._job?.seconds ?? null;
    if (result._job) e.shards = result._job.shards;
  }
  emit();
}

export const list = () => entries;
export const running = () => entries.filter((e) => e.ended === null);
export function subscribe(f) { listeners.add(f); return () => listeners.delete(f); }
function emit() { for (const f of listeners) f(); }

// A short line of what a request asked: its genome and its first few parameters.
function describe(p) {
  const parts = [];
  if (p.genome) parts.push(String(p.genome).replace(/^ncbi:/, ''));
  for (const [k, v] of Object.entries(p)) {
    if (k === 'genome' || k === 'op' || v === null || v === undefined) continue;
    if (parts.length >= 4) break;
    const s = typeof v === 'string' ? (v.length > 18 ? `${v.slice(0, 16)}…` : v) : Array.isArray(v) ? `${v.length} items` : typeof v === 'object' ? '{…}' : String(v);
    parts.push(`${k} ${s}`);
  }
  return parts.join(' · ');
}
