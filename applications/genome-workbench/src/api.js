// The server's interface (serve.py; the protocol is the README's §3): the engine and its
// workers, the genomes there are, the user's files, one operation, a job cut into shards, and the
// recipe. Every request to the engine is reported to the page's activity (core/activity.js).

import * as activity from './core/activity.js';

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

async function tracked(kind, path, op, params) {
  const e = activity.begin(kind, op, params);
  try {
    const r = await post(path, { ...params, op });
    activity.end(e, r);
    return r;
  } catch (error) {
    activity.end(e, null, error);
    throw error;
  }
}

export const engineState = () => fetch('/api/engine').then(json);
export const genomes = () => fetch('/api/genomes').then(json);
export const setup = () => fetch('/setup').then(json);
export const runSetup = (again = false) => fetch(`/setup/run${again ? '?again=1' : ''}`, { method: 'POST' }).then(json);
export const restart = () => post('/api/restart', {});

// One operation on one worker: run('genome.tracks', {genome, begin, end, bins}).
export const run = (op, params = {}) => tracked('run', '/api/run', op, params);
// A job cut into shards across the workers, answered with its result and its trace (_job).
export const job = (op, params = {}) => tracked('job', '/api/job', op, params);

// The user's files.
export const files = () => fetch('/api/files').then(json);
export const fileHead = (id) => fetch(`/api/files/${encodeURIComponent(id.replace(/^file:/, ''))}/head`).then(json);
export const fileText = (id) => fetch(`/api/files/${encodeURIComponent(id.replace(/^file:/, ''))}/text`).then(json);
export const removeFile = (id) => fetch(`/api/files/${encodeURIComponent(id.replace(/^file:/, ''))}`, { method: 'DELETE' }).then(json);
// A file sent to the library, with its progress (0 to 1) reported as it goes.
export function upload(file, onProgress = () => {}) {
  return new Promise((resolve, reject) => {
    const x = new XMLHttpRequest();
    x.open('POST', `/api/files?name=${encodeURIComponent(file.name)}`);
    x.upload.onprogress = (e) => { if (e.lengthComputable) onProgress(e.loaded / e.total); };
    x.onload = () => {
      let body;
      try { body = JSON.parse(x.responseText); } catch { body = { message: x.responseText.slice(0, 200) }; }
      if (x.status >= 200 && x.status < 300) resolve(body); else reject(new Error(body.message || x.statusText));
    };
    x.onerror = () => reject(new Error('the upload failed'));
    x.send(file);
  });
}

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
