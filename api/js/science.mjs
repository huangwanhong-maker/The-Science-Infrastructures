// The client of the API (api/README.md, api/openapi.json) for JavaScript: Node 18
// or later, and browsers, with nothing to install. It has the surface of the
// Python client (api/python/science/client.py), in JavaScript's own spelling.
//
//   import { Client } from './science.mjs';
//   const api = new Client('http://127.0.0.1:8800');          // or an atlas's server, which mounts the API
//   const gfs = await api.field('gfs-forecast');
//   await gfs.sample('tmp2m', 139.69, 35.69, { time: 'now+1d' });   // one place: a number (NaN where there is none)
//   const { values, grid } = await gfs.grid('tmp2m', { time: 0 });   // values: a Float32Array, rows north to south
//   await api.compute('astronomy', 'sunPosition', 2460000.5);      // a domain's function, run by the service
//   for await (const event of api.stream('garden')) console.log(event.readings);
//
// Streams are Server-Sent Events: through EventSource in a browser, which
// reconnects by itself, and through a reader of fetch's body in Node, which
// reconnects after a dropped connection or 45 s without a heartbeat, asking for
// what it missed by the last event's id. Times are UTC: a Date, milliseconds since
// 1970, ISO 8601 text, 'now' or 'now+6h'.

export const DEFAULT_URL = 'http://127.0.0.1:8800';

export class ApiError extends Error {
  constructor(status, message, url = '', job = null) {
    super(`${status ?? 'unreachable'}: ${message}${url ? ` (${url})` : ''}`);
    this.status = status;
    this.detail = message;
    this.url = url;
    this.job = job;
  }
}

const LITTLE_ENDIAN = new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;
const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
const enc = encodeURIComponent;

// A query value: dates and milliseconds as ISO 8601, lists comma-separated.
function text(value, key) {
  if (value instanceof Date) return value.toISOString();
  if (typeof value === 'number' && (key === 'time' || key === 'since' || key === 'until')) return new Date(value).toISOString();
  if (Array.isArray(value) || ArrayBuffer.isView(value)) return Array.from(value).join(',');
  return String(value);
}

const numbers = (list) => list.map((v) => (v === null ? NaN : v));

async function errorOf(response, url) {
  let message = response.statusText;
  try {
    const body = await response.text();
    try { message = JSON.parse(body).error || body; } catch { message = body || message; }
  } catch { /* no body */ }
  return new ApiError(response.status, message, url);
}

export class Client {
  constructor(url, { token, fetch: fetcher } = {}) {
    let base = (url || globalThis.process?.env?.SCIENCE_API_URL || DEFAULT_URL).replace(/\/+$/, '');
    if (base.endsWith('/api')) base += '/v1';
    else if (!base.endsWith('/api/v1')) base += '/api/v1';
    this.url = base;
    this.token = token ?? globalThis.process?.env?.SCIENCE_API_TOKEN ?? null;
    this.fetch = fetcher || globalThis.fetch.bind(globalThis);
  }

  address(path, query = {}) {
    const q = new URLSearchParams();
    for (const [k, v] of Object.entries(query)) if (v !== undefined && v !== null && v !== '') q.set(k, text(v, k));
    const s = q.toString();
    return `${this.url}${path}${s ? `?${s}` : ''}`;
  }

  // A request of the API, by its path under /api/v1; the response, or an ApiError.
  async request(method, path, { query, body, headers = {}, signal } = {}) {
    const url = this.address(path, query);
    const h = { Accept: 'application/json', ...headers };
    if (this.token) h.Authorization = `Bearer ${this.token}`;
    if (body !== undefined) h['Content-Type'] = 'application/json';
    let response;
    try {
      response = await this.fetch(url, { method, headers: h, body: body === undefined ? undefined : JSON.stringify(body), signal });
    } catch (error) {
      if (error.name === 'AbortError') throw error;
      throw new ApiError(null, `the service cannot be reached: ${error.cause?.code || error.message}`, url);
    }
    if (!response.ok) throw await errorOf(response, url);
    return response;
  }

  async json(method, path, options) {
    return (await this.request(method, path, options)).json();
  }

  get(path, query) { return this.json('GET', path, { query }); }

  // ---- the service and its catalogue ---------------------------------------------------
  info() { return this.get('/'); }
  catalog() { return this.get('/catalog'); }

  // ---- fields ----------------------------------------------------------------------------
  async fields() { return (await this.get('/fields')).fields; }
  async field(name) { return new Field(this, await this.get(`/fields/${enc(name)}`)); }

  // ---- observations and their stream ------------------------------------------------------
  observations(source, { since, until, bbox, limit } = {}) {
    return source ? this.get(`/observations/${enc(source)}`, { since, until, bbox, limit }) : this.get('/observations');
  }

  // Readings [{ station, lon, lat, time, values: { variable: value } }] for a sensor
  // source, defined first when `define` is { title, variables, stations? }.
  postObservations(source, readings, { define } = {}) {
    const body = { readings: readings.map((r) => (r.time instanceof Date ? { ...r, time: r.time.toISOString() } : r)) };
    if (define) body.define = define;
    return this.json('POST', `/observations/${enc(source)}`, { body });
  }

  // New readings as they arrive, as { id, source, readings }: an async iterator.
  // Stop it with `break`, or with options.signal (an AbortSignal).
  stream(source, options = {}) {
    const url = this.address(source ? `/observations/${enc(source)}/stream` : '/stream');
    if (typeof EventSource !== 'undefined' && !options.reader) return eventSourceStream(url, options);
    return readerStream(this.fetch, url, options);
  }

  // ---- events, tracks, regions, verification ------------------------------------------------
  async events(name, { since, until, bbox, limit } = {}) {
    return name ? this.get(`/events/${enc(name)}`, { since, until, bbox, limit }) : (await this.get('/events')).events;
  }
  async tracks(name) { return name ? this.get(`/tracks/${enc(name)}`) : (await this.get('/tracks')).tracks; }
  async regions(name) { return name ? this.get(`/regions/${enc(name)}`) : (await this.get('/regions')).regions; }
  async verification(name) { return name ? this.get(`/verification/${enc(name)}`) : (await this.get('/verification')).verification; }

  // ---- the compute bridge ---------------------------------------------------------------------
  async units(unit) { return unit ? this.get(`/compute/${enc(unit)}`) : (await this.get('/compute')).units; }

  // A pure exported function of a domain or foundation, run by the service. Typed
  // arrays are sent as lists; wrap one as { $typed: 'Float64Array', values } where the
  // function needs a typed array.
  async compute(unit, fn, ...args) {
    const body = { args: args.map((a) => (ArrayBuffer.isView(a) ? Array.from(a) : a)) };
    return (await this.json('POST', `/compute/${enc(unit)}/${enc(fn)}`, { body })).result;
  }

  // ---- plug-ins and jobs ---------------------------------------------------------------------
  async plugins() { return (await this.get('/plugins')).applications; }
  describePlugin(app, file) { return this.get('/plugins/describe', { app, file }); }
  job(id) { return this.get(`/jobs/${enc(id)}`); }

  // Runs a server plug-in and returns its job; with `wait` (the default) once it is
  // done. A failed run throws an ApiError carrying the job.
  async runPlugin(app, file, { chosen = {}, params = {}, context = {}, wait = true, poll = 1000 } = {}) {
    let job = await this.json('POST', '/plugins/run', { body: { app, file, chosen, params, context } });
    if (!wait) return job;
    while (job.state === 'starting' || job.state === 'running') {
      await sleep(poll);
      job = await this.job(job.id);
    }
    if (job.state === 'failed') throw new ApiError(null, job.message || 'the run failed', '', job);
    return job;
  }
}

// A field product of a service: its description, and its values at places, over
// time and on its grid, computed by the service as the Earth Atlas computes them.
export class Field {
  constructor(client, meta) {
    this.client = client;
    this.meta = meta;
    this.name = meta.name;
    this.variables = new Map(meta.variables.map((v) => [v.name, v]));
    this.times = (meta.times || []).map((t) => new Date(t));
    this.path = `/fields/${enc(meta.name)}`;
  }

  // Values at places, bilinear in space and linear in time: a number for one place,
  // an array for many; NaN where there is none.
  async sample(variable, lon, lat, { time } = {}) {
    const many = Array.isArray(lon) || Array.isArray(lat) || ArrayBuffer.isView(lon) || ArrayBuffer.isView(lat);
    const r = await this.client.get(`${this.path}/sample`, { var: variable, lon, lat, time });
    const values = numbers(r.values);
    return many ? values : values[0];
  }

  // The value at a place at every time: { times: [Date], values: [number], units }.
  async series(variable, lon, lat) {
    const r = await this.client.get(`${this.path}/series`, { var: variable, lon, lat });
    return { times: r.times.map((t) => (t ? new Date(t) : null)), values: numbers(r.values), units: r.units, variable };
  }

  // The variable on its grid, or a box of it ([west, south, east, north]), every
  // stride-th point, at a stored time (an index) or an instant:
  // { values: Float32Array (row-major, rows north to south, NaN for none), grid, time, units }.
  async grid(variable, { time, bbox, stride } = {}) {
    const response = await this.client.request('GET', `${this.path}/grid`, {
      query: { var: variable, time, bbox, stride: stride === 1 ? undefined : stride, format: 'f32' },
      headers: { Accept: 'application/octet-stream' },
    });
    const described = JSON.parse(response.headers.get('X-Grid') || '{}');
    const buffer = await response.arrayBuffer();
    let values;
    if (LITTLE_ENDIAN) values = new Float32Array(buffer);
    else {
      const view = new DataView(buffer);
      values = new Float32Array(buffer.byteLength / 4);
      for (let i = 0; i < values.length; i++) values[i] = view.getFloat32(4 * i, true);
    }
    const g = described.grid;
    return { values, grid: g, time: described.time, units: described.units, variable, at: (r, c) => values[r * g.nlon + c] };
  }

  // The variable at a stored time as a file (GeoTIFF 'tif', netCDF 'nc', 'csv' or
  // 'json'): { data: ArrayBuffer, filename, contentType }.
  async export(variable, { time, format = 'tif', bbox, stride } = {}) {
    const response = await this.client.request('GET', `${this.path}/export`, {
      query: { var: variable, time, format, bbox, stride: stride === 1 ? undefined : stride }, headers: { Accept: '*/*' },
    });
    const m = /filename="?([^";]+)"?/.exec(response.headers.get('Content-Disposition') || '');
    return { data: await response.arrayBuffer(), filename: m ? m[1] : `${this.name}-${variable}.${format}`, contentType: response.headers.get('Content-Type') };
  }
}

// ---- Server-Sent Events -----------------------------------------------------------------------

// In a browser: EventSource, which reconnects by itself and sends Last-Event-ID.
async function* eventSourceStream(url, { signal } = {}) {
  const source = new EventSource(url);
  const queue = [];
  let wake = null, failure = null;
  const push = (item) => { queue.push(item); if (wake) { wake(); wake = null; } };
  source.addEventListener('readings', (e) => push({ id: Number(e.lastEventId) || e.lastEventId, ...JSON.parse(e.data) }));
  source.onerror = () => {
    if (source.readyState === EventSource.CLOSED) { failure = new ApiError(null, 'the service closed the stream', url); push(null); }
  };
  signal?.addEventListener('abort', () => push(null), { once: true });
  try {
    for (;;) {
      while (queue.length) {
        const item = queue.shift();
        if (item === null) { if (failure) throw failure; return; }
        yield item;
      }
      await new Promise((resolve) => { wake = resolve; });
    }
  } finally {
    source.close();
  }
}

// In Node (and anywhere fetch streams): a reader of the response's body, reconnecting.
async function* readerStream(fetcher, url, { lastEventId = null, signal, retry = 3000, silence = 45000, reconnect = true } = {}) {
  let last = lastEventId;
  while (!signal?.aborted) {
    const control = new AbortController();
    const stop = () => control.abort();
    signal?.addEventListener('abort', stop, { once: true });
    let timer = null;
    const watch = () => { clearTimeout(timer); timer = setTimeout(stop, silence); };
    try {
      const headers = { Accept: 'text/event-stream' };
      if (last !== null) headers['Last-Event-ID'] = String(last);
      let response;
      try {
        response = await fetcher(url, { headers, signal: control.signal });
      } catch (error) {
        if (signal?.aborted) return;
        if (!reconnect) throw new ApiError(null, `the service cannot be reached: ${error.cause?.code || error.message}`, url);
        response = null;
      }
      if (response && !response.ok) throw await errorOf(response, url);      // 404, 503: not worth retrying
      if (response) {
        watch();
        const reader = response.body.getReader();
        const decoder = new TextDecoder();
        let buffer = '', id = null, kind = 'message', data = [];
        for (;;) {
          let chunk;
          try { chunk = await reader.read(); } catch (error) { if (signal?.aborted) return; break; }
          if (chunk.done) break;
          watch();
          buffer += decoder.decode(chunk.value, { stream: true });
          let i;
          while ((i = buffer.indexOf('\n')) >= 0) {
            const line = buffer.slice(0, i).replace(/\r$/, '');
            buffer = buffer.slice(i + 1);
            if (line === '') {
              if (data.length && kind === 'readings') yield { id, ...JSON.parse(data.join('\n')) };
              kind = 'message'; data = [];
            } else if (!line.startsWith(':')) {
              const colon = line.indexOf(':');
              const name = colon < 0 ? line : line.slice(0, colon);
              const value = colon < 0 ? '' : line.slice(colon + 1).replace(/^ /, '');
              if (name === 'data') data.push(value);
              else if (name === 'event') kind = value;
              else if (name === 'id') { id = /^\d+$/.test(value) ? Number(value) : value; last = id; }
              else if (name === 'retry' && /^\d+$/.test(value)) retry = Number(value);
            }
          }
        }
      }
    } finally {
      clearTimeout(timer);
      signal?.removeEventListener('abort', stop);
      control.abort();
    }
    if (!reconnect || signal?.aborted) return;
    await sleep(retry);
  }
}
