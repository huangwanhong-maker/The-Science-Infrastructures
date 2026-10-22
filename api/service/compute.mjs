// The compute bridge: a persistent Node process through which the API calls the
// pure functions of the domains and foundations while they are JavaScript, so that
// Python, C and C++ reach the same algorithms the pages run (api/README.md §2.2).
//
// The callable units are the public interfaces the dependency rules recognise:
// domains/<d>/js/index.js and foundations/<f>/js/index.js. Only their exported
// functions are callable (classes are not), by name, with JSON arguments; nothing
// is evaluated. The units are pure (development/tools/check-dependencies.mjs
// checks it: no network, no page, no storage), so a call reads and changes
// nothing outside itself.
//
// The protocol is JSON lines: one request per line on standard input, one answer
// per line on standard output, in order.
//
//   {"id": 1, "op": "units"}                                          -> {"id": 1, "ok": true, "result": [unit, ...]}
//   {"id": 2, "op": "call", "unit": "astronomy", "function": "sunPosition", "args": [2460000.5]}
//                                                                     -> {"id": 2, "ok": true, "result": {...}}
//   an error                                                          -> {"id": n, "ok": false, "kind": "not-found" | "failed", "error": "..."}
//
// Arguments are JSON; {"$typed": "Float64Array", "values": [...]} makes a typed
// array (null for NaN). Results are made JSON: typed arrays and sets as lists,
// maps and objects as objects, dates as ISO 8601, NaN and infinities as null.
// The worker has no time limit of its own: its host (science/bridge.py) kills
// and restarts it when a call takes too long.
//
// Run:  node api/service/compute.mjs        (the host starts it; a person need not)

import { existsSync, readFileSync, readdirSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { createInterface } from 'node:readline';
import { fileURLToPath, pathToFileURL } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');   // import.meta.dirname would need Node 20.11
const TYPED = { Float64Array, Float32Array, Int32Array, Uint32Array, Int16Array, Uint16Array, Int8Array, Uint8Array, Uint8ClampedArray };

// The units: each domain and foundation with a public interface in JavaScript.
function discover() {
  const out = new Map();
  for (const [kind, dir] of [['domain', 'domains'], ['foundation', 'foundations']]) {
    const base = join(ROOT, dir);
    if (!existsSync(base)) continue;
    for (const entry of readdirSync(base, { withFileTypes: true })) {
      const index = join(base, entry.name, 'js', 'index.js');
      if (entry.isDirectory() && existsSync(index)) out.set(entry.name, { name: entry.name, kind, index });
    }
  }
  return out;
}

const UNITS = discover();
const loaded = new Map();

function load(name) {
  const unit = UNITS.get(name);
  if (!unit) throw Object.assign(new Error(`no unit ${name}: the units are ${[...UNITS.keys()].join(', ')}`), { kind: 'not-found' });
  if (!loaded.has(name)) loaded.set(name, import(pathToFileURL(unit.index).href));
  return loaded.get(name);
}

const isClass = (f) => /^class[\s{]/.test(Function.prototype.toString.call(f));

// The parameters of a function as written in its source: ['jd'], ['snap', 'key', 'count = 256'].
function parameters(f) {
  const s = Function.prototype.toString.call(f);
  const open = s.indexOf('('), arrow = s.indexOf('=>');
  if (arrow >= 0 && (open < 0 || arrow < open)) return [s.slice(0, arrow).replace(/^async\s+/, '').trim()];
  if (open < 0) return [];
  let depth = 0, close = open, quote = null;
  for (; close < s.length; close++) {
    const ch = s[close];
    if (quote) { if (ch === quote && s[close - 1] !== '\\') quote = null; continue; }
    if (ch === '"' || ch === "'" || ch === '`') quote = ch;
    else if ('([{'.includes(ch)) depth++;
    else if (')]}'.includes(ch) && --depth === 0) break;
  }
  const parts = [];
  let part = '';
  depth = 0;
  for (const ch of s.slice(open + 1, close)) {
    if ('([{'.includes(ch)) depth++;
    else if (')]}'.includes(ch)) depth--;
    if (ch === ',' && depth === 0) { parts.push(part.trim()); part = ''; } else part += ch;
  }
  if (part.trim()) parts.push(part.trim());
  return parts.map((p) => p.replace(/\s+/g, ' '));
}

// The comments written above each export of a unit's modules, and the modules'
// own descriptions (the comments after `export * from` in index.js).
function documentation(unit) {
  const docs = new Map(), modules = [];
  const text = readFileSync(unit.index, 'utf8');
  for (const m of text.matchAll(/^export \* from '(\.\/[^']+)';\s*(?:\/\/\s*(.*))?$/gm)) {
    const file = resolve(dirname(unit.index), m[1]);
    modules.push({ module: m[1].slice(2), description: (m[2] || '').trim() });
    if (!existsSync(file)) continue;
    const lines = readFileSync(file, 'utf8').split(/\r?\n/);
    lines.forEach((line, i) => {
      const e = line.match(/^export\s+(?:async\s+)?(?:function\*?|const|let|class)\s+([A-Za-z_$][\w$]*)/);
      if (!e) return;
      const comment = [];
      for (let k = i - 1; k >= 0 && /^\s*\/\//.test(lines[k]); k--) comment.unshift(lines[k].replace(/^\s*\/\/\s?/, ''));
      const inline = line.match(/;\s*\/\/\s*(.+)$/);
      docs.set(e[1], { module: m[1].slice(2), doc: (comment.join(' ') || (inline ? inline[1] : '')).trim() });
    });
  }
  return { docs, modules };
}

// A JSON-able constant, small enough to list.
function constant(v) {
  if (typeof v === 'number') return Number.isFinite(v) ? v : null;
  if (typeof v === 'string' || typeof v === 'boolean') return v;
  try {
    const s = JSON.stringify(plain(v));
    return s && s.length <= 2000 ? JSON.parse(s) : undefined;
  } catch { return undefined; }
}

async function describe(name) {
  const unit = UNITS.get(name);
  const { docs, modules } = documentation(unit);
  try {
    const ns = await load(name);
    const functions = [], constants = {}, classes = [];
    for (const key of Object.keys(ns).sort()) {
      const v = ns[key];
      const d = docs.get(key) || {};
      if (typeof v === 'function') {
        if (isClass(v)) classes.push(key);
        else functions.push({ name: key, parameters: parameters(v), module: d.module || null, doc: d.doc || '' });
      } else {
        const c = constant(v);
        if (c !== undefined) constants[key] = c;
      }
    }
    return { unit: name, kind: unit.kind, modules, functions, constants, classes };
  } catch (error) {
    return { unit: name, kind: unit.kind, modules, functions: [], constants: {}, classes: [], error: `${error.name}: ${error.message}` };
  }
}

function revive(v) {
  if (Array.isArray(v)) return v.map(revive);
  if (v && typeof v === 'object') {
    if (typeof v.$typed === 'string') {
      const Type = TYPED[v.$typed];
      if (!Type || !Array.isArray(v.values)) throw Object.assign(new Error(`$typed is one of ${Object.keys(TYPED).join(', ')}, with values`), { kind: 'failed' });
      return Type.from(v.values, (x) => (x === null ? NaN : x));
    }
    return Object.fromEntries(Object.entries(v).map(([k, x]) => [k, revive(x)]));
  }
  return v;
}

// A value made JSON: see the header.
function plain(v, depth = 0, seen = new Set()) {
  if (v === null || v === undefined) return null;
  if (typeof v === 'number') return Number.isFinite(v) ? v : null;
  if (typeof v === 'string' || typeof v === 'boolean') return v;
  if (typeof v === 'bigint') return v.toString();
  if (typeof v === 'function' || typeof v === 'symbol') return undefined;
  if (v instanceof Date) return Number.isFinite(v.getTime()) ? v.toISOString() : null;
  if (ArrayBuffer.isView(v) && !(v instanceof DataView)) return Array.from(v, (x) => plain(x));
  if (v instanceof ArrayBuffer) return Array.from(new Uint8Array(v));
  if (depth > 64) throw new Error('the result is nested more than 64 deep');
  if (seen.has(v)) throw new Error('the result refers to itself');
  seen.add(v);
  let out;
  if (Array.isArray(v) || v instanceof Set) out = [...v].map((x) => plain(x, depth + 1, seen) ?? null);
  else if (v instanceof Map) out = Object.fromEntries([...v].map(([k, x]) => [String(k), plain(x, depth + 1, seen) ?? null]));
  else {
    out = {};
    for (const [k, x] of Object.entries(v)) {
      const p = plain(x, depth + 1, seen);
      if (p !== undefined) out[k] = p;
    }
  }
  seen.delete(v);
  return out;
}

async function call(unitName, fname, args) {
  const ns = await load(unitName);
  if (typeof fname !== 'string' || !Object.hasOwn(ns, fname)) {
    throw Object.assign(new Error(`${unitName} exports no function ${fname}`), { kind: 'not-found' });
  }
  const f = ns[fname];
  if (typeof f !== 'function' || isClass(f)) {
    throw Object.assign(new Error(`${unitName}.${fname} is not a function that can be called (a constant or a class)`), { kind: 'not-found' });
  }
  if (!Array.isArray(args)) throw Object.assign(new Error('args is a list'), { kind: 'failed' });
  return plain(await f(...revive(args)));
}

async function answer(request) {
  const { id = null, op } = request;
  try {
    if (op === 'units') return { id, ok: true, result: await Promise.all([...UNITS.keys()].sort().map(describe)) };
    if (op === 'unit') {
      if (!UNITS.has(request.unit)) throw Object.assign(new Error(`no unit ${request.unit}`), { kind: 'not-found' });
      return { id, ok: true, result: await describe(request.unit) };
    }
    if (op === 'call') return { id, ok: true, result: await call(request.unit, request.function, request.args ?? []) };
    if (op === 'ping') return { id, ok: true, result: 'pong' };
    throw Object.assign(new Error(`no operation ${op}`), { kind: 'failed' });
  } catch (error) {
    return { id, ok: false, kind: error.kind || 'failed', error: `${error.kind ? '' : `${error.name}: `}${error.message}` };
  }
}

// One request at a time, answered in order; the worker ends with its input.
const lines = createInterface({ input: process.stdin, crlfDelay: Infinity });
let chain = Promise.resolve();
lines.on('line', (line) => {
  if (!line.trim()) return;
  chain = chain.then(async () => {
    let request;
    try { request = JSON.parse(line); } catch { process.stdout.write(`${JSON.stringify({ id: null, ok: false, kind: 'failed', error: 'a request is one line of JSON' })}\n`); return; }
    let reply;
    try { reply = JSON.stringify(await answer(request)); } catch (error) { reply = JSON.stringify({ id: request.id ?? null, ok: false, kind: 'failed', error: `the result cannot be written as JSON: ${error.message}` }); }
    process.stdout.write(`${reply}\n`);
  });
});
lines.on('close', () => { chain.then(() => process.exit(0)); });
