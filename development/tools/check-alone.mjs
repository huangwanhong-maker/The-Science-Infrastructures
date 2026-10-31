// Each domain, run as a library on its own (ADR 0026).
//
// For each domain d, this copies into an empty directory only d's down-set: the
// domains d may depend on, transitively, the foundations they declare, and the
// harness of the checks. It then runs d's checks there. A domain that reached
// beyond its down-set would fail to load. The data library's cache is passed
// in as an input (SCIENCE_CACHE), so that the checks on real data run too.
//
// Run from the repository root:  node development/tools/check-alone.mjs [domain ...]

import { spawnSync } from 'node:child_process';
import { cpSync, existsSync, mkdtempSync, readFileSync, readdirSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const FOUNDATION_DEPENDS = { calx: [], morphe: [], physis: ['calx', 'morphe'], tyche: ['calx', 'morphe', 'physis'] };
const module = (d) => JSON.parse(readFileSync(join(ROOT, 'domains', d, 'module.json'), 'utf8'));
const all = readdirSync(join(ROOT, 'domains')).filter((d) => existsSync(join(ROOT, 'domains', d, 'module.json')));
const asked = process.argv.slice(2).length ? process.argv.slice(2) : all;

function downSet(d, seen = new Set()) {
  if (seen.has(d)) return seen;
  seen.add(d);
  for (const e of module(d).depends || []) downSet(e, seen);
  return seen;
}

let failures = 0;
for (const d of asked) {
  if (!existsSync(join(ROOT, 'domains', d, 'tests', 'run-tests.mjs'))) continue;
  const below = [...downSet(d)];
  const foundations = new Set();
  const addFoundation = (f) => { if (!foundations.has(f)) { foundations.add(f); (FOUNDATION_DEPENDS[f] || []).forEach(addFoundation); } };
  below.forEach((e) => (module(e).foundations || []).forEach(addFoundation));
  const dir = mkdtempSync(join(tmpdir(), `alone-${d}-`));
  try {
    // Each domain and foundation whole, its C++ with it, but not its history.
    const whole = (from, to) => cpSync(from, to, { recursive: true, filter: (src) => !/[\\/]\.git([\\/]|$)/.test(src) });
    for (const e of below) whole(join(ROOT, 'domains', e), join(dir, 'domains', e));
    for (const f of foundations) whole(join(ROOT, 'foundations', f), join(dir, 'foundations', f));
    cpSync(join(ROOT, 'development', 'tools', 'checks.mjs'), join(dir, 'development', 'tools', 'checks.mjs'));
    const r = spawnSync(process.execPath, [join(dir, 'domains', d, 'tests', 'run-tests.mjs')], {
      cwd: dir, encoding: 'utf8', maxBuffer: 1 << 26,
      // The copy's C++ is built in a directory removed afterwards (checks.mjs, cppChecks).
      env: { ...process.env, SCIENCE_CACHE: join(ROOT, 'data', 'cache'), SCIENCE_BUILD_EPHEMERAL: '1' },
    });
    const m = (r.stdout + r.stderr).match(/(\d+) of (\d+) checks passed\./);
    const ok = r.status === 0 && m && m[1] === m[2];
    if (!ok) failures++;
    const beneath = below.filter((e) => e !== d).sort().join(', ') || 'no other domain';
    console.log(`${ok ? 'ok  ' : 'FAIL'}  ${d.padEnd(13)} ${m ? `${m[1]} of ${m[2]}` : 'did not load'}, alone with ${beneath}${foundations.size ? `; ${[...foundations].sort().join(', ')}` : ''}`);
    if (!ok) console.log((r.stdout + r.stderr).replace(/^/gm, '        '));
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
}
process.exit(failures ? 1 : 0);
