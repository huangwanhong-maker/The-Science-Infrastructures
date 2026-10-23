// Every check of the infrastructure: the dependency rules first, then the suites
// of the foundations (their C++, built with CMake here or in WSL, then their
// JavaScript), of the domains in their order (each after those it stands on), of the
// shared visualization layer, and of the applications; last, each domain alone with
// its down-set.
// Run from the repository root:  node development/tools/run-tests.mjs [--verbose]

import { spawnSync } from 'node:child_process';
import { existsSync, readFileSync, readdirSync } from 'node:fs';
import { dirname, join, relative, resolve, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const verbose = process.argv.includes('--verbose');
const rel = (p) => relative(ROOT, p).split(sep).join('/');
const run = (file, ...args) => spawnSync(process.execPath, [file, ...args], { cwd: ROOT, encoding: 'utf8', maxBuffer: 1 << 26 });

const deps = run(join(ROOT, 'development/tools/check-dependencies.mjs'));
process.stdout.write(deps.stdout);
if (deps.status !== 0) { console.log('\nThe dependency rules are broken; the suites were not run.'); process.exit(1); }

// The domains in their order: by the longest chain of declared dependencies beneath each.
const declared = Object.fromEntries(readdirSync(join(ROOT, 'domains'), { withFileTypes: true })
  .filter((e) => e.isDirectory() && existsSync(join(ROOT, 'domains', e.name, 'module.json')))
  .map((e) => [e.name, JSON.parse(readFileSync(join(ROOT, 'domains', e.name, 'module.json'), 'utf8')).depends || []]));
const height = (n, memo = {}) => (memo[n] ??= Math.max(-1, ...declared[n].map((d) => height(d, memo))) + 1);
const ordered = Object.keys(declared).sort((a, b) => height(a) - height(b) || a.localeCompare(b));

// A suite is a script with its arguments, and the name it is reported under.
const script = (file) => ({ file, args: [], name: rel(file).replace(/\/(js\/)?tests\/run-tests\.mjs$/, '') });
const cpp = (dir) => ({ file: join(ROOT, 'development/tools/run-cpp.mjs'), args: [dir, rel(dir)], name: `${rel(dir)} (C++)` });
const foundations = readdirSync(join(ROOT, 'foundations'));
const suites = [
  ...foundations.map((f) => join(ROOT, 'foundations', f)).filter((d) => existsSync(join(d, 'CMakeLists.txt'))).map(cpp),
  ...[
  ...foundations.map((f) => join(ROOT, 'foundations', f, 'js/tests/run-tests.mjs')),
  ...ordered.map((d) => join(ROOT, 'domains', d, 'tests/run-tests.mjs')),
  join(ROOT, 'visualization/tests/run-tests.mjs'),
  join(ROOT, 'data/tests/run-tests.mjs'),
  join(ROOT, 'api/tests/run-tests.mjs'),
  ...readdirSync(join(ROOT, 'applications')).map((a) => join(ROOT, 'applications', a, 'tests/run-tests.mjs')),
  ].filter(existsSync).map(script),
];

console.log('\nSuites:');
let passed = 0, total = 0, skipped = 0, failedSuites = 0;
for (const { file, args, name } of suites) {
  const r = run(file, ...args);
  const out = r.stdout + r.stderr;
  const m = out.match(/(\d+) of (\d+) checks passed\./);
  const skips = (out.match(/^skip/gm) || []).length;
  if (m) { passed += Number(m[1]); total += Number(m[2]); }
  skipped += skips;
  const ok = r.status === 0 && m && m[1] === m[2];
  if (!ok) failedSuites++;
  console.log(`  ${ok ? 'ok  ' : 'FAIL'}  ${name.padEnd(28)} ${m ? `${m[1]} of ${m[2]}` : 'did not finish'}${skips ? `, ${skips} skipped` : ''}`);
  if (!ok || verbose) console.log(out.replace(/^/gm, '        '));
}
console.log(`\n${passed} of ${total} checks passed in ${suites.length} suites${skipped ? `; ${skipped} group(s) skipped for missing data` : ''}.`);

// Each domain again, alone with its down-set (check-alone.mjs).
console.log('\nEach domain alone, with only the domains and foundations beneath it:');
const alone = run(join(ROOT, 'development/tools/check-alone.mjs'));
process.stdout.write((alone.stdout + alone.stderr).replace(/^/gm, '  '));
process.exit(failedSuites || alone.status !== 0 ? 1 : 0);
