// The harness of the checks: shared by the tests of the domains, the
// visualization layer and the applications, so that every suite reports in the
// same way and the runner (run-tests.mjs) can add them up.
//
// A check compares a value with a limit; an assertion states a condition. A
// suite whose real data are missing says what it skipped, and why.

import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { cpus, homedir } from 'node:os';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

export const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
// The data library's cache is an input, not a dependency: SCIENCE_CACHE may
// point elsewhere, as when a domain's checks run alone (check-alone.mjs).
export const CACHE = process.env.SCIENCE_CACHE ? resolve(process.env.SCIENCE_CACHE) : resolve(ROOT, 'data/cache');
export const PRODUCTS = resolve(CACHE, 'products');

export function suite() {
  let failures = 0, count = 0;
  return {
    // Lines printed by a check program in another language (C++ on calx/verification), counted
    // as this suite's own.
    lines(text) {
      for (const line of text.split('\n')) {
        if (/^(ok  |FAIL)  /.test(line)) { count++; if (line.startsWith('FAIL')) failures++; console.log(line); }
        else if (line.startsWith('skip  ')) console.log(line);
      }
    },
    check(name, value, limit, unit = '') {
      count++;
      const ok = Number.isFinite(value) && Math.abs(value) <= limit;
      if (!ok) failures++;
      console.log(`${ok ? 'ok  ' : 'FAIL'}  ${name}: ${Number.isFinite(value) ? value.toExponential(3) : value} ${unit} (limit ${limit} ${unit})`);
    },
    assert(name, condition, detail = '') {
      count++;
      if (!condition) failures++;
      console.log(`${condition ? 'ok  ' : 'FAIL'}  ${name}${detail ? `: ${detail}` : ''}`);
    },
    skip(what) {
      console.log(`skip  ${what}`);
    },
    done() {
      console.log(`\n${count - failures} of ${count} checks passed.`);
      process.exit(failures ? 1 : 0);
    },
  };
}

// A reproducible sequence of uniform numbers in [0, 1) (a linear congruential
// generator, Numerical Recipes' constants).
export function seeded(seed) {
  let s = seed >>> 0;
  return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296);
}

// ---- C++ checks ------------------------------------------------------------------------------------
//
// A library's C++ is a CMake project whose checks are programs on calx/verification,
// registered with CTest, printing what this harness prints. cppChecks() configures and
// builds the project with this computer's CMake and compiler, or on Windows with WSL's,
// in a build directory kept between runs (~/.cache/science-build/<name>-<hash of the path>;
// removed afterwards when SCIENCE_BUILD_EPHEMERAL is set, as check-alone sets it),
// runs CTest, and adds the programs' lines to the suite. Without a toolchain it says what it
// skipped. A program that ends without printing its count is a failure.

const BUFFER = 1 << 26;

export function cppToolchain() {
  const ok = (cmd, args) => { try { return spawnSync(cmd, args, { encoding: 'utf8' }).status === 0; } catch { return false; } };
  if (process.platform !== 'win32' && ok('cmake', ['--version']) && (ok('c++', ['--version']) || ok('g++', ['--version']) || ok('clang++', ['--version']))) {
    return { side: 'local', run: (cmd, args) => spawnSync(cmd, args, { encoding: 'utf8', maxBuffer: BUFFER }), path: (q) => q, home: homedir() };
  }
  if (process.platform === 'win32' && ok('wsl', ['--exec', 'cmake', '--version']) && ok('wsl', ['--exec', 'g++', '--version'])) {
    const run = (cmd, args) => spawnSync('wsl', ['--exec', cmd, ...args], { encoding: 'utf8', maxBuffer: BUFFER });
    const home = run('sh', ['-c', 'echo $HOME']).stdout.trim();
    return { side: 'wsl', run, path: (q) => run('wslpath', ['-a', q]).stdout.trim(), home };
  }
  return null;
}

export function cppChecks(s, { project, name }) {
  const tools = cppToolchain();
  if (!tools) { s.skip(`${name}: no CMake and C++ compiler (here, or in WSL), so the C++ checks are not built`); return false; }
  const source = tools.path(resolve(project));
  const hash = createHash('sha1').update(source).digest('hex').slice(0, 10);
  const ephemeral = Boolean(process.env.SCIENCE_BUILD_EPHEMERAL);
  const build = `${tools.home}/.cache/science-build/${ephemeral ? 'ephemeral-' : ''}${name.replace(/[^A-Za-z0-9_-]+/g, '-')}-${hash}`;
  try {
    return buildAndRun(s, tools, source, build, name);
  } finally {
    if (ephemeral) tools.run('rm', ['-rf', build]);
  }
}

function buildAndRun(s, tools, source, build, name) {
  const fail = (what, r) => {
    s.assert(`${name}: ${what}`, false, ((r.stdout || '') + (r.stderr || '')).split('\n').filter((l) => /error|Error|FAILED/.test(l)).slice(0, 12).join(' | ') || `exit ${r.status}`);
    return false;
  };
  let r = tools.run('cmake', ['-S', source, '-B', build, '-DCMAKE_BUILD_TYPE=Release']);
  if (r.status !== 0) return fail('CMake configures the project', r);
  r = tools.run('cmake', ['--build', build, '-j', String(Math.max(2, cpus().length))]);
  if (r.status !== 0) return fail(`the C++ builds (${tools.side === 'wsl' ? 'GCC in WSL' : 'this computer\'s compiler'})`, r);
  r = tools.run('ctest', ['--test-dir', build, '-V', '--timeout', '1800']);
  const out = (r.stdout || '').split('\n').map((l) => l.replace(/^\d+: /, ''));
  s.lines(out.join('\n'));
  // Programs that ended without their count crashed or were killed.
  const finished = out.filter((l) => /^\d+ of \d+ checks passed\.$/.test(l)).length;
  const total = Number(((r.stdout || '').match(/tests passed, \d+ tests failed out of (\d+)/) || [])[1] || 0);
  if (finished < total) s.assert(`${name}: every check program finishes`, false, `${total - finished} of ${total} ended without their count`);
  return r.status === 0;
}
