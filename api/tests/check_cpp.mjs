// The checks of the C and C++ clients (api/README.md §3.5), for run-tests.mjs:
// where a compiler is found, api/Makefile builds the libraries, the examples
// and the check programs (the C++ one also as C++20), and the programs are run
// against a service, each holding its answers to the Python client's for the
// same requests (check_cpp.py, whose JSON they read on standard input).
//
// The compiler is this computer's own (c++, g++ or clang++, with make), or, on
// Windows, WSL's g++. A program built in WSL runs in WSL, and must reach the
// service: in WSL's mirrored network (networkingMode=mirrored) it reaches the
// Windows loopback, and so the service run-tests.mjs started; otherwise (WSL's
// NAT, where the Windows loopback is not reached) a service is started inside
// WSL, with WSL's own Python and numpy, on a free port and with its sensors'
// store in a directory of its own, and stopped and removed afterwards; the
// Python client's values then come from WSL's Python too. Setting
// SCIENCE_CHECK_CPP_SERVICE=wsl takes that way even where the loopback is
// reached, so that it is checked too. Without a compiler, or without a way to
// a service, the checks are skipped, saying why.

import { spawn, spawnSync } from 'node:child_process';
import { existsSync } from 'node:fs';
import { join, resolve } from 'node:path';

const HERE = import.meta.dirname;
const ROOT = resolve(HERE, '../..');
const BUFFER = 1 << 26;
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

// A C++ and a C compiler with make: this computer's, or WSL's.
function toolchain() {
  const ok = (cmd, args) => spawnSync(cmd, args, { stdio: 'ignore' }).status === 0;
  for (const [cxx, cc] of [['c++', 'cc'], ['g++', 'gcc'], ['clang++', 'clang']]) {
    if (ok(cxx, ['--version']) && ok(cc, ['--version']) && ok('make', ['--version'])) {
      return { side: 'native', cxx, cc, path: (p) => p, run: (cmd, args, options = {}) => spawnSync(cmd, args, { encoding: 'utf8', maxBuffer: BUFFER, ...options }) };
    }
  }
  if (process.platform === 'win32' && ['g++', 'gcc', 'make'].every((t) => ok('wsl', ['--exec', t, '--version']))) {
    // --exec runs the program itself, not a shell, so that paths with spaces pass whole.
    const run = (cmd, args, options = {}) => spawnSync('wsl', ['--exec', cmd, ...args], { encoding: 'utf8', maxBuffer: BUFFER, ...options });
    return { side: 'wsl', cxx: 'g++', cc: 'gcc', run, path: (p) => run('wslpath', ['-a', p]).stdout.trim() };
  }
  return null;
}

// A service inside WSL, on a port WSL's programs reach at its own loopback.
async function wslService(tools, port, reachable) {
  const { run, path } = tools;
  if (run('python3', ['-c', 'import numpy']).status !== 0) return null;
  const dir = run('mktemp', ['-d', '-t', 'science-api-store.XXXXXX']).stdout.trim();
  const store = existsSync(resolve(ROOT, 'observation/sensors/store.py')) ? [] : ['--sensor-store', path(resolve(HERE, 'stub_store.py'))];
  const child = spawn('wsl', ['--exec', 'env', `SCIENCE_OBSERVATIONS=${dir}/observations`, `SCIENCE_PLUGIN_OUTPUT=${dir}/plugins`,
    'sh', '-c', 'echo $$; exec "$@"', 'sh', 'python3', path(resolve(ROOT, 'api/service/server.py')),
    '--port', String(port), '--product-seconds', '0', '--heartbeat-seconds', '0.5', ...store], { stdio: ['ignore', 'pipe', 'pipe'] });
  let log = '';
  child.stdout.on('data', (d) => { log += d; });
  child.stderr.on('data', (d) => { log += d; });
  const url = `http://127.0.0.1:${port}`;
  const service = {
    url, log: () => log,
    stop() {
      const pid = /^\d+/.exec(log)?.[0];
      if (pid) {
        run('kill', ['-INT', pid]);                       // the service closes its compute workers
        for (let i = 0; i < 20 && run('kill', ['-0', pid]).status === 0; i++) spawnSync('wsl', ['--exec', 'sleep', '0.1']);
        run('kill', ['-KILL', pid]);
      }
      child.kill();
      run('rm', ['-rf', dir]);
    },
  };
  for (let i = 0; i < 100; i++) {
    if (child.exitCode !== null) break;
    if (/^\d+/.test(log) && reachable(url)) return service;
    await sleep(150);
  }
  service.stop();
  throw new Error(`the service in WSL did not start: ${log.slice(-600)}`);
}

// The JSON lines of a check program, reported through the harness; how many there were.
function forward(stdout, { check, assert, skip }) {
  let seen = 0;
  for (const line of (stdout || '').split('\n')) {
    if (!line.startsWith('{')) continue;
    let rec;
    try { rec = JSON.parse(line); } catch { continue; }
    seen++;
    if (rec.kind === 'check') check(rec.name, rec.value ?? NaN, rec.limit, rec.unit || '');
    else if (rec.kind === 'assert') assert(rec.name, rec.ok, rec.detail || '');
    else if (rec.kind === 'skip') skip(rec.name);
  }
  return seen;
}

const tail = (r) => `${r.error ? `${r.error.message} ` : ''}${(r.stderr || '').trim().split('\n').slice(-6).join(' | ')}`;

export async function native({ service, python, freePort, check, assert, skip, work }) {
  const tools = toolchain();
  if (!tools) {
    skip('api c and c++: no C and C++ compilers with make (c++, g++ or clang++ here, or g++ in WSL), so the C and C++ clients are not checked');
    return;
  }
  const { run, path, side } = tools;
  const build = side === 'wsl' ? run('mktemp', ['-d', '-t', 'science-api-cpp.XXXXXX']).stdout.trim() : join(work, 'native');
  const bin = (name) => `${build}/${name}`;
  let own = null;
  try {
    const made = run('make', ['-C', path(resolve(ROOT, 'api')), '-j4', `BUILD=${build}`, `CXX=${tools.cxx}`, `CC=${tools.cc}`, 'all', 'check'], { timeout: 900000 });
    const warnings = (made.stderr || '').split('\n').filter((l) => /warning:/.test(l));
    assert(`api c and c++: the libraries, the examples and the checks build with make (${tools.cxx} and ${tools.cc}${side === 'wsl' ? ' in WSL' : ''}; C++17, C++20 and C99)`,
      made.status === 0, tail(made));
    if (made.status !== 0) return;
    assert('api c and c++: ... without a warning (-Wall -Wextra -pedantic)', warnings.length === 0, warnings.slice(0, 3).join(' | '));
    const twenty = run(bin('check_cpp20'), ['--offline'], { timeout: 60000 });
    const lines = (twenty.stdout || '').split('\n').filter((l) => l.startsWith('{')).map((l) => JSON.parse(l));
    assert(`c++ client: built as C++20, its ${lines.length} offline checks hold`, twenty.status === 0 && lines.length > 20 && lines.every((r) => r.kind !== 'assert' || r.ok),
      lines.filter((r) => r.kind === 'assert' && !r.ok).map((r) => r.name).join(', '));

    // The service: the one run-tests.mjs started, where the programs reach it; else one inside WSL.
    let url = `http://127.0.0.1:${service.port}`;
    const reachable = (u) => run(bin('check_cpp'), ['--probe', u], { timeout: 20000 }).status === 0;
    let reference = (u) => spawnSync(python, [resolve(HERE, 'check_cpp.py'), u], {
      cwd: ROOT, encoding: 'utf8', maxBuffer: BUFFER, timeout: 300000, env: { ...process.env, PYTHONIOENCODING: 'utf-8' },
    });
    if (side === 'wsl' && (process.env.SCIENCE_CHECK_CPP_SERVICE === 'wsl' || !reachable(url))) {
      own = await wslService(tools, await freePort(), reachable);
      if (!own) {
        skip('api c and c++: the programs built in WSL do not reach the service on the Windows loopback, and WSL has no Python with numpy to run one of its own');
        return;
      }
      url = own.url;
      reference = (u) => run('python3', [path(resolve(HERE, 'check_cpp.py')), u], { timeout: 300000 });
    }
    const ref = reference(url);
    assert(`api c and c++: the Python client's values for the same requests (check_cpp.py${own ? ', in WSL' : ''})`, ref.status === 0 && ref.stdout.startsWith('{'), tail(ref));
    if (ref.status !== 0) return;
    for (const [program, source] of [['check_cpp', 'api-check-cpp'], ['check_c', 'api-check-c']]) {
      const r = run(bin(program), [url, source], { input: ref.stdout, timeout: 300000 });
      const seen = forward(r.stdout, { check, assert, skip });
      assert(`api c and c++: ${program} ran to the end${own ? ' (against a service in WSL)' : ''}`, r.status === 0 && seen > 0, tail(r));
    }
  } catch (error) {
    assert('api c and c++: the checks ran', false, `${error.stack || error}${own ? ` ${own.log().slice(-600)}` : ''}`);
  } finally {
    if (own) own.stop();
    if (side === 'wsl') run('rm', ['-rf', build]);
  }
}
