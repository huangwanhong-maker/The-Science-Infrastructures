// Checks of the data library: each data/tests/check_*.py, run with the Python of
// the repository's environment (.venv) or, failing that, the first Python with
// numpy, and read line by line into the shared harness (see checks.py).
// Run from the repository root:  node data/tests/run-tests.mjs

import { spawnSync } from 'node:child_process';
import { existsSync, readdirSync } from 'node:fs';
import { resolve } from 'node:path';

import { suite } from '../../development/tools/checks.mjs';

const { check, assert, skip, done } = suite();
const HERE = import.meta.dirname, ROOT = resolve(HERE, '../..');
const pythons = [resolve(ROOT, '.venv/Scripts/python.exe'), resolve(ROOT, '.venv/bin/python'), 'python3', 'python'];
const python = pythons.find((p) => (p.includes('.venv') ? existsSync(p) : true) && spawnSync(p, ['-c', 'import numpy'], { stdio: 'ignore' }).status === 0);

const scripts = readdirSync(HERE).filter((f) => /^check_.*\.py$/.test(f)).sort();
if (!python) {
  skip('data: no Python with numpy found');
} else {
  for (const script of scripts) {
    const r = spawnSync(python, [resolve(HERE, script)], { cwd: ROOT, encoding: 'utf8', timeout: 900000, maxBuffer: 1 << 26 });
    let seen = 0;
    for (const line of (r.stdout || '').split('\n')) {
      if (!line.startsWith('{')) continue;
      let rec;
      try { rec = JSON.parse(line); } catch { continue; }
      seen++;
      if (rec.kind === 'check') check(rec.name, rec.value, rec.limit, rec.unit || '');
      else if (rec.kind === 'assert') assert(rec.name, rec.ok, rec.detail || '');
      else if (rec.kind === 'skip') skip(rec.name);
    }
    assert(`${script}: ran to the end`, r.status === 0, (r.stderr || '').trim().split('\n').slice(-3).join(' | '));
    if (!seen) skip(`${script}: reported nothing`);
  }
}

done();
