// Checks of the evolution domain: its C++ (every module's checks, built with CMake and run on
// calx's verification, with the genomics domain beneath it), and that every module the README
// names as written has its directory, its README and its CMake target.
// Run from the repository root:  node domains/evolution/tests/run-tests.mjs

import { existsSync, readFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';

const { assert, skip, done, lines } = suite();
const here = dirname(fileURLToPath(import.meta.url));
const domain = resolve(here, '..');

{
  const readme = readFileSync(resolve(domain, 'README.md'), 'utf8');
  const written = [...readme.matchAll(/^\| `([a-z_]+)` \|[^\n]*\| (draft|stable) \|$/gm)].map((m) => m[1]);
  assert('the README names its written modules', written.length > 0, written.join(', '));
  for (const m of written) {
    const dir = resolve(domain, m);
    const cmake = existsSync(resolve(dir, 'CMakeLists.txt')) ? readFileSync(resolve(dir, 'CMakeLists.txt'), 'utf8') : '';
    assert(`module ${m}: its directory, README and target evolution::${m}`, existsSync(resolve(dir, 'README.md')) && cmake.includes(`evolution::${m}`) && existsSync(resolve(dir, 'include', 'evolution', m)));
  }
}

cppChecks({ assert, skip, lines }, { project: domain, name: 'evolution' });

done();
