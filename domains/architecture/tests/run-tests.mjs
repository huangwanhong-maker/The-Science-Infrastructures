// Checks of the architecture domain: its C++ (a building as a complex of rooms and its events; its heat as a Hodge
// Laplacian under the Sun), built with CMake and run on calx's verification.
// Run from the repository root:  node domains/architecture/tests/run-tests.mjs

import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';

const { assert, skip, done, lines } = suite();
const domain = resolve(dirname(fileURLToPath(import.meta.url)), '..');

cppChecks({ assert, skip, lines }, { project: domain, name: 'architecture' });

done();
