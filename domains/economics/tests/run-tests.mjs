// Checks of the economics domain: its C++ (transfer networks, circulation, production and the
// spacetime complex), built with CMake and run on calx's verification.
// Run from the repository root:  node domains/economics/tests/run-tests.mjs

import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';

const { assert, skip, done, lines } = suite();
const domain = resolve(dirname(fileURLToPath(import.meta.url)), '..');

cppChecks({ assert, skip, lines }, { project: domain, name: 'economics' });

done();
