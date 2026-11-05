// Checks of the physical-computing domain: its C++ (the formalism: tasks, schedules, the commuting square, the
// cost of a verified result and the materials' marks, the backend interface; and the fluids family's
// simulators), built with CMake and run on calx's verification. Run from the repository root:
//   node domains/physical-computing/tests/run-tests.mjs

import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';

const { assert, skip, done, lines } = suite();
const domain = resolve(dirname(fileURLToPath(import.meta.url)), '..');

cppChecks({ assert, skip, lines }, { project: domain, name: 'physical-computing' });

done();
