// Checks of the robotics domain: its C++ (vehicles on a sphere at two fidelities and the translation
// between them, the waypoint follower, survey planning, sensors on a vehicle), built with CMake and run
// on calx's verification.
// Run from the repository root:  node domains/robotics/tests/run-tests.mjs

import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';

const { assert, skip, done, lines } = suite();
const domain = resolve(dirname(fileURLToPath(import.meta.url)), '..');

cppChecks({ assert, skip, lines }, { project: domain, name: 'robotics' });

done();
