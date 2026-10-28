// Checks of the navigation bridge: its C++ (the catalogue of the Solar System's frames, joining
// astronomy, metrology and geodesy), built with CMake and run on calx's verification, against the
// ephemeris, NAIF's worked example of the lunar frames and closed forms of rotating frames.
// Run from the repository root:  node domains/navigation/tests/run-tests.mjs
// The kernels and the Earth orientation come through the data library (python
// data/acquisition/fetch.py jpl-de440s naif-lunar-frames iers-eop); the checks that need them are
// skipped, and say so, when they are missing.

import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';

const { assert, skip, lines, done } = suite();
cppChecks({ assert, skip, lines }, { project: resolve(dirname(fileURLToPath(import.meta.url)), '..'), name: 'navigation' });
done();
