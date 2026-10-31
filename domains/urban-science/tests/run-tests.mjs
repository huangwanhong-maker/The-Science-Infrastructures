// Checks of the urban-science domain: its C++ (the street system, its terrain, hours, mixing and visits; the
// rhythm of a record on the calendars' lattice; the response to forcing; the affinity of uses, catalysis and
// spacing; the crowd as a kinematic wave, person-time, transit and circulation; the classification), built with
// CMake and run on calx's verification.
// Run from the repository root:  node domains/urban-science/tests/run-tests.mjs

import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';

const { assert, skip, done, lines } = suite();
const domain = resolve(dirname(fileURLToPath(import.meta.url)), '..');

cppChecks({ assert, skip, lines }, { project: domain, name: 'urban-science' });

done();
