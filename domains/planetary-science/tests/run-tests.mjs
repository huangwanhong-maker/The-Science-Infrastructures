// Checks of the planetary-science domain: its C++ (magnetism: the Moon's crustal field against
// pyshtools on T2015_449, the law of its spectrum and its unresolved part, the interplanetary field
// carried from L1), built with CMake and run on calx's verification. The model comes through the
// data library (python data/acquisition/fetch.py moon-magnetic-t2015); without it the checks on it
// are skipped, and say so.
// Run from the repository root:  node domains/planetary-science/tests/run-tests.mjs

import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';

const { assert, skip, done, lines } = suite();
const domain = resolve(dirname(fileURLToPath(import.meta.url)), '..');

cppChecks({ assert, skip, lines }, { project: domain, name: 'planetary-science' });

done();
