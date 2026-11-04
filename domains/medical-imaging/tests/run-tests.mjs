// Checks of the medical-imaging domain: its C++ (DICOM series and the patient frame, the modalities' physics
// and forward models, the reconstruction of vessels from one image, and their integration across modalities),
// built with CMake and run on calx's verification.
// Run from the repository root:  node domains/medical-imaging/tests/run-tests.mjs

import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';

const { assert, skip, done, lines } = suite();
const domain = resolve(dirname(fileURLToPath(import.meta.url)), '..');

cppChecks({ assert, skip, lines }, { project: domain, name: 'medical-imaging' });

done();
