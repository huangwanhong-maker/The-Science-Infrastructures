// The C++ checks of one library as a suite: build its CMake project (here, or in WSL),
// run its check programs, and report as every suite does.
// Run from the repository root:  node development/tools/run-cpp.mjs <library directory> <name>

import { cppChecks, suite } from './checks.mjs';

const [dir, name] = process.argv.slice(2);
if (!dir || !name) {
  console.error('usage: node development/tools/run-cpp.mjs <library directory> <name>');
  process.exit(2);
}
const s = suite();
cppChecks(s, { project: dir, name });
s.done();
