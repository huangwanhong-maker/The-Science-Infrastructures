// The harness of the checks: shared by the tests of the domains, the
// visualization layer and the applications, so that every suite reports in the
// same way and the runner (run-tests.mjs) can add them up.
//
// A check compares a value with a limit; an assertion states a condition. A
// suite whose real data are missing says what it skipped, and why.

import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

export const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
// The data library's cache is an input, not a dependency: SCIENCE_CACHE may
// point elsewhere, as when a domain's checks run alone (check-alone.mjs).
export const CACHE = process.env.SCIENCE_CACHE ? resolve(process.env.SCIENCE_CACHE) : resolve(ROOT, 'data/cache');
export const PRODUCTS = resolve(CACHE, 'products');

export function suite() {
  let failures = 0, count = 0;
  return {
    check(name, value, limit, unit = '') {
      count++;
      const ok = Number.isFinite(value) && Math.abs(value) <= limit;
      if (!ok) failures++;
      console.log(`${ok ? 'ok  ' : 'FAIL'}  ${name}: ${Number.isFinite(value) ? value.toExponential(3) : value} ${unit} (limit ${limit} ${unit})`);
    },
    assert(name, condition, detail = '') {
      count++;
      if (!condition) failures++;
      console.log(`${condition ? 'ok  ' : 'FAIL'}  ${name}${detail ? `: ${detail}` : ''}`);
    },
    skip(what) {
      console.log(`skip  ${what}`);
    },
    done() {
      console.log(`\n${count - failures} of ${count} checks passed.`);
      process.exit(failures ? 1 : 0);
    },
  };
}

// A reproducible sequence of uniform numbers in [0, 1) (a linear congruential
// generator, Numerical Recipes' constants).
export function seeded(seed) {
  let s = seed >>> 0;
  return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296);
}
