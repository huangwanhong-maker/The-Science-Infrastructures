// Checks of the anatomy domain: its C++ (the frames of the body, its parts and systems, BodyParts3D's body
// read from the data library's product, the measures of its parts and its sections), built with CMake and
// run on calx's verification; and its knowledge base (knowledge/structures.json), whose form is checked here
// and whose keys are held to the release's concepts when the product has been built. The checks that read the
// body skip when the product has not been built.
// Run from the repository root:  node domains/anatomy/tests/run-tests.mjs

import { existsSync, readFileSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';

const { assert, skip, done, lines } = suite();
const domain = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const root = resolve(domain, '../..');

// ---- the knowledge base --------------------------------------------------------------------------------------
const kb = JSON.parse(readFileSync(join(domain, 'knowledge/structures.json'), 'utf8'));
const entries = Object.entries(kb.entries || {});
const FIELDS = ['summary', 'function', 'supply', 'drainage', 'innervation', 'clinical'];
assert('the knowledge base has entries, a status and its references', entries.length >= 40 && ['seed', 'draft', 'stable'].includes(kb.status) && Object.keys(kb.references || {}).length >= 3);
assert('every key is an FMA identifier and every entry names its concept and has a summary',
  entries.every(([k, e]) => /^FMA\d+$/.test(k) && typeof e.name === 'string' && e.name.length > 1 && typeof e.summary === 'string'));
const sentences = entries.flatMap(([k, e]) => FIELDS.filter((f) => e[f] !== undefined).map((f) => [k, f, e[f]]));
const ragged = sentences.filter(([, , s]) => typeof s !== 'string' || !/^[A-Z0-9"‘'(]/.test(s) || !/[.!?)]$/.test(s.trim()));
assert('every field is written as sentences: a capital first, a full stop last', ragged.length === 0, ragged.slice(0, 3).map(([k, f]) => `${k}.${f}`).join(', '));
const uncited = entries.filter(([, e]) => !(e.refs || []).length || e.refs.some((r) => !kb.references[r]));
assert('every entry cites references the base lists', uncited.length === 0, uncited.map(([k]) => k).join(', '));
const unused = Object.keys(kb.references).filter((r) => !entries.some(([, e]) => (e.refs || []).includes(r)));
assert('every reference listed is cited', unused.length === 0, unused.join(', '));

const concepts = join(root, 'data/cache/products/body-bp3d-concepts.tsv');
if (!existsSync(concepts)) {
  skip('the knowledge base against the release: the product is not built (python data/products/body_anatomy.py)');
} else {
  const names = new Map(readFileSync(concepts, 'utf8').split(/\r?\n/).slice(1).filter(Boolean).map((l) => { const [c, n] = l.split('\t'); return [c, n]; }));
  const absent = entries.filter(([k]) => !names.has(k)).map(([k]) => k);
  assert('every key is a concept of BodyParts3D 4.0', absent.length === 0, absent.join(', '));
  const misnamed = entries.filter(([k, e]) => names.has(k) && names.get(k).toLowerCase() !== e.name.toLowerCase() && !(k === 'FMA7131' && e.name === 'oesophagus'));
  assert('every entry names its concept as FMA does (oesophagus for FMA\'s esophagus)', misnamed.length === 0, misnamed.map(([k, e]) => `${k}: ${e.name} ≠ ${names.get(k)}`).join('; '));
}

// ---- the C++ ------------------------------------------------------------------------------------------------------
cppChecks({ assert, skip, lines }, { project: domain, name: 'anatomy' });

done();
