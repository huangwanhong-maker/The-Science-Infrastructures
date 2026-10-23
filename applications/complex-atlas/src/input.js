// A complex written by the user: Complex JSON, surface words ("abAB cdCD") or facets
// ("0 1 2; 0 2 3"), read into Complex JSON and checked before it is drawn or sent.

import { parseWords, parseFacets, checkComplex } from '../../../visualization/js/index.js';

// What kind of text it is, when the user has not said: JSON begins with a brace; facets are
// numbers with separators; anything else is read as words.
export function kindOf(text) {
  const t = String(text).trim();
  if (t.startsWith('{')) return 'json';
  if (/^[\d\s,;]+$/.test(t)) return 'facets';
  return 'words';
}

// The complex, or an Error whose message says what is wrong.
export function readComplex(text, kind = 'auto') {
  const k = kind === 'auto' ? kindOf(text) : kind;
  let complex;
  if (k === 'json') {
    try { complex = JSON.parse(text); } catch (error) { throw new Error(`not JSON: ${error.message}`); }
    if (complex && complex.complex && !complex.edges) complex = complex.complex;
  } else if (k === 'facets') complex = parseFacets(text);
  else complex = parseWords(text);
  const problems = checkComplex(complex);
  if (problems.length) throw new Error(problems.slice(0, 4).join('; ') + (problems.length > 4 ? `; and ${problems.length - 4} more` : ''));
  complex.name ??= k === 'json' ? 'yours' : k;
  complex.title ??= 'Your complex';
  return { complex, kind: k };
}
