// Checks of the history domain: that the codebook's vocabularies (colours, types of beings,
// occasions and ties) are the ones the C++ reads, so that the document and the code cannot
// drift apart; that every canon is well-formed JSON lines naming this codebook; and the
// domain's C++ (the chronology, the record, the construction and the readings), built with
// CMake and run on calx's verification.
// Run from the repository root:  node domains/history/tests/run-tests.mjs

import { readFileSync, readdirSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';

const { assert, skip, done, lines } = suite();
const here = dirname(fileURLToPath(import.meta.url));
const domain = resolve(here, '..');
const read = (path) => readFileSync(resolve(domain, path), 'utf8');

// ---- the codebook and the code ------------------------------------------------------------------
{
  const codebook = read('codebook.md');
  const record = read('record/include/history/record/record.hpp');
  const section = (n) => codebook.split(/\n## /).find((s) => s.startsWith(`${n}. `)) || '';
  // A vocabulary in the C++: the strings of the vector a function returns.
  const vocabulary = (fn) => {
    const m = record.match(new RegExp(`${fn}\\(\\) \\{[\\s\\S]*?\\{([\\s\\S]*?)\\};`));
    return new Set(m ? [...m[1].matchAll(/"([^"]+)"/g)].map((x) => x[1]) : []);
  };
  const same = (a, b) => a.size === b.size && [...a].every((x) => b.has(x));
  const show = (a) => [...a].sort().join(', ');

  const [colourTable, signTable] = section(6).split('**The sign of a tie');
  const colours = new Set([...colourTable.matchAll(/^\| `([a-z-]+)` \|/gm)].map((m) => m[1]));
  assert('the codebook\'s colours are the C++\'s', colours.size === 6 && same(colours, vocabulary('colours')), `${show(colours)} | ${show(vocabulary('colours'))}`);
  const occasionTypes = new Set([...(section(5).match(/\*\*Types:\*\*(.*)/) || ['', ''])[1].matchAll(/`([a-z]+)`/g)].map((m) => m[1]));
  assert('the codebook\'s types of occasion are the C++\'s', occasionTypes.size > 10 && same(occasionTypes, vocabulary('occasion_types')), show(occasionTypes));
  const beingTypes = new Set([...section(3).split('\n\n')[1].matchAll(/\*\*([a-z]+)\*\*/g)].map((m) => m[1]));
  assert('the codebook\'s types of being are the C++\'s', same(beingTypes, vocabulary('being_types')), `${show(beingTypes)} | ${show(vocabulary('being_types'))}`);
  const tieTypes = new Set([...(signTable || '').matchAll(/^\| `([a-z]+)`(?:, `([a-z]+)`)? \(?/gm)].flatMap((m) => [m[1], m[2]]).filter(Boolean).filter((t) => !colours.has(t)));
  assert('the codebook\'s types of tie are the C++\'s', same(tieTypes, vocabulary('tie_types')), `${show(tieTypes)} | ${show(vocabulary('tie_types'))}`);
  const version = Number((codebook.match(/\*Version (\d+),/) || [])[1]);
  const coded = Number((record.match(/CODEBOOK = (\d+);/) || [])[1]);
  assert('the codebook\'s version is the one the C++ reads', version === coded && version > 0, `${version} and ${coded}`);
}

// ---- the canons -------------------------------------------------------------------------------------
{
  const files = readdirSync(resolve(domain, 'canons')).filter((f) => f.endsWith('.jsonl'));
  assert('the domain has a canon', files.length > 0, files.join(', '));
  for (const file of files) {
    const text = read(`canons/${file}`);
    const rows = text.split('\n').filter((l) => l.trim());
    let parsed = [];
    try {
      parsed = rows.map((l) => JSON.parse(l));
    } catch (e) {
      assert(`canon ${file}: every line is JSON`, false, e.message);
      continue;
    }
    assert(`canon ${file}: every line is JSON with a kind and an id`, parsed.every((x) => x.kind && x.id), `${parsed.length} lines`);
    assert(`canon ${file}: its first line is the canon, coded by codebook 1`, parsed[0].kind === 'canon' && parsed[0].codebook === 1);
    assert(`canon ${file}: it says who drafted it`, Boolean(parsed[0].drafted_by), parsed[0].drafted_by);
  }
}

// ---- the domain's C++ ----------------------------------------------------------------------------------
cppChecks({ assert, skip, lines }, { project: domain, name: 'history' });

done();
