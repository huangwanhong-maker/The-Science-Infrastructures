// The dependency rules of the infrastructure, checked (ADR 0026).
//
// The domains form a partial order under "may depend on", declared in each
// domain's module.json; a domain's library is its down-set, and it must build
// and pass its tests with nothing else. This script checks that
//
//   1. the declared order is acyclic, and names only domains that exist;
//   2. every import in the JavaScript of the domains, the foundations, the shared
//      visualization layer and the applications, code and tests alike, reaches
//      another unit only through that unit's public interface (js/index.js), and
//      only a unit its module.json declares;
//   3. the code of the domains and foundations is pure: no network, no page,
//      no storage;
//   4. every domain has its README.md and module.json.
//
// Tests may also import the harness, development/tools/checks.mjs; a
// foundation's tests may not, since a foundation is a repository of its own.
//
// Run from the repository root:  node development/tools/check-dependencies.mjs

import { existsSync, readFileSync, readdirSync, statSync } from 'node:fs';
import { dirname, join, relative, resolve, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const rel = (p) => relative(ROOT, p).split(sep).join('/');
const problems = [];
const problem = (where, what) => problems.push(`${where}: ${what}`);

// Foundations may use one another as their own repositories declare.
const FOUNDATION_DEPENDS = { calx: [], morphe: [], physis: ['calx', 'morphe'], tyche: ['calx', 'morphe', 'physis'] };

// ---- the units and their declarations ----------------------------------------------
const units = new Map();       // key "domain:geodesy" -> { kind, name, dir, depends, foundations }
const read = (p) => JSON.parse(readFileSync(p, 'utf8'));
for (const name of readdirSync(join(ROOT, 'domains'))) {
  const dir = join(ROOT, 'domains', name);
  if (!statSync(dir).isDirectory()) continue;
  if (!existsSync(join(dir, 'README.md'))) problem(rel(dir), 'no README.md');
  if (!existsSync(join(dir, 'module.json'))) { problem(rel(dir), 'no module.json'); continue; }
  const m = read(join(dir, 'module.json'));
  if (m.name !== name) problem(rel(join(dir, 'module.json')), `name ${m.name} is not the directory's, ${name}`);
  if (existsSync(join(dir, 'js')) && !existsSync(join(dir, 'js', 'index.js'))) problem(rel(dir), 'js/ without its public interface, js/index.js');
  units.set(`domain:${name}`, { kind: 'domain', name, dir, depends: m.depends || [], foundations: m.foundations || [] });
}
{
  const m = read(join(ROOT, 'visualization', 'module.json'));
  units.set('edge:visualization', { kind: 'edge', name: 'visualization', dir: join(ROOT, 'visualization'), depends: m.depends || [], foundations: m.foundations || [] });
}
for (const name of readdirSync(join(ROOT, 'foundations'))) {
  const dir = join(ROOT, 'foundations', name);
  if (statSync(dir).isDirectory()) units.set(`foundation:${name}`, { kind: 'foundation', name, dir, depends: [], foundations: FOUNDATION_DEPENDS[name] || [] });
}
for (const name of readdirSync(join(ROOT, 'applications'))) {
  const dir = join(ROOT, 'applications', name);
  if (statSync(dir).isDirectory()) units.set(`application:${name}`, { kind: 'application', name, dir });
}
const domains = [...units.values()].filter((u) => u.kind === 'domain');
const domainNames = new Set(domains.map((d) => d.name));
for (const u of [...domains, units.get('edge:visualization')]) {
  for (const d of u.depends) if (!domainNames.has(d)) problem(`${u.name}/module.json`, `depends on ${d}, which is no domain`);
  for (const f of u.foundations) if (!units.has(`foundation:${f}`)) problem(`${u.name}/module.json`, `uses ${f}, which is no foundation`);
}

// ---- 1. the order is acyclic: its levels, by the longest chain below each domain ------------
const level = new Map();
const visiting = new Set();
function height(name, path = []) {
  if (level.has(name)) return level.get(name);
  if (visiting.has(name)) { problem('domains', `a cycle: ${[...path, name].join(' -> ')}`); return 0; }
  visiting.add(name);
  const d = units.get(`domain:${name}`);
  const h = d ? Math.max(-1, ...d.depends.filter((x) => domainNames.has(x)).map((x) => height(x, [...path, name]))) + 1 : 0;
  visiting.delete(name);
  level.set(name, h);
  return h;
}
domains.forEach((d) => height(d.name));

// ---- 2. imports ------------------------------------------------------------------------
function* files(dir) {
  if (!existsSync(dir)) return;
  for (const entry of readdirSync(dir, { withFileTypes: true })) {
    const p = join(dir, entry.name);
    if (entry.isDirectory()) { if (entry.name !== 'node_modules') yield* files(p); } else if (/\.(m?js)$/.test(entry.name)) yield p;
  }
}
function unitOf(path) {
  const r = rel(path).split('/');
  if (r[0] === 'domains' && r.length > 1) return units.get(`domain:${r[1]}`) || null;
  if (r[0] === 'visualization') return units.get('edge:visualization');
  if (r[0] === 'foundations' && r.length > 1) return units.get(`foundation:${r[1]}`) || null;
  if (r[0] === 'applications' && r.length > 1) return units.get(`application:${r[1]}`) || null;
  return null;
}
const interfaceOf = (u) => (u.kind === 'application' ? null : join(u.dir, 'js', 'index.js'));
const HARNESS = join(ROOT, 'development', 'tools', 'checks.mjs');
const SPECIFIERS = [
  /\b(?:import|export)\s[^'";]*?\bfrom\s*['"]([^'"]+)['"]/g,     // import ... from, export ... from
  /\bimport\s*['"]([^'"]+)['"]/g,                                 // import 'x'
  /\bimport\(\s*['"]([^'"]+)['"]\s*\)/g,                          // import('x')
  /\bnew URL\(\s*['"]([^'"]+)['"]\s*,\s*import\.meta\.url\s*\)/g, // new URL('x', import.meta.url), as for workers
];
let fileCount = 0, importCount = 0;
const sources = [
  ...domains.flatMap((d) => [...files(join(d.dir, 'js')), ...files(join(d.dir, 'tests'))]),
  ...files(join(ROOT, 'visualization', 'js')), ...files(join(ROOT, 'visualization', 'tests')),
  ...[...units.values()].filter((u) => u.kind === 'foundation').flatMap((u) => [...files(join(u.dir, 'js'))]),
  ...[...units.values()].filter((u) => u.kind === 'application').flatMap((u) => [...files(join(u.dir, 'src')), ...files(join(u.dir, 'plugins')), ...files(join(u.dir, 'tests'))]),
];
for (const file of sources) {
  fileCount++;
  const source = unitOf(file);
  const isTest = rel(file).split('/').includes('tests');
  const text = readFileSync(file, 'utf8').replace(/\/\*[\s\S]*?\*\//g, '').replace(/^\s*\/\/.*$/gm, '');
  for (const pattern of SPECIFIERS) {
    for (const [, spec] of text.matchAll(pattern)) {
      if (!spec.startsWith('.')) continue;          // packages ('three', 'node:fs') and URLs
      importCount++;
      const target = resolve(dirname(file), spec);
      const where = rel(file);
      if (target === HARNESS) {
        if (!isTest || source.kind === 'foundation') problem(where, 'only the tests of the infrastructure use the harness');
        continue;
      }
      const t = unitOf(target);
      if (!t) { problem(where, `imports ${rel(target)}, outside every unit`); continue; }
      if (t === source) continue;
      if (t.kind === 'application') { problem(where, `imports the application ${t.name}; nothing depends on an application`); continue; }
      if (target !== interfaceOf(t)) problem(where, `reaches past the public interface of ${t.name}: ${rel(target)}`);
      if (source.kind === 'application') continue;  // an application may use every domain, edge and foundation
      if (t.kind === 'domain' && !source.depends.includes(t.name)) problem(where, `${source.name} does not declare ${t.name} in its module.json`);
      if (t.kind === 'edge') problem(where, `${source.name} uses the visualization layer; only applications do`);
      if (t.kind === 'foundation' && !source.foundations.includes(t.name)) problem(where, `${source.name} does not declare the foundation ${t.name}`);
    }
  }
}

// ---- 3. purity of the domains and foundations ---------------------------------------------
const IMPURE = /\bfetch\s*\(|\bXMLHttpRequest\b|\bWebSocket\b|\bdocument\.|\bwindow\.|\blocalStorage\b|\bsessionStorage\b|\bfrom\s*['"]node:/;
for (const u of [...domains, ...[...units.values()].filter((x) => x.kind === 'foundation')]) {
  for (const file of files(join(u.dir, 'js'))) {
    if (rel(file).split('/').includes('tests')) continue;
    const text = readFileSync(file, 'utf8').replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/.*$/gm, '');
    const m = text.match(IMPURE);
    if (m) problem(rel(file), `is not pure: ${m[0].trim()}`);
  }
}

// ---- report ------------------------------------------------------------------------------
const byLevel = [];
for (const [name, h] of level) (byLevel[h] = byLevel[h] || []).push(name);
console.log('The order of the domains, by the longest chain beneath each:');
byLevel.forEach((names, h) => console.log(`  ${h}  ${names.sort().join(', ')}`));
console.log(`\n${fileCount} files, ${importCount} relative imports checked.`);
if (problems.length) {
  console.log(`\n${problems.length} problem(s):`);
  for (const p of problems) console.log(`  ${p}`);
  process.exit(1);
}
console.log('The dependency rules hold.');
