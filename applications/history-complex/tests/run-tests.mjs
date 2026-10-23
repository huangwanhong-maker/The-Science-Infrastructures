// Checks of what is the History Complex's own: how it writes dates (EDTF, BC years, the decimal year
// back to a date) and names; the foam's layout (deterministic, on the ring, heights in order of time);
// the fixtures it shows without the engine, each held to the protocol of the README (§3): indices in
// range, each lens's word a closed walk, χ = V − E + F, the Bianchi identity on every shell, attested
// signs kept by the reconstruction, spins in [0, 1], the foam's cells the lenses and the triangles,
// the persistence intervals in the order of time; and, when the engine is built, its answers on
// every canon of the domain. The science is checked with the history domain and morphe.
// Run from the repository root:  node applications/history-complex/tests/run-tests.mjs

import { existsSync, readdirSync, readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { spawnSync } from 'node:child_process';

import { suite, ROOT } from '../../../development/tools/checks.mjs';
import { dateText, isoOfYear, yearOfJdn, jdnOf, julianJdnOf, historicalJdnOf, yOfDate, shortNames, loopWords } from '../src/format.js';
import { placeBeings, heights, RING, HEIGHT } from '../src/foamlayout.js';

const { assert, check, skip, done } = suite();
const APP = resolve(ROOT, 'applications/history-complex');
const FIX = join(APP, 'tests/fixtures');
const read = (p) => JSON.parse(readFileSync(p, 'utf8'));

// ---- writing --------------------------------------------------------------------------------------
{
  assert('dates: a day, a month, a year, and an interval, written for reading',
    dateText('1785-02-12').includes('1785') && dateText('1797-03/1797-05').includes('1797') && dateText('1795?').includes('1795'), dateText('1797-03/1797-05'));
  assert('dates: a year BC', /BC/.test(dateText('-0398', 'julian')), dateText('-0398', 'julian'));
  assert('isoOfYear: 1808.75 is a day of October 1808', isoOfYear(1808.75).startsWith('1808-10'), isoOfYear(1808.75));
  assert('isoOfYear: a year BC keeps four digits and its sign', isoOfYear(-398) === '-0398-01-01', isoOfYear(-398));
  check('the decimal year of JDN 2451545 is 2000.0', yearOfJdn(jdnOf(2000, 1, 1)) - 2000, 1e-9);
  // The historians' convention, as the engine's chronology reads a reader's dates.
  check('the reform: 4 October 1582 (Julian) is followed by 15 October 1582 (Gregorian)', historicalJdnOf(1582, 10, 15) - historicalJdnOf(1582, 10, 4) - 1, 0);
  check('the reform begins at JDN 2299161', historicalJdnOf(1582, 10, 15) - 2299161, 0);
  check('before the reform a reader\'s date is Julian: 1 January 337 BC', historicalJdnOf(-336, 1, 1) - julianJdnOf(-336, 1, 1), 0);
  assert('the end of 338 BC is the Julian end of the year, so 337 BC begins after it', historicalJdnOf(-337, 12, 31) + 1 === julianJdnOf(-336, 1, 1) && jdnOf(-337, 12, 31) > julianJdnOf(-336, 1, 1), `${historicalJdnOf(-337, 12, 31)} ${julianJdnOf(-336, 1, 1)} ${jdnOf(-337, 12, 31)}`);
  check('yOfDate: the end of 338 BC is the last day of the year −337', yOfDate('-0337-12-31') - (-336), 1 / 300);
  assert('isoOfYear and yOfDate are inverse across the reform', isoOfYear(yOfDate('1582-10-20')) === '1582-10-20' && isoOfYear(yOfDate('1582-10-03')) === '1582-10-03', `${isoOfYear(yOfDate('1582-10-20'))} ${isoOfYear(yOfDate('1582-10-03'))}`);
  const names = shortNames(['Alexander von Humboldt', 'Wilhelm von Humboldt', 'Joseph Haydn']);
  assert('short names tell the two Humboldts apart', names[0] !== names[1] && names[2] === 'Haydn', names.join(', '));
  const kinds = shortNames(['Alexander the Great', 'Galileo Galilei', 'Académie royale des sciences', 'Society for Private Musical Performances', 'De revolutionibus orbium coelestium'],
    ['person', 'person', 'institution', 'group', 'work']);
  assert('short names: an epithet, a master by his first name, a Romance institution, a society, a work',
    kinds.join('|') === 'Alexander|Galileo|Académie|Musical Performances|De revolutionibus', kinds.join(', '));
}

// ---- the fixtures -----------------------------------------------------------------------------------
const builds = existsSync(FIX) ? readdirSync(FIX).filter((f) => /^build-.*\.json$/.test(f)) : [];
assert('fixtures: at least one canon built', builds.length > 0, builds.join(', '));
for (const file of builds) {
  const d = read(join(FIX, file));
  const name = file.replace(/^build-|\.json$/g, '');
  const nB = d.beings.length, nO = d.occasions.length, nT = d.ties.length;
  const inB = (i) => Number.isInteger(i) && i >= 0 && i < nB, inO = (i) => Number.isInteger(i) && i >= 0 && i < nO;
  assert(`${name}: occasions are the vertices, each once`, d.complex.vertices.length === nO && new Set(d.complex.vertices).size === nO);
  assert(`${name}: every edge is a being between two of its occasions`, d.complex.edges.every((e) => inB(e.being) && d.beings[e.being].occasions.includes(e.from) && d.beings[e.being].occasions.includes(e.to)));
  // Each lens's word is a closed walk on the edges.
  const vertexOf = new Map(d.complex.vertices.map((o, v) => [o, v]));
  const closed = d.complex.faces.every((f) => {
    let at = vertexOf.get(f.from);
    for (const x of f.boundary) {
      const e = d.complex.edges[Math.abs(x) - 1];
      const [tail, head] = [vertexOf.get(e.from), vertexOf.get(e.to)];
      if (x > 0 ? tail !== at : head !== at) return false;
      at = x > 0 ? head : tail;
    }
    return at === vertexOf.get(f.from);
  });
  assert(`${name}: every lens's word is a closed walk from its first meeting`, closed, `${d.complex.faces.length} lenses`);
  const [V, E, F] = d.invariants.counts;
  assert(`${name}: χ = V − E + F`, d.invariants.euler === V - E + F && V === nO && E === d.complex.edges.length && F === d.complex.faces.length, `${V} ${E} ${F} ${d.invariants.euler}`);
  const ranks = Object.fromEntries(d.invariants.homology.map((h) => [h.dimension, h.rank]));
  assert(`${name}: χ is the alternating sum of the Betti numbers`, d.invariants.euler === (ranks[0] || 0) - (ranks[1] || 0) + (ranks[2] || 0));
  assert(`${name}: the Bianchi identity on every shell`, d.invariants.shells.every((s) => s.product_to === s.product_from * (s.flips % 2 ? -1 : 1)), `${d.invariants.shells.length} shells`);
  assert(`${name}: attested rungs keep their signs`, d.ties.every((t) => t.rungs.every((g) => !g.attested || t.attested.some((a) => a.occasion === g.occasion && a.sign === g.sign))));
  assert(`${name}: ties name beings and occasions in range`, d.ties.every((t) => inB(t.a) && inB(t.b) && t.attested.every((a) => inO(a.occasion))));
  if (d.foam) {
    const f = d.foam;
    assert(`${name}: the foam's cells are the lenses and the triangles`, f.cells.length === F + f.triangles.length && f.cells.filter((c) => c.kind === 'lens').length === F);
    const j = f.cells.map((c) => c.spin);
    assert(`${name}: spins lie in [0, 1], within the sampling error`, Math.min(...j) > -0.06 && Math.max(...j) < 1.06, `${Math.min(...j).toFixed(3)} .. ${Math.max(...j).toFixed(3)}`);
    assert(`${name}: a cell outside the foams' support has spin 0`, f.cells.every((c) => c.support || c.spin === 0));
    const kept = f.ties.every((ft) => ft.posterior.every((p, k) => !d.ties[ft.tie].rungs[k].attested || p === d.ties[ft.tie].rungs[k].sign));
    assert(`${name}: the reconstruction keeps every attested sign`, kept);
    assert(`${name}: every reconstructed sign lies in [−1, 1]`, f.ties.every((ft) => ft.posterior.every((p) => p >= -1 && p <= 1)));
    // The band is the 5% and 95% points, which a mean need not lie between (a rare tail moves it).
    assert(`${name}: the posterior's frustration and its band lie between 0 and the triangles`,
      f.moments.every((m) => m.lo >= 0 && m.lo <= m.hi && m.hi <= m.triangles && m.mean >= 0 && m.mean <= m.triangles + 1e-9));
  }
  if (d.models) {
    const order = new Map(d.complex.vertices.map((o, k) => [o, k]));
    assert(`${name}: persistence: no class dies before it is born`, d.models.time.intervals.every((i) => i.death === null || order.get(i.death) >= order.get(i.birth)));
    const last = d.models.time.betti[d.models.time.betti.length - 1];
    assert(`${name}: persistence ends at the Betti numbers of the history`, last[1] === (ranks[1] || 0) && last[2] === (ranks[2] || 0), last.join(' '));
  }
}

// ---- the foam's layout ----------------------------------------------------------------------------
if (builds.length) {
  const d = read(join(FIX, builds[0]));
  const a = placeBeings(d), b = placeBeings(d);
  assert('foam: the layout of a canon is always the same', a.every((p, i) => p.x === b[i].x && p.z === b[i].z));
  check('foam: the farthest being stands on the ring', Math.max(...a.map((p) => Math.hypot(p.x, p.z))) - RING, 1e-9);
  const ys = heights(d, 'events'), yc = heights(d, 'calendar');
  const order = d.complex.vertices;
  assert('foam: heights by events rise with the order of time, from 0 to the top', order.every((o, k) => k === 0 || ys[o] > ys[order[k - 1]]) && ys[order[0]] === 0 && Math.abs(ys[order[order.length - 1]] - HEIGHT) < 1e-9);
  assert('foam: heights by calendar never fall', order.every((o, k) => k === 0 || yc[o] >= yc[order[k - 1]]));
  const w = loopWords(d.invariants.loops[0] || { beings: [], occasions: [], steps: [] }, d);
  assert('loops: told in words', typeof w.path === 'string');
}

// ---- the engine, when it is built ------------------------------------------------------------------
{
  const status = spawnSync('python', [join(APP, 'engine.py'), 'status'], { cwd: ROOT, encoding: 'utf8' });
  let state = null;
  try { state = JSON.parse(status.stdout); } catch { /* not built */ }
  // The engine reads the canon by the path its own side sees: under WSL, /mnt/c/....
  const sidePath = (p) => (state?.toolchain === 'wsl' ? p.replace(/^([A-Za-z]):[\\/]/, (_, c) => `/mnt/${c.toLowerCase()}/`).replace(/\\/g, '/') : p);
  if (!state?.current) skip('the engine is not built from the present sources: its answers on the canons are not checked');
  else {
    const canons = readdirSync(join(ROOT, 'domains/history/canons')).filter((f) => f.endsWith('.jsonl'));
    for (const c of canons) {
      const r = spawnSync('python', [join(APP, 'engine.py'), 'run', '-'], { cwd: ROOT, encoding: 'utf8', input: JSON.stringify({ op: 'build', canon: sidePath(join(ROOT, 'domains/history/canons', c)) }), maxBuffer: 1 << 28 });
      let answer = null;
      try { answer = JSON.parse(r.stdout); } catch { /* reported below */ }
      assert(`engine: ${c} built, its complex and foam complete`, answer && !answer.error && answer.foam && answer.models && answer.invariants.counts[0] === answer.occasions.length,
        answer?.error || (answer ? `${answer.invariants.counts.join(' ')}, ${answer.foam.foams} foams` : (r.stderr || '').split('\n').slice(-2).join(' ')));
    }
  }
}

done();
