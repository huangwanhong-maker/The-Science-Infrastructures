// Checks of the climate domain: the Köppen–Geiger classification on constructed
// climates (Beck et al. 2018, Table 2), and its symmetries.
// Run from the repository root:  node domains/climate/tests/run-tests.mjs
// The comparison with the published maps, which reads the data library's
// products, is the Earth Atlas's.

import { suite, seeded } from '../../../development/tools/checks.mjs';
import { koppen, koppenQuantities, KOPPEN_CODES, dryMonths } from '../js/index.js';

const { assert, done } = suite();

// Climates constructed to meet one rule each (Beck et al. 2018, Table 2):
// temperature a + b cos(2 pi (m - 6.5) / 12), warmest in July; precipitation by month.
const season = (a, b) => [...Array(12).keys()].map((m) => a + b * Math.cos((2 * Math.PI * (m - 6.5)) / 12));
const flat = (v) => Array(12).fill(v);
// Precipitation in summer (April-September) and in winter (October-March), month by month.
const halves = (summer, winter) => [...Array(12).keys()].map((m) => (m >= 3 && m <= 8 ? summer[m - 3] : winter[(m + 3) % 12]));
const cases = [
  ['Af', season(27, 1), flat(250), 'rainforest: every month at least 60 mm'],
  ['Am', season(27, 1), [40, ...Array(11).fill(300)], 'monsoon: a dry month, made up by the wet ones'],
  ['Aw', season(26, 2), halves(flat(200), flat(10)), 'savannah: a long dry winter'],
  ['BWh', season(28, 6), flat(5), 'hot desert'],
  ['BSk', season(10, 12), flat(25), 'cold steppe'],
  ['BSk', season(1, 6), flat(8), 'B before E: a dry climate with no month above 10 C is arid, not tundra'],
  ['Csa', season(17, 8), halves(flat(5), flat(100)), 'Mediterranean: dry, hot summer'],
  ['Cfb', season(11, 7), flat(60), 'oceanic'],
  ['Cwa', season(20, 7), halves(flat(250), flat(10)), 'dry winter, hot summer'],
  ['Dfb', season(5, 14), flat(50), 'humid continental'],
  ['Dwd', season(-12, 30), halves(flat(60), flat(2)), 'dry, very cold winter'],
  ['ET', season(-2, 8), flat(30), 'tundra'],
  ['EF', season(-30, 10), flat(10), 'ice cap'],
  ['Csa', season(15, 8), halves([3, 200, 3, 3, 3, 3], [15, 60, 60, 60, 60, 60]), 'both s and w met, more rain in winter: s'],
  ['Cwa', season(15, 8), halves([3, 400, 3, 3, 3, 3], [15, 60, 60, 60, 60, 60]), 'both s and w met, more rain in summer: w'],
];
const wrong = cases.filter(([code, T, P]) => koppen(T, P) !== code);
assert('köppen: fifteen constructed climates take their classes', !wrong.length,
  wrong.map(([code, T, P, what]) => `${what}: ${koppen(T, P)} for ${code}`).join('; ') || cases.map((c) => c[0]).join(' '));
assert('köppen: a missing month gives no class', koppen(season(10, 5), [NaN, ...flat(50).slice(1)]) === null);
const q = koppenQuantities(season(-10, 20), flat(30));
assert('köppen: summer is the warmer half-year', q.summer === 'April–September' && koppenQuantities(season(-10, -20), flat(30)).summer === 'October–March');

// The classes are invariant under the half-year shift of the months (Z/2), for
// climates whose half-years differ in temperature: summer and winter are
// defined by temperature, not by the calendar.
const shift = (a) => a.map((_, i) => a[(i + 6) % 12]);
let broken = 0, n = 0;
const rnd = seeded(7);
for (let k = 0; k < 4000; k++) {
  const T = season(-25 + 55 * rnd(), (rnd() - 0.5) * 50).map((t) => t + 2 * (rnd() - 0.5));
  const P = [...Array(12)].map(() => 300 * rnd() ** 2);
  if (Math.abs(T.slice(3, 9).reduce((a, b) => a + b) - shift(T).slice(3, 9).reduce((a, b) => a + b)) < 1e-9) continue;
  n++;
  if (koppen(T, P) !== koppen(shift(T), shift(P))) broken++;
}
assert('köppen: classes are invariant under the half-year shift', broken === 0, `${n - broken} of ${n} random climates`);
// More: within each half-year, the temperatures and the precipitations may be
// permuted independently, and the halves exchanged: the class is a function on
// the orbits of the group (S6 x S6) wr Z/2, of order 2 (6!)^4.
const HALF = [[3, 4, 5, 6, 7, 8], [9, 10, 11, 0, 1, 2]];
const shuffle = (r, xs) => { const a = [...xs]; for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(r() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; } return a; };
const act = (r, x, swap) => {
  const out = new Array(12);
  HALF.forEach((half, h) => { const to = HALF[swap ? 1 - h : h]; shuffle(r, half).forEach((m, k) => { out[to[k]] = x[m]; }); });
  return out;
};
let moved = 0, kept = 0;
const rnd3 = seeded(13);
for (let k = 0; k < 4000; k++) {
  const T = season(-25 + 55 * rnd3(), (rnd3() - 0.5) * 50).map((t) => t + 3 * (rnd3() - 0.5));
  const P = [...Array(12)].map(() => 300 * rnd3() ** 2);
  const swap = rnd3() < 0.5;
  const T2 = act(rnd3, T, swap), P2 = act(rnd3, P, swap);
  moved++;
  if (koppen(T, P) === koppen(T2, P2)) kept++;
}
assert('köppen: classes are invariant under the wreath product (S6 x S6) wr Z/2', kept === moved, `${kept} of ${moved} random climates and group elements`);
// Random climates with a wetter or drier half-year meet every class.
const classes = new Map();
const rnd2 = seeded(11);
for (let k = 0; k < 60000; k++) {
  const T = season(-40 + 75 * rnd2(), 35 * rnd2());
  const wet = [300 * rnd2() ** 2, 300 * rnd2() ** 2];
  const c = koppen(T, [...Array(12).keys()].map((m) => wet[m >= 3 && m <= 8 ? 0 : 1] * rnd2()));
  classes.set(c, (classes.get(c) || 0) + 1);
}
const rarest = [...classes].sort((x, y) => x[1] - y[1])[0];
assert('köppen: every one of the 30 classes is reachable', KOPPEN_CODES.every((c) => classes.has(c)), `${classes.size} classes met; the rarest, ${rarest[0]}, ${rarest[1]} times in 60 000`);

// Dry months after Gaussen (P <= 2 T).
const dry = dryMonths(season(17, 8), halves(flat(5), flat(100)));
assert('climate: the Mediterranean summer is dry, its winter wet', dry.slice(4, 9).every(Boolean) && !dry[0] && !dry[11]);

done();
