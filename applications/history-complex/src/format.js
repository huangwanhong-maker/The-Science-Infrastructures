// How the page writes dates, names, signs and numbers, and tells a loop of H₁ in words. Pure: the
// application's checks run these in Node.

const MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'];
const SHORT = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];

// ---- days -----------------------------------------------------------------------------------------

// The Julian day number of a Gregorian date (Meeus, Astronomical Algorithms, ch. 7), and back.
export function jdnOf(year, month, day) {
  const a = Math.floor((14 - month) / 12);
  const y = year + 4800 - a, m = month + 12 * a - 3;
  return day + Math.floor((153 * m + 2) / 5) + 365 * y + Math.floor(y / 4) - Math.floor(y / 100) + Math.floor(y / 400) - 32045;
}
export function gregorianOf(jdn) {
  const a = jdn + 32044, b = Math.floor((4 * a + 3) / 146097), c = a - Math.floor(146097 * b / 4);
  const d = Math.floor((4 * c + 3) / 1461), e = c - Math.floor(1461 * d / 4), m = Math.floor((5 * e + 2) / 153);
  return { year: 100 * b + d - 4800 + Math.floor(m / 10), month: m + 3 - 12 * Math.floor(m / 10), day: e - Math.floor((153 * m + 2) / 5) + 1 };
}
// The same for the Julian calendar (Meeus, ch. 7, without the Gregorian correction).
export function julianJdnOf(year, month, day) {
  const a = Math.floor((14 - month) / 12);
  const y = year + 4800 - a, m = month + 12 * a - 3;
  return day + Math.floor((153 * m + 2) / 5) + 365 * y + Math.floor(y / 4) - 32083;
}
export function julianOf(jdn) {
  const c = jdn + 32082, d = Math.floor((4 * c + 3) / 1461), e = c - Math.floor(1461 * d / 4), m = Math.floor((5 * e + 2) / 153);
  return { year: d - 4800 + Math.floor(m / 10), month: m + 3 - 12 * Math.floor(m / 10), day: e - Math.floor((153 * m + 2) / 5) + 1 };
}
// The historians' convention for a reader's dates, as the engine reads them (chronology.hpp): the
// Julian calendar before the Gregorian reform of 15 October 1582 (JDN 2299161), the Gregorian after.
const REFORM = 2299161;
export const historicalJdnOf = (year, month, day) => { const g = jdnOf(year, month, day); return g >= REFORM ? g : julianJdnOf(year, month, day); };
export const historicalOf = (jdn) => (jdn >= REFORM ? gregorianOf(jdn) : julianOf(jdn));
// The decimal year of a day in that convention, for drawing on a line of time.
export function yearOfJdn(jdn) {
  const { year } = historicalOf(jdn);
  const start = historicalJdnOf(year, 1, 1), end = historicalJdnOf(year + 1, 1, 1);
  return year + (jdn - start) / (end - start);
}
// The decimal year of a date the page holds (a year, a month or a day): its last day.
export function yOfDate(date) {
  const m = String(date).match(/^(-?\d{1,4})(?:-(\d\d))?(?:-(\d\d))?$/);
  if (!m) return null;
  const y = Number(m[1]), mo = Number(m[2] || 12);
  const d = m[3] ? Number(m[3]) : historicalJdnOf(mo === 12 ? y + 1 : y, mo === 12 ? 1 : mo + 1, 1) - historicalJdnOf(y, mo, 1);
  return yearOfJdn(historicalJdnOf(y, mo, d));
}
// The day a decimal year falls on, as an ISO date in that convention.
export function isoOfYear(y) {
  const year = Math.floor(y), start = historicalJdnOf(year, 1, 1), end = historicalJdnOf(year + 1, 1, 1);
  const g = historicalOf(Math.min(end - 1, start + Math.floor((y - year) * (end - start) + 1e-6)));
  const yy = `${g.year < 0 ? '-' : ''}${String(Math.abs(g.year)).padStart(4, '0')}`;
  return `${yy}-${String(g.month).padStart(2, '0')}-${String(g.day).padStart(2, '0')}`;
}

// ---- EDTF dates, written for reading ---------------------------------------------------------------

function part(p) {
  const q = p.match(/[~?%]$/)?.[0] || '';
  const m = p.replace(/[~?%]$/, '').match(/^(-?\d+)(?:-(\d\d))?(?:-(\d\d))?$/) || [];
  return { y: Number(m[1]), m: Number(m[2] || 0), d: Number(m[3] || 0), q };
}
// An astronomical year as it is read: year 0 is 1 BC, year −398 is 399 BC.
export const yearName = (y) => (y <= 0 ? `${1 - y} BC` : String(y));
// The round years of an axis from t0 to t1, every step: years BC where their names are round (440
// BC is the astronomical year −439), years AD where their numbers are.
export function roundYears(t0, t1, step) {
  const out = [];
  for (let y = Math.ceil(t0 / step) * step; y <= t1 + step; y += step) {
    const at = y <= 0 ? y + 1 : y;
    if (at >= t0 - 1e-9 && at <= t1 + 1e-9 && at !== out[out.length - 1]) out.push(at);
  }
  return out;
}
function onePart({ y, m, d }, short = false) {
  const months = short ? SHORT : MONTHS;
  if (d) return `${d} ${months[m - 1]} ${yearName(y)}`;
  if (m) return `${months[m - 1]} ${yearName(y)}`;
  return yearName(y);
}
const QUALIFIER = { '~': 'about ', '?': '', '%': 'about ' };
const DOUBT = { '?': ' (uncertain)', '%': ' (uncertain)', '~': '' };

// "1785-02-12" → "12 February 1785"; "1790-12~" → "about December 1790"; "1795?" → "1795 (uncertain)";
// "1797-03/1797-05" → "March–May 1797"; "1808-09-27/1808-10-14" → "27 Sep – 14 Oct 1808"; with the
// calendar when it is not the Gregorian.
export function dateText(date, calendar = 'gregorian') {
  if (!date) return '—';
  const pieces = String(date).split('/');
  let text;
  if (pieces.length === 2 && pieces[0] && pieces[1]) {
    const a = part(pieces[0]), b = part(pieces[1]);
    const q = a.q || b.q;
    if (a.y === b.y && a.m && b.m && !a.d && !b.d) text = `${MONTHS[a.m - 1]}–${MONTHS[b.m - 1]} ${yearName(a.y)}`;
    else if (a.y === b.y && a.m === b.m && a.d && b.d) text = `${a.d}–${b.d} ${MONTHS[a.m - 1]} ${yearName(a.y)}`;
    else if (a.y === b.y && a.d && b.d) text = `${a.d} ${SHORT[a.m - 1]} – ${b.d} ${SHORT[b.m - 1]} ${yearName(a.y)}`;
    else if (!a.m && !b.m) text = a.y <= 0 && b.y <= 0 ? `${1 - a.y}–${1 - b.y} BC` : `${yearName(a.y)}–${yearName(b.y)}`;
    else text = `${onePart(a, true)} – ${onePart(b, true)}`;
    text = (QUALIFIER[q] || '') + text + (DOUBT[q] || '');
  } else {
    const a = part(pieces[0] || pieces[1]);
    text = (QUALIFIER[a.q] || '') + onePart(a) + (DOUBT[a.q] || '');
  }
  // The Julian calendar is said where it could be mistaken for the Gregorian, after 1582.
  if (calendar === 'julian' && part(pieces[0] || pieces[1]).y >= 1582) text += ' (Julian)';
  else if (calendar === 'julian-march') text += ' (Old Style)';
  return text;
}
// The year an occasion is remembered by: its first.
export const yearText = (date) => { const m = String(date || '').match(/^-?\d{1,4}/); return m ? yearName(Number(m[0])) : '—'; };

// ---- names and words -------------------------------------------------------------------------------

const TITLES = new Set(['Duke', 'Earl', 'Viscount', 'Prince', 'Count', 'Marquess', 'Baron', 'Lord', 'Lady', 'King', 'Queen']);
// The name a being is called by in a sentence: "Joseph Haydn" → "Haydn", "Alexander I of Russia" →
// "Alexander I", "Duke of Wellington" → "Wellington", "G. W. F. Hegel" → "Hegel".
const NUMERAL = /^[IVXL]+$/;
const FIRST_NAMED = new Set(['Michelangelo', 'Leonardo', 'Raphael', 'Titian', 'Giotto', 'Galileo', 'Tycho']);
export function shortName(name) {
  if (!name) return '?';
  // What is in brackets or after a comma is a gloss ("Alma Mahler (born Schindler)", "Giuliano
  // de' Medici, Lorenzo's brother"), not the name.
  const s = String(name).replace(/\s*\([^)]*\)/g, '').split(', ')[0].trim();
  const tokens = s.split(/\s+/);
  const bare = TITLES.has(tokens[0]) || /^(Pope|Emperor|Empress|Saint|St\.)$/.test(tokens[0]) ? tokens.slice(1) : tokens;
  // A sovereign or a pope by the name and the numeral: "Julius II", "Cosimo I", "Franz Joseph I".
  const numeral = bare.findIndex((w, i) => i > 0 && NUMERAL.test(w));
  if (numeral > 0) return bare.slice(0, numeral + 1).join(' ');
  // A name with an epithet by the name: "Alexander the Great", "Lorenzo the Magnificent".
  const epithet = bare.findIndex((w, i) => i > 0 && w === 'the');
  if (epithet > 0) return bare.slice(0, epithet).join(' ');
  // Italian names by their custom: the Medici and Leonardo by the first name, Pico by "Pico",
  // Verrocchio by "Verrocchio"; the masters known by one name, by it.
  if (FIRST_NAMED.has(bare[0])) return bare[0];
  const particle = bare.findIndex((w, i) => i > 0 && /^(da|de'|de|di|del|della|dei|degli)$/.test(w));
  if (particle > 0) {
    if (/^(da|de')$/.test(bare[particle])) return bare[0];
    return particle >= 2 ? bare[particle - 1] : bare[bare.length - 1];
  }
  const of = s.split(' of ');
  if (of.length > 1) {
    const before = of[0].trim();
    if (TITLES.has(before)) return of.slice(1).join(' of ').split(/[,(]/)[0].trim();
    return before;
  }
  const words = s.replace(/^The\s+/, '').split(/\s+/);
  if (words.length > 1 && TITLES.has(words[0])) return words.slice(1).join(' ').split(' ').pop();
  // An institution or a group ("The Vienna Court Opera", "Society for Private Musical
  // Performances") keeps its last two words.
  if (/^The\s/.test(s) || words.includes('for')) return words.slice(-2).join(' ');
  return words[words.length - 1];
}

// An institution, a group, a state or a work by its proper name: without its article and what
// follows "of", kept whole up to three words ("Signoria", "Opera del Duomo", "Vienna Circle").
export function institutionName(name) {
  const s = String(name).replace(/\s*\([^)]*\)/g, '').split(', ')[0].trim().replace(/^The\s+/, '').split(' of ')[0].trim();
  const words = s.split(/\s+/);
  if (words.length <= 3) return s;
  // Named in the Romance way, the noun first and its qualifiers after in lower case ("Académie
  // royale des sciences"): by the noun.
  const lower = words.findIndex((w, i) => i > 0 && /^\p{Ll}/u.test(w));
  if (lower > 0 && !/^(for|and|of|in|on)$/.test(words[lower])) return words.slice(0, lower).join(' ');
  return words.slice(-2).join(' ');
}
// A work by the opening words of its title: "De revolutionibus orbium coelestium" is "De revolutionibus".
export function workName(name) {
  const words = String(name).replace(/\s*\([^)]*\)/g, '').trim().split(/\s+/);
  return words.length <= 3 ? words.join(' ') : words.slice(0, 2).join(' ');
}

// Short names that tell beings apart: where two share one ("Humboldt"), each takes the initial of
// its first name ("W. Humboldt", "A. Humboldt"); where that is not enough, its whole name. Types,
// when given, mark the institutions, groups, states and works.
export function shortNames(names, types = []) {
  const first = names.map((n, i) => (types[i] === 'work' ? workName(n) : types[i] && types[i] !== 'person' ? institutionName(n) : shortName(n)));
  const count = (xs) => xs.reduce((m, x) => m.set(x, (m.get(x) || 0) + 1), new Map());
  const c1 = count(first);
  // A person takes an initial; an institution, a group, a state or a work, its name without the gloss.
  const whole = (n) => String(n).replace(/\s*\([^)]*\)/g, '').split(', ')[0].trim().replace(/^The\s+/, '');
  const second = names.map((n, i) => (c1.get(first[i]) <= 1 ? first[i] : types[i] && types[i] !== 'person' ? whole(n) : `${String(n).trim()[0]}. ${first[i]}`));
  const c2 = count(second);
  return names.map((n, i) => (c2.get(second[i]) > 1 ? String(n) : second[i]));
}

// "A", "A and B", "A, B and C".
export function listText(items) {
  const a = items.filter((x) => x !== null && x !== undefined && x !== '');
  if (a.length <= 1) return a.join('');
  return `${a.slice(0, -1).join(', ')} and ${a[a.length - 1]}`;
}

export const signText = (s) => (s > 0 ? '+' : s < 0 ? '−' : '·');
export const signWord = (s) => (s > 0 ? 'positive' : s < 0 ? 'negative' : 'unknown');

// Numbers for reading: a true minus sign, a dash for nothing.
export function number(x, digits = 3) {
  if (x === null || x === undefined || Number.isNaN(x)) return '—';
  if (Number.isInteger(x)) return String(x).replace('-', '−');
  return Number(x).toFixed(digits).replace(/\.?0+$/, '').replace('-', '−');
}
export const plural = (n, one, many = `${one}s`) => `${n} ${n === 1 ? one : many}`;

// ---- a loop of H₁ in words ------------------------------------------------------------------------

// A loop passes beings and turns from one to the next at occasions: told as the beings in order,
// back to the first ("a → b → c → a"), and where it turns: the occasions, with their years.
export function loopWords(loop, data) {
  const names = loop.beings.map((b) => data.beings[b]?.short ?? shortName(data.beings[b]?.name));
  const path = names.length ? `${names.join(' → ')} → ${names[0]}` : '';
  const turns = loop.occasions.map((o) => {
    const occ = data.occasions[o];
    return occ ? `${occ.label} (${yearText(occ.date)})` : `occasion ${o}`;
  });
  return { path, turns, sentence: `${path}${turns.length ? `; where: ${listText(turns)}` : ''}.` };
}

// ---- the record's people ---------------------------------------------------------------------------

export function roleText(role, derived) {
  if (derived) return `${role || 'present'} (derived)`;
  return role || 'present';
}
