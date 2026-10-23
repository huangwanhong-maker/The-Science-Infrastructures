// The History Complex: a canon of the history domain, built by the engine into its complex, shown
// as a loom, a moment, a complex, a record and its readings. The state lives in the address
// (state.js); the engine answers by the protocol of the README, §3 (api.js).

import * as api from './api.js';
import { el, clear, $ } from './dom.js';
import { readHash, writeHash } from './state.js';
import { Loom } from './loom.js';
import { FoamView } from './foam3d.js';
import { MomentView } from './moment.js';
import { renderComplex } from './complexview.js';
import { renderRecord } from './record.js';
import { renderReadings } from './readings.js';
import { occasionCard, tieCard, beingCard, canonCard } from './card.js';
import { COLOURS } from './palette.js';
import { yearOfJdn, dateText, signText, shortNames, yearText, number, isoOfYear } from './format.js';

const state = { ...readHash(), engine: null, canons: [], data: null, check: null, checkTimer: null };
let loom = null, moment = null, foam = null, player = null;

// ---- the actions the views may take ---------------------------------------------------------------
// A choice made in the foam stays in the foam; from any other view it is followed on the loom.
const followIn = () => (state.view === 'foam' ? 'foam' : 'loom');
const act = {
  pick(sel) { state.focus = sel; state.loop = null; go(followIn()); },
  moment(date) { state.date = normaliseDate(date); go('moment'); },
  showLoop(i, where) { state.loop = i; state.focus = null; go(where || followIn()); },
  tab(id) { state.tab = id; save(); render(); },
  runCheck,
  filterText: '',
};

function go(view) { state.view = view; save(true); render(); }
function save(push = false) {
  const h = writeHash(state);
  if (h !== location.hash) (push ? history.pushState : history.replaceState).call(history, null, '', h || '#');
}
function normaliseDate(date) {
  const m = String(date || '').match(/^-?\d{1,4}(-\d\d){0,2}/);
  return m ? m[0] : null;
}

// ---- starting ---------------------------------------------------------------------------------------
async function start() {
  for (const b of document.querySelectorAll('.views button')) b.addEventListener('click', () => { state.loop = state.view === 'loom' ? state.loop : null; go(b.dataset.view); });
  window.addEventListener('hashchange', () => { Object.assign(state, readHash()); render(); });
  window.addEventListener('popstate', () => { Object.assign(state, readHash()); render(); });
  window.addEventListener('keydown', (e) => { if (e.key === 'Escape' && (state.focus || state.loop !== null)) { state.focus = null; state.loop = null; save(); render(); } });
  $('canon').addEventListener('change', (e) => { state.canon = e.target.value; state.focus = null; state.loop = null; save(true); load(); });
  await connect();
  try { TOURS = await (await fetch('tours.json')).json(); } catch { TOURS = {}; }
  try {
    const { canons } = await api.canons();
    state.canons = canons || [];
  } catch (error) {
    return notice('The canons could not be listed.', error.message);
  }
  const select = clear($('canon'));
  for (const c of state.canons) select.append(el('option', { value: c.id, text: c.title || c.id }));
  if (!state.canons.length) return notice('There is no canon yet.', 'A canon is a file of the history domain, domains/history/canons/<id>.jsonl, coded by its codebook.');
  if (!state.canons.some((c) => c.id === state.canon)) state.canon = state.canons[0].id;
  select.value = state.canon;
  save();
  if (state.engine && !state.engine.current) return offerBuild();
  load();
}

async function connect() {
  try {
    state.engine = await api.engineState();
  } catch (error) {
    state.engine = { built: false, current: false, error: error.message };
  }
  const e = state.engine, badge = $('engine');
  badge.className = `engine ${e.current ? 'ok' : 'stale'}`;
  badge.textContent = e.fixtures ? 'fixtures' : e.current ? `engine · ${e.toolchain}` : e.built ? 'engine stale' : 'no engine';
  badge.title = e.fixtures ? 'Answers from tests/fixtures/, without the engine' : e.built_at ? `Built ${e.built_at}` : 'The engine is not built';
}

async function load() {
  const title = state.canons.find((c) => c.id === state.canon)?.title || state.canon;
  status('Building the complex…');
  notice(`Building the complex of ${title}…`, 'The engine reads the canon, constructs the complex and computes its homology and fundamental group exactly; it takes a few seconds.');
  $('notice').classList.add('quiet');
  try {
    state.data = enrich(await api.build(state.canon));
  } catch (error) {
    status('');
    $('notice').classList.remove('quiet');
    if (error.status === 503) return offerBuild();
    return notice('The canon could not be built.', error.message);
  }
  status('');
  $('notice').hidden = true;
  $('notice').classList.remove('quiet');
  state.check = null;
  loadCheck();
  render();
}

// Lookups and the last days of occasions, which the views share.
function enrich(d) {
  d.sourceById = new Map(d.sources.map((s) => [s.id, s]));
  shortNames(d.beings.map((b) => b.name), d.beings.map((b) => b.type)).forEach((s, i) => { d.beings[i].short = s; });
  for (const o of d.occasions) o.yearLast = yearOfJdn(o.last);
  d.view = state;
  return d;
}

// ---- the views --------------------------------------------------------------------------------------
function render() {
  for (const b of document.querySelectorAll('.views button')) b.setAttribute('aria-selected', String(b.dataset.view === state.view));
  for (const v of ['foam', 'loom', 'moment', 'complex', 'record', 'readings']) $(`view-${v}`).hidden = v !== state.view;
  const d = state.data;
  if (!d) return;
  d.view = state;
  if (state.view !== 'foam') stopPlaying();
  if (state.view === 'foam') renderFoam(d);
  else if (state.view === 'loom') renderLoom(d);
  else if (state.view === 'moment') {
    moment ??= new MomentView($('view-moment'), { load: (date) => api.moment(state.canon, date), pick: act.pick, setDate: (date) => { state.date = date; save(); render(); } });
    moment.render(d, state);
  } else if (state.view === 'complex') renderComplex($('view-complex'), d, act);
  else if (state.view === 'record') renderRecord($('view-record'), d, state, act, state.check);
  else if (state.view === 'readings') renderReadings($('view-readings'), d, act, state);
}

function renderLoom(d) {
  loom ??= new Loom($('loom-stage'), {
    onPick: act.pick,
    onHover: hover,
    onDate: (date) => { state.date = date; save(); render(); },
  });
  loom.set(d, state);
  // The disciplines, shown or hidden; the marks.
  const tools = clear($('loom-tools'));
  const present = new Set([...d.ties.map((t) => t.colour), ...d.beings.map((b) => b.colour)]);
  tools.append(el('div', { class: 'toggles' }, ...COLOURS.filter((c) => present.has(c.id)).map((c) => {
    const on = !state.hide.has(c.id);
    return el('button', { class: on ? 'toggle on' : 'toggle', 'aria-pressed': String(on), title: on ? `Hide ${c.name.toLowerCase()}` : `Show ${c.name.toLowerCase()}`,
      onclick: () => { on ? state.hide.add(c.id) : state.hide.delete(c.id); save(); render(); } },
    el('span', { class: 'swatch', style: `background:${c.hex}` }), c.name);
  })),
  el('div', { class: 'marks' },
    el('span', { class: 'mark-key' }, el('span', { class: 'key lens' }), 'lens, positive'),
    el('span', { class: 'mark-key' }, el('span', { class: 'key lens neg' }), 'negative'),
    el('span', { class: 'mark-key' }, el('span', { class: 'key flip' }), 'flips'),
    el('span', { class: 'mark-key' }, el('span', { class: 'key tri', text: '▲' }), 'frustrated triangle'),
    el('span', { class: 'mark-key' }, el('span', { class: 'key hollow' }), 'derived presence')),
  el('span', { class: 'hint', text: 'wheel: zoom · drag: pan · double-click: whole span · click the years: set the date' }));
  renderCard($('card'), d);
}

// The card of what is followed, beside the loom and the foam.
function renderCard(panel, d) {
  const card = clear(panel);
  const f = state.focus;
  // On a tour, its step comes first, and what it follows beneath.
  if (state.view === 'foam' && state.tour !== null && tourOf(d)) {
    card.append(tourCard(d));
    if (f?.kind === 'being' && d.beings.some((b) => b.id === f.id)) card.append(beingCard(d, d.beings.findIndex((b) => b.id === f.id), act, state.check));
    else if (f?.kind === 'tie' && d.ties.some((t) => t.id === f.id)) card.append(tieCard(d, d.ties.findIndex((t) => t.id === f.id), act));
    else if (f?.kind === 'occasion' && d.occasions.some((o) => o.id === f.id)) card.append(occasionCard(d, d.occasions.findIndex((o) => o.id === f.id), act));
    else if (state.loop !== null && d.invariants.loops[state.loop]) card.append(loopCard(d, state.loop));
    return;
  }
  if (state.loop !== null && d.invariants.loops[state.loop]) {
    card.append(loopCard(d, state.loop));
  } else if (f?.kind === 'being' && d.beings.some((b) => b.id === f.id)) card.append(beingCard(d, d.beings.findIndex((b) => b.id === f.id), act, state.check));
  else if (f?.kind === 'tie' && d.ties.some((t) => t.id === f.id)) card.append(tieCard(d, d.ties.findIndex((t) => t.id === f.id), act));
  else if (f?.kind === 'occasion' && d.occasions.some((o) => o.id === f.id)) card.append(occasionCard(d, d.occasions.findIndex((o) => o.id === f.id), act));
  else card.append(canonCard(d), state.view === 'foam' ? foamCard(d) : dateCard());
  if (f || state.loop !== null) card.prepend(el('button', { class: 'link close', onclick: () => { state.focus = null; state.loop = null; save(); render(); } }, '← the canon'));
}

// ---- the foam ---------------------------------------------------------------------------------------
const FOAM_MODES = [
  ['discipline', 'Disciplines', 'Each lens in its tie\'s colour, darkening where the tie is negative'],
  ['spin', 'Spins', 'Each lens by its spin in the foam of the prior: how far the tie\'s change belongs to a closed configuration of changes (theory §6.2)'],
  ['flips', 'Flips', 'The lenses across which a tie changes sign: their holonomy is −1'],
  ['inference', 'Inferred', 'The lenses with an unattested end, brighter where the inferred sign is less certain'],
];

function renderFoam(d) {
  foam ??= new FoamView($('foam-stage'), { onPick: act.pick, onHover: hover });
  window.__foam = foam;  // for the gallery's pictures and for exploring from the console
  const order = d.complex.vertices;
  const first = d.occasions[order[0]].year, last = d.occasions[order[order.length - 1]].year;
  foam.set(d, { mode: state.mode, opened: state.opened, time: state.time, date: state.date, dateText: state.date ? dateText(state.date) : '',
    focus: state.focus, loop: state.loop, hide: state.hide });
  const tools = clear($('foam-tools'));
  const choice = (on, label, title, onclick) => el('button', { class: on ? 'toggle on' : 'toggle', 'aria-pressed': String(on), title, onclick }, label);
  const year = state.date ? Number((state.date.match(/^-?\d{1,4}/) || [first])[0]) : null;
  const slider = el('input', { type: 'range', min: String(Math.floor(first)), max: String(Math.ceil(last)), step: '1', value: String(year ?? Math.floor(first)),
    'aria-label': 'The date of the cut' });
  slider.addEventListener('input', () => { stopPlaying(); state.date = isoOfYear(Number(slider.value)); save(); render(); });
  tools.append(...[
    el('div', { class: 'toggles' }, el('span', { class: 'tool-label', text: 'Lenses' }), ...FOAM_MODES.map(([id, label, title]) => choice(state.mode === id, label, title, () => { state.mode = id; save(); render(); }))),
    el('div', { class: 'toggles' }, el('span', { class: 'tool-label', text: 'Occasions' }),
      choice(!state.opened, 'Knots', 'The history complex K: each occasion a point where the threads meet', () => { state.opened = false; save(); render(); }),
      choice(state.opened, 'Opened', 'The refined complex K̂: each occasion opened into its presences, with the rungs of its ties and their triangles (theory §5)', () => { state.opened = true; save(); render(); })),
    el('div', { class: 'toggles' }, el('span', { class: 'tool-label', text: 'Time' }),
      choice(state.time === 'calendar', 'Calendar', 'Heights by the calendar', () => { state.time = 'calendar'; save(); render(); }),
      choice(state.time === 'events', 'Events', 'Heights by the order of the occasions, each one step', () => { state.time = 'events'; save(); render(); })),
    el('div', { class: 'toggles cut' }, el('span', { class: 'tool-label', text: 'Cut' }), slider,
      choice(Boolean(player), player ? 'Pause' : 'Play', 'Raise the cut through the history', () => (player ? stopPlaying() : play(first, last))),
      state.date ? choice(false, 'No cut', 'Remove the cut', () => { stopPlaying(); state.date = null; save(); render(); }) : null,
      choice(false, 'Turn', 'Turn the history slowly round', () => { foam.spin(!foam.controls.autoRotate); }),
      choice(false, 'Reset', 'The first view', () => foam.reset())),
    tourOf(d) ? el('div', { class: 'toggles' }, el('span', { class: 'tool-label', text: 'Guide' }),
      choice(state.tour !== null, state.tour !== null ? 'On tour' : 'A tour', 'A guided tour of the canon, step by step', () => (state.tour !== null ? (state.tour = null, save(), render()) : takeStep(0)))) : null,
    el('div', { class: 'toggles' }, ...colourToggles(d))].filter(Boolean));
  renderCard($('foam-card'), d);
  foamLegend();
}

// The key inside the stage: what each mark of the foam is, under the chosen colouring.
function foamLegend() {
  const stage = $('foam-stage');
  let key = stage.querySelector('.foam-key');
  if (!key) { key = el('div', { class: 'foam-key' }); stage.append(key); }
  const row = (mark, text) => el('div', { class: 'foam-key-row' }, el('span', { class: `fk ${mark}` }), text);
  const lens = {
    discipline: row('film', 'a lens, in its discipline; darker where the tie is negative'),
    spin: row('film gold', 'a lens, brighter the likelier it is in a spin foam'),
    flips: row('film gold', 'a lens across which the tie changes sign'),
    inference: row('film mauve', 'a lens with an inferred end, brighter the less certain'),
  }[state.mode];
  key.replaceChildren(...[
    row('thread', 'a life'), row('knot', 'an occasion'), lens,
    state.opened ? row('rung', 'a tie at a meeting: gold +, violet −, dashed if inferred') : null,
    row('tri', state.opened ? 'a frustrated triangle of ties' : 'an occasion with a frustrated triangle'),
    state.date ? row('cut', 'the network of the moment at the cut') : null].filter(Boolean));
}

function colourToggles(d) {
  const present = new Set([...d.ties.map((t) => t.colour), ...d.beings.map((b) => b.colour)]);
  return COLOURS.filter((c) => present.has(c.id)).map((c) => {
    const on = !state.hide.has(c.id);
    return el('button', { class: on ? 'toggle on' : 'toggle', 'aria-pressed': String(on), title: on ? `Hide ${c.name.toLowerCase()}` : `Show ${c.name.toLowerCase()}`,
      onclick: () => { on ? state.hide.add(c.id) : state.hide.delete(c.id); save(); render(); } },
    el('span', { class: 'swatch', style: `background:${c.hex}` }), c.name);
  });
}

// The cut rises through the history, a year at a time.
function play(first, last) {
  stopPlaying();
  let y = state.date ? Number((state.date.match(/^-?\d{1,4}/) || [first])[0]) : Math.floor(first);
  if (y >= last) y = Math.floor(first);
  const span = Math.max(1, last - first);
  const step = Math.max(1, Math.round(span / 120));
  player = setInterval(() => {
    y += step;
    if (y > last) { stopPlaying(); render(); return; }
    state.date = isoOfYear(y);
    history.replaceState(null, '', writeHash(state));
    render();
  }, Math.max(60, Math.min(400, 24000 / (span / step))));
}
function stopPlaying() { if (player) { clearInterval(player); player = null; } }

// ---- the tours ---------------------------------------------------------------------------------------
let TOURS = {};
const tourOf = (d) => (Array.isArray(TOURS[d.canon.id]) && TOURS[d.canon.id].length ? TOURS[d.canon.id] : null);

// Takes a step of the canon's tour: its date, what it follows, its colouring, its occasions opened
// or not, its scale of time, and a loop named by its beings.
function takeStep(k) {
  const d = state.data, steps = tourOf(d);
  if (!steps) return;
  state.tour = Math.max(0, Math.min(steps.length - 1, k));
  const s = steps[state.tour].state || {};
  stopPlaying();
  state.date = s.date ?? null;
  state.focus = s.focus ?? null;
  state.mode = s.mode || 'discipline';
  state.opened = Boolean(s.opened);
  state.time = s.time === 'events' ? 'events' : 'calendar';
  state.loop = null;
  if (Array.isArray(s.loop)) {
    const want = new Set(s.loop);
    const i = d.invariants.loops.findIndex((l) => s.loop.every((id) => l.beings.some((b) => d.beings[b].id === id)) && l.beings.every((b) => want.has(d.beings[b].id)));
    const j = i >= 0 ? i : d.invariants.loops.findIndex((l) => s.loop.every((id) => l.beings.some((b) => d.beings[b].id === id)));
    state.loop = j >= 0 ? j : null;
  }
  state.view = 'foam';
  save();
  render();
}

function tourCard(d) {
  const steps = tourOf(d), k = state.tour, step = steps[k];
  return el('div', { class: 'card tour' },
    el('p', { class: 'eyebrow', text: `A tour of the canon · ${k + 1} of ${steps.length}` }),
    el('h2', { text: step.title }),
    el('p', { class: 'tour-text', text: step.text }),
    el('div', { class: 'row tour-steps' },
      el('button', { class: 'button', disabled: k === 0, onclick: () => takeStep(k - 1) }, '← Before'),
      k + 1 < steps.length ? el('button', { class: 'button primary', onclick: () => takeStep(k + 1) }, 'Next →') : null,
      el('button', { class: 'link', onclick: () => { state.tour = null; save(); render(); } }, 'End the tour')),
    el('ol', { class: 'tour-list' }, ...steps.map((s, i) => el('li', { class: i === k ? 'on' : '' }, el('button', { class: 'link', onclick: () => takeStep(i) }, s.title)))));
}

function foamCard(d) {
  const f = d.foam;
  const coupling = (of) => f?.couplings.find((c) => c.of === of);
  const bal = coupling('balance'), all = coupling('all');
  return el('div', { class: 'card date' },
    el('p', { class: 'eyebrow', text: 'Reading the foam' }),
    el('p', { class: 'note', text: 'Time rises. Each being is a thread, each occasion a knot where the threads of those present meet, and each tie between two of their meetings a ribbon: a lens, a 2-cell of the complex. Opened, an occasion shows its presences, the rungs of its ties (gold positive, violet negative, dashed where inferred) and their triangles, vermilion where frustrated. A cut at a date is the network of that moment.' }),
    f ? el('p', { class: 'note', text: `The theory fitted to this record: persistence tanh β = ${number(all?.t, 2)} over all ties, balance tanh β_Δ = ${number(bal?.t, 2)}; ${f.foams} independent foams.` }) : null,
    el('p', { class: 'note' }, 'How a record becomes a complex, and its signs a spin foam: ',
      el('a', { href: '/domains/history/theory/the-history-complex.md', target: '_blank', rel: 'noopener' }, 'the theory note'), '; how a history is recorded: ',
      el('a', { href: '/domains/history/codebook.md', target: '_blank', rel: 'noopener' }, 'the codebook'), '.'),
    state.date ? el('button', { class: 'button', onclick: () => go('moment') }, `The moment at ${dateText(state.date)}`) : el('p', { class: 'note', text: 'Drag the cut to a date, or play it upward through the history.' }));
}

function dateCard() {
  if (!state.date) return el('p', { class: 'note', text: 'Choose a being, a tie or an occasion to follow it; click the scale of years to choose a date.' });
  return el('div', { class: 'card date' }, el('p', { class: 'eyebrow', text: 'The chosen date' }), el('h3', { text: dateText(state.date) }),
    el('button', { class: 'button', onclick: () => go('moment') }, 'The moment at this date'));
}

function loopCard(d, i) {
  const loop = d.invariants.loops[i];
  return el('div', { class: 'card' }, el('p', { class: 'eyebrow', text: `A loop of H₁ · ${i + 1} of ${d.invariants.loops.length}` }),
    el('h2', { text: loop.beings.map((b) => d.beings[b].short).join(' → ') }),
    el('p', { class: 'note', text: 'A circle of lives meeting at occasions, across which the record holds no tie: the loop is drawn in mauve along the lanes it follows. A flat connection may wind around it; a historian may look here for a relation the canon is silent about.' }),
    el('ul', { class: 'occasions' }, ...loop.steps.map((s) => el('li', {}, el('span', { class: 'num', text: d.beings[s.being].short }),
      `${d.occasions[s.from].label} → ${d.occasions[s.to].label}`))));
}

// The note that follows the pointer on the loom.
function hover(h) {
  const box = $('hover'), d = state.data;
  if (!h || !d) { box.hidden = true; return; }
  let text;
  if (h.kind === 'occasion') { const o = d.occasions[h.occasion]; text = `${dateText(o.date, o.calendar)} · ${o.label}`; }
  else if (h.kind === 'being') text = d.beings[h.being].name;
  else if (h.kind === 'lens') {
    const f = d.complex.faces[h.face], t = d.ties[f.tie];
    const cell = d.foam?.cells.find((c) => c.kind === 'lens' && c.index === h.face);
    text = `${d.beings[t.a].short} — ${d.beings[t.b].short} (${t.type}): ${signText(f.signs[0])} → ${signText(f.signs[1])}${f.flip ? ', a flip' : ''}${cell && cell.support ? ` · spin ${number(cell.spin, 2)}` : ''}`;
  } else if (h.kind === 'triangle') {
    const T = d.foam.triangles[h.triangle];
    const product = T.ties.reduce((s, ti) => s * (d.ties[ti].rungs.find((g) => g.occasion === T.occasion)?.sign ?? 1), 1);
    text = `${T.beings.map((b) => d.beings[b].short).join(', ')} at ${d.occasions[T.occasion].label}: ${product < 0 ? 'frustrated' : 'balanced'}`;
  }
  box.textContent = text;
  box.hidden = false;
}
document.addEventListener('pointermove', (e) => {
  const box = $('hover');
  if (box.hidden) return;
  box.style.left = `${Math.min(e.clientX + 14, window.innerWidth - box.offsetWidth - 8)}px`;
  box.style.top = `${e.clientY + 16}px`;
});

// ---- the check against Wikidata ---------------------------------------------------------------------
async function loadCheck() {
  try {
    state.check = await api.check(state.canon);
  } catch {
    state.check = null;
  }
  clearTimeout(state.checkTimer);
  if (state.check?.status === 'running') state.checkTimer = setTimeout(async () => { await loadCheck(); render(); }, 3000);
}
async function runCheck() {
  try {
    state.check = await api.runCheck(state.canon);
  } catch (error) {
    state.check = { status: 'failed', message: error.message };
  }
  render();
  await loadCheck();
  render();
}

// ---- the engine, and what goes wrong ----------------------------------------------------------------
function status(text) { $('status').textContent = text; }

function notice(lead, ...rest) {
  const n = clear($('notice'));
  n.append(el('p', { class: 'lead', text: lead }), ...rest.map((r) => (typeof r === 'string' ? el('p', { class: 'note', text: r }) : r)));
  n.hidden = false;
}

function offerBuild() {
  const log = el('pre', { class: 'log' });
  const button = el('button', { class: 'button primary', onclick: () => build(button, log) }, state.engine?.built ? 'Build it again' : 'Build the engine');
  notice(state.engine?.built ? 'The engine was built from older sources.' : 'The engine has not been built on this computer.',
    el('p', { class: 'note', text: 'It is the history domain\'s C++, with morphe\'s topology, compiled once with g++ (natively, or in WSL on Windows); it takes a minute or two. To see the page without it, open it with ?fixtures in the address.' }),
    el('div', { class: 'row' }, button), log);
}

async function build(button, log) {
  button.disabled = true;
  try {
    await api.buildEngine();
  } catch (error) {
    log.textContent = error.message;
    button.disabled = false;
    return;
  }
  const poll = async () => {
    const s = await api.setup();
    log.textContent = s.log || '';
    if (s.running) return setTimeout(poll, 1000);
    if (s.state === 'done') { $('notice').hidden = true; await connect(); load(); return; }
    button.disabled = false;
  };
  poll();
}

start();
