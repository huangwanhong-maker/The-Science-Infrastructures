// The layer panel: what is drawn on the surface and by which model (or, for
// climate classes, which period, scenario and source), what moves, which lines
// and events are shown, the laboratory's own model, and how the base map looks.

import { cssGradient } from '../render/colormaps.js';
import { DOMAINS, byKey } from '../layers/catalogue.js';

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

export const LANGUAGES = [
  ['en', 'English'], ['zh', '中文（简体）'], ['zht', '中文（繁體）'], ['ja', '日本語'], ['ko', '한국어'],
  ['fr', 'Français'], ['de', 'Deutsch'], ['es', 'Español'], ['pt', 'Português'], ['it', 'Italiano'],
  ['ru', 'Русский'], ['ar', 'العربية'], ['hi', 'हिन्दी'],
];

function choice(name, value, label, checked, swatch) {
  const mark = swatch ? `<span class="swatch" style="background:${swatch}"></span>` : '<span class="dotmark"></span>';
  return `<label class="choice"><input type="radio" name="${name}" value="${esc(value)}"${checked ? ' checked' : ''}>${mark}<span>${esc(label)}</span></label>`;
}

// The Köppen–Geiger hues in a gradient, for the swatch of a class map.
const CLASS_SWATCH = 'linear-gradient(90deg, #6f86c9, #8ea9dc 18%, #d98d7f 30%, #e2c07c 42%, #cfd08a 55%, #93c08f 66%, #a98ab8 78%, #79b4c4 88%, #b9b3b3)';

function swatchOf(def) {
  if (def.kind === 'classes') return CLASS_SWATCH;
  const st = def.style;
  return st ? cssGradient(st.colormap, { reverse: st.reverse, clip: st.clip }) : null;
}

// Buttons choosing one model among the definitions of a key.
function models(kind, defs, current) {
  if (defs.length < 2) return '';
  return `<div class="models" role="group" aria-label="Model">${defs.map((d) =>
    `<button type="button" data-model="${kind}" data-id="${esc(d.id)}" aria-pressed="${d.id === current}">${esc(d.model)}</button>`).join('')}</div>`;
}

// Maps of classes: a row of periods, the scenarios of a future period, the
// sources of a period, and another map to compare with. Each button chooses a
// map by its id, keeping the scenario where it can.
const shortPeriod = (p) => `${p[0]}–${String(p[1]).slice(String(p[0]).slice(0, 2) === String(p[1]).slice(0, 2) ? 2 : 0)}`;
const samePeriod = (a, b) => a[0] === b[0] && a[1] === b[1];

function classChoices(kind, defs, active, compare) {
  const html = [];
  const button = (d, label, pressed) => `<button type="button" data-model="${kind}" data-id="${esc(d.id)}" aria-pressed="${pressed}">${esc(label)}</button>`;
  const library = defs.filter((d) => !d.computed);
  const periods = [];
  for (const d of library) if (!periods.some((p) => samePeriod(p, d.period))) periods.push(d.period);
  periods.sort((a, b) => a[0] - b[0]);
  const pick = (period, scenario) => library.find((d) => samePeriod(d.period, period) && d.scenario === scenario)
    || library.find((d) => samePeriod(d.period, period) && d.scenario === 'SSP2-4.5')
    || library.find((d) => samePeriod(d.period, period));
  html.push(`<div class="models periods" role="group" aria-label="Period">${periods.map((p) => button(pick(p, active.scenario), shortPeriod(p), samePeriod(p, active.period))).join('')}</div>`);
  const scenarios = library.filter((d) => samePeriod(d.period, active.period) && d.scenario);
  if (scenarios.length) html.push(`<div class="models" role="group" aria-label="Scenario">${scenarios.map((d) => button(d, d.scenario, d.id === active.id)).join('')}</div>`);
  const sources = defs.filter((d) => samePeriod(d.period, active.period) && !d.scenario);
  if (sources.length > 1) html.push(`<div class="models" role="group" aria-label="Source">${sources.map((d) => button(d, d.computed ? 'Ours, from normals' : 'Beck et al.', d.id === active.id)).join('')}</div>`);
  const options = [`<option value="">—</option>`, ...defs.filter((d) => d.id !== active.id)
    .map((d) => `<option value="${esc(d.id)}"${compare && d.id === compare.id ? ' selected' : ''}>${esc(d.model)}</option>`)];
  const word = compare && compare.period && compare.period[1] < active.period[0] ? 'changed since' : 'differs from';
  html.push(`<label class="compare"><span>${word}</span><select data-compare aria-label="Compare with">${options.join('')}</select></label>`);
  return html.join('');
}

// A list of keyed choices: one radio per key; the models of the chosen key below it.
function keyed(kind, groups, selected, { swatches = false, compare = null } = {}) {
  const html = [];
  for (const [key, defs] of groups) {
    const active = defs.find((d) => d.id === selected);
    const def = active || defs[0];
    html.push(choice(kind, def.id, def.label, Boolean(active), swatches ? swatchOf(def) : null));
    if (active && active.kind === 'classes') {
      html.push(classChoices(kind, defs, active, compare && compare.id ? compare : null));
    } else if (active && defs.length > 1) {
      html.push(models(kind, defs, active.id));
      if (compare) {
        const options = [`<option value="">—</option>`, ...defs.filter((d) => d.id !== active.id)
          .map((d) => `<option value="${esc(d.id)}"${d.id === compare.id ? ' selected' : ''}>${esc(d.model)}</option>`)];
        html.push(`<label class="compare"><span>minus</span><select data-compare aria-label="Compare with">${options.join('')}</select></label>`);
      }
    }
  }
  return html.join('');
}

export function buildPanel(el, catalogue, state, years, extra, on) {
  const html = [];
  // ---- surface ----
  html.push('<h3 style="margin-top:6px;border-top:0;padding-top:0">Surface</h3>');
  const fields = catalogue.surfaces.filter((s) => s.kind === 'field' || s.kind === 'classes');
  const indicators = catalogue.surfaces.filter((s) => s.kind === 'indicator');
  const compare = catalogue.surfaces.find((s) => s.id === state.compare) || null;
  const byDomain = new Map();
  for (const s of fields) {
    if (!byDomain.has(s.domain)) byDomain.set(s.domain, []);
    byDomain.get(s.domain).push(s);
  }
  for (const [domain, list] of byDomain) {
    html.push(`<div class="group-name">${esc(DOMAINS[domain] || domain)}</div>`);
    html.push(keyed('surface', byKey(list), state.surface, { swatches: true, compare: compare || {} }));
  }
  if (indicators.length) {
    const active = indicators.find((s) => s.id === state.surface);
    const current = active || indicators.find((s) => s.id === state.indicator) || indicators[0];
    html.push(`<div class="group-name">${esc(DOMAINS.geography)}, by country</div>`);
    const options = indicators.map((s) => `<option value="${esc(s.id)}"${s.id === current.id ? ' selected' : ''}>${esc(s.label)}</option>`).join('');
    html.push(`<label class="choice indicator"><input type="radio" name="surface" value="${esc(current.id)}"${active ? ' checked' : ''}><span class="swatch" style="background:${swatchOf(current)}"></span><select id="indicator" aria-label="Indicator">${options}</select></label>`);
    if (active && years) {
      html.push(`<div class="indicator-row"><span></span><span class="year"><button type="button" data-year="-1" aria-label="Earlier year">‹</button><span>${years.year}</span><button type="button" data-year="1" aria-label="Later year">›</button></span></div>`);
    }
  }
  html.push(choice('surface', 'none', 'Relief only', state.surface === 'none', null));

  // ---- motion, lines, events ----
  html.push('<h3>Motion</h3>');
  html.push(choice('flow', 'none', 'None', state.flow === 'none'));
  html.push(keyed('flow', byKey(catalogue.flows), state.flow));
  if (catalogue.lines.length) {
    html.push('<h3>Lines</h3>');
    html.push(choice('lines', 'none', 'None', state.lines === 'none'));
    html.push(keyed('lines', byKey(catalogue.lines), state.lines));
  }
  if (extra.plates) {
    html.push('<h3>Solid Earth</h3>');
    html.push(`<label class="check"><input type="checkbox" data-plates${state.plates ? ' checked' : ''}> Plate boundaries and names</label>`);
    if (extra.mechanisms) {
      html.push(`<label class="check"><input type="checkbox" data-mechanisms${state.mechanisms ? ' checked' : ''}> Earthquake mechanisms (Global CMT)</label>`);
      if (state.mechanisms) {
        html.push(`<div class="segmented sub-choice">${[['year', 'The past year'], ['large', 'Magnitude 6.5 and more since 1976']]
          .map(([k, l]) => `<button type="button" data-mechanism-mode="${k}" aria-pressed="${state.mechanismMode === k}">${l}</button>`).join('')}</div>`);
      }
    }
    if (extra.stations) {
      html.push(`<label class="check"><input type="checkbox" data-stations${state.stations ? ' checked' : ''}> GPS stations: how the ground moves</label>`);
      if (state.stations) {
        html.push(`<div class="segmented sub-choice">${[['plate', 'Relative to their plates'], ['igs', 'IGS20']]
          .map(([k, l]) => `<button type="button" data-station-frame="${k}" aria-pressed="${state.stationFrame === k}">${l}</button>`).join('')}</div>`);
      }
    }
  }
  if (catalogue.tracks.length) {
    html.push('<h3>Cyclones</h3>');
    for (const t of catalogue.tracks) {
      const on = state.tracks.includes(t.id);
      const swatch = t.kind === 'official' ? 'linear-gradient(90deg, #f6ece8 50%, transparent 50%) 0 0 / 8px 100%' : extra.modelColour(t.model);
      html.push(`<label class="check track"><input type="checkbox" data-track="${esc(t.id)}"${on ? ' checked' : ''}><span class="track-key" style="background:${swatch}"></span> ${esc(t.label)}</label>`);
    }
  }
  if (catalogue.events.length) {
    html.push('<h3>Events</h3>');
    const kinds = new Map();
    for (const e of catalogue.events) {
      if (!kinds.has(e.kind)) kinds.set(e.kind, []);
      kinds.get(e.kind).push(e);
    }
    for (const [kind, list] of kinds) {
      const checked = list.some((e) => e.id === state.events);
      const label = kind === 'earthquake' ? 'Earthquakes' : kind === 'disaster' ? 'Disaster alerts (GDACS)' : kind;
      html.push(`<label class="check"><input type="checkbox" data-events="${esc(list[0].id)}"${checked ? ' checked' : ''}> ${esc(label)}</label>`);
      if (list.length > 1) {
        const span = (e) => {
          const days = (Date.parse(e.window.end) - Date.parse(e.window.start)) / 86400000;
          return days > 300 ? 'past year' : days > 20 ? 'past month' : days > 5 ? 'past week' : 'past day';
        };
        html.push(`<div class="segmented" style="margin-left:26px">${list.map((e) => `<button type="button" data-events-window="${esc(e.id)}" aria-pressed="${e.id === state.events}">${span(e)}</button>`).join('')}</div>`);
      }
    }
  }

  // ---- the laboratory's model ----
  if (extra.lab) {
    const lab = extra.lab;
    const kinds = [['barotropic', 'Barotropic'], ['shallow-water', 'Shallow water']];
    html.push('<h3>Laboratory models</h3>');
    html.push(`<div class="lab"><p>The flow at 500 hPa, spectral T42, from ${esc(lab.input)}: the barotropic vorticity equation of 1950, or the shallow-water equations.</p>`);
    if (lab.running) {
      html.push(`<div class="progress lab-progress"><span style="width:${Math.round(lab.progress * 100)}%"></span></div>`);
    } else {
      html.push(`<div class="runs">${kinds.map(([k, label]) => `<button type="button" class="run" data-run-lab="${k}">Run ${label.toLowerCase()}</button>`).join('')}</div>`);
    }
    // Runs of the laboratory's models, and AI models run here (products with a
    // verification), compared at the same start when they share one.
    const columns = kinds.filter(([k]) => lab.runs[k] && lab.runs[k].length).map(([k, label]) => [label, lab.runs[k]]);
    const ai = (extra.ai || []).filter((a) => !columns.length || a.start === lab.start);
    if (!columns.length && ai.length) html.push(`<p class="note">${esc(ai[0].label)}, run here from the GFS analysis of ${esc(new Date(ai[0].start).toISOString().slice(0, 13).replace('T', ' '))} UTC:</p>`);
    columns.push(...ai.slice(0, columns.length ? 3 : 1).map((a) => [a.short, a.leads]));
    if (columns.length) {
      const leads = columns[0][1].filter((r) => r.hours % 24 === 0).map((r) => r.hours);
      const cell = (rows, h) => { const r = rows.find((x) => x.hours === h); return r ? `${Math.round(r.model)} m` : '—'; };
      const persistence = (h) => { const r = columns[0][1].find((x) => x.hours === h); return r ? `${Math.round(r.persistence)} m` : '—'; };
      html.push(`<table class="verify"><thead><tr><th>Lead</th>${columns.map(([label]) => `<th>${esc(label)}</th>`).join('')}<th>Persistence</th></tr></thead><tbody>${leads
        .map((h) => `<tr><th>+${h} h</th>${columns.map(([, rows]) => `<td>${cell(rows, h)}</td>`).join('')}<td>${persistence(h)}</td></tr>`).join('')}</tbody></table>`);
      const all = columns.flatMap(([, rows]) => rows);
      html.push(`<p class="note">Root-mean-square error of height, 20–80° N, against ${esc(lab.against || 'GFS')}${all.some((r) => r.against === 'forecast') ? ' (its forecast, beyond the present)' : ' analyses'}.</p>`);
    }
    if (lab.error) html.push(`<p class="note">${esc(lab.error)}</p>`);
    html.push('</div>');
  }

  // ---- map ----
  html.push('<h3>Map</h3>');
  const styles = [['relief', 'Relief'], ...extra.imagery.map((e) => [e.name, e.name === 'true-colour' ? `Satellite, ${e.date ? e.date.slice(5).replace('-', '/') : 'day'}` : e.name === 'blue-marble' ? 'Blue Marble' : e.title]), ['plain', 'Plain']];
  html.push(`<div class="segmented wrap">${styles.map(([s, label]) => `<button type="button" data-style="${esc(s)}" aria-pressed="${state.style === s}">${esc(label)}</button>`).join('')}</div>`);
  const checks = [['night', 'Day and night'], ...(extra.lights ? [['lights', 'Lights at night']] : []), ...(extra.clouds ? [['clouds', 'Forecast clouds']] : []),
    ['coast', 'Coastlines'], ['borders', 'Borders'], ['rivers', 'Rivers'], ['graticule', 'Graticule'], ['labels', 'Names'], ['reliefOnData', 'Relief on data']];
  html.push(`<div class="checks">${checks.map(([k, l]) => `<label class="check"><input type="checkbox" data-base="${k}"${state[k] ? ' checked' : ''}> ${l}</label>`).join('')}</div>`);
  html.push(`<label class="select"><span>Names in</span><select id="language">${LANGUAGES.map(([c, n]) => `<option value="${c}"${c === state.lang ? ' selected' : ''}>${n}</option>`).join('')}</select></label>`);
  el.innerHTML = html.join('');

  el.querySelectorAll('input[name="surface"]').forEach((i) => i.addEventListener('change', () => on.surface(i.value)));
  el.querySelectorAll('input[name="flow"]').forEach((i) => i.addEventListener('change', () => on.flow(i.value)));
  el.querySelectorAll('input[name="lines"]').forEach((i) => i.addEventListener('change', () => on.lines(i.value)));
  el.querySelectorAll('button[data-model]').forEach((b) => b.addEventListener('click', () => on[b.dataset.model](b.dataset.id)));
  const cmp = el.querySelector('select[data-compare]');
  if (cmp) cmp.addEventListener('change', () => on.compare(cmp.value || null));
  el.querySelectorAll('input[data-plates]').forEach((i) => i.addEventListener('change', () => on.plates(i.checked)));
  el.querySelectorAll('input[data-mechanisms]').forEach((i) => i.addEventListener('change', () => on.mechanisms(i.checked)));
  el.querySelectorAll('button[data-mechanism-mode]').forEach((b) => b.addEventListener('click', () => on.mechanismMode(b.dataset.mechanismMode)));
  el.querySelectorAll('input[data-stations]').forEach((i) => i.addEventListener('change', () => on.stations(i.checked)));
  el.querySelectorAll('button[data-station-frame]').forEach((b) => b.addEventListener('click', () => on.stationFrame(b.dataset.stationFrame)));
  el.querySelectorAll('input[data-track]').forEach((i) => i.addEventListener('change', () => on.track(i.dataset.track, i.checked)));
  el.querySelectorAll('input[data-events]').forEach((i) => i.addEventListener('change', () => on.events(i.checked ? i.dataset.events : null)));
  el.querySelectorAll('button[data-events-window]').forEach((b) => b.addEventListener('click', () => on.events(b.dataset.eventsWindow)));
  el.querySelectorAll('button[data-style]').forEach((b) => b.addEventListener('click', () => on.base('style', b.dataset.style)));
  el.querySelectorAll('input[data-base]').forEach((i) => i.addEventListener('change', () => on.base(i.dataset.base, i.checked)));
  el.querySelectorAll('button[data-year]').forEach((b) => b.addEventListener('click', () => on.year(Number(b.dataset.year))));
  el.querySelectorAll('[data-run-lab]').forEach((b) => b.addEventListener('click', () => on.runLab(b.dataset.runLab)));
  const ind = el.querySelector('#indicator');
  if (ind) ind.addEventListener('change', () => on.surface(ind.value));
  el.querySelector('#language').addEventListener('change', (e) => on.base('lang', e.target.value));
}
