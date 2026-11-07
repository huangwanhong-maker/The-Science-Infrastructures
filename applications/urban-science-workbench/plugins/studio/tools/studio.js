// The studio: one's own data in space and time, and one's own code over it. A file dropped anywhere on the page, or a record of the library, goes on the map and on the line
// of time; a script, in JavaScript (run in a worker of this page), Python or C++ (run by the server, as a process with a time limit), is given what is shown and answers what to
// draw, with tables, curves and notes here. Every script obeys one contract in its three languages: scripts/README.md, and the examples.

import { spacetime as ST } from '/visualization/js/gis/index.js';
import { SOURCES } from '../library.js';
import { TEMPLATES } from './templates.js';

const API = '/studio';
const STUDIO_HEADER = { 'Content-Type': 'application/json', 'X-Workbench': 'studio' };
async function post(path, body) {
  const r = await fetch(`${API}${path}`, { method: 'POST', headers: STUDIO_HEADER, body: JSON.stringify(body) });
  const a = await r.json().catch(() => ({ message: `the server answered ${r.status}` }));
  if (!r.ok) throw new Error(a.message || `the server answered ${r.status}`);
  return a;
}
const LANGUAGE = { js: 'JavaScript', py: 'Python', cpp: 'C++' };
const languageOf = (file) => file.split('.').pop();

// A script in JavaScript, run off the page's thread: its own text and a small runner, with the contract's module at the same address the page has it.
function runInWorker(code, message, { onProgress, seconds = 120 } = {}) {
  const head = `import * as api from ${JSON.stringify(`${location.origin}/visualization/js/gis/spacetime.js`)};\n`;
  const tail = `
;self.onmessage = async ({ data }) => {
  try {
    const plugin = create(api);
    if (data.describe) { self.postMessage({ description: plugin.describe() }); return; }
    const context = { ...data.context, progress: (fraction, text = '') => self.postMessage({ progress: fraction, text }) };
    const params = {};
    for (const p of plugin.params || []) params[p.name] = data.params && data.params[p.name] !== undefined ? data.params[p.name] : p.default;
    const results = await plugin.run(data.records, params, context);
    self.postMessage({ results });
  } catch (error) { self.postMessage({ error: String(error && error.stack || error) }); }
};`;
  const url = URL.createObjectURL(new Blob([head, code, tail], { type: 'text/javascript' }));
  return new Promise((resolve, reject) => {
    const worker = new Worker(url, { type: 'module' });
    const stop = setTimeout(() => { worker.terminate(); reject(new Error(`stopped after ${seconds} s`)); }, seconds * 1000);
    const end = () => { clearTimeout(stop); worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = (e) => { end(); reject(new Error(`line ${(e.lineno || 2) - 1}: ${e.message}`)); };
    worker.onmessage = ({ data }) => {
      if (data.progress !== undefined) return onProgress?.(data.progress, data.text);
      end();
      if (data.error) reject(new Error(data.error)); else resolve(data);
    };
    worker.postMessage(message);
  });
}

export async function render(root, ctx) {
  const { dom: { h }, chart, earth, st } = ctx;
  const say = {};
  // ---- data -----------------------------------------------------------------------------------------------------------------------------------------
  const held = h('div', { class: 'st-held' });
  const mapping = h('div', {});
  const picker = h('input', { type: 'file', accept: '.csv,.tsv,.txt,.json,.geojson,.jsonl', hidden: true, onchange: () => { if (picker.files[0]) takeFile(picker.files[0]); picker.value = ''; } });
  const drop = h('div', { class: 'st-drop', onclick: () => picker.click() }, h('b', {}, 'Drop a file anywhere on the page'), h('span', {}, 'A table with a place and a time (CSV, TSV, JSON lines), GeoJSON, or a record of the contract'), picker);
  const library = h('div', { class: 'st-library' }, Object.entries(SOURCES).filter(([key]) => key !== 'walkers-week' || ctx.area?.id === 'melbourne').map(([key, s]) => h('button', { class: 'choice', onclick: () => fromLibrary(key) }, h('b', {}, s.title), h('span', {}, s.credit))));
  say.data = h('p', { class: 'muted' });

  function drawHeld() {
    const list = st.entries();
    held.replaceChildren(...list.map((e) => h('div', { class: 'st-record' },
      h('b', {}, e.rec.title),
      h('span', {}, `${e.rec.kind}, ${e.rec.count.toLocaleString('en')}`),
      h('span', { class: 'mono' }, `${ST.formatTime(e.rec.span[0], 'full', 3 * ST.YEAR)} – ${ST.formatTime(e.rec.span[1], 'full', 3 * ST.YEAR)}`),
      h('button', { class: 'st-go', title: 'go to it', onclick: () => st.fit(e.id) }, 'go to it'))));
    if (!list.length) held.replaceChildren(h('p', { class: 'faint' }, 'Nothing is shown yet.'));
  }
  async function fromLibrary(key) {
    const s = SOURCES[key];
    say.data.textContent = `reading ${s.title.toLowerCase()}…`;
    try {
      const raw = await s.read();
      st.show(raw, { ...s.style, credit: s.credit });
      st.fit();
      st.timeline.play();
      say.data.textContent = '';
    } catch (e) { say.data.textContent = String(e.message || e); }
    drawHeld();
  }
  function showRecords(records, name) {
    const left = {};
    for (const r of records) for (const [why, n] of Object.entries(st.show(r).record.left || {})) left[why] = (left[why] || 0) + n;
    st.fit();
    const leftOut = Object.entries(left).map(([why, n]) => ` ${n.toLocaleString('en')} left out: ${why}.`).join('');
    say.data.textContent = `${name}: ${records.length === 1 ? 'on the map' : `${records.length} records on the map`}.${leftOut}`;
    drawHeld();
  }
  // A table is asked what its columns are, only if they cannot be told.
  function askColumns(table, title) {
    const g = ST.guessColumns(table);
    const choice = (label, key) => {
      const s = h('select', { class: 'st-select', onchange: () => { g[key] = s.value || null; } }, h('option', { value: '' }, '—'), table.columns.map((c) => h('option', { value: c, selected: c === g[key] ? true : null }, c)));
      return h('label', { class: 'st-map' }, h('span', {}, label), s);
    };
    const draw = () => {
      try {
        const record = ST.recordFromTable(table, g, title);
        mapping.replaceChildren();
        showRecords([record], title);
      } catch (e) { say.data.textContent = String(e.message || e); }
    };
    if (g.lat && g.lon && g.time) return draw();
    say.data.textContent = `${title}: ${table.rows.length.toLocaleString('en')} rows. Which columns say where and when?`;
    mapping.replaceChildren(h('div', { class: 'st-mapping' }, choice('latitude', 'lat'), choice('longitude', 'lon'), choice('time', 'time'), choice('which thing (for tracks)', 'track'), h('button', { class: 'primary', onclick: draw }, 'Draw it')));
  }
  async function takeFile(file) {
    say.data.textContent = `reading ${file.name}…`;
    try {
      const out = ST.readRecords(await file.text(), file.name);
      if (out.table) askColumns(out.table, out.title);
      else showRecords(out.records, out.title);
    } catch (e) { say.data.textContent = String(e.message || e); }
  }
  // Files dropped on the page come here, also one dropped before this panel was open.
  ctx.onDispose(st.acceptFiles(takeFile));
  st.onChange(drawHeld);

  // ---- code -----------------------------------------------------------------------------------------------------------------------------------------
  const lines = h('pre', { class: 'st-lines', 'aria-hidden': 'true' });
  const editor = h('textarea', { class: 'st-code', spellcheck: 'false', 'aria-label': 'the script', wrap: 'off' });
  const countLines = () => { lines.textContent = Array.from({ length: editor.value.split('\n').length }, (_, k) => k + 1).join('\n'); lines.scrollTop = editor.scrollTop; };
  editor.addEventListener('input', countLines);
  editor.addEventListener('scroll', () => { lines.scrollTop = editor.scrollTop; });
  editor.addEventListener('keydown', (e) => {
    if (e.key === 'Tab') { e.preventDefault(); const { selectionStart: a, selectionEnd: b } = editor; editor.setRangeText('  ', a, b, 'end'); countLines(); }
    if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) { e.preventDefault(); run(); }
  });
  const name = h('input', { class: 'st-name', 'aria-label': 'the script\'s name', placeholder: 'its name' });
  const langKeys = h('div', { class: 'keys-tight' }, Object.entries(LANGUAGE).map(([k, n]) => h('button', { 'data-lang': k, onclick: () => setLanguage(k, true) }, n)));
  const scripts = h('select', { class: 'st-select st-scripts', 'aria-label': 'a script', onchange: () => open(scripts.value) });
  const form = h('div', { class: 'st-form' });
  const output = h('div', { class: 'st-output' });
  const status = h('span', { class: 'st-status' });
  const runKey = h('button', { class: 'primary', title: 'Ctrl Enter', onclick: () => run() }, 'Run');
  const saveKey = h('button', { title: 'keep it in plugins/studio/scripts', onclick: () => save() }, 'Keep');
  let language = 'py', values = {}, description = null, shownByRun = [];

  function setLanguage(k, fresh) {
    language = k;
    for (const b of langKeys.children) b.classList.toggle('on', b.dataset.lang === k);
    if (fresh) {
      editor.value = TEMPLATES[k];
      name.value = `mine.${k}`;
      scripts.value = '';
      countLines();
      describe();
    }
  }
  async function refreshScripts(selected) {
    try {
      const a = await (await fetch(`${API}/scripts`)).json();
      scripts.replaceChildren(h('option', { value: '' }, 'a new script'), ...a.scripts.map((s) => h('option', { value: s.file, selected: s.file === selected ? true : null }, `${s.file} (${LANGUAGE[languageOf(s.file)]})`)));
      if (a.refusal) say.code = a.refusal;
    } catch { /* the page still has JavaScript */ }
  }
  async function open(file) {
    if (!file) return setLanguage(language, true);
    const r = await fetch(`${API}/script?file=${encodeURIComponent(file)}`);
    const a = await r.json();
    if (!r.ok) { status.textContent = a.message; return; }
    editor.value = a.code;
    name.value = file;
    setLanguage(languageOf(file), false);
    countLines();
    describe();
  }
  function drawForm() {
    form.replaceChildren();
    values = {};
    for (const p of description?.params || []) {
      values[p.name] = p.default;
      const input = p.kind === 'choice'
        ? h('select', { class: 'st-select', onchange: () => { values[p.name] = input.value; } }, (p.options || []).map((o) => h('option', { value: o, selected: o === p.default ? true : null }, o)))
        : p.kind === 'boolean'
          ? h('input', { type: 'checkbox', checked: p.default ? true : null, onchange: () => { values[p.name] = input.checked; } })
          : h('input', { type: 'number', value: p.default, min: p.min ?? undefined, max: p.max ?? undefined, step: p.step ?? (p.kind === 'integer' ? 1 : 'any'), onchange: () => { values[p.name] = Number(input.value); } });
      form.append(h('label', { class: 'st-param', title: p.help || '' }, h('span', {}, p.label), input));
    }
  }
  async function describe() {
    description = null;
    try {
      if (language === 'js') { description = (await runInWorker(editor.value, { describe: true }, { seconds: 15 })).description; drawForm(); status.textContent = ''; } else if (name.value && scripts.value) {
        if (language === 'cpp') status.textContent = 'compiling it, the first time…';
        const r = await fetch(`${API}/describe?file=${encodeURIComponent(scripts.value)}`);
        const a = await r.json();
        if (r.ok) { description = a; drawForm(); status.textContent = ''; } else status.textContent = a.message;
      } else { form.replaceChildren(); }
    } catch (e) { status.textContent = String(e.message || e); }
  }
  const contextNow = () => {
    const v = earth.view, DEG = 180 / Math.PI;
    return { lat: v.lat * DEG, lon: v.lon * DEG, range_m: v.range, now_ms: st.now(), window_ms: Number.isFinite(st.window()) ? st.window() : null, area: ctx.area ? { id: ctx.area.id, title: ctx.area.title || ctx.area.name, lat: ctx.area.lat, lon: ctx.area.lon } : null };
  };
  async function save() {
    const file = name.value.trim() || `mine.${language}`;
    try {
      const a = await post('/script', { file, code: editor.value });
      name.value = a.file;
      await refreshScripts(a.file);
      status.textContent = `kept as ${a.file}`;
      return a.file;
    } catch (e) { status.textContent = String(e.message || e); output.replaceChildren(h('pre', { class: 'error-box st-error' }, String(e.message || e))); return null; }
  }
  function showResults(results) {
    output.replaceChildren();
    for (const id of shownByRun) st.remove(id);
    shownByRun = [];
    for (const r of results) {
      try {
        if (ST.KINDS.includes(r.kind)) shownByRun.push(st.show(r).id);
        else if (r.kind === 'note') output.append(h('p', {}, r.text));
        else if (r.kind === 'table') output.append(h('h2', {}, r.title), h('table', { class: 'sources' }, h('thead', {}, h('tr', {}, r.columns.map((c) => h('th', {}, c)))), h('tbody', {}, r.rows.slice(0, 200).map((row) => h('tr', {}, row.map((c) => h('td', {}, String(c ?? '—'))))))), r.note ? h('p', { class: 'note' }, r.note) : null);
        else if (r.kind === 'series') {
          const lines_ = r.lines.map((l, k) => ({ name: l.name, x: l.x, y: l.y, color: chart.category(k) }));
          output.append(h('h2', {}, r.title), chart.lines(lines_, { xlabel: r.xlabel, ylabel: r.ylabel, width: 460, height: 200, xfmt: r.time_axis ? (t) => ST.formatTime(t, 'full', 3 * ST.YEAR) : undefined }), r.lines.length > 1 ? chart.legend(lines_) : null, r.note ? h('p', { class: 'note' }, r.note) : null);
        }
      } catch (e) { output.append(h('p', { class: 'error-box' }, `${r.title || r.kind}: ${e.message || e}`)); }
    }
    if (shownByRun.length) st.timeline.play();
    drawHeld();
  }
  async function run() {
    const t0 = performance.now();
    runKey.disabled = true;
    status.textContent = 'running…';
    output.replaceChildren();
    const records = st.raw();
    try {
      let results;
      if (language === 'js') {
        results = (await runInWorker(editor.value, { records, params: values, context: contextNow() }, { onProgress: (f, t) => { status.textContent = `${Math.round(f * 100)} %${t ? `: ${t}` : ''}`; } })).results;
      } else {
        const file = await save();
        if (!file) return;
        scripts.value = file;
        status.textContent = language === 'cpp' ? 'compiling, the first time, and running…' : 'running…';
        const { job } = await post('/run', { file, records, params: values, context: contextNow() });
        results = await new Promise((resolve, reject) => {
          const poll = async () => {
            const a = await (await fetch(`${API}/job/${job}`)).json();
            if (a.state === 'done') return resolve(a.results);
            if (a.state === 'failed' || a.state === 'cancelled') return reject(new Error(a.message || a.state));
            status.textContent = `${Math.round((a.progress || 0) * 100)} %${a.text ? `: ${a.text}` : ''}`;
            setTimeout(poll, 300);
          };
          poll();
        });
        describe();
      }
      showResults(results || []);
      status.textContent = `done in ${((performance.now() - t0) / 1000).toFixed(1)} s`;
    } catch (e) {
      status.textContent = 'it failed';
      output.replaceChildren(h('pre', { class: 'error-box st-error' }, String(e.message || e)));
    } finally {
      runKey.disabled = false;
    }
  }

  // ---- the sheet -------------------------------------------------------------------------------------------------------------------------------------
  root.append(
    h('h1', {}, 'Your data and your code'),
    h('p', { class: 'muted' }, 'Anything with a place and a time goes on this map, and runs through time above the clock. A script is given what is shown, and answers what to draw.'),
    h('h2', {}, 'What is shown'), held, drop, say.data, mapping,
    h('div', { class: 'lens-group' }, 'Or begin from a record of the library'), library,
    h('h2', {}, 'A script'),
    h('div', { class: 'st-bar' }, scripts, langKeys),
    h('div', { class: 'st-editor' }, lines, editor),
    form,
    h('div', { class: 'st-bar' }, runKey, saveKey, name, status),
    output,
    h('details', { class: 'st-contract' }, h('summary', {}, 'What a script is given and answers'),
      h('p', {}, 'Given: the records shown, the parameters it declares, and the view (where it looks, the moment and its window). Answers: records, tables, curves and notes.'),
      h('p', {}, 'Space is degrees on WGS 84, and time is milliseconds since 1970 in UTC. A record is one of:'),
      h('ul', {},
        h('li', {}, h('b', {}, 'events'), ', points in space and time; ', h('code', {}, 'time, lat, lon, values')),
        h('li', {}, h('b', {}, 'tracks'), ', moving things; ', h('code', {}, 'tracks: [{id, time, lat, lon, values}]')),
        h('li', {}, h('b', {}, 'cells'), ', places with a series; ', h('code', {}, 'time, cells: [{id, ring}], values: {name: [[…] at each time]}')),
        h('li', {}, h('b', {}, 'grid'), ', a raster through time; ', h('code', {}, 'time, box, width, height, values')),
        h('li', {}, h('b', {}, 'flows'), ', between places; ', h('code', {}, 'places, flows: [{from, to, time}]'))),
      h('p', {}, 'A plug-in of the workbench itself reaches the same line of time by ', h('code', {}, 'ctx.st.show(record)'), '. plugins/studio/scripts/README.md has the whole of it.')),
  );
  drawHeld();
  await refreshScripts('density.py');
  await open('density.py').catch(() => setLanguage('js', true));
}
