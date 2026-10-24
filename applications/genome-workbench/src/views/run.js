// The Engine view: where the computation happens. The machine and what the engine was compiled
// for; the workers the cluster keeps, here and (when configured) on other machines; a job cut
// into shards and pulled by free workers, drawn on a timeline with each shard's threads; and
// the same work on 1 to n threads. The recipe that builds the engine is here too.
//
// A calm page of the system: four tiles of its state along the top, the workers and the build
// beneath them, then the measurements, each a card that keeps its last run. The tiles and the
// workers refresh themselves while the page is open.

import { el, clear } from '../dom.js';
import * as api from '../api.js';
import { chart, timeline } from '../draw/chart.js';
import * as fmt from '../format.js';
import { CATEGORIES, token } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, seg, select, fld, card, metrics, callout, empty, toast } from '../ui/kit.js';

let view = null;
// The measurements' choices and last results, kept when the view is mounted again (a change of theme).
const kept = { job: null, scaling: null, devices: null, poll: null, jobParams: { reads: 100000, shards: 16, genome: 'sim:repeats-200k' }, targets: 4000 };

export function mount(section, ctx) {
  view = { ctx, section, e: null, target: null, setup: null, tiles: el('div', { class: 'sys-tiles' }) };
  section.append(
    el('p', { class: 'lede' }, 'The engine is the genomics and evolution domains in C++, compiled for this processor: its vector lanes as wide as the processor allows (calx’s CPU backend), its threads in calx’s scheduler. The server keeps several engines alive as workers, and can keep more on other machines; a job is cut into shards that any worker can compute from a few numbers, and merged in the order of the shards, so that its answer never depends on who computed what.'),
    view.tiles,
    el('div', { class: 'sys-cols' },
      el('div', { class: 'sys-col' }, workersCard(), opsCard()),
      el('div', { class: 'sys-col' }, buildCard())),
    el('div', { class: 'section-head' }, el('h3', {}, 'Measurements'), el('span', { class: 'sub' }, 'the same work across the workers, the threads and the devices')),
    el('div', { class: 'stack', style: { gap: '16px' } }, jobCard(), scalingCard(), devicesCard()));
  drawTiles();
  if (kept.job) drawJob();
  if (kept.scaling) drawScaling();
  if (kept.devices) drawDevices();
  clearInterval(kept.poll);
  kept.poll = setInterval(() => { if (!view.section.hidden && document.visibilityState === 'visible') quietRefresh(); }, 10000);
}

export async function show(ctx) {
  view.ctx = ctx;
  await refresh();
}

async function refresh() {
  let e;
  try { e = await api.engineState(); } catch (err) { drawTiles(err); return; }
  view.e = e;
  drawAll();
  // Ask a worker for its target, which also starts the workers.
  if (e.current) {
    const target = await api.run('target').catch(() => null);
    if (target) { view.target = target; view.e = await api.engineState().catch(() => e); }
  }
  drawAll();
}

// The state alone, without asking a worker anything: for the page's own refreshing.
async function quietRefresh() {
  const e = await api.engineState().catch(() => null);
  if (!e) return;
  view.e = e;
  drawTiles();
  drawWorkers();
}

function drawAll() { drawTiles(); drawWorkers(); drawOps(); if (view.setup) drawBuild(view.setup); }

// ---- the tiles ----------------------------------------------------------------------------------------
const brand = (s) => (s || '').replace(/\(R\)|\(TM\)/g, '').replace(/ CPU\b/, '').replace(/\s*@.*$/, '').replace(/\s+/g, ' ').trim();
const clock = (s) => (s || '').match(/@\s*(.*)$/)?.[1] || '';

function tile(ic, label, value, sub) {
  return el('div', { class: 'sys-tile' }, el('div', { class: 'glyph' }, icon(ic)), el('div', { class: 'k' }, label), el('div', { class: 'v' }, value), el('div', { class: 's' }, sub));
}
const flag = (text, on) => el('span', { class: `badge ${on ? 'good' : ''}`, title: on ? 'this processor has it' : 'this processor lacks it' }, on ? icon('check', 'sm') : null, text);

function drawTiles(err = null) {
  const box = clear(view.tiles), e = view.e, t = view.target;
  if (err) { box.append(callout(`The server did not answer: ${err.message}`, 'bad')); return; }
  if (!e) { for (let i = 0; i < 4; i++) box.append(el('div', { class: 'skeleton', style: { height: '104px', borderRadius: 'var(--r-3)' } })); return; }
  const workers = e.cluster?.workers || [];
  const alive = workers.filter((w) => w.alive);
  const tg = t || e.target || alive[0] || {};
  const best = tg.best || tg.built_for || '';
  const [dot, state, why] = !e.built ? ['off', 'Not built', 'build it with the recipe, below'] : e.current ? ['ok', 'Current', 'built from the sources as they stand'] : ['warn', 'Sources changed', 'build it again to bring the changes in'];
  const built = e.built_at ? `built ${e.built_at.slice(0, 16).replace('T', ' ')} UTC${e.toolchain ? `, g++ in ${e.toolchain === 'wsl' ? 'WSL' : e.toolchain}` : ''}` : why;
  const served = workers.reduce((a, w) => a + (w.served || 0), 0), busy = workers.reduce((a, w) => a + (w.seconds || 0), 0);
  // The processor's model as the value, its maker and generation beside the threads.
  const full = brand(tg.brand), model = full.replace(/^\d+(st|nd|rd|th) Gen\s+/i, '').replace(/^(Intel|AMD|Apple)\s+/i, ''), maker = full.slice(0, full.length - model.length).trim();
  box.append(
    tile('run', 'Engine', [el('span', { class: `dot ${dot}` }), el('span', {}, state)], e.built ? (e.current ? built : `${why}; ${built}`) : why),
    tile('gpu', 'Processor', model || tg.architecture || '–',
      el('span', { class: 'flags' }, el('span', {}, [maker, clock(tg.brand), `${tg.threads ?? '–'} threads a worker`, t?.hardware_threads ? `${t.hardware_threads} on the machine` : null].filter(Boolean).join(' · ')),
        flag('SSE4.1', true), flag('AVX2', t ? t.avx2 : /avx2|avx512/.test(best)), flag('AVX-512', t ? t.avx512 : /avx512/.test(best)))),
    tile('layers', 'Vector lanes', tg.lane_bytes ? `${tg.lane_bytes} bytes` : '–', tg.built_for ? `compiled for ${tg.built_for}: ${tg.lane_bytes} lanes of 8 bits, ${tg.lane_bytes / 2} of 16` : 'the width the engine was compiled for'),
    tile('activity', 'Workers', `${alive.length} alive${workers.length > alive.length ? ` of ${workers.length}` : ''}`, workers.length ? `${fmt.count(served)} requests served · ${fmt.number(busy, 3)} s of work` : 'they start with the first request'));
}

// ---- the workers ----------------------------------------------------------------------------------------
function workersCard() {
  view.workerTable = el('div', { class: 'table-wrap' });
  view.config = el('code', {});
  const howto = el('div', { class: 'howto', hidden: true },
    el('p', { class: 'note' }, 'Build the engine there from the same sources (the recipe), then write the command that starts it into the cluster file, and restart the workers. A remote worker takes the shards of jobs on simulated genomes; genomes read from this computer’s files stay with the local workers.'),
    el('pre', { class: 'file' }, '{\n  "workers": [\n    {"name": "eu_region_server", "command": ["ssh", "eu_region_server",\n      "~/The-Science-Infrastructures/data/cache/engines/genome/genome --serve"]}\n  ]\n}'));
  const more = btn('Adding a worker on another machine', { icon: 'chevron', kind: 'ghost', size: 'sm', onClick: () => { howto.hidden = !howto.hidden; more.classList.toggle('open', !howto.hidden); } });
  const restart = btn('Restart', { icon: 'refresh', size: 'sm', title: 'stop the workers; they start again with the next request', onClick: async () => {
    try { await api.restart(); toast('The workers were stopped; they start again with the next request.', 'good'); refresh(); } catch (err) { toast(err.message, 'error'); }
  } });
  return card({ title: 'The workers', sub: 'engines kept alive by the server, each with its pool of threads', actions: [restart],
    body: [view.workerTable, el('div', { class: 'card-foot' }, el('div', { class: 'small faint' }, view.local = el('span', {}), ' remote ones from ', view.config), more, howto)] });
}

function drawWorkers() {
  const e = view.e, box = clear(view.workerTable);
  const workers = e?.cluster?.workers || [];
  view.local.textContent = `${e?.cluster?.local ?? 2} local;`;
  view.config.textContent = e?.cluster?.config || 'cluster.json';
  if (!workers.length) { box.append(el('p', { class: 'note', style: { padding: '12px 16px' } }, e ? 'No worker has started yet; they start with the first request.' : '…')); return; }
  const table = el('table', { class: 'data wide' }, el('tr', {}, ...[['worker'], ['where'], ['threads', 1], ['lanes'], ['requests', 1], ['busy', 1], ['state']].map(([h, n]) => el('th', { class: n ? 'num' : '' }, h))));
  workers.forEach((x, k) => {
    table.append(el('tr', {},
      el('td', {}, el('span', { class: 'swatch', style: { background: CATEGORIES[k % CATEGORIES.length] } }), ' ', x.name),
      el('td', {}, x.remote ? 'remote' : 'this computer'),
      el('td', { class: 'num' }, x.threads ?? '–'),
      el('td', {}, x.built_for ? `${x.built_for}, ${x.lane_bytes} B` : '–'),
      el('td', { class: 'num' }, fmt.count(x.served)),
      el('td', { class: 'num' }, `${fmt.number(x.seconds, 3)} s`),
      el('td', {}, !x.alive ? el('span', { class: 'badge bad', title: x.error || '' }, 'stopped') : x.busy ? el('span', { class: 'badge accent' }, el('span', { class: 'spin' }), 'working') : el('span', { class: 'badge good' }, 'alive'))));
  });
  box.append(table);
}

// ---- the operations ----------------------------------------------------------------------------------------
const DOMAINS = { align: 'alignment', assembly: 'assembly', expression: 'expression', genome: 'genomes', index: 'indexes', pop: 'populations', reads: 'reads', run: 'the machine', target: 'the machine', trees: 'trees', tool: 'the Toolbox' };

function opsCard() {
  view.opsBody = el('div');
  view.opsCard = card({ title: 'What the engine does', sub: 'its operations, by the part of the work they serve', body: view.opsBody });
  return view.opsCard;
}

function drawOps() {
  const t = view.target, box = clear(view.opsBody);
  view.opsCard.hidden = !t?.operations;
  if (!t?.operations) {
    box.append(el('p', { class: 'note', style: { padding: '12px 16px 14px' } }, view.e?.built ? 'Listed when the engine is current: a worker is asked for them then.' : 'Nothing yet: the engine is not built.'));
    return;
  }
  const groups = new Map();
  for (const op of t.operations) {
    const d = DOMAINS[op.split('.')[0]] || op.split('.')[0];
    if (!groups.has(d)) groups.set(d, []);
    groups.get(d).push(op);
  }
  view.opsCard.querySelector('.card-head .sub').textContent = `${t.operations.length} operations, by the part of the work they serve`;
  const list = el('div', { class: 'ops' });
  for (const [d, ops] of groups) list.append(el('div', { class: 'k' }, d), el('div', { class: 'chips' }, ...ops.map((op) => el('span', { class: 'kind', title: op }, op.includes('.') ? op.slice(op.indexOf('.') + 1) : op))));
  box.append(list);
}

// ---- the build ----------------------------------------------------------------------------------------
// The recipe (recipe.json through the server's /setup): its steps, a button that runs them, and the
// log while they run.
function buildCard() {
  view.buildBody = el('div', { class: 'build' }, el('div', { class: 'skeleton', style: { height: '120px', margin: '14px 16px' } }));
  api.setup().then((v) => { drawBuild(v); if (v.running) pollSetup(); }).catch((err) => clear(view.buildBody).append(el('div', { style: { padding: '14px 16px' } }, callout(err.message, 'bad'))));
  return card({ title: 'The build', sub: 'the recipe: the engine compiled for this processor, and the genomes fetched', body: view.buildBody });
}

function pollSetup() {
  clearTimeout(view.setupTimer);
  view.setupTimer = setTimeout(async () => {
    const v = await api.setup().catch(() => null);
    if (!v) return;
    drawBuild(v);
    if (v.running) pollSetup();
    else if (v.state === 'done') refresh();
  }, 1200);
}

function drawBuild(v) {
  view.setup = v;
  const box = clear(view.buildBody), e = view.e;
  const steps = v.steps || [];
  const allDone = steps.length && steps.every((s) => s.done);
  if (e?.built && !e.current && !v.running) box.append(el('div', { class: 'pad' }, callout('The engine’s sources have changed since it was built. Building it again brings the changes in, in about two minutes; the workers then start on the new engine.', 'warn')));
  box.append(el('ol', { class: 'steps' }, ...steps.map((s, k) => el('li', { class: s.done ? 'done' : '' },
    el('span', { class: 'mark' }, v.running && !s.done && steps.slice(0, k).every((x) => x.done) ? el('span', { class: 'spin' }) : s.done ? icon('check', 'sm') : String(k + 1)),
    el('span', { class: 't' }, s.title),
    el('span', { class: 'when' }, s.done ? (s.updated ? new Date(s.updated * 1000).toLocaleString([], { dateStyle: 'medium', timeStyle: 'short' }) : 'done') : 'to do')))));
  const go = btn(v.running ? 'Running the recipe…' : allDone ? 'Build again' : 'Build the engine and fetch the genomes', { icon: v.running ? null : allDone ? 'refresh' : 'play', kind: !allDone || (e?.built && !e.current) ? 'primary' : '', disabled: v.running || !v.writable,
    onClick: async () => { try { drawBuild(await api.runSetup(allDone)); pollSetup(); } catch (err) { toast(err.message, 'error'); } } });
  if (v.running) go.prepend(el('span', { class: 'spin' }));
  if (v.log) {
    const log = el('pre', { class: 'file log' }, v.log);
    box.append(el('div', { class: 'pad last' }, log));
    requestAnimationFrame(() => { log.scrollTop = log.scrollHeight; });
  }
  if (v.state === 'failed') box.append(el('div', { class: 'pad last' }, callout(`The recipe failed at: ${(v.failed || []).join(', ')}.`, 'bad')));
  box.append(el('div', { class: 'run-bar' }, go, !v.writable && v.refusal ? el('span', { class: 'small faint', style: { alignSelf: 'center' } }, v.refusal) : null));
}

// ---- running a measurement ----------------------------------------------------------------------------------
// The button disabled and turning, a line with the seconds since the run began at the top of the
// card's results, and the last results dimmed until the new ones come. Returns the function that ends it.
function working(button, body, text) {
  const keep = [...button.childNodes];
  button.disabled = true;
  button.replaceChildren(el('span', { class: 'spin' }), el('span', {}, 'Running…'));
  const elapsed = el('span', { class: 'elapsed' });
  const line = el('div', { class: 'wait pad-x' }, el('span', { class: 'spin' }), el('span', {}, text), elapsed);
  const old = [...body.children];
  if (old.length && !old[0].classList.contains('empty')) old.forEach((c) => c.classList.add('is-busy'));
  else clear(body);
  body.prepend(line);
  const t0 = performance.now();
  const tick = () => { elapsed.textContent = `${((performance.now() - t0) / 1000).toFixed(1)} s`; };
  tick();
  const timer = setInterval(tick, 100);
  return (error = null) => {
    clearInterval(timer);
    button.replaceChildren(...keep);
    button.disabled = false;
    line.remove();
    [...body.children].forEach((c) => c.classList.remove('is-busy'));
    if (error) body.prepend(el('div', { class: 'pad' }, callout(error.message, 'bad')));
  };
}

const width = (body) => Math.max(480, (body.clientWidth || view.section.clientWidth - 60) - 34);

// ---- a job ----------------------------------------------------------------------------------------
function jobCard() {
  const p = kept.jobParams;
  view.jobBtn = btn('Run the job', { icon: 'play', kind: 'primary', onClick: runJob });
  view.jobBody = el('div', {}, empty({ icon: 'activity', title: 'No job run yet', text: 'Short reads simulated from a genome, each read from its own random stream, so that a shard is just which reads; their qualities and errors are measured. Each shard is a task graph of chunks on a worker’s threads.' }));
  const strip = el('div', { class: 'ctl-strip' },
    fld('reads', seg([20000, 100000, 300000, 1000000].map((n) => [n, n >= 1e6 ? `${n / 1e6} M` : `${n / 1000} k`, `${fmt.count(n)} reads`]), p.reads, (v) => { p.reads = v; })),
    fld('shards', seg([1, 4, 8, 16, 32, 64].map((n) => [n, String(n)]), p.shards, (v) => { p.shards = v; })),
    fld('genome', select([['sim:repeats-200k', 'a simulated genome (200 kb)'], ['sim:bacterium-5m', 'a simulated bacterium (5 Mb)'], ['current', 'the genome in hand']], p.genome, (v) => { p.genome = v; })),
    el('span', { class: 'grow' }), view.jobBtn);
  return card({ title: 'A job across the workers', sub: 'reads simulated in shards, pulled by whichever worker is free, merged in the order of the shards', body: [strip, view.jobBody] });
}

async function runJob() {
  const p = kept.jobParams;
  const genome = p.genome === 'current' ? view.ctx.genome : p.genome;
  const stop = working(view.jobBtn, view.jobBody, `Simulating ${fmt.count(p.reads)} reads in ${p.shards} shards across the workers…`);
  const t0 = performance.now();
  try {
    const r = await api.job('reads.sample', { genome, reads: p.reads, shards: p.shards, seed: 7, length: 150, error_first: 0.001, error_last: 0.012 });
    kept.job = { r, wall: (performance.now() - t0) / 1000 };
    stop();
    drawJob();
    refresh();
  } catch (e) { stop(e); }
}

function drawJob() {
  const { r, wall } = kept.job, box = clear(view.jobBody), W = Math.min(1100, width(view.jobBody));
  const job = r._job;
  const workers = [...new Set(job.trace.map((t) => t.worker).concat(job.merge.worker ? [job.merge.worker] : []))].sort();
  const colour = (name) => CATEGORIES[workers.indexOf(name) % CATEGORIES.length];
  const engineSeconds = job.trace.reduce((s, t) => s + (t.engine_seconds || 0), 0);
  box.append(el('div', { class: 'pad' },
    metrics([[fmt.count(r.reads), 'reads'], [fmt.count(r.bases), 'bases'], [`${job.shards}`, `shards on ${workers.length} worker${workers.length > 1 ? 's' : ''}`], [fmt.seconds(job.seconds), 'from plan to merge'],
      [fmt.seconds(engineSeconds), 'of engine time in the shards'], [fmt.percent(r.error_per_base, 3), 'of bases in error'], [fmt.count(r.minimizers), 'minimizers sketch them']]),
    el('p', { class: 'note', style: { margin: '-4px 0 0' } }, `${fmt.seconds(wall)} in all, seen from this page.`)));
  // The timeline: one row per worker, the shards as bars, the merge last.
  const rows = workers;
  const bars = job.trace.map((t) => ({ row: rows.indexOf(t.worker), start: t.start, end: t.end, colour: colour(t.worker), opacity: 0.55 + 0.4 * ((t.shard % 2) === 0), title: `shard ${t.shard}: ${fmt.count(t.size)} reads, ${fmt.seconds(t.engine_seconds)} in the engine` }));
  if (job.merge.worker) bars.push({ row: rows.indexOf(job.merge.worker), start: job.merge.start, end: Math.max(job.merge.end, job.merge.start + job.seconds * 0.004), colour: token('--ink'), title: 'the merge' });
  const t = el('div', { class: 'fig' }); t.innerHTML = timeline({ width: W, rows, bars, seconds: job.seconds });
  box.append(el('div', { class: 'block' }, el('div', { class: 'eyebrow' }, 'shards on the workers, in seconds from the plan'), t));
  // One shard's threads: the busiest shard's task graph.
  const big = job.trace.filter((x) => x.threads?.tasks?.length).sort((a, b) => b.threads.tasks.length - a.threads.tasks.length)[0];
  if (big) {
    const th = big.threads;
    const trows = [...Array(th.threads).keys()];
    const tbars = th.tasks.map(([task, thread, s, e]) => ({ row: thread, start: s, end: e, colour: colour(big.worker), opacity: 0.5 + 0.45 * (task % 2), title: `chunk ${task}` }));
    const tt = el('div', { class: 'fig' }); tt.innerHTML = timeline({ width: W, rows: trows, bars: tbars, seconds: th.seconds, rowHeight: 14, label: (k) => `thread ${k}` });
    box.append(el('div', { class: 'block' }, el('div', { class: 'eyebrow' }, `inside shard ${big.shard} on ${big.worker}: its chunks on ${th.threads} threads, ${fmt.percent(th.utilisation, 0)} busy`), tt));
  }
  // What was measured: the error rate and mean quality along the read.
  const n = r.error_rate.length;
  const a = el('div', { class: 'fig' }); a.innerHTML = chart({ width: 440, height: 220, x: [1, n], y: [0, Math.max(...r.error_rate) * 1.15 || 0.01], xlabel: 'position in the read', ylabel: 'errors per base', series: [{ points: r.error_rate.map((v, k) => [k + 1, v]), colour: token('--accent'), area: true }] });
  const b = el('div', { class: 'fig' }); b.innerHTML = chart({ width: 440, height: 220, x: [1, n], y: [Math.min(...r.mean_quality) - 2, Math.max(...r.mean_quality) + 2], xlabel: 'position in the read', ylabel: 'mean Phred quality', series: [{ points: r.mean_quality.map((v, k) => [k + 1, v]), colour: token('--C') }] });
  box.append(el('div', { class: 'block' }, el('div', { class: 'eyebrow' }, 'the simulated sequencer, measured'), el('div', { class: 'fig-row' }, a, b)));
}

// ---- threads ----------------------------------------------------------------------------------------
function scalingCard() {
  view.scaleBtn = btn('Measure', { icon: 'play', kind: 'primary', onClick: runScaling });
  view.scaleBody = el('div', {}, empty({ icon: 'activity', title: 'Not measured yet', text: 'One computation on one worker with 1, 2, 4 and all its threads: its speed-up, and the trace of its task graph.' }));
  const strip = el('div', { class: 'ctl-strip' }, el('span', { class: 'small faint strip-note' }, '60,000 reads of a simulated genome (200 kb), on one worker'), el('span', { class: 'grow' }), view.scaleBtn);
  return card({ title: 'Threads', sub: 'one computation on 1, 2, 4 … threads, and its speed-up', body: [strip, view.scaleBody] });
}

async function runScaling() {
  const stop = working(view.scaleBtn, view.scaleBody, 'Running on 1, 2, 4 … threads…');
  try {
    kept.scaling = await api.run('run.scaling', { genome: 'sim:repeats-200k', reads: 60000, seed: 3 });
    stop();
    drawScaling();
  } catch (e) { stop(e); }
}

function drawScaling() {
  const r = kept.scaling, box = clear(view.scaleBody);
  const runs = r.runs;
  const maxT = Math.max(...runs.map((x) => x.threads));
  const c = el('div', { class: 'fig' });
  c.innerHTML = chart({ width: 440, height: 240, x: [1, maxT], y: [0, maxT], xlabel: 'threads', ylabel: 'speed-up over one thread',
    series: [{ points: [[1, 1], [maxT, maxT]], colour: token('--ink-faint'), dash: '4 4', label: 'linear' }, { points: runs.map((x) => [x.threads, x.speedup]), colour: token('--accent'), dots: true, label: 'measured' }] });
  const last = runs[runs.length - 1].trace;
  const tt = el('div', { class: 'fig' });
  tt.innerHTML = timeline({ width: Math.max(420, Math.min(640, width(view.scaleBody) - 480)), rows: [...Array(last.threads).keys()], seconds: last.seconds, rowHeight: 14, label: (k) => `thread ${k}`,
    bars: last.tasks.map(([task, thread, s, e]) => ({ row: thread, start: s, end: e, colour: CATEGORIES[task % 3], opacity: 0.8, title: `chunk ${task}` })) });
  const best = runs.reduce((a, x) => (x.speedup > a.speedup ? x : a), runs[0]);
  box.append(el('div', { class: 'pad' },
    metrics(runs.map((x) => [fmt.seconds(x.seconds), `${x.threads} thread${x.threads > 1 ? 's' : ''}${x.threads > 1 ? `, ×${fmt.number(x.speedup, 2)}` : ''}`, x === best && x.threads > 1 ? 'good' : ''])),
    el('p', { class: 'note', style: { margin: '-4px 0 0' } }, `${fmt.count(r.reads)} reads: ${r.work}. Hyperthreads share a core's lanes, so the curve bends at the number of cores.`)),
    el('div', { class: 'block' }, el('div', { class: 'fig-row' }, c, el('div', {}, el('div', { class: 'eyebrow' }, `the task graph on ${last.threads} threads`), tt))));
}

// ---- devices ----------------------------------------------------------------------------------------
// The same Smith–Waterman scores of one query against many targets, on the GPU by WebGPU in this
// page (src/gpu/smithwaterman.js) and on the CPU by the engine; the scores must agree exactly.
function devicesCard() {
  view.devBtn = btn('Score on both', { icon: 'play', kind: 'primary', onClick: runDevices });
  view.devBody = el('div', {}, empty({ icon: 'gpu', title: 'Not scored yet', text: 'A GPU is a device of its own: here, the browser\'s, through WebGPU. The kernel is written in WGSL in the page: one invocation per target, Gotoh\'s local recurrences down it, as SWIPE and CUDASW++ schedule them. The engine scores the same pairs on the CPU three ways: Gotoh\'s scalar rows, Farrar\'s stripes along the query, and one target to each lane of a vector register (Rognes 2011), on every thread. Since the scores are integers, all of them must agree exactly.' }));
  const strip = el('div', { class: 'ctl-strip' },
    fld('targets', seg([1000, 4000, 16000].map((n) => [n, fmt.count(n)]), kept.targets, (v) => { kept.targets = v; })),
    el('span', { class: 'small faint strip-note' }, 'reads of 150 bases against a query of 120, Smith–Waterman with affine gaps'),
    el('span', { class: 'grow' }), view.devBtn);
  return card({ title: 'Devices: the GPU and the CPU', sub: 'one kernel of integers on every device and schedule, which must agree exactly', body: [strip, view.devBody] });
}

async function runDevices() {
  const stop = working(view.devBtn, view.devBody, 'Drawing targets, then scoring them on both devices…');
  try {
    const gpu = await import('../gpu/smithwaterman.js');
    const dev = await gpu.device();
    // The targets: reads of a simulated genome, the query one of them shortened, from the engine.
    const batch = await api.run('align.batch_input', { genome: 'sim:repeats-200k', targets: kept.targets, length: 150, query: 120, seed: 11 });
    const decode = (s) => Uint8Array.from(s, (c) => 'ACGTN'.indexOf(c));
    const query = decode(batch.query), targets = batch.targets.map(decode);
    const M = [];
    for (let a = 0; a < 5; a++) for (let b = 0; b < 5; b++) M.push(a === 4 || b === 4 ? -1 : a === b ? 2 : -4);
    // Three schedules of the same scores on the CPU: scalar rows, stripes along the query, a target to a lane.
    const kernels = [];
    for (const kernel of ['gotoh', 'striped', 'inter']) kernels.push(await api.run('align.batch', { query: batch.query, targets: batch.targets, gap_open: 4, gap_extend: 2, kernel }));
    const cpu = kernels[0];
    let g = null;
    if (dev) {
      await gpu.scores(query, targets.slice(0, 64), M, 4, 2);   // the first run compiles the pipeline
      g = await gpu.scores(query, targets, M, 4, 2);
    }
    const agrees = (scores) => scores.every((v, k) => v === cpu.scores[k]);
    const rows = kernels.map((x) => ({ device: 'the CPU engine', kernel: x.kernel, seconds: x.seconds, gcups: x.gcups, ok: agrees(x.scores) }));
    if (g) rows.push({ device: 'the GPU, in this page', kernel: g.adapter, seconds: g.seconds, gcups: g.gcups, ok: agrees(g.scores) });
    kept.devices = { rows, gpu: !!g, same: g ? agrees(g.scores) : null, targets: targets.length, query: query.length, cells: cpu.cells };
    stop();
    drawDevices();
  } catch (e) { stop(e); }
}

function drawDevices() {
  const d = kept.devices, box = clear(view.devBody);
  const fastest = d.rows.reduce((a, x) => (x.gcups > a.gcups ? x : a), d.rows[0]);
  const table = el('table', { class: 'data wide' }, el('tr', {}, ...[['device'], ['kernel'], ['time', 1], ['GCUPS', 1], ['scores']].map(([h, n]) => el('th', { class: n ? 'num' : '' }, h))),
    ...d.rows.map((x, i) => el('tr', {}, el('td', {}, x.device), el('td', { class: 'faint' }, x.kernel), el('td', { class: 'num' }, fmt.seconds(x.seconds)),
      el('td', { class: 'num' }, x === fastest ? el('b', { style: { fontWeight: 600 } }, fmt.number(x.gcups, 3)) : fmt.number(x.gcups, 3)),
      el('td', {}, i === 0 ? el('span', { class: 'faint small' }, 'the reference') : el('span', { class: `badge ${x.ok ? 'good' : 'bad'}` }, x.ok ? 'equal' : 'differ')))));
  const allSame = d.rows.every((x) => x.ok);
  box.append(el('div', { class: 'pad' },
    metrics([[fmt.count(d.targets), 'targets'], [fmt.count(d.query), 'bases of query'], [fmt.count(d.cells), 'cells'], [fmt.number(fastest.gcups, 3), `GCUPS at best, ${fastest.device === 'the CPU engine' ? 'on the CPU' : 'on the GPU'}`, 'good']])),
    el('div', { class: 'table-wrap' }, table),
    el('div', { class: 'pad last' }, d.gpu ? callout(allSame ? 'Every score agrees: the devices and the schedules computed the same integers, as they must.' : 'The scores differ. An integer kernel must agree exactly; this is a fault to report.', allSame ? '' : 'bad', allSame ? 'check' : null)
      : callout('This browser offers no WebGPU adapter, so only the CPU ran. In Chrome or Edge on a machine with a GPU the kernel runs there.', 'warn')));
}
