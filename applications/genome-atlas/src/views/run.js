// The Run view: where the computation happens. The machine and what the engine was compiled
// for; the workers the cluster keeps, here and (when configured) on other machines; a job cut
// into shards and pulled by free workers, drawn on a timeline with each shard's threads; and
// the same work on 1 to n threads. The recipe that builds the engine is here too.

import { el, clear, notice, field } from '../dom.js';
import * as api from '../api.js';
import { setupPanel } from '../setup.js';
import { chart, timeline } from '../draw/chart.js';
import * as fmt from '../format.js';
import { CATEGORIES, token } from '../palette.js';

let view = null;

export function mount(section, ctx) {
  view = { ctx, section, machine: el('div'), workers: el('div'), job: el('div'), scaling: el('div'), devices: null, setup: el('div') };
  section.append(
    el('h2', { class: 'view-title' }, 'Run'),
    el('p', { class: 'lede' }, 'The engine is the genomics and evolution domains in C++, compiled for this processor: its vector lanes as wide as the processor allows (calx’s CPU backend), its threads in calx’s scheduler. The server keeps several engines alive as workers, and can keep more on other machines; a job is cut into shards that any worker can compute from a few numbers, and merged in the order of the shards, so that its answer never depends on who computed what.'),
    el('h3', { class: 'section-title' }, 'The machine'), view.machine,
    el('h3', { class: 'section-title' }, 'The workers'), view.workers,
    el('h3', { class: 'section-title' }, 'A job across the workers'), view.job,
    el('h3', { class: 'section-title' }, 'Threads'), view.scaling,
    el('h3', { class: 'section-title' }, 'Devices: the GPU and the CPU'), view.devices = el('div'),
    el('h3', { class: 'section-title' }, 'The engine’s recipe'), view.setup,
  );
  buildJob();
  buildScaling();
  buildDevices();
  view.setup.append(setupPanel());
}

export async function show(ctx) {
  view.ctx = ctx;
  await refresh();
}

async function refresh() {
  let e;
  try { e = await api.engineState(); } catch (err) { clear(view.machine).append(el('p', { class: 'error' }, err.message)); return; }
  // Ask a worker for its target, which also starts the workers.
  let target = null;
  if (e.current) target = await api.run('target').catch(() => null);
  if (target) e = await api.engineState().catch(() => e);
  drawMachine(e, target);
  drawWorkers(e);
}

function chip(text, on, title) { return el('span', { class: `badge ${on ? 'good' : ''}`, title: title || null }, text); }

function drawMachine(e, t) {
  const m = clear(view.machine);
  if (!t) { m.append(el('p', { class: 'note' }, e.built ? 'The engine is built but stale: build it again below.' : 'The engine is not built yet: build it below.')); return; }
  const kv = el('dl', { class: 'kv' });
  const row = (k, ...v) => kv.append(el('dt', {}, k), el('dd', {}, ...v));
  row('processor', t.brand || t.architecture);
  row('threads', `${t.threads} in each worker’s pool (${t.hardware_threads} on the machine)`);
  row('instruction sets', el('span', { class: 'row' }, chip('SSE4.1', true), chip('AVX2', t.avx2), chip('AVX-512', t.avx512)));
  row('compiled for', `${t.built_for}: vector lanes of ${t.lane_bytes} bytes (${t.lane_bytes} lanes of 8 bits, ${t.lane_bytes / 2} of 16)`);
  if (e.built_at) row('built', `${e.built_at.replace('T', ' ').replace('Z', ' UTC')} with g++ (${e.toolchain})`);
  row('operations', `${t.operations.length}: ${t.operations.join(', ')}`);
  m.append(kv);
}

function drawWorkers(e) {
  const w = clear(view.workers);
  const workers = e.cluster?.workers || [];
  const table = el('table', { class: 'data' }, el('tr', {}, ...['worker', 'where', 'threads', 'lanes', 'requests', 'busy (s)', 'state'].map((h) => el('th', {}, h))));
  for (const x of workers) {
    table.append(el('tr', {},
      el('td', {}, el('span', { class: 'swatch', style: { background: CATEGORIES[workers.indexOf(x) % CATEGORIES.length] } }), ' ', x.name),
      el('td', {}, x.remote ? 'remote' : 'this computer'),
      el('td', { class: 'num' }, x.threads ?? '–'),
      el('td', {}, x.built_for ? `${x.built_for}, ${x.lane_bytes} B` : '–'),
      el('td', { class: 'num' }, fmt.count(x.served)),
      el('td', { class: 'num' }, fmt.number(x.seconds, 3)),
      el('td', {}, x.alive ? el('span', { class: 'badge good' }, 'alive') : el('span', { class: 'badge bad', title: x.error || '' }, 'stopped'))));
  }
  if (!workers.length) w.append(el('p', { class: 'note' }, 'No worker has started yet; they start with the first request.'));
  else w.append(table);
  w.append(el('div', { class: 'row', style: { marginTop: '10px' } },
    el('button', { class: 'button small', onclick: async () => { try { await api.restart(); notice('The workers were stopped; they start again with the next request.'); refresh(); } catch (err) { notice(err.message, true); } } }, 'Restart the workers'),
    el('span', { class: 'faint small' }, `${e.cluster?.local ?? 2} local; remote workers from ${e.cluster?.config || 'cluster.json'}`)));
  w.append(el('details', { class: 'small', style: { marginTop: '8px', maxWidth: '820px' } },
    el('summary', { class: 'faint' }, 'Adding a worker on another machine'),
    el('p', { class: 'note' }, 'Build the engine there from the same sources (the recipe), then write the command that starts it into the cluster file, and restart the workers. A remote worker takes the shards of jobs on simulated genomes; genomes read from this computer’s files stay with the local workers.'),
    el('pre', { class: 'mono', style: { background: 'var(--paper-3)', padding: '8px 10px', borderRadius: '6px' } },
      '{\n  "workers": [\n    {"name": "eu_region_server", "command": ["ssh", "eu_region_server",\n      "~/The-Science-Infrastructures/data/cache/engines/genome/genome --serve"]}\n  ]\n}')));
}

// ---- a job ----------------------------------------------------------------------------------------
function buildJob() {
  const box = view.job;
  const reads = el('select', {});
  for (const n of [20000, 100000, 300000, 1000000]) reads.append(el('option', { value: n, selected: n === 100000 }, fmt.count(n)));
  const shards = el('select', {});
  for (const n of [1, 4, 8, 16, 32, 64]) shards.append(el('option', { value: n, selected: n === 16 }, String(n)));
  const genomeSel = el('select', {},
    el('option', { value: 'sim:repeats-200k' }, 'a simulated genome (200 kb)'),
    el('option', { value: 'sim:bacterium-5m' }, 'a simulated bacterium (5 Mb)'),
    el('option', { value: 'current' }, 'the genome chosen above'));
  const result = el('div');
  const go = el('button', { class: 'button primary', onclick: async () => {
    go.disabled = true;
    clear(result).append(el('p', { class: 'wait' }, 'Simulating reads across the workers…'));
    const genome = genomeSel.value === 'current' ? view.ctx.genome : genomeSel.value;
    const t0 = performance.now();
    try {
      const r = await api.job('reads.sample', { genome, reads: Number(reads.value), shards: Number(shards.value), seed: 7, length: 150, error_first: 0.001, error_last: 0.012 });
      drawJob(result, r, (performance.now() - t0) / 1000);
      refresh();
    } catch (e) { clear(result).append(el('p', { class: 'error' }, e.message)); }
    go.disabled = false;
  } }, 'Run');
  box.append(
    el('p', { class: 'note' }, 'Short reads simulated from a genome (each read from its own random stream, so a shard is just which reads), with their qualities and errors measured; each shard is a task graph of chunks on a worker’s threads.'),
    el('div', { class: 'row' }, field('reads', reads), field('shards', shards), field('genome', genomeSel), go),
    result);
}

function drawJob(box, r, wall) {
  clear(box);
  const job = r._job;
  const workers = [...new Set(job.trace.map((t) => t.worker).concat(job.merge.worker ? [job.merge.worker] : []))].sort();
  const colour = (name) => CATEGORIES[workers.indexOf(name) % CATEGORIES.length];
  const engineSeconds = job.trace.reduce((s, t) => s + (t.engine_seconds || 0), 0);
  box.append(el('p', { class: 'note' },
    `${fmt.count(r.reads)} reads (${fmt.count(r.bases)} bases) in ${job.shards} shards on ${workers.length} worker${workers.length > 1 ? 's' : ''}: ${fmt.seconds(job.seconds)} from plan to merge, ${fmt.seconds(engineSeconds)} of engine time in the shards. ${fmt.number(r.error_per_base * 100, 3)}% of bases are errors; ${fmt.count(r.minimizers)} minimizers sketch them.`));
  // The timeline: one row per worker, the shards as bars, the merge last.
  const rows = workers;
  const bars = job.trace.map((t) => ({ row: rows.indexOf(t.worker), start: t.start, end: t.end, colour: colour(t.worker), opacity: 0.55 + 0.4 * ((t.shard % 2) === 0), title: `shard ${t.shard}: ${fmt.count(t.size)} reads, ${fmt.seconds(t.engine_seconds)} in the engine` }));
  if (job.merge.worker) bars.push({ row: rows.indexOf(job.merge.worker), start: job.merge.start, end: Math.max(job.merge.end, job.merge.start + job.seconds * 0.004), colour: token('--ink'), title: 'the merge' });
  const width = Math.min(980, view.section.clientWidth - 70);
  const t = el('div'); t.innerHTML = timeline({ width, rows, bars, seconds: job.seconds });
  box.append(el('p', { class: 'eyebrow' }, 'shards on the workers, in seconds from the plan'), t);
  // One shard's threads: the busiest shard's task graph.
  const big = job.trace.filter((x) => x.threads?.tasks?.length).sort((a, b) => b.threads.tasks.length - a.threads.tasks.length)[0];
  if (big) {
    const th = big.threads;
    const trows = [...Array(th.threads).keys()];
    const tbars = th.tasks.map(([task, thread, s, e]) => ({ row: thread, start: s, end: e, colour: colour(big.worker), opacity: 0.5 + 0.45 * (task % 2), title: `chunk ${task}` }));
    const tt = el('div'); tt.innerHTML = timeline({ width, rows: trows, bars: tbars, seconds: th.seconds, rowHeight: 14, label: (k) => `thread ${k}` });
    box.append(el('p', { class: 'eyebrow' }, `inside shard ${big.shard} on ${big.worker}: its chunks on ${th.threads} threads, ${fmt.percent(th.utilisation, 0)} busy`), tt);
  }
  // What was measured: the error rate and mean quality along the read.
  const n = r.error_rate.length;
  const q = el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '20px' } });
  const a = el('div'); a.innerHTML = chart({ width: 440, height: 220, x: [1, n], y: [0, Math.max(...r.error_rate) * 1.15 || 0.01], xlabel: 'position in the read', ylabel: 'errors per base', series: [{ points: r.error_rate.map((v, k) => [k + 1, v]), colour: token('--accent'), area: true }] });
  const b = el('div'); b.innerHTML = chart({ width: 440, height: 220, x: [1, n], y: [Math.min(...r.mean_quality) - 2, Math.max(...r.mean_quality) + 2], xlabel: 'position in the read', ylabel: 'mean Phred quality', series: [{ points: r.mean_quality.map((v, k) => [k + 1, v]), colour: token('--C') }] });
  q.append(a, b);
  box.append(el('p', { class: 'eyebrow' }, 'the simulated sequencer, measured'), q);
  void wall;
}

// ---- threads ----------------------------------------------------------------------------------------
function buildScaling() {
  const box = view.scaling;
  const result = el('div');
  const go = el('button', { class: 'button', onclick: async () => {
    go.disabled = true;
    clear(result).append(el('p', { class: 'wait' }, 'Running on 1, 2, 4 … threads…'));
    try {
      const r = await api.run('run.scaling', { genome: 'sim:repeats-200k', reads: 60000, seed: 3 });
      const runs = r.runs;
      const maxT = Math.max(...runs.map((x) => x.threads));
      const c = el('div');
      c.innerHTML = chart({ width: 440, height: 240, x: [1, maxT], y: [0, maxT], xlabel: 'threads', ylabel: 'speed-up over one thread',
        series: [{ points: [[1, 1], [maxT, maxT]], colour: token('--ink-faint'), dash: '4 4', label: 'linear' }, { points: runs.map((x) => [x.threads, x.speedup]), colour: token('--accent'), dots: true, label: 'measured' }] });
      const last = runs[runs.length - 1].trace;
      const tt = el('div');
      tt.innerHTML = timeline({ width: Math.min(560, view.section.clientWidth - 520), rows: [...Array(last.threads).keys()], seconds: last.seconds, rowHeight: 14, label: (k) => `thread ${k}`,
        bars: last.tasks.map(([task, thread, s, e]) => ({ row: thread, start: s, end: e, colour: CATEGORIES[task % 3], opacity: 0.8, title: `chunk ${task}` })) });
      clear(result).append(el('p', { class: 'note' }, `${fmt.count(r.reads)} reads: ${r.work}. ${runs.map((x) => `${x.threads} thread${x.threads > 1 ? 's' : ''}: ${fmt.seconds(x.seconds)}`).join(' · ')}. Hyperthreads share a core's lanes, so the curve bends at the number of cores.`),
        el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '24px' } }, c, el('div', {}, el('p', { class: 'eyebrow' }, `the task graph on ${last.threads} threads`), tt)));
    } catch (e) { clear(result).append(el('p', { class: 'error' }, e.message)); }
    go.disabled = false;
  } }, 'Measure');
  box.append(el('p', { class: 'note' }, 'One computation on one worker with 1, 2, 4 and all its threads: its speed-up, and the trace of its task graph.'), go, result);
}

// ---- devices ----------------------------------------------------------------------------------------
// The same Smith–Waterman scores of one query against many targets, on the GPU by WebGPU in this
// page (src/gpu/smithwaterman.js) and on the CPU by the engine; the scores must agree exactly.
function buildDevices() {
  const box = view.devices;
  const result = el('div');
  const count = el('select', {});
  for (const n of [1000, 4000, 16000]) count.append(el('option', { value: n, selected: n === 4000 }, `${fmt.count(n)} targets`));
  const go = el('button', { class: 'button', onclick: async () => {
    go.disabled = true;
    clear(result).append(el('p', { class: 'wait' }, 'Drawing targets, then scoring them on both devices…'));
    try {
      const gpu = await import('../gpu/smithwaterman.js');
      const dev = await gpu.device();
      // The targets: reads of a simulated genome, the query one of them shortened, from the engine.
      const batch = await api.run('align.batch_input', { genome: 'sim:repeats-200k', targets: Number(count.value), length: 150, query: 120, seed: 11 });
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
      const same = g ? agrees(g.scores) : null;
      const rows = kernels.map((x) => ['the CPU engine', x.kernel, x.seconds, x.gcups, agrees(x.scores)]);
      if (g) rows.push(['the GPU, in this page', g.adapter, g.seconds, g.gcups, same]);
      const table = el('table', { class: 'data' }, el('tr', {}, ...['device', 'kernel', 'time', 'GCUPS', 'scores'].map((h) => el('th', {}, h))),
        ...rows.map(([d, k, s, c, ok], i) => el('tr', {}, el('td', {}, d), el('td', { class: 'faint' }, k), el('td', { class: 'num' }, fmt.seconds(s)), el('td', { class: 'num' }, fmt.number(c, 3)),
          el('td', {}, i === 0 ? el('span', { class: 'faint small' }, 'the reference') : el('span', { class: `badge ${ok ? 'good' : 'bad'}` }, ok ? 'equal' : 'differ')))));
      clear(result).append(
        el('p', { class: 'note' }, `${fmt.count(targets.length)} reads of 150 bases against a query of ${query.length}: ${fmt.count(cpu.cells)} cells.`),
        table,
        g ? el('p', {}, el('span', { class: `badge ${same ? 'good' : 'bad'}` }, same ? 'every score agrees' : 'the scores differ'), ' ', el('span', { class: 'faint small' }, same ? 'The devices and the schedules computed the same integers, as they must.' : 'An integer kernel must agree exactly; this is a fault to report.'))
          : el('p', { class: 'note' }, 'This browser offers no WebGPU adapter, so only the CPU ran. In Chrome or Edge on a machine with a GPU the kernel runs there.'));
    } catch (e) { clear(result).append(el('p', { class: 'error' }, e.message)); }
    go.disabled = false;
  } }, 'Score on both');
  box.append(
    el('p', { class: 'note' }, 'A GPU is a device of its own: here, the browser\'s, through WebGPU. The kernel is written in WGSL in the page: one invocation per target, Gotoh\'s local recurrences down it, as SWIPE and CUDASW++ schedule them. The engine scores the same pairs on the CPU three ways: Gotoh\'s scalar rows, Farrar\'s stripes along the query, and one target to each lane of a vector register (Rognes 2011), on every thread. Since the scores are integers, all of them must agree exactly.'),
    el('div', { class: 'row' }, count, go), result);
}
