// Smith–Waterman scores of one query against many targets on the GPU, through WebGPU: the
// atlas's first accelerator. The kernel is the inter-sequence schedule of Rognes (2011, BMC
// Bioinformatics 12, 221, SWIPE) carried to a GPU (as CUDASW++ does: Liu, Schmidt and Maskell
// 2009): one invocation per target, running Gotoh's local recurrences down the target with a
// column of H and E for the query in its private memory, so that every invocation is
// independent and none waits for another.
//
//   H(i, j) = max(0, H(i−1, j−1) + s(q_i, t_j), E(i, j), F(i, j))
//   E(i, j) = max(H(i, j−1) − o − e, E(i, j−1) − e)        a gap in the query
//   F(i, j) = max(H(i−1, j) − o − e, F(i−1, j) − e)        a gap in the target
//
// with a gap of length L costing o + L·e, the convention of the genomics domain's align module.
// The scores are integers, so the GPU's and the CPU's must agree exactly: the page checks them
// against the engine's for the same pairs. Queries up to MAX_QUERY codes.

export const MAX_QUERY = 256;

const WGSL = /* wgsl */ `
struct Params { qlen: u32, ntargets: u32, gap_open: i32, gap_extend: i32 };
@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> query: array<u32>;
@group(0) @binding(2) var<storage, read> targets: array<u32>;
@group(0) @binding(3) var<storage, read> offsets: array<u32>;
@group(0) @binding(4) var<storage, read> matrix: array<i32>;
@group(0) @binding(5) var<storage, read_write> scores: array<i32>;

const NEG: i32 = -1000000;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let t = gid.x;
  if (t >= params.ntargets) { return; }
  var H: array<i32, ${MAX_QUERY}>;
  var E: array<i32, ${MAX_QUERY}>;
  for (var i = 0u; i < params.qlen; i = i + 1u) { H[i] = 0; E[i] = NEG; }
  let oe = params.gap_open + params.gap_extend;
  var best = 0;
  for (var j = offsets[t]; j < offsets[t + 1u]; j = j + 1u) {
    let b = targets[j];
    var diag = 0;
    var up = 0;
    var f = NEG;
    for (var i = 0u; i < params.qlen; i = i + 1u) {
      let e = max(H[i] - oe, E[i] - params.gap_extend);
      f = max(up - oe, f - params.gap_extend);
      let h = max(max(0, diag + matrix[query[i] * 5u + b]), max(e, f));
      diag = H[i];
      H[i] = h;
      E[i] = e;
      up = h;
      best = max(best, h);
    }
  }
  scores[t] = best;
}`;

let devicePromise = null;

// The GPU, or null where the browser has no WebGPU or no adapter.
export function device() {
  if (!devicePromise) {
    devicePromise = (async () => {
      if (typeof navigator === 'undefined' || !navigator.gpu) return null;
      const adapter = await navigator.gpu.requestAdapter({ powerPreference: 'high-performance' });
      if (!adapter) return null;
      const dev = await adapter.requestDevice();
      let info = {};
      try { info = adapter.info || (await adapter.requestAdapterInfo?.()) || {}; } catch { /* older browsers */ }
      const pipeline = dev.createComputePipeline({ layout: 'auto', compute: { module: dev.createShaderModule({ code: WGSL }), entryPoint: 'main' } });
      return { dev, pipeline, name: [info.vendor, info.architecture, info.description, info.device].filter(Boolean).join(' ') || 'a WebGPU adapter' };
    })().catch(() => null);
  }
  return devicePromise;
}

// Scores of the query against each target (codes 0–4), with the 5×5 matrix and the gaps.
export async function scores(query, targets, matrix, gapOpen, gapExtend) {
  const gpu = await device();
  if (!gpu) throw new Error('this browser has no WebGPU adapter');
  if (query.length > MAX_QUERY) throw new Error(`a query has at most ${MAX_QUERY} bases on the GPU`);
  const { dev, pipeline } = gpu;
  const n = targets.length;
  const offsets = new Uint32Array(n + 1);
  for (let k = 0; k < n; k++) offsets[k + 1] = offsets[k] + targets[k].length;
  const flat = new Uint32Array(Math.max(1, offsets[n]));
  for (let k = 0; k < n; k++) flat.set(targets[k], offsets[k]);
  const buffer = (data, usage) => {
    const b = dev.createBuffer({ size: Math.max(16, data.byteLength), usage, mappedAtCreation: true });
    new data.constructor(b.getMappedRange()).set(data);
    b.unmap();
    return b;
  };
  const S = GPUBufferUsage.STORAGE;
  const params = new Int32Array([query.length, n, gapOpen, gapExtend]);
  const bufs = [
    buffer(new Uint32Array(params.buffer), GPUBufferUsage.UNIFORM),
    buffer(Uint32Array.from(query), S), buffer(flat, S), buffer(offsets, S), buffer(Int32Array.from(matrix), S),
  ];
  const out = dev.createBuffer({ size: Math.max(16, 4 * n), usage: S | GPUBufferUsage.COPY_SRC });
  const read = dev.createBuffer({ size: Math.max(16, 4 * n), usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
  const group = dev.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [...bufs, out].map((b, k) => ({ binding: k, resource: { buffer: b } })) });
  const t0 = performance.now();
  const enc = dev.createCommandEncoder();
  const pass = enc.beginComputePass();
  pass.setPipeline(pipeline);
  pass.setBindGroup(0, group);
  pass.dispatchWorkgroups(Math.ceil(n / 64));
  pass.end();
  enc.copyBufferToBuffer(out, 0, read, 0, 4 * n);
  dev.queue.submit([enc.finish()]);
  await read.mapAsync(GPUMapMode.READ);
  const seconds = (performance.now() - t0) / 1000;
  const result = new Int32Array(read.getMappedRange().slice(0, 4 * n));
  read.unmap();
  for (const b of [...bufs, out, read]) b.destroy();
  const cells = query.length * offsets[n];
  return { scores: result, seconds, cells, gcups: cells / seconds / 1e9, adapter: gpu.name };
}

// The same recurrences in JavaScript, for a few pairs: the page's own check of the kernel.
export function scoreOne(query, target, matrix, gapOpen, gapExtend) {
  const m = query.length, NEG = -1e9, oe = gapOpen + gapExtend;
  const H = new Int32Array(m), E = new Int32Array(m).fill(NEG);
  let best = 0;
  for (let j = 0; j < target.length; j++) {
    let diag = 0, up = 0, f = NEG;
    for (let i = 0; i < m; i++) {
      const e = Math.max(H[i] - oe, E[i] - gapExtend);
      f = Math.max(up - oe, f - gapExtend);
      const h = Math.max(0, diag + matrix[query[i] * 5 + target[j]], e, f);
      diag = H[i];
      H[i] = h;
      E[i] = e;
      up = h;
      if (h > best) best = h;
    }
  }
  return best;
}
