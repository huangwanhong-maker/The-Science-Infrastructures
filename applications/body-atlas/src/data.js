// The body's data: the product's catalogue and packs (data/products/body_anatomy.py), the engine's catalogue
// of measures and trees (served by serve.py), and the anatomy domain's knowledge base.
//
// Frames. The product and the engine speak LPS millimetres (DICOM's patient frame: +x the subject's left, +y
// posterior, +z superior; the anatomy domain's frames.hpp). The scene is three.js's, y up, right-handed:
// world = (x_L, z_S − zc, −y_P), so that +y is superior, +x the subject's left and +z anterior (L × S = A);
// zc centres the body's height. Seen from +z, the default view, the subject faces the viewer, the subject's
// left on the viewer's right. The two maps below are the only places the frames meet.

export let ZC = 0;
export function setCentre(zc) { ZC = zc; }
export const lpsToWorld = (p) => [p[0], p[2] - ZC, -p[1]];
export const worldToLps = (w) => [w[0], -w[2], w[1] + ZC];

async function json(url) {
  const r = await fetch(url);
  if (!r.ok) {
    let message = `${url}: ${r.status}`;
    try { message = (await r.json()).message || message; } catch (_) { /* not JSON */ }
    throw new Error(message);
  }
  return r.json();
}

export const loadStatus = () => json('/api/status');
export const loadProduct = () => json('/data/products/body-bp3d.json');
export const loadCatalogue = () => json('/api/catalogue');
export async function loadKnowledge() {
  try { return await json('/domains/anatomy/knowledge/structures.json'); } catch (_) { return { entries: {}, references: {} }; }
}

export async function ask(op, body) {
  const r = await fetch(`/api/${op}`, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) });
  const out = await r.json();
  if (!r.ok) throw new Error(out.message || `${op} failed`);
  return out;
}

async function gunzip(url, onBytes) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  const total = Number(r.headers.get('Content-Length')) || 0;
  const reader = r.body.getReader();
  const chunks = [];
  let got = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    got += value.length;
    if (onBytes) onBytes(got, total);
  }
  const stream = new Blob(chunks).stream().pipeThrough(new DecompressionStream('gzip'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

// Each element's mesh in the scene's frame: positions (Float32, world mm), normals (Int8, ×127), indices.
export function decodeElement(e, buffer) {
  const n = e.vertices;
  const q = new Uint16Array(buffer.buffer, buffer.byteOffset + e.positionOffset, n * 3);
  const nq = new Int8Array(buffer.buffer, buffer.byteOffset + e.normalOffset, n * 4);
  const [lx, ly, lz, hx, hy, hz] = e.box;
  const sx = (hx - lx || 1e-3) / 65535, sy = (hy - ly || 1e-3) / 65535, sz = (hz - lz || 1e-3) / 65535;
  const position = new Float32Array(n * 3);
  const normal = new Int8Array(n * 3);
  for (let v = 0; v < n; v++) {
    const x = lx + q[3 * v] * sx, y = ly + q[3 * v + 1] * sy, z = lz + q[3 * v + 2] * sz;
    position[3 * v] = x;
    position[3 * v + 1] = z - ZC;
    position[3 * v + 2] = -y;
    normal[3 * v] = nq[4 * v];
    normal[3 * v + 1] = nq[4 * v + 2];
    normal[3 * v + 2] = -nq[4 * v + 1];
  }
  const count = e.triangles * 3;
  const index = e.indexBytes === 2
    ? new Uint16Array(buffer.buffer.slice(buffer.byteOffset + e.indexOffset, buffer.byteOffset + e.indexOffset + count * 2))
    : new Uint32Array(buffer.buffer.slice(buffer.byteOffset + e.indexOffset, buffer.byteOffset + e.indexOffset + count * 4));
  return { position, normal, index };
}

// The packs, fetched in parallel; onPack(k, bytes) as each arrives.
export async function loadPacks(product, onPack, onProgress) {
  const sizes = product.packs.map(() => [0, 0]);
  const report = () => {
    if (!onProgress) return;
    const got = sizes.reduce((s, x) => s + x[0], 0), total = sizes.reduce((s, x) => s + x[1], 0);
    onProgress(got, total);
  };
  await Promise.all(product.packs.map(async (p, k) => {
    const bytes = await gunzip(`/data/products/${p.file}`, (got, total) => { sizes[k] = [got, total]; report(); });
    onPack(k, bytes);
  }));
}
