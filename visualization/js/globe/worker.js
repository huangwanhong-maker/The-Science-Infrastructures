// A worker that fetches terrain tiles, reads them and builds their meshes (terrain.js), so that
// the page's thread only uploads the result to the graphics card.
//
// Messages in: { id, url, figure, exaggeration, surface } to build a tile (url null for a flat one,
// with z, x, y; surface 'visible' or 'floor'), { cancel: id } to abandon one. Messages out: { id, mesh } with the mesh's arrays
// transferred, { id, error }, or { id, aborted: true }.

import { buildMesh, decodeTile, flatTile } from './terrain.js';
import { figure } from './figure.js';

const pending = new Map();

self.onmessage = async (event) => {
  const msg = event.data;
  if (msg.cancel !== undefined) {
    pending.get(msg.cancel)?.abort();
    pending.delete(msg.cancel);
    return;
  }
  const { id, url, z, x, y } = msg;
  const controller = new AbortController();
  pending.set(id, controller);
  try {
    const fig = figure(msg.figure);
    let tile;
    if (url) {
      const r = await fetch(url, { signal: controller.signal });
      if (!r.ok) throw new Error(`terrain ${z}/${x}/${y}: HTTP ${r.status}`);
      tile = decodeTile(await r.arrayBuffer());
    } else {
      tile = flatTile(z, x, y);
    }
    const mesh = buildMesh(fig, tile, { exaggeration: msg.exaggeration ?? 1, surface: msg.surface || 'visible' });
    mesh.header = tile.header;
    pending.delete(id);
    self.postMessage({ id, mesh }, [mesh.position.buffer, mesh.normal.buffer, mesh.up.buffer, mesh.uv.buffer, mesh.provenance.buffer, mesh.height.buffer, mesh.index.buffer]);
  } catch (err) {
    pending.delete(id);
    if (err && err.name === 'AbortError') self.postMessage({ id, aborted: true });
    else self.postMessage({ id, error: String(err && err.message || err) });
  }
};
