// A pool of terrain workers (worker.js): a loader for TileTree that builds meshes off the page's
// thread, with requests that can be abandoned.

export class TerrainWorkers {
  constructor(count = Math.max(2, Math.min(6, (globalThis.navigator?.hardwareConcurrency || 4) - 1))) {
    this.workers = Array.from({ length: count }, () => new Worker(new URL('./worker.js', import.meta.url), { type: 'module' }));
    this.next = 0;
    this.seq = 0;
    this.waiting = new Map();
    for (const w of this.workers) {
      w.onmessage = (event) => {
        const { id, mesh, error, aborted } = event.data;
        const job = this.waiting.get(id);
        if (!job) return;
        this.waiting.delete(id);
        if (mesh) job.resolve(mesh);
        else if (aborted) job.reject(new DOMException('aborted', 'AbortError'));
        else job.reject(new Error(error));
      };
    }
  }

  // A mesh for tile (z, x, y) of a figure ('earth', 'moon' or a figure object), from `url`, or flat
  // when url is null; of the tile's visible surface or its floor.
  build({ url, figure, z, x, y, exaggeration = 1, surface = 'visible' }, signal) {
    const id = ++this.seq;
    const worker = this.workers[this.next++ % this.workers.length];
    return new Promise((resolve, reject) => {
      this.waiting.set(id, { resolve, reject });
      if (signal) {
        if (signal.aborted) { this.waiting.delete(id); reject(new DOMException('aborted', 'AbortError')); return; }
        signal.addEventListener('abort', () => worker.postMessage({ cancel: id }), { once: true });
      }
      worker.postMessage({ id, url, figure, z, x, y, exaggeration, surface });
    });
  }

  terminate() {
    for (const w of this.workers) w.terminate();
  }
}
