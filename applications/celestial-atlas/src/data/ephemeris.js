// Fetching of the planetary ephemeris (JPL DE440s) prepared by data/, with
// progress. The astronomy domain reads the kernel (SPK); fetching it is the
// application's part, since the domains touch no network.

import { SPK } from '../../../../domains/astronomy/js/index.js';

export async function loadEphemeris(url, onProgress) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`${url}: HTTP ${response.status}`);
  const total = Number(response.headers.get('Content-Length')) || 0;
  if (!onProgress || !response.body || !total) {
    return new SPK(await response.arrayBuffer());
  }
  const buffer = new Uint8Array(total);
  const reader = response.body.getReader();
  let done = 0;
  for (;;) {
    const { value, done: finished } = await reader.read();
    if (finished) break;
    buffer.set(value, done);
    done += value.length;
    onProgress(done / total);
  }
  return new SPK(buffer.buffer);
}
