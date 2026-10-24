// The server's protocol (README §3): each call returns the engine's or the server's answer, or throws
// with its message.

async function call(path, body, method) {
  const init = body === undefined ? { method: method || 'GET' } : {
    method: method || 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body),
  };
  const response = await fetch(path, init);
  const text = await response.text();
  let data;
  try {
    data = JSON.parse(text);
  } catch {
    throw new Error(`the server answered ${response.status} without JSON`);
  }
  if (!response.ok) throw new Error(data.message || `the server answered ${response.status}`);
  return data;
}

const withDataset = (dataset) => (dataset ? { dataset } : {});

export const engineState = () => call('/api/engine');
export const dataState = () => call('/api/data');
export const map = (year, dataset) => call('/api/map', { ...(year !== undefined && year !== null ? { year } : {}), ...withDataset(dataset) });
export const timeline = (dataset) => call('/api/timeline', withDataset(dataset));
export const spacetime = (first, last, dataset, currencies) => call('/api/spacetime', { first, last, ...withDataset(dataset), ...(currencies ? { currencies } : {}) });
export const production = (table, extra = {}) => call('/api/production', { table, ...extra });
export const circuits = (network) => call('/api/circuits', { network });
export const setup = () => call('/setup');
export const runSetup = (again = false) => call(`/setup/run${again ? '?again=1' : ''}`, {});

export const datasets = () => call('/api/datasets');
export const addDataset = (kind, title, files, base) => call('/api/datasets', { kind, title, files, ...(base ? { base } : {}) });
export const removeDataset = (id) => call(`/api/datasets/${encodeURIComponent(id)}`, undefined, 'DELETE');

export const plugins = () => call('/api/plugins');
export const describePlugin = (file) => call(`/api/plugins/describe?file=${encodeURIComponent(file)}`);
export const runPlugin = (file, params, dataset, year) => call('/api/plugins/run', { file, params, ...withDataset(dataset), ...(Number.isInteger(year) ? { year } : {}) });
export const job = (id) => call(`/api/plugins/job/${encodeURIComponent(id)}`);
export const cancelJob = (id) => call(`/api/plugins/job/${encodeURIComponent(id)}/cancel`, {});
