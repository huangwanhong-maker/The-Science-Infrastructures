// The server's protocol (README §3): each call returns the engine's answer or throws with its message.

async function call(path, body) {
  const response = await fetch(path, body === undefined ? {} : {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body),
  });
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

export const engineState = () => call('/api/engine');
export const dataState = () => call('/api/data');
export const map = (year) => call('/api/map', { year });
export const spacetime = (first, last, currencies) => call('/api/spacetime', { first, last, currencies });
export const production = (table, extra = {}) => call('/api/production', { table, ...extra });
export const circuits = (network) => call('/api/circuits', { network });
export const setup = () => call('/setup');
export const runSetup = (again = false) => call(`/setup/run${again ? '?again=1' : ''}`, {});
