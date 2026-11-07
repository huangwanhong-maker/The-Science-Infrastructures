// The server's API: the plug-ins, one request to a plug-in's engine, and the recipe that builds them.

export async function plugins() {
  const r = await fetch('/api/plugins');
  if (!r.ok) throw new Error(`the plug-ins could not be listed (${r.status})`);
  return (await r.json()).plugins;
}

export function asker(plugin) {
  return async (op, body = {}) => {
    // High priority: the tiles of the city and the globe hold most of a host's six connections while the camera moves, and are
    // asked at low priority, so that what a person asked for takes the next connection to free.
    const r = await fetch(`/api/${plugin}/${op}`, { method: 'POST', priority: 'high', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) });
    const answer = await r.json().catch(() => ({ message: `the server answered ${r.status}` }));
    if (!r.ok) throw new Error(answer.message || `the server answered ${r.status}`);
    return answer;
  };
}

export async function setup() {
  const r = await fetch('/setup');
  return r.json();
}
export async function runSetup(again = false) {
  const r = await fetch(`/setup/run${again ? '?again=1' : ''}`, { method: 'POST' });
  return r.json();
}
