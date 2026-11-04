// The server's API: the plug-ins, one request to a plug-in's engine, and the recipe that builds them.

export async function plugins() {
  const r = await fetch('/api/plugins');
  if (!r.ok) throw new Error(`the plug-ins could not be listed (${r.status})`);
  return (await r.json()).plugins;
}

export function asker(plugin) {
  return async (op, body = {}) => {
    const r = await fetch(`/api/${plugin}/${op}`, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) });
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
