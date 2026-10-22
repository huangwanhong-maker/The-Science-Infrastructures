// The data on this computer. The application's server reports its recipe at
// /setup: each step (a dataset fetched or a product built), whether it is done,
// and a run in progress, if any. The panel shows this, and lets the user fetch
// and build what is missing, or any step again, from the page.

const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const when = (t) => (t ? `${new Date(t * 1000).toISOString().slice(0, 16).replace('T', ' ')} UTC` : '');

export function setupPanel(container) {
  let timer = null;
  let ranHere = false;          // a run was started from this page: new data need a reload

  async function refresh() {
    clearTimeout(timer);
    let view = null;
    try {
      const r = await fetch('setup', { cache: 'no-store' });
      if (r.ok) view = await r.json();
    } catch { /* no server: reported below */ }
    if (!view) {
      container.innerHTML = '<h3>Data on this computer</h3><p class="setup-note">This page is not served by its own server, so its data cannot be fetched from here; the README gives the commands.</p>';
      return;
    }
    render(view);
    if (view.running) timer = setTimeout(refresh, 2000);
  }

  function render(view) {
    const missing = view.steps.filter((s) => !s.done);
    let note = missing.length ? `${missing.length} of ${view.steps.length} steps are missing.` : `All ${view.steps.length} steps are done.`;
    if (view.running) note += ' Fetching and building…';
    else if (view.state === 'failed') note += ` The last run failed at: ${view.failed.join(', ')}.`;
    else if (view.state === 'done' && ranHere) note += ' Reload the page to see the new data.';
    const rows = view.steps.map((s) => `<tr class="${s.done ? 'done' : 'missing'}">
      <th>${esc(s.title)}</th>
      <td>${s.done ? when(s.updated) : 'missing'}</td>
      <td>${view.running ? '' : `<button type="button" class="link" data-step="${esc(s.id)}">${s.done ? 'again' : 'fetch'}</button>`}</td></tr>`).join('');
    container.innerHTML = `
      <h3>Data on this computer</h3>
      <p class="setup-note">${esc(note)}</p>
      ${missing.length && !view.running ? '<p><button type="button" class="link" data-run="missing">Fetch and build what is missing</button></p>' : ''}
      <table class="setup"><tbody>${rows}</tbody></table>
      ${view.log ? `<pre class="setup-log">${esc(view.log)}</pre>` : ''}`;
    const log = container.querySelector('.setup-log');
    if (log && view.running) log.scrollTop = log.scrollHeight;
  }

  container.addEventListener('click', async (e) => {
    const all = e.target.closest('[data-run]'), one = e.target.closest('[data-step]');
    if (!all && !one) return;
    const query = one ? `?steps=${encodeURIComponent(one.dataset.step)}&again=1` : '';
    try { await fetch(`setup/run${query}`, { method: 'POST' }); ranHere = true; } catch { /* reported by refresh */ }
    refresh();
  });

  return { refresh, stop: () => clearTimeout(timer) };
}
