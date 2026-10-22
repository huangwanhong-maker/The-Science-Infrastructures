// The scorecard: how well each model forecast what happened, from a product of
// kind `verification` (data/products/verification.py). One small chart per
// quantity, measure and region: the score against the lead time, a line per
// model in the model's colour, persistence (the forecast that nothing changes)
// as a dashed grey reference beneath them, and each point as solid as its
// number of cases warrants.

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

const METRICS = {
  rmse: 'root-mean-square error', acc: 'anomaly correlation', bias: 'mean error (bias)',
  track_km: 'track error', pressure_hpa: 'central-pressure error',
};
// For these, larger is better; for the others, nearer zero.
const HIGHER_IS_BETTER = new Set(['acc']);
const PERSISTENCE = '#9a8a96';

function nice(lo, hi, n = 4) {
  const raw = (hi - lo) / n || 1, mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw);
  const out = [];
  for (let v = Math.floor(lo / step) * step; v <= hi + step * 1e-6; v += step) out.push(+v.toPrecision(6));
  return out;
}

function chart(group, colourOf) {
  const W = 300, H = 150, L = 40, R = 8, T = 10, B = 26;
  const metric = group[0].metric;
  const leads = [...new Set(group.flatMap((s) => s.leads.map((l) => l.hours)))].sort((a, b) => a - b);
  const values = group.flatMap((s) => s.leads.map((l) => l.value)).filter(Number.isFinite);
  if (!leads.length || !values.length) return '';
  let lo = Math.min(0, ...values), hi = Math.max(...values);
  if (metric === 'acc') { lo = Math.min(0, lo); hi = 1; }
  if (metric === 'bias') { const m = Math.max(Math.abs(lo), Math.abs(hi)) || 1; lo = -m; hi = m; }
  const yTicks = nice(lo, hi);
  lo = Math.min(lo, yTicks[0]); hi = Math.max(hi, yTicks[yTicks.length - 1]);
  const x0 = leads[0], x1 = leads[leads.length - 1] || 1;
  const X = (h) => L + ((h - x0) / (x1 - x0 || 1)) * (W - L - R), Y = (v) => T + (1 - (v - lo) / (hi - lo || 1)) * (H - T - B);
  const step = x1 - x0 > 96 ? 48 : x1 - x0 > 24 ? 24 : 6;
  const xTicks = leads.filter((h) => h % step === 0);
  const maxCount = Math.max(1, ...group.flatMap((s) => s.leads.map((l) => l.count || 1)));
  // Reference lines: zero for a bias, 0.6 for the anomaly correlation (the usual limit of useful skill).
  const refs = [metric === 'bias' ? 0 : null, metric === 'acc' ? 0.6 : null].filter((v) => v !== null && v >= lo && v <= hi)
    .map((v) => `<line x1="${L}" x2="${W - R}" y1="${Y(v)}" y2="${Y(v)}" class="ref"/>`).join('');
  const grid = yTicks.map((v) => `<line x1="${L}" x2="${W - R}" y1="${Y(v)}" y2="${Y(v)}" class="grid"/><text x="${L - 5}" y="${Y(v) + 3.5}" class="yt">${esc(v)}</text>`).join('')
    + xTicks.map((h) => `<text x="${X(h)}" y="${H - B + 14}" class="xt">${h}</text>`).join('');
  // Persistence first, beneath the models.
  const ordered = [...group].sort((a, b) => (b.model === 'Persistence') - (a.model === 'Persistence'));
  const lines = ordered.map((s) => {
    const pts = s.leads.filter((l) => Number.isFinite(l.value)).sort((a, b) => a.hours - b.hours);
    const d = pts.map((l, i) => `${i ? 'L' : 'M'}${X(l.hours).toFixed(1)},${Y(l.value).toFixed(1)}`).join('');
    const reference = s.model === 'Persistence', c = reference ? PERSISTENCE : colourOf(s.model);
    const dots = pts.map((l) => {
      const solid = 0.35 + 0.65 * Math.min(1, Math.log1p(l.count || 1) / Math.log1p(maxCount));
      return `<circle cx="${X(l.hours).toFixed(1)}" cy="${Y(l.value).toFixed(1)}" r="${reference ? 1.6 : 2.3}" fill="${c}" fill-opacity="${solid.toFixed(2)}"><title>${esc(s.model)}, +${l.hours} h: ${esc(+l.value.toPrecision(4))} ${esc(s.units === '1' ? '' : s.units || '')} (${l.count ?? '?'} cases)</title></circle>`;
    }).join('');
    return `<path d="${d}" stroke="${c}" class="line${reference ? ' reference' : ''}"/>${dots}`;
  }).join('');
  return `<svg viewBox="0 0 ${W} ${H}" class="score-chart" role="img">${grid}${refs}<text x="${W - R}" y="${H - 4}" class="xt end">lead (h)</text>${lines}</svg>`;
}

// verification: the product; colourOf(model) -> colour.
export function renderScores(el, verification, colourOf) {
  if (!verification || !verification.scores || !verification.scores.length) {
    el.innerHTML = '<p class="note">No scores yet. They are built by <code>data/products/verification.py</code> from forecasts kept by <code>data/products/archive.py</code>; the schedule (<code>data/acquisition/schedule.py</code>) runs both every six hours.</p>';
    return;
  }
  const groups = new Map();
  for (const s of verification.scores) {
    const key = `${s.quantity}|${s.metric}|${s.region || ''}`;
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(s);
  }
  const blocks = [...groups.values()].map((group) => {
    const s = group[0];
    const models = group.map((g) => `<span><i class="${g.model === 'Persistence' ? 'reference' : ''}" style="background:${g.model === 'Persistence' ? PERSISTENCE : colourOf(g.model)}"></i>${esc(g.model)}</span>`).join('');
    const better = HIGHER_IS_BETTER.has(s.metric) ? 'higher is better' : 'nearer zero is better';
    const against = String(s.against || '').replace(/\s*\(.*$/, '');
    const unit = s.units && s.units !== '1' ? `${s.units} · ` : '';
    return `<section class="score">
      <h4>${esc(s.title || s.quantity)}: ${esc(s.quantity === 'cyclones' && s.metric === 'bias' ? 'central-pressure bias' : METRICS[s.metric] || s.metric)}${s.region ? `, ${esc(s.region)}` : ''}</h4>
      <p class="note" title="${esc([s.definition, s.against].filter(Boolean).join(' — '))}">${esc(unit)}${better}${against ? ` · against ${esc(against)}` : ''}</p>
      ${chart(group, colourOf)}
      <div class="keys">${models}</div>
    </section>`;
  }).join('');
  const notes = (verification.notes || []).map((n) => `<li>${esc(n)}</li>`).join('');
  el.innerHTML = `<p class="note">Updated ${esc((verification.updated || '').slice(0, 16).replace('T', ' '))} UTC. A faint point rests on few cases; hover a point for its count.</p>
    <div class="scores">${blocks}</div>${notes ? `<details class="score-notes"><summary>Notes</summary><ul>${notes}</ul></details>` : ''}`;
}
