// The card about a place (or about the Earth as a whole when no place is
// chosen): where it is, the Sun there, and every value the loaded data give
// for it, with the course of each value over time for a chosen place. Where
// several models give the same quantity, their courses are drawn together.

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

// lines: [{ times, values, cursor, colour }]; the first is drawn strongest.
function sparkline(lines) {
  const W = 280, H = 34, pad = 3;
  const all = lines.flatMap((l) => l.values.map((v, i) => [l.times[i], v])).filter(([, v]) => Number.isFinite(v));
  if (all.length < 2) return '';
  const t0 = Math.min(...lines.map((l) => l.times[0])), t1 = Math.max(...lines.map((l) => l.times[l.times.length - 1]));
  let lo = Math.min(...all.map((p) => p[1])), hi = Math.max(...all.map((p) => p[1]));
  if (hi - lo < 1e-9) { lo -= 1; hi += 1; }
  const x = (t) => pad + ((t - t0) / (t1 - t0)) * (W - 2 * pad);
  const y = (v) => H - pad - ((v - lo) / (hi - lo)) * (H - 2 * pad);
  let svg = '';
  lines.forEach((l, k) => {
    const ok = l.values.map((v, i) => [l.times[i], v]).filter(([, v]) => Number.isFinite(v));
    if (ok.length < 2) return;
    const d = ok.map(([t, v], i) => `${i ? 'L' : 'M'}${x(t).toFixed(1)},${y(v).toFixed(1)}`).join('');
    if (k === 0) {
      const area = `${d}L${x(ok[ok.length - 1][0]).toFixed(1)},${H}L${x(ok[0][0]).toFixed(1)},${H}Z`;
      svg += `<path d="${area}" fill="${l.colour}" opacity=".14"/>`;
    }
    svg += `<path d="${d}" fill="none" stroke="${l.colour}" stroke-width="${k === 0 ? 1.5 : 1.1}" opacity="${k === 0 ? 1 : 0.85}" vector-effect="non-scaling-stroke"/>`;
  });
  const c = lines[0].cursor;
  if (c && c.time >= t0 && c.time <= t1) {
    const cx = x(c.time);
    svg += `<line x1="${cx}" x2="${cx}" y1="0" y2="${H}" stroke="#f4dbe7" stroke-width="1" vector-effect="non-scaling-stroke" opacity=".6"/>`;
    if (Number.isFinite(c.value)) svg += `<circle cx="${cx}" cy="${y(c.value)}" r="2.6" fill="#f4dbe7"/>`;
  }
  return `<svg viewBox="0 0 ${W} ${H}" preserveAspectRatio="none" aria-hidden="true">${svg}</svg>`;
}

// info: { title, sub, pinned, rows: [[term, html]], sections: [{ heading, rows, figure }],
//         series: [{ label, now, lines }], models: [{ model, colour }], hint }
export function renderPlace(el, info) {
  const rows = (list) => `<dl>${list.map(([dt, dd, kind]) => `<dt${kind ? ` class="${kind}"` : ''}>${esc(dt)}</dt><dd${kind ? ` class="${kind}"` : ''}>${dd}</dd>`).join('')}</dl>`;
  let html = '';
  if (info.pinned) html += '<button type="button" class="unpin" data-action="unpin" aria-label="Forget this place">✕</button>';
  html += `<h2 class="${info.title.length > 22 ? 'long' : ''}">${esc(info.title)}</h2>`;
  if (info.sub) html += `<p class="sub">${esc(info.sub)}</p>`;
  if (info.rows && info.rows.length) html += rows(info.rows);
  for (const s of info.sections || []) {
    if (!s.rows.length && !s.figure) continue;
    html += `<h4>${esc(s.heading)}</h4>${rows(s.rows)}`;
    if (s.figure) html += `<figure>${s.figure.svg}${s.figure.caption ? `<figcaption>${esc(s.figure.caption)}</figcaption>` : ''}</figure>`;
  }
  if (info.series && info.series.length) {
    html += '<h4>Over time</h4>';
    if (info.models) html += `<div class="model-keys">${info.models.map((m) => `<span><i style="background:${m.colour}"></i>${esc(m.model)}</span>`).join('')}</div>`;
    html += '<div class="series">';
    for (const s of info.series) {
      const line = sparkline(s.lines);
      if (!line) continue;
      html += `<div class="row"><span>${esc(s.label)}</span><span class="now">${s.now}</span>${line}</div>`;
    }
    html += '</div>';
  }
  if (info.hint) html += `<p class="hint">${esc(info.hint)}</p>`;
  el.innerHTML = html;
}
