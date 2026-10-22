// The orientation compass, the toolbar beneath the stage, and the search.

import { esc } from './card.js';

// ---- the compass ------------------------------------------------------------------------------

const AXES = [
  ['R', [1, 0, 0], '#e7a4b4'], ['L', [-1, 0, 0], '#e7a4b4'],
  ['A', [0, 1, 0], '#a5d6bb'], ['P', [0, -1, 0], '#a5d6bb'],
  ['S', [0, 0, 1], '#c8b4ea'], ['I', [0, 0, -1], '#c8b4ea'],
];
export const VIEW_UPS = { S: [0, 1, 0], I: [0, 1, 0] };

export function drawCompass(svg, basis) {
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const items = AXES.map(([letter, v, colour]) => ({ letter, v, colour, x: dot(v, basis.right) * 29, y: -dot(v, basis.up) * 29, z: dot(v, basis.back) }));
  items.sort((a, b) => a.z - b.z);
  const parts = ['<circle class="ring" r="44"/>'];
  for (const it of items) {
    const back = it.z < -0.05;
    parts.push(`<line x1="0" y1="0" x2="${it.x.toFixed(1)}" y2="${it.y.toFixed(1)}" stroke="${it.colour}" stroke-opacity="${back ? 0.3 : 0.85}"/>`);
  }
  for (const it of items) {
    const back = it.z < -0.05;
    parts.push(`<text class="${back ? 'back' : ''}" x="${(it.x * 1.32).toFixed(1)}" y="${(it.y * 1.32).toFixed(1)}" fill="${it.colour}" data-view="${it.letter}"><title>View from the ${{ R: 'right', L: 'left', A: 'front', P: 'back', S: 'top', I: 'bottom' }[it.letter]}</title>${it.letter}</text>`);
  }
  svg.innerHTML = parts.join('');
}
export const VIEW_DIRECTIONS = Object.fromEntries(AXES.map(([l, v]) => [l, v]));

// ---- the toolbar ------------------------------------------------------------------------------

const opt = (attr, value, label, pressed, title = '') => `<button type="button" class="opt" ${attr}="${esc(value)}" aria-pressed="${Boolean(pressed)}"${title ? ` title="${esc(title)}"` : ''}>${label}</button>`;

export function toolbarHtml(state) {
  const g = [];
  if (state.mode === 'sections') {
    g.push(`<div class="group"><span class="caption">Main section</span><div class="row">${['axial', 'coronal', 'sagittal'].map((p) => opt('data-main', p, p[0].toUpperCase() + p.slice(1), state.mainPlane === p)).join('')}</div></div>`);
    g.push(`<div class="group"><span class="caption">Colour</span><div class="row"><input type="range" min="0" max="1" step="0.05" value="${state.overlay}" data-range="overlay" aria-label="Overlay strength"><span class="value">${Math.round(state.overlay * 100)}%</span></div></div>`);
    return g.join('');
  }
  if (state.mode === 'atlases' || state.mode === 'data') {
    g.push(`<div class="group"><span class="caption">Show</span><div class="row">${opt('data-show', '3d', 'In the round', state.show !== 'sections')}${opt('data-show', 'sections', 'In sections', state.show === 'sections')}</div></div>`);
    if (state.show === 'sections') {
      g.push(`<div class="group"><span class="caption">Main section</span><div class="row">${['axial', 'coronal', 'sagittal'].map((p) => opt('data-main', p, p[0].toUpperCase() + p.slice(1), state.mainPlane === p)).join('')}</div></div>`);
      return g.join('');
    }
  }
  g.push(`<div class="group"><span class="caption">Cortex</span><div class="row">${opt('data-cortex', 'solid', 'Solid', state.cortexStyle === 'solid')}${opt('data-cortex', 'glass', 'Glass', state.cortexStyle === 'glass')}${opt('data-cortex', 'hidden', 'Hidden', state.cortexStyle === 'hidden')}</div></div>`);
  const cut = state.cut;
  g.push(`<div class="group"><span class="caption">Cut</span><div class="row">${opt('data-cut', 'off', 'Off', !cut)}${opt('data-cut', 'sagittal', 'Sagittal', cut && cut.axis === 'sagittal')}${opt('data-cut', 'coronal', 'Coronal', cut && cut.axis === 'coronal')}${opt('data-cut', 'axial', 'Axial', cut && cut.axis === 'axial')}</div>${cut ? `<div class="row"><input type="range" min="${cut.min}" max="${cut.max}" step="1" value="${cut.position}" data-range="cut" aria-label="Position of the cut"><span class="value">${Math.round(cut.position)} mm</span>${opt('data-flip', '1', '⇄', false, 'Keep the other side')}${opt('data-cap', '1', 'Face', cut.cap, 'Show the section on the cut face')}</div>` : ''}</div>`);
  if (state.mode !== 'anatomy') { g.push(`<div class="group"><span class="caption">Names</span><div class="row">${opt('data-toggle', 'labels', 'Labels', state.labels)}</div></div>`); return g.join(''); }
  g.push(`<div class="group"><span class="caption">Explode</span><div class="row"><input type="range" min="0" max="1" step="0.02" value="${state.explode}" data-range="explode" aria-label="Explode" style="width:92px"></div></div>`);
  g.push(`<div class="group"><span class="caption">Names</span><div class="row">${opt('data-toggle', 'labels', 'Labels', state.labels)}</div></div>`);
  return g.join('');
}

// ---- the search -------------------------------------------------------------------------------

export function searchResults(query, entries, limit = 14) {
  const q = query.trim().toLowerCase();
  if (!q) return [];
  const scored = [];
  for (const e of entries) {
    const names = [e.name, ...(e.also || [])].map((s) => s.toLowerCase());
    let score = 0;
    for (const n of names) {
      if (n === q) score = Math.max(score, 4);
      else if (n.startsWith(q)) score = Math.max(score, 3);
      else if (n.split(/[\s,()'’-]+/).some((w) => w.startsWith(q))) score = Math.max(score, 2);
      else if (n.includes(q)) score = Math.max(score, 1);
    }
    if (score) scored.push([score, e]);
  }
  scored.sort((a, b) => b[0] - a[0] || a[1].name.length - b[1].name.length);
  return scored.slice(0, limit).map((s) => s[1]);
}

export function resultsHtml(list) {
  if (!list.length) return '<div class="note" style="padding:6px 8px">Nothing by that name.</div>';
  return list.map((e, i) => `<button type="button" data-select="${esc(e.id)}" class="${i === 0 ? 'active' : ''}"><i style="background:${esc(e.colour || 'transparent')}"></i><span>${esc(e.name)}</span><span class="path">${esc(e.path || '')}</span></button>`).join('');
}
