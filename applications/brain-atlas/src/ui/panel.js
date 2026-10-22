// The panel's HTML: the point card, the choice of parcellation, the data
// mapping, the head and its electrodes, the legend, and the dialog about the
// data and method.

import { cssGradient, COLORMAP_NAMES } from '../render/colormaps.js';
import { tidy } from './data.js';
import { CITATIONS } from '../../../../domains/neuroscience/js/index.js';

export const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const rgb = (c) => (c ? `rgb(${c[0]}, ${c[1]}, ${c[2]})` : 'transparent');
const mm = (p) => p.map((x) => (Math.round(x * 10) / 10).toFixed(1)).join(', ');

// The card of the current point: the region there, its place in the tree, the
// coordinates in every frame that applies, and the selected region's measures.
export function pointCard({ atlas, point, voxel, region, path, selected, stats, talairach, ccf, value, scaled, pointsHere, hemisphere }) {
  const html = [];
  if (region) {
    const colour = region.colour || (region.parent && region.parent.colour);
    html.push(`<h2><span class="swatch" style="background:${rgb(colour)}"></span>${esc(region.name)}</h2>`);
    if (path.length > 1) html.push(`<p class="path">${path.slice(0, -1).map((r) => `<b>${esc(r.name)}</b>`).join('<span class="sep">›</span>')}</p>`);
  } else {
    html.push('<h2>Outside the parcellation</h2>');
    html.push(`<p class="sub">${atlas === 'human' ? 'No labelled region at this point: white matter, cerebrospinal fluid, or outside the brain.' : 'No structure at this point.'}</p>`);
  }
  html.push('<dl>');
  html.push(`<dt>${atlas === 'human' ? 'MNI 152' : 'RAS'} <small>mm</small></dt><dd>${mm(point)}</dd>`);
  if (scaled) html.push(`<dt>Subject <small>mm</small></dt><dd>${mm(scaled)}</dd>`);
  if (voxel) html.push(`<dt>Voxel</dt><dd>${voxel.map((x) => Math.round(x)).join(', ')}</dd>`);
  if (talairach) {
    html.push(`<dt>Talairach <small>Brett</small></dt><dd>${mm(talairach.brett)}</dd>`);
    html.push(`<dt>Talairach <small>Lancaster, SPM</small></dt><dd>${mm(talairach.spm)}</dd>`);
  }
  if (ccf) html.push(`<dt>CCF <small>µm, posterior · inferior · right</small></dt><dd>${ccf.map((x) => Math.round(x)).join(', ')}</dd>`);
  if (hemisphere) html.push(`<dt>Hemisphere</dt><dd>${esc(hemisphere)}</dd>`);
  if (value !== null && value !== undefined) html.push(`<dt>Mapped value</dt><dd>${tidy(value)}</dd>`);
  html.push('</dl>');
  if (pointsHere && pointsHere.length) html.push(`<p class="note"><b>Points near here:</b> ${pointsHere.map((p) => esc(p.name || mm(p.position)) + (p.value !== null ? ` (${tidy(p.value)})` : '')).join('; ')}</p>`);
  if (selected) {
    html.push(`<h3>Selected: ${esc(selected.name)}</h3>`);
    html.push('<dl>');
    if (stats) {
      html.push(`<dt>Volume</dt><dd>${stats.volume >= 1000 ? `${tidy(stats.volume / 1000)} cm³` : `${tidy(stats.volume)} mm³`} <small>${stats.voxels.toLocaleString()} voxels</small></dd>`);
      html.push(`<dt>Centroid <small>mm</small></dt><dd>${mm(stats.centroid)}</dd>`);
      html.push(`<dt>Extent <small>mm</small></dt><dd>${[0, 1, 2].map((a) => tidy(stats.bounds.max[a] - stats.bounds.min[a])).join(' × ')}</dd>`);
    } else {
      html.push('<dt>Volume</dt><dd><small>no voxels in this map</small></dd>');
    }
    html.push(`<dt>Leaves</dt><dd>${selected.leaves().length}</dd>`);
    html.push('</dl>');
    html.push('<div class="actions"><button type="button" data-action="centre">Go to its centroid</button><button type="button" data-action="deselect">Deselect</button></div>');
  } else {
    html.push('<p class="hint">Click a slice to move the point; scroll to step through slices; shift-drag to pan, ctrl-scroll to zoom. Drag the round view to turn it, click it to pick. Choose a region in the tree for its surface and measures.</p>');
  }
  return html.join('');
}

export function parcellationPanel(atlas, entries, current, overlay) {
  const html = ['<h3>Parcellation</h3>'];
  for (const p of entries) {
    html.push(`<label class="choice"><input type="radio" name="parcellation" value="${esc(p.id)}"${p.id === current ? ' checked' : ''}${p.reader === 'built-in' || !p.reader ? ' disabled' : ''}><span>${esc(p.name)}<span class="meta">${esc(p.modality)}${p.reader === 'built-in' ? ' · built in, no volume here' : ''}</span></span></label>`);
  }
  html.push(`<div class="slider"><span>Overlay</span><input type="range" id="overlay" min="0" max="1" step="0.05" value="${overlay}" aria-label="Overlay opacity"><span id="overlay-value">${Math.round(overlay * 100)}%</span></div>`);
  html.push('<div class="actions"><button type="button" data-action="toggle-brain" aria-pressed="true">Brain surface</button><button type="button" data-action="toggle-region" aria-pressed="true">Region surface</button></div>');
  return html.join('');
}

export function dataPanel(state, parcellationName) {
  const maps = COLORMAP_NAMES.map((n) => `<option value="${n}"${n === state.colormap ? ' selected' : ''}>${n}</option>`).join('');
  return `<h3>Data on the brain</h3>
<p class="note">Paste a table: a region and a value on each line (name, acronym, label or FreeSurfer key; separated by commas or tabs). Regions of <b>${esc(parcellationName)}</b> are coloured by value.</p>
<div class="data-form">
<textarea id="values-text" placeholder="region, value&#10;Area hOc1 (V1, 17, CalcS) - left hemisphere, 0.82&#10;left precentral, 1.5">${esc(state.valuesText)}</textarea>
<div class="row"><button type="button" class="pill" data-action="map-values">Map the values</button><label>colours <select id="colormap">${maps}</select></label><button type="button" class="pill" data-action="clear-values">Clear</button></div>
<div class="row">built in: <button type="button" class="pill" data-builtin="volume">volume of each region</button><button type="button" class="pill" data-builtin="hemisphere">hemisphere</button><button type="button" class="pill" data-builtin="depth">depth in the tree</button><button type="button" class="pill" data-builtin="group">group</button></div>
<div id="values-report" class="report">${state.valuesReport || ''}</div>
</div>
<p class="note" style="margin-top:12px">Or points: x, y, z on each line, then a value and a name if you like, in ${state.atlas === 'human' ? 'MNI or Talairach coordinates' : 'RAS millimetres or CCF micrometres'}.</p>
<div class="data-form">
<textarea id="points-text" placeholder="x, y, z, value, name&#10;-38, -22, 56, 3.1, left hand area&#10;0, -60, 20">${esc(state.pointsText)}</textarea>
<div class="row"><label>given in <select id="points-space">${state.atlas === 'human'
    ? ['mni:MNI 152', 'brett:Talairach (Brett)', 'spm:Talairach (Lancaster, SPM)', 'fsl:Talairach (Lancaster, FSL)', 'other:Talairach (Lancaster, other)'].map((o) => { const [v, t] = o.split(':'); return `<option value="${v}"${v === state.pointSpace ? ' selected' : ''}>${t}</option>`; }).join('')
    : ['ras:RAS mm', 'ccf:CCF µm'].map((o) => { const [v, t] = o.split(':'); return `<option value="${v}"${v === state.pointSpace ? ' selected' : ''}>${t}</option>`; }).join('')}</select></label>
<button type="button" class="pill" data-action="show-points">Show the points</button><button type="button" class="pill" data-action="clear-points">Clear</button></div>
<div id="points-report" class="report">${state.pointsReport || ''}</div>
</div>`;
}

export function headPanel(head) {
  if (!head.available) return '';
  const field = (key, label, value, was) => `<label class="field"><span>${label} <span class="was">template ${tidy(was)} mm</span></span><input type="number" step="1" min="60" max="400" id="head-${key}" value="${Math.round(value)}"></label>`;
  return `<h3>Head and electrodes</h3>
<p class="note">The scalp of the template, with the electrodes of the 10–20 system placed by the arcs between the nasion, inion and preauricular points. Set a subject's head size and the brain is scaled to it; paste channel values and they become a map on the scalp and the brain.</p>
<div class="actions"><button type="button" data-action="toggle-head" aria-pressed="${head.on}">Show the head</button>
<button type="button" data-system="10-20" aria-pressed="${head.system === '10-20'}">10–20</button><button type="button" data-system="10-10" aria-pressed="${head.system === '10-10'}">10–10</button></div>
<div style="margin-top:8px">
${field('length', 'Length, nasion to inion', head.measurements.length, head.template.length)}
${field('breadth', 'Breadth, ear to ear', head.measurements.breadth, head.template.breadth)}
${field('height', 'Height, ears to vertex', head.measurements.height, head.template.height)}
<div class="actions"><button type="button" data-action="apply-head">Scale the brain to this head</button><button type="button" data-action="reset-head">Template size</button></div>
<p class="note">Scales: ${head.scales.map((s) => s.toFixed(3)).join(' × ')}. Straight distances, not tape measures.</p>
</div>
<p class="note" style="margin-top:12px">Channel values, one per line: <b>Cz, 12.5</b>. Interpolated over the scalp by spherical splines (Perrin et al. 1989) and shown on the brain beneath.</p>
<div class="data-form">
<textarea id="channels-text" placeholder="Fz, 2.1&#10;Cz, 4.0&#10;Pz, 1.2&#10;C3, -1.0&#10;C4, 3.5">${esc(head.channelsText)}</textarea>
<div class="row"><button type="button" class="pill" data-action="map-channels">Map the channels</button><button type="button" class="pill" data-action="clear-channels">Clear</button></div>
<div id="channels-report" class="report">${head.channelsReport || ''}</div>
</div>
<p class="note">Landmarks on the template head are provisional: scalp points in the direction of the usual positions. The 10–20 positions follow Jasper (1958) and Oostenveld and Praamstra (2001); 65 positions of the 10–10 system above the landmark level.</p>`;
}

export function legend({ title, unit, colormap, range, keys, note }) {
  if (!title) return '';
  const html = [`<div><div class="title"><span>${esc(title)}</span><span class="unit">${esc(unit || '')}</span></div>`];
  if (colormap && range) {
    html.push(`<div class="bar"><span style="background:${cssGradient(colormap)}"></span></div>`);
    html.push(`<div class="ticks">${[0, 0.5, 1].map((t) => `<span style="left:${t * 100}%">${tidy(range[0] + t * (range[1] - range[0]))}</span>`).join('')}</div>`);
  }
  if (keys) html.push(`<div class="keys">${keys.map((k) => `<span><i style="background:${k.colour}"></i>${esc(k.name)}</span>`).join('')}</div>`);
  if (note) html.push(`<div class="note">${esc(note)}</div>`);
  html.push('</div>');
  return html.join('');
}

export function aboutHtml(atlasId, entries) {
  const cite = (k) => `<li>${esc(CITATIONS[k] || k)}</li>`;
  return `<p>The Brain Atlas shows a template brain in three slices and in the round, divides it by a parcellation, and maps data onto it. Every number is computed by the neuroscience domain of The Science Infrastructures from files the data library fetched; nothing is computed on a server.</p>
<h4>The template</h4>
<p>${atlasId === 'human' ? 'MNI 152 ICBM 2009c nonlinear asymmetric, the average of 152 adult brains and the reference space of the EBRAINS human brain atlas, at 1 mm (Fonov et al. 2011), from TemplateFlow. Coordinates are RAS millimetres: x to the right, y to the front, z to the top, with the origin near the anterior commissure. Talairach coordinates are approximations by Brett (1999) and Lancaster et al. (2007).' : 'The Allen Mouse Brain Common Coordinate Framework v3, the average of 1,675 mouse brains at 50 µm (Wang et al. 2020). Its own axes run posterior, inferior and right from the volume’s corner, in micrometres; the atlas shows it in RAS millimetres.'}</p>
<h4>The parcellations</h4>
<ul>${entries.map((p) => `<li><b>${esc(p.name)}</b>: ${esc(p.description)} <i>${esc(p.licence)}</i></li>`).join('')}</ul>
<h4>Method</h4>
<ul>
<li>Slices sample the template and the label volume on the graphics card by the point’s coordinates; the views follow the neurological convention, the subject’s left on the viewer’s left, the sagittal view from the left.</li>
<li>Surfaces are extracted from masks by surface nets (Gibson 1998) in a worker: the brain from the template’s mask, a region from its labels, the scalp from the head.</li>
<li>Values by region colour every voxel of the region’s leaves; a value given for a parent colours all beneath it.</li>
<li>Electrodes are placed by the arcs of the 10–20 system on the template’s scalp (Jasper 1958; Oostenveld and Praamstra 2001); channel values are interpolated by spherical splines (Perrin et al. 1989) and shown on the scalp and on the brain beneath.</li>
<li>A subject’s head size scales the template about the ears; the electrodes are placed again on the scaled scalp.</li>
</ul>
<h4>Sources</h4>
<ul>${['Fonov2011', 'Amunts2020', 'Schaefer2018', 'Yeo2011', 'Desikan2006', 'Wang2020', 'Brett2001', 'Lancaster2007', 'Jasper1958', 'Oostenveld2001', 'Perrin1989', 'Dickscheid2026'].map(cite).join('')}</ul>`;
}
