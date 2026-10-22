// The anatomy card: what is known of a structure, and what the atlas measures of it.

export const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmt1 = (x) => (Math.round(x * 10) / 10).toFixed(1);

export function volumeText(mm3) {
  if (!Number.isFinite(mm3)) return '—';
  if (mm3 >= 1000) return `${(mm3 / 1000).toFixed(mm3 >= 100000 ? 0 : 1)} <small>cm³</small>`;
  return `${Math.round(mm3)} <small>mm³</small>`;
}

const SECTIONS = [
  ['location', 'Where it lies'],
  ['function', 'What it does'],
  ['connections', 'Connections'],
  ['supply', 'Blood supply'],
  ['clinical', 'In the clinic'],
  ['development', 'Development'],
];

function block(value) {
  if (!value) return '';
  if (Array.isArray(value)) return `<ul>${value.map((v) => `<li>${v}</li>`).join('')}</ul>`;
  return `<p>${value}</p>`;
}

// model: {kind: 'structure' | 'system' | 'region' | 'welcome', ...}
export function cardHtml(model) {
  if (model.kind === 'welcome') return welcome(model);
  if (model.kind === 'system') return systemCard(model);
  const k = model.knowledge || {};
  const inherited = model.inherited || {};
  const html = [];
  html.push('<div class="card-head">');
  if (model.path && model.path.length) {
    html.push(`<div class="over">${model.path.map((p, i) => `${i ? '<span class="sep">›</span>' : ''}<button type="button" data-select="${esc(p.id)}">${esc(p.name)}</button>`).join('')}</div>`);
  }
  html.push(`<h2><span class="sw" style="background:${esc(model.colour || 'transparent')}"></span><span>${esc(model.name)}</span></h2>`);
  if (k.latin) html.push(`<div class="latin">${esc(k.latin)}</div>`);
  if (model.sideText) html.push(`<div class="side">${esc(model.sideText)}</div>`);
  html.push('<button type="button" class="close" data-action="deselect" aria-label="Close"><svg><use href="#i-close"/></svg></button>');
  html.push('<div class="card-actions">');
  html.push('<button type="button" data-action="focus"><svg><use href="#i-focus"/></svg>Focus</button>');
  html.push(`<button type="button" data-action="isolate" aria-pressed="${Boolean(model.isolated)}"><svg><use href="#i-isolate"/></svg>Isolate</button>`);
  if (model.hasSections) html.push('<button type="button" data-action="sections"><svg><use href="#i-slices"/></svg>In sections</button>');
  if (model.canHide) html.push('<button type="button" data-action="hide"><svg><use href="#i-hide"/></svg>Hide</button>');
  html.push('</div></div>');

  html.push('<div class="card-body">');
  if (k.summary) html.push(`<p class="summary">${k.summary}</p>`);
  else if (model.fallback) html.push(`<p class="summary">${model.fallback}</p>`);
  if (model.facts && model.facts.length) {
    html.push(`<div class="facts">${model.facts.map(([key, v]) => `<div><div class="k">${esc(key)}</div><div class="v">${v}</div></div>`).join('')}</div>`);
  }
  if (k.also && k.also.length) html.push(`<p class="note">Also called ${k.also.map((a) => `<b>${esc(a)}</b>`).join(', ')}.</p>`);
  for (const [field, title] of SECTIONS) {
    const own = k[field];
    const from = !own && inherited[field];
    if (!own && !from) continue;
    html.push(`<div class="sect${field === 'clinical' ? ' clinical' : ''}"><h4>${title}</h4>`);
    if (from) html.push(`<div class="inherited">Of the ${esc(from.name.toLowerCase())} as a whole:</div>`);
    html.push(block(own || from.value));
    html.push('</div>');
  }
  if (k.note) html.push(`<div class="sect"><h4>Remark</h4><p>${k.note}</p></div>`);
  if (model.near && model.near.length) {
    html.push(`<div class="sect"><h4>Runs beneath</h4><div class="chips">${model.near.map(chip).join('')}</div></div>`);
  }
  if (model.contacts && model.contacts.length) {
    html.push(`<div class="sect"><h4>Borders, in this atlas</h4><div class="chips">${model.contacts.map(chip).join('')}</div></div>`);
  }
  if (model.children && model.children.length) {
    html.push(`<div class="sect"><h4>Contains</h4><div class="chips">${model.children.map(chip).join('')}</div></div>`);
  }
  if (model.systems && model.systems.length) {
    html.push(`<div class="sect"><h4>Part of</h4><div class="chips">${model.systems.map((s) => `<button type="button" data-system="${esc(s.id)}">${esc(s.name)}</button>`).join('')}</div></div>`);
  }
  if (model.hasSections) {
    html.push('<div class="sect"><h4>In the sections</h4><div class="mini-sections"><div><canvas data-mini="axial"></canvas><span>axial</span></div><div><canvas data-mini="coronal"></canvas><span>coronal</span></div><div><canvas data-mini="sagittal"></canvas><span>sagittal</span></div></div></div>');
  }
  const refs = [...(model.sources || []), ...(k.refs || [])];
  if (refs.length) html.push(`<div class="card-refs">${refs.map((r) => `<p>${r}</p>`).join('')}</div>`);
  html.push('</div>');
  return html.join('');
}

function chip(c) {
  return `<button type="button" data-select="${esc(c.id)}"${c.side ? ` data-side="${esc(c.side)}"` : ''}><i style="background:${esc(c.colour || 'transparent')}"></i>${esc(c.name)}${c.note ? ` <small>${esc(c.note)}</small>` : ''}</button>`;
}

function systemCard(model) {
  return `<div class="card-head"><div class="over">A functional system</div><h2><span>${esc(model.name)}</span></h2>
<button type="button" class="close" data-action="deselect" aria-label="Close"><svg><use href="#i-close"/></svg></button></div>
<div class="card-body"><p class="summary">${model.summary}</p>
<div class="sect"><h4>Its parts in this atlas</h4><div class="chips">${model.members.map(chip).join('')}</div></div>
<p class="note">The other structures are drawn as glass while a system is shown. Choose a part to read about it.</p></div>`;
}

function welcome(model) {
  return `<div class="card-empty">
<h2>The human brain</h2>
<p>${esc(model.count)} structures of one brain, the average of 152 adults in the MNI 152 ICBM 2009c space, drawn from the segmentations of FreeSurfer, MASSP 2.0, the Harvard–Oxford and Julich-Brain atlases, and the fibre bundles of the EBRAINS atlas.</p>
<p>Click any part to read its anatomy: where it lies, what it does, what it connects, what supplies it and what follows when it fails. Double-click to fly to it.</p>
<p class="hint">Drag to turn the brain, scroll to come closer, right-drag to move it. Try the systems on the left, a cut through the brain below, or a guided tour above.</p>
</div>`;
}

export { fmt1 };
