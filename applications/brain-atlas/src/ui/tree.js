// The tree of regions: collapsible, searchable, with each region's colour and
// the number of leaves beneath it; a click selects a region.

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const rgb = (c) => (c ? `rgb(${c[0]}, ${c[1]}, ${c[2]})` : null);

export class RegionTree {
  constructor(container, { onSelect, colourOf = (r) => r.colour, maxRoots = 400 }) {
    this.el = container;
    this.onSelect = onSelect;
    this.colourOf = colourOf;
    this.maxRoots = maxRoots;
    this.parcellation = null;
    this.selected = null;
    this.query = '';
    this.open = new Set();
    this.el.addEventListener('click', (e) => {
      const row = e.target.closest('.row');
      if (!row) return;
      const index = Number(row.dataset.index);
      const region = this.parcellation.all[index];
      if (e.target.classList.contains('twist')) { this.toggle(region); return; }
      this.onSelect(region);
    });
  }

  setParcellation(p) {
    this.parcellation = p;
    this.open = new Set();
    // Open the first two levels.
    for (const r of p.all) if (r.depth < 1) this.open.add(r.index);
    this.selected = null;
    this.render();
  }

  select(region) {
    this.selected = region;
    if (region) for (const a of region.ancestors()) this.open.add(a.index);
    this.render();
    const row = this.el.querySelector('.row.selected');
    if (row) row.scrollIntoView({ block: 'nearest' });
  }

  toggle(region) {
    if (this.open.has(region.index)) this.open.delete(region.index); else this.open.add(region.index);
    this.render();
  }

  search(query) {
    this.query = query.trim();
    this.render();
  }

  render() {
    const p = this.parcellation;
    if (!p) { this.el.innerHTML = ''; return; }
    if (this.query) {
      const found = p.find(this.query, { limit: 80 });
      this.el.innerHTML = `<ul>${found.map((r) => this.row(r, true)).join('')}${found.length === 80 ? '<li class="more">… more; refine the search</li>' : ''}${found.length ? '' : '<li class="more">nothing found</li>'}</ul>`;
      return;
    }
    const node = (r) => {
      const open = this.open.has(r.index);
      const children = r.children.length && open ? `<ul>${r.children.map(node).join('')}</ul>` : '';
      return `<li class="${open ? '' : 'collapsed'}">${this.row(r, false, open)}${children}</li>`;
    };
    const roots = p.roots.slice(0, this.maxRoots);
    this.el.innerHTML = `<ul>${roots.map(node).join('')}${p.roots.length > roots.length ? `<li class="more">and ${p.roots.length - roots.length} more</li>` : ''}</ul>`;
  }

  row(r, flat, open = false) {
    const colour = rgb(this.colourOf(r));
    const twist = r.children.length ? `<span class="twist">${open ? '▾' : '▸'}</span>` : '<span class="twist leaf">·</span>';
    const count = r.children.length ? `<span class="count">${r.leaves().length}</span>` : '';
    const name = flat ? r.path().map((x) => x.name).join(' › ') : r.name;
    const dim = !colour && r.isLeaf ? ' dim' : '';
    return `<div class="row${r === this.selected ? ' selected' : ''}${dim}" data-index="${r.index}" title="${esc(r.path().map((x) => x.name).join(' › '))}">${flat ? '' : twist}<span class="swatch${colour ? '' : ' none'}" style="${colour ? `background:${colour}` : ''}"></span><span class="text">${esc(name)}${r.acronym && r.acronym !== r.name ? ` <small style="color:var(--text-faint)">${esc(r.acronym)}</small>` : ''}</span>${count}</div>`;
  }
}
