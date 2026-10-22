// The tree of structures: collapsible, with each structure's colour, an eye to
// show or hide it, and a click to choose it. It reads the catalogue of the domain.

import { esc } from './card.js';

export class StructureTree {
  // nodes: [{id, name, parent, colour, drawable}]; isVisible(id) -> true | false | null (no surface);
  // events: onSelect(id), onToggle(id)
  constructor(el, { nodes, isVisible, onSelect, onToggle }) {
    this.el = el;
    this.nodes = nodes;
    this.byId = new Map(nodes.map((n) => [n.id, n]));
    this.children = new Map();
    for (const n of nodes) {
      if (!this.children.has(n.parent)) this.children.set(n.parent, []);
      this.children.get(n.parent).push(n);
    }
    this.isVisible = isVisible;
    this.onSelect = onSelect;
    this.onToggle = onToggle;
    this.open = new Set(['brain', 'cerebrum', 'diencephalon', 'brainstem', 'cerebellum', 'ventricular-system']);
    this.selected = null;
    this.matches = null;
    el.addEventListener('click', (e) => {
      const row = e.target.closest('.row');
      if (!row) return;
      const id = row.dataset.id;
      if (e.target.closest('.eye')) { this.onToggle(id); return; }
      if (e.target.closest('.twist')) { this.toggleOpen(id); return; }
      this.onSelect(id);
    });
  }

  toggleOpen(id) {
    if (this.open.has(id)) this.open.delete(id); else this.open.add(id);
    this.render();
  }

  select(id) {
    this.selected = id;
    for (let n = this.byId.get(id); n; n = this.byId.get(n.parent)) if (n.id !== id) this.open.add(n.id);
    this.render();
    const row = this.el.querySelector('.row.selected');
    if (row) row.scrollIntoView({ block: 'nearest' });
  }

  // Show only the structures matching a text, with their ancestors open.
  filter(text) {
    const q = text.trim().toLowerCase();
    if (!q) { this.matches = null; this.render(); return; }
    this.matches = new Set();
    const keep = new Set();
    for (const n of this.nodes) {
      if (n.name.toLowerCase().includes(q) || (n.also || []).some((a) => a.toLowerCase().includes(q))) {
        this.matches.add(n.id);
        for (let p = n; p; p = this.byId.get(p.parent)) keep.add(p.id);
      }
    }
    this.keep = keep;
    this.render();
  }

  render() {
    const node = (n) => {
      const kids = (this.children.get(n.id) || []).filter((c) => !this.matches || this.keep.has(c.id));
      const open = this.matches ? true : this.open.has(n.id);
      const vis = this.isVisible(n.id);
      const eye = vis === null ? '<span class="eye" style="visibility:hidden"></span>'
        : `<span class="eye${vis ? '' : ' off'}" title="${vis ? 'Hide' : 'Show'}"><svg><use href="#${vis ? 'i-eye' : 'i-eye-off'}"/></svg></span>`;
      const sw = n.colour ? `<span class="sw" style="background:${esc(n.colour)}"></span>` : '<span class="sw group"></span>';
      const twist = `<span class="twist${kids.length ? ' has' : ' leaf'}"><svg><use href="#i-chev"/></svg></span>`;
      const cls = ['row', n.id === this.selected ? 'selected' : '', this.matches && this.matches.has(n.id) ? 'match' : '', vis === false ? 'dim' : ''].join(' ');
      const inner = kids.length && open ? `<ul>${kids.map(node).join('')}</ul>` : '';
      return `<li class="${open ? '' : 'collapsed'}"><div class="${cls}" data-id="${esc(n.id)}">${twist}${sw}<span class="name">${esc(n.name)}</span>${eye}</div>${inner}</li>`;
    };
    const roots = (this.children.get('brain') || []).filter((c) => !this.matches || this.keep.has(c.id));
    this.el.innerHTML = `<ul>${roots.map(node).join('')}</ul>`;
  }
}
