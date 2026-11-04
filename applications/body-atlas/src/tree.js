// The rail: the systems (shown, glass or hidden), the parts as FMA's part-of or is-a tree, and the search.

import { systemColour } from './palette.js';

const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
// FMA's names as a reader expects them: a capital first, and a closing roman numeral in capitals (segment IV).
export const capital = (s) => (s ? (s[0].toUpperCase() + s.slice(1)).replace(/ (i{1,3}|iv|v|vi{1,3}|ix|x)$/, (m) => m.toUpperCase()) : s);

// ---- systems ----------------------------------------------------------------------------------------------

export function renderSystems(ul, systems, styles, { onToggle, onGlass }) {
  ul.innerHTML = '';
  for (const s of systems) {
    const li = document.createElement('li');
    li.dataset.key = s.key;
    li.style.setProperty('--sw', systemColour(s.key));
    li.innerHTML = `<span class="sw" style="background:${systemColour(s.key)}"></span><span class="nm">${esc(s.name)}</span>`
      + `<span class="n">${s.count}</span><button class="gl" title="Draw as glass"><svg class="i"><use href="#i-glass"/></svg></button>`;
    li.addEventListener('click', (ev) => {
      if (ev.target.closest('.gl')) { ev.stopPropagation(); onGlass(s.key); } else onToggle(s.key);
    });
    ul.appendChild(li);
  }
  updateSystems(ul, styles);
}

export function updateSystems(ul, styles) {
  for (const li of ul.children) {
    const style = styles[li.dataset.key];
    li.classList.toggle('off', style === 'off');
    li.classList.toggle('glass', style === 'glass');
  }
}

// ---- the parts tree ---------------------------------------------------------------------------------------

export class PartsTree {
  constructor(el, atlas, { onChoose }) {
    this.el = el;
    this.atlas = atlas;
    this.onChoose = onChoose;
    this.relation = 'partof';
    this.open = new Set();
    this.current = null;
  }

  setRelation(r) {
    this.relation = r;
    this.render();
    if (this.current) this.reveal(this.current);
  }

  children(id) {
    const list = (this.relation === 'isa' ? this.atlas.childrenIsa : this.atlas.childrenPartof).get(id) || [];
    return list;
  }

  render() {
    this.el.innerHTML = '';
    const roots = this.relation === 'isa' ? this.atlas.rootsIsa : this.atlas.rootsPartof;
    const ul = document.createElement('ul');
    for (const r of roots) ul.appendChild(this.node(r));
    this.el.appendChild(ul);
  }

  node(id) {
    const c = this.atlas.concepts.get(id);
    const kids = this.children(id);
    const li = document.createElement('li');
    li.dataset.id = id;
    if (!kids.length) li.classList.add('leaf');
    const n = (this.relation === 'isa' ? c.elementsIsa : c.elementsPartof).length;
    const colour = c.system ? systemColour(c.system) : 'transparent';
    li.innerHTML = `<div class="row${this.current === id ? ' on' : ''}"><svg class="tw"><use href="#i-chev"/></svg>`
      + `<span class="dot" style="background:${colour}"></span><span>${esc(capital(c.name))}</span>${n > 1 ? `<span class="cnt">${n}</span>` : ''}</div>`;
    const row = li.firstChild;
    row.querySelector('.tw').addEventListener('click', (ev) => { ev.stopPropagation(); this.toggle(li); });
    row.addEventListener('click', () => this.onChoose(id, this.relation));
    row.addEventListener('dblclick', () => this.toggle(li));
    if (this.open.has(id)) this.expand(li);
    return li;
  }

  toggle(li) {
    if (li.classList.contains('open')) { li.classList.remove('open'); this.open.delete(li.dataset.id); li.querySelector('ul')?.remove(); }
    else this.expand(li);
  }

  expand(li) {
    const id = li.dataset.id;
    const kids = this.children(id);
    if (!kids.length) return;
    this.open.add(id);
    li.classList.add('open');
    if (li.querySelector('ul')) return;
    const ul = document.createElement('ul');
    for (const k of kids) ul.appendChild(this.node(k));
    li.appendChild(ul);
  }

  // Opens the path to a concept in the present tree, marks it and scrolls to it.
  reveal(id) {
    this.current = id;
    for (const r of this.el.querySelectorAll('.row.on')) r.classList.remove('on');
    const parent = this.relation === 'isa' ? (x) => this.atlas.concepts.get(x)?.isa : (x) => this.atlas.concepts.get(x)?.partof;
    const path = [];
    for (let x = id; x; x = parent(x)) path.unshift(x);
    if (!path.length || !this.atlas.concepts.has(path[0])) return;
    let scope = this.el;
    for (const x of path) {
      const li = [...scope.querySelectorAll(':scope > ul > li')].find((l) => l.dataset.id === x);
      if (!li) return;
      if (x !== id) this.expand(li);
      else {
        li.firstChild.classList.add('on');
        li.firstChild.scrollIntoView({ block: 'nearest' });
      }
      scope = li;
    }
  }
}

// ---- search -----------------------------------------------------------------------------------------------

export function searchIndex(atlas) {
  const items = [];
  for (const c of atlas.concepts.values()) items.push({ id: c.id, name: c.name, low: c.name.toLowerCase(), n: Math.max(c.elementsIsa.length, c.elementsPartof.length), system: c.system });
  return items;
}

export function search(items, text, limit = 14) {
  const q = text.trim().toLowerCase();
  if (q.length < 2) return [];
  const words = q.split(/\s+/);
  const scored = [];
  for (const it of items) {
    if (!words.every((w) => it.low.includes(w))) continue;
    const score = (it.low === q ? 0 : it.low.startsWith(q) ? 1 : it.low.split(/\s+/).some((w) => w.startsWith(words[0])) ? 2 : 3) + it.low.length / 200 - Math.min(it.n, 50) / 400;
    scored.push([score, it]);
  }
  scored.sort((a, b) => a[0] - b[0]);
  return scored.slice(0, limit).map((s) => s[1]);
}

export function renderFound(ol, results, systemsByKey, onChoose) {
  ol.innerHTML = '';
  ol.hidden = !results.length;
  results.forEach((r, k) => {
    const li = document.createElement('li');
    if (k === 0) li.classList.add('on');
    const kind = r.n > 1 ? `${r.n} parts` : (systemsByKey[r.system]?.name || '');
    li.innerHTML = `<span>${esc(capital(r.name))}</span><span class="k">${esc(kind)}</span>`;
    li.addEventListener('mousedown', (ev) => { ev.preventDefault(); onChoose(r.id); });
    ol.appendChild(li);
  });
}
