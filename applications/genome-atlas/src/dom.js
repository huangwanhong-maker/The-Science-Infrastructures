// Making elements. Text from genomes, files and the engine is set as text, never as markup; the
// large drawings are built as strings with esc() on every piece of text in them.

export function el(tag, attrs = {}, ...children) {
  return fill(document.createElement(tag), attrs, children);
}

const SVG = 'http://www.w3.org/2000/svg';
export function svg(tag, attrs = {}, ...children) {
  return fill(document.createElementNS(SVG, tag), attrs, children);
}

function fill(node, attrs, children) {
  for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === 'class') node.setAttribute('class', v);
    else if (k === 'text') node.textContent = v;
    else if (k.startsWith('on')) node.addEventListener(k.slice(2), v);
    else if (k === 'dataset') Object.assign(node.dataset, v);
    else if (k === 'style' && typeof v === 'object') Object.assign(node.style, v);
    else node.setAttribute(k, v === true ? '' : v);
  }
  for (const c of children.flat()) if (c !== null && c !== undefined && c !== false) node.append(c instanceof Node ? c : document.createTextNode(String(c)));
  return node;
}

// Escapes text for markup built as a string.
const ESC = { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' };
export const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ESC[c]);

export const clear = (node) => { while (node.firstChild) node.firstChild.remove(); return node; };
export const $ = (id) => document.getElementById(id);

// A tooltip that follows the pointer; html is built by the caller with esc().
export function hover(html, event) {
  const h = $('hover');
  if (!html) { h.hidden = true; return; }
  h.innerHTML = html;
  h.hidden = false;
  const pad = 14, w = h.offsetWidth, ht = h.offsetHeight;
  let x = event.clientX + pad, y = event.clientY + pad;
  if (x + w > window.innerWidth - 8) x = event.clientX - w - pad;
  if (y + ht > window.innerHeight - 8) y = event.clientY - ht - pad;
  h.style.left = `${x}px`;
  h.style.top = `${y}px`;
}

// A notice at the foot of the page, for a few seconds.
let noticeTimer = null;
export function notice(text, error = false, ms = 5200) {
  const n = $('notice');
  n.textContent = text;
  n.className = `notice${error ? ' error' : ''}`;
  n.hidden = false;
  clearTimeout(noticeTimer);
  noticeTimer = setTimeout(() => { n.hidden = true; }, ms);
}

// A row of toggle buttons: options are [value, label, title]; onPick(value) on a choice.
export function toggles(options, value, onPick) {
  const box = el('div', { class: 'toggles', role: 'group' });
  options.forEach(([v, label, title], k) => {
    box.append(el('button', { class: v === value ? 'on' : '', title: title || null, onclick: () => {
      for (const b of box.children) b.classList.remove('on');
      box.children[k].classList.add('on');
      onPick(v);
    } }, label));
  });
  return box;
}

// A labelled field: field('k', input).
export const field = (label, control) => el('label', { class: 'field' }, el('span', {}, label), control);
