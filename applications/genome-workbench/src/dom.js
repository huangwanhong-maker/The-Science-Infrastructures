import { icon } from './ui/icons.js';

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

// A notice for a few seconds, as a toast in the corner (the kit's toast; views written before the
// kit call it notice).
let toastBox = null;
export function toast(text, kind = '', ms = 5200) {
  if (!toastBox) { toastBox = el('div', { class: 'toasts', 'aria-live': 'polite' }); document.body.append(toastBox); }
  const t = el('div', { class: `toast ${kind}`.trim() }, icon(kind === 'error' ? 'warn' : kind === 'good' ? 'check' : 'info', 'sm'), el('div', {}, text));
  toastBox.append(t);
  setTimeout(() => { t.style.transition = 'opacity .3s'; t.style.opacity = '0'; setTimeout(() => t.remove(), 320); }, ms);
  return t;
}
export const notice = (text, error = false, ms = 5200) => toast(text, error ? 'error' : '', ms);

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
