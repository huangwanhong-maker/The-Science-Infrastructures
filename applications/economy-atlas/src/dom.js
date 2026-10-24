// Small helpers for building the page: elements, SVG, and numbers for reading.

const SVG = 'http://www.w3.org/2000/svg';

function build(node, attrs, children) {
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === undefined || v === null || v === false) continue;
    if (k === 'text') node.textContent = v;
    else if (k === 'class') node.setAttribute('class', v);
    else if (k.startsWith('on')) node.addEventListener(k.slice(2), v);
    else if (k === 'style' && typeof v === 'object') Object.assign(node.style, v);
    else node.setAttribute(k, v === true ? '' : v);
  }
  for (const c of children.flat(Infinity)) {
    if (c === null || c === undefined || c === false) continue;
    node.append(c instanceof Node ? c : document.createTextNode(String(c)));
  }
  return node;
}

export const el = (tag, attrs = {}, ...children) => build(document.createElement(tag), attrs, children);
export const svg = (tag, attrs = {}, ...children) => build(document.createElementNS(SVG, tag), attrs, children);

export function clear(node) {
  while (node.firstChild) node.removeChild(node.firstChild);
  return node;
}

// Children appended as el() appends them: arrays flattened, null and false left out.
export const put = (node, ...children) => build(node, {}, children);
export const fill = (node, ...children) => put(clear(node), ...children);

// Millions of US dollars, read as thousands of millions or millions of millions.
export function money(millions) {
  const a = Math.abs(millions);
  if (a >= 1e6) return `${(millions / 1e6).toFixed(2)} tn $`;
  if (a >= 1e3) return `${(millions / 1e3).toFixed(a >= 1e5 ? 0 : 1)} bn $`;
  return `${millions.toFixed(0)} m $`;
}

// A quantity in the data's own unit: thousands, millions and thousands of millions shortened.
export function quantity(x) {
  if (!Number.isFinite(x)) return '—';
  const a = Math.abs(x);
  if (a >= 1e9) return `${(x / 1e9).toFixed(2)} bn`;
  if (a >= 1e6) return `${(x / 1e6).toFixed(2)} m`;
  if (a >= 1e4) return `${(x / 1e3).toFixed(1)} k`;
  return x.toLocaleString('en-GB', { maximumFractionDigits: a >= 100 ? 0 : 2 });
}

export const percent = (x, digits = 1) => `${(100 * x).toFixed(digits)}%`;
export const signed = (x, digits = 3) => `${x >= 0 ? '+' : '−'}${Math.abs(x).toFixed(digits)}`;
export const fixed = (x, digits = 3) => (Number.isFinite(x) ? x.toFixed(digits) : '—');

// A bar of shares with its key: parts [{label, value, colour}].
export function shareBar(parts) {
  const total = parts.reduce((s, p) => s + Math.max(0, p.value), 0) || 1;
  return el('div', {},
    el('div', { class: 'bar' }, parts.map((p) => el('span', { style: { width: `${(100 * Math.max(0, p.value)) / total}%`, background: p.colour }, title: `${p.label}: ${percent(p.value / total)}` }))),
    el('div', { class: 'bar-key' }, parts.map((p) => el('span', {}, el('i', { style: { background: p.colour } }), `${p.label} ${percent(p.value / total)}`))));
}
