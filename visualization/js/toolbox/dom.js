// Small helpers for a toolbox's page: elements and SVG, with a text attribute, and tables saved as CSV.

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

// A file handed to the reader: text made in the page, saved under a name.
export function download(name, text, type = 'text/csv') {
  const url = URL.createObjectURL(new Blob([text], { type: `${type};charset=utf-8` }));
  const a = el('a', { href: url, download: name });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

// Rows as CSV: a header, then each row, cells quoted where they must be.
export function csv(header, rows) {
  const cell = (v) => {
    if (v === null || v === undefined || (typeof v === 'number' && !Number.isFinite(v))) return '';
    const s = String(v);
    return /[",\n;]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  return [header, ...rows].map((r) => r.map(cell).join(',')).join('\n') + '\n';
}
