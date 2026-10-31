// Small helpers for building the page: h() for elements, svg() for SVG elements, and the hover note.

export function h(tag, attrs = {}, ...children) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === undefined || v === null || v === false) continue;
    if (k === 'class') el.className = v;
    else if (k === 'html') el.innerHTML = v;
    else if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v === true ? '' : v);
  }
  for (const c of children.flat()) if (c !== null && c !== undefined && c !== false) el.append(c instanceof Node ? c : String(c));
  return el;
}

const SVGNS = 'http://www.w3.org/2000/svg';
export function svg(tag, attrs = {}, ...children) {
  const el = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === undefined || v === null) continue;
    if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v);
  }
  for (const c of children.flat()) if (c !== null && c !== undefined) el.append(c instanceof Node ? c : document.createTextNode(String(c)));
  return el;
}

export function hover(el, text) {
  const box = document.getElementById('hover');
  el.addEventListener('mousemove', (e) => {
    box.textContent = typeof text === 'function' ? text() : text;
    box.hidden = false;
    box.style.left = `${Math.min(e.clientX + 14, window.innerWidth - 330)}px`;
    box.style.top = `${e.clientY + 14}px`;
  });
  el.addEventListener('mouseleave', () => { box.hidden = true; });
}

export function css(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

export function fmt(x, digits = 4) {
  if (typeof x !== 'number') return String(x);
  if (Number.isInteger(x)) return x.toLocaleString('en');
  const a = Math.abs(x);
  if (a !== 0 && (a < 1e-3 || a >= 1e6)) return x.toExponential(2);
  return x.toPrecision(digits).replace(/\.?0+$/, '');
}

export function busy(root, text = 'computing…') {
  root.replaceChildren(h('p', { class: 'busy' }, text));
}
export function failure(root, error) {
  root.replaceChildren(h('p', { class: 'error-box' }, String(error.message || error)));
}
