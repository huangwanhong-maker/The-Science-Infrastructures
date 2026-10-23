// Making elements. Text from canons and from the engine is set as text, never as markup.

export function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  return fill(node, attrs, children);
}

const SVG = 'http://www.w3.org/2000/svg';
export function svg(tag, attrs = {}, ...children) {
  const node = document.createElementNS(SVG, tag);
  return fill(node, attrs, children);
}

function fill(node, attrs, children) {
  for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === 'class') node.setAttribute('class', v);
    else if (k === 'text') node.textContent = v;
    else if (k.startsWith('on')) node.addEventListener(k.slice(2), v);
    else if (k === 'dataset') Object.assign(node.dataset, v);
    else node.setAttribute(k, v === true ? '' : v);
  }
  for (const c of children.flat()) if (c !== null && c !== undefined && c !== false) node.append(c instanceof Node ? c : document.createTextNode(String(c)));
  return node;
}

export const clear = (node) => { while (node.firstChild) node.firstChild.remove(); return node; };
export const $ = (id) => document.getElementById(id);
