// A description's derivation drawn over its text: every node of the named parse tree a bracket
// over the span it derives, leaves just above the letters and the start at the top, so that one
// reads the grammar's analysis of the text as a constituency tree laid flat.
//
//   derivationSvg({ text, tree, width?, highlight? })
//
// tree is morphe's named parse tree ({ name, begin, end, text, children }); a node's row is its
// height (a node without children stands on the text). A rejected text is drawn with its viable
// prefix underlined and the place of the fault marked: derivationSvg({ text, furthest, expected }).

import { esc } from '../dom.js';
import { token, CATEGORIES } from '../palette.js';

const CW = 8.2;      // a letter's width at 13.5 px of the monospaced face
const ROW = 21;      // a row of brackets

function height(n) { return n.children?.length ? 1 + Math.max(...n.children.map(height)) : 1; }

// A rule's colour, the same for one name wherever it is drawn.
const colourCache = new Map();
export function ruleColour(name) {
  if (!colourCache.has(name)) {
    let h = 0;
    for (const c of name) h = (h * 31 + c.charCodeAt(0)) >>> 0;
    colourCache.set(name, CATEGORIES[h % CATEGORIES.length]);
  }
  return colourCache.get(name);
}

export function derivationSvg({ text, tree = null, furthest = null, expected = [], highlight = null }) {
  const n = text.length;
  const pad = 14;
  const textW = Math.max(1, n) * CW;
  const H = tree ? height(tree) : 0;
  const top = 8, textY = top + H * ROW + 18;
  const width = Math.max(240, textW + 2 * pad);
  const out = [];
  const ink = token('--ink'), faint = token('--ink-faint');
  const xOf = (i) => pad + i * CW;
  // The nodes, by height.
  const nodes = [];
  const walk = (node, depth) => {
    const h = height(node);
    nodes.push({ node, h, depth });
    for (const c of node.children || []) walk(c, depth + 1);
  };
  if (tree) walk(tree, 0);
  for (const { node, h, depth } of nodes) {
    const x0 = xOf(node.begin), x1 = xOf(node.end);
    const y = top + (H - h) * ROW;
    const col = ruleColour(node.name);
    const lit = highlight && highlight.begin === node.begin && highlight.end === node.end && highlight.name === node.name;
    const wide = x1 - x0;
    const id = `${node.name}:${node.begin}:${node.end}`;
    out.push(`<g class="dnode" data-node="${esc(id)}" data-begin="${node.begin}" data-end="${node.end}" style="cursor:pointer">`);
    out.push(`<title>${esc(`${node.name}  ·  ${node.text}`)}</title>`);
    out.push(`<rect x="${(x0 + 0.8).toFixed(1)}" y="${(y + 2).toFixed(1)}" width="${Math.max(2, wide - 1.6).toFixed(1)}" height="${ROW - 5}" rx="4" fill="${col}" fill-opacity="${lit ? 0.32 : depth === 0 ? 0.07 : 0.13}" stroke="${col}" stroke-opacity="${lit ? 0.95 : 0.55}" stroke-width="${lit ? 1.4 : 0.8}" />`);
    // The rule's name where it fits, shortened where it nearly does.
    const room = Math.floor((wide - 8) / 6.4);
    if (room >= 3) {
      const label = node.name.length <= room ? node.name : `${node.name.slice(0, Math.max(2, room - 1))}…`;
      out.push(`<text x="${((x0 + x1) / 2).toFixed(1)}" y="${(y + ROW / 2 + 2.5).toFixed(1)}" text-anchor="middle" style="font: italic 11.5px var(--serif); fill:${ink}; opacity:.82">${esc(label)}</text>`);
    }
    out.push('</g>');
  }
  // The text, letter by letter, so that each sits under its brackets.
  for (let i = 0; i < n; i++) {
    const bad = furthest != null && !tree && i >= furthest;
    out.push(`<text x="${(xOf(i) + CW / 2).toFixed(1)}" y="${textY}" text-anchor="middle" style="font: 13.5px var(--mono); fill:${bad ? token('--bad') : ink}${bad && i > furthest ? ';opacity:.45' : ''}">${esc(text[i] === ' ' ? '·' : text[i])}</text>`);
  }
  let height_ = textY + 12;
  if (furthest != null && !tree) {
    // The viable prefix underlined; a caret at the fault, with what could have come.
    out.push(`<line x1="${xOf(0)}" y1="${textY + 6}" x2="${xOf(furthest)}" y2="${textY + 6}" stroke="${token('--good')}" stroke-width="2" stroke-linecap="round" />`);
    const cx = xOf(furthest) + (furthest < n ? CW / 2 : 2);
    out.push(`<path d="M${cx - 5},${textY + 18} L${cx},${textY + 9} L${cx + 5},${textY + 18}" fill="none" stroke="${token('--bad')}" stroke-width="1.6" stroke-linejoin="round" />`);
    if (expected.length) {
      const shown = expected.slice(0, 8).join('  ') + (expected.length > 8 ? `  … ${expected.length - 8} more` : '');
      out.push(`<text x="${Math.min(cx, width - 10)}" y="${textY + 34}" text-anchor="${cx > width * 0.6 ? 'end' : 'start'}" style="font: italic 12px var(--serif); fill:${faint}">expected ${esc(shown)}</text>`);
    }
    height_ = textY + 44;
  }
  return `<svg class="derivation" width="${width}" height="${height_}" viewBox="0 0 ${width} ${height_}" xmlns="http://www.w3.org/2000/svg">${out.join('')}</svg>`;
}
