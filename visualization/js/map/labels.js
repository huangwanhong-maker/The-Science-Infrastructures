// The names on a map: placed from the most important down, each where it overlaps none placed before it
// (the greedy rule of label placement; Imhof 1975 on priority, Christensen, Marks and Shieber 1995 on
// its quality), and only those a naming policy lets the map draw. Pure arithmetic, checked in Node.
//
// A naming policy is the data library's (data/products/naming.json, ADR 0069): names for places by
// their code, and the codes of territories (and seas) whose names a map does not draw, their status
// being in dispute. The map draws a place's name from the policy where it has one.

// A policy from its JSON ({names: {code: name|null}, unlabelled: [codes], seas: [names]}); without one,
// every name is drawn as given.
export function namingPolicy({ names = {}, unlabelled = [], seas = [] } = {}) {
  const off = new Set(unlabelled), offSeas = new Set(seas);
  return {
    name: (code, given) => (code in names ? names[code] : given),
    labelled: (code) => !off.has(code) && !(code in names && names[code] === null),
    seaLabelled: (name) => !offSeas.has(name),
  };
}

// items: [{id, text, x, y}] in order of priority, in drawing units; px the drawing units of a screen
// pixel; size the type's size in pixels; the box of each name is estimated from its length. Returns the
// names placed, each with its box [x0, y0, x1, y1].
export function placeLabels(items, { px = 1, size = 12, charWidth = 0.55, lift = 6, limit = Infinity, policy = null } = {}) {
  const placed = [];
  for (const it of items) {
    if (placed.length >= limit) break;
    if (!it.text || (policy && !policy.labelled(it.id))) continue;
    const w = charWidth * size * px * it.text.length, h = (size + 1) * px;
    const box = [it.x - w / 2, it.y - lift * px - h, it.x + w / 2, it.y - lift * px + 2 * px];
    if (placed.some((p) => box[0] < p.box[2] && p.box[0] < box[2] && box[1] < p.box[3] && p.box[1] < box[3])) continue;
    placed.push({ ...it, box, y: it.y - lift * px });
  }
  return placed;
}
