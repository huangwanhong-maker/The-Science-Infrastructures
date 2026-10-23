// The disciplines of the codebook (§6) and their colours, in their fixed order: the order of the
// loom's groups of lanes and of every legend. The six were validated on the page's paper (#F6EFE7)
// by the dataviz method (OKLCH lightness band and chroma floor; adjacent pairs in this order:
// worst CVD ΔE 12.2, normal vision 17.6; every one at least 3:1 against the paper). All pairs cannot
// be told apart under colour-blindness at once, so a colour never stands alone: lanes are grouped
// and labelled by discipline, and every legend names it.

export const COLOURS = [
  { id: 'music', name: 'Music', hex: '#C2457E' },
  { id: 'art', name: 'Art', hex: '#B5711C' },
  { id: 'international-relations', name: 'International relations', hex: '#2F6FB5' },
  { id: 'science', name: 'Science', hex: '#2F8F46' },
  { id: 'letters', name: 'Letters', hex: '#7656C2' },
  { id: 'philosophy', name: 'Philosophy', hex: '#77801C' },
];

// Reserved: frustration (a triangle of ties whose signs multiply to −1) is marked in vermilion, and
// in no other place; the page's own accent is a faded mauve, used for choices and cursors only.
export const FRUSTRATION = '#C8402B';
export const ACCENT = '#9C6B85';
export const INK = '#3A2C33';
export const PAPER = '#F6EFE7';

const UNKNOWN = { id: 'unknown', name: 'Unassigned', hex: '#8E7F86' };
export const colourOf = (id) => COLOURS.find((c) => c.id === id) || (id && id !== 'unknown' ? { ...UNKNOWN, id, name: id } : UNKNOWN);
export const colourIndex = (id) => { const i = COLOURS.findIndex((c) => c.id === id); return i < 0 ? COLOURS.length : i; };
