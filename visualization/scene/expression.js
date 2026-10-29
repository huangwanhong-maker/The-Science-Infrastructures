// Visual expressions of fields: the device-neutral part of the scene model (ADR 0005; the
// infrastructure's ADR 0053; Cosmopoiesis's formalism, Definition 3.10).
//
// A drawing of a field is not a kind of phenomenon but a composite map from the field to marks,
// V = M ∘ ρ ∘ r: a restriction r of the section to a sub-base it is drawn on, a reduction ρ of its fibre
// to what is shown, and a mark M whose geometry the sub-base decides and whose channels the reduced
// fibre decides. An expression here is that composite as data: serialisable, without a renderer, and
// checked by the formalism's rules before anything is drawn.
//
//   field        the field's descriptor as its engine gives it (physis fields: name, title, frame, time,
//                base, fibre {kind, unit, basis}, representation, provenance, parts, source)
//   restriction  the sub-base: a sampling {kind: sphere | plane | volume | points, …} (grid.js), or
//                {kind: 'curves', seeds, within?} for integral curves through the base
//   reduction    {kind: none | magnitude | component | radial | horizontal | on_axes, axis?, basis?}
//   mark         {kind: drape | cells | slices | shell | isosurface | arrows | curves | glyphs, …}
//   colour       {map, scale: auto | linear | log | diverging, centre?, floor?, opacity: constant | rising}
//
// The defaults follow the formalism's tables: a vector over a region is drawn as its integral curves,
// or reduced (to its magnitude, or a component) and coloured on a surface or in a volume; a magnitude is
// sequential, logarithmic where it spans decades; a signed component is diverging about zero; a label
// is categorical, never ordered. The rules (validate) are the formalism's too, among them its remark on
// invariance: a drawing is the same in every frame exactly when its reduction is an invariant or names
// its frame, so a component is never drawn without its frame.

import { DEFAULT_MAP, MAPS } from '../styles/palettes.js';

export const RESTRICTIONS = ['sphere', 'plane', 'volume', 'points', 'curves'];
export const REDUCTIONS = ['none', 'magnitude', 'component', 'radial', 'horizontal', 'on_axes'];
export const MARKS = { sphere: ['drape'], plane: ['drape'], volume: ['cells', 'slices', 'shell', 'isosurface', 'arrows'], points: ['glyphs'], curves: ['curves'] };

// Whether a reduction gives a quantity that is the same in every frame.
export const invariant = (kind) => ['none', 'magnitude', 'radial', 'horizontal'].includes(kind);
// Whether a reduction has a sign about a natural zero.
export const signed = (kind) => ['component', 'radial'].includes(kind);

const AXES = ['x', 'y', 'z'];
const isVector = (field) => field?.fibre?.kind === 'vector' || field?.fibre?.kind === 'covector';

// The name of what is shown: "|B|", "B_r", "B_x (ITRS)".
export function quantityName(field, reduction) {
  const symbol = field.symbol || ({ T: 'B', 'm s^-2': 'g', K: 'T', Pa: 'p', 'kg m^-3': 'ρ', 'm s^-1': 'v' }[field.fibre?.unit] ?? 'v');
  switch (reduction.kind) {
    case 'magnitude': return `|${symbol}|`;
    case 'radial': return `${symbol}_r`;
    case 'horizontal': return `${symbol}_h`;
    case 'component': return `${symbol}_${AXES[reduction.axis] ?? '?'} (${reduction.basis || field.frame})`;
    case 'on_axes': return `${symbol} on ${reduction.basis || 'the given axes'}`;
    default: return symbol;
  }
}

// An expression of a field on a restriction, its reduction, mark and colour chosen by the formalism's
// tables unless the options choose them.
export function express(field, restriction, options = {}) {
  if (!field || !field.fibre) throw new Error('express: a field is its descriptor, with its fibre');
  if (!RESTRICTIONS.includes(restriction?.kind)) throw new Error(`express: the restriction is one of ${RESTRICTIONS.join(', ')}`);
  const vector = isVector(field);
  const label = field.fibre.kind === 'label';
  let reduction = options.reduction || null;
  if (!reduction) {
    if (!vector || restriction.kind === 'curves') reduction = { kind: 'none' };
    else if (restriction.kind === 'points') reduction = { kind: 'none' };
    else reduction = { kind: 'magnitude' };
  }
  if (reduction.kind === 'component' && !reduction.basis) reduction = { ...reduction, basis: field.frame };
  const mark = options.mark || { kind: restriction.kind === 'curves' ? 'curves' : MARKS[restriction.kind][0] };
  // The colour: of the scalar drawn, or along a curve of its |v|.
  const scale = options.scale || (label ? 'categorical' : signed(reduction.kind) ? 'diverging' : options.about !== undefined ? 'diverging' : 'auto');
  const kind = label ? 'categorical' : scale === 'diverging' ? 'diverging' : 'sequential';
  const colour = {
    map: options.map || DEFAULT_MAP[kind],
    scale,
    centre: scale === 'diverging' ? (options.about ?? 0) : null,
    floor: options.floor ?? null,
    opacity: options.opacity || 'constant',
  };
  const expression = {
    field, restriction, reduction, mark, colour,
    // Curves are coloured by the vector's magnitude along them.
    shows: restriction.kind === 'curves' ? `${quantityName(field, { kind: 'magnitude' })} along its curves` : quantityName(field, reduction),
    unit: field.fibre.unit,
    frame: invariant(reduction.kind) ? null : reduction.basis || field.frame,
    title: options.title || `${quantityName(field, reduction)}: ${field.title}`,
  };
  const problems = validate(expression);
  if (problems.length) throw new Error(`express: ${problems.join('; ')}`);
  return expression;
}

// The formalism's rules for an expression: the problems that keep it from being drawn, none if it may be.
export function validate(e) {
  const out = [];
  const vector = isVector(e.field);
  const { restriction: r, reduction: rho, mark, colour } = e;
  if (!REDUCTIONS.includes(rho.kind)) out.push(`no reduction "${rho.kind}"`);
  if (!vector && rho.kind !== 'none') out.push(`a scalar is drawn as it is; "${rho.kind}" reduces a vector`);
  if (r.kind === 'curves') {
    if (!vector) out.push('integral curves are of a vector field');
    if (mark.kind !== 'curves') out.push('a restriction to curves is drawn as curves');
    if (rho.kind !== 'none') out.push('a curve follows the vector itself, not a reduction of it');
  } else {
    if (!(MARKS[r.kind] || []).includes(mark.kind)) out.push(`a ${r.kind} is drawn as ${MARKS[r.kind].join(' or ')}, not as ${mark.kind}`);
    const scalarShown = !vector || !['none', 'on_axes'].includes(rho.kind);
    if (!scalarShown && !['arrows', 'glyphs'].includes(mark.kind)) out.push("a vector's colour needs a reduction to a scalar: its magnitude, a component, its radial or horizontal part");
    if (['arrows', 'glyphs'].includes(mark.kind) && vector && !['none', 'on_axes'].includes(rho.kind)) out.push('arrows draw a vector, not a reduction of it to a scalar');
  }
  if (!invariant(rho.kind) && !(e.frame || rho.basis)) out.push('a component is never drawn without its frame');
  if (colour) {
    const m = MAPS[colour.map];
    if (!m) out.push(`no colour map "${colour.map}"`);
    if (e.field.fibre.kind === 'label' && m && m.kind !== 'categorical') out.push('a label is drawn in categorical hues, never on an ordered map');
    if (colour.scale === 'diverging') {
      if (colour.centre === null || colour.centre === undefined) out.push('a diverging scale needs its centre');
      if (m && m.kind !== 'diverging') out.push(`a diverging scale needs a diverging map, not ${m.kind} "${colour.map}"`);
      if (['magnitude', 'horizontal'].includes(rho.kind) && colour.centre === 0) out.push('a magnitude has no sign: about zero its scale is sequential');
      if (e.field.fibre.unit === 'K' && rho.kind === 'none' && colour.centre === 0) out.push('an absolute temperature has no natural zero to diverge about: give a reference');
    } else if (m && m.kind === 'diverging') out.push('a diverging map needs a diverging scale about a centre');
    if (colour.scale === 'log' && signed(rho.kind)) out.push('a signed quantity has no logarithmic scale: take a diverging one with a floor');
  }
  return out;
}

// The expressions a field offers by its base and fibre, for a menu: each a label and the options of
// express() with a restriction to be completed by the caller (a sphere's radius, a volume's centre).
export function offers(field) {
  if (!field?.fibre || field.available === false) return [];
  const vector = isVector(field);
  const out = [];
  if (vector) {
    out.push({ label: 'field lines', restriction: 'curves', options: {} });
    out.push({ label: 'its magnitude on a sphere', restriction: 'sphere', options: { reduction: { kind: 'magnitude' } } });
    out.push({ label: 'its radial part on a sphere', restriction: 'sphere', options: { reduction: { kind: 'radial' } } });
    out.push({ label: 'in a volume about a point', restriction: 'volume', options: { reduction: { kind: 'magnitude' } } });
  } else {
    out.push({ label: 'on a sphere', restriction: 'sphere', options: {} });
    out.push({ label: 'on a slice', restriction: 'plane', options: {} });
    out.push({ label: 'in a volume about a point', restriction: 'volume', options: {} });
  }
  return out;
}
