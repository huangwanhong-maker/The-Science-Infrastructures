// The public interface of the scene model (ADR 0005; the infrastructure's ADR 0053): device-neutral
// descriptions of what is drawn, which every renderer reads.

export * from './grid.js';       // the samplings of physis's fields: places and orders of their points
export * from './samples.js';    // a field's answers decoded: values, qualities, reasons, notes; curves
export * from './scale.js';      // where a value falls on a colour map: linear, logarithmic, diverging
export * from './expression.js'; // visual expressions M ∘ ρ ∘ r, their defaults and the formalism's rules
export * from './legend.js';     // what an expression shows, in what unit, frame and provenance
export * from './volume.js';     // the arithmetic of a field in a volume: cells shown, gradient, centre, level
