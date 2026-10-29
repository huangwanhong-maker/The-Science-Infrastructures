// The public interface of the shared styles (ADR 0005; the infrastructure's ADR 0053): colour in OKLab,
// the warm family's colour maps, and colour-vision deficiencies simulated to validate them.

export * from './colour.js';   // sRGB and OKLab, ramps, tables for textures, gradients for legends
export * from './palettes.js'; // the named colour maps, by kind
export * from './cvd.js';      // Machado et al.'s simulation of protanopia, deuteranopia and tritanopia
