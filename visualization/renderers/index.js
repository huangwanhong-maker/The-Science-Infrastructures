// The public interface of the renderers' arithmetic, without a device (ADR 0005; the infrastructure's
// ADR 0053). Each technique's own renderer is a directory beside it: three/ (three.js) today.

export * from './raster.js'; // a value's colour under an expression, drapes as images, curves as coloured segments
