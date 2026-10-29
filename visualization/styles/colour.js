// Colour: sRGB and OKLab, and colour maps made in OKLab (ADR 0005; the infrastructure's ADR 0053).
//
// Colours are [r, g, b] in [0, 1] of sRGB. OKLab (Ottosson 2020, "A perceptual color space for image
// processing") is a space whose Euclidean distances approximate perceived differences and whose L is a
// perceived lightness from 0 (black) to 1 (white); a map interpolated in it changes by even steps, and a
// map whose L is monotone reads in order without its hue (for a reader who cannot tell the hues apart,
// or on a grey print). The conversions are Ottosson's, through linear sRGB (IEC 61966-2-1's transfer
// function).

export const toLinear = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
export const toGamma = (c) => (c <= 0.0031308 ? 12.92 * c : 1.055 * c ** (1 / 2.4) - 0.055);

// sRGB → OKLab [L, a, b].
export function oklab([r, g, b]) {
  const R = toLinear(r), G = toLinear(g), B = toLinear(b);
  const l = Math.cbrt(0.4122214708 * R + 0.5363325363 * G + 0.0514459929 * B);
  const m = Math.cbrt(0.2119034982 * R + 0.6806995451 * G + 0.1073969566 * B);
  const s = Math.cbrt(0.0883024619 * R + 0.2817188376 * G + 0.6299787005 * B);
  return [0.2104542553 * l + 0.793617785 * m - 0.0040720468 * s, 1.9779984951 * l - 2.428592205 * m + 0.4505937099 * s, 0.0259040371 * l + 0.7827717662 * m - 0.808675766 * s];
}

// OKLab → sRGB, clamped to the gamut's cube.
export function rgbOfLab([L, a, b]) {
  const l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3, m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3, s = (L - 0.0894841775 * a - 1.291485548 * b) ** 3;
  const clamp = (x) => Math.min(1, Math.max(0, toGamma(x)));
  return [clamp(4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s), clamp(-1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s), clamp(-0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s)];
}

export const hexRgb = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16) / 255);
export const cssColour = (rgb) => `#${rgb.map((c) => Math.round(Math.min(1, Math.max(0, c)) * 255).toString(16).padStart(2, '0')).join('')}`;

// The perceived difference of two colours: their distance in OKLab (a just noticeable difference is
// some 0.02; the dataviz convention of ×100 makes it 2).
export function deltaE(x, y) {
  const a = oklab(x), b = oklab(y);
  return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

// The OKLab lightness of a colour.
export const lightness = (rgb) => oklab(rgb)[0];

// A ramp through stops (hex), interpolated in OKLab: t in [0, 1] → [r, g, b]; cyclic ramps wrap.
export function ramp(stops, { cyclic = false } = {}) {
  const labs = stops.map((h) => oklab(hexRgb(h)));
  const n = cyclic ? labs.length : labs.length - 1;
  return (t) => {
    t = cyclic ? ((t % 1) + 1) % 1 : Math.min(1, Math.max(0, Number.isFinite(t) ? t : 0));
    const x = t * n, i = Math.min(Math.floor(x), n - 1), f = x - i;
    const a = labs[i], b = labs[(i + 1) % labs.length];
    return rgbOfLab([a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f]);
  };
}

// A colour function sampled into an RGBA table of `steps` entries over u in [0, 1] (for a texture's
// lookup): colour(u) → [r, g, b], alpha(u) → [0, 1].
export function table(colour, steps = 256, alpha = () => 1) {
  const out = new Uint8ClampedArray(4 * steps);
  for (let k = 0; k < steps; k++) {
    const u = steps > 1 ? k / (steps - 1) : 0, c = colour(u);
    out[4 * k] = Math.round(255 * c[0]);
    out[4 * k + 1] = Math.round(255 * c[1]);
    out[4 * k + 2] = Math.round(255 * c[2]);
    out[4 * k + 3] = Math.round(255 * alpha(u));
  }
  return out;
}

// A colour function over u in [0, 1] as a CSS gradient, for a legend's bar.
export function gradientCss(colour, { stops = 9, direction = 'to right' } = {}) {
  const parts = Array.from({ length: stops }, (_, k) => {
    const u = k / (stops - 1);
    return `${cssColour(colour(u))} ${(100 * u).toFixed(1)}%`;
  });
  return `linear-gradient(${direction}, ${parts.join(', ')})`;
}
