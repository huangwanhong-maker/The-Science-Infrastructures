// Colours for maps of classes.
//
// A class map is read by its conventional colours: for the Köppen–Geiger
// classes, the scheme of Peel et al. (2007), which Beck et al. give with their
// maps (blue for the tropics, red and orange for the arid lands, yellow and
// green for the temperate, violet and cyan for the cold, grey for the polar).
// Those colours are pure primaries. Here each keeps its hue, but its chroma is
// reduced and its lightness drawn toward the middle, in the OKLab space
// (Ottosson 2020), so that the map sits with the atlas's gentle palette while
// remaining the map that readers know.

const toLinear = (c) => { c /= 255; return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4; };
const toGamma = (c) => 255 * (c <= 0.0031308 ? 12.92 * c : 1.055 * c ** (1 / 2.4) - 0.055);

function oklab([r, g, b]) {
  const [R, G, B] = [r, g, b].map(toLinear);
  const l = Math.cbrt(0.4122214708 * R + 0.5363325363 * G + 0.0514459929 * B);
  const m = Math.cbrt(0.2119034982 * R + 0.6806995451 * G + 0.1073969566 * B);
  const s = Math.cbrt(0.0883024619 * R + 0.2817188376 * G + 0.6299787005 * B);
  return [0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
    1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
    0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s];
}

function srgb([L, a, b]) {
  const l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3;
  const m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3;
  const s = (L - 0.0894841775 * a - 1.2914855480 * b) ** 3;
  return [4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
    -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
    -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s]
    .map((c) => Math.round(Math.max(0, Math.min(255, toGamma(Math.max(0, c))))));
}

// A colour with its chroma scaled by `chroma` and its lightness moved a
// fraction `pull` of the way toward `target`.
export function soften(rgb, { chroma = 0.6, pull = 0.34, target = 0.74 } = {}) {
  const [L, a, b] = oklab(rgb);
  return srgb([L + (target - L) * pull, a * chroma, b * chroma]);
}

// RGBA by class value (256 entries) for the shader, from a product's legend.
export function classColours(classes, alpha = 236) {
  const out = new Uint8Array(256 * 4);
  for (const c of classes) out.set([...soften(c.colour), alpha], 4 * c.value);
  return out;
}

export const cssColour = (rgb) => `rgb(${rgb.join(',')})`;
