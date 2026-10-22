// Colour maps: perceptually uniform, ordered from the low end to the high end.
// The same maps as the Earth Atlas draws with, until they move to the shared
// visualization layer.
//
// thermal ... rain: cmocean (Thyng, K. M., Greene, C. A., Hetland, R. D., Zimmerle,
// H. M., and DiMarco, S. F. (2016), "True colors of oceanography: guidelines for
// effective and accurate colormap selection", Oceanography 29(3), 9-13),
// MIT licence, Copyright (c) 2015 Kristen M. Thyng; 33 of the 256 entries of
// each table, interpolated linearly in sRGB between them.
// orchid: our own, from plum through Dusty Orchid (#CA8DBE) to cream.

const STOPS = {
  orchid: '24142b 3d1f48 5c2d66 7d4180 a2629f ca8dbe e3b5d3 f4dbe7 fdf1ec fff8e4',
  thermal: '042333 052943 082e54 0e3167 19337c 273491 39339e 49369f 563b9c 624198 6d4694 784b91 83508f 8e548c 995889 a55c86 b15f82 bb627d c76676 d36a6e de7065 e7765b ef7e52 f58749 f99242 fb9e3e fcaa3c fbb73d f9c441 f7d246 f3df4c eeed53 e8fa5b',
  balance: '181c43 20275b 263276 293c93 2548b0 1059be 0f6abd 237aba 3888ba 4e95ba 65a2bc 7eaebf 98bac5 b0c5cd c7d2d7 dee0e1 f1eceb eadcd8 e3c9c1 ddb6aa d8a493 d2927c cd7f65 c76d50 c05a3c b9452c af3024 a11d25 901029 7b0e29 650f24 510d1c 3c0912',
  matter: 'feedb0 fde2a5 fbd79b facc91 f9c187 f8b67e f6ab76 f4a06e f29567 ef8a60 ed7f5b e97457 e66a54 e15f53 dc5553 d54c54 ce4356 c73d59 bf365b b62f5e ac2a60 a32562 992263 8e1e63 841c63 791a61 6e195f 63185b 591757 4e1651 44144b 391245 2f0f3e',
  haline: '2a186c 2d1a7f 2e1d93 2a25a1 2034a2 16409d 0f4a98 0c5394 0f5b90 15638e 1c6a8c 22718a 297789 2e7e89 338588 388c88 3c9387 3f9986 44a184 49a882 4faf7e 56b77a 5fbe75 6bc46e 78cb68 89d061 9bd55c afd95d c1dd64 d2e16f e1e57c f0ea8a fdef9a',
  deep: 'fdfecc edf7c3 ddf2ba ccecb3 bbe6ac aae0a8 98daa4 87d4a3 78cea3 6bc6a3 62bfa4 5ab7a4 55aea3 51a6a2 4e9ea1 4b96a0 488e9e 46879c 437f9a 417799 3f6f97 3e6795 3e5f93 3e5691 404d8c 414583 413e76 3f3869 3b325c 372c4f 322643 2d2037 281a2c',
  ice: '040613 0b0d1e 13132a 1a1a35 212041 28264e 2e2c5a 333267 383975 3c3f82 3e478f 3f4f9a 3f57a3 3e60ab 3e69b0 4072b4 427bb7 4682ba 4b8bbd 5194c0 579cc3 5fa4c6 67adc9 6fb5cd 79bed0 85c6d3 92ced7 a0d6dc afdde1 bee5e7 cdecef dcf4f6 eafdfd',
};

const parse = (s) => s.split(' ').map((h) => [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)]);
export const COLORMAPS = Object.fromEntries(Object.entries(STOPS).map(([k, v]) => [k, parse(v)]));
export const COLORMAP_NAMES = Object.keys(STOPS);

// The colour (r, g, b in 0..255) at t in [0, 1].
export function colourAt(name, t) {
  const stops = COLORMAPS[name] || COLORMAPS.orchid;
  const x = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(x)), f = x - i;
  return stops[i].map((c, k) => c + (stops[i + 1][k] - c) * f);
}

// A CSS linear-gradient for a legend or a swatch.
export function cssGradient(name, direction = 'to right') {
  const parts = [];
  for (let i = 0; i <= 16; i++) {
    const [r, g, b] = colourAt(name, i / 16);
    parts.push(`rgb(${Math.round(r)}, ${Math.round(g)}, ${Math.round(b)}) ${(i / 16) * 100}%`);
  }
  return `linear-gradient(${direction}, ${parts.join(', ')})`;
}

// Hues for categories, warm and distinct, cycling.
export const CATEGORY_COLOURS = ['#CA8DBE', '#f6c47e', '#9fb6e8', '#8fd1b4', '#e8a0a0', '#d9b8ff', '#f0d9a0', '#a9d8e0', '#e3b5d3', '#c9d99a', '#f3a6c2', '#b8c4f0'];
export const categoryColour = (i) => CATEGORY_COLOURS[i % CATEGORY_COLOURS.length];
export const hexToRgb = (hex) => [1, 3, 5].map((k) => parseInt(hex.slice(k, k + 2), 16));
