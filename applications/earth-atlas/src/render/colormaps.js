// Colour maps: perceptually uniform, ordered from the low end to the high end.
//
// thermal ... rain: cmocean (Thyng, K. M., Greene, C. A., Hetland, R. D., Zimmerle,
// H. M., and DiMarco, S. F. (2016), "True colors of oceanography: guidelines for
// effective and accurate colormap selection", Oceanography 29(3), 9-13),
// MIT licence, Copyright (c) 2015 Kristen M. Thyng; 33 of the 256 entries of
// each table, interpolated linearly in sRGB between them.
// orchid: the laboratory's own, from plum through Dusty Orchid (#CA8DBE) to
// cream, with OKLab lightness rising evenly from 0.22 to 0.98.

const STOPS = {
  thermal: '042333 052943 082e54 0e3167 19337c 273491 39339e 49369f 563b9c 624198 6d4694 784b91 83508f 8e548c 995889 a55c86 b15f82 bb627d c76676 d36a6e de7065 e7765b ef7e52 f58749 f99242 fb9e3e fcaa3c fbb73d f9c441 f7d246 f3df4c eeed53 e8fa5b',
  balance: '181c43 20275b 263276 293c93 2548b0 1059be 0f6abd 237aba 3888ba 4e95ba 65a2bc 7eaebf 98bac5 b0c5cd c7d2d7 dee0e1 f1eceb eadcd8 e3c9c1 ddb6aa d8a493 d2927c cd7f65 c76d50 c05a3c b9452c af3024 a11d25 901029 7b0e29 650f24 510d1c 3c0912',
  matter: 'feedb0 fde2a5 fbd79b facc91 f9c187 f8b67e f6ab76 f4a06e f29567 ef8a60 ed7f5b e97457 e66a54 e15f53 dc5553 d54c54 ce4356 c73d59 bf365b b62f5e ac2a60 a32562 992263 8e1e63 841c63 791a61 6e195f 63185b 591757 4e1651 44144b 391245 2f0f3e',
  haline: '2a186c 2d1a7f 2e1d93 2a25a1 2034a2 16409d 0f4a98 0c5394 0f5b90 15638e 1c6a8c 22718a 297789 2e7e89 338588 388c88 3c9387 3f9986 44a184 49a882 4faf7e 56b77a 5fbe75 6bc46e 78cb68 89d061 9bd55c afd95d c1dd64 d2e16f e1e57c f0ea8a fdef9a',
  dense: 'e6f1f1 d8eaec cbe4e8 bedde5 b1d6e3 a5cfe2 9ac8e2 8fc1e2 86b9e3 7fb1e4 79a9e4 75a0e5 7397e4 738ee2 7584df 767bdb 7871d5 7969cf 7960c7 7957bd 784eb3 7746a9 753d9d 723692 6f2e85 6b2879 66216c 601c5e 5a1851 521544 491238 40102e 360e24',
  tempo: 'fff6f4 f4eee8 e9e7dd dde0d1 d2d9c7 c6d3bc bacdb3 aec7aa a1c1a1 94bc9a 86b693 78b18e 69ab89 59a685 49a083 3a9a81 2a937f 1f8d7e 16867c 117f7a 117777 147074 176870 19616c 1b5968 1c5263 1c4b5f 1c435a 1b3c56 1a3551 182d4c 172548 151d44',
  deep: 'fdfecc edf7c3 ddf2ba ccecb3 bbe6ac aae0a8 98daa4 87d4a3 78cea3 6bc6a3 62bfa4 5ab7a4 55aea3 51a6a2 4e9ea1 4b96a0 488e9e 46879c 437f9a 417799 3f6f97 3e6795 3e5f93 3e5691 404d8c 414583 413e76 3f3869 3b325c 372c4f 322643 2d2037 281a2c',
  speed: 'fffdcd f9f3bc f4eaab efe19a e9d889 e3cf77 dbc866 d3c155 c9ba45 beb436 b2af28 a5aa1c 98a512 8aa009 7c9c06 6e9707 5f920c 528d12 448818 36821e 287c23 1c7527 126e2a 0c672c 0b602d 0e582d 11512c 15492a 174127 193a23 19321e 192b19 172313',
  ice: '040613 0b0d1e 13132a 1a1a35 212041 28264e 2e2c5a 333267 383975 3c3f82 3e478f 3f4f9a 3f57a3 3e60ab 3e69b0 4072b4 427bb7 4682ba 4b8bbd 5194c0 579cc3 5fa4c6 67adc9 6fb5cd 79bed0 85c6d3 92ced7 a0d6dc afdde1 bee5e7 cdecef dcf4f6 eafdfd',
  rain: 'eeedf3 e8e5e6 e3dcda ded4cb dbcbbc d6c3ac cfbc9c c3b791 b5b38a a8ae85 9aaa80 8ca57c 7ea177 6f9c74 5f9871 4f936e 3d8e6e 2d896e 1c836e 0f7c6e 06756e 036e6d 05676c 0c5f6a 145867 1a5063 1f485d 234157 253951 25324a 252a44 23233e 221b38',
  orchid: '24142b 3d1f48 5c2d66 7d4180 a2629f ca8dbe e3b5d3 f4dbe7 fdf1ec fff8e4',
};

const parse = (s) => s.split(' ').map((h) => [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)]);
export const COLORMAPS = Object.fromEntries(Object.entries(STOPS).map(([k, v]) => [k, parse(v)]));

// The colour (r, g, b in 0..255) at t in [0, 1].
export function colourAt(name, t) {
  const stops = COLORMAPS[name];
  const x = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(x)), f = x - i;
  return stops[i].map((c, k) => c + (stops[i + 1][k] - c) * f);
}

// A 256-entry RGBA lookup table. Options: reverse; clip [t0, t1] of the map to
// use; fade, the fraction of the low end over which opacity rises from 0 to 1.
export function lookupTable(name, { reverse = false, clip = [0, 1], fade = 0 } = {}) {
  const out = new Uint8Array(256 * 4);
  for (let i = 0; i < 256; i++) {
    const s = i / 255;
    const t = clip[0] + (clip[1] - clip[0]) * (reverse ? 1 - s : s);
    const [r, g, b] = colourAt(name, t);
    out.set([r, g, b, fade > 0 ? Math.round(255 * Math.min(1, s / fade)) : 255], 4 * i);
  }
  return out;
}

// A CSS linear-gradient for a legend or a swatch.
export function cssGradient(name, options = {}, direction = 'to right') {
  const lut = lookupTable(name, options);
  const parts = [];
  for (let i = 0; i <= 16; i++) {
    const k = Math.round((i / 16) * 255) * 4;
    parts.push(`rgba(${lut[k]}, ${lut[k + 1]}, ${lut[k + 2]}, ${(lut[k + 3] / 255).toFixed(2)}) ${(i / 16) * 100}%`);
  }
  return `linear-gradient(${direction}, ${parts.join(', ')})`;
}
