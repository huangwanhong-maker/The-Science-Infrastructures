// Tiles of an equirectangular map service (NASA's Trek WMTS): which tiles,
// at which level, cover a region of a body at a given resolution.
//
// Level z has 2^(z+1) columns and 2^z rows of square tiles of `tile` pixels,
// from 180 W and 90 N, in east longitude; so a pixel of level z spans
// 360 / (2^(z+1) tile) degrees.

// box: [west, south, width, height] in degrees (west + width <= 180);
// ppd: screen pixels per degree wanted; live: { tile, levels }.
// Returns { z, c0, c1, r0, r1, box: the region the tiles cover }.
export function tilesFor(box, ppd, live, maxTiles = 48) {
  const [west, south, span, tall] = box;
  let z = Math.max(0, Math.min(live.levels, Math.ceil(Math.log2((ppd * 360) / (2 * live.tile)))));
  for (;;) {
    const cols = 2 ** (z + 1), rows = 2 ** z;
    const c0 = Math.max(0, Math.floor(((west + 180) / 360) * cols));
    const c1 = Math.min(cols - 1, Math.ceil(((west + span + 180) / 360) * cols) - 1);
    const r0 = Math.max(0, Math.floor(((90 - (south + tall)) / 180) * rows));
    const r1 = Math.min(rows - 1, Math.ceil(((90 - south) / 180) * rows) - 1);
    if ((c1 - c0 + 1) * (r1 - r0 + 1) <= maxTiles || z === 0) {
      return {
        z, c0, c1, r0, r1,
        box: [-180 + (c0 * 360) / cols, 90 - ((r1 + 1) * 180) / rows, ((c1 - c0 + 1) * 360) / cols, ((r1 - r0 + 1) * 180) / rows],
      };
    }
    z--;
  }
}
