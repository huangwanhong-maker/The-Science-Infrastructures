// Paints the base map (src/render/relief.js) away from the page's main thread:
// at 8192 x 4096 pixels it takes about two seconds.
//
// Message in:  { layers: the geography product's land, lakes, glaciers and ice
//                shelves; coords: Float32Array; elevation: { meta, codes } or
//                null; width, height, style }
// Message out: { rgba: Uint8Array (width x height x 4) } or { error }

import { baseImage } from '../render/relief.js';
import { FieldProduct } from '../core/grid.js';

self.onmessage = (event) => {
  try {
    const { layers, coords, elevation, width, height, style } = event.data;
    const field = elevation ? new FieldProduct(elevation.meta, elevation.codes.buffer) : null;
    const { data } = baseImage({ elevation: field, geography: { layers }, coords, width, height, style });
    const rgba = new Uint8Array(data.buffer);
    self.postMessage({ rgba }, [rgba.buffer]);
  } catch (e) {
    self.postMessage({ error: e.message });
  }
};
