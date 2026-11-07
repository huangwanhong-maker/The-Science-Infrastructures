// The Earth by night, whatever the hour: NASA's Black Marble, which the globe already shows on the side the Sun has left, here over
// all of it. Nothing is computed; the Sun is only taken from the ground for as long as this is on.

import { card } from '../shared.js';

export async function layer(ctx) {
  const { earth, dom: { h } } = ctx;
  earth.setLook({ nightAll: true });
  return {
    legend: card(h, 'The Earth by night', { note: 'Annual composites of night-time light, whatever the clock says.', source: 'NASA Black Marble (VIIRS), through GIBS' }),
    dispose() { earth.setLook({ nightAll: false }); },
  };
}
