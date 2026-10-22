// The Earth's crust beneath a place, from CRUST1.0 (Laske et al. 2013): its
// layers from the surface to the Moho, with their P- and S-wave speeds and
// densities, and the mantle beneath.
//
// Provisional stand-in for core/domains/geophysics (laboratory ADR 0024).
//
// CRUST1.0 averages each 1-degree cell, so a place takes the values of the
// cell that contains it, not an interpolation between cells (the model's
// readme: "please use the midpoint of this cell"). A layer's base is the top
// of the next layer present; layers of zero thickness are absent.

export const LAYERS = ['water', 'ice', 'upper sediments', 'middle sediments', 'lower sediments',
  'upper crystalline crust', 'middle crystalline crust', 'lower crystalline crust', 'mantle'];

// The index of the cell containing (lon, lat) on the product's grid.
export function cellIndex(grid, lon, lat) {
  const r = Math.min(grid.nlat - 1, Math.max(0, Math.floor((grid.lat0 + grid.dlat / 2 - lat) / grid.dlat)));
  const c = (((Math.floor((lon - (grid.lon0 - grid.dlon / 2)) / grid.dlon)) % grid.nlon) + grid.nlon) % grid.nlon;
  return r * grid.nlon + c;
}

// The column at a cell: { layers: [{ name, top, bottom, thickness, vp, vs, rho }],
// moho (km below sea level), crust (km, solid surface to Moho), sediments (km),
// water (km), mantle: { vp, vs, rho }, oceanic }.
// get(quantity, layer) -> the value at the cell (tops in km above sea level).
export function column(get) {
  const tops = LAYERS.map((_, k) => get('top', k));
  const layers = [];
  for (let k = 0; k < 8; k++) {
    const bottom = tops[k + 1];
    const thickness = tops[k] - bottom;
    if (thickness <= 1e-6) continue;
    layers.push({ name: LAYERS[k], top: tops[k], bottom, thickness, vp: get('vp', k), vs: get('vs', k), rho: get('rho', k) });
  }
  const water = tops[0] - tops[1], sediments = tops[2] - tops[5];
  return {
    layers, moho: -tops[8], crust: tops[2] - tops[8], sediments, water,
    surface: tops[0], solid: tops[2],
    mantle: { vp: get('vp', 8), vs: get('vs', 8), rho: get('rho', 8) },
    // Oceanic crust: under the sea and thinner than 15 km of rock.
    oceanic: water > 0.5 && tops[2] - tops[8] < 15,
  };
}
