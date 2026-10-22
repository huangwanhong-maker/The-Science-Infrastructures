// Terrain and water from the relief: a Plugin. From the elevation chosen
// (ETOPO 2022), the slope of the surface and the area draining through each
// cell, by the geography domain's hydrology: depressions filled by
// Priority-Flood (Barnes, Lehman and Mulla 2014), flow to the steepest of the
// eight neighbours with distances on the sphere (D8), and the areas of the
// cells carried down to the sea. With a place chosen on the map, its
// watershed (every cell that drains to it), outlined, with the main stem above
// it and the path of its water down to the sea, as a sketch, and a table.
//
// A source finer than the grid asked for is averaged onto it in blocks, so
// that the 5′ relief (9.3 million cells) is taken at 10′ (2.3 million) unless
// asked otherwise. Caveats: on a grid of 10′ (18 km) a river is a chain of
// cells, a valley narrower than a cell may be missed (a river may leave its
// basin through a pass the cells cannot see), and the floors of lakes and dry
// depressions are filled (theory/drainage-from-elevation.md).

import {
  slope, fillDepressions, flowDirections, flowAccumulation, watershed, downstream, mainStem, snapToStream,
  cellOutline, pathOf, cellOf, centreOf,
} from '../../../domains/geography/js/index.js';

export function create(api) {
  class Terrain extends api.Plugin {
    id = 'terrain';
    title = 'Terrain and water';
    domain = 'geography';
    description = 'Slope and the area draining through each cell, from the relief; with a place chosen, its watershed, the main stem above it and its way to the sea.';
    parameters = {
      cells: { label: 'Grid', default: 10, choices: [[10, '10′ (2.3 million cells)'], [20, '20′ (0.6 million cells)'], [0, 'the source’s own']] },
      snap: { label: 'Take the place onto a river within', default: 1, min: 0, max: 5, step: 1, unit: 'cells' },
    };

    needs() { return [{ name: 'relief', label: 'Relief', standard_name: 'altitude' }]; }

    // The elevation on a grid of `minutes` (0: the source's), averaged in blocks from a finer one.
    elevation(input, minutes) {
      const { product, names: [name] } = input, g = product.grid, z = product.slice(name, 0);
      const k = minutes ? Math.max(1, Math.round(minutes / (g.dlon * 60))) : 1;
      if (k === 1 || g.nlon % k || g.nlat % k) return { z, grid: g };
      const grid = { lon0: g.lon0 + ((k - 1) * g.dlon) / 2, dlon: g.dlon * k, nlon: g.nlon / k, lat0: g.lat0 - ((k - 1) * g.dlat) / 2, dlat: g.dlat * k, nlat: g.nlat / k, periodic: g.periodic };
      const out = new Float32Array(grid.nlon * grid.nlat);
      for (let r = 0; r < grid.nlat; r++) {
        for (let c = 0; c < grid.nlon; c++) {
          let s = 0, m = 0;
          for (let a = 0; a < k; a++) for (let b = 0; b < k; b++) { const v = z[(r * k + a) * g.nlon + c * k + b]; if (v === v) { s += v; m++; } }
          out[r * grid.nlon + c] = m ? s / m : NaN;
        }
      }
      return { z: out, grid };
    }

    async run(inputs, params, context) {
      const time = {}, mark = (what, t) => { time[what] = Date.now() - t; return Date.now(); };
      let t = Date.now();
      const { z, grid } = this.elevation(inputs.relief, Number(params.cells ?? 10));
      t = mark('grid', t);
      await this.breathe(context, 0.1, 'slope');
      const steep = slope(z, grid);
      t = mark('slope', t);
      await this.breathe(context, 0.25, 'filling depressions');
      const { filled, outlet, raised } = fillDepressions(z, grid);
      t = mark('fill', t);
      await this.breathe(context, 0.6, 'flow');
      const dirs = flowDirections(filled, grid, outlet), acc = flowAccumulation(dirs, grid);
      t = mark('flow', t);
      const upstream = new Float32Array(acc.length);
      for (let i = 0; i < acc.length; i++) upstream[i] = outlet[i] ? NaN : acc[i];
      const cells = grid.nlon * grid.nlat, total = Object.values(time).reduce((a, b) => a + b, 0);
      const made = api.makeProduct({
        name: this.productName, title: this.title, model: `${this.title} (plug-in)`, domain: this.domain, grid, times: [],
        variables: [
          { name: 'slope', long_name: 'Slope of the surface (rise over run)', units: '1', slices: [steep] },
          { name: 'upstream_area', long_name: 'Area draining through the cell (flow accumulation)', units: 'km2', slices: [upstream] },
        ],
        source: { dataset: inputs.relief.meta.name, provider: `computed in the page from ${inputs.relief.model}`, credit: `${this.title} (plug-in), from ${inputs.relief.model}` },
      });
      const minutes = Math.round(grid.dlon * 60);
      const results = [api.field(made, { surface: 'upstream_area' })];
      const timing = `${(cells / 1e6).toFixed(2)} million cells of ${minutes}′: grid ${time.grid} ms, slope ${time.slope} ms, filling ${time.fill} ms (${raised.toLocaleString('en-GB')} cells raised), flow ${time.flow} ms; ${(total / 1000).toFixed(1)} s in all`;
      if (!context.pinned) {
        results.push(api.note(`${timing}. Choose a place on the map and run again for its watershed.`));
        return results;
      }
      // The place's watershed, its main stem and its way to the sea.
      const [lon, lat] = context.pinned, chosen = cellOf(grid, lon, lat);
      if (chosen < 0 || outlet[chosen]) {
        results.push(api.note(`${timing}. The place chosen is at sea (or below sea level by a sea), where water leaves the land: choose a place on land.`));
        return results;
      }
      const cell = snapToStream(upstream, grid, chosen, Math.round(Number(params.snap ?? 1)));   // onto land only: NaN at sea
      const shed = watershed(dirs, grid, cell), rings = cellOutline(grid, shed);
      const down = pathOf(grid, downstream(dirs, grid, cell)), stem = pathOf(grid, mainStem(dirs, acc, grid, cell));
      const [plon, plat] = centreOf(grid, cell), end = down.coords.length / 2 - 1;
      const sea = [down.coords[2 * end], down.coords[2 * end + 1]];
      const deg = (x, pos, neg) => `${Math.abs(x).toFixed(2)}° ${x >= 0 ? pos : neg}`;
      const km = (m) => `${Math.round(m / 1000).toLocaleString('en-GB')} km`;
      results.push(api.sketch({
        title: `Watershed of ${deg(plat, 'N', 'S')}, ${deg(plon, 'E', 'W')}: ${Math.round(acc[cell]).toLocaleString('en-GB')} km²`,
        polylines: [
          ...rings.map((coords) => ({ coords, colour: '#e8c3dc', width: 1.6, alpha: 0.95 })),
          { coords: Array.from(stem.coords), colour: '#a9cdf0', width: 1.6, alpha: 0.9 },
          { coords: Array.from(down.coords), colour: '#8fb8e8', width: 2.2, alpha: 0.95 },
        ],
        markers: [{ lon: plon, lat: plat, colour: '#f6c47e', size: 4 }, { lon: sea[0], lat: sea[1], colour: '#8fb8e8', size: 3 }],
      }));
      results.push(api.table(`The watershed of the place, on ${inputs.relief.model} at ${minutes}′`, ['', ''], [
        ['Area draining to the place', `${Math.round(acc[cell]).toLocaleString('en-GB')} km²`],
        ['Cells', shed.length.toLocaleString('en-GB')],
        ['The place, on the river', `${deg(plat, 'N', 'S')}, ${deg(plon, 'E', 'W')}${cell !== chosen ? ' (moved onto the river)' : ''}`],
        ['Elevation of its cell', `${Math.round(z[cell])} m`],
        ['Main stem above the place', km(stem.length)],
        ['Down to the sea', km(down.length)],
        ['Reaching the sea at', `${deg(sea[1], 'N', 'S')}, ${deg(sea[0], 'E', 'W')}`],
      ], timing));
      return results;
    }
  }
  return new Terrain();
}
