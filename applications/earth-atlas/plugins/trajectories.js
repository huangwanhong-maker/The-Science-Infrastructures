// Air parcels carried through the wind: a Trajectories plug-in. With a place
// chosen on the map, the parcels start on a ring around it; otherwise on a
// lattice over the globe. The base class advances each through the wind's own
// times by a Runge–Kutta step on the sphere and returns a sketch: the paths,
// and a point on each that moves with the time shown. This is the Lagrangian
// view of a flow, where the particles of the Motion layer are its Eulerian one.

import { move } from '../../../domains/geodesy/js/index.js';

export function create(api) {
  class AirParcels extends api.Trajectories {
    id = 'trajectories';
    title = 'Air parcels';
    description = 'Parcels carried through the wind chosen, from a ring around the chosen place, or from a lattice over the globe when no place is chosen.';
    parameters = {
      hours: { label: 'Hours ahead', default: 48, min: 6, max: 240, step: 6, unit: 'h' },
      step: { label: 'Time step', default: 1, min: 0.25, max: 6, step: 0.25, unit: 'h' },
      count: { label: 'Parcels', default: 24, min: 1, max: 200, step: 1 },
      radius: { label: 'Ring radius', default: 300, min: 0, max: 3000, step: 10, unit: 'km' },
    };

    needs() { return [{ name: 'wind', label: 'Wind', standard_names: ['eastward_wind', 'northward_wind'], times: 2 }]; }

    async seeds(inputs, params, context) {
      const n = Math.max(1, Math.round(params.count));
      if (!context.pinned) return api.lattice(n);
      const [lon, lat] = context.pinned, radius = params.radius * 1000;
      return Array.from({ length: n }, (_, k) => {
        const azimuth = (2 * Math.PI * k) / n;
        const [plon, plat] = move(lon, lat, radius * Math.sin(azimuth), radius * Math.cos(azimuth));
        return { lon: plon, lat: plat };
      });
    }
  }
  return new AirParcels();
}
