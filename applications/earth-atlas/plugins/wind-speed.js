// The smallest plug-in: the speed of the wind, a Diagnostic with one output.
// It recomputes a quantity the atlas already shows, so that the two can be set
// side by side, and shows the whole contract in a page. It is not listed in
// index.json: load it from the panel's URL field as plugins/wind-speed.js, or
// from the console, atlas.plugins.load('plugins/wind-speed.js').

export function create(api) {
  class WindSpeed extends api.Diagnostic {
    // Identity: the product is named plugin-wind-speed, and the model "Wind speed (plug-in)".
    id = 'wind-speed';
    title = 'Wind speed';
    description = 'The speed of the wind chosen, from its two components: the smallest plug-in, to be read first.';

    // What it needs: a wind, at whatever level the user chooses; run() receives
    // inputs.wind = { product, names: [u, v], ... }.
    needs() { return [{ name: 'wind', label: 'Wind', standard_names: ['eastward_wind', 'northward_wind'] }]; }

    // What it returns: one variable, on the wind's grid, at each of its times.
    outputs() { return [{ name: 'speed', long_name: 'Wind speed', standard_name: 'wind_speed', units: 'm s-1' }]; }

    // The algorithm, for time `index`: Float32Arrays in, Float32Arrays out; NaN stays NaN.
    compute(inputs, params, index) {
      const { product, names: [u, v] } = inputs.wind;
      const a = product.slice(u, index), b = product.slice(v, index), speed = new Float32Array(a.length);
      for (let i = 0; i < a.length; i++) speed[i] = Math.hypot(a[i], b[i]);
      return { speed };
    }
  }
  return new WindSpeed();
}
