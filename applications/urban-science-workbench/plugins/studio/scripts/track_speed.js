// How fast each thing moves: the speed along every track shown, by the great-circle distance between its samples, given back as the same tracks coloured by speed, with the fastest in a table
// and the distribution of speeds as a curve. An example of a script in JavaScript: it runs in a worker of the page, and is given and answers what a Python script is.
//
// Limits: a speed between two samples is a mean over the time between them, so that a track sampled sparsely reads slower than it went.

const EARTH_RADIUS_KM = 6371.0088;
const RAD = Math.PI / 180;
const distanceKm = (la0, lo0, la1, lo1) => {
  const a = Math.sin(((la1 - la0) * RAD) / 2) ** 2 + Math.cos(la0 * RAD) * Math.cos(la1 * RAD) * Math.sin(((lo1 - lo0) * RAD) / 2) ** 2;
  return 2 * EARTH_RADIUS_KM * Math.asin(Math.min(1, Math.sqrt(a)));
};

export function create(api) {
  return new class extends api.SpaceTimeAnalysis {
    id = 'track-speed';
    title = 'How fast each thing moves';
    description = 'The speed along every track shown, by the great-circle distance between its samples.';
    params = [api.param('fastest', 'Fastest listed', 8, { kind: 'integer', min: 1, max: 100 })];

    run(records, params, context) {
      const found = records.find((r) => r.kind === 'tracks');
      if (!found) throw new Error('This needs tracks: show a record of tracks first.');
      const all = [];
      const items = found.tracks.map((t, n) => {
        context.progress(n / found.tracks.length, t.name || String(t.id));
        const speed = t.time.map((_, k) => {
          const a = Math.max(0, k - 1), b = Math.min(t.time.length - 1, k + 1), hours = (t.time[b] - t.time[a]) / 3600000;
          return hours > 0 ? distanceKm(t.lat[a], t.lon[a], t.lat[b], t.lon[b]) / hours : null;
        });
        const known = speed.filter((v) => v !== null);
        all.push(...known);
        return { id: t.id, name: t.name, time: t.time, lat: t.lat, lon: t.lon, values: { speed }, mean: known.reduce((p, q) => p + q, 0) / Math.max(1, known.length) };
      });
      const ranked = [...items].sort((p, q) => q.mean - p.mean).slice(0, params.fastest);
      const top = Math.max(...all), bins = new Array(24).fill(0);
      for (const v of all) bins[Math.min(23, Math.floor((v / top) * 24))]++;
      return [
        api.tracks(`How fast: ${found.title}`, items, { units: { speed: 'km/h' }, show: { value: 'speed', map: 'plum' } }),
        api.table('The fastest on average', ['track', 'km/h'], ranked.map((t) => [t.name || t.id, Math.round(t.mean * 10) / 10])),
        api.series('Speeds at every sample', [{ name: 'samples', x: bins.map((_, k) => ((k + 0.5) / 24) * top), y: bins }], { xlabel: 'km/h', ylabel: 'samples' }),
      ];
    }
  }();
}
