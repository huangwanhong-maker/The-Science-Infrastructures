// Records of the data library, as records in space and time. Names of regions in the catalogues are not carried, since they are the catalogues' own and not under the
// data library's policy of names (ADR 0069): a hover says where by the map itself.

const kept = new Map();
const json = (file) => {
  if (!kept.has(file)) kept.set(file, fetch(`/data/products/${file}`).then(async (r) => { if (!r.ok) throw new Error(`The data library has no ${file} yet: its recipe fetches it.`); return r.json(); }).catch((e) => { kept.delete(file); throw e; }));
  return kept.get(file);
};
const columnsOf = (d) => Object.fromEntries(d.fields.map((f, k) => [f, k]));

// A small hexagon about a point, for a place that has a position and no shape of its own.
const hexAbout = (lat, lon, radius) => Array.from({ length: 6 }, (_, k) => {
  const a = (k * Math.PI) / 3;
  return [lon + ((radius * Math.cos(a)) / 111000) / Math.cos((lat * Math.PI) / 180), lat + (radius * Math.sin(a)) / 111000];
});

export const SOURCES = {
  'walkers-week': {
    title: "Melbourne's walkers through a mean week", credit: 'City of Melbourne, pedestrian counting system', style: { value: 'walkers', map: 'plum', scale: 'linear', top: 340, zone: 'Australia/Melbourne', weekly: true }, kind: 'cells',
    async read() {
      const r = await fetch('/api/street-oceanography/profiles', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: '{}' });
      const answer = await r.json();
      if (!r.ok) throw new Error(answer.message || 'The counters could not be read.');
      const sensors = answer.filter((s) => s.profile?.length === 168);
      // A mean week has no date of its own: it is laid on a nominal one, a Monday, at Melbourne's midnight in summer time (UTC+11).
      const monday = Date.UTC(2026, 2, 1, 13);
      return {
        kind: 'cells', title: this.title, units: { walkers: 'an hour' },
        note: 'Each counter\'s mean count at each hour of the week: its Monday is a nominal date.',
        time: Array.from({ length: 168 }, (_, h) => monday + h * 3600000),
        cells: sensors.map((s) => ({ id: String(s.id), name: s.name, ring: hexAbout(s.lat, s.lon, 34) })),
        values: { walkers: Array.from({ length: 168 }, (_, h) => sensors.map((s) => s.profile[h])) },
      };
    },
  },
  gcmt: {
    title: 'Earthquakes since 1976', credit: 'Global CMT Project (Ekström et al. 2012)', style: { value: 'magnitude', map: 'plum', scale: 'linear', pixel: 6 }, kind: 'events',
    async read() {
      const d = await json('events-gcmt.json'), c = columnsOf(d), e = d.events;
      return { kind: 'events', title: this.title, time: e.map((x) => x[c.time_ms]), lat: e.map((x) => x[c.lat]), lon: e.map((x) => x[c.lon]), values: { magnitude: e.map((x) => x[c.mw]), depth: e.map((x) => x[c.depth_km]) }, units: { magnitude: 'Mw', depth: 'km' } };
    },
  },
  usgs: {
    title: 'Earthquakes this month', credit: 'USGS ComCat', style: { value: 'magnitude', map: 'plum', pixel: 8 }, kind: 'events',
    async read() {
      const d = await json('events-usgs-earthquakes-month.json'), c = columnsOf(d), e = d.events;
      return { kind: 'events', title: this.title, time: e.map((x) => x[c.time_ms]), lat: e.map((x) => x[c.lat]), lon: e.map((x) => x[c.lon]), values: { magnitude: e.map((x) => x[c.mag]), depth: e.map((x) => x[c.depth_km]) }, units: { magnitude: 'M', depth: 'km' } };
    },
  },
  ibtracs: {
    title: 'Tropical cyclones, the last three years', credit: 'IBTrACS (NOAA NCEI)', style: { value: 'wind', map: 'plum', window_ms: 30 * 86400000 }, kind: 'tracks',
    async read() {
      const d = await json('tracks-ibtracs-last3years.json'), c = columnsOf(d);
      return { kind: 'tracks', title: this.title, units: { wind: 'm/s', pressure: 'hPa' }, tracks: d.tracks.map((t) => ({ id: t.id, name: t.name, time: t.points.map((p) => p[c.time_ms]), lat: t.points.map((p) => p[c.lat]), lon: t.points.map((p) => p[c.lon]), values: { wind: t.points.map((p) => p[c.wind_ms]), pressure: t.points.map((p) => p[c.pressure_hpa]) } })) };
    },
  },
};
