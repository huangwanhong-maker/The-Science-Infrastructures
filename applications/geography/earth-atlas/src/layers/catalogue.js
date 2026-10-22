// What can be drawn: layers built from the products in the library.
//
// Layers are typed by the shape of their data, not by subject:
//   surface  a quantity over the whole surface: a field (one variable, or the
//            speed of a vector field), a value per country (an indicator), or
//            a class per cell (climate classes, by period and scenario);
//   flow     a vector field shown by moving particles;
//   lines    isolines of a field, with its highs and lows;
//   events   things that happened at places and instants;
//   tracks   tropical cyclones: found here in each model's sea-level pressure,
//            or observed and forecast officially (a tracks product).
// A field variable is recognised by its CF standard name and its level, so
// that any product carrying, say, geopotential height at 500 hPa offers the
// same layer, drawn the same way. Products of different models for the same
// quantity share one *key*: the panel offers them as a choice of model, and
// two of them can be compared as a difference. A variable that is not
// recognised is still offered, with a range taken from its own values: new
// data, such as our own model output, can be seen as soon as it is a product.
//
// A climatology (monthly normals) is offered as surfaces that follow the month
// of the time shown, and feeds the place card's climate diagram.
//
// Styles (colour map, range, scale) live here, in the application, and not in
// the data. Ranges are in the units shown (src/core/units.js); `ranges` gives
// a range per level where the quantity's values depend on the level.

const QUANTITIES = [
  { id: 'air-temperature', standard: 'air_temperature', label: 'Air temperature', colormap: 'thermal', range: [-35, 40], ranges: { '850 hPa': [-35, 30], '500 hPa': [-50, 0] } },
  { id: 'surface-temperature', standard: 'surface_temperature', label: 'Surface temperature', colormap: 'thermal', range: [-40, 50] },
  { id: 'pressure', standard: 'air_pressure_at_sea_level', label: 'Pressure at sea level', colormap: 'dense', reverse: true, range: [975, 1040] },
  { id: 'height', standard: 'geopotential_height', label: 'Geopotential height', colormap: 'haline', range: [4900, 5950], ranges: { '850 hPa': [1150, 1600], '500 hPa': [4900, 5950], '250 hPa': [9600, 11000] } },
  { id: 'precipitation', standard: 'precipitation_flux', label: 'Precipitation', colormap: 'haline', range: [0.1, 20], scale: 'log', fade: 0.2 },
  { id: 'humidity', standard: 'relative_humidity', label: 'Relative humidity', colormap: 'tempo', reverse: true, range: [0, 100] },
  { id: 'cloud', standard: 'cloud_area_fraction', label: 'Cloud cover', colormap: 'ice', range: [0, 100], fade: 0.12 },
  { id: 'water-vapour', standard: 'atmosphere_mass_content_of_water_vapor', label: 'Precipitable water', colormap: 'rain', reverse: true, range: [0, 70] },
  { id: 'cape', standard: 'atmosphere_convective_available_potential_energy_wrt_surface', label: 'Convective energy (CAPE)', colormap: 'thermal', range: [0, 3500], scale: 'sqrt', fade: 0.12 },
  { id: 'vorticity', standard: 'atmosphere_absolute_vorticity', label: 'Absolute vorticity', colormap: 'balance', range: [-30, 30] },
  { id: 'relative-vorticity', standard: 'atmosphere_relative_vorticity', label: 'Relative vorticity', colormap: 'balance', range: [-20, 20] },
  { id: 'shortwave-down', standard: 'surface_downwelling_shortwave_flux', label: 'Sunlight reaching the surface', colormap: 'thermal', range: [0, 1000] },
  { id: 'shortwave-up', standard: 'surface_upwelling_shortwave_flux', label: 'Sunlight reflected by the surface', colormap: 'thermal', range: [0, 500] },
  { id: 'longwave-down', standard: 'surface_downwelling_longwave_flux', label: 'Infrared from the sky', colormap: 'thermal', range: [100, 450] },
  { id: 'longwave-up', standard: 'surface_upwelling_longwave_flux', label: 'Infrared from the surface', colormap: 'thermal', range: [150, 550] },
  { id: 'sst', standard: 'sea_surface_temperature', label: 'Sea-surface temperature', colormap: 'thermal', range: [-2, 32] },
  { id: 'sst-anomaly', standard: 'surface_temperature_anomaly', label: 'Sea-surface temperature anomaly', colormap: 'balance', range: [-4, 4] },
  { id: 'sea-ice', match: (v) => v.name === 'ice' || /sea ice concentration/i.test(v.long_name || ''), label: 'Sea-ice concentration', colormap: 'ice', range: [0, 1], fade: 0.12, percent: true },
  { id: 'wave-height', standard: 'sea_surface_wave_significant_height', label: 'Significant wave height', colormap: 'deep', reverse: true, range: [0, 9] },
  { id: 'swell-height', standard: 'sea_surface_swell_wave_significant_height', label: 'Swell height', colormap: 'deep', reverse: true, range: [0, 7] },
  { id: 'wind-sea-height', standard: 'sea_surface_wind_wave_significant_height', label: 'Wind-sea height', colormap: 'deep', reverse: true, range: [0, 7] },
  { id: 'wave-period', standard: 'sea_surface_wave_period_at_variance_spectral_density_maximum', label: 'Peak wave period', colormap: 'speed', reverse: true, range: [2, 18] },
  { id: 'swell-period', standard: 'sea_surface_swell_wave_period', label: 'Swell period', colormap: 'speed', reverse: true, range: [2, 18] },
  { id: 'wind-sea-period', standard: 'sea_surface_wind_wave_period', label: 'Wind-sea period', colormap: 'speed', reverse: true, range: [1, 12] },
  { id: 'elevation', standard: 'altitude', label: 'Elevation', base: true },
  // The crust (CRUST1.0): maps only; the place card draws the whole column.
  { id: 'moho', match: (v) => v.name === 'moho_depth', label: 'Depth of the Moho', colormap: 'deep', range: [0, 70], noProbe: true },
  { id: 'crust-thickness', match: (v) => v.name === 'crust_thickness', label: 'Thickness of the crust', colormap: 'matter', reverse: true, range: [0, 70], noProbe: true },
  { id: 'sediment-thickness', match: (v) => v.name === 'sediment_thickness', label: 'Thickness of sediments', colormap: 'speed', reverse: true, range: [0, 12], noProbe: true },
];

// Monthly normals of a climatology, recognised by standard name and cell methods.
const NORMALS = [
  { id: 'tmax', match: (v) => v.standard_name === 'air_temperature' && /maximum within days/.test(v.cell_methods || ''), label: 'Daily maximum temperature', colormap: 'thermal', range: [-35, 45] },
  { id: 'tmin', match: (v) => v.standard_name === 'air_temperature' && /minimum within days/.test(v.cell_methods || ''), label: 'Daily minimum temperature', colormap: 'thermal', range: [-45, 30] },
  { id: 'ppt', match: (v) => v.standard_name === 'precipitation_amount', label: 'Precipitation in the month', colormap: 'haline', range: [2, 500], scale: 'log', fade: 0.15 },
];

// Names of series of class maps.
const CLASS_SERIES = { 'koppen-geiger': 'Climate classes (Köppen–Geiger)' };

// Vector quantities. `toward`: named by where they flow (currents, plates), not
// where they come from (wind). `pace` scales how fast particles move on screen,
// and `siPerUnit` converts the stored units to m/s, so that the legend can say
// how much faster than real the particles go.
const VECTORS = [
  { id: 'wind', standard: 'wind', label: 'Wind', colormap: 'orchid', range: [0, 25], ranges: { '850 hPa': [0, 35], '500 hPa': [0, 50], '250 hPa': [0, 80] } },
  { id: 'current', standard: 'sea_water_velocity', label: 'Ocean current', colormap: 'orchid', range: [0, 1.5], toward: true },
  { id: 'plate', standard: 'plate_velocity', label: 'Plate motion', colormap: 'matter', reverse: true, range: [0, 100], toward: true,
    pace: 0.12, siPerUnit: 1 / (1000 * 365.25 * 86400) },
];

// Short names and units for World Bank indicators, by code.
export const INDICATORS = {
  'SP.POP.TOTL': ['Population', ''],
  'EN.POP.DNST': ['Population density', 'per km²'],
  'NY.GDP.PCAP.CD': ['GDP per person', 'US$'],
  'SP.DYN.LE00.IN': ['Life expectancy', 'years'],
  'SP.DYN.TFRT.IN': ['Births per woman', ''],
  'SP.URB.TOTL.IN.ZS': ['Urban population', '%'],
  'AG.LND.FRST.ZS': ['Forest area', '% of land'],
  'EG.ELC.ACCS.ZS': ['Access to electricity', '%'],
  'EN.GHG.CO2.PC.CE.AR5': ['CO₂ emitted per person', 't'],
  'EN.GHG.ALL.PC.CE.AR5': ['Greenhouse gases per person', 't CO₂e'],
  'SI.POV.GINI': ['Gini index of income', ''],
};

export const DOMAINS = {
  meteorology: 'Meteorology', oceanography: 'Oceanography', climate: 'Climate',
  geography: 'Geography', geophysics: 'Geophysics',
};

// Models in the order they are offered.
const MODEL_ORDER = ['GFS', 'AIGFS', 'AIGEFS mean', 'Pangu-Weather (laboratory)', 'Pangu-Weather hindcast (laboratory)',
  'Shallow water (laboratory)', 'Barotropic (laboratory)'];
export const modelRank = (m) => { const i = MODEL_ORDER.indexOf(m); return i < 0 ? MODEL_ORDER.length : i; };

const isDirection = (v) => /direction/.test(v.standard_name || '') || v.units === 'degree';
const isSpread = (v) => /realization:\s*standard_deviation/.test(v.cell_methods || '');
const isMean = (v) => /realization:\s*mean/.test(v.cell_methods || '');

// "air temperature at 2m" -> " at 2 m"; a stated level takes precedence.
function height(v) {
  if (v.level) return ` at ${v.level}`;
  const m = /\bat\s*(\d+(?:\.\d+)?)\s*m\b/i.exec(v.long_name || '');
  return m ? ` at ${m[1]} m` : '';
}

const levelOf = (v) => v.level || height(v).replace(/^ at /, '');
const slug = (s) => (s ? `.${s.replace(/\s+/g, '')}` : '');
const capital = (s) => s.charAt(0).toUpperCase() + s.slice(1);

function source(meta) {
  const s = meta.source || {};
  return { dataset: s.dataset, profile: s.profile, provider: s.provider, licence: s.licence, citation: s.citation, retrieved: s.retrieved_at, credit: s.credit };
}

const periodText = (p) => `${p[0]}–${p[1]}`;

// Class maps of one series: the present first, then by period and scenario.
function classSurfaces(library) {
  const out = [];
  for (const [name, meta] of library.classMeta) {
    const period = periodText(meta.period);
    out.push({ id: `classes:${name}`, key: `classes:${meta.series}`, kind: 'classes', label: CLASS_SERIES[meta.series] || meta.title,
      model: meta.computed ? `${period}, ${meta.model}` : meta.scenario ? `${period}, ${meta.scenario}` : period,
      period: meta.period, scenario: meta.scenario || null, computed: Boolean(meta.computed), product: name,
      domain: meta.domain || 'climate', sourceName: meta.computed ? meta.model : `the map of ${(meta.source && meta.source.credit || '').split(',').slice(1, 2).join('').trim() || meta.title}` });
  }
  const present = (d) => (d.period[1] === 2020 && !d.computed ? 0 : 1);
  return out.sort((a, b) => present(a) - present(b) || a.period[0] - b.period[0] || (a.scenario || '').localeCompare(b.scenario || '') || a.computed - b.computed);
}

export function buildCatalogue(library) {
  const surfaces = [], flows = [], lines = [], events = [], probes = [], clouds = [], normals = [], tracks = [];
  surfaces.push(...classSurfaces(library));
  for (const [name, meta] of library.fieldMeta) {
    const domain = meta.domain || 'geography';
    const model = meta.model || meta.title;
    if (meta.climatology) {
      const found = {};
      for (const v of meta.variables) {
        const known = NORMALS.find((q) => q.match(v));
        if (!known) continue;
        found[known.id] = v.name;
        surfaces.push({ id: `normal:${known.id}@${name}`, key: `normal:${known.id}`, kind: 'field', label: `${known.label}, normal`, model, product: name,
          variables: [v.name], domain: 'climate', units: v.units, level: '', climatology: meta.climatology, source: source(meta),
          style: { colormap: known.colormap, range: known.range, scale: known.scale || 'linear', fade: known.fade } });
      }
      if (found.tmax && found.tmin && found.ppt) normals.push({ product: name, model, period: meta.climatology.period, variables: found });
      continue;
    }
    const taken = new Set();
    for (const vec of meta.vectors || []) {
      const known = VECTORS.find((k) => k.standard === vec.standard_name) || { id: vec.name, label: capital(vec.name.replace(/_/g, ' ')), colormap: 'orchid', range: [0, 20] };
      const u = meta.variables.find((v) => v.name === vec.u);
      if (isSpread(u)) continue;
      const level = levelOf(u);
      const id = `${known.id}${slug(level)}@${name}`;
      const label = `${known.label}${height(u)}`;
      const range = (known.ranges && known.ranges[level]) || known.range;
      flows.push({ id, key: `${known.id}|${level}`, kind: 'vector', label, model, product: name, variables: [vec.u, vec.v], domain, units: u.units, source: source(meta),
        pace: known.pace || 1, siPerUnit: known.siPerUnit || 1 });
      surfaces.push({ id: `${id}:speed`, key: `speed:${known.id}|${level}`, kind: 'field', label: `${known.label} speed${height(u)}`, model, product: name, variables: [vec.u, vec.v], domain, units: u.units, level,
        style: { colormap: known.colormap, reverse: known.reverse, range, scale: 'linear' }, source: source(meta) });
      taken.add(vec.u); taken.add(vec.v);
      // The plates' own rows in the place card describe their motion.
      if (known.id !== 'plate') {
        probes.push({ kind: 'vector', label: `${known.label}${height(u)}`, model, product: name, variables: [vec.u, vec.v], units: u.units, domain, key: `${known.id}|${level}`, toward: Boolean(known.toward) });
      }
    }
    const dir = meta.variables.find((v) => v.standard_name === 'sea_surface_wave_from_direction_at_variance_spectral_density_maximum');
    const per = meta.variables.find((v) => v.standard_name === 'sea_surface_wave_period_at_variance_spectral_density_maximum');
    if (dir && per) {
      flows.push({ id: `waves@${name}`, key: 'waves|', kind: 'waves', label: 'Wave energy', model, product: name, variables: [dir.name, per.name], domain, units: 'm s-1', source: source(meta) });
    }
    for (const v of meta.variables) {
      if (taken.has(v.name)) continue;
      const known = QUANTITIES.find((q) => (q.standard && q.standard === v.standard_name) || (q.match && q.match(v)));
      const level = levelOf(v);
      if (isDirection(v)) {
        probes.push({ kind: 'direction', label: capital(v.long_name || v.name), model, product: name, variables: [v.name], units: v.units, domain });
        continue;
      }
      const spread = isSpread(v);
      const qid = known ? known.id : v.name;
      const label = known ? `${known.label}${height(v)}${spread ? ', ensemble spread' : ''}` : capital(v.long_name || v.name);
      const shownModel = isMean(v) ? `${model} mean` : model;
      if (!(known && known.noProbe)) probes.push({ kind: 'scalar', id: known && known.id, label, model: shownModel, product: name, variables: [v.name], units: v.units, domain, percent: known && known.percent, spread, key: `${qid}|${level}|${spread ? 'spread' : ''}` });
      if (known && known.base) continue;
      if (v.standard_name === 'cloud_area_fraction' && !spread) clouds.push({ product: name, variable: v.name, model });
      let style;
      if (spread) {
        style = { colormap: 'orchid', range: null, scale: 'linear', auto: [0, v.p98] };
      } else if (known) {
        style = { colormap: known.colormap, reverse: known.reverse, range: (known.ranges && known.ranges[level]) || known.range, scale: known.scale || 'linear', fade: known.fade };
      } else {
        style = { colormap: 'matter', reverse: true, clip: [0, 0.9], range: null, scale: 'linear', auto: [v.p02, v.p98] };
      }
      surfaces.push({ id: `${qid}${slug(level)}${spread ? '-spread' : ''}@${name}`, key: `${qid}|${level}|${spread ? 'spread' : ''}`, kind: 'field', label, model: shownModel, product: name, variables: [v.name], domain, units: v.units, level, style, source: source(meta), percent: known && known.percent, spread });
      if (spread) continue;
      if (v.standard_name === 'air_pressure_at_sea_level') {
        if (meta.times && meta.times.length > 4) {
          tracks.push({ id: `cyclones@${name}`, kind: 'model', label: `Found in ${shownModel}`, model: shownModel, product: name, variable: v.name });
        }
        lines.push({ id: `isobars@${name}`, key: 'isobars|', label: 'Isobars', model: shownModel, product: name, variable: v.name, step: 400, base: 101300, extrema: { radius: 10, prominence: 500, polarRows: 8 }, units: v.units, labelUnits: 'hPa' });
      }
      if (v.standard_name === 'air_temperature') {
        lines.push({ id: `isotherms${slug(level)}@${name}`, key: `isotherms|${level}`, label: `Isotherms${height(v)}`, model: shownModel, product: name, variable: v.name, step: 5, base: 273.15, units: v.units, labelUnits: '°C' });
      }
      if (v.standard_name === 'geopotential_height') {
        const step = /250/.test(level) ? 120 : 60;
        lines.push({ id: `heights${slug(level)}@${name}`, key: `heights|${level}`, label: `Height contours${height(v)}`, model: shownModel, product: name, variable: v.name, step, base: 0, extrema: { radius: 10, prominence: 40, polarRows: 8, minAbsLat: 22 }, units: v.units, labelUnits: 'm' });
      }
    }
  }
  for (const e of library.list('indicators')) {
    const [label, unit] = INDICATORS[e.code] || [e.quantity || e.title, e.unit];
    surfaces.push({ id: `indicator:${e.name}`, key: `indicator:${e.name}`, kind: 'indicator', label, product: e.name, code: e.code, domain: e.domain || 'geography',
      units: unit, title: e.title, style: { colormap: 'matter', reverse: true, clip: [0.04, 0.96] }, source: e.source });
  }
  for (const e of library.list('tracks')) {
    tracks.push({ id: `tracks:${e.name}`, kind: 'official', label: 'Observed and official forecast', model: (e.source && e.source.credit) || e.title, product: e.name, source: e.source });
  }
  for (const e of library.list('events')) {
    // Mechanisms have a layer of their own (src/layers/mechanisms.js).
    if (e.kind === 'mechanism') continue;
    events.push({ id: e.name, label: e.title, kind: e.kind, product: e.name, count: e.count, window: e.window, domain: e.domain });
  }
  const order = ['meteorology', 'oceanography', 'climate', 'geography', 'geophysics'];
  const levelRank = (level) => { const m = /(\d+)\s*hPa/.exec(level || ''); return m ? 2000 - Number(m[1]) : 0; };
  const within = (s) => (s.key.startsWith('speed:wind') ? 0 : QUANTITIES.findIndex((q) => s.key.startsWith(`${q.id}|`)) + 1 || 500);
  const rank = (s) => (order.indexOf(s.domain) + 1 || 99) * 1e6 + within(s) * 1e4 + levelRank(s.level) + (s.key.endsWith('|spread') ? 5000 : 0);
  surfaces.sort((a, b) => rank(a) - rank(b) || modelRank(a.model) - modelRank(b.model));
  const flowRank = (f) => (f.key.startsWith('wind|') ? levelRank(f.key.split('|')[1]) : 5000);
  flows.sort((a, b) => flowRank(a) - flowRank(b) || modelRank(a.model) - modelRank(b.model));
  lines.sort((a, b) => a.key.localeCompare(b.key) || modelRank(a.model) - modelRank(b.model));
  clouds.sort((a, b) => modelRank(a.model) - modelRank(b.model));
  tracks.sort((a, b) => (a.kind === 'official' ? -1 : 0) - (b.kind === 'official' ? -1 : 0) || modelRank(a.model) - modelRank(b.model));
  return { surfaces, flows, lines, events, probes, clouds, normals, tracks };
}

// Definitions grouped by key, in catalogue order: [[key, [defs by model]]].
export function byKey(defs) {
  const groups = new Map();
  for (const d of defs) {
    if (!groups.has(d.key)) groups.set(d.key, []);
    groups.get(d.key).push(d);
  }
  return [...groups];
}
