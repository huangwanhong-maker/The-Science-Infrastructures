// Imagery for the globes: map services whose tiles are cut to the globe's pyramid.
//
// A service is asked for exactly the part of the body a tile covers, 256 pixels square, in
// geographic coordinates: a WMS by its bounding box (OGC 06-042; version 1.3.0 orders an
// EPSG:4326 box latitude first), a WMTS whose matrix is the pyramid itself (NASA Trek) by its
// row and column. Each service states its body, its role (a base, the night's lights, a field of
// data laid over the base, or a reference of lines and names), its deepest useful level, whether
// it changes from day to day, and the credit it requires; the page shows the credits of what it
// draws. The services named here allow browsers to fetch from them (CORS) and their terms allow
// this use:
//
//   EOX Sentinel-2 cloudless 2016   Earth, ~10 m, CC BY 4.0 (EOX IT Services GmbH)
//   NASA GIBS                       Earth: the Blue Marble, VIIRS true colour of a day (250 m), the
//                                   Black Marble at night; the sea's temperature of a day (GHRSST
//                                   MUR, 1 km) and its anomaly, rain of a day (GPM IMERG), and
//                                   coastlines, borders and place names; NASA's data, without
//                                   restriction (the references from OpenStreetMap, ODbL)
//   NASA Moon Trek                  the Moon: LRO's wide-angle camera mosaic at 303 px/°
//
// A service that changes daily is asked for a date (YYYY-MM-DD) or, without one, for its latest.
//
// The table is the default; an application may add services of its own of the same form.

export const IMAGERY = {
  'eox-s2cloudless-2016': {
    body: 'earth', title: 'Sentinel-2 cloudless 2016 (EOX)', kind: 'wmts', role: 'base',
    template: 'https://tiles.maps.eox.at/wmts/1.0.0/s2cloudless/default/WGS84/{z}/{row}/{col}.jpg', maxLevel: 13,
    credit: 'Sentinel-2 cloudless – https://s2maps.eu by EOX IT Services GmbH (Contains modified Copernicus Sentinel data 2016), CC BY 4.0',
  },
  'gibs-bluemarble': {
    body: 'earth', title: 'Blue Marble, relief and bathymetry (NASA)', kind: 'wms', version: '1.3.0', role: 'base',
    url: 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi', layer: 'BlueMarble_ShadedRelief_Bathymetry', format: 'image/jpeg', maxLevel: 7,
    credit: 'Blue Marble: NASA Earth Observatory, via NASA GIBS',
  },
  'gibs-viirs-day': {
    body: 'earth', title: 'The Earth on a day: VIIRS true colour (NASA)', kind: 'wms', version: '1.3.0', daily: true, role: 'base',
    url: 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi', layer: 'VIIRS_SNPP_CorrectedReflectance_TrueColor', format: 'image/jpeg', maxLevel: 8,
    credit: 'VIIRS (Suomi NPP) corrected reflectance: NASA EOSDIS LANCE and GIBS',
  },
  'gibs-blackmarble': {
    body: 'earth', title: 'Night lights: Black Marble 2016 (NASA)', kind: 'wms', version: '1.3.0', night: true, role: 'night',
    url: 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi', layer: 'VIIRS_Black_Marble', format: 'image/jpeg', maxLevel: 7,
    credit: 'Black Marble 2016: NASA Earth Observatory, via NASA GIBS',
  },
  'gibs-sst': {
    body: 'earth', title: 'The sea\'s temperature of a day (GHRSST MUR, 1 km)', kind: 'wms', version: '1.3.0', daily: true, role: 'data',
    url: 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi', layer: 'GHRSST_L4_MUR_Sea_Surface_Temperature', format: 'image/png', transparent: true, maxLevel: 7,
    credit: 'GHRSST Level 4 MUR global foundation sea surface temperature (NASA JPL, PO.DAAC), via NASA GIBS',
  },
  'gibs-sst-anomaly': {
    body: 'earth', title: 'The sea\'s temperature anomaly of a day (GHRSST MUR)', kind: 'wms', version: '1.3.0', daily: true, role: 'data',
    url: 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi', layer: 'GHRSST_L4_MUR_Sea_Surface_Temperature_Anomalies', format: 'image/png', transparent: true, maxLevel: 7,
    credit: 'GHRSST Level 4 MUR sea surface temperature anomaly (NASA JPL, PO.DAAC), via NASA GIBS',
  },
  'gibs-imerg': {
    body: 'earth', title: 'Rain of a day (GPM IMERG)', kind: 'wms', version: '1.3.0', daily: true, role: 'data',
    url: 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi', layer: 'IMERG_Precipitation_Rate', format: 'image/png', transparent: true, maxLevel: 5,
    credit: 'GPM IMERG precipitation rate (NASA GSFC), via NASA GIBS',
  },
  'gibs-coastlines': {
    body: 'earth', title: 'Coastlines', kind: 'wms', version: '1.3.0', role: 'reference',
    url: 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi', layer: 'Coastlines_15m', format: 'image/png', transparent: true, maxLevel: 12,
    credit: 'Coastlines: © OpenStreetMap contributors (ODbL), via NASA GIBS',
  },
  'gibs-reference': {
    body: 'earth', title: 'Borders and roads', kind: 'wms', version: '1.3.0', role: 'reference',
    url: 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi', layer: 'Reference_Features_15m', format: 'image/png', transparent: true, maxLevel: 12,
    credit: 'Reference features: © OpenStreetMap contributors (ODbL), via NASA GIBS',
  },
  'gibs-labels': {
    body: 'earth', title: 'Place names', kind: 'wms', version: '1.3.0', role: 'reference',
    url: 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi', layer: 'Reference_Labels_15m', format: 'image/png', transparent: true, maxLevel: 12,
    credit: 'Place labels: © OpenStreetMap contributors (ODbL), via NASA GIBS',
  },
  'trek-lro-wac': {
    body: 'moon', title: 'LRO wide-angle camera mosaic (NASA Moon Trek)', kind: 'wmts', role: 'base',
    template: 'https://trek.nasa.gov/tiles/Moon/EQ/LRO_WAC_Mosaic_Global_303ppd_v02/1.0.0/default/default028mm/{z}/{row}/{col}.jpg', maxLevel: 7,
    credit: 'LRO WAC global mosaic, 303 px/° (NASA/GSFC/Arizona State University), via NASA Moon Trek',
  },
};

// The address of a tile's image. `date` (YYYY-MM-DD) for services that change daily.
export function imageryUrl(service, z, x, y, { date, size = 256 } = {}) {
  const span = 180 / 2 ** z;
  const west = -180 + x * span, north = 90 - y * span, east = west + span, south = north - span;
  if (service.kind === 'wmts') return service.template.replace('{z}', z).replace('{row}', y).replace('{col}', x);
  const box = service.version === '1.3.0' ? `${south},${west},${north},${east}` : `${west},${south},${east},${north}`;
  const crs = service.version === '1.3.0' ? 'CRS' : 'SRS';
  const params = [
    'SERVICE=WMS', `VERSION=${service.version}`, 'REQUEST=GetMap', `LAYERS=${service.layer}`, 'STYLES=',
    `${crs}=EPSG:4326`, `BBOX=${box}`, `WIDTH=${size}`, `HEIGHT=${size}`, `FORMAT=${encodeURIComponent(service.format)}`,
  ];
  if (service.transparent) params.push('TRANSPARENT=TRUE');
  if (service.daily && date) params.push(`TIME=${date}`);
  return `${service.url}?${params.join('&')}`;
}

// A body's map as a layer of TileTree: one equirectangular image of the whole body (east longitude
// from −180° at the left edge to 180° at the right, north up), cut into the pyramid's two root tiles
// (west and east hemispheres), the deeper tiles drawing parts of them; and, where a tile service of
// the same pyramid exists (`live`: NASA Trek's WMTS for the Moon, Mars and Mercury, as the data
// library's solar-system-maps names it), its tiles from level 1 to its deepest.
//   map: { url, body, title, credit, live?: { template, levels, credit } }
export function mapLayer(id, map) {
  let whole = null;
  const image = (signal) => {
    if (!whole) {
      whole = fetch(map.url, { signal }).then((r) => (r.ok ? r.blob() : null)).then((b) => (b ? createImageBitmap(b) : null))
        .catch((e) => { whole = null; throw e; });
    }
    return whole;
  };
  const live = map.live ? { kind: 'wmts', template: map.live.template, maxLevel: (map.live.levels || 1) - 1 } : null;
  return {
    id,
    maxLevel: live ? live.maxLevel : 0,
    service: { body: map.body, title: map.title, credit: map.live ? `${map.credit}; at close range ${map.live.credit}` : map.credit, kind: 'map' },
    load: async (z, x, y, signal) => {
      if (z >= 1 && live) {
        const r = await fetch(imageryUrl(live, z, x, y), { signal, mode: 'cors' });
        if (!r.ok) return null;
        const blob = await r.blob();
        return blob.type.startsWith('image/') ? createImageBitmap(blob, { imageOrientation: 'flipY' }) : null;
      }
      if (z !== 0) return null;
      const img = await image(signal);
      if (!img) return null;
      const half = img.width / 2;
      return createImageBitmap(img, x * half, 0, half, img.height, { imageOrientation: 'flipY' });
    },
  };
}

// A loader of a service's images as ImageBitmaps (decoded off the main thread), for TileTree's
// layers; a failure resolves to null so that the tile keeps its ancestor's image.
export function imageryLayer(id, service, options = {}) {
  return {
    id,
    maxLevel: service.maxLevel,
    service,
    load: async (z, x, y, signal) => {
      const url = imageryUrl(service, z, x, y, { date: options.date?.() });
      const r = await fetch(url, { signal, mode: 'cors' });
      if (!r.ok) return null;
      const blob = await r.blob();
      if (!blob.type.startsWith('image/')) return null; // a service's XML error report
      return createImageBitmap(blob, { imageOrientation: 'flipY' });
    },
  };
}
