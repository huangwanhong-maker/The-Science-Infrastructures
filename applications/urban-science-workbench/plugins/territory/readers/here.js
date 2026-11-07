// What the world's grids know of any point on land: its H3 level 6 cell in Space2Stats, and the hazards named for its
// administrative unit by ThinkHazard!, which come only from the name of a place and are said to.

import { compact, nearest, tileAt } from '../shared.js';

const CACHE = new Map();
const HAZARD_LEVEL = { HIG: 'high', MED: 'medium', LOW: 'low', VLO: 'very low', NPR: 'no data' };

async function cellAt(lat, lon) {
  const t = await tileAt(lat, lon);
  if (!CACHE.has(t)) CACHE.set(t, nearest(t.centres, 1));
  const k = CACHE.get(t)(lon, lat, 0.08);
  return k >= 0 ? { t, k } : null;
}

export async function read(ctx, at) {
  if (Math.abs(at.lat) > 80) return null;
  const c = await cellAt(at.lat, at.lon);
  if (!c) return null;
  const { t: { fields: f }, k } = c;
  const rows = [];
  const add = (label, v, text) => { if (v !== null && v !== undefined && Number.isFinite(v)) rows.push({ label, value: text ? text(v) : compact(v) }); };
  add('people, 2025', f.sum_pop_2025[k]);
  add('of them in flood zones', f.pop_flood_pct[k], (v) => `${(v * (v <= 1 ? 100 : 1)).toFixed(0)} %`);
  add('built-up land, 2020', f.sum_built_area_m_2020[k], (v) => `${(v / 1e6).toFixed(2)} km²`);
  add('night light, 2024', f.sum_viirs_ntl_2024[k]);
  return rows.length ? { title: 'This 36 km² cell (Space2Stats)', rows } : null;
}

// For the risk tools.
export const hazardLevel = (code) => HAZARD_LEVEL[code] || code;
