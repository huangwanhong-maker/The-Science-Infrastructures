// What the plug-in's layer and tool share: the district (the shell's box, network and buildings), the record of
// its events and its event complex as the engine built them, how the spacetime is drawn, and listeners, so that the
// map and the drawer show one record. A module is loaded once by its URL, so the layer and the tool hold this object.
//
// Time is a Modified Julian Date in UTC (the engine's): the clock's instant is MJD = ms / 86 400 000 + 40 587.
// Height is time: z = scale · (years since the window's start) + lift, scale metres a year, from the ground or from
// above the district's tallest building.

import { districtOf, boxOf, networkOf, footprintsIn } from '/src/district.js';

export const state = {
  area: null, box: null, record: null, request: null, error: null, cut: null,
  view: { scale: 10, mode: 'refined', lift: 'ground', show: { building: true, street: true, junction: false, 'party-wall': true, frontage: false, 'cross-section': false }, censuses: true },
};
const listeners = new Set();
export function on(f) { listeners.add(f); return () => listeners.delete(f); }
export function emit(what) { for (const f of [...listeners]) try { f(what); } catch (e) { console.error(e); } }

export const mjdOf = (clock) => clock.instant / 86400000 + 40587;
export const dateOf = (mjd) => new Date((mjd - 40587) * 86400000).toISOString().slice(0, 10);
export const yearOf = (mjd) => 1970 + (mjd - 40587) / 365.2425;

export function setView(changes) { Object.assign(state.view, changes); emit('view'); }

// The height of the tallest building of the district: where the spacetime starts when lifted above it.
export function skyline() {
  const R = state.record;
  if (!R) return 0;
  return R.couplings.reduce((m, c) => (c.type === 'building' ? Math.max(m, c.top || 0) : m), 0);
}
export function z(mjd) {
  const R = state.record;
  const base = state.view.lift === 'above' ? skyline() + 20 : 0;
  return base + (state.view.scale * (mjd - R.window.from)) / 365.2425;
}

let building = null;
export function ensureRecord(ctx) {
  if (state.area !== ctx.area.id) Object.assign(state, { area: ctx.area.id, box: null, record: null, error: null, cut: null }), building = null;
  if (state.record) return Promise.resolve(state.record);
  if (building) return building;
  building = (async () => {
    if (!ctx.model) throw new Error('the city model is not built yet');
    state.box = boxOf(ctx.frame, districtOf(ctx.area));
    const net = networkOf(ctx.model, state.box);
    const buildings = footprintsIn(ctx.model, state.box);
    state.request = { area: ctx.area.id, ...net, buildings };
    state.record = await ctx.ask('record', state.request);
    emit('record');
    return state.record;
  })();
  building.catch((e) => { building = null; state.error = String(e.message || e); emit('error'); });
  return building;
}

// The network of the moment at an instant, from the engine (the couplings and interfaces whose lives cross it).
export async function cutAt(ctx, mjd) {
  if (!state.record) return null;
  const R = state.record;
  const t = Math.min(R.window.to - 1e-3, Math.max(R.window.from, mjd));
  // The engine keeps the record it last built; should it have been restarted since, the district goes again.
  try { state.cut = await ctx.ask('cut', { t }); } catch { state.cut = await ctx.ask('cut', { t, ...state.request }); }
  emit('cut');
  return state.cut;
}

// Colours (sRGB, 0–1): the cells by type in the workbench's sea and sand; the events by type in the visualization
// layer's warm categorical hues, ordered so that neighbours stay apart in every vision (palettes.js).
const hex = (s) => [1, 3, 5].map((i) => parseInt(s.slice(i, i + 2), 16) / 255);
export const EVENT_TYPES = [
  { type: 'birth', name: 'built', hue: '#c9d38a' },
  { type: 'alteration', name: 'refurbished or completed', hue: '#8fb3dd' },
  { type: 'institution', name: 'registered (heritage)', hue: '#b9a2e0' },
  { type: 'regime', name: 'a street event (permit)', hue: '#e48f78' },
  { type: 'death', name: 'demolished', hue: '#e0a3b6' },
  { type: 'occurrence', name: 'an occurrence', hue: '#f0cf9a' },
];
export const eventColour = (type) => hex((EVENT_TYPES.find((e) => e.type === type) || EVENT_TYPES[5]).hue);
export function cellColour(type, night) {
  const day = { building: '#5b8fb9', street: '#c98b6b', laneway: '#c98b6b', mall: '#c98b6b', lane: '#c98b6b', junction: '#8a99a6',
    'party-wall': '#4f8a8b', frontage: '#d69a8e', 'cross-section': '#a9b4bc' };
  const dark = { building: '#8fb9da', street: '#e0ae90', laneway: '#e0ae90', mall: '#e0ae90', lane: '#e0ae90', junction: '#7f8d99',
    'party-wall': '#8fbdb6', frontage: '#e2b0a5', 'cross-section': '#5d6a74' };
  return hex((night ? dark : day)[type] || (night ? '#7f8d99' : '#8a99a6'));
}
export const kindOfCell = (type) => (['street', 'laneway', 'mall', 'lane'].includes(type) ? 'street' : type);
export const css = (rgb) => `rgb(${rgb.map((c) => Math.round(c * 255)).join(',')})`;
