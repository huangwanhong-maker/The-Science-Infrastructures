// A raster on the ground about the view, asked for any place as the camera moves: the year each cell was first built (WSF
// Evolution, with the scale of years above the clock, so that the city grows on the ground), the degree of urbanisation (GHS-SMOD)
// or the depth of a river flood (GloFAS, at a return period one chooses). It is laid on the Earth's own surface by its shader.

import { card, gradient, tone } from '../shared.js';

const URBAN = [null, 'very low density', 'low density', 'rural cluster', 'suburban', 'semi-dense', 'dense cluster', 'urban centre'];

function kinds(theme) {
  const orchid = tone('orchid', theme, { from: 0.08 });
  const sea = tone('sea', theme, { from: 0.15, keep: theme === 'day' });
  return {
    settlement: {
      title: 'Settlement spreading', nearest: true, floor: 0.002,
      ramp: (u) => { if (u < 0.02) return [0, 0, 0, 0]; const c = orchid(Math.min(1, u)); return [c[0], c[1], c[2], 0.92]; },
      opacity: 0.9, source: 'World Settlement Footprint Evolution, DLR, CC BY 4.0', note: 'The year each 30 m cell first stood built, seen from Landsat: the older the deeper.',
    },
    urbanisation: {
      title: 'Degree of urbanisation', nearest: true, floor: 0.002,
      ramp: (u) => { const k = Math.floor(u * 8); if (k < 1 || k > 7) return [0, 0, 0, 0]; const c = orchid(k / 7); return [c[0], c[1], c[2], 0.25 + 0.09 * k]; },
      opacity: 0.85, source: 'GHS-SMOD R2023A for 2025, European Commission JRC, CC BY 4.0',
    },
    flood: {
      title: 'River floods', nearest: false, floor: 0.002,
      ramp: (u) => { if (u <= 0) return [0, 0, 0, 0]; const c = sea(Math.min(1, Math.sqrt(u / 0.6))); return [c[0], c[1], c[2], Math.min(0.92, 0.4 + 4 * u)]; },
      opacity: 0.92, source: 'GloFAS river flood hazard maps, CEMS, European Union, CC BY 4.0', note: 'The depth a river flood of this return period reaches, on a grid of about 90 m. Small streams and the sea are not in it.',
    },
  };
}

export async function layer(ctx) {
  const { earth, dom: { h } } = ctx;
  const kind = ctx.spec.raster;
  const K = kinds(ctx.theme)[kind];
  let dead = false, shown = null, asking = false, rp = 100, cut = 1;
  const status = h('div', { class: 'ld-quiet' }, 'asking for the ground here…');

  async function refresh(force = false) {
    if (dead || asking) return;
    const v = earth.view, DEG = 180 / Math.PI, lat = v.lat * DEG, lon = v.lon * DEG;
    const reach = kind === 'urbanisation' ? 1.2 : 0.3;
    if (!force && shown && Math.abs(lat - shown.lat) < reach && Math.abs(lon - shown.lon) < reach / Math.max(0.3, Math.cos(v.lat))) return;
    asking = true;
    status.textContent = 'asking for the ground here…';
    try {
      const q = `kind=${kind}&lat=${lat.toFixed(3)}&lon=${lon.toFixed(3)}${kind === 'flood' ? `&rp=${rp}` : ''}`;
      const r = await fetch(`/world/raster?${q}`);
      const meta = await r.json();
      if (!r.ok) throw new Error(meta.message || `the service answered ${r.status}`);
      const bytes = new Uint8Array(await (await fetch(meta.url)).arrayBuffer());
      if (dead) return;
      earth.setDrape({ data: bytes, width: meta.width, height: meta.height, box: meta.box, ramp: K.ramp, nearest: K.nearest, opacity: K.opacity, floor: K.floor, cut });
      shown = { lat, lon, key: q };
      status.textContent = `${meta.width.toLocaleString('en')} × ${meta.height.toLocaleString('en')} cells about this view`;
    } catch (e) {
      status.textContent = String(e.message || e);
      shown = { lat, lon };   // not asked again until the view has moved on
    } finally {
      asking = false;
    }
  }
  const timer = setInterval(refresh, 1500);
  refresh();

  let years = null;
  if (kind === 'settlement') {
    const span = Array.from({ length: 31 }, (_, i) => 1985 + i);
    years = ctx.years({ label: 'Standing by', years: span, value: 2015, onChange: (y) => { cut = ((y - 1984) * 8 + 2) / 255; earth.setDrapeCut(cut); } });
  }
  const extra = [];
  if (kind === 'flood') {
    const keys = h('div', { class: 'keys-tight', title: 'a flood expected once in this many years' }, [10, 20, 50, 100, 200, 500].map((p) => h('button', { class: p === rp ? 'on' : '', onclick: (e) => {
      rp = p;
      for (const b of e.currentTarget.parentElement.children) b.classList.toggle('on', b === e.currentTarget);
      refresh(true);
    } }, String(p))));
    extra.push(h('div', { class: 'ld-line' }, 'once in so many years'), keys);
  }
  if (kind === 'urbanisation') {
    extra.push(h('div', { class: 'legend ld-swatches' }, ...URBAN.map((name, k) => (name ? h('span', {}, h('i', { class: 'dot', style: `background:${gradientStop(K, k)}` }), name) : null)).filter(Boolean)));
  }
  const hasRamp = kind !== 'urbanisation';
  const colour = (u) => { const c = K.ramp(kind === 'flood' ? 0.02 + u * 0.6 : 0.03 + 0.97 * u); return [c[0], c[1], c[2]]; };
  const legend = card(h, K.title, { extra, colour: hasRamp ? colour : null, low: kind === 'settlement' ? '1985' : hasRamp ? '0' : undefined, high: kind === 'settlement' ? '2015' : '7.6 m', note: K.note, source: K.source });
  legend.append(status);
  return {
    legend,
    dispose() {
      dead = true;
      clearInterval(timer);
      earth.setDrape(null);
      if (years) years.close();
    },
  };
}

function gradientStop(K, k) {
  const c = K.ramp((k * 32 + 16) / 256);
  return `rgb(${c.slice(0, 3).map((v) => Math.round(v * 255)).join(',')})`;
}
export { gradient };
