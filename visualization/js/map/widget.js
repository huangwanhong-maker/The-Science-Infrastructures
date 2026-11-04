// An interactive 2D map in an SVG element: the base (the sphere, the graticule, the countries of a
// GeoJSON), an overlay group in which an application draws its layers, and the instruments of a GIS:
// zoom by the wheel and buttons, pan by dragging, a pick when the pointer is pressed and released in
// place, and a readout of the projection, the place under the pointer and the zoom. What lies in the
// overlay is drawn by the application, which `onSettle` calls again when the view has changed, so that
// sizes given in screen pixels stay so (px() is the drawing's units of one pixel at the present zoom).
//
// Areas in dispute (a GeoJSON such as Natural Earth's breakaway and disputed areas) are drawn over the
// countries in no claimant's colour, hatched on the paper, and without a name (ADR 0069): whatever
// colours the countries does not colour them. Their colours are the page's --map-disputed-ground and
// --map-disputed-hatch; the hatch keeps its spacing on the screen at every zoom.

import { planar } from './projection.js';
import { identity, transform, zoomAt, fit, centreOn } from './view.js';

const SVG = 'http://www.w3.org/2000/svg';
const make = (tag, attrs = {}) => {
  const n = document.createElementNS(SVG, tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== null && v !== undefined) n.setAttribute(k, v);
  return n;
};

// The code a Natural Earth country is known by: ISO_A3, else ISO_A3_EH, else ADM0_A3.
export function naturalEarthCode(p) {
  for (const k of ['ISO_A3', 'ISO_A3_EH', 'ADM0_A3']) if (p[k] && p[k] !== '-99' && p[k].length === 3) return p[k];
  return null;
}

export function createMap(host, { projection = 'equal-earth', scale = 300, countries = null, disputed = null, codeOf = naturalEarthCode, nameOf = (p) => p.NAME || '',
  maxZoom = 18, onPointer = null, onLeave = null, onPick = null, onSettle = null, readout = null } = {}) {
  const P = planar(projection, scale);
  const state = { view: { ...identity }, built: false };
  let world, over, hatch = null, settle = null;

  function build() {
    host.setAttribute('viewBox', P.viewBox());
    while (host.firstChild) host.removeChild(host.firstChild);
    world = make('g', { class: 'world' });
    const sphere = P.outline();
    const defs = make('defs');
    const clip = make('clipPath', { id: `${host.id || 'map'}-sphere` });
    clip.append(make('path', { d: sphere }));
    defs.append(clip);
    world.append(defs, make('path', { class: 'sphere', d: sphere }), make('path', { class: 'graticule', d: P.graticule(), 'clip-path': `url(#${host.id || 'map'}-sphere)` }));
    const g = make('g', { class: 'countries' });
    for (const f of countries?.features || []) {
      g.append(make('path', { class: 'country', d: P.feature(f.geometry), 'data-code': codeOf(f.properties) || '', 'data-name': nameOf(f.properties) || '' }));
    }
    const areas = make('g', { class: 'disputed' });
    if (disputed?.features?.length) {
      const id = `${host.id || 'map'}-disputed`;
      hatch = make('pattern', { id, patternUnits: 'userSpaceOnUse', width: 8, height: 8 });
      hatch.append(make('rect', { width: 8, height: 8, style: 'fill: var(--map-disputed-ground, #EEEDE8); stroke: none' }),
        make('path', { d: 'M0 0V8', style: 'stroke: var(--map-disputed-hatch, #A9A69C); stroke-width: 1.6; fill: none' }));
      defs.append(hatch);
      // (Names are not kept: an area in dispute is not named.)
      for (const f of disputed.features) areas.append(make('path', { class: 'area-disputed', d: P.feature(f.geometry), fill: `url(#${id})`, 'data-disputed': '1' }));
    }
    over = make('g', { class: 'over' });
    world.append(g, areas, over);
    host.append(world);
    apply();
    state.built = true;
  }

  const apply = () => {
    world?.setAttribute('transform', transform(state.view));
    hatch?.setAttribute('patternTransform', `rotate(45) scale(${(1 / state.view.k).toFixed(4)})`);
  };
  const pxPerUnit = () => Math.max(1, host.clientWidth) / host.viewBox.baseVal.width;
  const point = (e) => {
    const p = host.createSVGPoint();
    p.x = e.clientX;
    p.y = e.clientY;
    return p.matrixTransform(host.getScreenCTM().inverse());
  };
  // The place under the pointer, [lon, lat] in degrees, or null off the map.
  const place = (e) => { const p = point(e); return P.unproject((p.x - state.view.x) / state.view.k, (p.y - state.view.y) / state.view.k); };

  function show(e) {
    if (!readout) return;
    let where = '';
    const ll = e ? place(e) : null;
    if (ll) where = ` · ${Math.abs(ll[1]).toFixed(1)}° ${ll[1] >= 0 ? 'N' : 'S'}, ${Math.abs(ll[0]).toFixed(1)}° ${ll[0] >= 0 ? 'E' : 'W'}`;
    readout.replaceChildren(Object.assign(document.createElement('b'), { textContent: P.name }), `${where} · ×${state.view.k.toFixed(1)}`);
  }

  function setView(v, now = false) {
    state.view = v;
    apply();
    show(null);
    clearTimeout(settle);
    if (now) onSettle?.(state.view);
    else settle = setTimeout(() => onSettle?.(state.view), 140);
  }

  function wire() {
    host.addEventListener('wheel', (e) => {
      e.preventDefault();
      const p = point(e);
      setView(zoomAt(state.view, state.view.k * Math.exp(-e.deltaY * 0.0016), p.x, p.y, { max: maxZoom }));
    }, { passive: false });
    let drag = null, moved = false;
    host.addEventListener('pointerdown', (e) => {
      if (e.button !== 0) return;
      drag = { p: point(e), x: state.view.x, y: state.view.y };
      moved = false;
    });
    host.addEventListener('pointermove', (e) => {
      if (drag) {
        const p = point(e);
        const dx = p.x - drag.p.x, dy = p.y - drag.p.y;
        if (!moved && Math.hypot(dx, dy) > 4) { moved = true; host.setPointerCapture(e.pointerId); host.classList.add('dragging'); onLeave?.(); }
        if (moved) { state.view = { k: state.view.k, x: drag.x + dx, y: drag.y + dy }; apply(); return; }
      }
      onPointer?.(e);
      show(e);
    });
    host.addEventListener('pointerup', (e) => {
      if (drag && !moved) onPick?.(e);
      if (drag && moved) setView(state.view);
      drag = null;
      host.classList.remove('dragging');
    });
    host.addEventListener('pointerleave', () => { onLeave?.(); show(null); });
    host.addEventListener('dblclick', () => setView({ ...identity }, true));
  }

  build();
  wire();
  show(null);
  return {
    planar: P, host,
    get view() { return state.view; },
    get world() { return world; },
    get over() { return over; },
    countries: () => host.querySelectorAll('.countries .country'),
    // The drawing's units of one screen pixel at the present zoom.
    px: () => 1 / pxPerUnit() / state.view.k,
    setView,
    zoomBy: (f) => setView(zoomAt(state.view, state.view.k * f, 0, 0, { max: maxZoom })),
    reset: () => setView({ ...identity }, true),
    fit: (points, opts) => setView(fit(points, P.W, P.H, { max: maxZoom, ...opts }), true),
    centreOn: (x, y) => setView(centreOn(state.view, x, y), true),
    point, place, rebuild: build,
  };
}
