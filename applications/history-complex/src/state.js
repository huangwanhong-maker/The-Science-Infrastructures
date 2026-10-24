// The page's state in the address: which view, which canon, which date, what is followed, which
// disciplines are hidden, which loop of H₁ is drawn; so that a view can be linked to and returned to.
//
//   #view=loom&canon=europe-1785-1828&date=1814-11-15&focus=being:beethoven&hide=art,science&loop=2
//   #canon=europe-1785-1828&mode=spin&open=1&time=events     (the foam, the default view)

export const VIEWS = ['foam', 'loom', 'moment', 'complex', 'record', 'readings', 'field'];
export const FOAM_MODES = ['discipline', 'spin', 'flips', 'inference'];

export function readHash(hash = typeof location !== 'undefined' ? location.hash : '') {
  const p = new URLSearchParams(String(hash).replace(/^#/, ''));
  const focus = (p.get('focus') || '').match(/^(being|tie|occasion):(.+)$/);
  const loop = p.get('loop');
  return {
    view: VIEWS.includes(p.get('view')) ? p.get('view') : 'foam',
    canon: p.get('canon') || null,
    date: /^-?\d{1,4}(-\d\d){0,2}$/.test(p.get('date') || '') ? p.get('date') : null,
    focus: focus ? { kind: focus[1], id: focus[2] } : null,
    hide: new Set((p.get('hide') || '').split(',').filter(Boolean)),
    loop: loop !== null && /^\d+$/.test(loop) ? Number(loop) : null,
    tab: p.get('tab') || null,
    mode: FOAM_MODES.includes(p.get('mode')) ? p.get('mode') : 'discipline',
    opened: p.get('open') === '1',
    time: p.get('time') === 'events' ? 'events' : 'calendar',
    tour: /^\d+$/.test(p.get('tour') || '') ? Number(p.get('tour')) : null,
  };
}

export function writeHash(s) {
  const p = new URLSearchParams();
  if (s.view && s.view !== 'foam') p.set('view', s.view);
  if (s.canon) p.set('canon', s.canon);
  if (s.date) p.set('date', s.date);
  if (s.focus) p.set('focus', `${s.focus.kind}:${s.focus.id}`);
  if (s.hide && s.hide.size) p.set('hide', [...s.hide].join(','));
  if (s.loop !== null && s.loop !== undefined) p.set('loop', String(s.loop));
  if (s.tab) p.set('tab', s.tab);
  if (s.mode && s.mode !== 'discipline') p.set('mode', s.mode);
  if (s.opened) p.set('open', '1');
  if (s.time === 'events') p.set('time', 'events');
  if (s.tour !== null && s.tour !== undefined) p.set('tour', String(s.tour));
  return `#${p.toString().replace(/%2C/g, ',').replace(/%3A/g, ':')}`;
}
