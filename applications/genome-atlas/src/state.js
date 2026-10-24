// The page's state in the address, so that a view can be linked to and returned to:
//
//   #view=genome&genome=ncbi:NC_000913.3&begin=3920000&end=3932000&feature=812
//   #view=align&a=HEAGAWGHEE&b=PAWHEAE&alphabet=protein&matrix=BLOSUM50

export const VIEWS = ['genome', 'align', 'index', 'reads', 'assembly', 'expression', 'trees', 'populations', 'run'];

export function readHash(hash = typeof location !== 'undefined' ? location.hash : '') {
  const p = new URLSearchParams(String(hash).replace(/^#/, ''));
  const out = { view: VIEWS.includes(p.get('view')) ? p.get('view') : 'genome', genome: p.get('genome') || null, params: {} };
  for (const [k, v] of p) if (k !== 'view' && k !== 'genome') out.params[k] = v;
  return out;
}

export function writeHash(s) {
  const p = new URLSearchParams();
  if (s.view && s.view !== 'genome') p.set('view', s.view);
  if (s.genome) p.set('genome', s.genome);
  for (const [k, v] of Object.entries(s.params || {})) if (v !== null && v !== undefined && v !== '') p.set(k, String(v));
  return `#${p.toString().replace(/%3A/g, ':').replace(/%2C/g, ',')}`;
}

// A number from the hash, within bounds, or the fallback.
export function hashNumber(params, key, fallback, lo = -Infinity, hi = Infinity) {
  const raw = params[key];
  const v = Number(raw);
  return raw !== undefined && raw !== '' && Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : fallback;
}
