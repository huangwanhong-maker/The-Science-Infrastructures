// What the views share: the datasets known, the one each view draws, and the layers plug-ins have made.
// A change is announced to whoever listens, so that the map, the dock and the panel follow it.

export const shared = {
  datasets: [],                       // the library's products and one's own, from the server
  flows: 'library:trade',             // the flows the map draws
  currencies: 'library:currencies',   // the record spacetime draws
  layers: [],                         // [{key, file, title, language, year, dataset, result}]: plug-ins' results
  colour: 'attraction',               // what colours the places: 'attraction' or a layer's key
  arcs: 'all',                        // which arcs: 'all', 'driven', 'circulating', 'none' or a layer's key
};

const listeners = new Map();

export function on(event, f) {
  if (!listeners.has(event)) listeners.set(event, new Set());
  listeners.get(event).add(f);
}

export function emit(event, detail) {
  for (const f of listeners.get(event) || []) f(detail);
}

export const datasetsOf = (kind) => shared.datasets.filter((d) => d.kind === kind);
export const dataset = (id) => shared.datasets.find((d) => d.id === id) || null;

export function setDatasets(list) {
  shared.datasets = list || [];
  for (const kind of ['flows', 'currencies']) {
    const key = kind === 'flows' ? 'flows' : 'currencies';
    if (!dataset(shared[key])) shared[key] = datasetsOf(kind)[0]?.id || null;
  }
  emit('datasets');
}
