// The network being composed: a specification of the neuroscience domain's engine (snn), held
// as the very JSON the engine reads (domains/neuroscience/snn/include/neuroscience/snn/spec.hpp),
// so that what the page shows, saves and runs are one object.
//
// A specification has its settings (duration, resolution, seed, integrator), its populations
// (a model, a size, parameters, initial states, a bias), its projections (a source, a target,
// a receptor, a rule of connection, a weight, a delay, a synapse), its devices and its
// recorders. Quantities are written with their units ("20 ms"); the engine checks their
// dimensions against the models'.

export const RULES = {
  all_to_all: { label: 'all to all', value: null },
  one_to_one: { label: 'one to one', value: null },
  fixed_indegree: { label: 'fixed in-degree', value: 'inputs per neuron' },
  fixed_outdegree: { label: 'fixed out-degree', value: 'outputs per neuron' },
  fixed_probability: { label: 'fixed probability', value: 'probability' },
  uniform: { label: 'uniform (mean field)', value: null },
};

// The synapses and the fields each takes, with the unit of each ('w' is the weight's).
export const SYNAPSES = {
  static: [],
  stdp: [['w_max', 'w'], ['lambda', '1'], ['alpha', '1'], ['mu_plus', '1'], ['mu_minus', '1'], ['tau_plus', 's'], ['tau_minus', 's']],
  triplet: [['w_max', 'w'], ['A2_plus', '1'], ['A3_plus', '1'], ['A2_minus', '1'], ['A3_minus', '1'], ['tau_plus', 's'], ['tau_minus', 's'], ['tau_x', 's'], ['tau_y', 's']],
  tsodyks: [['U', '1'], ['tau_rec', 's'], ['tau_fac', 's']],
};

// The devices and their fields: 'in' is the target's input unit, 'w' a receptor's weight.
export const DEVICES = {
  poisson: { label: 'Poisson spikes', fields: [['rate', 'Hz'], ['receptor', 'receptor'], ['weight', 'w']] },
  dc: { label: 'constant current', fields: [['amplitude', 'in']] },
  step_current: { label: 'step current', fields: [['times', 'list s'], ['amplitudes', 'list in']] },
  noise: { label: 'noise current', fields: [['mean', 'in'], ['sd', 'in'], ['interval', 's']] },
  sine: { label: 'sine current', fields: [['offset', 'in'], ['amplitude', 'in'], ['frequency', 'Hz'], ['phase', '1']] },
  spike_times: { label: 'spikes at times', fields: [['times', 'list s'], ['receptor', 'receptor'], ['weight', 'w']] },
};
export const DEVICE_COMMON = [['start', 's'], ['stop', 's'], ['first', 'index'], ['count', 'index']];

export const RECORDERS = {
  spikes: { label: 'spikes', fields: [['population', 'population'], ['first', 'index'], ['count', 'index']] },
  rate: { label: 'population rate', fields: [['population', 'population'], ['bin', 's'], ['first', 'index'], ['count', 'index']] },
  state: { label: 'state', fields: [['population', 'population'], ['neurons', 'list index'], ['components', 'components'], ['interval', 's']] },
  weights: { label: 'weights', fields: [['projection', 'projection'], ['interval', 's']] },
};

export const INTEGRATORS = ['default', 'dopri5', 'rk4', 'heun', 'euler', 'rush_larsen'];

export function empty() {
  return { duration: '1 s', resolution: '0.1 ms', seed: 1, populations: [], projections: [], devices: [], recorders: [] };
}

// A specification with every list present (a file or an example may leave some out).
export function normalise(spec) {
  for (const k of ['populations', 'projections', 'devices', 'recorders']) if (!Array.isArray(spec[k])) spec[k] = [];
  return spec;
}

// A copy, so that presets and the page never share an object.
export const clone = (x) => JSON.parse(JSON.stringify(x));

// The rule's kind, and its value, of a projection.
export function ruleOf(p) {
  const r = p.rule ?? {};
  const kind = Object.keys(RULES).find((k) => k in r) ?? (('pairs' in r) ? 'pairs' : 'all_to_all');
  return { kind, value: r[kind], autapses: r.autapses, multapses: r.multapses };
}

export function setRule(p, kind, value) {
  const old = ruleOf(p);
  const rule = { [kind]: RULES[kind]?.value ? (value ?? old.value ?? 1) : true };
  if (old.autapses !== undefined) rule.autapses = old.autapses;
  if (old.multapses !== undefined) rule.multapses = old.multapses;
  p.rule = rule;
}

export function unique(names, base) {
  if (!names.includes(base)) return base;
  for (let i = 2; ; i++) if (!names.includes(`${base} ${i}`)) return `${base} ${i}`;
}

// A new population of a model, named after it.
export function addPopulation(spec, model, size = 100) {
  const name = unique(spec.populations.map((p) => p.name), model.replace(/_/g, ' '));
  spec.populations.push({ name, model, size });
  return name;
}

// Renames a population everywhere it is named.
export function renamePopulation(spec, from, to) {
  for (const p of spec.populations) if (p.name === from) p.name = to;
  for (const p of spec.projections) {
    if (p.source === from) p.source = to;
    if (p.target === from) p.target = to;
  }
  for (const d of spec.devices) if (d.target === from) d.target = to;
  for (const r of spec.recorders) if (r.population === from) r.population = to;
}

// Removes a population and everything that names it; weight recorders follow their projections.
export function removePopulation(spec, name) {
  const kept = spec.projections.map((p) => p.source !== name && p.target !== name);
  const index = [];
  let n = 0;
  for (const k of kept) index.push(k ? n++ : -1);
  spec.populations = spec.populations.filter((p) => p.name !== name);
  spec.projections = spec.projections.filter((_, i) => kept[i]);
  spec.devices = spec.devices.filter((d) => d.target !== name);
  spec.recorders = spec.recorders
    .filter((r) => r.population !== name && !(r.model === 'weights' && index[r.projection] < 0))
    .map((r) => (r.model === 'weights' ? { ...r, projection: index[r.projection] } : r));
}

export function removeProjection(spec, i) {
  spec.projections.splice(i, 1);
  spec.recorders = spec.recorders
    .filter((r) => !(r.model === 'weights' && r.projection === i))
    .map((r) => (r.model === 'weights' && r.projection > i ? { ...r, projection: r.projection - 1 } : r));
}

// What is wrong with a specification before the engine sees it: references, sizes, names.
// The engine checks the rest (models, parameters, dimensions) and says what it refuses.
export function problems(spec, catalogue = null) {
  normalise(spec);
  const out = [];
  const names = spec.populations.map((p) => p.name);
  const model = (n) => catalogue?.get(spec.populations.find((p) => p.name === n)?.model);
  if (!spec.populations.length) out.push('the network has no population');
  names.forEach((n, i) => {
    if (!n) out.push(`population ${i + 1} has no name`);
    else if (names.indexOf(n) !== i) out.push(`two populations are named "${n}"`);
  });
  for (const p of spec.populations) {
    if (!(Number.isInteger(p.size) && p.size > 0)) out.push(`${p.name}: its size is a whole number above zero`);
    if (catalogue && !catalogue.get(p.model)) out.push(`${p.name}: no model "${p.model}"`);
  }
  spec.projections.forEach((p, i) => {
    const what = `projection ${i + 1} (${p.source} → ${p.target})`;
    if (!names.includes(p.source)) out.push(`${what}: no population "${p.source}"`);
    if (!names.includes(p.target)) out.push(`${what}: no population "${p.target}"`);
    const m = model(p.target);
    if (m && !m.receptors.some((r) => r.name === p.receptor)) out.push(`${what}: ${m.name} has no receptor "${p.receptor}"`);
  });
  for (const d of spec.devices) if (!names.includes(d.target)) out.push(`a ${d.model} device names no population ("${d.target}")`);
  for (const r of spec.recorders) {
    if (r.model === 'weights') {
      if (!(r.projection >= 0 && r.projection < spec.projections.length)) out.push('a weight recorder names no projection');
    } else if (!names.includes(r.population)) out.push(`a ${r.model} recorder names no population ("${r.population}")`);
  }
  return out;
}

// The specification as the engine is given it: empty fields dropped.
export function forEngine(spec) {
  const strip = (o) => {
    if (Array.isArray(o)) return o.map(strip);
    if (o && typeof o === 'object') {
      const r = {};
      for (const [k, v] of Object.entries(o)) {
        if (v === null || v === undefined || v === '') continue;
        if (typeof v === 'object' && !Array.isArray(v) && Object.keys(v).length === 0 && k !== 'rule') continue;
        r[k] = strip(v);
      }
      return r;
    }
    return o;
  };
  return strip(spec);
}

// The neurons in a specification.
export const neurons = (spec) => spec.populations.reduce((n, p) => n + (Number(p.size) || 0), 0);
