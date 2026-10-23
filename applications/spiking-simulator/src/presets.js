// Networks to start from: each a specification of the engine with what it shows and the
// papers it comes from. They are the simulator's worked examples, and its gallery.

// n values from a to b, evenly, each written with its unit.
const ramp = (n, a, b, unit, digits = 4) => Array.from({ length: n }, (_, i) => `${Number((a + ((b - a) * i) / (n - 1)).toPrecision(digits))} ${unit}`);

// Brunel's sparse network of excitatory and inhibitory integrate-and-fire neurons (Brunel 2000,
// J. Comput. Neurosci. 8, 183-208): N_E = 4 N_I, each neuron receiving C_E = 0.1 N_E excitatory
// and C_I = C_E / 4 inhibitory inputs of J = 0.1 mV and -g J, after D = 1.5 ms, and external
// Poisson input at eta nu_thr per synapse, nu_thr = theta / (J C_E tau) the rate that alone
// would bring the mean potential to threshold. Its four states (Fig. 8) differ in g and eta.
function brunel({ g, eta, duration = '1 s', scale = 1, bin = '0.5 ms', analysis = { from: '100 ms' } }) {
  const NE = Math.round(10000 * scale), NI = Math.round(2500 * scale), CE = Math.round(NE / 10), CI = Math.round(CE / 4);
  const J = 0.1, theta = 20, tau = 20e-3;
  const nuThr = theta / (J * CE * tau);   // Hz
  const neuron = { tau_m: '20 ms', t_ref: '2 ms', E_L: '0 mV', V_reset: '10 mV', V_th: '20 mV' };
  const pop = (name, size) => ({ name, model: 'lif_psc_delta', size, parameters: neuron, initial: { V: { uniform: ['0 mV', '20 mV'] } } });
  const proj = (source, target, C, w) => ({ source, target, receptor: 'delta', rule: { fixed_indegree: C, autapses: false }, weight: `${Number(w.toPrecision(6))} mV`, delay: '1.5 ms' });
  return {
    duration, resolution: '0.1 ms', seed: 1, analysis,
    populations: [pop('E', NE), pop('I', NI)],
    projections: [proj('E', 'E', CE, J), proj('I', 'E', CI, -g * J), proj('E', 'I', CE, J), proj('I', 'I', CI, -g * J)],
    devices: ['E', 'I'].map((t) => ({ model: 'poisson', target: t, rate: `${Number((eta * nuThr * CE).toPrecision(6))} Hz`, receptor: 'delta', weight: `${J} mV` })),
    recorders: [
      { model: 'spikes', population: 'E', first: 0, count: Math.min(NE, 400) },
      { model: 'spikes', population: 'I', first: 0, count: Math.min(NI, 100) },
      { model: 'rate', population: 'E', bin },
      { model: 'rate', population: 'I', bin },
      { model: 'state', population: 'E', neurons: [0, 1, 2], components: ['V'], interval: '0.1 ms' },
    ],
  };
}

export const PRESETS = [
  {
    id: 'adex',
    title: 'One neuron and its adaptation',
    about: 'An adaptive exponential integrate-and-fire neuron under a step of current: its spikes come closer at first, then spread as the adaptation current w builds up with each spike (Brette and Gerstner 2005).',
    spec: () => ({
      duration: '600 ms', resolution: '0.1 ms', seed: 1,
      populations: [{ name: 'neuron', model: 'adex', size: 1 }],
      devices: [{ model: 'dc', target: 'neuron', amplitude: '800 pA', start: '100 ms', stop: '500 ms' }],
      recorders: [
        { model: 'spikes', population: 'neuron' },
        { model: 'state', population: 'neuron', neurons: [0], components: ['V', 'w'], interval: '0.1 ms' },
      ],
    }),
  },
  {
    id: 'brunel-ai',
    title: "Brunel's network: asynchronous and irregular",
    about: 'Twelve thousand five hundred integrate-and-fire neurons, inhibition five times as strong as excitation (g = 5) and external drive twice threshold (eta = 2): each neuron fires irregularly at a low rate, and the population almost steadily (Brunel 2000, Fig. 8C).',
    spec: () => brunel({ g: 5, eta: 2 }),
  },
  {
    id: 'brunel-sr',
    title: "Brunel's network: synchronous and regular",
    about: 'Excitation dominant (g = 3): the neurons fire near their highest rate, regularly and almost in step (Brunel 2000, Fig. 8A). The busiest of the four states.',
    spec: () => brunel({ g: 3, eta: 2, duration: '300 ms' }),
  },
  {
    id: 'brunel-si-fast',
    title: "Brunel's network: fast oscillation",
    about: 'Strong inhibition and strong drive (g = 6, eta = 4): the population oscillates fast, near 180 Hz, while each neuron fires irregularly and seldom, skipping most cycles (Brunel 2000, Fig. 8B).',
    spec: () => brunel({ g: 6, eta: 4, duration: '500 ms', bin: '0.25 ms', analysis: { from: '100 ms', bin: '0.25 ms', sync_bin: '1 ms' } }),
  },
  {
    id: 'brunel-si-slow',
    title: "Brunel's network: slow oscillation",
    about: 'Inhibition dominant and weak drive (g = 4.5, eta = 0.9): a slow oscillation of the population, the neurons firing in irregular bursts (Brunel 2000, Fig. 8D).',
    spec: () => brunel({ g: 4.5, eta: 0.9, duration: '2 s', bin: '1 ms', analysis: { from: '200 ms', bin: '2 ms', sync_bin: '5 ms' } }),
  },
  {
    id: 'excitability',
    title: 'Classes of excitability, in one run',
    about: 'Four populations, each neuron given its own constant current, rising along the population: the firing rate against the current is each model\'s f-I curve. Wang and Buzsáki\'s interneuron and Morris and Lecar\'s SNIC set begin to fire at any low rate (class 1); Hodgkin and Huxley\'s axon and Morris and Lecar\'s Hopf set begin at a finite rate (class 2) (Hodgkin 1948; Rinzel and Ermentrout 1998).',
    spec: () => ({
      duration: '3 s', resolution: '0.05 ms', seed: 1, analysis: { from: '500 ms' },
      populations: [
        { name: 'Wang-Buzsáki', model: 'wang_buzsaki', size: 120, bias: ramp(120, 0, 100, 'pA') },
        { name: 'Morris-Lecar SNIC', model: 'morris_lecar', size: 120, parameters: { g_Ca: '400 nS', V_3: '12 mV', V_4: '17.4 mV', phi: '66.6667 Hz' }, bias: ramp(120, 3, 6, 'nA') },
        { name: 'Hodgkin-Huxley', model: 'hodgkin_huxley', size: 120, bias: ramp(120, 0, 1500, 'pA') },
        { name: 'Morris-Lecar Hopf', model: 'morris_lecar', size: 120, bias: ramp(120, 6, 14, 'nA') },
      ],
      recorders: [
        { model: 'spikes', population: 'Wang-Buzsáki' },
        { model: 'spikes', population: 'Morris-Lecar SNIC' },
        { model: 'spikes', population: 'Hodgkin-Huxley' },
        { model: 'spikes', population: 'Morris-Lecar Hopf' },
      ],
    }),
  },
  {
    id: 'quantum',
    title: 'Quantum emitters',
    about: 'Two-level emitters followed along exact quantum trajectories; each detected photon is a spike. The driven ones, starting together in their ground state, emit with a damped Rabi oscillation that settles to the master equation\'s rate; the undriven ones are excited by a pi pulse every 100 ms and then decay at gamma. No emitter ever gives two photons at once (antibunching).',
    spec: () => ({
      duration: '600 ms', resolution: '0.1 ms', seed: 1, analysis: { bin: '1 ms', window: '50 ms', sync_bin: '1 ms' },
      populations: [
        { name: 'driven', model: 'quantum_emitter', size: 2000, parameters: { Omega: '400 Hz', gamma: '100 Hz' } },
        { name: 'pulsed', model: 'quantum_emitter', size: 2000, parameters: { Omega: '0 Hz', gamma: '100 Hz' } },
      ],
      devices: [{ model: 'spike_times', target: 'pulsed', receptor: 'pulse', weight: 3.141592653589793, times: ['50 ms', '150 ms', '250 ms', '350 ms', '450 ms', '550 ms'] }],
      recorders: [
        { model: 'spikes', population: 'driven', first: 0, count: 200 },
        { model: 'spikes', population: 'pulsed', first: 0, count: 200 },
        { model: 'rate', population: 'driven', bin: '1 ms' },
        { model: 'rate', population: 'pulsed', bin: '1 ms' },
      ],
    }),
  },
  {
    id: 'stdp',
    title: 'A neuron learns from its inputs',
    about: 'A thousand Poisson inputs at 10 Hz reach one neuron through synapses with additive STDP, depression slightly stronger than potentiation. The weights drift until the neuron fires at a moderate rate, whatever their start: STDP balances the input (Song, Miller and Abbott 2000).',
    spec: () => ({
      duration: '20 s', resolution: '0.1 ms', seed: 1,
      populations: [
        { name: 'inputs', model: 'parrot', size: 1000 },
        { name: 'neuron', model: 'lif_psc_exp', size: 1 },
      ],
      projections: [{
        source: 'inputs', target: 'neuron', receptor: 'excitatory', rule: { all_to_all: true },
        weight: { uniform: ['0 pA', '40 pA'] }, delay: '1 ms',
        synapse: { model: 'stdp', w_max: '40 pA', lambda: 0.005, alpha: 1.05, tau_plus: '20 ms', tau_minus: '20 ms' },
      }],
      devices: [{ model: 'poisson', target: 'inputs', rate: '10 Hz', receptor: 'spike', weight: 1 }],
      recorders: [
        { model: 'spikes', population: 'neuron' },
        { model: 'spikes', population: 'inputs', first: 0, count: 100 },
        { model: 'rate', population: 'neuron', bin: '500 ms' },
        { model: 'weights', projection: 0, interval: '100 ms' },
      ],
    }),
  },
];

export const preset = (id) => PRESETS.find((p) => p.id === id);
