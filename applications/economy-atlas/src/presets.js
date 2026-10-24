// Examples to start from: the transfer networks of the user's paper on loop deficits (Huang 2026,
// The Geometry of Value) and input–output tables from the textbook (Miller and Blair 2009).

export const NETWORKS = [
  {
    id: 'triangle',
    title: 'The paper\'s triangle',
    about: 'Three positions passing value round: 1.2, 0.5 and 1.0 registered per unit surrendered. Round the loop 0.60 of each unit comes back: a deficit of −ln 0.6 = 0.511, whoever describes the value how.',
    network: {
      positions: ['a', 'b', 'c'],
      legs: [{ from: 0, to: 1, ratio: 1.2, flow: 10, kind: 'value' }, { from: 1, to: 2, ratio: 0.5, flow: 12, kind: 'value' }, { from: 2, to: 0, ratio: 1.0, flow: 6, kind: 'value' }],
    },
  },
  {
    id: 'five',
    title: 'Five positions, six legs',
    about: 'The paper\'s example of the tree gauge: a spanning tree r→1, r→2, 1→3, 2→4, and two chords, 1→2 and 3→4, which carry the holonomies of their cycles, 0.975 and 1.0725. Two numbers describe the network up to re-description.',
    network: {
      positions: ['r', '1', '2', '3', '4'],
      legs: [{ from: 0, to: 1, ratio: 1.3 }, { from: 0, to: 2, ratio: 0.8 }, { from: 1, to: 3, ratio: 0.9 }, { from: 2, to: 4, ratio: 1.5 }, { from: 1, to: 2, ratio: 0.6 }, { from: 3, to: 4, ratio: 1.1 }],
    },
  },
  {
    id: 'fee',
    title: 'A fee round a circuit',
    about: 'A worker sells to a firm, the firm to a platform, and the platform pays the worker after keeping a tenth: every leg but one registers value at par, and the loop has a deficit of −ln 0.9 = 0.105. No price system can describe it.',
    network: {
      positions: ['worker', 'firm', 'platform'],
      legs: [{ from: 0, to: 1, ratio: 1.0, flow: 100, kind: 'labour' }, { from: 1, to: 2, ratio: 1.0, flow: 100, kind: 'goods' }, { from: 2, to: 0, ratio: 0.9, flow: 100, kind: 'money' }],
    },
  },
  {
    id: 'prices',
    title: 'A market with prices',
    about: 'Four goods traded at the ratios of their prices (1, 2, 5 and 8): every loop returns exactly what went in, the deficits all vanish, and the potential recovers the log prices.',
    network: {
      positions: ['bread', 'cloth', 'tools', 'land'],
      legs: [{ from: 0, to: 1, ratio: 0.5 }, { from: 1, to: 2, ratio: 0.4 }, { from: 2, to: 3, ratio: 0.625 }, { from: 0, to: 2, ratio: 0.2 }, { from: 1, to: 3, ratio: 0.25 }, { from: 3, to: 0, ratio: 8 }],
    },
  },
  {
    id: 'round-trip',
    title: 'A round trip',
    about: 'Two positions sending each other the same amount: the flow is pure circulation round the hole it makes, with no position gaining or losing; the transport 0.95 out and 1.0 back is a deficit of 0.051 on each round.',
    network: {
      positions: ['A', 'B'],
      legs: [{ from: 0, to: 1, ratio: 0.95, flow: 30, kind: 'goods' }, { from: 1, to: 0, ratio: 1.0, flow: 30, kind: 'money' }],
    },
  },
];

export const TABLES = [
  {
    id: 'miller-blair',
    title: 'Miller and Blair\'s two sectors',
    about: 'The textbook\'s example (2009, §2.3): agriculture and manufacturing, with A = [[0.15, 0.25], [0.20, 0.05]] and L = (1/0.7575) [[0.95, 0.25], [0.20, 0.85]].',
    table: { sectors: ['agriculture', 'manufacturing'], z: [[150, 500], [200, 100]], final_demand: [350, 1700] },
  },
  {
    id: 'three',
    title: 'Three sectors with a cycle',
    about: 'Agriculture feeds industry, industry equips services, services serve agriculture, and agriculture also feeds services directly: a circular flow, of which the decomposition separates what runs round from what is driven.',
    table: { sectors: ['agriculture', 'industry', 'services'], z: [[40, 120, 30], [25, 180, 140], [60, 90, 110]], final_demand: [210, 455, 540] },
  },
];
