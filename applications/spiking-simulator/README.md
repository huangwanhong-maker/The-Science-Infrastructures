# Spiking Simulator

*Status: draft (2026-10-23). Serves use and learning. Draws on: neuroscience (its `snn` engine), through calx and physis.* Decision: [ADR 0038](../../development/decisions/0038-the-spiking-simulator.md).

## 1. Purpose

A place to compose networks of spiking neurons and run them: a population of any of the neuroscience domain's twenty-two models, classical or quantum, connected by projections with their rules, delays and plasticity, driven by devices, and recorded. The simulation is done by the domain's C++ engine, exact where the dynamics are linear and giving the same spikes on any number of threads; the page composes, runs and reports, and every network is a specification (JSON, with its quantities in units) that the command-line simulator runs as well.

It is apart from the Brain Atlas, which draws the brain's anatomy and signals; this application is for the dynamics of networks.

## 2. Concepts

- **A specification is the network.** The page holds one specification, the very JSON the engine reads (`domains/neuroscience/snn/include/neuroscience/snn/spec.hpp`): settings, populations, projections, devices, recorders. The editor changes it field by field; the Specification tab shows it as text, to edit, save or open. A network saved here runs unchanged with `snn run`.
- **Quantities keep their units.** A field takes "20 ms", "0.1 mV" or "8 kHz", a bare number where the quantity has no dimension, or a distribution, `{"uniform": ["0 mV", "20 mV"]}`. The engine checks every dimension against the model's and says what it refuses.
- **The engine runs on the server.** The page sends the specification to its server, which runs the engine as a separate program and returns its results: spike times and senders, binned rates, recorded states with their units, weights, and the engine's own figures (run time, spikes, synaptic events).
- **A run is shown and profiled.** Its spikes are drawn as a raster over the population rates; its activity is described population by population (rates, the irregularity of intervals as CV and LV, the Fano factor of counts, Golomb's synchrony χ, correlations, and the rhythm from the spectrum of the rate), with the distributions behind them and, where each neuron had its own current, the f–I curve; recorded states and weights are traced; the network is drawn. The statistics are the engine's own, computed and checked in the domain's C++ (`snn/analysis.hpp`), and so is the profile: where the engine's time went, phase by phase and thread by thread, with advice when a network is too small for its threads.
- **The engine is built here.** It is compiled once with g++ (11 or later), natively or in WSL on Windows, into `data/cache/engines/snn/`, and compiled again when its sources change; the application's recipe does it, from the command line or from the page.

## 3. Use

```
python data/acquisition/recipe.py applications/spiking-simulator/recipe.json   # build the engine (a minute)
python applications/spiking-simulator/serve.py --open
```

It opens at http://127.0.0.1:8830/. Choose an example, or compose a network: add populations from the catalogue on the right, connect them, drive them and record them on the left, and run (Ctrl+Enter). The run's views are tabs: Activity (drag across a chart to look closer in time, double-click to see all of it), Statistics, States, Weights, Profile, Network, Summary (and its files), Specification. A link `#preset=brunel-ai&run=1` opens an example and runs it; `&tab=statistics` opens a tab of its own, and `&from=300&to=350` the activity on a span of time, in milliseconds.

The examples:

| Example | What it shows | Source |
|---|---|---|
| One neuron and its adaptation | An AdEx neuron's spikes spreading as its adaptation current builds | Brette and Gerstner 2005 |
| Brunel's network, four states | Asynchronous irregular, synchronous regular, fast and slow oscillations of 12,500 neurons | Brunel 2000, Fig. 8 |
| Classes of excitability, in one run | The f-I curves of four conductance models, each neuron with its own current | Hodgkin 1948; Rinzel and Ermentrout 1998 |
| Quantum emitters | Photons of driven and pulsed two-level emitters along quantum trajectories: Rabi oscillation, decay, antibunching | Dalibard, Castin and Mølmer 1992 |
| A neuron learns from its inputs | STDP balancing a thousand Poisson inputs onto one neuron | Song, Miller and Abbott 2000 |

## 4. Structure

| Path | Content |
|---|---|
| `engine.py` | Builds the engine (finding g++ natively or in WSL), keeps its catalogue of models, runs a specification and can cancel a run |
| `serve.py` | The page and its interface: `/api/engine`, `/api/models`, `/api/run`, `/api/cancel`, and the recipe at `/setup` |
| `recipe.json` | The one step that builds the engine |
| `src/spec.js` | The specification: its parts, references kept whole on renaming and removal, the problems found before a run |
| `src/units.js` | Quantities shown with their prefixes and written as the engine reads them |
| `src/editor.js`, `src/catalogue.js`, `src/results.js` | The editor, the catalogue of models, the summary of a run |
| `src/views.js` | The visualizer and profiler: activity, statistics, states, weights, the engine's profile, the network |
| `src/charts.js` | Charts on canvas: scales, ticks, hairline axes, lines, bands, bars, dots, a crosshair with its tooltip, a range dragged to look closer |
| `src/presets.js` | The examples |
| `src/palette.js` | The populations' colours, validated for colour-blind readers on the page's surface |
| `tests/run-tests.mjs` | Its checks |

## 5. Limits

On this computer a run may have 250,000 neurons, 60 s of simulated time and 15 minutes of computing, on every core. On a public address the server allows 50,000 neurons, 10 s and one minute, on at most four threads, and builds the engine only for the holder of its token (`--token`); each limit can be set (`--max-neurons`, `--max-duration`, `--max-seconds`). One run at a time; Stop ends it.

## 6. Checks

`node applications/spiking-simulator/tests/run-tests.mjs`: quantities shown and read; the charts' ticks and scales; the specification's references kept whole when populations and projections are renamed or removed; every example whole, and every model, parameter, component and receptor it names present in the engine's catalogue; and, when the engine is built, every example run by the engine with its neurons firing, its activity described neuron by neuron and its profile given thread by thread. The science is checked with the neuroscience domain (its `snn` module).

## 7. Open Questions

1. **A run that reports as it goes,** for long simulations, and results too large to send whole (rasters thinned on the server).
2. **Comparisons:** two runs side by side, or a parameter swept across runs.
3. **Python bindings** of the engine, beside the command line and this page.
