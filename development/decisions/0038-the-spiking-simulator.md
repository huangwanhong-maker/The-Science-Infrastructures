# 0038. The Spiking Simulator

- **Status:** Accepted
- **Date:** 2026-10-23
- **Decided by:** the user asked for the simulator of spiking networks to be built as an application, after the engine and its models, and then a visualizer and profiler, with simulations profiled, visualized and added to the gallery; and set the scope: simulators now, classical and quantum, with devices and the compiler that would reach them only planned. Claude designed the application and wrote this record.

## 1. Context

The neuroscience domain's `snn` module (ADR 0036) is an engine: twenty-two neuron models, plasticity, a parallel and deterministic network, and a command-line simulator that reads a specification in JSON. It has no place where a network is composed, run and read. ADR 0036 placed such an application apart from the Brain Atlas.

## 2. Decision

1. **`applications/spiking-simulator/`**, a web application on port 8830, in the pattern of the atlases: ES modules without a build step, a server in Python's standard library, a recipe, its own checks.
2. **The engine runs on the server, as a program.** The page sends a specification; the server runs the domain's command-line simulator on it and returns its results. The C++ is not compiled to the browser: the engine's parallelism and its speed are the point, and the domain's library stays pure (it reads only its standard input).
3. **The engine is built by the recipe** with g++, natively or in WSL, into `data/cache/engines/snn/`, with the hash of its sources recorded so that a changed header makes it stale; the server refuses to run a stale engine and the page offers to build it.
4. **The specification is the document.** The page edits the engine's own JSON, so that what is shown, saved and run is one object, and a network made in the page runs from the command line.
5. **Limits on a public address**: fewer neurons, shorter runs, fewer threads, and building only with a token.
6. **In stages, as the user ordered:** the simulator first; then the visualizer and profiler, whose analysis of spike trains and profile of the engine belong to the domain's C++, checked there, with the page drawing their results; then simulations profiled and drawn for the gallery.

## 3. Alternatives Considered

- **WebAssembly in the page:** no toolchain for it here, a single thread in most browsers, and a second build of the engine to keep checked.
- **Python bindings first:** they serve notebooks, not a page; they remain planned.
- **The Brain Atlas as host:** ADR 0036 and the user placed the networks apart from it.

## 4. Consequences

- The applications are four; CLAUDE.md §5 and §7 and the README name the new one and its port.
- A machine that runs the page needs g++ (or WSL on Windows), as the domain's checks already do.
