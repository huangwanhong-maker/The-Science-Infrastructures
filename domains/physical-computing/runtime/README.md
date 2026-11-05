# runtime: Programs Run on Backends

*Status: planned (ADR 0076 §2.3).* The runtime runs a compiled program under its schedule on a backend, and returns the readings with their provenance.

## 1. The Backend Interface

Every backend, of every family, answers four requests:

- `reset(state)`;
- `apply(control, input, duration)`;
- `read(sensors)`, which returns readings in their units;
- `describe()`, which returns the target description the compiler reads.

There are three kinds of backend:

1. **The simulator**, written first for every family. It is the substrate's model, deterministic under its seed. It declares its integrator, its step and its own error, which is kept apart from the task's tolerance.
2. **The device**, later. A real substrate is an instrument behind `observation/`, and a work drives it; the domain's runtime drives simulators only, so that the core stays pure (CLAUDE.md §3.3). A device is held to its family's simulator.
3. **The surrogate**, optional. A learned digital twin of a device, held to the device's readings.

## 2. Provenance

Every reading records what made it:

- the kind of backend;
- the model and its version;
- the seed;
- the conditions (temperature, drive).

It records too whether it was observed or modelled, so that a reading from a device is never taken for one from a simulator. This follows the distinction Cosmopoiesis draws among its instruments.
