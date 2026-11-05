# ir: The Intermediate Representation

*Status: optional, future (ADR 0076 §2.3, §2.7).* At the user's request, it is not implemented now. Until it exists, tasks, compiled programs and schedules are C++ values. The pipeline is built so that it can be added later without changing the stages around it.

## 1. What it Must Carry

- **Tasks** of the six kinds ([`tasks/`](../tasks/)), with their metric, tolerance and check.
- **Programs** compiled for a target: the encoding, control and readout, with their parameters and units.
- **Schedules** ([`schedule/`](../schedule/)), kept apart from programs.
- **Continuous time** and **tolerances**, where calx's IR has discrete operations and declared precisions.
- **A version**, so that a text written today can still be read when the format changes.

## 2. When it is Built

It is built once the compilers of two or more families show what it must carry, so that it is designed from instances and not ahead of them (CLAUDE.md §2.2: an abstraction after its pattern has appeared twice).

## 3. Open Question

Whether it is calx's IR, extended with continuous time and a tolerance in place of a precision; or a dialect beside it, in the manner of MLIR; or the domain's own. The first is calx's decision.
