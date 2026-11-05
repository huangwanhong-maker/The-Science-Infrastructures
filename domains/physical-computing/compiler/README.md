# compiler: Tasks Lowered onto Substrates

*Status: planned (ADR 0076 §2.3).* The compiler carries a task onto a target substrate. It chooses the encoding $E$, the control $u$ and the readout $R$ of the square (theory note §1.1). The generic passes live here; the lowering specific to a target lives in its family's directory under [`substrates/`](../substrates/).

## 1. The Passes

| Pass | What it does | Stands on |
|---|---|---|
| Capabilities | Matches a task to the targets that accept it, by each family's *target description*: the kinds of task, dimensions, ranges, precision, speed, energy per operation, materials | — |
| Scaling | Scales time and amplitude to the substrate's units and ranges: an analogue computer's machine units; a deep substrate's dimensionless groups (theory note §1.6) | physis `units` |
| Encoding | Dual rails for signed quantities; masks for time multiplexing; codes by rate, phase or amplitude; Ising embeddings of QUBO problems | physis `chemistry` (dual rails) |
| Readout | Linear readouts trained by ridge regression; thresholds; the signs of amplitudes; the binarisation of phases | tyche `learning` |
| Fitting and search | The substrate's free parameters found by search (evolution in materio) | heuresis |
| Lowering | The target-specific step: a netlist of integrators, a reaction network, DNA gates, a mesh's phases, a hysteron design | each family |

## 2. Invariants

- A pass never loosens a task's tolerance. If the target cannot meet it, compilation fails and says why, rather than returning a looser program.
- A compiled program states the tolerance it is expected to meet, and the runtime's readings are held to it.
