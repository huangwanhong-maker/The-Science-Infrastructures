# tasks: What is Computed

*Status: draft (2026-11-06): the kinds and tolerances, `include/physical_computing/tasks/task.hpp`.* The first stage of the pipeline: a task says what is to be computed, on no particular substrate, with the set $K$ it must hold on, its metric $d$, its tolerance $\varepsilon$ and its check (theory note §1.1). Until the [IR](../ir/) exists, a task is a C++ value.

## 1. The Kinds of Task

| Kind | What | Checked by | Examples |
|---|---|---|---|
| Filter | a causal map from input series to output series, with fading memory | the normalised error on held-out inputs | delayed recall; NARMA10 (Atiya and Parlos 2000); prediction of Mackey–Glass; spoken digits (the Free Spoken Digit Dataset, CC BY-SA 4.0) |
| Flow | a dynamical system to realise: a polynomial ODE, a PDE by the method of lines, a delay system | invariants, not trajectories (theory note §1.4): the Lyapunov spectrum, the invariant measure by the energy distance, the bifurcations; or a trajectory over a horizon in Lyapunov times | Lorenz 63 and 96; Kuramoto–Sivashinsky (Pathak et al. 2018); polynomial food webs |
| Energy | an Ising, QUBO or Potts problem, to minimise or to sample at an inverse temperature β | exact enumeration where small; known optima; the divergence from the Boltzmann law | Max-Cut; invertible logic (Camsari et al. 2017) |
| Function | a function on counts or bits, or a linear map | exactly, or to a declared precision | semilinear functions of reaction networks; Boolean circuits; matrix–vector products |
| Phenomenon | a behaviour to produce | its normal form: the place and kind of a bifurcation, its scaling, the early warnings before it | a Hopf point crossed; a pattern of a given wavelength; synchrony at Kuramoto's threshold |
| Learning | a filter the substrate learns by its own rule (mode 4) | the rule's convergence, and the filter learnt | coupled learning of a linear map |

## 2. Each Task Carries

- its input set and its law (for filters, the input process);
- its metric, tolerance and budget (of time, energy and repetitions);
- its check, a function of the readings and nothing else;
- its units, through physis's `units`, so that the compiler can scale it to a substrate.

## 3. References

- Atiya, A. F. and Parlos, A. G. (2000), New results on recurrent network training: unifying the algorithms and accelerating convergence, *IEEE Transactions on Neural Networks* 11, 697–709.
- The rest are in the [theory note](../theory/physical-computation.md).
