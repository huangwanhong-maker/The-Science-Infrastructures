# analogue-electronics: Analogue Computers, Programmable Arrays and Oscillators

*Status: planned (ADR 0076).* Modes: analogue, reservoir, programmed. **Breadth:** the broadest of all families, since Shannon's general-purpose analogue computer generates exactly the solutions of polynomial ODEs (theory note §1.6). Priority: first.

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| Analogue computers: integrators, summers, multipliers, coefficients | any polynomial dynamics within ±1 machine unit, chaos among it | silicon and copper; reusable, programmed by patching; low power |
| Field-programmable analogue arrays (Hasler 2020) | the same netlists, routed by switches | as above; reprogrammed electrically |
| CMOS oscillators; oscillator-based Ising machines (Wang and Roychowdhury 2019) | synchronisation; phases binarised by second-harmonic injection to solve Ising problems | standard CMOS |

## 2. The Simulators

1. **The netlist.** A polynomial ODE is compiled to integrators, summers, multipliers and coefficients, with time scaled by $\tau$ and amplitudes to machine units. Non-idealities: each integrator's finite gain, as a leak $\delta$; offsets; multiplier error; saturation at ±1; thermal noise. It stands on physis `dynamics`, calx's integrators and random streams. Checks:
   - the ideal netlist of Lorenz 63, scaled, equals the direct integration to rounding over a short horizon;
   - a leak $\delta$ in each of the $n$ integrators lowers the divergence, and so the sum of the Lyapunov exponents, by exactly $n\delta$;
   - the harmonic oscillator $\ddot x = -x$ decays as $e^{-\delta t}$;
   - saturation is reported, never silent.
2. **The programmable array.** The netlist on a routing fabric, each switch a series resistance and a parasitic capacitance, so that each connection is a first-order lag. Check: as the lag tends to 0, the array's run tends to the netlist's.
3. **Oscillators.** A ring oscillator of $N$ stages has frequency $1/(2Nt_d)$. A network of them follows the phase model $\dot\theta_i = \omega_i - K\sum_j J_{ij}\sin(\theta_i - \theta_j) - K_s \sin 2\theta_i$. Checks:
   - above a threshold in $K_s$ the phases settle in $\{0, \pi\}$;
   - on small instances, the Ising states found are compared with exact enumeration (physis);
   - Kuramoto's threshold $K_c = 2/(\pi g(0))$ holds when $K_s = 0$.

## 3. Tasks and their Lowering

- **Flows** (simulators 1 and 2): the netlist; the scaling is to machine units.
- **Filters** (simulator 1 as a reservoir): an input injected at a summer; a linear readout of the integrators.
- **Energy** (simulator 3): $J$ and $h$ to couplings and biases; the readout is the phase binarised.
- **Phenomena:** a bifurcation crossed by a coefficient potentiometer.

## 4. What the Foundations Must Gain

- physis `dynamics`: flows of runtime dimension, Lyapunov spectra;
- physis `statistical_mechanics`: Ising on graphs, exactly for small $n$.

## 5. References

- Hasler, J. (2020), Large-scale field-programmable analog arrays, *Proceedings of the IEEE* 108, 1283–1302.
- Shannon (1941); Graça and Costa (2003); Ulmann (2022); Mead (1989); Huang et al. (2017): in the theory note.
- Wang, T. and Roychowdhury, J. (2019), OIM: oscillator-based Ising machines for solving combinatorial optimisation problems, *Unconventional Computation and Natural Computation* (UCNC 2019), LNCS 11493, 232–256.
