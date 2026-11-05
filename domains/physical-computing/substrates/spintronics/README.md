# spintronics: Spin-Torque Oscillators, Magnetic Networks and Stochastic Junctions

*Status: planned (ADR 0076).* Modes: reservoir, programmed. **Breadth:** networks of one kind of unit. Fast, solid state, at room temperature.

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| Spin-torque nano-oscillators, alone (Torrejon et al. 2017) and coupled (Romera et al. 2018) | amplitude relaxation as memory; synchronisation | CoFeB, MgO, Pt, Ta, W; cobalt, platinum, tantalum and tungsten are critical raw materials; made on chip lines; reusable |
| Magnetic networks and artificial spin ices (Gartside et al. 2022) | collective switching, many metastable states, complex attractors | permalloy (Ni, Fe), lithographed; reusable |
| Stochastic magnetic tunnel junctions | telegraph switching: the random source of p-bits ([`probabilistic`](../probabilistic/)) | as the oscillators |

## 2. The Simulators

1. **The auto-oscillator** (Slavin and Tiberkevich 2009). The complex amplitude $c$, with power $p = |c|^2$, obeys
   $$\dot c = -i\omega(p)c - \Gamma_+(p)c + \Gamma_-(p)c,$$
   with $\omega = \omega_0 + Np$, $\Gamma_+ = \Gamma_G(1 + Qp)$ and $\Gamma_- = \sigma I(1 - p)$. Thermal noise is complex white noise. Checks:
   - the steady power $p_0 = (\zeta - 1)/(\zeta + Q)$, where $\zeta = I/I_{th}$ and $I_{th} = \Gamma_G/\sigma$;
   - the frequency $\omega_0 + Np_0$;
   - two oscillators coupled electrically lock within Adler's range, widened by the nonlinear frequency shift by $\sqrt{1 + \nu^2}$.
2. **Magnetic networks.** Each nanomagnet is a macrospin, and its dynamics follow the Landau–Lifshitz–Gilbert equation; the magnets are coupled by their dipolar fields. Checks:
   - quasi-static switching on the Stoner–Wohlfarth astroid, $h_{sw}(\theta) = (\cos^{2/3}\theta + \sin^{2/3}\theta)^{-3/2}$;
   - precession at $\gamma\mu_0 H$;
   - energy that never increases under Gilbert damping.
3. **The stochastic junction.** A two-state telegraph process with dwell times $\tau_\pm = \tau_0 \exp(\Delta(1 \mp I/I_c))$, after Néel and Brown. Checks:
   - dwell times exponentially distributed, with the Néel–Brown mean;
   - the mean state as a function of the current.

## 3. Tasks and their Lowering

- **Filters** (simulator 1): the input modulates the current, and the readout is the amplitude sampled at virtual nodes, as Torrejon et al. did.
- **Phenomena:** synchronisation, the switching of networks.
- **Energy** (simulator 2): the ground states of spin ices.
- **A random source** (simulator 3) for the probabilistic family.

## 4. What the Foundations Must Gain

- physis `dynamics`: stochastic flows;
- the Landau–Lifshitz–Gilbert equation, a general law, in physis (its `electromagnetism`, planned, or a module of magnetism; physis's decision).

## 5. References

- Slavin and Tiberkevich (2009); Torrejon et al. (2017); Romera et al. (2018): in the theory note.
- Adler, R. (1946), A study of locking phenomena in oscillators, *Proceedings of the IRE* 34, 351–357.
- Brown, W. F. (1963), Thermal fluctuations of a single-domain particle, *Physical Review* 130, 1677–1686.
- Gartside, J. C. et al. (2022), Reconfigurable training and reservoir computing in an artificial spin-vortex ice via spin-wave fingerprinting, *Nature Nanotechnology* 17, 460–469.
- Gilbert, T. L. (2004), A phenomenological theory of damping in ferromagnetic materials, *IEEE Transactions on Magnetics* 40, 3443–3449.
- Stoner, E. C. and Wohlfarth, E. P. (1948), A mechanism of magnetic hysteresis in heterogeneous alloys, *Philosophical Transactions of the Royal Society A* 240, 599–642.
