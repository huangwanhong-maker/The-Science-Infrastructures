# chemical: The BZ Medium, Reaction Networks as Programs, and Electrochemical Oscillators

*Status: planned (ADR 0076).* Modes: analogue, programmed, reservoir. **Breadth:** the broadest with analogue electronics. Continuous reaction networks realise any polynomial ODE after a dual-rail encoding (Fages et al. 2017), and their exact theory says what they compute (theory note §1.2). Priority: first, on physis's `chemistry`, already written.

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| The Belousov–Zhabotinsky medium | oscillations, waves, spirals, Turing patterns, chaos | malonic acid, bromate (a strong oxidiser), a cerium or ferroin catalyst, sulfuric acid; consumed |
| Programmable BZ arrays (Parrilla-Gutierrez et al. 2020) | 5 × 5 cells coupled and controlled; memory and recognition of patterns | as above |
| Reaction networks as programs | stable computation of the semilinear functions; any polynomial ODE by dual rails | ideal networks, realised by [`molecular`](../molecular/) |
| Electrochemical oscillators (Kiss, Zhai and Hudson 2002) | Kuramoto's transition seen on 64 electrodes | nickel in acid; the electrodes dissolve slowly |

## 2. The Simulators

1. **The Oregonator.** Tyson and Fife's two variables, $\varepsilon\dot u = u(1 - u) - fv(u - q)/(u + q)$ and $\dot v = u - v$, and Field and Noyes's three. They run on grids and graphs by physis's `reaction_diffusion`, with a photosensitive term for control by light. Checks:
   - the Hopf points in $f$ found by continuation agree with the zero of the trace at the steady state, solved independently;
   - Turing's condition on a graph agrees with physis;
   - the threshold of excitability under light.
2. **BZ arrays.** Coupled Oregonators on the graph of cells. Check: two coupled cells synchronise in phase or in anti-phase as the phase reduction predicts from the coupling's sign.
3. **Reaction networks as programs.** On physis's networks:
   - semilinear functions (min, max, sums, multiples), checked exactly by the network's Markov chain;
   - dual-rail compilations of polynomial ODEs (Lorenz 63), checked against the ODE.
4. **Electrochemical oscillators.** A phase model with the electrodes' spread of natural frequencies, coupled through a shared resistance; a two-variable model of nickel's electrodissolution comes later, provisional. Check: $K_c = 2/(\pi g(0))$, and $r(K) = \sqrt{1 - K_c/K}$ for Lorentzian frequencies.

## 3. Tasks and their Lowering

- **Flows and functions** (simulator 3): the dual rails, and the network.
- **Phenomena** (simulators 1 and 2): a Hopf point, waves.
- **Filters:** a BZ reservoir, its input as light, its readout from a camera's image.

## 4. What the Foundations Must Gain

- physis `chemistry`: the dual-rail compiler;
- physis `dynamics`: continuation and its Hopf points.

## 5. References

- Field and Noyes (1974); Parrilla-Gutierrez et al. (2020); Chen, Doty and Soloveichik (2014); Fages et al. (2017); Kiss, Zhai and Hudson (2002): in the theory note.
- Field, R. J., Körös, E. and Noyes, R. M. (1972), Oscillations in chemical systems. II. Thorough analysis of temporal oscillation in the bromate–cerium–malonic acid system, *Journal of the American Chemical Society* 94, 8649–8664.
- Tyson, J. J. and Fife, P. C. (1980), Target patterns in a realistic model of the Belousov–Zhabotinskii reaction, *Journal of Chemical Physics* 73, 2224–2237.
