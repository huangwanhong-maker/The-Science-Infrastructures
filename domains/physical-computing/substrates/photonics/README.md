# photonics: Interferometer Meshes, Delay Loops and Coherent Ising Machines

*Status: planned (ADR 0076).* Modes: programmed, reservoir. **Breadth:** any unitary, by a mesh; delay dynamics, by a loop. Fast; linear operations cheap, strong nonlinearity costly.

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| Meshes of Mach–Zehnder interferometers (Shen et al. 2017) | any linear unitary map | silicon or silicon nitride photonics (silicon metal is a critical raw material); room temperature |
| Delay-loop reservoirs (Appeltant et al. 2011; Larger et al. 2012) | one nonlinear node with delayed feedback; time-multiplexed virtual nodes; delay-induced chaos | fibre, modulators (lithium niobate: niobium is a critical raw material) |
| Coherent Ising machines (Inagaki et al. 2016; McMahon et al. 2016) | degenerate parametric oscillators that settle into Ising states | fibre loops, pump lasers |

## 2. The Simulators

1. **The mesh.** A target unitary is factorised into the phases $(\theta, \phi)$ of interferometers by Clements's scheme (physis `symmetry`). Each interferometer has a loss, phase errors and a finite extinction. Checks:
   - the reconstruction equals the target to rounding;
   - a uniform loss $l$ scales the norm by $(1 - l)^{\text{depth}}$.
2. **The delay loop.** $\tau_R \dot x = -x + \beta \sin^2(x(t - \tau_D) + \phi + \gamma\, m(t)u(t))$ (after Ikeda 1979), with a mask $m(t)$ setting $N$ virtual nodes $\theta$ apart. Checks:
   - the fixed point loses stability when $|\beta \sin 2(x^* + \phi)| > 1$, the map limit of a long delay;
   - the memory capacity is at most $N$;
   - the echo state property holds below the instability.
3. **The coherent Ising machine.** $\dot x_i = (p - 1 - x_i^2)x_i + \xi\sum_j J_{ij}x_j + \text{noise}$. Checks:
   - the threshold of oscillation $p_{th} = 1 - \xi\lambda_{\max}(J)$, in closed form;
   - on small instances, the ground states are compared with exact enumeration.

## 3. Tasks and their Lowering

- **Functions:** a linear map is factorised by its singular values, $U\Sigma V^\dagger$, into two meshes and attenuators.
- **Filters** (simulator 2): the mask and the virtual nodes come from the schedule.
- **Energy** (simulator 3).

## 4. What the Foundations Must Gain

- physis `symmetry`: the factorisation of a unitary into two-mode rotations;
- physis `dynamics`: delay equations;
- physis `statistical_mechanics`: Ising on graphs.

## 5. References

- Reck et al. (1994); Clements et al. (2016); Shen et al. (2017); Appeltant et al. (2011); Larger et al. (2012); Inagaki et al. (2016): in the theory note.
- Ikeda, K. (1979), Multiple-valued stationary state and its instability of the transmitted light by a ring cavity system, *Optics Communications* 30, 257–261.
- McMahon, P. L. et al. (2016), A fully programmable 100-spin coherent Ising machine with all-to-all connections, *Science* 354, 614–617.
