# probabilistic: p-bits, Probabilistic Computers and Thermodynamic Computing

*Status: planned (ADR 0076).* Mode: programmed (sampling). **Breadth:** any Boltzmann law on spins with pairwise couplings, and Gaussian laws; strong for criticality and for sampling. Priority: first.

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| p-bits and probabilistic computers (Camsari et al. 2017; Borders et al. 2019) | sampling of Boltzmann laws; invertible logic; phase transitions | stochastic tunnel junctions ([`spintronics`](../spintronics/)) or CMOS; room temperature |
| Thermodynamic linear algebra (Aifer et al. 2024) | linear systems solved by the stationary law of coupled noisy units | analogue electronics |
| Thermodynamic sampling units, as Extropic reports them | Gibbs sampling of graphical models in hardware | CMOS. The claims of efficiency are projections, to be read in their primary form |

## 2. The Simulators

1. **A network of p-bits.** $m_i = \mathrm{sgn}(\tanh(\beta I_i) - r_i)$ with $I_i = \sum_j J_{ij}m_j + h_i$ and $r_i$ uniform on $(-1, 1)$. Updates are sequential or by colours of the graph. The device's delays and correlated noise are declared. Checks:
   - the stationary law is Boltzmann's, against exact enumeration for small $n$ (physis);
   - the AND gate of Camsari et al. has the truth table as its ground states;
   - with its output clamped, it samples consistent inputs.
2. **Thermodynamic linear algebra.** $dx = -Ax\,dt + \sqrt{2/\beta}\,dW$ has stationary covariance $A^{-1}/\beta$, and its time averages solve linear systems. Checks:
   - the Ornstein–Uhlenbeck stationary covariance, in closed form;
   - an error that falls as the inverse square root of the time;
   - the dissipation against the thermodynamic uncertainty relation's floor (physis).
3. **Sampling of graphical models.** Heat-bath updates on factor graphs, the units of a thermodynamic sampler. Provisional; checked against exact marginals.

## 3. Tasks and their Lowering

- **Energy:** a QUBO problem is mapped to $J$ and $h$; the schedule anneals $\beta(t)$; the readout is the empirical law or the best state.
- **Functions:** invertible logic; linear algebra (simulator 2).

## 4. What the Foundations Must Gain

- physis `statistical_mechanics`: Ising on graphs (exact and by Glauber's dynamics); the thermodynamic uncertainty relation;
- calx: streams, already written.

## 5. References

- Camsari et al. (2017): in the theory note.
- Aifer, M. et al. (2024), Thermodynamic linear algebra, *npj Unconventional Computing* 1, 13. To be checked.
- Borders, W. A. et al. (2019), Integer factorization using stochastic magnetic tunnel junctions, *Nature* 573, 390–393.
- Glauber, R. J. (1963), Time-dependent statistics of the Ising model, *Journal of Mathematical Physics* 4, 294–307.
