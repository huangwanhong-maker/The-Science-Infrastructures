# learning-networks: Networks that Learn by Local Rules

*Status: planned (ADR 0076).* Mode: learning, by the substrate's own rules (theory note §1.2, mode 4). Each element adjusts itself from what it alone sees, so no digital computer trains the network, and the cost of training, often the largest, falls away. It fits physis closely: a resistor or flow network is a Kirchhoff network with conductances. Priority: first.

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| Resistor networks (Dillavou et al. 2022) | supervised learning by coupled learning, in hardware | resistors and transistors; reusable |
| Flow networks | conductances learnt through the radii of tubes | tubing, water |
| Elastic networks (Stern et al. 2021; Pashine et al. 2019) | stiffnesses learnt, allosteric responses, directed ageing | polymer or rubber networks |

## 2. The Simulators

1. **Coupled learning in a resistor network.** In the free state the inputs are fixed. In the clamped state the outputs are nudged a fraction $\eta$ towards their targets. Each conductance changes by $\Delta k_e = -(\alpha/\eta)\,[(\Delta V_e^C)^2 - (\Delta V_e^F)^2]$, clipped to its range, with the steps and limits of transistor resistors (Dillavou et al.). Checks:
   - Kirchhoff's laws: conductances in series and in parallel;
   - as $\eta \to 0$, the rule descends the true cost to first order;
   - a realisable linear map is learnt to a small error.
2. **Flow networks.** The same rule on conductances $\pi r^4/(8\mu L)$, learnt through the radii. Check: it agrees with simulator 1 under that change of parameters.
3. **Elastic networks.** Linear springs in two dimensions; the equilibria by linear solves; the stiffnesses learnt. Check: a realisable target ratio of strains is learnt.
4. **Equilibrium propagation** (Scellier and Bengio 2017), an alternative rule. Check: its gradient agrees with differentiation through the network's equilibrium.

## 3. Tasks and their Lowering

- **Learning** (all four simulators): the inputs and outputs are placed on nodes; the schedule alternates free and clamped phases.
- **Functions:** the static map once it is learnt.

## 4. What the Foundations Must Gain

physis `discretisation`: Kirchhoff networks on graphs whose conductances adapt by local rules. The same structure serves *Physarum* ([`living`](../living/)) and crossbars ([`memristive`](../memristive/)).

## 5. References

- Stern et al. (2021); Scellier and Bengio (2017); Dillavou et al. (2022): in the theory note.
- Pashine, N., Hexner, D., Liu, A. J. and Nagel, S. R. (2019), Directed aging, memory, and nature's greed, *Science Advances* 5, eaax4215.
