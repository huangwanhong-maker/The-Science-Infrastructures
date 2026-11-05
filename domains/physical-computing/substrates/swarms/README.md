# swarms: Robot Swarms

*Status: planned (ADR 0076).* Modes: programmed, analogue. **Breadth:** the network dynamics of a programmable unit; the rule each robot runs is a program, and the swarm is its many-body realisation. The robots are robotics's objects. When this family's simulator is written, the domain declares robotics, which gains a module `swarm` (ADR 0076 §2.5).

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| Kilobots and their kin (Rubenstein, Ahler and Nagpal 2012; Rubenstein, Cornejo and Nagpal 2014) | self-assembly of a thousand robots; gradients; consensus | batteries, plastics, electronics; reusable |
| Swarm morphogenesis (Slavkov et al. 2018) | Turing patterns grown by a swarm | as above |

## 2. The Simulators

1. **A swarm on a plane**, on robotics's `swarm`: units moved by vibration or by differential drive, exchanging infrared messages within a radius, each running a local rule. Checks:
   - a gradient formed from a seed equals the hop count on the communication graph (breadth-first search, morphe);
   - self-assembly as reported, qualitatively.
2. **Morphogenesis.** A reaction–diffusion system computed over the communication graph, each robot moving by its morphogen. Checks:
   - the unstable modes are those physis's `turing_modes` predicts on the graph's Laplacian;
   - the pattern's wavelength.

## 3. Tasks and their Lowering

- **Phenomena:** patterns.
- **Functions:** gradients, consensus.

A task's rule is compiled into each robot's program.

## 4. What the Foundations Must Gain

Nothing new: physis's `chemistry` (Turing's condition on graphs) and morphe's graphs are written. Robotics gains `swarm`.

## 5. References

- Rubenstein, M., Ahler, C. and Nagpal, R. (2012), Kilobot: a low cost scalable robot system for collective behaviors, *IEEE International Conference on Robotics and Automation*, 3293–3298.
- Rubenstein, M., Cornejo, A. and Nagpal, R. (2014), Programmable self-assembly in a thousand-robot swarm, *Science* 345, 795–799.
- Slavkov, I. et al. (2018), Morphogenesis in robot swarms, *Science Robotics* 3, eaau9178.
- Turing, A. M. (1952), The chemical basis of morphogenesis, *Philosophical Transactions of the Royal Society B* 237, 37–72.
