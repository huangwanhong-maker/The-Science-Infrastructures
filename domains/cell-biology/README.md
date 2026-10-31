# cell-biology

*Status: planned. Field, height 1, on genomics and molecular science.* The cell as a system: its genes expressed and regulated, the pathways that carry its signals and its metabolism, its membranes and compartments, its cycle of growth and division, and the states it can be in. Systems biology and synthetic biology are its modules. Placed in the order with the other domains of life at the user's request (2026-10-31, [ADR 0061](../../development/decisions/0061-the-domains-of-life.md)); its work is tier 3.3 of the [README's order](../../README.md#46-the-order-of-implementation).

## 1. Purpose

- **A regulatory network is a dynamical system on a signed graph.** Genes are vertices and regulation is a signed edge. The dynamics are ordinary differential equations with Hill functions, Boolean maps on {0,1}ⁿ (Kauffman), or a chemical master equation when the molecules are few. The network's attractors are the cell's stable states. Thomas's rules, now theorems, tie the graph to the dynamics: multistationarity needs a positive circuit in the interaction graph, and a stable oscillation a negative one (Gouzé 1998; Soulé 2003). A positive loop of two repressions can make a switch (Gardner et al. 2000), and a negative loop of three a clock (Elowitz and Leibler 2000).
- **A pathway is a reaction network.** Its stoichiometric matrix S, its kinetics, its conservation laws (the left null space of S) and its steady fluxes (Sv = 0) are physis's `chemistry`; what names a kinase, a receptor or glucose is here. Flux balance analysis is a linear programme over the polytope {v : Sv = 0, l ≤ v ≤ u}. Elementary flux modes are the support-minimal elements of the flux cone {v : Sv = 0, v_i ≥ 0 for every irreversible i}, its extreme rays when every reaction is irreversible.
- **A cell is a point in a space of states.** Its expression profile, read one cell at a time, is a point of a space of many dimensions. Types are clusters, differentiation is a path, and RNA velocity is a vector field on the set the cells occupy (La Manno et al. 2018). The counting of transcripts stays in genomics (`expression`); the cell's states are here.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `regulation` | Gene regulatory networks: transcription factors bound at genomics's motifs; Hill functions; switches and oscillators (the toggle, the repressilator); Boolean networks and their attractors; stochastic expression (the telegraph model, bursts) on physis's master equations; Thomas's circuit rules checked on each network | planned |
| `signalling` | Signalling pathways as reaction networks: receptors and ligands, cascades and their ultrasensitivity (Huang and Ferrell 1996), robust adaptation (Barkai and Leibler 1997), network motifs (Milo et al. 2002) | planned |
| `metabolism` | Metabolic networks: stoichiometric matrices read from SBML, conservation laws, flux balance analysis, elementary flux modes, enzyme kinetics from molecular science | planned |
| `membranes` | Membranes and compartments: diffusion and transport between compartments; the Nernst and Goldman–Hodgkin–Katz potentials; osmosis and the cell's volume | planned |
| `cycle` | Growth and division: the cell cycle as a sequence of switches (Tyson and Novák), checkpoints, the control of size (sizer, timer, adder) | planned |
| `states` | Cells as points of expression space: types by clustering, trajectories and pseudotime, RNA velocity, Waddington's landscape where a potential exists | planned |
| `mechanics` | The cytoskeleton, molecular motors as ratchets, the shape and motility of cells (later) | planned |

## 3. Dependencies

| Source | What cell biology uses from it |
|---|---|
| genomics | Genes and their sequences, motifs and their position weight matrices (`annotation`), measured expression (`expression`) |
| molecular-science | The molecules of the pathways, enzyme kinetics, the structures of proteins |
| physis | `chemistry` (reaction networks, mass action, their exact Markov chains, reaction–diffusion), `statistical_mechanics` (master equations, entropy production), `dynamics` |
| calx | Integrators (the θ-method for stiff systems), linear algebra, optimisation (linear programming, to be written there for flux balance) |
| tyche | Regression; clustering and principal components (`multivariate`, planned) |
| morphe | Graphs and their canonical forms (motifs), signed graphs and their circuits |

## 4. Used By

Microbiology, physiology and developmental biology stand on it. No work uses it yet.

## 5. Open Questions

- Where genomics's `expression` ends and `states` begins: counting molecules per gene there, the cell as a point here. Counting single cells is probably genomics's.
- Whether Boolean networks and their attractors, finite dynamical systems that name nothing in nature, belong in physis's `dynamics` (discrete maps, planned) or in morphe; in either case not here.
- Whether SBML is read by the data library (an edge) or by the domain (from a path it is given).

## 6. References

- Alberts, B. et al. (2022), *Molecular Biology of the Cell*, 7th ed., W. W. Norton.
- Alon, U. (2019), *An Introduction to Systems Biology: Design Principles of Biological Circuits*, 2nd ed., CRC Press.
- Phillips, R., Kondev, J., Theriot, J. and Garcia, H. (2012), *Physical Biology of the Cell*, 2nd ed., Garland Science.
- Feinberg, M. (2019), *Foundations of Chemical Reaction Network Theory*, Springer.
- Kauffman, S. A. (1969), Metabolic stability and epigenesis in randomly constructed genetic nets, *Journal of Theoretical Biology* 22, 437–467.
- Gouzé, J.-L. (1998), Positive and negative circuits in dynamical systems, *Journal of Biological Systems* 6, 11–15.
- Soulé, C. (2003), Graphic requirements for multistationarity, *ComPlexUs* 1, 123–133.
- Gardner, T. S., Cantor, C. R. and Collins, J. J. (2000), Construction of a genetic toggle switch in *Escherichia coli*, *Nature* 403, 339–342.
- Elowitz, M. B. and Leibler, S. (2000), A synthetic oscillatory network of transcriptional regulators, *Nature* 403, 335–338.
- Huang, C.-Y. F. and Ferrell, J. E. (1996), Ultrasensitivity in the mitogen-activated protein kinase cascade, *PNAS* 93, 10078–10083.
- Barkai, N. and Leibler, S. (1997), Robustness in simple biochemical networks, *Nature* 387, 913–917.
- Milo, R. et al. (2002), Network motifs: simple building blocks of complex networks, *Science* 298, 824–827.
- Orth, J. D., Thiele, I. and Palsson, B. Ø. (2010), What is flux balance analysis?, *Nature Biotechnology* 28, 245–248.
- Schuster, S., Fell, D. A. and Dandekar, T. (2000), A general definition of metabolic pathways useful for systematic organization and analysis of complex metabolic networks, *Nature Biotechnology* 18, 326–332.
- Tyson, J. J., Chen, K. C. and Novák, B. (2003), Sniffers, buzzers, toggles and blinkers: dynamics of regulatory and signaling pathways in the cell, *Current Opinion in Cell Biology* 15, 221–231.
- Taheri-Araghi, S. et al. (2015), Cell-size control and homeostasis in bacteria, *Current Biology* 25, 385–391.
- La Manno, G. et al. (2018), RNA velocity of single cells, *Nature* 560, 494–498.
