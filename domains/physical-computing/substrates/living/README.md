# living: Mycelium, *Physarum*, Microcosms, Gene Circuits, Biofilms and Neural Cultures

*Status: planned (ADR 0076).* Modes: reservoir, programmed, analogue, learning. **Depth:** the deepest family for populations and cells. These are live systems that have been seen to cross a Hopf bifurcation, show chaos and give early warnings before collapse. They grow themselves and biodegrade. They are slow, and they vary from one culture to the next.

## 1. Placement

Each substrate is this domain's as a device; the science of its organism belongs to its own domain (ADR 0076 §2.5).

- **Neural cultures and organoids** use neuroscience's `snn`, which is built; the domain declares neuroscience when that simulator is written.
- **The rest** hold the least model their simulators need, marked provisional, until their domain is founded. Microcosms and food webs go to ecology, gene circuits to cell-biology, biofilms and *Physarum* to microbiology, and the electrophysiology of mycelium to physiology.

Engineered organisms need containment. The ethics of computing with living neural tissue are unsettled (Smirnova et al. 2023), and it is an open question of the domain whether that tissue is in scope beyond its simulation.

## 2. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| Mycelium (Adamatzky 2022; LaRocco et al. 2025) | electrical spiking; mycelial memristors | grown, biodegradable |
| *Physarum* (Tero, Kobayashi and Nakagaki 2007; Tero et al. 2010) | adaptive transport networks; shortest paths | grown on oats |
| Ecological microcosms: chemostats, plankton, microbial food webs (Fussmann et al. 2000; Benincà et al. 2008; Becks et al. 2005; Dai et al. 2012; Carpenter et al. 2011) | a Hopf point crossed, chaos, critical slowing down before collapse | grows itself, biodegradable; days per run |
| Synthetic gene circuits (Elowitz and Leibler 2000; Gardner, Cantor and Collins 2000; Danino et al. 2010; Basu et al. 2005) | oscillators, switches, synchronised clocks, patterns | fed on sugar; containment |
| Biofilms (Prindle et al. 2015; Liu et al. 2015) | electrical waves of potassium; collective oscillation | as above |
| Neural cultures and organoids (Kagan et al. 2022; Cai et al. 2023) | spiking, plasticity, reservoirs | living tissue; ethics unsettled |

## 3. The Simulators

1. **Mycelium.** FitzHugh–Nagumo units on the graph of a mycelium, after Adamatzky's models, and the mycelial memristor as a volatile memristor (shared with [`soft-matter`](../soft-matter/)). Checks:
   - the units' threshold of excitability and their refractory time;
   - the speed of spikes on a chain against the travelling-wave solution;
   - the memristor's pinched loop.
2. ***Physarum*.** Tero's adaptation, $\dot D_e = |Q_e| - D_e$, on physis's Kirchhoff networks. Checks:
   - convergence to the shortest path where it is unique (Bonifaci, Mehlhorn and Varma 2012);
   - of two parallel paths, the shorter wins.
3. **The chemostat.** The Rosenzweig–MacArthur predator and prey, logistic and Holling II, and its form in a chemostat with a nutrient. Checks:
   - the Hopf point where the prey's isocline peaks, $x^* = (K - 1/(ah))/2$: the paradox of enrichment (Rosenzweig 1971);
   - tyche's early warnings on a noisy approach to a fold.
4. **Food chains.** Hastings and Powell's three species. Checks: the largest Lyapunov exponent is positive at their parameters ($b_1 = 3$), and the sum of the exponents is negative.
5. **Gene circuits.** Noise in each is intrinsic, by Gillespie's method (physis `chemistry`).
   - The repressilator, with its mRNA and protein. Check: it oscillates exactly when $(\beta + 1)^2/\beta < 3X^2/(4 + 2X)$ (Elowitz and Leibler 2000), in agreement with continuation's Hopf point.
   - The toggle switch. Check: with Hill exponents 2 and symmetric strengths, it is bistable exactly when $\alpha > 2$.
   - Repressilators coupled by quorum sensing (Garcia-Ojalvo, Elowitz and Strogatz 2004). Check: synchrony grows with the coupling.
6. **Biofilms.** The least model is a bistable chain, $u_t = Du_{xx} + u(1 - u)(u - a)$, for the potassium wave. Checks:
   - the front's speed $\sqrt{2D}\,(1/2 - a)$;
   - pinning on a coarse chain (physis's McKean).
7. **Neural cultures and organoids.** Networks of neuroscience's `snn`, stimulated through electrode arrays and read as rates. Checks:
   - the echo state property and separation;
   - memory capacity;
   - the reported tasks, qualitatively.

## 4. Tasks and their Lowering

- **Phenomena** (simulators 3 to 6).
- **Filters** (simulators 1 and 7).
- **Functions** (simulator 2): shortest paths.
- **Learning** (simulator 7): plasticity.

## 5. What the Foundations Must Gain

- physis `discretisation`: adaptive Kirchhoff networks;
- physis `dynamics`: continuation, Lyapunov spectra;
- tyche `time_series`: early warnings.

## 6. References

- Adamatzky (2022); Tero, Kobayashi and Nakagaki (2007); Bonifaci, Mehlhorn and Varma (2012); Fussmann et al. (2000); Benincà et al. (2008); Becks et al. (2005); Dai et al. (2012); Carpenter et al. (2011); Elowitz and Leibler (2000); Gardner, Cantor and Collins (2000); Danino et al. (2010); Basu et al. (2005); Prindle et al. (2015); Liu et al. (2015); Cai et al. (2023); Smirnova et al. (2023): in the theory note.
- FitzHugh, R. (1961), Impulses and physiological states in theoretical models of nerve membrane, *Biophysical Journal* 1, 445–466.
- Garcia-Ojalvo, J., Elowitz, M. B. and Strogatz, S. H. (2004), Modeling a synthetic multicellular clock: repressilators coupled by quorum sensing, *Proceedings of the National Academy of Sciences* 101, 10955–10960.
- Hastings, A. and Powell, T. (1991), Chaos in a three-species food chain, *Ecology* 72, 896–903.
- Kagan, B. J. et al. (2022), In vitro neurons learn and exhibit sentience when embodied in a simulated game-world, *Neuron* 110, 3952–3969.
- LaRocco, J. et al. (2025), Sustainable memristors from shiitake mycelium for high-frequency bioelectronics, *PLOS ONE*. To be checked.
- Nagumo, J., Arimoto, S. and Yoshizawa, S. (1962), An active pulse transmission line simulating nerve axon, *Proceedings of the IRE* 50, 2061–2070.
- Rosenzweig, M. L. (1971), Paradox of enrichment: destabilization of exploitation ecosystems in ecological time, *Science* 171, 385–387.
- Rosenzweig, M. L. and MacArthur, R. H. (1963), Graphical representation and stability conditions of predator–prey interactions, *American Naturalist* 97, 209–223.
- Tero, A. et al. (2010), Rules for biologically inspired adaptive network design, *Science* 327, 439–442.
