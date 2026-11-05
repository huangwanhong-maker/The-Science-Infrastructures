# mechanical: Metamaterials, Oscillators, Soft Bodies and Granular Media

*Status: planned (ADR 0076).* Modes: programmed, reservoir, analogue. **Breadth:** multistable elements on a designed graph. **Depth:** granular criticality. Made of polymers, silicone, steel, sand and rice: abundant and reusable.

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| Mechanical metamaterials and hysterons (Bertoldi et al. 2017; Kwakernaak and van Hecke 2023) | multistability, hysteresis, memory, counting | polymers, silicone; reusable |
| Mechanical oscillators: metronomes (Pantaleone 2002; Martens et al. 2013) | synchronisation, chimeras | wood, steel, springs |
| Soft bodies (Hauser et al. 2011; Nakajima et al. 2013) | a body that computes as it moves | silicone, cables |
| Granular media and experimental landscapes (Frette et al. 1996; Paola et al. 2009) | avalanches, jamming, criticality in rice piles but rarely in sand; landscapes and strata | sand, rice, water |

## 2. The Simulators

1. **Interacting hysterons** (van Hecke 2021). Each hysteron switches up at $h_i^+$ and down at $h_i^-$, its thresholds shifted by the others' states, $h_i^\pm(S) = h_i^\pm - \sum_j c_{ij}s_j$. Their transition graph is traced under driving. Checks:
   - non-interacting hysterons have return-point memory (Preisach 1935), exactly;
   - a design that counts shows its period-$k$ orbit.
2. **Metronomes.** Small-angle pendula with an escapement of van der Pol type on a moving platform (Pantaleone 2002). Checks:
   - two metronomes on a light free platform synchronise in phase;
   - two platforms coupled by a spring show chimeras (Martens et al. 2013), qualitatively.
3. **The soft body.** A network of masses with nonlinear springs and dampers, its inputs forces and its readouts lengths (Hauser et al. 2011). Checks:
   - the echo state property: a negative conditional exponent;
   - a Volterra filter emulated better than by a linear readout of the input alone.
4. **Granular media.** The abelian sandpile and the Oslo model (physis `statistical_mechanics`), the rice pile's random thresholds. Checks:
   - the recurrent configurations number the spanning trees, the determinant of the reduced Laplacian (Dhar 1990), exactly;
   - the Oslo model's exponent of avalanche sizes, near 1.55 (Christensen et al. 1996).

   Experimental landscapes come later, after geography.

## 3. Tasks and their Lowering

- **Functions** (simulator 1): counting and sequences; the design is searched by heuresis.
- **Phenomena** (simulators 2 and 4).
- **Filters** (simulator 3).

## 4. What the Foundations Must Gain

- physis `statistical_mechanics`: sandpiles;
- physis `mechanics`: already written.

## 5. References

- Bertoldi et al. (2017); Kwakernaak and van Hecke (2023); Martens et al. (2013); Hauser et al. (2011); Nakajima et al. (2013); Frette et al. (1996); Paola et al. (2009): in the theory note.
- Bak, P., Tang, C. and Wiesenfeld, K. (1987), Self-organized criticality: an explanation of 1/f noise, *Physical Review Letters* 59, 381–384.
- Christensen, K., Corral, Á., Frette, V., Feder, J. and Jøssang, T. (1996), Tracer dispersion in a self-organized critical system, *Physical Review Letters* 77, 107–110.
- Dhar, D. (1990), Self-organized critical state of sandpile automaton models, *Physical Review Letters* 64, 1613–1616.
- Pantaleone, J. (2002), Synchronization of metronomes, *American Journal of Physics* 70, 992–1000.
- Preisach, F. (1935), Über die magnetische Nachwirkung, *Zeitschrift für Physik* 94, 277–302.
- van Hecke, M. (2021), Profusion of transition pathways for interacting hysterons, *Physical Review E* 104, 054608.
