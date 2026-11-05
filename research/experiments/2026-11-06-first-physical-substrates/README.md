# The First Physical Substrates

*Status: draft (2026-11-06). Learning and research.* One example of each simulator of the physical-computing domain's first three families, fluids, chemical and molecular ([ADR 0076](../../../development/decisions/0076-physical-computing.md)), run and drawn for the gallery's §15. The simulators and the laws beneath them (physis's decision 0022) are checked in their libraries; this record runs them on examples and draws them. It fits nothing.

## 1. The Question and Why

The user asked for the domain's "formalism and simulators", the fluids, chemical and molecular families first, "and example visualization". Each example shows one phenomenon the substrate can produce or compute, with the check that holds it to a closed form or a published result:

| Figure | Substrate | What it shows | Held to |
|---|---|---|---|
| 1 | convection cell | onset, steady rolls, chaos | Chandrasekhar's 1707.762, 1100.65, 27π⁴/4; Lorenz's r_H |
| 2 | rotating tank | baroclinic waves from noise | Phillips's growth rates; Θ_m = 16/(π² + (mL/R)²) |
| 3 | Faraday layer | the first tongue, onset against drive | 4γω₀/(k tanh kh) at the tip |
| 4 | walking droplet | walking past β = 2, orbits in a well | u = √(2β − 1 − √(1 + 4β))/√2 |
| 5 | wave tank | a reservoir's state; fading memory | e^(−γt) |
| 6 | capillary network | a maze solved by flow | the unique path, dead ends without flow |
| 7 | foam | coarsening | von Neumann's law, statistically |
| 8 | bubble logic; Phillips's machine | AND and OR in one junction; the multiplier | the truth table; I/(1 − c) |
| 9 | BZ medium, stirred | two Hopf points, relaxation oscillation | the trace and continuation agree |
| 10 | BZ medium, in a gel | a broken wave turning in a dish with no-flux walls | (a picture; excitability and pinning are checked on a chain) |
| 11 | BZ array | locking | the gap's spread below 1% of a period |
| 12 | reaction networks | Lorenz on dual rails; max with probability one | the ODE; exact absorption |
| 13 | electrodes | Kuramoto's transition | √(1 − 2γ/K) |
| 14 | strand displacement | convergence to the formal network | error ∝ 1/C_max |
| 15 | seesaw gates | OR and AND by thresholds | truth tables |
| 16 | PEN predator–prey | the paradox of enrichment; a front | K_H in closed form; Bramson's correction |

## 2. How to Run It

```
cmake -S . -B ~/.cache/science-build/exp-physical-substrates -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-physical-substrates -j4 --target substrates
~/.cache/science-build/exp-physical-substrates/substrates results            # all sections, about 5 minutes
~/.cache/science-build/exp-physical-substrates/substrates results spiral     # one section
.venv/Scripts/python.exe figures.py [--only 1,5]                              # from the repository's root
```

`substrates.cpp` writes one JSON file per section to `results/`; `figures.py` draws `assets/gallery/physcomp-01-…` to `physcomp-16-…` at 1600 × 1000.

## 3. What Was Found on the Way

Each of these was found by a check failing, and each was fixed at its cause, not by loosening the check:

1. **A tank never forgets its volume.** The uniform mode of a basin is neutral, so the echo state property holds for every departure of the surface but its mean. The motors' forcing was made volume-conserving, as a plunger's is.
2. **A periodic box is not a channel.** It admits modes uniform across the gap, which walls forbid and which grew fastest. physis's two-layer model now keeps a channel by images. The annulus's transition depends on its whole azimuthal wavenumbers, Θ_m = 16/(π² + (mL/R)²), not on the long-channel 16/π².
3. **A pole in the kinetics.** The Oregonator's (u − q)/(u + q) let an implicit stage's Newton iteration converge beyond u = −q. physis's split step now damps Newton, halves failing steps and keeps concentrations non-negative.
4. **Fronts must be resolved.** A pulse of the BZ medium propagates on lattices that resolve its front (spacing 0.1 against a width near 0.2), and is pinned on coarser ones: the first runs showed no propagation at a spacing of 0.25.
5. **A model that did not oscillate was not kept.** The PEN activator–inhibitor model first written had no Hopf point at any parameters tried. The predator–prey of Fujii and Rondelez, in Rosenzweig and MacArthur's form with its Hopf point in closed form, replaced it.
6. **Stiff cascades.** Strand displacement converges to its formal network only as fuel and cascade rates grow together; integrated implicitly, its error falls exactly as 1/C_max.
7. **A torus mimics a spiral.** A spiral's wavelength spans some T/ε front widths: at ε = 0.04 no affordable torus held one. At ε = 0.12 a plane wave circulated round the torus and looked sustained without turning. In a dish with no-flux walls (physis's `SparseLaplacian::grid`, the domain's `BZMedium::dish`) a plane wave leaves; the broken wave keeps turning about its free end to t = 40, its excitation passing from quadrant to quadrant. At ε = 0.16 the medium no longer carries a wave at all.

## 4. Limits

These are simulations of models, each with its declared limits (in its header and its family's README), not measurements. The Faraday thresholds count bulk damping only. The electrodes are phases only. The bubble gate's junction rule is a model. The seesaw rates are declared. The PEN parameters are scaled and illustrative. A spiral's picture is not a check.
