# developmental-biology

*Status: planned. Field, height 2, on cell biology.* How an organism forms from one cell: positions read from gradients, patterns that make themselves, clocks that segment a body, fates chosen and kept, tissues that grow and rearrange, and plants that grow by rewriting. Embryology, morphogenesis and regeneration are its modules. Placed in the order with the other domains of life at the user's request (2026-10-31, [ADR 0061](../../development/decisions/0061-the-domains-of-life.md)); its work is in tier 4 of the [README's order](../../README.md#46-the-order-of-implementation).

## 1. Purpose

- **Position is information.** A morphogen made at a source, diffusing with coefficient D and decaying at rate k, settles into the gradient c(x) = c₀ e^{−x/λ} with λ = √(D/k). A cell reads its position by thresholds on c (Wolpert's French flag); how precisely is bounded by the noise of the reading.
- **Patterns can make themselves.** Two species that react and diffuse at different rates can destabilise a uniform state (Turing 1952). physis's `chemistry` already holds Turing's condition, in closed form and on any graph Laplacian. This domain names the activators and inhibitors (Gierer and Meinhardt 1972) and the tissues they pattern.
- **A body is segmented by a clock and a wave.** The cells of the presomitic mesoderm oscillate by delayed negative feedback (Lewis 2003) while a wavefront of maturation sweeps past; each period freezes one segment (Cooke and Zeeman 1976).
- **A fate is an attractor.** Cell biology's regulatory networks have several stable states, and differentiation is a path among their basins: Waddington's landscape, where a potential exists. Lineage trees record which cell came from which.
- **Growth is rewriting.** A plant's architecture is a word rewritten in parallel: a Lindenmayer system (Lindenmayer 1968). A D0L system is a morphism h of the free monoid iterated from an axiom w; its growth function n ↦ |hⁿ(w)| = π Mⁿ η, with M the incidence matrix of h, has a rational generating series. An epithelium is a 2-complex of cells, and its vertex model an energy on the complex (Farhadifar et al. 2007). Its rearrangements (T1 transitions) and divisions are local moves on the complex.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `patterning` | Morphogen gradients by synthesis, diffusion and decay, their scaling and precision; positional information; Turing patterns of named systems (Gierer–Meinhardt) on physis's reaction–diffusion | planned |
| `clocks` | The segmentation clock: oscillators with delayed autorepression, their coupling and synchrony, the clock and wavefront, and the size of segments | planned |
| `fate` | Fates as attractors of regulatory networks from cell biology, the landscape, decisions under noise; lineage trees from recorded barcodes | planned |
| `tissues` | Epithelia as 2-complexes of cells: the vertex model's energy, T1 transitions and divisions as moves on the complex, the statistics of packing | planned |
| `growth` | Growth and form: L-systems and their growth functions, phyllotaxis and the golden angle (Vogel 1979), allometry | planned |

## 3. Dependencies

| Source | What developmental biology uses from it |
|---|---|
| cell-biology | Regulatory networks, signalling, the cell cycle, the states of cells |
| morphe | `formal` (rewriting, to be extended there to parallel rewriting for L-systems), `topology` (cell complexes and local moves on them) |
| physis | `chemistry` (reaction–diffusion, Turing's condition), `continuum` (diffusion), `mechanics` (the forces of the vertex model), `dynamics` (delay equations, planned there) |
| calx | Integrators, optimisation, linear algebra |

## 4. Used By

No work uses it yet.

## 5. Open Questions

- Lineage trees from barcodes are trees inferred from characters. Either developmental biology declares evolution and uses its phylogeny, or the inference of trees from characters, which names no organism, moves down to a foundation.
- Evolutionary developmental biology: a module here once this domain declares evolution, or a bridge over the two.
- The mechanics of tissues in three dimensions (folding, gastrulation), once physis's `continuum` has elasticity.

## 6. References

- Gilbert, S. F. and Barresi, M. J. F. (2016), *Developmental Biology*, 11th ed., Sinauer.
- Wolpert, L. (1969), Positional information and the spatial pattern of cellular differentiation, *Journal of Theoretical Biology* 25, 1–47.
- Turing, A. M. (1952), The chemical basis of morphogenesis, *Philosophical Transactions of the Royal Society of London B* 237, 37–72.
- Gierer, A. and Meinhardt, H. (1972), A theory of biological pattern formation, *Kybernetik* 12, 30–39.
- Meinhardt, H. (1982), *Models of Biological Pattern Formation*, Academic Press.
- Cooke, J. and Zeeman, E. C. (1976), A clock and wavefront model for control of the number of repeated structures during animal morphogenesis, *Journal of Theoretical Biology* 58, 455–476.
- Lewis, J. (2003), Autoinhibition with transcriptional delay: a mechanism for the zebrafish somitogenesis oscillator, *Current Biology* 13, 1398–1408.
- Waddington, C. H. (1957), *The Strategy of the Genes*, George Allen & Unwin.
- Lindenmayer, A. (1968), Mathematical models for cellular interactions in development, I and II, *Journal of Theoretical Biology* 18, 280–315.
- Rozenberg, G. and Salomaa, A. (1980), *The Mathematical Theory of L Systems*, Academic Press.
- Prusinkiewicz, P. and Lindenmayer, A. (1990), *The Algorithmic Beauty of Plants*, Springer.
- Vogel, H. (1979), A better way to construct the sunflower head, *Mathematical Biosciences* 44, 179–189.
- Farhadifar, R., Röper, J.-C., Aigouy, B., Eaton, S. and Jülicher, F. (2007), The influence of cell mechanics, cell-cell interactions, and proliferation on epithelial packing, *Current Biology* 17, 2095–2104.
