# Kepler's Third Law Rediscovered from the Ephemeris

*Status: draft (2026-11-05).* Research: a case study of the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)), its module `programs`, on the astronomy domain's data.

## 1. The Question

Can a search over formulas recover Kepler's third law, T = 2π √(a³/GM), from the motions of the Solar System as JPL's ephemeris gives them, and what does typing the formulas by physical dimension change?

Kepler found T² ∝ a³ in 1619 from Tycho's tables, in units (the astronomical unit and the year) in which the constant is 1. That is the trap such a search falls into. In those units T = a^{3/2} fits the planets to 6 × 10⁻⁴, and it is not a law: it holds only because the Sun's GM is nearly 4π² au³ yr⁻². A formula is a law only if it survives a change of units. Buckingham's Π theorem states that condition for formulas whose symbols carry dimensions. heuresis types its programs by physis's dimensions, so the typed search can only write formulas that survive a change of units. The study measures what that buys on real data, and where real data depart from the two-body law.

## 2. Method

### 2.1 The Orbits from the Ephemeris

The program [kepler.cpp](kepler.cpp) reads DE440s through the astronomy domain's SPK reader (`astronomy::ephemerides::Spk`). It measures ten orbits over the kernel's span less a day at each end, 300.06 Julian years (1849-12-26 to 2150-01-21 TDB):

- **The planets** (the barycentres of Mercury's to Pluto's systems; the Earth as the Earth–Moon barycentre). Each is taken in **Jacobi coordinates**: the planet relative to the barycentre of the Sun and the planets inside its orbit, with μ = G(M_⊙ + Σ_inside mᵢ + m) (Murray and Dermott 1999, §6.2). For comparison each is also taken heliocentric, with μ = G(M_⊙ + m).
- **The Moon** about the Earth, with μ = G(M_E + M_M).
- **GM** are DE440's (Park et al. 2021), as the domain's `bodies.hpp` holds them.

For each orbit the program computes:

- **The mean semi-major axis**: the osculating a = 1/(2/r − v²/μ) at every sample (0.25 d for the Moon, 0.5 d to 30 d for the planets), averaged over the span.
- **The mean eccentricity**, from the eccentricity vector, averaged the same way.
- **The sidereal period, from the ephemeris alone.** The ecliptic longitude (the J2000 ecliptic, mean obliquity 84381.406″, IAU 2006) is unwrapped. Where the span holds 20 revolutions or more, the period is 2π over the slope of its least-squares line. Otherwise it is the mean time for the longitude to advance by 2π, from up to 400 starts spread over the span, each return found by calx's Brent method. The returns are biased where the apsides turn fast: the Moon's turn in 8.85 years, and its mean return falls 242 s short of its sidereal month. The line is biased where there are few revolutions: for Pluto, which completes 1.14 in the span, it is off by 5.7 %. Each method is used where it is not biased.

### 2.2 The Search

The cases (a, GM, e, T) are given in two systems of units: Kepler's own (the astronomical unit, 149 597 870 700 m by IAU 2012 Resolution B2, and the Julian year of 365.25 d), and SI. There are two problems:

- **Two centres**: the nine planets and the Moon. The variables are a (a length), GM (L³ T⁻²) and e (dimensionless), and the target is T (a time).
- **One centre**: the nine planets alone. The variables are a and e, with an ephemeral constant k of the sort of a gravitational parameter, drawn in [1, 100] au³ yr⁻² (its range converted with its unit). Without k, no typed term of the sort of a time could be formed at all.

The eccentricity is offered as a variable the law does not need. Without it, a and GM admit no dimensionless group, so every typed term of the sort of a time is C √(a³/GM), and the typed search would only choose C.

Each problem, in each system of units, is searched by `heuresis::programs::GeneticProgramming` twice: typed by dimension, and with every symbol dimensionless. In the untyped search k is a bare number in [1, 100] in either system, since nothing tells it that it is a GM. Both searches share the primitives (+, −, ×, ÷, square, cube, √), the constants (2, π, and an ephemeral c in [−1, 1]) and the options. Each runs 30 seeds of 30 000 evaluations, driven by `heuresis::search::run`. The objective is the RMS relative error of T over the cases. A term whose value is not finite is a `search::failure` with its reason.

### 2.3 The Test of a Law

A term is counted as the law when three things hold:

- **It has the law's form.** For two centres it equals 2π √(a³/GM) to 10⁻⁹ at six probe points that vary a, GM and e. For one centre it is proportional to a^{3/2} to 10⁻⁹ and does not depend on e.
- **It fits.** Its error on the data is below 1 %.
- **It survives a change of units.** Its error stays below 1 % after the units change (Kepler's to SI, or back), with its constants converted by their declared dimension. An untyped constant is a bare number and keeps its value.

### 2.4 The Front of Formulas

Every term evaluated is kept by size, and the best of each size over all 30 seeds is retained. The non-dominated set of (size, error), by `heuresis::multiobjective::non_dominated`, is the front of formulas.

## 3. Code Version

- The infrastructure at 4153cb9 (heuresis is a submodule there).
- heuresis 6b3cc88, calx 90aa710, morphe 14b3fb8, physis e2b98d8, tyche b4034b8.
- The astronomy domain as committed at 4153cb9.

To build and run (in WSL, g++ 11.5 and CMake 3.28), from this directory:

```
cmake -S . -B ~/.cache/science-build/exp-kepler -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-kepler -j4
~/.cache/science-build/exp-kepler/kepler ../../../data/cache/jpl-de440s/de440s.bsp results
.venv/Scripts/python.exe figures.py      # from Windows, the repository's environment
```

The 240 searches take 62 s on eight threads, and the orbits 3.4 s.

## 4. Parameters

| | |
|---|---|
| Population | 500, ramped half-and-half of depths 2 to 5 |
| Limits | depth 10, size 40 nodes; subtree mutation grows depth 3 |
| Variation | crossover 0.8, subtree mutation 0.1, point mutation 0.1 (rate 0.15) |
| Selection | tournaments of 7 under lexicographic parsimony; 2 elites |
| Dimension bound | exponents within ±4 |
| Budget | 30 000 evaluations a run; seeds 1000 to 1029 |

## 5. Data

| | |
|---|---|
| Source | JPL planetary ephemeris DE440s (Park, Folkner, Williams and Boggs 2021, AJ 161, 105), from NAIF |
| File | `data/cache/jpl-de440s/de440s.bsp`, sha256 c1c7feea…49f2, retrieved 2026-10-21 |
| Licence | NAIF's generic kernels are distributed by NASA/JPL without restriction |
| Constants | GM: DE440 (the astronomy domain's `bodies.hpp`); au: IAU 2012 B2; obliquity: IAU 2006 (Capitaine, Wallace and Chapront 2003) |

The kernel stays in the data library's cache. The tables in `results/` are this study's own: `orbits.csv` holds the ten orbits, `gp_runs.json` the 240 runs with their progress, `gp_fronts.json` the fronts with each formula's error in both systems of units, and `run.log` the program's report.

## 6. Results

### 6.1 The Law in the Ephemeris

| Body | Frame | a (au) | e | T (yr) | Method | T / 2π√(a³/μ) − 1 | Heliocentric |
|---|---|---|---|---|---|---|---|
| Mercury | Jacobi | 0.387098 | 0.2056 | 0.240847 | line | +1.4 × 10⁻⁶ | +1.4 × 10⁻⁶ |
| Venus | Jacobi | 0.723330 | 0.0068 | 0.615197 | line | +5.3 × 10⁻⁶ | +5.0 × 10⁻⁶ |
| Earth | Jacobi | 0.999998 | 0.0167 | 1.000018 | line | +3.9 × 10⁻⁶ | −1.3 × 10⁻⁶ |
| Mars | Jacobi | 1.523671 | 0.0934 | 1.880849 | line | +2.4 × 10⁻⁵ | +1.2 × 10⁻⁵ |
| Jupiter | Jacobi | 5.202784 | 0.0484 | 11.86246 | line | +4.9 × 10⁻⁵ | +3.7 × 10⁻⁵ |
| Saturn | Jacobi | 9.543280 | 0.0543 | 29.44624 | returns | −5.9 × 10⁻⁴ | −2.5 × 10⁻³ |
| Uranus | Jacobi | 19.19264 | 0.0473 | 84.01061 | returns | −2.2 × 10⁻⁴ | −2.7 × 10⁻³ |
| Neptune | Jacobi | 30.07446 | 0.0086 | 164.7947 | returns | −1.6 × 10⁻⁴ | −2.7 × 10⁻³ |
| Pluto | Jacobi | 39.48520 | 0.2490 | 247.9822 | returns | +1.2 × 10⁻⁴ | −2.5 × 10⁻³ |
| Moon | geocentric | 0.002563 | 0.0555 | 0.07480263 | line | +5.3 × 10⁻³ | |

What these residuals show:

- **The inner planets** keep the two-body law to a few parts in 10⁶.
- **The heliocentric axes of the outer planets** are 0.10 to 0.13 % longer than their Jacobi axes. Their T falls 2.5 × 10⁻³ short of the law's, because an outer planet circles the Sun together with the planets inside its orbit. Jacobi coordinates bring the residuals down to 1 to 6 × 10⁻⁴.
- **Saturn's remaining −5.9 × 10⁻⁴** is the scale of its great inequality with Jupiter: a near 5:2 resonance of period about 900 years, of which the span holds a third.
- **The Moon's period, 27.321662 d**, is the sidereal month, but it runs 0.53 % longer than the two-body law gives for its mean osculating axis. The Sun's tide lengthens it. That excess is about as large as the two-body correction (1 + m/M)^{−1/2} = −0.61 %, and opposite in sign. On these data a search could not tell Newton's form G(M + m) from Kepler's GM alone.

The law's own RMS relative error on the ten orbits is 1.68 × 10⁻³, the Moon's departure dominating. On the planets alone, with the Sun's GM in place of each μ, it is 6.1 × 10⁻⁴.

### 6.2 Typed and Untyped Searches

| Problem | Units | Search | Seeds meeting the law | Evaluations to it (median) | Final best is the law | Median final error | …after the change of units |
|---|---|---|---|---|---|---|---|
| two centres | au, yr | typed | 30 / 30 | 541 | 21 | 1.68 × 10⁻³ | 1.68 × 10⁻³ |
| two centres | au, yr | untyped | 0 / 30 | | 0 | 1.03 × 10⁻¹ | 1.7 × 10⁹ |
| two centres | SI | typed | 30 / 30 | 541 | 19 | 1.68 × 10⁻³ | 1.68 × 10⁻³ |
| two centres | SI | untyped | 2 / 30 | 8 254 | 0 | 5.24 × 10⁻¹ | 7.2 × 10⁶ |
| one centre | au, yr | typed | 30 / 30 | 803 | 27 | 4.41 × 10⁻⁴ | 4.41 × 10⁻⁴ |
| one centre | au, yr | untyped | 0 / 30 | | 0 | 6.01 × 10⁻⁴ | 1.8 × 10⁹ |
| one centre | SI | typed | 30 / 30 | 803 | 28 | 4.41 × 10⁻⁴ | 4.41 × 10⁻⁴ |
| one centre | SI | untyped | 0 / 30 | | 0 | 4.68 × 10⁻¹ | 1.3 × 10⁷ |

- **Typed by dimension, the law is met on every seed**: within 15 135 evaluations at worst for two centres, and 1 772 for one. Both systems of units meet it at the same evaluation on 59 of 60 seeds; the typed terms are covariant, and rounding breaks the ties that remain. The forms found are 2π √(a³/GM) written as (a + a)π √(a/GM), at 10 nodes, or as 2π a²/√(GM a), at 12. The canonical 9-node form 2π √(a³/GM) was never met.
- **Past the law, the typed front buys 4 % of error** with constants that nudge 2π to 2π(1 + 4 × 10⁻⁴). That fits the Moon's solar excess, which is overfitting. For two centres a formula past the law ended as the best on 9 of 30 seeds in Kepler's units, and on 11 in SI.
- **Untyped, in Kepler's units, the search finds the units.** a^{3/2} fits the planets at once. The Moon is patched with additive constants ((√a)³ + 0.0668, at 5 nodes) and with the eccentricity (e³/a + (√a)³). After the change to SI every one of these formulas is wrong by a factor of 10⁹.
- **Untyped, in SI, there is no such shortcut.** The search met the law on 2 seeds in 30, at 5 279 and 11 228 evaluations. On both seeds it went on to formulas fitting the ten cases better (errors 7.7 and 8.7 × 10⁻⁴), which are wrong by 10⁶ or more in au and years.
- **With one centre, typing makes the constant k a GM.** The front's law formulas imply GM = 39.508 to 39.528 au³ yr⁻² (a size-7 formula, 39.29). That is not the Sun's 39.4769 but lies within the planets' Jacobi μ, 39.4769 (Mercury) to 39.5299 (Neptune, Pluto): the outer planets are swung by the Sun and the planets inside them. The constant is found to about 10⁻³ only, because heuresis draws constants and does not fit them (§8).

### 6.3 The Fronts

In Kepler's units, two centres:

| Size | Typed | Error | Untyped | Error | In SI |
|---|---|---|---|---|---|
| 1 | | | a | 6.4 × 10⁻¹ | 3.7 × 10³ |
| 3 | | | (√a)³ | 3.2 × 10⁻¹ | 1.7 × 10⁹ |
| 5 | √(a³/GM) | 8.4 × 10⁻¹ | (√a)³ + 0.06675 | 1.0 × 10⁻¹ | 1.7 × 10⁹ |
| 8 | √(a/GM) a/0.1665 | 4.5 × 10⁻² | e³/a + (√a)³ | 4.4 × 10⁻² | 1.7 × 10⁹ |
| 10 | √(a/GM) π(a + a): **the law** | 1.68 × 10⁻³ | e³/(GM/e) + √(a³) | 2.0 × 10⁻² | 1.7 × 10⁹ |
| 22 | (a few constants on the law) | 1.62 × 10⁻³ | a formula of a, GM and e | 6.7 × 10⁻⁴ | 1.7 × 10⁹ |

The untyped front falls below the law's own error from size 22 on. On ten cases, an untyped formula of 20 nodes and five constants can fit the perturbations that the law leaves.

## 7. Figures

See [gallery.md](gallery.md):

- `search-10-kepler-law`: the law on log–log axes, and the residuals body by body.
- `search-11-kepler-front`: the fronts of formulas, typed and untyped.
- `search-12-kepler-success`: success and error over evaluations.
- `search-13-kepler-units`: the change of units, and the GM read from the one-centre front.

## 8. Caveats

- **Constants are drawn, not fitted.** `heuresis::programs` draws ephemeral constants uniformly and redraws them by point mutation (its README, open question 1). The one-centre GM is found to about 10⁻³, and the fronts carry constants of no meaning. Local fitting (calx's minimisers) or linear scaling would sharpen both; it is a gap of the library, not of this study.
- **Formulas are not simplified.** (a + a)π √(a/GM) is printed as the search found it (heuresis, open question 2). The figures fold constant subterms for reading and nothing more.
- **The axes are means of osculating axes.** For the Moon, a is ambiguous at the level of 0.1 to 0.4 % because of the Sun's perturbation. Its departure of 0.53 % is a property of that choice as much as of the Moon.
- **Saturn's period** is measured over a third of the great inequality, so its value depends on the window. Neptune's and Pluto's are measured over one or two revolutions.
- **Ten cases are few.** Both searches overfit them once the law is met. A larger set of orbits (satellites of the giant planets, asteroids) would separate the law from its corrections; DE440s holds only the planets and the Moon.
- **The success test is strict.** For two centres it asks for 2π exactly. A formula with C within 10⁻⁴ of 2π is counted as not the law, although it fits as well.
- **The typed search is easy here by construction.** Once a, GM and e are given, dimensional analysis leaves only functions of e times √(a³/GM). The typed search's speed measures how little remains to search, not the strength of the search.

## 9. References

- Barenblatt, G. I. (1996), *Scaling, Self-Similarity, and Intermediate Asymptotics*, Cambridge (Buckingham's Π theorem, ch. 1).
- Capitaine, N., Wallace, P. T. and Chapront, J. (2003), Expressions for IAU 2000 precession quantities, *A&A* 412, 567–586.
- IAU (2012), Resolution B2, on the re-definition of the astronomical unit of length.
- Koza, J. R. (1992), *Genetic Programming*, MIT Press.
- Montana, D. J. (1995), Strongly typed genetic programming, *Evol. Comput.* 3, 199–230.
- Murray, C. D. and Dermott, S. F. (1999), *Solar System Dynamics*, Cambridge (eqs. 2.133–2.134; §6.2, Jacobi coordinates).
- Park, R. S., Folkner, W. M., Williams, J. G. and Boggs, D. H. (2021), The JPL planetary and lunar ephemerides DE440 and DE441, *AJ* 161, 105.
- Schmidt, M. and Lipson, H. (2009), Distilling free-form natural laws from experimental data, *Science* 324, 81–85.
- Udrescu, S.-M. and Tegmark, M. (2020), AI Feynman: a physics-inspired method for symbolic regression, *Sci. Adv.* 6, eaay2631.
