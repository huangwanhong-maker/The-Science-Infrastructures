# planetary-science

*Status: draft. Bridge, height 3.* The planets and moons as worlds: their interiors and crusts, their surfaces, their atmospheres, and the fields about them ([ADR 0052](../../development/decisions/0052-fields-known-in-part.md)). Founded at the user's request, earlier than its place in tier 3 of the order (README §4.6), for the Moon's magnetic field in [Cosmopoiesis](../../applications/cosmopoiesis-lab/). Its first module is the Moon's crustal magnetism and the interplanetary field about the Moon. Serves use and research.

## 1. Concepts

- **A body without a dynamo has a crustal field.** The Moon's field at its surface is the remanence of its crust, some nanotesla over most of it and hundreds over its strongest anomalies (Reiner Gamma; Gerasimovich, on the far side). It is known from orbit to a degree: T2015_449, the expansion of Tsunakawa et al. (2015)'s surface vector mapping of Kaguya's and Lunar Prospector's magnetometers ([theory §1](theory/magnetism.md)).
- **What a model does not resolve has a law.** The model's spectrum follows $A\,n^\beta\rho^{2n}$ over the degrees it resolves. The power the model lacks below that law, at its own highest degrees and above them, is drawn as a Gaussian field, labelled synthetic ([theory §2](theory/magnetism.md)).
- **An expansion is not continued downward.** Below its reference sphere, on the low ground of the maria, the model is read at the sphere above the mount ([theory §3](theory/magnetism.md)).
- **The solar wind brings a field.** The interplanetary field measured at L1 reaches the Moon about an hour later, except inside the Earth's magnetotail ([theory §4](theory/magnetism.md)).

The methods are the foundations': potential fields, their spectra and random fields are physis's `analysis`; the Legendre functions are calx's. This domain names the bodies, the models, their bands and their ground.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `magnetism` | The crustal field of a body without a dynamo: T2015_449 read, the law of its spectrum on degrees 30–250 and its effective depth, the unresolved part drawn from a seed, the reading on low ground; the interplanetary field carried from L1 to a body near the Earth, the magnetotail, GSE's axes; the Moon's crust as one of physis's fields (`LunarCrustalField`, `field.hpp`: the model or its unresolved part, low ground read at the sphere above with a note, the expansion shared; [ADR 0053](../../development/decisions/0053-fields-and-their-visual-expression.md)) ([`magnetism/`](magnetism/include/planetary_science/magnetism/)) | draft |
| `gravity` | The bodies' gravity fields (GRAIL's GRGM1200A for the Moon) on physis's expansions, with gravitation's normalisation | planned |
| `surfaces` | The bodies' figures and terrain beyond the data library's tiles: craters, slopes, illumination | planned |
| `interiors` | Density, pressure and temperature with depth; tides of solid bodies | planned |
| `atmospheres` | The atmospheres of Mars, Venus, Titan and the giants, on the methods meteorology shares through physis | planned |

## 3. Conventions

- **Frames.** The Moon's fields are in its mean-Earth frame (MOON_ME, the navigation domain's catalogue); positions are planetocentric, the latitude $\psi$ and the east longitude $\lambda$; a field's components are north, east and down in the local horizon.
- **Units.** Tesla for fields, metres for radii and heights, seconds for travel times (SI, CLAUDE.md §4.1).
- **Expansions.** Schmidt's semi-normalisation without the Condon–Shortley phase (physis decision 0010).

## 4. Dependencies

It may use astronomy, geophysics and meteorology (`module.json`). `magnetism` uses only the foundations: physis `analysis` and `fields`, and calx.

## 5. Checks

`node domains/planetary-science/tests/run-tests.mjs` builds the C++ checks with CMake and runs them on calx's verification (`magnetism/tests/check_magnetism.cpp`):

- the carriage from L1 against the geometry of the Sun–Earth line;
- GSE's axes, orthonormal and right-handed;
- the nearest sample of a series, and the local horizon;
- low ground, on a dipole;
- T2015_449 at five mounts on and above the Moon against pyshtools (`tests/fixtures/lunar-field.txt`, made by `tests/make_fixtures.py`);
- the law: bounded, and fitted to its spectrum;
- the unresolved part's power above degree 449 against the law, and its determinism;
- the crust as a field: pyshtools's mounts on the frame's axes, low ground read at the sphere above with its note, the base's bounds with their reasons, a sphere read by rows against point by point, and the composite of the model and its unresolved part, synthetic, with its parts.

The model comes through the data library (`python data/acquisition/fetch.py moon-magnetic-t2015`); without it those checks are skipped, and say so.

## 6. References

- Tsunakawa, H., Takahashi, F., Shimizu, H., Shibuya, H. and Matsushima, M. (2015), Surface vector mapping of magnetic anomalies over the Moon using Kaguya and Lunar Prospector observations, *Journal of Geophysical Research: Planets* 120, 1160–1185.
- Wieczorek, M. A. (2018), Strength, depth, and geometry of magnetic sources in the crust of the Moon from localized power spectrum analysis, *Journal of Geophysical Research: Planets* 123, 291–316; T2015_449, Zenodo, doi:10.5281/zenodo.3873648 (CC BY 4.0).
- Voorhies, C. V., Sabaka, T. J. and Purucker, M. (2002), On magnetic spectra of Earth and Mars, *Journal of Geophysical Research* 107, 5034.
- Smith, D. E. et al. (2010), Initial observations from the Lunar Orbiter Laser Altimeter (LOLA), *Geophysical Research Letters* 37, L18204.
- Mailyan, B., Munteanu, C. and Haaland, S. (2008), What is the best method to calculate the solar wind propagation delay?, *Annales Geophysicae* 26, 2383–2394.
- Lue, C. et al. (2011), Strong influence of lunar crustal fields on the solar wind flow, *Geophysical Research Letters* 38, L03202.

## 7. Open Questions

- **Space weather.** The carriage of the interplanetary field and the Earth's magnetotail belong to the space-weather bridge (README §4.4, tier 4) when it is founded; they are here until then.
- **Other models.** Ravat, Purucker and Olsen (2020)'s lunar model, from along-track gradients, could not be fetched (its host challenges automatic clients); Mars's crustal field (Langlais et al. 2019) is the next body's.
- **A non-stationary prior.** The unresolved part is drawn from a stationary, isotropic law; the Moon's anomalies cluster.
- **The lobes' field.** Inside the magnetotail the lobes' field (some 10 nT) is not given.
