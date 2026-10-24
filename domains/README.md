# domains: The Sciences, One Library Each

*Status: draft.* The knowledge of the infrastructure, organised by subject ([ADR 0026](../development/decisions/0026-domains-as-libraries.md)).

## 1. Purpose

A domain holds one science's objects and their models: the planets and their ephemerides, the ocean column, the plates and their motion, the climates of the Earth. It is a library: pure, tested, and usable on its own.

The foundations below it hold laws and methods that name no object in nature; the applications and research above it are organised by purpose and draw on as many domains as they need. The test of placement is simple: content that names a thing in nature belongs to a domain.

## 2. The Order of the Domains

Let D be the set of domains, and let d → e mean that domain d may depend on domain e, as d's `module.json` declares. The reflexive and transitive closure ≤ of → is a partial order on D; the checker in [`development/tools/`](../development/tools/) verifies that → has no cycle, so that ≤ is antisymmetric.

- **A domain's library is its down-set.** ↓d = { e ∈ D : e ≤ d } is all that d needs, together with the foundations it declares; d builds and passes its tests with ↓d alone.
- **Every domain is integrative to some degree.** How integrative is its height: the length of the longest chain beneath it. Metrology has height 0; climate, with chains climate → meteorology → geodesy → metrology, has height 4.
- **A bridge is a domain with several fields beneath it.** Tides stand on astronomy, oceanography, geodesy and meteorology. "Bridge" is a label in `module.json`, not a separate kind of thing.
- **Works compose domains.** An application or a study uses a set S of domains, and so needs the order ideal ↓S = ⋃_{d∈S} ↓d. Two works may share domains freely: the domains are covered by the works, not partitioned among them.

```
height 0   metrology · neuroscience · history · genomics · economics
height 1   geodesy · evolution
height 2   astronomy · geophysics · geography · oceanography · meteorology
height 3   tides · navigation · air-sea                                  (bridges)
height 4   climate                                                       (bridge)
```

*Remark.* A domain may live in a repository of its own and be attached here as a git submodule (the graduation rule of ADR 0026); neuroscience is the first ([ADR 0028](../development/decisions/0028-neuroscience-and-the-brain-atlas.md)). The order and the checks do not see the difference.

## 3. Index

| Domain | Kind | Scope | May depend on | Status |
|---|---|---|---|---|
| [metrology](metrology/) | field | Time scales, reference frames, Earth orientation | — | draft |
| [neuroscience](neuroscience/) | field | The brain: reference spaces, parcellations and maps, the head and its electrodes; later surfaces, connectivity, neurons, signals (a submodule) | — | draft |
| [history](history/) | field | Records of histories by a codebook; chronology and calendars; the history complex, its topology and its signed ties as a ℤ₂ field, their spin foam computed exactly on physis; the history as a field theory (ADR 0041, 0042) | — | draft |
| [genomics](genomics/) | field | Genomes and their sequences: alphabets, formats and genetic codes; alignment as a semiring sum and its fast schedules; indexes and sketches; k-mers and assembly; mapping and variants; HMMs and annotation; multiple alignment; expression; simulation with its truth (ADR 0043) | — | draft |
| [economics](economics/) | field | Value carried between positions and through time: transfer networks, their loop deficits and prices; the circulation of flows on a network with faces; production as a circular flow; the spacetime complex of a record of transactions; places and the faces of a map (ADR 0044) | — | draft |
| [geodesy](geodesy/) | field | The figure of the Earth, its gravity field and rotation, reference systems; places, distances, directions | metrology | draft |
| [evolution](evolution/) | field | Phylogeny (distances, trees, likelihood, parsimony, the bootstrap, dating), population genetics (the coalescent, Wright–Fisher, diversity, linkage, structure), selection (dN/dS, McDonald–Kreitman) (ADR 0043) | genomics | draft |
| [astronomy](astronomy/) | field | Positional astronomy, ephemerides, the rotation of the bodies, galaxies; later astrophysics and cosmology | metrology, geodesy | draft |
| [geophysics](geophysics/) | field | The solid Earth: plate motion, earthquake mechanisms, the crust, sections; later seismology, gravity and magnetism | geodesy | draft |
| [geography](geography/) | field | All branches: physical geography, cartography and GIS, remote sensing, human geography | geodesy | draft |
| [oceanography](oceanography/) | field | Physical, chemical and biological oceanography; seawater; waves; the sea floor | geodesy | draft |
| [meteorology](meteorology/) | field | Atmospheric science: dynamics, moist processes, radiation, chemistry, prediction | geodesy | draft |
| [tides](tides/) | bridge | The tide-generating potential, ocean and atmospheric tides, tidal analysis and prediction | astronomy, geodesy, oceanography, meteorology | planned |
| [navigation](navigation/) | bridge | Celestial navigation, charts and routes; later GNSS | astronomy, geodesy, geography | planned |
| [air-sea](air-sea/) | bridge | Exchange of momentum, heat, water and gases at the sea surface; wind waves | meteorology, oceanography | planned |
| [climate](climate/) | bridge | The coupled system, energy balance, climatic zones, variability, orbital forcing | meteorology, oceanography, geography, astronomy, air-sea | draft |

Future domains join the same order, each beside or above the domains it stands on. The root README lists all fifty-five, with what each would stand on and its priority ([§4.4](../README.md#44-the-domains)), and the disciplines held as modules of them ([§4.5](../README.md#45-disciplines-held-as-modules)). Neuroscience, the first beyond the Earth and the sky, stands beside metrology; history, the first of the humanities, beside both ([ADR 0041](../development/decisions/0041-the-history-domain.md)); genomics, the first of the life sciences' molecular half, beside them, and evolution above it ([ADR 0043](../development/decisions/0043-the-genome-infrastructure.md)).

## 4. The Contract of a Domain

Each domain directory holds:

| Path | Content |
|---|---|
| `README.md` | Purpose, status, scope by module, dependencies, references and open questions |
| `module.json` | `name`; `kind` (`field` or `bridge`); `status`; `depends`, the domains it may use; `foundations`, the foundation libraries it may use |
| `js/` | Provisional JavaScript, until the C++ exists (ADR 0003, 0011). `js/index.js` is the public interface; nothing outside the domain imports anything else. |
| `tests/` | Its checks, run by `node domains/<name>/tests/run-tests.mjs`, with their reference values in `tests/fixtures/`, each with its source and retrieval date |
| `theory/` | Its notes (§5) |
| later | The C++: headers, sources and their tests, after [architecture §6](../docs/architecture.md) |

The rules:

- **Pure.** A domain has no network access, no instrument I/O, no user interface and no hidden global state, and reads only files whose paths it is given. Data reach it through `data/` and an application.
- **Through interfaces.** A domain imports another domain, or a foundation, only through its public interface, and only if its `module.json` declares it.
- **Tests stay within the down-set.** A domain's tests import what the domain may import, the shared harness (`development/tools/checks.mjs`) and Node's built-in modules. A check that needs two domains that do not stand on each other, or an application's code, belongs to that application.
- **Real data are optional.** Checks against the data library's products are skipped, and say so, when a product is missing.

All the rules are checked by `node development/tools/check-dependencies.mjs`, and all the checks run by `node development/tools/run-tests.mjs`.

## 5. Notes

A domain's notes are expositions, derivations and proofs, written to be understood. Notes on mathematics, numerics and general physics live in the foundation that implements them.

Notes are written in Markdown with LaTeX mathematics; the format is provisional. Each note begins with front matter:

```yaml
---
title: The Equilibrium Tide
status: draft              # planned | seed | draft | stable | deprecated
implemented-by:            # code that implements results from this note
  - domains/tides/js/potential.js
references: [Pugh2014]     # keys in docs/references.bib
---
```

When writing a note:

- state the hypotheses;
- separate definitions from results;
- number equations, so that code can cite them;
- follow the glossary.

## 6. Reference Models and Systems

Mature, open models and systems show how physics, data and observation are combined in practice.

- **Ocean:** the models MITgcm, MOM6, NEMO, ROMS, FESOM and Oceananigans.jl, and the ECCO state estimate.
- **Meteorology and the coupled system:** WRF, MPAS and ICON; the Earth-system models CESM and E3SM; the ERA5 reanalysis; and machine-learning weather models such as GraphCast and ECMWF's AIFS.
- **Geography:**
  - PROJ and GDAL (permissive licences), which may be linked;
  - QGIS and GRASS GIS (GPL), which we study but do not link.

We study these systems, reimplement their algorithms in this framework where useful, and validate against them.

- Reimplementing an algorithm is always allowed.
- Copying code requires a compatible licence.

See also the [idea on learning from ocean models](../development/ideas/2026-10-21-learning-from-ocean-models.md).
