# The Science Infrastructures

An integrative infrastructure for learning, using and researching all domains of science. Each domain, or each group of closely connected domains, is a laboratory. Shared libraries carry mathematics, computation and physics, and common conventions let the laboratories be read and combined as one body of knowledge.

<p align="center">
<a href="https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/blob/main/docs/gallery.md"><img src="https://raw.githubusercontent.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/main/assets/gallery/earth-20-plates.jpg" width="32%" alt="The tectonic plates in motion, in the Earth Atlas"></a>
<a href="https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/blob/main/docs/gallery.md"><img src="https://raw.githubusercontent.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/main/assets/gallery/earth-15-cyclone-forecasts.jpg" width="32%" alt="Forecasts of a typhoon from several models, in the Earth Atlas"></a>
<a href="https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/blob/main/docs/gallery.md"><img src="https://raw.githubusercontent.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/main/assets/gallery/celestial-02-saturn.jpg" width="32%" alt="Saturn, in the Celestial Atlas"></a>
</p>
<p align="center"><i>From the first laboratory: the Earth Atlas and the Celestial Atlas (<a href="https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/blob/main/docs/gallery.md">gallery</a>).</i></p>

## Abstract

**English.** This programme treats the sciences as one connected body of knowledge, to be learned, used and extended in theory and in practice. They share their foundations: mathematical structure (symmetry, algebra, geometry and topology), computation, physical law, and the logic of inference from data. They also meet in shared forms: fields on curved spaces, networks, sequences, populations and rates of change. The programme is organised as laboratories, one for each domain or group of connected domains. Each laboratory is its own repository, built on shared libraries and common conventions, and each can grow, change or be replaced without disturbing the others. The first laboratory studies astronomy, geography, oceanography and meteorology as one body of knowledge, with geometry as their foundation.

**中文.** 本计划把各门科学视为一个相互联系的知识整体，从理论与实践两方面加以学习、应用和拓展。它们共享基础——数学结构（对称、代数、几何与拓扑）、计算、物理定律以及从数据进行推断的逻辑——并在共同的形式中相遇：弯曲空间上的场、网络、序列、种群与变化率。本计划以"实验室"为单位组织，每个实验室对应一个领域或一组相互联系的领域，各自为独立的代码仓库，建立在共享的库与共同的约定之上，可以独立地成长、改变或被替换。第一个实验室以几何学为基础，把天文学、地理学、海洋学与气象学作为一个整体来研究。

**日本語.** 本プログラムは、諸科学を一つのつながった知識体系として捉え、理論と実践の両面から学び、使い、発展させる。諸科学は基盤を共有している。すなわち、数学的構造（対称性・代数・幾何・位相）、計算、物理法則、そしてデータからの推論の論理である。また、曲がった空間上の場、ネットワーク、系列、個体群、変化率といった共通の形式において出会う。本プログラムは「研究室」を単位として組織される。各研究室は一つの分野、または互いに結びついた分野の一群に対応する独立したリポジトリである。共有ライブラリと共通の規約の上に築かれ、他を乱すことなく成長し、変化し、置き換えられる。最初の研究室は、幾何学を基盤として天文学・地理学・海洋学・気象学を一体として研究する。

## 1. Motivation and Scope

The sciences are taught and practised apart, yet the same structures recur across them.

- **Symmetry.** Groups act on molecules and crystals, on the sphere of the sky, on quantum states and on sequences. Their representations give spectra, selection rules and harmonics.
- **Conservation and change.** The same partial differential equations transport heat in the ocean, chemicals in a cell and populations across a landscape.
- **Networks.** Food webs, neural circuits, reaction networks and economies are graphs with dynamics on them.
- **Sequences and grammars.** Genomes, languages and programs are strings with structure, studied by formal language theory.
- **Inference.** Every field estimates hidden states from noisy observations, and quantifies its uncertainty.

An infrastructure that holds the sciences together can make these connections explicit. It can also share one implementation of each method, test it once, and reuse it everywhere.

Every artefact serves one of three activities, and says which:

- **Learning:** derivations, worked examples and expositions, written to be understood;
- **Use:** tested tools whose numbers can be relied on;
- **Research:** open questions, experiments and new results, built on the first two.

## 2. Design Philosophy

The programme will grow over decades and change direction as it learns. Its governing principle is the one its first laboratory adopted: **modular, extensible, revisable**.

- **Laboratories as modules.** A laboratory states its scope and exposes its products: data, tools and notes. It depends on others only through the shared libraries and those published products.
- **Extension by addition.** A new domain is a new laboratory, and a new method is added to a shared library. Neither requires editing what already exists.
- **General before particular.** Mathematics, computation and physics are written once, in shared libraries, and used by every laboratory.
- **Revisable by record.** Decisions are recorded with their reasons, and are superseded, never rewritten. Every module and note carries a status: planned, seed, draft, stable or deprecated.
- **Verifiable.** Numbers are checked against authoritative references, and theory is checked by computation where it can be.

## 3. Architecture

| Layer | What it holds | Where | Status |
|---|---|---|---|
| Shared libraries | calx (computation), morphe (mathematical structures), physis (geometry and physics) | Today inside the first laboratory (`astronomy and geometry/lib/`); to be raised here (ADR 0001, open) | draft |
| Laboratories | One per domain or group of domains, each its own repository, attached here as a submodule | §4 | one working |
| Bridges | Topics that join laboratories, such as biogeochemistry or climate and health | `bridges/` | planned |
| This repository | The map of the sciences, the shared conventions, the decisions that cross laboratories, and the roadmap | `.` | draft |

Dependencies point one way: shared libraries ← laboratories ← bridges ← this repository. No library depends on a laboratory. No laboratory depends on this repository, or on another laboratory except through shared libraries and published products.

## 4. Laboratories

| Domains | Laboratory | Path | Status |
|---|---|---|---|
| Astronomy, geography, oceanography, meteorology, geodesy; geophysics proposed | [Integrative Astronomy, Geography, Oceanography Laboratory](https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory) | `astronomy and geometry/` | draft: two applications working, the C++ core planned |
| Mathematics, logic, computation, statistics | (the shared libraries for now) | | planned |
| Physics and chemistry; materials | | | planned |
| Life sciences: molecular and cell biology, genomics, ecology, evolution, neuroscience | | | planned |
| Health and medicine | | | planned |
| Social and behavioural sciences: demography, economics, linguistics, cognition | | | planned |
| Engineering and technology | | | planned |

Rows are a map of intent, not commitments. A laboratory is created when its first work begins, and its scope is settled in a decision record.

## 5. Organisation of the Repository

```
.
├── astronomy and geometry/   the first laboratory (git submodule)
├── development/              decisions · ideas · roadmap
├── CLAUDE.md                 working conventions for this repository
├── LICENSE                   MIT
└── README.md
```

Clone with every laboratory and library:

```
git clone --recurse-submodules https://github.com/huangwanhong-maker/The-Science-Infrastructures.git
```

## 6. Status

**Founded on 2026-10-22.**

- **This repository** holds the first laboratory as a submodule, and records the founding decision ([ADR 0001](development/decisions/0025-an-umbrella-of-laboratories.md)).
- **The first laboratory** has two applications:
  - the Celestial Atlas, which shows the Sun, the planets, the stars and the galaxies at any instant from 1849 to 2150;
  - the Earth Atlas, which shows the atmosphere, the ocean, the land, its people and the solid Earth on one globe, with forecasts from an AI model the laboratory runs itself.

  Both are tested against published references (66 and 151 checks).
- **Next:**
  - raise the shared libraries to this repository;
  - choose the second laboratory;
  - the first laboratory's C++ core.

  See the [roadmap](development/roadmap.md).

## 7. License

The software is released under the MIT License. See [`LICENSE`](LICENSE). Each laboratory and library carries the same licence, and datasets keep their own.
