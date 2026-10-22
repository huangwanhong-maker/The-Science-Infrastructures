# The Science Infrastructures

An integrative infrastructure for learning, using and researching all domains of science, as one connected body of knowledge. Foundation libraries carry mathematics, computation and physics. Each science is a domain: a library of its own objects and models, pure and tested, that runs on its own with only the domains beneath it. Applications and studies, organised by purpose, compose the domains they need.

<p align="center">
<a href="docs/gallery.md"><img src="assets/gallery/earth-20-plates.jpg" width="32%" alt="The tectonic plates in motion, in the Earth Atlas"></a>
<a href="docs/gallery.md"><img src="assets/gallery/earth-15-cyclone-forecasts.jpg" width="32%" alt="Forecast tracks of a typhoon from several models, in the Earth Atlas"></a>
<a href="docs/gallery.md"><img src="assets/gallery/celestial-02-saturn.jpg" width="32%" alt="Saturn, in the Celestial Atlas"></a>
</p>
<p align="center"><i>The Earth Atlas and the Celestial Atlas, the first two applications: see the <a href="docs/gallery.md">gallery</a>.</i></p>

## Abstract

**English.** This programme treats the sciences as one connected body of knowledge, to be learned, used and extended in theory and in practice. They share their foundations: mathematical structure (symmetry, algebra, geometry and topology), computation, physical law, and the logic of inference from data. They also meet in shared forms: fields on curved spaces, networks, sequences, populations and rates of change. The programme has three tiers. Foundation libraries hold mathematics, computation and physics. Each science is a domain: a library of its own objects and models, pure and tested, which runs on its own with only the domains beneath it. Applications and studies, organised by purpose, compose the domains they need. Each part can grow, change or be replaced without disturbing the others. The first domains are those of the Earth and the sky: metrology, geodesy, astronomy, geophysics, geography, oceanography and meteorology, with climate, tides, navigation and air–sea exchange as bridges between them.

**中文.** 本计划把各门科学视为一个相互联系的知识整体，从理论与实践两方面加以学习、应用和拓展。它们共享基础——数学结构（对称、代数、几何与拓扑）、计算、物理定律以及从数据进行推断的逻辑——并在共同的形式中相遇：弯曲空间上的场、网络、序列、种群与变化率。本计划分为三层：基础库承载数学、计算与物理；每门科学是一个领域，即一个关于其自身对象与模型的库，纯粹且经过检验，仅凭其下的领域即可独立运行；应用与研究按目的组织，组合各自所需的领域。各部分都可以成长、改变或被替换，而不影响其他部分。首批领域是关于地球与天空的：计量学、大地测量学、天文学、地球物理学、地理学、海洋学与气象学，并以气候、潮汐、导航与海气交换作为它们之间的桥梁。

**日本語.** 本プログラムは、諸科学を一つのつながった知識体系として捉え、理論と実践の両面から学び、使い、発展させる。諸科学は基盤を共有している。すなわち、数学的構造（対称性・代数・幾何・位相）、計算、物理法則、そしてデータからの推論の論理である。また、曲がった空間上の場、ネットワーク、系列、個体群、変化率といった共通の形式において出会う。本プログラムは三つの層からなる。基盤ライブラリが数学・計算・物理を担う。各科学は一つの分野であり、その対象とモデルを収めた純粋で検証済みのライブラリとして、その下にある分野だけで単独に動く。目的ごとに組織されたアプリケーションと研究が、必要な分野を組み合わせる。各部分は、他を乱すことなく成長し、変化し、置き換えられる。最初の分野は地球と空に関わるものである。計量学・測地学・天文学・地球物理学・地理学・海洋学・気象学であり、気候・潮汐・航法・大気海洋間の交換がそれらを結ぶ橋となる。

## 1. Motivation and Scope

### 1.1 Structures Shared by the Sciences

The sciences are taught and practised apart, yet the same structures recur across them.

- **Symmetry.** Groups act on molecules and crystals, on the sphere of the sky, on quantum states and on sequences. Their representations give spectra, selection rules and harmonics.
- **Conservation and change.** The same partial differential equations transport heat in the ocean, chemicals in a cell and populations across a landscape.
- **Networks.** Food webs, neural circuits, reaction networks and economies are graphs with dynamics on them.
- **Sequences and grammars.** Genomes, languages and programs are strings with structure, studied by formal language theory.
- **Inference.** Every field estimates hidden states from noisy observations, and quantifies its uncertainty.

An infrastructure that holds the sciences together can make these connections explicit. It can also share one implementation of each method, test it once, and reuse it everywhere.

### 1.2 The First Domains: the Earth and the Sky

The first domains share a geometry, and meet in phenomena.

- **Spheres and rotations.** The celestial sphere, the figure of the Earth, and the ocean and atmosphere, thin fluid shells on a rotating sphere, are all described in the language of the sphere and of SO(3). Changes of reference frame are elements of SO(3), and fields on the sphere expand in spherical harmonics, which carry its irreducible representations.
- **Tides join astronomy to oceanography.** The Moon and the Sun set up a tide-generating potential to which the ocean responds; each tidal constituent oscillates at an integer combination of six fundamental astronomical frequencies, so the tidal spectrum carries a lattice.
- **Geography maps the surface on which the others meet,** on geodesy and map projections; its climatic zones follow from how sunlight falls on the sphere.
- **Navigation joins astronomy to geography:** a position at sea is fixed from star sights by solving a spherical triangle.
- **Ocean and atmosphere form one coupled fluid system,** exchanging momentum, heat, water and gases at the sea surface, both driven by sunlight whose distribution astronomy sets.
- **The solid Earth moves beneath them all:** plates turning about their poles, earthquakes on their boundaries, the crust above the mantle.

### 1.3 Three Activities

Every artefact serves one of three activities, and says which:

- **Learning:** derivations, worked examples and expositions, written to be understood;
- **Use:** tested tools whose numbers can be relied on;
- **Research:** open questions, experiments and new results, built on the first two.

## 2. Design Philosophy

The programme will grow over decades and change direction as it learns. Its governing principle is **modular, extensible, revisable**.

- **Knowledge by subject, work by purpose.** A domain holds one science's objects and models; an application or a study composes as many domains as its purpose needs. Works cover the domains; none owns one.
- **Independent, and checked to be.** Each domain declares what it may depend on; the declared order is acyclic, and each domain passes its checks alone with what lies beneath it.
- **Extension by addition.** A new science is a new domain, and a new general method is added to a foundation. Neither requires editing what already exists.
- **General before particular.** Mathematics, computation and physics are written once, in the foundations, and used by every domain.
- **Revisable by record.** Decisions are recorded with their reasons, and are superseded, never rewritten. Every module and note carries a status: planned, seed, draft, stable or deprecated.
- **Verifiable.** Numbers are checked against authoritative references, and theory is checked by computation where it can be.

## 3. Architecture

| Tier | What it holds | Where |
|---|---|---|
| Works | Applications and studies, organised by purpose, each composing domains | `applications/`, `research/`, `notebooks/` |
| Shared edges | Data (registry, acquisition, products, requests on demand), visualization, observation | `data/`, `visualization/`, `observation/` |
| Domains | One library per science, pure and tested; fields, and bridges between them | `domains/` |
| Foundations | [calx](https://github.com/huangwanhong-maker/calx) (computation), [morphe](https://github.com/huangwanhong-maker/morphe) (mathematical structures), [physis](https://github.com/huangwanhong-maker/physis) (geometry and physics); inference planned | `foundations/` (submodules) |

The domains are ordered by what each may depend on:

```
height 4   climate
height 3   tides · navigation · air-sea
height 2   astronomy · geophysics · geography · oceanography · meteorology
height 1   geodesy
height 0   metrology
           ─────────────────────────────── foundations: physis ─► calx, morphe
```

Dependencies point down. No foundation depends on a domain; no domain depends on an edge or a work; nothing depends on an application. The full design is in [docs/architecture.md](docs/architecture.md), and its reasons in [ADR 0026](development/decisions/0026-domains-as-libraries.md).

## 4. Domains

| Domain | Kind | Scope | Status |
|---|---|---|---|
| [metrology](domains/metrology/) | field | Time scales, reference frames, Earth orientation | draft |
| [geodesy](domains/geodesy/) | field | The figure of the Earth, its gravity and rotation, reference systems | draft |
| [astronomy](domains/astronomy/) | field | Ephemerides, apparent places, the rotation of the bodies, galaxies | draft |
| [geophysics](domains/geophysics/) | field | Plate motion, earthquake mechanisms, the crust | draft |
| [geography](domains/geography/) | field | Cartography and GIS, physical and human geography, remote sensing | draft |
| [oceanography](domains/oceanography/) | field | The ocean column, waves, seawater, the sea floor | draft |
| [meteorology](domains/meteorology/) | field | Atmospheric dynamics, our own models, cyclones | draft |
| [climate](domains/climate/) | bridge | Climate classes; energy balance and orbital forcing planned | draft |
| [tides](domains/tides/), [navigation](domains/navigation/), [air-sea](domains/air-sea/) | bridges | Tides; celestial navigation and routes; exchange at the sea surface | planned |

Further domains join the same order as the work reaches them. As a map of intent, not commitments:

- chemistry and materials;
- molecular and cell biology, genomics, ecology, evolution and neuroscience;
- health and medicine;
- demography, economics, linguistics and cognition;
- engineering.

Mathematics, logic, computation and statistics are foundations.

## 5. Applications

- The [Celestial Atlas](applications/celestial-atlas/README.md) shows the Sun, the Moon, the planets, the stars and the galaxies, at any instant from 1849 to 2150, and the sky from any place on the Earth.
- The [Earth Atlas](applications/earth-atlas/README.md) shows the atmosphere, the ocean, the land, its people and the solid Earth on one globe:
  - forecasts from NOAA's models, from an AI model we run ourselves (Pangu-Weather) and from two simple models of our own;
  - tropical cyclones tracked in every forecast;
  - the ocean beneath any place, in depth;
  - the plates and their motion, 20,000 GPS stations and 70,000 earthquake mechanisms;
  - the crust, and sections through subduction zones.

The [index of applications](applications/README.md) lists those planned.

## 6. Organisation of the Repository

```
.
├── foundations/     calx · morphe · physis (git submodules); inference planned
├── domains/         metrology · geodesy · astronomy · geophysics · geography · oceanography ·
│                    meteorology · tides · navigation · air-sea · climate
├── data/            registry · acquisition · formats · products · requests on demand
├── visualization/   scene model · renderers · styles; provisional drawing helpers
├── observation/     instruments · field work
├── applications/    celestial-atlas · earth-atlas
├── research/        questions · experiments · papers
├── notebooks/       exploration and teaching
├── docs/            architecture · gallery
├── development/     decisions · ideas · roadmap · tools
└── assets/          images · fonts · style tokens
```

To clone it with the foundations:

```
git clone --recurse-submodules https://github.com/huangwanhong-maker/The-Science-Infrastructures.git
```

To run every check (the dependency rules, the suites of the foundations, domains, visualization and applications, and each domain alone):

```
node development/tools/run-tests.mjs
```

To run an application: `python applications/earth-atlas/serve.py --open`, or `python applications/celestial-atlas/serve.py --open`. Their data come through the data library ([data/README.md](data/README.md)).

## 7. Status

**As of 2026-10-22.** Two applications work. Their science is the provisional JavaScript of eight domains, of physis and of the visualization layer, tested against published references in 218 checks. Each domain passes its checks on its own. There is no C++ yet. The next steps are in the [roadmap](development/roadmap.md):

- Phase 0, the toolchain and the build;
- the choice of the next domains.

## 8. License

The software is released under the MIT License. See [`LICENSE`](LICENSE). The foundation libraries carry the same licence, and datasets keep their own.
