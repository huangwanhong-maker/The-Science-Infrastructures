# The Science Infrastructures

An integrative infrastructure for learning, using and researching all domains of science, as one connected body of knowledge. Foundation libraries carry mathematics, computation and physics. Each science is a domain: a library of its own objects and models, pure and tested, that runs on its own with only the domains beneath it. Applications and studies, organised by purpose, compose the domains they need.

<p align="center">
<a href="docs/gallery.md"><img src="assets/gallery/earth-20-plates.jpg" width="32%" alt="The tectonic plates in motion, in the Earth Atlas"></a>
<a href="docs/gallery.md"><img src="assets/gallery/earth-15-cyclone-forecasts.jpg" width="32%" alt="Forecast tracks of a typhoon from several models, in the Earth Atlas"></a>
<a href="docs/gallery.md"><img src="assets/gallery/celestial-02-saturn.jpg" width="32%" alt="Saturn, in the Celestial Atlas"></a>
</p>
<p align="center"><i>The Earth Atlas and the Celestial Atlas, the first two applications: see the <a href="docs/gallery.md">gallery</a>.</i></p>

<p align="center">
<a href="docs/gallery.md"><img src="assets/gallery/brain-02-hippocampus.jpg" width="32%" alt="The deep structures of the brain through a glass cortex, with the hippocampus and its card, in the Brain Atlas"></a>
<a href="docs/gallery.md"><img src="assets/gallery/brain-03-coronal-section.jpg" width="32%" alt="A coronal cut through the brain, its face washed in the colours of the structures, in the Brain Atlas"></a>
<a href="docs/gallery.md"><img src="assets/gallery/brain-08-scalp-map.jpg" width="32%" alt="The electrodes of the 10-20 system on the scalp, with a map of channel values on the cortex, in the Brain Atlas"></a>
</p>
<p align="center">
<a href="docs/gallery.md"><img src="assets/gallery/brain-01-gyri.jpg" width="32%" alt="The cortex by gyrus and lobe, named as on an engraved plate, in the Brain Atlas"></a>
<a href="docs/gallery.md"><img src="assets/gallery/brain-07-networks.jpg" width="32%" alt="Yeo's seven networks on the cortex, in the Brain Atlas"></a>
<a href="docs/gallery.md"><img src="assets/gallery/brain-09-mouse.jpg" width="32%" alt="The Allen mouse brain and its major divisions, in the Brain Atlas"></a>
</p>
<p align="center"><i>The Brain Atlas: every part of the brain clickable for its anatomy, with sections, atlases, data, the head and its electrodes, and the mouse.</i></p>

## Quick Start

**Online, with nothing to install.** The three applications run on the project's server:

- **Earth Atlas:** http://173.212.213.207:9201/ — the atmosphere, the ocean, the land, the people and the solid Earth on one globe.
- **Brain Atlas:** http://173.212.213.207:9202/ — the human brain in the round, every part with its anatomy; sections, atlases, data, the head and its electrodes, and the mouse brain.
- **Celestial Atlas:** http://173.212.213.207:9200/ — the Sun, the Moon, the planets, the stars and the galaxies, at any instant from 1849 to 2150.

The Earth and Celestial Atlases link to each other: a place on the Earth leads to its sky, and a sky to its weather.

**On your own computer.** You need [Git](https://git-scm.com/), [Python 3.12 or later](https://www.python.org/downloads/) (on Windows, tick *Add python.exe to PATH* when installing) and a current browser. Nothing else: no Node, no build step, no account. Paste the block for your system into a terminal, in the folder where you keep code. The first run downloads about 450 MB of public data and takes ten to thirty minutes, most of it waiting for data servers.

Windows (PowerShell or Command Prompt):

```
git clone --recurse-submodules https://github.com/huangwanhong-maker/The-Science-Infrastructures.git
cd The-Science-Infrastructures
python -m venv .venv
.venv\Scripts\python -m pip install -r data\requirements.txt
.venv\Scripts\python data\acquisition\recipe.py applications\earth-atlas\recipe.json
.venv\Scripts\python data\acquisition\recipe.py applications\celestial-atlas\recipe.json
.venv\Scripts\python applications\earth-atlas\serve.py --open
```

macOS and Linux:

```
git clone --recurse-submodules https://github.com/huangwanhong-maker/The-Science-Infrastructures.git
cd The-Science-Infrastructures
python3 -m venv .venv
.venv/bin/python -m pip install -r data/requirements.txt
.venv/bin/python data/acquisition/recipe.py applications/earth-atlas/recipe.json
.venv/bin/python data/acquisition/recipe.py applications/celestial-atlas/recipe.json
.venv/bin/python applications/earth-atlas/serve.py --open
```

The Earth Atlas opens in your browser at http://127.0.0.1:8790/. For the Celestial Atlas, run the last line with `applications/celestial-atlas/serve.py` in a second terminal; it opens at http://127.0.0.1:8770/, and the two link to each other. For the Brain Atlas, fetch its data (about 80 MB) and serve it:

```
python data/acquisition/recipe.py applications/brain-atlas/recipe.json
python applications/brain-atlas/serve.py --open
```

It opens at http://127.0.0.1:8810/.

*If something goes wrong.*

- **A data step failed** (a data server was down, or a network blocked it). Everything else still works: open the page, click **Data** at the top right, and fetch the missing steps there with one click, now or later. Running the recipe again does the same; it skips what is already there, and `--status` lists the steps.
- **`python` is not found,** or Windows opens the Store: install Python from python.org with *Add python.exe to PATH*, open a new terminal, and start again at the third line.
- **`foundations/` is empty:** the clone was made without submodules. Run `git submodule update --init`.
- **The port is taken:** add `--port 8791` (or any free port) to the last line.
- **Serving other computers:** add `--host 0.0.0.0`, and `--sibling URL` with the other atlas's address so that the links between them work. Each atlas's README has the details.

The checks (`node development/tools/run-tests.mjs`) need Node 22, and nothing else does.

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
| Shared edges | Data (registry, acquisition, products, requests on demand), visualization, observation and sensors, and the API for other programs | `data/`, `visualization/`, `observation/`, `api/` |
| Domains | One library per science, pure and tested; fields, and bridges between them | `domains/` |
| Foundations | [calx](https://github.com/huangwanhong-maker/calx) (computation), [morphe](https://github.com/huangwanhong-maker/morphe) (mathematical structures), [physis](https://github.com/huangwanhong-maker/physis) (geometry and physics); inference, and the quantitative social sciences, planned (§4.3) | `foundations/` (submodules) |

The domains are ordered by what each may depend on:

```
height 4   climate
height 3   tides · navigation · air-sea
height 2   astronomy · geophysics · geography · oceanography · meteorology
height 1   geodesy
height 0   metrology · neuroscience
           ─────────────────────────────── foundations: physis ─► calx, morphe
```

Dependencies point down. No foundation depends on a domain; no domain depends on an edge or a work; nothing depends on an application. The full design is in [docs/architecture.md](docs/architecture.md), and its reasons in [ADR 0026](development/decisions/0026-domains-as-libraries.md).

## 4. Domains

Every science is meant to become a domain here: a library of its own objects and models, placed in the order by what it may depend on, tested against published references, and drawn by the works that need it. The sciences are taken in the broad sense of the German *Wissenschaft*, with the social sciences and the humanities among them, for they too measure, model and infer, and they share structures with the rest (§1.1): a language is a grammar, a chronology is a reference frame, a market is a network, a melody is a signal.

### 4.1 The Domains Built

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
| [neuroscience](domains/neuroscience/) | field | The brain: reference spaces, parcellations and maps; the head and its electrodes (a repository of its own, as a submodule) | draft |

### 4.2 The Next Domains

The next five, in this order, each reusing what the one before made. *Status: planned.*

| # | Domain or module | Kind and place | What it builds on | Its first content |
|---|---|---|---|---|
| 1 | **Spiking neural networks**, a module of [neuroscience](domains/neuroscience/) | module of a field, height 0 | The neuroscience domain and the Brain Atlas: neurons placed in its structures, activity shown on its surfaces | Leaky integrate-and-fire neurons, synapses and their delays, network simulation, spike-timing-dependent plasticity |
| 2 | **Computational linguistics** (`linguistics`) | field, height 0 | morphe's formal structures: languages, automata, grammars | Automata, context-free grammars, parsing, grammar inference; syntax, semantics and generative representations |
| 3 | **Genomics** (`genomics`) | field, height 0 | The sequence, grammar, graph and inference algorithms of the two before | Sequence alignment, k-mers, genome assembly on de Bruijn graphs, variant analysis |
| 4 | **Molecular science** (`molecular-science`) | field, height 0 | physis's chemistry and thermodynamics; graphs | Atoms, molecules, proteins and their interactions: molecular graphs, conformers, reaction networks and their dynamics |
| 5 | **Drug discovery** (`drug-discovery`) | bridge, height 1, over genomics and molecular science | Genomics, molecular science, molecular biology and pharmacology | Evidence joining targets and compounds, molecular properties, pharmacological models |

```
height 1   drug-discovery ─► genomics, molecular-science
height 0   linguistics · genomics · molecular-science · neuroscience (+ snn)
```

### 4.3 A Foundation for the Quantitative Social Sciences

Beside calx, morphe, physis and the planned inference foundation, a fifth foundation is planned for the methods that the social sciences share and that name no object in nature. *Status: planned; its name and its boundary with the inference foundation are open decisions.*

- **Measurement:** scales and indices, reliability and validity, item response theory.
- **Samples and surveys:** designs, weights, estimation for populations.
- **Causal inference:** potential outcomes, directed acyclic graphs, matching, instruments, differences in differences, regression discontinuity.
- **Panels and time:** panel models, event histories, durations.
- **Networks:** centrality, communities, diffusion on graphs.
- **Strategy and choice:** games and their equilibria, social choice and voting, mechanism design.
- **Simulation:** agent-based models of many interacting people.

### 4.4 The Map of the Sciences

A map of intent, not of commitments: each science that may join, by branch, with the structures through which it meets the rest. Each joins the order when the work reaches it. *Status of every row: planned,* except where a built domain already covers part of it.

**The formal sciences** are foundations, not domains: mathematics in morphe; computation and numerical methods in calx; geometry and the laws of physics in physis; statistics and inference in the planned inference foundation; the quantitative methods of the social sciences in the foundation of §4.3. Logic, the theory of computation and information theory belong with morphe and calx.

**The physical sciences.**

| Domain | Scope | Meets |
|---|---|---|
| Condensed matter and materials | Crystals, defects, electronic and mechanical properties | Space groups and their representations; band structure by Fourier analysis |
| Particle and nuclear physics | Particles, decays, nuclear structure and reactions | Lie groups and their representations; scattering; radioactive dating for geology and archaeology |
| Optics and photonics | Light in media, instruments, spectra | Fourier optics; the telescopes of astronomy; remote sensing |
| Plasma physics | Ionised gases, the Sun's corona, the ionosphere | Fluid equations shared with meteorology; space weather |
| Chemistry | Physical, inorganic, organic and analytical chemistry, beyond molecular science (§4.2) | Thermodynamics and kinetics; reaction networks; spectra |

**The Earth and planetary sciences**, beyond the domains built.

| Domain | Scope | Meets |
|---|---|---|
| Cryosphere | Ice sheets, glaciers, sea ice, permafrost | Climate, oceanography, geophysics (isostasy) |
| Hydrology | Rivers, groundwater, the water balance of catchments | Meteorology's precipitation; graphs of drainage; geography |
| Soil science | Soils, their formation and their carbon | Climate classes; chemistry; ecology |
| Volcanology | Eruptions, magma, volcanic hazards | Geophysics; the atmosphere's aerosols |
| Planetary science | The planets and moons as worlds: surfaces, interiors, atmospheres | Astronomy's ephemerides and rotation models; geophysics; meteorology |
| Space weather | The Sun's activity and the Earth's magnetosphere | Plasma physics; geophysics (the magnetic field); navigation (satellite positioning) |

**The life sciences**, with neuroscience built (§4.1) and genomics next (§4.2).

| Domain | Scope | Meets |
|---|---|---|
| Molecular and cell biology | Genes, proteins, pathways, the cell | Genomics; molecular science; reaction networks |
| Structural biology | The shapes of proteins and nucleic acids | Molecular science; the rotation group; inference from diffraction and microscopy |
| Developmental biology | How an embryo forms | Reaction–diffusion equations; growth on surfaces |
| Physiology | Organs and their regulation | Control and dynamical systems; neuroscience |
| Microbiology and immunology | Microbes, infection, immune responses | Genomics; population dynamics; epidemiology |
| Ecology | Populations, communities, ecosystems | Networks (food webs); the climate; geography |
| Evolution | Phylogeny, selection, population genetics | Trees and their inference; genomics; stochastic processes |
| Botany and zoology | Plants and animals, their forms and distributions | Geography's biomes; evolution |

**The health sciences.**

| Domain | Scope | Meets |
|---|---|---|
| Anatomy and medical imaging | The body's structures, as images and as models | The Brain Atlas's methods: templates, segmentations, surfaces |
| Epidemiology and public health | Disease in populations, its causes and its spread | Causal inference (§4.3); compartmental models; geography |
| Pharmacology | Drugs in the body: kinetics and effects | Drug discovery (§4.2); compartmental models |
| Clinical medicine and nutrition | Diagnosis, treatment, diet | Inference under uncertainty; the life sciences |

**The cognitive and behavioural sciences.**

| Domain | Scope | Meets |
|---|---|---|
| Psychology | Perception, learning, memory, emotion, personality | Measurement (§4.3); neuroscience |
| Cognitive science | Minds as systems that compute | Neuroscience; linguistics; learning theory |
| Linguistics | Beyond computational linguistics (§4.2): phonetics, typology, historical linguistics | Signals and spectra; trees of descent shared with evolution |

**The social sciences**, on the foundation of §4.3.

| Domain | Scope | Meets |
|---|---|---|
| Demography | Births, deaths, migration, populations by age | Geography's population grids; survival analysis; projections |
| Economics | Production, markets, money, growth, inequality | Games and equilibria; time series; networks of trade |
| Sociology | Groups, institutions, social networks, stratification | Networks; measurement; agent-based models |
| Political science | States, elections, conflict, policy | Social choice and voting; causal inference; geography's divisions |
| Anthropology | Cultures, kinship, human variation | Kinship as graphs and algebras; genomics |
| Human geography and urban studies | Cities, mobility, land use | Geography (built); networks; remote sensing |
| Education | Learning and teaching, assessment | Item response theory; the learning activity of the infrastructure itself |
| Law | Legal systems, rules and their reasoning | Formal logic; networks of citation |

**The humanities.**

| Domain | Scope | Meets |
|---|---|---|
| History and chronology | Events, sources and their dating | Calendars and time scales (metrology); astronomical dating of eclipses; inference from sources |
| Archaeology | Material remains of the past | Radiocarbon and dendrochronology (the climate's archives); geophysical survey; geography |
| Philology and literature | Texts, their transmission and their style | Grammars and sequences (linguistics); stemmata as trees; stylometry |
| Philosophy | Knowledge, reasoning, ethics, the philosophy of each science | Logic and formal semantics (morphe) |
| Musicology | Music, its scales, harmony and performance | Groups of pitch classes; spectra and acoustics |
| Art history and visual culture | Images, their colours, techniques and dating | Colour science; image analysis; spectra of pigments |
| Religious studies and classics | Traditions, their texts and their calendars | Philology; chronology |

**Engineering and technology.**

| Domain | Scope | Meets |
|---|---|---|
| Computing and artificial intelligence | Learning machines, their theory and their use | calx; the inference foundation; neuroscience's networks |
| Electrical and electronic engineering | Circuits, signals, power systems | Signals and networks; control |
| Mechanical, civil and aerospace engineering | Structures, machines, flight | physis's continuum mechanics; navigation |
| Chemical and materials engineering | Processes and materials at scale | Molecular science; thermodynamics |
| Energy systems | Generation, storage, grids | Meteorology (wind and sun); networks; economics |
| Control and robotics | Systems that sense and act | Dynamical systems; inference; navigation |

## 5. Applications

- The [Celestial Atlas](applications/celestial-atlas/README.md) shows the Sun, the Moon, the planets, the stars and the galaxies, at any instant from 1849 to 2150, and the sky from any place on the Earth.
- The [Earth Atlas](applications/earth-atlas/README.md) shows the atmosphere, the ocean, the land, its people and the solid Earth on one globe:
  - forecasts from NOAA's models, from an AI model we run ourselves (Pangu-Weather) and from two simple models of our own;
  - tropical cyclones tracked in every forecast;
  - the ocean beneath any place, in depth;
  - the plates and their motion, 20,000 GPS stations and 70,000 earthquake mechanisms;
  - the crust, and sections through subduction zones.

- The [Brain Atlas](applications/brain-atlas/README.md) shows the human brain in the round, its cortex by gyrus and its deep nuclei, ventricles and fibre bundles drawn as surfaces, every part clickable for its anatomy: where it lies, what it does, what it connects, what supplies it and what follows when it fails. It cuts the brain and shows it in sections, colours the cortex by Julich-Brain's areas or Schaefer's networks, maps data by structure and by coordinate, places the electrodes of the 10–20 system on the scalp and carries their values onto the cortex, scales the brain to a subject's head, and shows the mouse brain in the Allen framework.

The first two run online on the project's server, the [Celestial Atlas](http://173.212.213.207:9200/) and the [Earth Atlas](http://173.212.213.207:9201/), and all three on your own computer in a few lines (Quick Start, above). The [index of applications](applications/README.md) lists those planned.

## 6. Organisation of the Repository

```
.
├── foundations/     calx · morphe · physis (git submodules); inference planned
├── domains/         metrology · geodesy · astronomy · geophysics · geography · oceanography ·
│                    meteorology · tides · navigation · air-sea · climate · neuroscience (submodule)
├── data/            registry · acquisition · formats · products · requests on demand
├── visualization/   scene model · renderers · styles; provisional drawing helpers
├── observation/     instruments · field work · sensors
├── api/             the service and its clients in Python, C and C++, JavaScript
├── applications/    celestial-atlas · earth-atlas · brain-atlas
├── research/        questions · experiments · papers
├── notebooks/       exploration and teaching
├── docs/            architecture · gallery
├── development/     decisions · ideas · roadmap · tools
└── assets/          images · fonts · style tokens
```

To clone it with the foundations and the neuroscience domain:

```
git clone --recurse-submodules https://github.com/huangwanhong-maker/The-Science-Infrastructures.git
```

To run every check (the dependency rules, the suites of the foundations, domains, visualization and applications, and each domain alone):

```
node development/tools/run-tests.mjs
```

To run an application, follow the Quick Start at the top of this page; each application's README has the details, and their data come through the data library ([data/README.md](data/README.md)).

## 7. Status

**As of 2026-10-22.** Three applications work, two of them online, all on any computer (Quick Start). Their science is the provisional JavaScript of nine domains, of physis and of the visualization layer, tested against published references and closed-form answers in 380 checks. Each domain passes its checks on its own. Neuroscience, the first domain beyond the Earth and the sky, is a repository of its own attached as a submodule, with the Brain Atlas drawing it ([ADR 0028](development/decisions/0028-neuroscience-and-the-brain-atlas.md)). The Earth Atlas takes plug-ins: scripts for prediction, assimilation and analysis, loaded on request ([ADR 0027](development/decisions/0027-earth-atlas-plugins.md)). There is no C++ yet. The next steps are in the [roadmap](development/roadmap.md):

- Phase 0, the toolchain and the build;
- the next domains, in the order of §4.2: spiking neural networks, computational linguistics, genomics, molecular science, drug discovery; and the foundation for the quantitative social sciences (§4.3).

## 8. License

The software is released under the MIT License. See [`LICENSE`](LICENSE). The foundation libraries carry the same licence, and datasets keep their own.
