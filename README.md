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
<p align="center">
<a href="docs/gallery.md"><img src="assets/gallery/brain-10-alpha-sources.jpg" width="32%" alt="The alpha rhythm of closed eyes estimated on the cortex from 64 EEG channels, in the Brain Atlas"></a>
<a href="docs/gallery.md"><img src="assets/gallery/brain-13-eeg-and-meg.jpg" width="32%" alt="One current in the cortex with its EEG potential and MEG field, in the Brain Atlas"></a>
<a href="docs/gallery.md"><img src="assets/gallery/brain-15-fmri-language.jpg" width="32%" alt="An fMRI map of listening to stories on the cortex, in the Brain Atlas"></a>
</p>
<p align="center"><i>The Brain Atlas: every part of the brain clickable for its anatomy, with sections, atlases, data, the head and its electrodes, and the mouse; and EEG, MEG and fMRI carried onto the brain.</i></p>

<p align="center">
<a href="docs/gallery.md"><img src="assets/gallery/snn-02-fast-oscillation.jpg" width="32%" alt="The fast oscillation of Brunel's network, spikes in narrow volleys, in the Spiking Simulator"></a>
<a href="docs/gallery.md"><img src="assets/gallery/snn-06-classes-of-excitability.jpg" width="32%" alt="f-I curves of three neuron models, two rising continuously and one jumping, in the Spiking Simulator"></a>
<a href="docs/gallery.md"><img src="assets/gallery/snn-07-quantum-emitters.jpg" width="32%" alt="Photons of quantum emitters as spikes, in the Spiking Simulator"></a>
</p>
<p align="center"><i>The Spiking Simulator: networks of classical and quantum spiking neurons, run by the neuroscience domain's C++ engine and shown with their statistics and the engine's profile.</i></p>

## Quick Start

**Online, with nothing to install.** The three atlases run on the project's server; the fourth application does not:

- **Earth Atlas:** http://173.212.213.207:9201/ — the atmosphere, the ocean, the land, the people and the solid Earth on one globe.
- **Brain Atlas:** http://173.212.213.207:9202/ — the human brain in the round, every part with its anatomy; sections, atlases, data, the head and its electrodes, the mouse brain, and EEG, MEG and fMRI mapped onto the brain.
- **Celestial Atlas:** http://173.212.213.207:9200/ — the Sun, the Moon, the planets, the stars and the galaxies, at any instant from 1849 to 2150.
- **Spiking Simulator:** not run on our online server, whose resources are too limited for it — networks of spiking neurons, classical and quantum, composed in the page and run by the neuroscience domain's C++ engine. Each run is a simulation that takes seconds on every core of the machine serving it, so it runs on your own computer: see [the Spiking Simulator on your computer](#the-spiking-simulator-on-your-computer), four lines.

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

<a id="the-spiking-simulator-on-your-computer"></a>**The Spiking Simulator on your computer.** It needs [Git](https://git-scm.com/), [Python 3.12 or later](https://www.python.org/downloads/) and a C++ compiler, `g++` (GCC 11 or later); no data and no Python packages. First the compiler, once:

- **Linux:** `sudo apt install g++` (or your distribution's GCC).
- **Windows:** WSL with g++ in it. In PowerShell as administrator, `wsl --install`, restart, open the *Ubuntu* app and run `sudo apt install g++`. The simulator finds it there by itself.
- **macOS:** untried; the `g++` of Apple's command-line tools (`xcode-select --install`) may serve.

Then, in a terminal (on macOS and Linux, write `python3`; skip the first two lines if you have the repository already):

```
git clone --recurse-submodules https://github.com/huangwanhong-maker/The-Science-Infrastructures.git
cd The-Science-Infrastructures
python data/acquisition/recipe.py applications/spiking-simulator/recipe.json
python applications/spiking-simulator/serve.py --open
```

The third line builds the engine (a minute or two, once; the page offers to do it too), and the fourth opens the simulator at http://127.0.0.1:8830/. Choose one of the eight examples (Brunel's balanced network, the classes of excitability, quantum emitters, STDP and more) and press **Run**; its activity, statistics and the engine's profile follow in the tabs. Its [README](applications/spiking-simulator/README.md) tells the rest.

*If something goes wrong.*

- **A data step failed** (a data server was down, or a network blocked it). Everything else still works: open the page, click **Data** at the top right, and fetch the missing steps there with one click, now or later. Running the recipe again does the same; it skips what is already there, and `--status` lists the steps.
- **`python` is not found,** or Windows opens the Store: install Python from python.org with *Add python.exe to PATH*, open a new terminal, and start again at the third line.
- **`foundations/` is empty:** the clone was made without submodules. Run `git submodule update --init`.
- **The Spiking Simulator says there is no engine:** it found no g++. Install one as above (on Windows, in WSL), then build again, from the page or with the recipe.
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

## 4. The Sciences and the Order of their Implementation

Every science is meant to become a domain here: a library of its own objects and models, placed in the order by what it may depend on, tested against published references, and drawn by the works that need it. The sciences are taken in the broad sense of the German *Wissenschaft*, with the social sciences and the humanities among them, for they too measure, model and infer, and they share structures with the rest (§1.1): a language is a grammar, a chronology is a reference frame, a market is a network, a melody is a signal.

This section lists every science once, and gives each a priority. *Status: the order of tier 1 was set by the user on 2026-10-23; the rest is proposed, and revisable.*

### 4.1 Foundations, Domains and Modules

Each discipline appears in the list exactly once, as one of three things.

- **A foundation** holds laws and methods that name no object in nature: mathematics, computation, physics, inference. A discipline made of such methods, such as logic, statistics, optics or control, is part of a foundation, not a domain.
- **A domain** holds the objects of one science and their models, in a library of its own. It is a *field*, or a *bridge* when it exists to join several fields beneath it, as tides join astronomy and oceanography. Its *height* is the length of the longest chain of domains beneath it.
- **A module** is a part of a domain. A discipline whose objects a domain already holds is a module of that domain, not a domain of its own: hydrology is a module of geography, volcanology of geophysics. A module may graduate into a domain when it grows large enough to stand alone (ADR 0026).

The branches under which the domains are grouped below, the Earth and the sky, matter, life and the rest, are headings for the reader. Like a programme, they own nothing.

### 4.2 The Tiers of Priority

| Tier | Meaning |
|---|---|
| **built** | Has code and checks against published references; every built domain is a `draft` |
| **1** | The next work, in the order the user set |
| **2** | The work the roadmap already plans for the built domains (Phases 1–3): the three planned bridges and the inference foundation |
| **3** | One step beyond: domains that stand directly on what is built or in tier 1, and serve a work that exists |
| **4** | The rest of the map, each when a work reaches it |

A priority such as 1.3 is the third place in tier 1. Within a tier, a domain comes after the domains it stands on. Between two candidates, the one comes first on which more later domains stand, which reuses more of what exists, or whose results can be checked against open references.

### 4.3 The Foundations

| Priority | Foundation | What it holds | Status |
|---|---|---|---|
| built | [**calx**](https://github.com/huangwanhong-maker/calx) | Computation: the retargetable stack, and numerical methods on floating-point numbers and arrays; with morphe, the theory of computation and information theory. Its C++ stack is roadmap Phase 0. | draft |
| 1.1 | [**morphe**](https://github.com/huangwanhong-maker/morphe), the algebraic foundation | Mathematical structures and their equivalences. Next: arithmetic (ℤ, ℚ, finite fields, polynomial rings); algebra (groups, rings, modules, representations, lattices, homological algebra); category theory (functors, natural transformations, limits, monoidal structure, algebras and coalgebras of functors); topology (simplicial and CW complexes, homology and cohomology, persistent homology, fundamental groups); geometry (algebraic geometry with Gröbner bases and varieties; differential geometry with physis); equivalence (canonical forms, bisimulation); formal structures (words, trees, automata, grammars, rewriting) | draft |
| 1.2 | [**physis**](https://github.com/huangwanhong-maker/physis), the foundation of theoretical physics | Geometry and the laws of physics, on morphe and calx. Next: units as a free abelian group; differential and Riemannian geometry and the rotation groups; Lagrangian and Hamiltonian mechanics; gravitation and general relativity; thermodynamics; continuum mechanics; radiative transfer; and, to be proposed in its own decision log, electromagnetism and optics, acoustics, plasma physics, quantum mechanics, statistical mechanics and field theory | draft |
| 2.1 | **inference** (name open) | Statistics and estimation, least squares and harmonic analysis, data assimilation, uncertainty, control, machine learning | planned |
| 3.8 | **the quantitative social methods** (name and boundary with inference open) | Measurement (scales, reliability, item response theory); samples and surveys; causal inference (potential outcomes, directed acyclic graphs, matching, instruments, differences in differences, regression discontinuity); panels, event histories and durations; networks (centrality, communities, diffusion); strategy and choice (games and equilibria, social choice and voting, mechanism design); agent-based simulation | planned |

### 4.4 The Domains

Fifty-five domains, nine of them built. A name links to the domain's directory where it exists; the others are placed in the order only when their work begins.

**The Earth and the sky.**

| Domain | Kind | Stands on | Scope | Priority |
|---|---|---|---|---|
| [metrology](domains/metrology/) | field | — | Time scales, reference frames, Earth orientation | built |
| [geodesy](domains/geodesy/) | field | metrology | The figure of the Earth, its gravity and rotation, reference systems; places, distances, directions | built |
| [astronomy](domains/astronomy/) | field | metrology, geodesy | Ephemerides, apparent places, the rotation of the bodies, galaxies; later astrophysics and cosmology | built |
| [geophysics](domains/geophysics/) | field | geodesy | Plate motion, earthquake mechanisms and aftershocks, the crust, sections | built |
| [geography](domains/geography/) | field | geodesy | Cartography and GIS, remote sensing, physical and human geography, terrain and drainage, zonal statistics | built |
| [oceanography](domains/oceanography/) | field | geodesy | The ocean column, waves, seawater, the sea floor, tsunami travel times | built |
| [meteorology](domains/meteorology/) | field | geodesy | Atmospheric dynamics, our own models, cyclones, fire weather, heat, settling of ash | built |
| [climate](domains/climate/) | bridge | meteorology, oceanography, geography, astronomy, air-sea | Climate classes; later energy balance and orbital forcing | built |
| [tides](domains/tides/) | bridge | astronomy, geodesy, oceanography, meteorology | The tide-generating potential, ocean and atmospheric tides, tidal analysis and prediction | 2.2 |
| [air-sea](domains/air-sea/) | bridge | meteorology, oceanography | Exchange of momentum, heat, water and gases at the sea surface | 2.3 |
| [navigation](domains/navigation/) | bridge | astronomy, geodesy, geography | Celestial navigation, charts and routes; later satellite positioning | 2.4 |
| cryosphere | field | geodesy | Ice sheets, glaciers, sea ice, snow, permafrost | 3.4 |
| geology | field | geophysics, chemistry | Minerals, rocks, strata, the geological time scale | 3.5 |
| planetary-science | bridge | astronomy, geophysics, meteorology | The planets and moons as worlds: surfaces, interiors, atmospheres | 3.6 |
| space-weather | bridge | astronomy, geophysics, meteorology | The Sun's activity, the magnetosphere, the ionosphere | 4 |

**Matter.**

| Domain | Kind | Stands on | Scope | Priority |
|---|---|---|---|---|
| molecular-science | field | — | Atoms, molecules, proteins and their interactions: molecular graphs, conformers, reaction networks and their dynamics | 1.6 |
| chemistry | field | molecular-science | Substances in bulk: equilibria, solutions, electrochemistry; inorganic, organic and analytical chemistry | 3.1 |
| condensed-matter | field | — | Crystals and their defects; electronic, magnetic and mechanical properties; materials | 4 |
| nuclear-physics | field | — | Nuclides, their decays and reactions; radiometric dating | 4 |
| particle-physics | field | — | The particles and their interactions; scattering | 4 |

**Life.**

| Domain | Kind | Stands on | Scope | Priority |
|---|---|---|---|---|
| [neuroscience](domains/neuroscience/) | field | — | The brain: reference spaces, parcellations and maps, the head and its electrodes; spiking neural networks next | built; its module of spiking networks 1.3 |
| genomics | field | — | Sequence alignment, k-mers, genome assembly on de Bruijn graphs, variant analysis | 1.5 |
| evolution | field | genomics | Phylogeny, population genetics, selection; the systematics of plants, animals and fungi | 3.2 |
| cell-biology | field | genomics, molecular-science | Genes and their expression, pathways, the cell | 3.3 |
| ecology | field | geography, climate, evolution | Populations, communities, ecosystems; biodiversity and biogeography | 3.7 |
| microbiology | field | genomics, cell-biology | Bacteria, archaea and viruses; infection | 4 |
| physiology | field | cell-biology, neuroscience | Organs and their regulation; immunology, endocrinology, nutrition | 4 |
| developmental-biology | field | cell-biology | How an embryo forms | 4 |

**Health.**

| Domain | Kind | Stands on | Scope | Priority |
|---|---|---|---|---|
| drug-discovery | bridge | genomics, molecular-science | Evidence joining targets and compounds, molecular properties, pharmacological models | 1.7 |
| epidemiology | field | demography | Disease in populations, its causes and its spread; public health | 3.10 |
| anatomy | field | — | The structures of the body, as images and as models | 4 |
| clinical-medicine | bridge | anatomy, physiology, drug-discovery, epidemiology | Diagnosis and treatment, under uncertainty | 4 |

**Mind and language.**

| Domain | Kind | Stands on | Scope | Priority |
|---|---|---|---|---|
| linguistics | field | — | Computational linguistics first: automata, context-free grammars, parsing, grammar inference, syntax, semantics and generative representations | 1.4 |
| psychology | field | neuroscience | Perception, learning, memory, emotion, personality | 4 |
| cognitive-science | bridge | neuroscience, psychology, linguistics | Minds as systems that compute | 4 |

**Society**, on the foundation of the quantitative social methods.

| Domain | Kind | Stands on | Scope | Priority |
|---|---|---|---|---|
| demography | field | — | Births, deaths, migration, populations by age | 3.9 |
| economics | field | — | Production, markets, money, growth, inequality | 4 |
| sociology | field | — | Groups, institutions, social networks, stratification | 4 |
| political-science | field | — | States, elections, conflict, policy | 4 |
| anthropology | field | evolution | Cultures, kinship, human variation | 4 |
| education | field | — | Learning and teaching, assessment | 4 |
| law | field | — | Legal systems, their rules and their reasoning | 4 |

**The humanities.**

| Domain | Kind | Stands on | Scope | Priority |
|---|---|---|---|---|
| history | field | metrology, astronomy | Events, sources and their dating; chronology and calendars | 4 |
| archaeology | bridge | history, geology | The material remains of the past, and their dating | 4 |
| philology | field | linguistics | Texts, their transmission and their style; literature; the classics | 4 |
| philosophy | field | — | Knowledge, reasoning and ethics; the philosophy of each science | 4 |
| musicology | field | — | Music: its scales, harmony and performance | 4 |
| art-history | field | — | Images: their colours, techniques and dating | 4 |
| religious-studies | field | philology, history | Traditions, their texts and their calendars | 4 |

**Engineering.**

| Domain | Kind | Stands on | Scope | Priority |
|---|---|---|---|---|
| electrical-engineering | field | — | Circuits, signals, power systems | 4 |
| mechanical-engineering | field | — | Machines and mechanisms; flight | 4 |
| civil-engineering | field | — | Structures, hydraulics, transport | 4 |
| chemical-engineering | field | chemistry | Processes and reactors at scale | 4 |
| robotics | field | mechanical-engineering, electrical-engineering | Systems that sense and act | 4 |
| energy-systems | bridge | meteorology, electrical-engineering, economics | Power from wind and sun, storage, grids | 4 |

### 4.5 Disciplines Held as Modules

These disciplines are not domains here. Each is held by the foundation or domain named, and may graduate into a domain of its own when it outgrows it.

| Discipline | Held by |
|---|---|
| Logic, the theory of computation, information theory | morphe and calx |
| Optics, acoustics, plasma physics, electromagnetism, quantum and statistical mechanics | physis |
| Statistics, data assimilation, control, machine learning and artificial intelligence | the inference foundation |
| Econometrics, psychometrics, survey methodology, game theory | the quantitative social methods |
| Cartography and GIS, remote sensing, hydrology, soils, human, economic and urban geography | geography |
| Seismology, volcanology, geomagnetism, gravimetry | geophysics |
| Astrophysics, cosmology | astronomy |
| Atmospheric chemistry, aerosols, atmospheric radiation | meteorology |
| Chemical and biological oceanography | oceanography |
| Palaeoclimate | climate |
| Mineralogy, petrology, stratigraphy, geochemistry, palaeontology | geology |
| Spiking neural networks (1.3), connectomics, the brain's signals | neuroscience |
| Structural biology, biochemistry, quantum chemistry of molecules | molecular-science |
| Pharmacology, until it graduates onto physiology | drug-discovery |
| Botany, zoology and mycology as systematics; population genetics | evolution |
| Biodiversity, biogeography, fisheries | ecology |
| Immunology, endocrinology, nutrition | physiology |
| Medical imaging, histology | anatomy |
| Public health | epidemiology |
| Computational linguistics (1.4), phonetics, typology, historical linguistics | linguistics |
| Chronology and calendars, beside metrology's time scales | history |
| Literature, stylometry, the classics | philology |
| Materials science and engineering | condensed-matter |
| Aerospace engineering | mechanical-engineering |

### 4.6 The Order of Implementation

```
first    the toolchain and the build (roadmap Phase 0)
tier 1   morphe → physis → spiking neural networks → linguistics → genomics → molecular-science
         → drug-discovery
tier 2   inference → tides → air-sea → navigation
tier 3   chemistry → evolution → cell-biology → cryosphere → geology → planetary-science → ecology
         → the quantitative social methods → demography → epidemiology
tier 4   the rest of §4.4, each when a work reaches it
```

The domains built and those of tiers 1 to 3 in the order, by height; a superscript is the tier, and a name without one is built:

```
height 5   ecology³
height 4   climate
height 3   tides² · navigation² · air-sea² · geology³ · planetary-science³
height 2   astronomy · geophysics · geography · oceanography · meteorology · cryosphere³
height 1   geodesy · drug-discovery¹ · chemistry³ · evolution³ · cell-biology³ · epidemiology³
height 0   metrology · neuroscience · linguistics¹ · genomics¹ · molecular-science¹ · demography³
```

*Rationale.*

- **Tier 1** is the user's order. Each step reuses what the one before made:
  - spiking networks are placed in the Brain Atlas's structures;
  - linguistics stands on morphe's formal structures;
  - genomics reuses the sequences, grammars and graphs;
  - molecular science stands on physis's chemistry and thermodynamics;
  - drug discovery joins genomics and molecular science.
- **Tier 2** finishes what the built domains began, as the roadmap plans:
  - the inference foundation and tides form the first vertical slice (Phase 1), with NOAA's predictions as the reference;
  - air–sea is the bridge that climate already declares (Phase 3);
  - navigation joins the sky of the Celestial Atlas to the charts of the Earth Atlas.
- **Tier 3** opens the next branches where the existing works already reach:
  - chemistry, for the chemistry of the ocean and the air, and for geology;
  - evolution and cell biology, on genomics;
  - the cryosphere, geology and planetary science, for sea ice, the rocks beneath the plates, and the bodies of the Celestial Atlas;
  - ecology, over the Earth Atlas's vegetation and climate classes;
  - the quantitative social methods, and on them demography and epidemiology, over the Earth Atlas's population grids and its measures of exposure.

## 5. Applications

- The [Celestial Atlas](applications/celestial-atlas/README.md) shows the Sun, the Moon, the planets, the stars and the galaxies, at any instant from 1849 to 2150, and the sky from any place on the Earth.
- The [Earth Atlas](applications/earth-atlas/README.md) shows the atmosphere, the ocean, the land, its people and the solid Earth on one globe:
  - forecasts from NOAA's models, from an AI model we run ourselves (Pangu-Weather) and from two simple models of our own;
  - tropical cyclones tracked in every forecast;
  - the ocean beneath any place, in depth;
  - the plates and their motion, 20,000 GPS stations and 70,000 earthquake mechanisms;
  - the crust, and sections through subduction zones.

- The [Brain Atlas](applications/brain-atlas/README.md) shows the human brain in the round, its cortex by gyrus and its deep nuclei, ventricles and fibre bundles drawn as surfaces, every part clickable for its anatomy: where it lies, what it does, what it connects, what supplies it and what follows when it fails. It cuts the brain and shows it in sections, colours the cortex by Julich-Brain's areas or Schaefer's networks, maps data by structure and by coordinate, places the electrodes of the 10–20 system on the scalp and carries their values onto the cortex, scales the brain to a subject's head, and shows the mouse brain in the Allen framework.

- The [Spiking Simulator](applications/spiking-simulator/README.md) composes networks of spiking neurons from twenty-two models, classical and quantum, and runs them in the neuroscience domain's C++ engine, exactly where the dynamics are linear and with the same spikes on any number of threads. Its examples include Brunel's balanced network in its four states, the classes of excitability, quantum emitters and STDP. It is not run on our online server, whose resources are too limited for it; it runs on your own computer, with g++ (in WSL on Windows), as [the Spiking Simulator on your computer](#the-spiking-simulator-on-your-computer) shows in four lines.

The three atlases run online on the project's server, the [Celestial Atlas](http://173.212.213.207:9200/), the [Earth Atlas](http://173.212.213.207:9201/) and the [Brain Atlas](http://173.212.213.207:9202/), and on your own computer in a few lines (Quick Start, above). The [index of applications](applications/README.md) lists those planned.

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
├── applications/    celestial-atlas · earth-atlas · brain-atlas · spiking-simulator
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

**As of 2026-10-23.** Four applications work: the three atlases, online and on any computer (Quick Start), and the Spiking Simulator ([ADR 0038](development/decisions/0038-the-spiking-simulator.md)), on any computer with g++. Their science is the provisional JavaScript of nine domains, of physis and of the visualization layer, tested against published references and closed-form answers in 422 checks. Each domain passes its checks on its own. Neuroscience, the first domain beyond the Earth and the sky, is a repository of its own attached as a submodule, with the Brain Atlas drawing it ([ADR 0028](development/decisions/0028-neuroscience-and-the-brain-atlas.md)). The Earth Atlas takes plug-ins: scripts for prediction, assimilation and analysis, loaded on request ([ADR 0027](development/decisions/0027-earth-atlas-plugins.md)). The first C++ is written: calx's and physis's first modules, on a [map of the mathematics](docs/mathematics.md) that places every area of mathematics in a foundation ([ADR 0035](development/decisions/0035-the-mathematical-foundations.md)), and the neuroscience domain's networks of spiking neurons, which simulate Brunel's network of 12,500 neurons and 15.6 million synapses at 2.6 times real time on a laptop, with the same spikes on any number of threads ([ADR 0036](development/decisions/0036-spiking-networks.md)); twenty-two neuron models, the first quantum one among them, are composed and run in the Spiking Simulator. The next steps are in the [roadmap](development/roadmap.md), in the order of §4.6:

- Phase 0, the toolchain and the build;
- tier 1: the foundations morphe, the algebraic foundation (category theory, algebra, topology, geometry), and physis, the foundation of theoretical physics; then spiking neural networks, computational linguistics, genomics, molecular science and drug discovery;
- tier 2: the inference foundation, and the three planned bridges, tides, air–sea and navigation;
- tier 3: chemistry, evolution, cell biology, the cryosphere, geology, planetary science, ecology, the quantitative social methods, demography and epidemiology.

## 8. License

The software is released under the MIT License. See [`LICENSE`](LICENSE). The foundation libraries carry the same licence, and datasets keep their own.
