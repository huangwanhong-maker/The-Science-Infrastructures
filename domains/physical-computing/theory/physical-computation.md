---
title: Physical Computation — Definitions, Modes, Transfer, Cost and Generativity
status: draft
date: 2026-11-05
related: [ADR 0076, calx, physis (units, dynamics, control, statistical_mechanics, chemistry, continuum, discretisation, quantum), tyche (learning, time_series, information), heuresis, the neuroscience, robotics, ecology, microbiology, meteorology and oceanography domains, Cosmopoiesis]
---

# Physical Computation: Definitions, Modes, Transfer, Cost and Generativity

*Status: draft (2026-11-06). The theory note of the physical-computing domain ([ADR 0076](../../../development/decisions/0076-physical-computing.md)), promoted from the idea of 2026-11-05, which was written at the user's suggestion ("I think emerging computation tech worth a domain") after two assessments of candidate substrates that the user passed on, and with the user's three criteria: complex systems, a high capacity to generate phenomena, and kindness to the environment ("important for climate, Earth environment, ecosystem").*

*§1 is the domain's theory. §§2–6 keep the proposal as it was written. The decision amends two of them: §2.2, because a living substrate is now the domain's own as a device, its organism's science belonging to its domain (ADR 0076 §2.5); and §4, because the first work is now a simulator for every family of substrates, not a slice of three (ADR 0076 §2.7). The plan of each family is in [`../substrates/`](../substrates/).*

## 1. The Claim Made Precise

### 1.1 The Definition of a Physical Computer

The slogan "any dynamical system can compute" is true only vacuously unless the encoding and the reading are restricted. By Putnam's argument (1988), any open physical system, read through a suitably chosen map of its states, realises any finite automaton; the map then does the computing. A precise definition puts the restriction where it belongs (Horsman, Stepney, Wagner and Kendon 2014).

**Definition.** A *physical computer* is a tuple $(\mathcal{S}, \Phi_u, \mathcal{M}, C, E, R, \varepsilon)$: a substrate with states $\mathcal{S}$ and its evolution $\Phi_u : \mathcal{S} \to \mathcal{S}$ under a control $u$; an abstract model $C : \mathcal{M} \to \mathcal{M}$ with a metric $d$; an encoding $E : \mathcal{M} \to \mathcal{S}$ and a readout $R : \mathcal{S} \to \mathcal{M}$; and a tolerance $\varepsilon$. It *computes $C$ on $K \subseteq \mathcal{M}$* when the square commutes within $\varepsilon$:
$$ d\big(R(\Phi_u(E(m))),\ C(m)\big) \le \varepsilon \qquad \text{for all } m \in K. $$
It is *non-trivial* when $E$ and $R$ belong to a class fixed in advance, cheaper than $C$: linear readouts for reservoirs, a fixed change of variables for analogues.

So the questions asked of any substrate are: which class of $E$ and $R$, what $\varepsilon$ on what $K$, at what cost (§1.5), and how $\varepsilon$ was verified. A substrate alone is not a computer.

### 1.2 Four Modes of Physical Computation

1. **Analogue, by conjugacy.** The substrate obeys the target's law up to a change of variables $h$: $h \circ \Phi^t_{\mathcal{S}} = \Phi^t_{\text{target}} \circ h$ (a conjugacy; a semiconjugacy when $h$ forgets). Phillips's hydraulic machine of the economy (1950), in which water computes the circular flow of income, electrical analogues of mechanical systems, and excitable chemical media as models of cardiac tissue are of this kind. It is checked by the laws through $h$, and then holds for every input at once.
2. **Reservoir, by training.** A generic driven substrate $s_{t+1} = F(s_t, u_t)$ with a trained readout $y_t = W \psi(s_t)$. If $F$ has the *echo state property* (the state is asymptotically a function of the input's past alone, whatever the initial state) and the family of reservoirs separates points, linear readouts approximate uniformly any causal, time-invariant filter with fading memory (Boyd and Chua 1985; Maass, Natschläger and Markram 2002; Grigoryeva and Ortega 2018). Water in a bucket (Fernando and Sojakka 2003), a spin-torque oscillator (Torrejon et al. 2017), a delayed optoelectronic loop (Appeltant et al. 2011), memristors, mycelium and neural organoids (Cai et al. 2023) have been used so.
3. **Programmed, by design.** The substrate's dynamics are designed to carry out a given computation. Deterministic chemical reaction networks stably compute exactly the semilinear functions (Chen, Doty and Soloveichik 2014); stochastic ones simulate Turing machines with an error that can be made as small as wished (Soloveichik, Cook, Winfree and Bruck 2008); DNA strand displacement realises them (Qian and Winfree 2011). Ising machines relax towards a Hamiltonian's ground states (Inagaki et al. 2016), p-bits sample its Boltzmann law (Camsari, Faria, Sutton and Datta 2017), and memristor crossbars multiply a matrix by a vector by Ohm's and Kirchhoff's laws.
4. **Learning, by the substrate's own rules.** The substrate is trained by its physics rather than by a digital computer. A network of resistors, pipes or springs adjusts its own conductances or stiffnesses by a local rule that compares a free state with a gently clamped one (coupled learning, Stern, Hexner, Rocks and Liu 2021; equilibrium propagation, Scellier and Bengio 2017), as a self-adjusting resistor network has shown (Dillavou, Stern, Liu and Durian 2022); or it is trained by backpropagation through a digital twin of itself (Wright et al. 2022). The first way removes the digital training, often the largest energy cost.

*Remark.* The modes differ in what is verified: the law (1), the filter on held-out inputs (2), the program's specification (3), the learning rule's convergence and the filter learnt (4).

### 1.3 Simulating and Producing Need Opposite Properties

The user's first goal asks for a substrate that can *simulate* complex phenomena (chaos, nonlinearity, changing attractors, critical transitions) or *really produce* them. These ask for opposite properties.

- A reservoir must forget its initial state, so that it is repeatable: the echo state property is a uniform contraction of the input-driven state. It *computes* chaos only when its output is fed back as its input, which makes the closed loop an autonomous system whose attractor imitates the target's (Pathak, Hunt, Girvan, Lu and Ott 2018, on the Kuramoto–Sivashinsky equation, for several Lyapunov times).
- A substrate that *produces* chaos is sensitive to its initial state. It is repeatable only in its statistics, or when its input drives it into generalised synchronisation.

The boundary between the two regimes, the "edge of chaos", is where much of the literature looks for computational power (Bertschinger and Natschläger 2004). That it is optimal in general is disputed (Mitchell, Hraber and Crutchfield 1993); the optimum depends on the task. The domain makes the regime a declared, tunable parameter of every substrate, measured by its conditional Lyapunov exponent under its input.

### 1.4 What Transfers from a Substrate to its Target

Chaos makes trajectories incomparable after a few Lyapunov times, so a substrate is held to a target by invariants: the invariant measure (a Wasserstein or energy distance between long-run empirical laws), the Lyapunov spectrum and its Kaplan–Yorke dimension, power spectra, and the places and kinds of its bifurcations.

Knowledge that an AI finds on a substrate transfers to another system by universality. Near a codimension-one bifurcation every smooth system is locally equivalent to the bifurcation's normal form (saddle-node, transcritical, pitchfork, Hopf). So critical slowing down, the lag-one autocorrelation rising towards one and the variance growing as a fold is approached (Scheffer et al. 2009), is shared, and so are the scaling exponents; the rest is not. A substrate is a world for discovery exactly as far as its normal forms and universality classes are the target's.

### 1.5 The Cost of a Verified Result

The user's second criterion is measured per verified result, not per operation:
$$ \text{cost} = \frac{E_{\text{fabrication}} + E_{\text{materials}}}{N_{\text{results over the lifetime}}} + E_{\text{operation}} + E_{\text{conversion}} + E_{\text{verification}}, $$
in joules, each term with its source declared, and those not measured declared as assumptions.

Two floors are physical:

- erasing one bit at temperature $T$ dissipates at least $k_B T \ln 2$ (Landauer 1961), $2.87 \times 10^{-21}$ J at 300 K;
- a current of a Markov jump process in its steady state, measured to relative precision $\epsilon$ ($\epsilon^2 = \mathrm{Var}\,J / \langle J\rangle^2$), costs a dissipation of at least $2 k_B T / \epsilon^2$ (the thermodynamic uncertainty relation; Barato and Seifert 2015; Gingrich, Horowitz, Perunov and England 2016).

Both are computable exactly on physis's master equations, which already give entropy production (physis decision 0016). The term most often left out is the conversion between the analogue substrate and the digital world, an analogue-to-digital conversion per reading; Murmann's survey of converters gives its energy against resolution and rate, and it often outweighs the substrate's own.

Energy is not the whole of it. The materials are counted by five marks:

1. **Supply.** Whether the substrate needs a raw material on the EU's list of critical raw materials (Regulation (EU) 2024/1252, Annex II): vanadium, niobium, germanium and antimony for the phase-change oxides and chalcogenides; cobalt, tantalum, tungsten and platinum for most spintronics. *Caveat.* The list measures economic importance and the risk to supply, not scarcity or harm alone; silicon metal is on it too.
2. **Toxicity and persistence.** Liquid-crystal monomers are persistent, bioaccumulative and toxic (Su et al. 2019); acrylamide is toxic where alginate and agarose are not.
3. **Reuse.** Whether a run consumes the substrate: a foam coarsens, a DNA circuit burns its fuel, a reaction its reagents, an active colloid its peroxide; water in a tank, a circuit and a metamaterial are used again.
4. **Conditions.** Room temperature, or cryogenics and vacuum.
5. **End of life.** Biodegradable, recyclable, or neither.

### 1.6 Generativity: Breadth and Depth

How many phenomena a substrate can produce is two questions, not one.

- **Breadth** is the class of laws the substrate can be made to obey, up to conjugacy. For three classes it is greatest, and there it is a theorem.
  - The electronic analogue computer (integrators, adders and multipliers) is Shannon's general-purpose analogue computer, which generates exactly the solutions of polynomial ordinary differential equations (Shannon 1941; Graça and Costa 2003); polynomial time is characterised by polynomial ODEs of polynomial length (Bournez, Graça and Pouly 2017).
  - Continuous chemical reaction networks realise any polynomial ODE after a dual-rail encoding $x = x^+ - x^-$, with $x^+ + x^- \to \varnothing$, and are Turing complete over the reals (Fages, Le Guludec, Bournez and Pouly 2017).
  - DNA strand displacement realises any reaction network up to its rate constants (Soloveichik, Seelig and Winfree 2010).

  Every polynomial chaos (Lorenz, Rössler), every polynomial food web, every reaction–diffusion by the method of lines is then in reach, within the substrate's precision and dynamic range. Below these come networks of one kind of unit on a programmable graph (oscillators, p-bits, a swarm of robots), whose breadth is that of the network dynamics of that unit.
- **Depth** is whether the substrate *is* an instance of its law rather than an approximation to it. Water obeys the Navier–Stokes equations; a rotating, differentially heated annulus of water has the atmosphere's baroclinic instability, its waves, vacillation and geostrophic turbulence, in a regime diagram set by two dimensionless numbers (Fultz et al. 1959; Hide and Mason 1975). A deep substrate produces what every model of its law omits.

Simulating asks for breadth and producing for depth: this is §1.3 again, seen from the substrate.

*Remark.* A deep substrate stands for its target by similarity. The dimensionless groups of a law form the kernel of its dimension matrix, a sublattice of $\mathbb{Z}^n$ (Buckingham 1914); physis's `units` already holds dimensions as a free abelian group. Two systems with the same law in dimensionless form and equal groups are conjugate by a rescaling. A laboratory cannot match every group (no annulus has the atmosphere's Reynolds number), and the transfer then holds only where the unmatched groups are asymptotically irrelevant (incomplete similarity; Barenblatt 1996). That is why small experimental landscapes reproduce real ones as well as they do (Paola, Straub, Mohrig and Reinhardt 2009).

## 2. Placement

### 2.1 A Domain, and what goes to the Foundations

The substrates name objects in nature (a memristor, a spin-torque oscillator, a Belousov–Zhabotinsky medium, a rotating annulus), so they belong to a domain. The theory of computation by dynamics names none; it goes to the foundations, written once, where every domain can use it:

| Where | What | Status there |
|---|---|---|
| physis `dynamics` | discrete-time and delay dynamics; Lyapunov spectra (Benettin et al. 1980), conditional exponents of driven systems; invariant measures; bifurcations of equilibria and cycles, their normal forms and their continuation; test systems with known invariants (Lorenz 63 and 96, Kuramoto–Sivashinsky, Mackey–Glass, the logistic map) | planned in its README |
| physis `control` | systems with inputs and outputs; the echo state property and its sufficient conditions; system identification | planned |
| physis `units` | dimensionless groups as the kernel of the dimension matrix; similarity of a substrate and its target | dimensions written |
| physis `statistical_mechanics` | Landauer's bound and the thermodynamic uncertainty relation; spin models and their transitions, for Ising machines and p-bits | entropy production written; spin models planned |
| physis `chemistry` | the Oregonator and excitable media on its reaction–diffusion; polynomial ODEs compiled to reaction networks by dual rails | kinetics and reaction–diffusion written |
| physis `discretisation` | flows whose conductances adapt by a local rule (coupled learning, *Physarum*'s adaptation) | flows with conductances written |
| tyche `learning` | readouts; memory and information-processing capacity (Dambre, Verstraeten, Schrauwen and Massar 2012); next-generation reservoir computing (Gauthier, Bollt, Griffith and Barbosa 2021), the digital baseline every substrate must beat at equal cost | networks written |
| tyche `time_series` | delay embeddings (Takens 1981; false nearest neighbours), the correlation dimension, early warning signals of transitions (Dakos et al. 2012) | surrogates written |
| tyche `information` | transfer entropy, active information storage | planned |
| heuresis | configuring a substrate by search (evolution in materio; Miller and Downing 2002) | its executors and samplers serve as they are |

### 2.2 Boundaries with the Other Libraries

- **calx.** A device that executes calx's IR at a declared precision is a calx backend: the GPUs, neuromorphic chips and quantum processors already planned there. A substrate that computes by its own physics, with an error that is part of its physics, is this domain's. Whether calx will one day accept analogue backends with a declared precision (an analogue matrix–vector product) is calx's decision, under its open precision policy (CLAUDE.md §8.5).
- **Substrates another domain holds.** A substrate whose objects another domain already names stays that domain's model: living neural tissue is neuroscience's, a plankton community or a chemostat ecology's, a bacterial biofilm microbiology's, a swarm of robots robotics's. Used as a computer, it is composed with this domain's contract (§1.1) in a work, so that physical computing need not stand above ecology (height 5) to use one. Networks of spiking neurons as models of the brain stay in neuroscience.
- **robotics.** A compliant body as a reservoir (morphological computation; Hauser et al. 2011; Nakajima et al. 2013) is a substrate here, on physis's mechanics; robotics's bodies may use it.
- **heuresis.** A quantum device as a sampler of a law stays heuresis's seam. An Ising machine or an array of p-bits modelled from its physics is this domain's, and may serve as one of heuresis's samplers.

### 2.3 The Domain's Place in the Order

A field at height 0, on calx, morphe, physis and tyche (and heuresis when a module configures a substrate by search), C++ from the start, listed in the README's §4.4 under engineering beside electrical-engineering. By §2.2 it composes other domains' substrates in works rather than declaring them, so it rises only if a module of its own needs another domain's objects.

### 2.4 Real Substrates

The domain simulates substrates and holds their models; it touches no device. A real substrate (a dish of the BZ medium under a camera, a rotating tank, a memristor on a probe station) is an instrument behind `observation/`, its readings held to the domain's model as observed beside modelled, as Cosmopoiesis's instruments are.

## 3. The Substrates

### 3.1 The First Candidates

Assessed against the user's criteria and against what the infrastructure already holds. *Phenomena* is what the substrate is known to show; *beneath*, what it would stand on; *checks*, the closed forms or published results it would be held to.

| Substrate | Modes | Phenomena | Beneath | Checks | Maturity and cost |
|---|---|---|---|---|---|
| Coupled nonlinear oscillators: spin-torque, Mott (VO₂, NbO₂), electronic | reservoir, programmed | synchronisation transitions, multistability, chimeras, chaos | physis `dynamics` (Stuart–Landau, Kuramoto); the auto-oscillator model (Slavin and Tiberkevich 2009) | Kuramoto's threshold $K_c = 2\gamma$ and order $r = \sqrt{1 - 2\gamma/K}$ for Lorentzian frequencies (Kuramoto 1984; Ott and Antonsen 2008) | vowels recognised by four coupled oscillators (Romera et al. 2018); spin-torque and Mott devices need critical raw materials (§1.5), electronic ones do not |
| Memristors and crossbars | programmed, reservoir | hysteresis, negative differential resistance, chaos (Kumar, Strachan and Williams 2017) | a state-dependent Ohm's law; circuits | the pinched loop (Chua 1971, 2014); the HP model (Strukov et al. 2008) | nearest to practice; analogue precision of a few bits, drift; conversion dominates the cost |
| Ising machines, p-bits | programmed | phase transitions, sampling, critical slowing down | physis `statistical_mechanics` (spin models) | Onsager's critical coupling $\beta J_c = \tfrac12 \ln(1 + \sqrt 2)$ on the square lattice | p-bits from stochastic magnetic tunnel junctions at room temperature; coherent Ising machines optical |
| Chemical media and reaction networks (BZ, CRNs, DNA) | analogue, programmed, reservoir | oscillations, waves, spirals, Turing patterns, chaos (Roux, Simoyi and Swinney 1983) | physis `chemistry`, written | the Oregonator's Hopf point (Field and Noyes 1974); the semilinear functions of CRNs | slow (seconds to minutes), consumes reagents; best as a scientific world (Parrilla-Gutierrez et al. 2020) |
| Photonics: meshes, delay loops | programmed, reservoir | interference (linear), delay-induced chaos | physis `symmetry` (any unitary as a mesh; Reck et al. 1994; Clements et al. 2016); `radiation`, planned | unitary decompositions, exact | fast; linear operations cheap, nonlinearity costly (Shen et al. 2017; Larger et al. 2012) |
| Water and fluids | analogue, reservoir | waves, turbulence | physis `continuum` (fluids planned) | Phillips's machine against the model of income it was built for (economics) | the substrate is free; the reading is by camera (§3.2 for the rotating annulus) |
| Soft bodies | reservoir | nonlinear elasticity, contact | physis `mechanics`, written | Nakajima et al.'s arm | the body computes as it moves |
| Mycelium, *Physarum* | reservoir, programmed | electrical spiking (Adamatzky 2022), adaptive networks | physis `discretisation` (flows with conductances) | *Physarum*'s dynamics converge to the shortest path (Tero, Kobayashi and Nakagaki 2007; Bonifaci, Mehlhorn and Varma 2012) | grown and biodegradable; slow and variable |
| Neural cultures, organoids | reservoir | spiking, plasticity | neuroscience `snn`, composed in a work (§2.2) | — | living tissue; its ethics unsettled (Smirnova et al. 2023) |
| Quantum reservoirs | reservoir | unitary and open dynamics | physis `quantum`, written | small systems exactly (Fujii and Nakajima 2017) | mostly cryogenic |

*Remark.* *Physarum* fits physis closely: Tero, Kobayashi and Nakagaki's model is a Kirchhoff flow whose conductances grow with their flux and decay without it, and physis already holds flows with conductances on complexes (the economics domain uses them). Coupled learning (§1.2, mode 4) is the same structure with a different local rule.

### 3.2 Further Candidates, for Generativity and the Environment

The second list passed on, and those it leaves out, assessed by §1.5's marks and §1.6's two axes.

| Substrate | Breadth, depth | Phenomena | Materials and reuse | Beneath | Checks |
|---|---|---|---|---|---|
| Electronic analogue computers; field-programmable analogue arrays | broadest: polynomial ODEs (§1.6) | any polynomial dynamics, chaos among them, within a precision near $10^{-3}$ | silicon and copper; reusable, rewired by patching or switches; low power in subthreshold CMOS (Mead 1989; Hasler 2020) | physis `dynamics` | the compiled ODE's invariants; components' gain, offset, saturation and noise (Ulmann 2022); PDEs by hybrid analogue and digital solution (Huang et al. 2017) |
| Water in rotating, stratified and convecting layers; Faraday waves | deepest for the atmosphere and ocean | baroclinic waves, vacillation, geostrophic turbulence (Hide and Mason 1975); convection rolls and spiral defect chaos (Morris, Bodenschatz, Cannell and Ahlers 1993); walking droplets (Couder and Fort 2006; Bush 2015) | water, heat and a turntable; wholly reusable; read by camera | physis `continuum` (fluids planned); meteorology and oceanography as targets | the annulus's regime diagram; the onset of convection at $Ra_c = 1707.76$ between rigid plates (Chandrasekhar 1961) |
| Ecological microcosms: chemostats, plankton | deep for populations | a Hopf bifurcation crossed in a live predator–prey system (Fussmann, Ellner, Shertzer and Hairston 2000); chaos in a plankton community (Benincà et al. 2008) and a microbial food web (Becks et al. 2005); critical slowing down before collapse (Dai, Vorselen, Korolev and Gore 2012); early warnings in a whole lake (Carpenter et al. 2011) | grows itself; biodegradable; slow (days) and variable | ecology's models, composed in a work (§2.2); tyche's early warnings | the published bifurcations and indicators |
| Synthetic gene circuits; biofilms | broad within designed circuits; deep for bacterial excitable media | oscillators (Elowitz and Leibler 2000), switches (Gardner, Cantor and Collins 2000), synchronised clocks (Danino et al. 2010), patterns (Basu et al. 2005), electrical waves in biofilms (Prindle et al. 2015; Liu et al. 2015) | room temperature, fed on sugar; engineered organisms need containment | cell-biology and microbiology's models, composed in a work | the circuits' bifurcations |
| Enzymatic DNA networks (the PEN toolbox) | broad: any reaction network (§1.6) | oscillators, bistability, travelling waves, a French flag (Montagne et al. 2011; Zadorin et al. 2017) | biodegradable; its fuel consumed; synthesis costly | physis `chemistry` | the compiled network against its ODE |
| Mechanical metamaterials; mechanical oscillators | multistable elements on a designed graph | multistability, hysteresis, counting (Kwakernaak and van Hecke 2023; Bertoldi et al. 2017); chimeras among coupled metronomes (Martens, Thutupalli, Fourrière and Hallatschek 2013) | polymers, silicone, steel; reusable | physis `mechanics` | the counted sequences; Kuramoto's threshold |
| Physical learning networks | broad for the filters they learn | resistor, flow and elastic networks that learn by local rules (§1.2, mode 4) | resistors and transistors; reusable; no digital training | physis `discretisation` | the coupled-learning fixed point; the filter learnt |
| Robot swarms | network dynamics of a programmable unit | self-assembly of a thousand robots (Rubenstein, Cornejo and Nagpal 2014), Turing patterns grown by a swarm (Slavkov et al. 2018) | reusable; batteries | robotics's models, composed in a work | the patterns' wavelengths |
| Liquid crystals; active nematics | deep: nematohydrodynamics | electroconvection and defect turbulence (Kramer and Pesch 1995); active turbulence (Sanchez et al. 2012) | low voltage, reusable; the monomers persistent and toxic | physis `continuum` | the onset of electroconvection |
| Granular media; experimental landscapes | deep, irreproducible in detail | avalanches and jamming; criticality in rice piles (Frette et al. 1996) but rarely in sand (Nagel 1992); landscapes and strata (Paola et al. 2009) | sand, rice, water; reusable | physis `mechanics`; geography and geophysics as targets | the laws of avalanche sizes |
| Active colloids | deep: many-body motion | flocking of Quincke rollers (Bricard et al. 2013); living crystals (Palacci et al. 2013) | hydrogen peroxide as fuel, or hexadecane as solvent | physis `statistical_mechanics` | the flocking transition |
| Foams; capillary networks | deep | coarsening by von Neumann's law $\dot A \propto (n - 6)$ (von Neumann 1952), topological rearrangements | soap, water, air; a run coarsens for good, but is cheaply remade | physis `computational_geometry` | von Neumann's law |
| Ionic membranes; hydrogels | moderate | excitability, memory, self-oscillating gels (Yoshida, Takahashi, Yamaguchi and Ichijo 1996) | alginate and agarose biodegradable; polyacrylamide's monomer toxic | physis `chemistry` | the gel's oscillation |
| Electrochemical oscillators | network of one unit | Kuramoto's transition seen on 64 nickel electrodes (Kiss, Zhai and Hudson 2002) | nickel in acid; electrodes dissolve slowly | physis `chemistry`, `dynamics` | Kuramoto's threshold |
| Phase-transition materials (VO₂, NbO₂, GeSbTe) | network of one unit | threshold switching, hysteresis, metastable states | vanadium, niobium, germanium and antimony, all critical raw materials | circuits | the measured switching curves |

*Caveat.* Generativity claimed is not generativity shown. The walking droplets' analogue of single-particle interference, for one, was not reproduced in a later, more careful experiment (Andersen et al. 2015). Every row's phenomena are to be held to reproducible checks before the domain relies on them.

### 3.3 The Earth as its Own Substrate

The greenest substrate is one that already runs. The atmosphere, the ocean, lakes and forests compute their own futures at no added cost, and the data library already reads them; for knowledge of the Earth and its ecosystems they are the substrate of last appeal. A built substrate earns its place by what nature gives only once: control, repetition and speed. Each candidate is therefore asked which class of phenomena it lets us repeat under control, at what cost per verified result, and what it shows that the observed Earth cannot.

## 4. A First Slice

Proposed so that the domain begins where it is broadest, most rooted and most checkable, and serves all three criteria.

1. **`computer`**: the contract of §1.1 (substrate, control, encoding, readout, model, tolerance) and the account of §1.5, energy and materials, with physis's physical floors; every substrate a model behind one interface.
2. **`analogue`**: Shannon's general-purpose analogue computer. A polynomial ODE is compiled to integrators, adders and multipliers, with the non-idealities of real components (finite gain, offset, saturation, noise), and each phenomenon is asked what precision it tolerates (chaos in its statistics, a fold in its place). Networks of phase and amplitude oscillators on a programmable graph belong here, with spin-torque oscillators by the auto-oscillator model as a realisation to compare, not preferred (§1.5); they serve both as reservoirs (spoken digits from the Free Spoken Digit Dataset, CC BY-SA 4.0) and as worlds (synchronisation, multistability, chaos).
3. **`chemical`**: the same polynomial ODEs compiled to reaction networks by dual rails and to DNA strand displacement, run by mass action and as exact Markov chains on physis's `chemistry`; and the BZ medium by the Oregonator.

Breadth is then checkable: Lorenz's attractor is realised in a circuit, a reaction network and a DNA network, each held to the ODE's invariants (its Lyapunov spectrum, its invariant measure).

Beneath them, each written first in its own library: physis's planned `dynamics` (Lyapunov spectra, bifurcations, delay, the test systems) and `control` (driven systems), and tyche's readouts, capacities and early warnings. Next come `fluids` (the rotating annulus and convection, once physis's fluids are written: depth for the atmosphere and the ocean), memristor crossbars, p-bits and physical learning networks; and works that compose ecology's microcosms with the contract.

A case study for the first goal, in `research/experiments/`: a substrate tuned towards a fold or a Hopf point, and an AI (heuresis's programs typed by dimension, as with Kepler's third law, or a language model as an external proposer) asked to find its law and to warn of its transition from its readings alone, scored against the known normal form.

## 5. Remarks on the Assessments Passed On

### 5.1 The First Assessment

The first assessment is a fair survey, and its distinction between simulating and producing is right. This proposal departs from it in four places.

1. Its condition ("encode, control, read") is necessary but not sufficient: without a restriction on the class and cost of $E$ and $R$ the claim is vacuous (§1.1).
2. Simulating and producing need opposite properties, not only different uses (§1.3); what transfers from a substrate is what is universal (§1.4).
3. It leaves out Ising machines and p-bits (the strongest for critical transitions and sampling), chemical reaction networks and DNA (the substrates with an exact theory of what they compute), soft bodies, *Physarum*, and the digital baseline, next-generation reservoir computing.
4. More than half its links are news summaries or aggregators rather than the papers, and several figures it reports (Extropic's factor of 10⁴, the memristor twin's gains) are, as it says itself, projections. Before the domain relies on any of them, a survey goes to the primary papers (CLAUDE.md §4.3), as was done for the reconstruction of vessels.

### 5.2 The Second List

The second list (mechanical metamaterials, ionic membranes and hydrogels, magnetic and spintronic networks, electrochemical oscillators, active colloids, phase-transition materials, granular media, capillary networks and foams) names media for their dynamics, and adds good candidates (§3.2). But interesting dynamics is not yet computing. It states for none of them an encoding, a readout or a tolerance. Three lose or consume their state in a run (foams coarsen, active colloids burn fuel, granular piles keep their history), which makes repetition costly. Two depend on critical raw materials (the phase-transition oxides; most spintronics). Its magnetic networks repeat the first assessment's spintronics. It also leaves out the three broadest substrates (analogue electronics, reaction networks, DNA), the deepest for the Earth (rotating and convecting water, ecological microcosms), and learning by the substrate's own rules.

## 6. Open Questions

1. **The domain's name.** `physical-computing` (proposed; Horsman et al. 2014 and Nakajima's "physical reservoir computing" use it), or `unconventional-computing` (Adamatzky's term, and Jaeger, Noheda and van der Wiel's 2023), which is defined against the digital computer and so ages.
2. **Living substrates.** Whether living neural tissue and engineered organisms are in scope, and under what ethical rule and rule of containment.
3. **A work to draw it.** A workbench of substrates; Cosmopoiesis, which already constructs worlds; or case studies and the gallery at first.
4. **Analogue backends.** Whether calx accepts them with a declared precision; calx's decision.

## References

- Adamatzky, A. (2022), Language of fungi derived from their electrical spiking activity, *Royal Society Open Science* 9, 211926.
- Andersen, A. et al. (2015), Double-slit experiment with single wave-driven particles and its relation to quantum mechanics, *Physical Review E* 92, 013006.
- Appeltant, L. et al. (2011), Information processing using a single dynamical node as complex system, *Nature Communications* 2, 468.
- Barato, A. C. and Seifert, U. (2015), Thermodynamic uncertainty relation for biomolecular processes, *Physical Review Letters* 114, 158101.
- Barenblatt, G. I. (1996), *Scaling, Self-Similarity, and Intermediate Asymptotics*, Cambridge University Press.
- Basu, S., Gerchman, Y., Collins, C. H., Arnold, F. H. and Weiss, R. (2005), A synthetic multicellular system for programmed pattern formation, *Nature* 434, 1130–1134.
- Becks, L., Hilker, F. M., Malchow, H., Jürgens, K. and Arndt, H. (2005), Experimental demonstration of chaos in a microbial food web, *Nature* 435, 1226–1229.
- Benettin, G., Galgani, L., Giorgilli, A. and Strelcyn, J.-M. (1980), Lyapunov characteristic exponents for smooth dynamical systems and for Hamiltonian systems; a method for computing all of them, *Meccanica* 15, 9–20 and 21–30.
- Benincà, E. et al. (2008), Chaos in a long-term experiment with a plankton community, *Nature* 451, 822–825.
- Bertoldi, K., Vitelli, V., Christensen, J. and van Hecke, M. (2017), Flexible mechanical metamaterials, *Nature Reviews Materials* 2, 17066.
- Bertschinger, N. and Natschläger, T. (2004), Real-time computation at the edge of chaos in recurrent neural networks, *Neural Computation* 16, 1413–1436.
- Bonifaci, V., Mehlhorn, K. and Varma, G. (2012), Physarum can compute shortest paths, *Journal of Theoretical Biology* 309, 121–133.
- Bournez, O., Graça, D. S. and Pouly, A. (2017), Polynomial time corresponds to solutions of polynomial ordinary differential equations of polynomial length, *Journal of the ACM* 64, 38.
- Boyd, S. and Chua, L. O. (1985), Fading memory and the problem of approximating nonlinear operators with Volterra series, *IEEE Transactions on Circuits and Systems* 32, 1150–1161.
- Bricard, A., Caussin, J.-B., Desreumaux, N., Dauchot, O. and Bartolo, D. (2013), Emergence of macroscopic directed motion in populations of motile colloids, *Nature* 503, 95–98.
- Buckingham, E. (1914), On physically similar systems; illustrations of the use of dimensional equations, *Physical Review* 4, 345–376.
- Bush, J. W. M. (2015), Pilot-wave hydrodynamics, *Annual Review of Fluid Mechanics* 47, 269–292.
- Cai, H. et al. (2023), Brain organoid reservoir computing for artificial intelligence, *Nature Electronics* 6, 1032–1039.
- Camsari, K. Y., Faria, R., Sutton, B. M. and Datta, S. (2017), Stochastic p-bits for invertible logic, *Physical Review X* 7, 031014.
- Carpenter, S. R. et al. (2011), Early warnings of regime shifts: a whole-ecosystem experiment, *Science* 332, 1079–1082.
- Chandrasekhar, S. (1961), *Hydrodynamic and Hydromagnetic Stability*, Oxford University Press.
- Chen, H.-L., Doty, D. and Soloveichik, D. (2014), Deterministic function computation with chemical reaction networks, *Natural Computing* 13, 517–534.
- Chua, L. O. (1971), Memristor — the missing circuit element, *IEEE Transactions on Circuit Theory* 18, 507–519.
- Chua, L. O. (2014), If it's pinched it's a memristor, *Semiconductor Science and Technology* 29, 104001.
- Clements, W. R. et al. (2016), Optimal design for universal multiport interferometers, *Optica* 3, 1460–1465.
- Couder, Y. and Fort, E. (2006), Single-particle diffraction and interference at a macroscopic scale, *Physical Review Letters* 97, 154101.
- Dai, L., Vorselen, D., Korolev, K. S. and Gore, J. (2012), Generic indicators for loss of resilience before a tipping point leading to population collapse, *Science* 336, 1175–1177.
- Dakos, V. et al. (2012), Methods for detecting early warnings of critical transitions in time series illustrated using simulated ecological data, *PLoS ONE* 7, e41010.
- Dambre, J., Verstraeten, D., Schrauwen, B. and Massar, S. (2012), Information processing capacity of dynamical systems, *Scientific Reports* 2, 514.
- Danino, T., Mondragón-Palomino, O., Tsimring, L. and Hasty, J. (2010), A synchronized quorum of genetic clocks, *Nature* 463, 326–330.
- Dillavou, S., Stern, M., Liu, A. J. and Durian, D. J. (2022), Demonstration of decentralized physics-driven learning, *Physical Review Applied* 18, 014040.
- Elowitz, M. B. and Leibler, S. (2000), A synthetic oscillatory network of transcriptional regulators, *Nature* 403, 335–338.
- European Union (2024), Regulation (EU) 2024/1252 establishing a framework for ensuring a secure and sustainable supply of critical raw materials (the Critical Raw Materials Act), Annex II.
- Fages, F., Le Guludec, G., Bournez, O. and Pouly, A. (2017), Strong Turing completeness of continuous chemical reaction networks and compilation of mixed analog-digital programs, *Computational Methods in Systems Biology* (CMSB 2017), LNCS 10545, 108–127.
- Fernando, C. and Sojakka, S. (2003), Pattern recognition in a bucket, *Advances in Artificial Life* (ECAL 2003), LNCS 2801, 588–597.
- Field, R. J. and Noyes, R. M. (1974), Oscillations in chemical systems. IV. Limit cycle behavior in a model of a real chemical reaction, *Journal of Chemical Physics* 60, 1877–1884.
- Frette, V. et al. (1996), Avalanche dynamics in a pile of rice, *Nature* 379, 49–52.
- Fujii, K. and Nakajima, K. (2017), Harnessing disordered-ensemble quantum dynamics for machine learning, *Physical Review Applied* 8, 024030.
- Fultz, D., Long, R. R., Owens, G. V., Bohan, W., Kaylor, R. and Weil, J. (1959), Studies of thermal convection in a rotating cylinder with some implications for large-scale atmospheric motions, *Meteorological Monographs* 4 (21), American Meteorological Society.
- Fussmann, G. F., Ellner, S. P., Shertzer, K. W. and Hairston, N. G. (2000), Crossing the Hopf bifurcation in a live predator–prey system, *Science* 290, 1358–1360.
- Gardner, T. S., Cantor, C. R. and Collins, J. J. (2000), Construction of a genetic toggle switch in *Escherichia coli*, *Nature* 403, 339–342.
- Gauthier, D. J., Bollt, E., Griffith, A. and Barbosa, W. A. S. (2021), Next generation reservoir computing, *Nature Communications* 12, 5564.
- Gingrich, T. R., Horowitz, J. M., Perunov, N. and England, J. L. (2016), Dissipation bounds all steady-state current fluctuations, *Physical Review Letters* 116, 120601.
- Graça, D. S. and Costa, J. F. (2003), Analog computers and recursive functions over the reals, *Journal of Complexity* 19, 644–664.
- Grigoryeva, L. and Ortega, J.-P. (2018), Echo state networks are universal, *Neural Networks* 108, 495–508.
- Hasler, J. (2020), Large-scale field-programmable analog arrays, *Proceedings of the IEEE* 108, 1283–1302.
- Hauser, H., Ijspeert, A. J., Füchslin, R. M., Pfeifer, R. and Maass, W. (2011), Towards a theoretical foundation for morphological computation with compliant bodies, *Biological Cybernetics* 105, 355–370.
- Hide, R. and Mason, P. J. (1975), Sloping convection in a rotating fluid, *Advances in Physics* 24, 47–100.
- Horsman, C., Stepney, S., Wagner, R. C. and Kendon, V. (2014), When does a physical system compute?, *Proceedings of the Royal Society A* 470, 20140182.
- Huang, Y., Guo, N., Seok, M., Tsividis, Y., Mandli, K. and Sethumadhavan, S. (2017), Hybrid analog-digital solution of nonlinear partial differential equations, *Proceedings of the 50th IEEE/ACM International Symposium on Microarchitecture* (MICRO-50).
- Inagaki, T. et al. (2016), A coherent Ising machine for 2000-node optimization problems, *Science* 354, 603–606.
- Jaeger, H., Noheda, B. and van der Wiel, W. G. (2023), Toward a formal theory for computing machines made out of whatever physics offers, *Nature Communications* 14, 4911.
- Kiss, I. Z., Zhai, Y. and Hudson, J. L. (2002), Emerging coherence in a population of chemical oscillators, *Science* 296, 1676–1678.
- Kramer, L. and Pesch, W. (1995), Convection instabilities in nematic liquid crystals, *Annual Review of Fluid Mechanics* 27, 515–541.
- Kumar, S., Strachan, J. P. and Williams, R. S. (2017), Chaotic dynamics in nanoscale NbO₂ Mott memristors for analogue computing, *Nature* 548, 318–321.
- Kuramoto, Y. (1984), *Chemical Oscillations, Waves, and Turbulence*, Springer.
- Kwakernaak, L. J. and van Hecke, M. (2023), Counting and sequential information processing in mechanical metamaterials, *Physical Review Letters* 130, 268204.
- Landauer, R. (1961), Irreversibility and heat generation in the computing process, *IBM Journal of Research and Development* 5, 183–191.
- Larger, L. et al. (2012), Photonic information processing beyond Turing: an optoelectronic implementation of reservoir computing, *Optics Express* 20, 3241–3249.
- Liu, J. et al. (2015), Metabolic co-dependence gives rise to collective oscillations within biofilms, *Nature* 523, 550–554.
- Maass, W., Natschläger, T. and Markram, H. (2002), Real-time computing without stable states: a new framework for neural computation based on perturbations, *Neural Computation* 14, 2531–2560.
- Martens, E. A., Thutupalli, S., Fourrière, A. and Hallatschek, O. (2013), Chimera states in mechanical oscillator networks, *Proceedings of the National Academy of Sciences* 110, 10563–10567.
- Mead, C. (1989), *Analog VLSI and Neural Systems*, Addison-Wesley.
- Miller, J. F. and Downing, K. (2002), Evolution in materio: looking beyond the silicon box, *Proceedings of the NASA/DoD Conference on Evolvable Hardware*, 167–176.
- Mitchell, M., Hraber, P. T. and Crutchfield, J. P. (1993), Revisiting the edge of chaos: evolving cellular automata to perform computations, *Complex Systems* 7, 89–130.
- Montagne, K., Plasson, R., Sakai, Y., Fujii, T. and Rondelez, Y. (2011), Programming an in vitro DNA oscillator using a molecular networking strategy, *Molecular Systems Biology* 7, 466.
- Morris, S. W., Bodenschatz, E., Cannell, D. S. and Ahlers, G. (1993), Spiral defect chaos in large aspect ratio Rayleigh–Bénard convection, *Physical Review Letters* 71, 2026–2029.
- Murmann, B., *ADC Performance Survey 1997–present*, online, updated yearly.
- Nagel, S. R. (1992), Instabilities in a sandpile, *Reviews of Modern Physics* 64, 321–325.
- Nakajima, K. et al. (2013), A soft body as a reservoir: case studies in a dynamic model of octopus-inspired soft robotic arm, *Frontiers in Computational Neuroscience* 7, 91.
- Ott, E. and Antonsen, T. M. (2008), Low dimensional behavior of large systems of globally coupled oscillators, *Chaos* 18, 037113.
- Palacci, J., Sacanna, S., Steinberg, A. P., Pine, D. J. and Chaikin, P. M. (2013), Living crystals of light-activated colloidal surfers, *Science* 339, 936–940.
- Paola, C., Straub, K., Mohrig, D. and Reinhardt, L. (2009), The "unreasonable effectiveness" of stratigraphic and geomorphic experiments, *Earth-Science Reviews* 97, 1–43.
- Parrilla-Gutierrez, J. M. et al. (2020), A programmable chemical computer with memory and pattern recognition, *Nature Communications* 11, 1442.
- Pathak, J., Hunt, B., Girvan, M., Lu, Z. and Ott, E. (2018), Model-free prediction of large spatiotemporally chaotic systems from data: a reservoir computing approach, *Physical Review Letters* 120, 024102.
- Phillips, A. W. (1950), Mechanical models in economic dynamics, *Economica* 17, 283–305.
- Prindle, A. et al. (2015), Ion channels enable electrical communication in bacterial communities, *Nature* 527, 59–63.
- Putnam, H. (1988), *Representation and Reality*, MIT Press, appendix.
- Qian, L. and Winfree, E. (2011), Scaling up digital circuit computation with DNA strand displacement cascades, *Science* 332, 1196–1201.
- Reck, M., Zeilinger, A., Bernstein, H. J. and Bertani, P. (1994), Experimental realization of any discrete unitary operator, *Physical Review Letters* 73, 58–61.
- Romera, M. et al. (2018), Vowel recognition with four coupled spin-torque nano-oscillators, *Nature* 563, 230–234.
- Roux, J.-C., Simoyi, R. H. and Swinney, H. L. (1983), Observation of a strange attractor, *Physica D* 8, 257–266.
- Rubenstein, M., Cornejo, A. and Nagpal, R. (2014), Programmable self-assembly in a thousand-robot swarm, *Science* 345, 795–799.
- Sanchez, T., Chen, D. T. N., DeCamp, S. J., Heymann, M. and Dogic, Z. (2012), Spontaneous motion in hierarchically assembled active matter, *Nature* 491, 431–434.
- Scellier, B. and Bengio, Y. (2017), Equilibrium propagation: bridging the gap between energy-based models and backpropagation, *Frontiers in Computational Neuroscience* 11, 24.
- Scheffer, M. et al. (2009), Early-warning signals for critical transitions, *Nature* 461, 53–59.
- Shannon, C. E. (1941), Mathematical theory of the differential analyzer, *Journal of Mathematics and Physics* 20, 337–354.
- Shen, Y. et al. (2017), Deep learning with coherent nanophotonic circuits, *Nature Photonics* 11, 441–446.
- Slavin, A. and Tiberkevich, V. (2009), Nonlinear auto-oscillator theory of microwave generation by spin-polarized current, *IEEE Transactions on Magnetics* 45, 1875–1918.
- Slavkov, I. et al. (2018), Morphogenesis in robot swarms, *Science Robotics* 3, eaau9178.
- Smirnova, L. et al. (2023), Organoid intelligence (OI): the new frontier in biocomputing and intelligence-in-a-dish, *Frontiers in Science* 1, 1017235.
- Soloveichik, D., Cook, M., Winfree, E. and Bruck, J. (2008), Computation with finite stochastic chemical reaction networks, *Natural Computing* 7, 615–633.
- Soloveichik, D., Seelig, G. and Winfree, E. (2010), DNA as a universal substrate for chemical kinetics, *Proceedings of the National Academy of Sciences* 107, 5393–5398.
- Stern, M., Hexner, D., Rocks, J. W. and Liu, A. J. (2021), Supervised learning in physical networks: from machine learning to learning machines, *Physical Review X* 11, 021045.
- Strukov, D. B., Snider, G. S., Stewart, D. R. and Williams, R. S. (2008), The missing memristor found, *Nature* 453, 80–83.
- Su, H. et al. (2019), Persistent, bioaccumulative, and toxic properties of liquid crystal monomers and their detection in indoor residential dust, *Proceedings of the National Academy of Sciences* 116, 26450–26458.
- Takens, F. (1981), Detecting strange attractors in turbulence, *Dynamical Systems and Turbulence*, Lecture Notes in Mathematics 898, 366–381.
- Tero, A., Kobayashi, R. and Nakagaki, T. (2007), A mathematical model for adaptive transport network in path finding by true slime mold, *Journal of Theoretical Biology* 244, 553–564.
- Torrejon, J. et al. (2017), Neuromorphic computing with nanoscale spintronic oscillators, *Nature* 547, 428–431.
- Ulmann, B. (2022), *Analog Computing*, second edition, De Gruyter.
- von Neumann, J. (1952), Discussion, in *Metal Interfaces*, American Society for Metals, 108–110.
- Wright, L. G. et al. (2022), Deep physical neural networks trained with backpropagation, *Nature* 601, 549–555.
- Yoshida, R., Takahashi, T., Yamaguchi, T. and Ichijo, H. (1996), Self-oscillating gel, *Journal of the American Chemical Society* 118, 5134–5135.
- Zadorin, A. S. et al. (2017), Synthesis and materialization of a reaction–diffusion French flag pattern, *Nature Chemistry* 9, 990–996.

*Caveat.* These references were written from memory for this seed; they are to be checked against the papers before any of them is cited in a module or a decision.
