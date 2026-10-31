# 0060. Streets as Open Systems: Street Oceanography Made Exact, the Urban-Science Domain and the Urban Building Workbench

- **Status:** Accepted, provisionally. The domain's name and the workbench's colour are provisional; the workbench's other plug-ins are planned.
- **Date:** 2026-10-31
- **Decided by:**
  - **The request.** The user asked Claude to "build the urban building workbench. And street oceanography lab as one part of it", from their papers on street oceanography. The workbench was to "utilize dynamical system formalism, especially fluid dynamics, chemical dynamics, etc to analyze street". It was also to equip general street design ("human flow analysis, architecture design, cultural analysis"), beginning "from what we have": "As street oceanography already a formalism, we begin from it."
  - **The physics and mathematics asked for.** "From this opportunity we can also implement relative fluid dynamics, chemical dynamics, thermodynamics, tide, etc's physical library. And frequency domain, wave, spectrum analysis libraries." The user's examples: a street holds activities of different frequencies (some stores close earlier, some later); kinds of people flow to particular stores at particular times; a store coupled with another "may like chemical structure, exhibit special functions … Like marriage store and diamond stores coupled at space time together".
  - **The user's correction, during the work.** "Not necessary mention the paper details in it. But implement a 3D GIS System.. And Visualization, analysis, etc can on top of it." The workbench was rebuilt on a three-dimensional GIS of cities (§2.5). Its analyses speak plainly of what was measured, proved or modelled, and how. The papers' apparatus stays in the domain's theory note.
  - **What Claude did.** It designed the work and wrote this record. CLAUDE.md §1 retires "laboratory", so the user's "street oceanography lab" is the workbench's plug-in *Street Oceanography*.

## 1. Context

The user's working papers, *Street Oceanography: An Oceanographic Epistemology of Socio-Cultural Dynamical Systems* ([SO]) and *Towards Street Oceanography: Ginza Street as a Case Study* ([G]), propose a field organised around the street as an open system. Heterogeneous flows pass through a heterogeneous terrain, forced on many timescales; each alters what passes through it and is altered in turn. The paper borrows from oceanography a disciplinary form, not its equations, and it withholds dynamical laws (Withholdings 2.1–2.5). It names:
- nine domains of question and eight classes of instrument;
- the standing of each borrowed feature: carried, carried in part, held in reserve, or set aside;
- five null accounts a finding must defeat;
- the classification of segments as the field's first test;
- an order of work that begins with declared opening hours.

The infrastructure had no urban science. Its foundations had:
- the FFT and Welch's method (calx);
- regression and ARIMA (tyche);
- reaction kinetics, the Hodge decomposition of flows and heat conduction as a balance law (physis).

They had no multitaper, coherence, bispectrum, harmonic (tidal) analysis or Lomb–Scargle. They had no kinematic waves, no transport, and no thermodynamics of Markov processes.

## 2. Decision

### 2.1 The Standing of the Formalism

The domain keeps the paper's standing ([theory/street-oceanography.md](../../domains/urban-science/theory/street-oceanography.md)):
- what the paper carries is computed as a definition;
- what holds of any street is proved as a theorem;
- every mechanism (fluid, chemical, thermodynamic) is a representation in the sense of [SO] §24, offered for sufficiency, consequence or exclusion, never prediction.

Every finding the workbench shows names its standing, its instrument and blind spots, the null it competes with, and the uses its result makes available ([SO] Requirements 28.1 and 33.2).

### 2.2 The General Methods Go to the Foundations

- **calx, decision 0013:**
  - discrete prolate spheroidal tapers;
  - the tapered periodogram;
  - Thomson's multitaper estimate with adaptive weights;
  - cross spectra and coherence;
  - the bispectrum and bicoherence.
- **tyche, decision 0005:**
  - spectral inference: χ² intervals, the levels of coherence and bicoherence, the red-noise (AR(1)) null, surrogates, aliasing;
  - harmonic analysis at given frequencies, the method of tidal analysis, with the Rayleigh criterion and errors from the residual spectrum. The tides domain will put its astronomical constituents on it;
  - the Lomb–Scargle periodogram with Baluev's false-alarm probability;
  - in `spatial`, the co-location quotient with stratified permutation nulls, and Clark and Evans's test.
- **physis, decision 0016:**
  - in `continuum`, kinematic waves on networks of links (Lighthill–Whitham–Richards, Godunov's scheme in demand–supply form, junctions with or without FIFO, scheduled changes of a link's law) and the transport of a scalar on a graph;
  - a new module, `statistical_mechanics`: master equations, stationary distributions, entropy production, cycle affinities and Kolmogorov's criterion, the relaxation of relative entropy, phase-type first passage and absorption, and the entropy of mixing with its decomposition across scales;
  - in `chemistry`, reaction–diffusion on graphs, Turing's condition in closed form and on any Laplacian, and the pinning of fronts on a discrete chain (McKean's caricature in closed form).

The user's "thermodynamics" is here the statistical mechanics of Markov processes and mixing. It defines entropy production and free-energy-like departures from rates and compositions alone. The paper sets aside the energetic vocabulary for streets as figurative (Claim 27.2), and this needs no conservation law. Their "tide" is harmonic analysis at the calendars' frequencies. Their "wave" is the kinematic wave of a crowd and the spectral analysis of records.

### 2.3 The Urban-Science Domain

`domains/urban-science/` is a field at height 0, standing on calx, morphe, physis and tyche, in C++ from the start. Its namespace is `urban::`. Its modules:
- **`street`:** holdings, terrain and its dominance order, the terrain length and compression ratio, the departure from mixing, the effective terrain in time, and the visits of kinds of people.
- **`rhythm`:** the calendars as a lattice, records, spectra with classified lines, harmonic constants, selectivity, the civil share of the clock, bicoherence, and the paper's Figure 1.
- **`response`:** susceptibility, and the comparison across a repeating alteration.
- **`affinity`:** co-location, catalysis in space and time, and spacing.
- **`flow`:** kinematic streets, person-time, transit, and circulation.
- **`classification`:** the plane of segments.

Architecture, culture and movement beyond counts are planned modules.

The name *urban science* is provisional. README §4.5 keeps urban geography within geography: this domain holds streets, buildings and their flows as objects; geography keeps places on the Earth. The two will meet when the domain declares geography.

### 2.4 The Formalism Made Exact, and Twelve Amendments

The theory note proves or cites every result. Twelve amendments to the papers ([theory note §10](../../domains/urban-science/theory/street-oceanography.md)):
1. Claim 9.1 is a theorem, and stronger: no weighting is invariant under changes of unit. The dominance order is the weighting-free content, and with signs of effect unfixed nothing is comparable.
2. L_grad is defined as an integral scale.
3. The departure from mixing is constructed: the mutual information of place and use, decomposed exactly across scales.
4. Exchange blind to kind never raises the departure, which gives the persistence criterion its form. Fronts pin on a discrete street.
5. The day and the week are one calendar.
6. The record-length constraint is the Rayleigh criterion for separating lines.
7. Discriminator 29.1 needs separability.
8. Discriminator 29.1 needs an order of explanation, with its chance rate.
9. Bicoherence discriminates only the irregular part.
10. Selectivity of response is defined.
11. Person-time's exact budget is Little's identity.
12. Circulation is topological: only the non-exact part of the log rate ratios drives it. This joins *Street Spin Foam*'s prediction to the Hodge decomposition.

### 2.5 The Urban Building Workbench: a 3D GIS with Analyses on Top

`applications/urban-building-workbench/` serves on port 8970. Its base is a three-dimensional GIS of a city. The general part is the visualization layer's `gis`:
- the local frame of a city, exactly on WGS 84;
- rings, ribbons and colours;
- a three.js scene of extruded buildings, streets, spaces, columns and bars, with picking, map controls and the Sun's light and shadows.

On that base, the workbench draws the data library's city model: OpenStreetMap's buildings with their heights (a default where none is mapped, drawn paler), ways with their widths, parks, water and squares, for Melbourne's centre. A clock runs in the city's own time zone and lights it by the real Sun of the instant (the astronomy domain's subsolar point). The buildings stand as built, lowered to a fifth of their height, or of glass, so that what moves in the streets shows between towers. Ginza was drawn too at first; the user had it removed ("Remove Ginza."), and the city model, the holdings and the workbench now hold Melbourne alone.

Plug-ins put layers and analyses on top, as the Chemistry Workbench's plug-ins put views (ADR 0051):
- a plug-in is a directory with a manifest of its *layers* and *tools*, and the shell names none;
- a layer draws on the map, follows the clock, and answers for what is picked there with actions that open tools;
- a tool opens in a drawer beside the map and can act on it: select, fly, take two picked stores;
- a manifest of status *planned* names what is to come: *Human flow*, *Architecture* and *Cultural analysis* are planned so;
- the engines are built by `engine.py` into `data/cache/engines/urban/<plug-in>/` under the sessions' WSL lock, and kept alive by the server, which tells each where the data library's products are; the engine reads nothing else.

**Street Oceanography** draws the street as a sea on the map (§2.7), and its records:
- the tide on the streets, the currents of walkers each way, and the waves of the crowd on Swanston Street (a model);
- pedestrians at Melbourne's sensors, a column for the mean count at the clock's hour of the week;
- shops and services along the named streets, coloured by use, standing tall when open by their declared hours;
- businesses by use from Melbourne's census of land use.

It also has eight tools:
- **Tides and currents of the city:** every street's tide and current together (§2.7), and the surge on the clock's date.
- **Rhythm of a street:** a sensor's spectrum against red noise, its lines read on the calendars, harmonic constants, selectivity, whose clock it keeps at the changes of daylight saving, and bicoherence.
- **Weather and holidays:** susceptibility to rain, warmth and holidays, and every street under one weather.
- **Mix of uses:** the departure from an even mix and its scales, the distance over which uses change, the census years, and why no single score ranks streets.
- **Opening hours:** where and when a street is open, the rhythm of its hours, the order of opening and closing, and which kinds of people its open shops draw.
- **Stores that work together:** two stores picked on the map as a pathway at their real separation and hours, its visitors spreading along the street on the map, and which uses keep company beyond shared blocks.
- **Crowds and circulation:** Swanston Street ordinary and with its carriageway opened, and circulation in a district of blocks.
- **Kinds of streets:** Melbourne's streets in the plane of compression and selectivity, with Clark and Evans's test; a point picked flies the map there.

Its colour is the user's: "UI main Color please utilize gentle sea-like blue". It is a gentle sea blue (#5B8FB9), with sand (#C98B6B), the shore beside it, as its minor colour, on misty paper by day and a deep-sea night. Type is serif and lines are hairlines; data colours come from the visualization layer's validated maps, which gained a sequential *sea* and a diverging *sea-brick*.

### 2.7 The Street as a Sea, on the Map

The user found the first build's dynamics only in the drawer, as charts ("I cannot find dynamical system analysis, especially oceanography one, and visualization on 3d GIS"). They are now drawn on the map, and measured wherever the records allow. Melbourne's live record counts each way at 102 sensors (October 2024 to October 2026, the two directions summing to the total at every hour). Its sensors are read as an oceanographer reads tide gauges and current meters:
- **The tide** (the domain's `rhythm/tide.hpp`). Each street's counts are fitted at the calendars' constituents in civil time. Its block of street is drawn as a channel of water standing to the predicted level, a metre for every hundred people an hour, rising and falling as the clock runs. On a date in the record the water stands at that hour's count, coloured by its surge, with a tide mark at the prediction. Two more readings colour it: the hour of high water, and the kind of tide by Courtier's form number, with the day index of the lattice as the tidal species.
- **The currents** (`flow/currents.hpp`). Each counting line carries two lanes of walkers. Each lane holds as many walkers as walk there on average at that hour (N = λL/v), moving in real time, those walking forward on the left. They are coloured by whether they walk toward the centre of the counts. An arrow gives the net current. The current meters give the commuters' flood and ebb, and a lower bound on the walkers' entropy production from the measured edges (Proposition 6.7 of the theory note). They also give a sign test of handedness: on 29 of 41 footpaths the residual current runs as left-hand traffic would, P = 0.0058.
- **The waves.** The kinematic street of §2.4 runs on Swanston Street. Its blocks and cross streets come from the map, with a signal at each street and none at a lane. Walkers enter at the counted rate. The engine keeps four minutes second by second, and the page plays them: the density as water, walkers moving at the model's speed, queues at red, the lamps of the crossings.
- **The plume.** Two stores picked on the map as a pathway are drawn on it. The first store's visitors spread along the street minute by minute, rose where the need has been raised; an arc joins the stores, as thick as their coupling; and a column rises at the second with those who have bought.

The fits of a hundred tides took two minutes, because calx's QR applied its reflections column by column to a matrix stored by rows. Applied row by row, the factors are the same to the bit (calx's checks hold this) and a fit takes a quarter of the time. The fits now run on every core, and the server warms them when it starts.

### 2.6 The Records

The data library gains manifests and a product builder (`data/products/street_records.py`), for:
- the City of Melbourne's hourly pedestrian counts, from May 2009 to December 2022 and from October 2024, with the sensors' locations;
- its census of land use (CLUE), businesses by year, use and place;
- Open-Meteo's hourly weather at the city's centre;
- Victoria's public holidays (Nager.Date);
- OpenStreetMap's holdings and declared opening hours along five streets of Melbourne and a residential street of Tokyo (the paper's uniform case, Claim 26.1), read by a parser of the `opening_hours` grammar (`data/formats/opening_hours.py`);
- each sensor's counts each way, from the live record.

Melbourne gives the long, real records the paper's order of work needs first. Ginza, the companion paper's case, was read through its declared hours and its rule in the first build, and removed at the user's word. The domain's theory note still discusses that paper, and its crowd model keeps the paper's street as its default.

### 2.8 Architecture, a Domain of Its Own

Buildings were first planned as a module of urban science. At the user's word ("Only architecture independent in domain", after "I think architecture deed can independent. But street/city building? It is highly interdisciplinary domain"), architecture is a domain of its own (`domains/architecture`, planned, field, height 0 on calx, morphe and physis). It holds buildings' form, plans and districts as cell complexes, frontages, the rules that bound them, and their environment. A building is an object with models no other domain holds.

Street and city building is not a domain of its own. The street and the city are urban science's object already: "streets and buildings as open systems". A second domain would hold them twice. The infrastructure's form for an interdisciplinary subject is the bridge, a domain with several fields beneath it. Urban science becomes the bridge of the street and the city when its modules declare architecture, geography or economics. The city composed for a purpose, designing a street or reading a district, is a work's: the Urban Building Workbench, independent among the applications. Its planned plug-in *Architecture* begins when the architecture domain's first module is written.

## 3. Alternatives Considered

- **Implementing the paper as dynamical laws.** Navier–Stokes for crowds, mass action for stores, an energy for streets: this is what "utilize fluid dynamics, chemical dynamics" reads as at first. The paper itself withholds every such law (§2) and sets aside the energetic vocabulary (Claim 27.2). Implementing them as claims would have contradicted the formalism being implemented. They are implemented as representations, labelled so, beside the measurables they would have to answer to.
- **Street oceanography as a module of geography.** Geography holds urban geography (README §4.5) and places on the Earth. But the paper's object is a system of terrain, flows, carried properties and record, not a place. The workbench is also to hold architecture and culture. A domain of its own keeps geography to places and lets the two meet by declaration.
- **The spectral methods in the domain.** Multitaper estimation, harmonic analysis and co-location name no object in nature, and the tides, climate and neuroscience domains need them too. They went to calx and tyche (§2.2), as the user's principle requires (ADR 0042).
- **Only the constructed example of the paper.** Real records exist and are open: Melbourne's counts since 2009, its census of land use since 2002, OpenStreetMap's declared hours. The paper's order of work begins with exactly these cheap instruments (§36.3), so the workbench begins with them.

## 4. Consequences

- The user's examples become exact objects, each named in a view:
  - stores that close at different hours are lines of the week;
  - kinds of people sent to stores at hours are the visits of a contraction;
  - the bridal and ring stores are a catalytic pathway with a coupling length and window.
- The foundations gain spectral and harmonic analysis that the tides domain, climate and neuroscience can use, kinematic waves and transport that traffic and hydrology can use, and the statistical mechanics of Markov processes that chemistry and economics can use.
- The paper's first test is begun on real segments, with its caveat: the sensors were not chosen by its design. Its decisive instrument, carried properties at entry and departure, remains unbuilt, and the workbench says so wherever composition is drawn.
- What the records hold, as found:
  - **Heights.** OpenStreetMap maps the height of 14% of Melbourne's buildings. The rest stand at defaults by kind, drawn paler. The City of Melbourne's building outlines, with measured heights, would replace them.
  - **Bridal stores.** OpenStreetMap tags no bridal store on the streets, and CLUE names 11, by heuristic. The bridal–jeweller pathway therefore runs on stores the user picks, not on a found cluster.
  - **Currents.** On the mean weekday the counting lines take in about 5,200 more people than they let out by 17:00, and by midnight let out more than they took in: arrivals come partly by ways no sensor counts (the trains beneath). The lines are a sample, not a cordon.
  - **Gaps.** Melbourne's counts have a gap from November 2022 to October 2024, and the archive's September 2010 timestamps were repaired by order of record.
- The domain's name, and whether urban geography moves here, are left to the user. The colour is the user's (§2.5).
