# applications: Works That Compose the Domains

*Status: draft.* This directory is an edge at the top of the architecture. It may depend on everything, and nothing depends on it.

## 1. Purpose

An application is a complete program that people use. Where a domain holds the knowledge of one science, an application is organised by purpose, and draws on as many domains as its purpose needs. The Earth Atlas, for instance, shows the work of geography, meteorology, oceanography, geophysics, geodesy and climate on one globe.

Every application says which of the three activities it serves:

- **learning:** it shows how something works, and invites the reader to change it and watch;
- **use:** it gives numbers and pictures that can be relied on;
- **research:** it is the instrument of an open question.

## 2. Organisation

- **One directory per application,** named for what it is, not for a field.
- **Dependencies.** An application uses the domains, the foundations and the shared edges (`data/`, `visualization/`) only through their public interfaces. Its README names the domains it draws on.
- **Its own server and tests.** Each web application has its own local server and port, and its own tests of what is particular to it; the checks of the science live with the domains.
- **Its own data recipe.** `recipe.json` lists the datasets the application fetches and the products it builds, in order. `data/acquisition/recipe.py` runs it from the command line and is the quick start; the server reports it to the page at `/setup` and runs it on request, so that the page can fetch what is missing.
- **Applications do not depend on one another.** Two applications may link to each other's views, as the two atlases do, but neither imports the other.

## 3. Index

| Application | Serves | Draws on | Status |
|---|---|---|---|
| [Celestial Atlas](celestial-atlas/) — [online](http://173.212.213.207:9200/) | use, learning | metrology, astronomy, geodesy | draft |
| [Earth Atlas](earth-atlas/) — [online](http://173.212.213.207:9201/) | use, learning | geography, meteorology, oceanography, geophysics, geodesy, astronomy, climate | draft |
| [Brain Atlas](brain-atlas/) | use, learning | neuroscience | draft |
| [Spiking Simulator](spiking-simulator/) | use, learning | neuroscience (its C++ engine) | draft |
| [Complex Atlas](complex-atlas/) | learning, research | the foundations morphe and physis (its C++ engine) | draft |
| [History Complex](history-complex/) | research, learning | history (its C++ engine), the foundation morphe | draft |
| [Genome Atlas](genome-atlas/) | use, learning, research | genomics and evolution (its C++ engine, kept alive as workers of a cluster), the foundations calx and morphe | draft |
| Orbital sandbox | learning | astronomy | planned |
| Ephemeris calculator | use | metrology, astronomy | planned |
| Map projections | learning | geography | planned |
| Terrain and watersheds | learning | geography | planned |
| Remote sensing | learning | geography | planned |
| Land-use atlas | use | geography | planned |
| Rotating tank | learning | oceanography, meteorology | planned |
| Ocean-profile explorer | use | oceanography | planned |
| Weather maps | learning | meteorology | planned |
| Atmospheric soundings | learning | meteorology | planned |
| Foucault pendulum | learning | geodesy | planned |
| Tide-predicting machine | learning | tides | planned |
| Tide predictor | use | tides | planned |
| Sight-reduction trainer | learning | navigation | planned |
| Route planner | use | navigation, geography | planned |
| Air–sea exchange | learning | air-sea | planned |
| Energy balance | learning | climate | planned |
| Milankovitch cycles | learning | climate, astronomy | planned |

The two atlases run online on the project's server (the links above), and on any computer with Git and Python in a few lines: the Quick Start of the [repository's README](../README.md#quick-start), or of each atlas's own.

The planned applications in brief:

- **Orbital sandbox:** two-body and n-body motion, Kepler's laws, perturbations.
- **Ephemeris calculator:** tables of positions in any time scale and reference frame.
- **Map projections:** projections and their distortions.
- **Terrain and watersheds:** terrain models, slopes, drainage networks and basins.
- **Remote sensing:** satellite imagery, its correction, band combinations and land-cover classification.
- **Land-use atlas:** population, settlement and land use, mapped and analysed.
- **Rotating tank:** shallow water on a rotating sphere, with Kelvin and Rossby waves; the same dynamics hold in the atmosphere.
- **Ocean-profile explorer:** temperature and salinity profiles from Argo floats.
- **Weather maps:** geostrophic wind, fronts and cyclones on real analyses.
- **Atmospheric soundings:** radiosonde profiles, thermodynamic diagrams and stability.
- **Foucault pendulum:** the rotation of the Earth made visible.
- **Tide-predicting machine:** a digital version of Kelvin's machine, each constituent a wheel and the tide their sum.
- **Tide predictor:** predictions for any tide gauge with a harmonic analysis.
- **Sight-reduction trainer:** a position from star sights by spherical trigonometry.
- **Route planner:** great-circle and rhumb-line routes on charts.
- **Air–sea exchange:** wind stress, heat and gas fluxes, and the growth of wind waves.
- **Energy balance:** the Earth's radiation budget and the greenhouse effect, from a single-layer model.
- **Milankovitch cycles:** the Earth's orbital elements, insolation and the ice ages.
