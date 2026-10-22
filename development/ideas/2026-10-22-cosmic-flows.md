---
title: Cosmic Flows from Peculiar Velocities
status: seed
date: 2026-10-22
related: [applications/astronomy/celestial-atlas, core/domains/astronomy, core/inference]
---

# Cosmic Flows from Peculiar Velocities

## 1. Idea

A galaxy's velocity has two parts:

- the Hubble flow, H₀d, set by its distance;
- a peculiar velocity, set by the gravity of the structures around it.

Where distances are measured independently of redshift, the difference v − H₀d maps the motions that gravity induces: the infall towards Virgo, and the flow towards the Great Attractor and Laniakea.

## 2. Data and Methods

- **Data.** Cosmicflows-4 (Tully et al. 2023; VizieR J/ApJ/944/94) gives distances for 56 000 galaxies and peculiar velocities for groups. It can be registered as a further `tap` manifest.
- **Display.** Only the radial component of the velocity is observed. The first view colours galaxies by approaching or receding relative to the Hubble flow.
- **Reconstruction.** The full three-dimensional velocity field and the density field behind it come from reconstruction: Wiener filtering, or Bayesian inference with the gravitational potential as a latent field. This belongs in `core/inference`.
- **Representation.** Streamlines of the reconstructed flow, which connect to the atlas's 3D map.

## 3. Questions

- Which reconstruction is simplest to verify? Linear theory with a Wiener filter is a candidate.
- How should the large distance errors (about 20 %) be shown honestly?
