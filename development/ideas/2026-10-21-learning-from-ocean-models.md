---
title: Learning from Operational Ocean Models
status: seed
date: 2026-10-21
related: [lib/physis/continuum, lib/physis/chemistry, core/domains/oceanography, core/inference/assimilation]
---

# Learning from Operational Ocean Models

## 1. Context

Fluid and chemical dynamics are central to oceanography. Global observing systems, such as Argo floats and satellite altimetry, do not simulate anything themselves. Their data are combined with numerical ocean models through data assimilation. Examples of the resulting products are:

- the ECCO state estimate;
- Copernicus Marine;
- NOAA's operational forecasts.

## 2. Idea

Study the open models and reimplement their core algorithms in our own layered framework, validating against the originals. Candidates:

- **Ocean models:**
  - MITgcm, which has an adjoint, used by ECCO;
  - MOM6, NEMO, ROMS and FESOM;
  - Oceananigans.jl, which runs on GPUs.
- **The TEOS-10 equation of state for seawater**, from the GSW library.
- **Simple biogeochemical models**, such as NPZD (nutrient, phytoplankton, zooplankton, detritus). These would be the first coupling of chemistry and biology with transport.

Reimplementing an algorithm from its published description is always allowed. Copying code requires that the model's licence be compatible with ours.

## 3. Path

Each step has standard test cases to validate against.

1. Shallow water on the sphere.
2. Tides with real bathymetry.
3. A primitive-equation model on a small domain.
4. Biogeochemical tracers.
