---
title: The Research Arc of Prediction, a Plan
status: seed
date: 2026-10-23
related: [development/decisions/0031-observations-scores-and-hazards.md, development/ideas/2026-10-22-primitive-equation-model.md, foundations/calx, foundations/physis]
---

# The Research Arc of Prediction, a Plan

Deferred by the user's decision of 2026-10-23 ([ADR 0031](../decisions/0031-observations-scores-and-hazards.md)): it rests on two decisions CLAUDE.md §8 lists as the user's, the name of the inference foundation and the C++ toolchain. Written down so that it can start when they are taken.

## 1. Idea

The Earth Atlas now observes, keeps and scores. What it cannot yet do is the loop at the heart of numerical weather prediction: a model of its own, started from an analysis that it makes itself from observations, run forward, and scored, cycle after cycle. Each part exists in a small form; the arc joins them and moves them where they belong.

## 2. Steps

1. **The inference foundation** (name to be chosen; the fourth foundation, beside calx, morphe and physis). Its first modules are those assimilation needs, written once for every science:
   - linear-Gaussian estimation: optimal interpolation, the Kalman filter, 3D-Var as the minimisation of a quadratic cost;
   - ensembles: the ensemble Kalman filter (stochastic and square-root forms), localisation, inflation;
   - diagnostics: innovation statistics, the Desroziers relations for estimating error variances.

   The optimal interpolation of the Earth Atlas's plug-in moves here, and the plug-in calls it.
2. **The C++ toolchain** (Phase 0): on the server, g++ 13 and CMake; on Windows, WSL, the MSVC Build Tools, or LLVM. A first C++ module, with its JavaScript stand-in kept as the check, under the same tests.
3. **Upper-air observations**: radiosondes (IGRA2, NOAA NCEI), and aircraft reports if an open source is found. Without them, a 500 hPa model cannot be assimilated: surface reports constrain the surface.
4. **A primitive-equation model** ([the idea](2026-10-22-primitive-equation-model.md)): dry first, in σ coordinates, spectral in the horizontal, semi-implicit; the baroclinic wave of Jablonowski and Williamson (2006) as its first test.
5. **The cycle**: every six hours, the model's forecast as background, observations assimilated by 3D-Var or an ensemble Kalman filter, the analysis as the next start; scored against GFS's and ECMWF's analyses by the existing verification, and run by the existing schedule.
6. **Learned corrections**, studied against the cycle: whether a correction learned from the archive's errors improves the model, as NOAA's AIGFS suggests, and what the physics gives that data alone do not.

## 3. Open Questions

- The foundation's name, and whether it lives in its own repository (as the other foundations do).
- Whether the first C++ core is calx's spectral transforms, or the model's dynamics directly.
- How far an ensemble can go on this hardware (8 cores and 29 GB on the server; an RTX 3050 Ti with 4 GB on the laptop).
