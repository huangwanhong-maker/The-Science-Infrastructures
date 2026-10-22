# inference: Where Models Meet Observations

*Status: planned.* Layer 3, the top of `core/`.

## 1. Purpose

This layer estimates states and parameters from observations, together with their uncertainty. Much of the research in this project will happen here, for example:

- orbit determination;
- tidal analysis;
- ocean state estimation;
- numerical weather prediction.

## 2. Planned Modules

| Module | Scope |
|---|---|
| `estimation` | Least squares, maximum likelihood, Bayesian estimation |
| `filtering` | Kalman filters and their nonlinear variants; particle filters |
| `assimilation` | Variational assimilation (3D-Var, 4D-Var), ensemble assimilation, adjoint models |
| `uncertainty` | Uncertainty propagation, sensitivity analysis |
| `learning` | Surrogate models, neural operators, equation discovery. Python frameworks are allowed in this module. |

## 3. Open Question

Variational assimilation needs gradients. Automatic differentiation is an open decision in calx, which would build it as a transformation on its IR.
