# 0048. The Inference Foundation, tyche

- **Status:** Accepted
- **Date:** 2026-11-04 (the library founded 2026-10-24)
- **Decided by:** the user, who asked for "a quantitative economics toolbox" in the Economy Atlas and, offered one foundation or two for the methods it needs, chose two: "If something better independent to the economics package, then make new foundation or some else. Low level things better in cpp or python." The user approved a foundation for inference and its name, *tyche* (Greek *τύχη*, chance). Claude drew its boundaries and wrote its first modules and this record; tyche's own decisions 0001–0003 record its founding in detail.

## 1. Context

The map of the mathematics ([docs/mathematics.md](../../docs/mathematics.md), §6; [ADR 0035](0035-the-mathematical-foundations.md)) designed a fourth foundation for chance, *inference*: probability, stochastic processes, statistics, inference and learning, on calx, morphe and physis. Its name was open decision 3 of CLAUDE.md, and its first module was to be the one a domain first needed. That domain was economics: the Economy Atlas's toolbox needed distributions with accurate tails, the description of samples, least squares with the covariances applied work uses, generalised linear models, and the time series of macroeconomics. The genomics domain had meanwhile written the statistics it needed itself (mixtures, negative binomial models, false discovery rates), for want of the foundation.

## 2. Decision

### 2.1 The Library

tyche is a foundation library in its own repository, [`huangwanhong-maker/tyche`](https://github.com/huangwanhong-maker/tyche), attached as the submodule `foundations/tyche`. It holds what is general about chance: laws and methods that name no object in nature. It depends on calx (linear algebra, special functions, optimisation, random streams), morphe (exact and combinatorial structures) and physis (dynamical systems, for assimilation and the uncertainty of models); none of them depends on it.

```
calx    morphe
   \    /   \
   physis    |
      \      |
       tyche
```

### 2.2 Its First Modules

| Module | Content |
|---|---|
| `probability` | distributions (normal, Student's t, χ², F, gamma, beta, binomial and Poisson) with their densities, distribution functions, survival functions and quantiles accurate in the tails, and p-values |
| `descriptive` | moments, quantiles of every Hyndman–Fan type, weighted statistics, covariance and correlation matrices, ranks with Spearman's and Kendall's correlations, summaries |
| `regression` | least squares, weighted or not, by Householder QR, with classical, robust (HC0–HC3) and clustered covariances, Wald tests and predictions; generalised linear models by IRLS (Gaussian, binomial with logit or probit, Poisson) with model-based, robust or clustered covariances |
| `time_series` | autocovariance, the ACF and PACF by Durbin–Levinson, Ljung–Box; AR by least squares and Yule–Walker; ARIMA by conditional least squares, with forecasts; the augmented Dickey–Fuller test with MacKinnon's p-values, KPSS, Engle and Granger's cointegration; the Hodrick–Prescott filter; VARs, their order, impulse responses and forecasts |

`estimation` and `testing`, which the map places between `probability` and `regression`, were not separate modules at first: the likelihoods and tests these four needed were written where they were used.

### 2.3 Its Boundary with the Methods of the Social Sciences

The methods that belong to the practice of the social sciences (measurement and index numbers, games, matching, social choice, and the econometric designs: panels, instruments, differences in differences, Poisson regression with absorbed effects) are a foundation of their own, nomos, on tyche ([ADR 0056](0056-nomos-and-a-toolbox-of-quantitative-economics.md)). tyche estimates and tests in general; nomos holds the designs and measures of the social sciences and estimates on tyche's regressions.

## 3. Alternatives Considered

- **One foundation, with the social methods as its modules.** Offered, and declined by the user: the social methods have their own literature, their own checks (axioms of indices and of voting rules, equilibria) and a reach beyond statistics (games and choice are combinatorial), and the two libraries can grow at their own pace.
- **The statistics in the economics domain.** Rejected by the principle of reuse: a regression names no economy, and the genomics, urban-science and planetary-science domains needed the same methods.
- **Borrowing statsmodels or R.** For the page, rejected: the infrastructure's computation is C++ (ADR 0003), its numbers must be checked in its own suites, and a statistics written here is how its methods are learnt. Python work may still use SciPy and statsmodels directly (docs/mathematics.md §10.3).

## 4. Consequences

- Open decision 3 of CLAUDE.md is closed for inference; with nomos (ADR 0056) it is closed whole.
- tyche grew after its founding by other sessions' work, each recorded in its own decisions and in the infrastructure's records: Gaussian conditioning (its decision 0004, [ADR 0054](0054-bodies-dynamics-and-knowledge.md), when it was published at the user's word); spectral inference, harmonic analysis and point patterns (0005, [ADR 0060](0060-streets-as-open-systems.md)); change points, penalised regression and factorisations, fine-mapping, variance components and multiple testing, and learning (0006–0008, [ADR 0064](0064-abnormal-genomes-their-description-and-expression.md)).
- The statistics the genomics domain wrote for itself before tyche existed may now move down to it, as their second users appear.
