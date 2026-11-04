# 0056. The Methods of the Social Sciences, nomos, and a Toolbox of Quantitative Economics

- **Status:** Accepted, provisionally
- **Date:** 2026-11-04 (the work begun 2026-10-24)
- **Amends:** [ADR 0044](0044-value-in-space-and-time.md), [ADR 0046](0046-the-economy-atlas-data-and-plugins.md)
- **Decided by:** the user, after ADR 0046: "I think function is so limite. And color can be lighter, and more positive. And I think can add quantantative economics toolbox." Asked which foundations and which groups of tools, the user chose two new foundations and "Econometrics, Time series, Inequality and indices, Trade, games and models, All you think is needed. If something better independent to the economics package, then make new foundation or some else. Low level things better in cpp or python." The user approved the name *nomos* (Greek *νόμος*, law, custom; the *-nomy* of economy). Claude designed the library, the domain's modules, the toolbox and the map's new functions, and wrote this record; two subagents wrote tyche's first modules ([ADR 0048](0048-the-inference-foundation-tyche.md)) and nomos's to the boundaries stated here.

## 1. Context

The Economy Atlas drew value carried between places and through time (ADR 0044), took one's own data and plug-ins (ADR 0046), but answered no ordinary question of quantitative economics: how unequal incomes are, whether poorer economies grow faster, what a regression of life expectancy on income says, whether the pound follows a random walk, how far distance holds trade back. The methods those questions need are of two kinds. Some are general statistics (distributions, regression, time series): tyche, the inference foundation (ADR 0048). Others belong to the practice of the social sciences but name no society: measures of inequality, poverty and concentration, index numbers, games, matching, social choice, and the econometric designs (panels, instruments, differences in differences, Poisson regression with absorbed effects). The map of the mathematics placed these in a planned foundation for the quantitative social sciences (§7), whose name and boundary with inference were open decision 3.

## 2. Decision

### 2.1 nomos, the Foundation for the Methods of the Social Sciences

nomos is a foundation library in its own repository, [`huangwanhong-maker/nomos`](https://github.com/huangwanhong-maker/nomos), attached as the submodule `foundations/nomos`. It depends on tyche (for regression and probability), calx (linear algebra, the simplex method, random streams) and morphe (for the graphs and orders its planned modules need); nothing in the foundations depends on it. Its first six modules (its decision 0003), 242 checks:

| Module | Content |
|---|---|
| `measurement` | inequality (the weighted Gini coefficient, the Lorenz curve, the generalised entropy family, Atkinson's indices, quantile shares, the Palma ratio), poverty (Foster–Greer–Thorbecke, Watts), concentration (Herfindahl–Hirschman, concentration ratios, numbers equivalents), specialisation (location quotients, dissimilarity, Grubel–Lloyd), and their decompositions by groups |
| `indices` | the index numbers of the CPI Manual, bilateral and elementary, direct and chained |
| `games` | finite games: best responses, pure and mixed equilibria by support enumeration, iterated dominance, zero-sum games by linear programming |
| `matching` | deferred acceptance with capacities from either side, stability, top trading cycles |
| `choice` | profiles of ballots, scoring rules, instant runoff, Condorcet methods (Copeland, Schulze, Kemeny–Young) |
| `econometrics` | fixed effects absorbed by alternating projections; panels (within, between, random effects, Hausman); two-stage least squares with the first stage and overidentification; differences in differences; average marginal effects; Poisson pseudo-maximum likelihood with two sets of absorbed effects |

**The boundary with tyche.** tyche estimates and tests in general; nomos holds the measures and designs of the social sciences and estimates on tyche's regressions. Absorbing fixed effects is general, but stays in nomos, its established home, until another science needs it (nomos's decision 0003).

*Remark.* Two checks of the panel estimators failed when the work stopped, and the fault was in one estimator and one test. Random effects took its classical covariance from the transformed regression's residual variance, which the effects inflate when they are correlated with the regressors, so that Hausman's statistic could come out negative; it now carries σ̂²_ε (ZᵀZ)⁻¹ on Swamy and Arora's scale, the within estimator's, so that V_FE − V_RE is positive semi-definite (Baltagi 2013, ch. 2). The test of the variance components simulated a common trend the one-way model leaves out; it now simulates the model it tests.

### 2.2 The Economics Domain Gains Growth, the Short Run and Trade

Three modules, header-only C++, on the foundations (the domain declares tyche and nomos):

| Module | Content |
|---|---|
| `growth` | the Solow–Swan model (Cobb–Douglas in closed form, CES by calx's Runge–Kutta and Brent's method), its steady state, convergence and golden rule; growth accounting with Törnqvist shares; the perpetual inventory of capital |
| `macro` | IS–LM: the equilibrium, its multipliers, the curves and aggregate demand |
| `trade` | the partners of trade (partner concentration, the intensity index, which is the location quotient of the bilateral matrix, and two-way trade, whose world share is circulation's reciprocal exchange exactly); the gravity of trade by PPML with exporter and importer effects (nomos) or with sizes (tyche), on great-circle distances (physis) |

The domain's suite runs 109 checks (68 before).

### 2.3 The Toolbox

A fifth view of the Economy Atlas, the Toolbox, holds 27 tools in seven groups:

| Group | Tools | On |
|---|---|---|
| statistics | summary and correlation | tyche |
| econometrics | linear regression; logit, probit and Poisson; panels; instrumental variables; differences in differences; the convergence of incomes | tyche, nomos |
| time series | autocorrelation; unit roots; ARIMA and forecasts; trend and cycle; vector autoregressions; cointegration | tyche |
| inequality and indices | inequality; poverty; concentration; index numbers | nomos |
| trade | the gravity of trade; the concentration of partners; trade intensity | economics, on nomos and tyche |
| growth and macroeconomics | the Solow model; IS–LM; growth accounting | economics |
| games and choice | two-player games; zero-sum games; stable matching; voting rules | nomos |

- **A tool is a declaration** (`src/tools.js`): the roles it gives the columns of a table, its parameters, and an example on the library's data or a table pasted, so that every tool runs at once. Its form is built from the declaration.
- **The engine answers it** (`engine/tools/`, a header for each group): it reads the columns the roles name, calls the foundations and the domain, and answers blocks: figures, tables, lines with bands, bars, scatters with a line, matrices, notes, and values for the economies, which become a layer of the map as a plug-in's do.
- **Tables.** The library's tables (`tables.py`): the economies by year, 2000–2024, their trade beside the World Bank's indicators (GDP per head, population, life expectancy, fertility, urban share, forest share, access to electricity, CO₂ and greenhouse gases per person, the Gini index, density), with the World Bank's regions and income groups; the currencies by month, long and wide. One's own tables are the datasets of kind `table` (CSV or JSON); a table may also be pasted into the form. The recipe fetches the World Bank's indicators and makes the tables.
- **Examples tell the truth of their data.** Where the library's data cannot carry a valid design (an instrument, a treatment), the example is simulated and says so: sixty markets whose cost of supply moves price but not demand, twenty regions half treated from 2020. Least squares is shown beside two-stage least squares so that the bias the instruments remove is seen.

### 2.4 The Map Does More

- A place is found by its name or code.
- **Play** goes through the years, each year's decomposition fetched while the one before is shown.
- The economies are coloured by the flows' own measures (net inflow, exports, imports, the share circulating) or by any of the World Bank's indicators (in the map's year, or the latest of the five before it; on a log scale where the values span decades), besides α and the plug-ins' and tools' values.
- A place's card shows its trade over the years, and the indicator the map is coloured by.
- The year's places and net flows are saved as CSV, and the map as SVG with its styles.

### 2.5 Lighter and More Positive Colours

The palette of ADR 0046, lightened at the user's word and validated again with the dataviz method (OKLab ΔE ≥ 8 between neighbours under protanopia and deuteranopia; text ≥ 4.5:1 on the paper). Day is the page's own theme; night is chosen, no longer taken from the system.

| Job | ADR 0046 | Now |
|---|---|---|
| paper | ivory | `#FFFDF7`, a lighter ivory |
| gold, the main colour | `#C9A85C` | `#DDB65A` |
| value sent out | steel blue `#3A6EA5` | sky blue `#4A86C5` |
| value drawn in, rose | `#C4506E` | `#E0607E` |
| driven flow | indigo `#4B4F8F` | `#5B63B0` |
| circulation round the holes | `#A98A2C` | `#C39A2E`, its contrast relieved by labels and tables |
| a magnitude | the neutral to deep gold | the neutral through light gold to a warm amber (`#A8701A`) |

## 3. Alternatives Considered

- **One foundation, the social methods as modules of inference.** Offered and declined (ADR 0048, §3).
- **The methods in the economics domain.** Rejected: a Gini coefficient or a voting rule names no economy, and history and the planned domains of political science and sociology need games, choice and matching.
- **The toolbox in Python, on statsmodels and SciPy.** Rejected for the engine: the infrastructure computes in C++ (ADR 0003), and its numbers are checked in its own suites. The tables are made in Python, as the data library's products are.
- **Dense dummies for gravity's effects.** Rejected: 400 columns for 40,000 pairs in every iteration, where the projections cost a few passes over the cells (nomos's decision 0003).

## 4. Consequences

- Open decision 3 of CLAUDE.md is closed: the inference foundation is tyche (ADR 0048), the foundation of the quantitative social sciences nomos, and their boundary is §2.1's.
- The order of the foundations: calx and morphe, physis on both, tyche on the three, nomos on tyche; commits follow it. `development/tools/check-dependencies.mjs` and `check-alone.mjs` know nomos.
- The economics domain declares tyche and nomos; its checks run alone with them.
- Not yet: nomos's planned modules (psychometrics, surveys, networks, mechanisms, agents); GDP for places of one's own that are not the library's economies (gravity's sizes then fall back to total trade); the tools' answers kept between visits.
