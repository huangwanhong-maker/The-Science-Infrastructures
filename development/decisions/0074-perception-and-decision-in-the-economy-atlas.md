# 0074. Perception and Decision in the Economy Atlas

- **Status:** Accepted, provisionally. Built: tyche's hidden Markov models (its decision 0011), nomos's Kirman model (its decision 0005), the economics domain's `regimes` and `policy`, and the Economy Atlas's Decisions view.
- **Date:** 2026-11-05
- **Decided by:**
  - **The request.** After ADR 0071 the user asked to "visualize such algorithms in economics state perception and decision making, and save galley images. Comprehensive as possible."
  - **What Claude did.** It chose what is perceived (the regime of an economy, from the World Bank's growth of real GDP) and what is decided (a policymaker's support or restraint, a fund's support among economies, a budget on support, two policymakers, and two markets of agents). It wrote the methods these need into the foundations and the domain, the view, its gallery pictures and this record.

## 1. Context

ADR 0071 put POMDPs of every family into tyche's `decision` and agent-based models into nomos's `agents`, and gave them their first use in a city, the maintenance of its surfaces. The Economy Atlas drew none of them.

A POMDP needs a state that is not seen, observations that bear on it, and a model of both. The economics domain had none:
- **The state.** The regime of an economy, read from its growth: Hamilton's (1989) reading of the business cycle as a hidden Markov chain.
- **The fit.** tyche had the panel likelihood of a chain observed at intervals (its decision 0010), but no model whose observations are continuous and noisy, fitted to many sequences.
- **Agents.** nomos's agents held one model, Schelling's.

## 2. Decision

### 2.1 The Foundations

- **tyche `processes`** gained `hidden_markov.hpp` (its decision 0011):
  - Gaussian emissions;
  - the filter, the smoother and Viterbi's path, with missing years allowed;
  - Baum–Welch over many sequences, optionally with one σ for every state.
- **tyche `decision`.** HSVI records its interval before each trial.
- **nomos `agents`** gained Kirman's herding, with the exact law of its birth–death chain (its decision 0005).

### 2.2 The Domain

**`regimes`:** `perceive_regimes` fits three regimes to every economy's growth at once.
- **Data.** The World Bank's NY.GDP.MKTP.KD.ZG, 1961–2024, a profile of the data library's `worldbank-indicators`; aggregates are excluded.
- **Standardisation.** Each economy is read in its own deviations: its median and 1.4826 times its median absolute deviation, clipped at ±6.
- **The model.** One σ for all three regimes, so that they differ by their means alone.
- **The fit.** On 212 economies and 11,208 years:
  - contraction at −3.86σ, normal growth at 0.01σ and a boom at 3.39σ, σ = 1.01;
  - mean durations of 1.5, 22.2 and 1.8 years;
  - 4.5 %, 92.0 % and 3.5 % of the years.
- **The readings.** Each economy's filtered and smoothed laws, its path, and the belief carried a year ahead.

**`policy`:** decisions on that belief, every assumption declared in `PolicyAssumptions` and carried into every answer.
- **The policymaker's POMDP.** Its state is the regime, its observation the year's growth in six bins, and its actions wait, support and restrain.
  - Support moves half of a contraction's persistence to normal growth.
  - Restraint moves 35 % of a boom's.
  - Each has its side effect: overheating a boom, or deepening a contraction.
  - A year is worth −4, 1 or −0.5 by regime, less the action's cost (0.15 for support, 0.2 for restraint), discounted at 0.9.
- **The fund's model.** Two actions, each economy one component of a weakly coupled problem, support costing one unit of a year's budget.
- **Support's cost.** A cost of one for each year of support, for the constrained problem.

### 2.3 The View

**The Decisions view**, beside the Toolbox, in eight sections, each computed by the engine (`engine/decide.hpp`, the request `decide`):
- **Perception.**
  - An economy's growth year by year, with its filtered and smoothed chance of contraction.
  - The regimes and how they follow one another.
  - The world's contractions year by year.
  - The chance of contraction in a year, as a layer of the map.
- **Decision.**
  - The policy over the simplex of beliefs, solved by HSVI, with an economy's beliefs drawn over it.
  - The action on each year's belief.
  - HSVI's interval closing.
  - The value along an edge of the simplex, against the bounds from the states.
  - The assumptions, editable.
- **Search.** POMCP's action values at one belief against the number of simulations, beside HSVI's intervals.
- **A fund.** Support among every economy with a record in a year, under a budget:
  - the Lagrangian bound against its multiplier;
  - the economies chosen, as a layer of the map.
- **A budget on support itself.**
  - The plans for each price of support.
  - The two the budget mixes.
  - The mixture, whose welfare equals the dual bound.
- **Two policymakers.** A Dec-POMDP of two neighbouring economies, each seeing only its own growth:
  - the best joint policy's trees;
  - JESP's equilibrium;
  - the value of acting alone.
- **Herding.** Kirman's market and its law against the exact one.
- **Sorting.** Schelling's lattice, tick by tick.

## 3. Alternatives Considered

- **Regimes from raw growth, pooled.** Rejected after the fit. Pooled in per cent, the regimes separated kinds of economy, fast growers and slow, rather than phases of the cycle. Each economy is now read against its own normal.
- **A deviation for each regime.** Rejected after the fit: the worst regime became turbulence, a mean near zero with a large σ, not contraction.
- **The policy's first defaults** (a contraction worth −2, support costing 0.3 and ending a quarter of its persistence). Rejected: support was never worth its cost, and every decision was to wait. The values are assumptions; the page says so and lets them be changed.
- **The fund's ties broken by the order of the codes.** Rejected after its pictures. In 2020, 42 economies were as sure of contraction as the model can be, and the greedy allocation keeps the order of equal gains, so the budget went alphabetically. Beliefs are now kept to six places, so that equal ones tie exactly. The economies are ordered by the depth of their fall in their own deviations, unclipped. The page states the rule and the number tied.
- **The budget's plans at prices 0, 0.5, 1, ….** Rejected after its picture: the budget binds at λ = 0.088, between the first two, so the two plans mixed were not among those drawn. The prices are now spread geometrically from 0.001, and the two plans are drawn.
- **Kirman's law from the run drawn (4,000 units of time).** Rejected after its picture: the market switches seldom, so one run leaves the two ends unequal by chance. The law is read from a run fifty times as long, within total variation 0.01 of the exact law.
- **The view as tools of the Toolbox.** Rejected: a tool answers one form, while a decision is read against the perception it acts on, across sections that share an economy and a year.

## 4. Consequences

- **Checks.**
  - tyche's `check_tyche_hidden_markov` (7).
  - nomos's `check_nomos_agents`, now 16.
  - economics' `check_economics_policy` (9), on simulated economies:
    - the pooled fit recovers three regimes' means and persistence;
    - every reading's laws sum to 1, and a contraction's mean duration is 1/(1 − A₀₀);
    - standardised, the most probable paths recover each economy's own regimes in 90.7 % of its years (85 % required);
    - the POMDP is well formed, and support and restraint act as declared;
    - certain of contraction it supports, of a boom it restrains, of normal growth it waits;
    - HSVI's interval closes to 0.01;
    - the fund's model keeps wait and support, support costing one a year.

    The fund's bound and the budget's mixture rest on tyche's checks of `structured` and `constrained` (its decision 0009).
- **What the pictures show** (the gallery, §8.1).
  - The United States is in contraction with chance 0.52 in 1982, 0.88 in 2009 and 0.67 in 2020.
  - 41 % of the economies with a record were more likely than not in contraction in 2020, and 18 % in 2009.
  - In 2020 a budget of 20 supports Curaçao, Mauritius, Montenegro and Kosovo first, by the depth of their falls.
  - Alone, neither policymaker finds support worth its cost; together, they support at once after the shock.
- **The assumptions stay assumptions.** What support and restraint do, and what a year is worth, are not estimated from the data. Estimating them would need a record of policies and their effects, which the data library does not hold.
- **Planned.**
  - Policy's effects estimated from a record of policy rates (BIS) and fiscal balances.
  - Regimes of more than one series (growth and inflation).
  - Agents that plan.
