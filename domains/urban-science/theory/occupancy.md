# Owners and Occupiers: the Uses of a District's Buildings, Their Vacancy, and Interventions on It

- **Status:** draft
- **Activity:** research
- **Module:** `occupancy` (`occupancy/include/urban/occupancy/`: `model.hpp`, `calibration.hpp`, `planning.hpp`), checked by `occupancy/tests/check_occupancy.cpp`
- **Decisions:** [ADR 0071](../../../development/decisions/0071-agents-and-decisions-for-urban-planning.md) (planned), [ADR 0075](../../../development/decisions/0075-agents-in-the-city.md) (built); the design, [development/ideas/2026-11-04-agents-and-decisions-for-urban-planning.md](../../../development/ideas/2026-11-04-agents-and-decisions-for-urban-planning.md), §1, §2, §10, §11
- **Foundations:** physis (the nested chain and its fit, decision 0020), nomos (`agents`: Gillespie's schedule, logit choice, the ODD+D description), tyche (`decision`: POMDPs, HSVI, components joined by a budget and the Lagrangian bound; `processes`: transition rows by uniformisation), calx (random streams, categorical draws, Nelder and Mead, linear solves); within the domain, `events` (the record, the classifications of uses, the memory) and `affinity` (co-location)
- **Case study:** [research/experiments/2026-11-05-owners-and-occupiers](../../../research/experiments/2026-11-05-owners-and-occupiers/README.md)

This note states the model the module implements, with its definitions and results and their proofs or sources, and says where a statement is a theorem about the model and where it is a finding on the City of Melbourne's record. Names in `code` are those of the headers.

## 1. The Question and the Record

A building's use is seen once a year, at the City of Melbourne's Census of Land Use and Employment (CLUE): its predominant space use, one of CLUE's uses (Office, Retail - Shop, Unoccupied - Unused, …), which fall into families by CLUE's own prefixes (`urban::events::prefix_classification`). The questions are how the uses change, why buildings fall vacant and stay so, and which interventions on vacant space reduce vacancy, when a building's use between censuses is not seen.

The module `events` already reads the record as a *nested chain* on the classification (physis decision 0020): a building half keeps its family for 7.3 years and its use for 2.9 in Event Spacetime's district, 7.7 ± 0.2 and 2.9 ± 0.1 over the Hoddle Grid's 1,535 buildings. This module puts agents on that law, calibrates them to the record, and plans on them.

## 2. The Model

### 2.1 States and Rules

**Definition 2.1** (classification). States $x \in U$ are uses, classified at levels $0 \le \ell \le L$ (the whole, the family, the use: $L = 2$), with a law $\pi$ on $U$, the uses' frequency over the record. $B_\ell(x)$ is the level-$\ell$ block of $x$. The *vacant family* $F_\varnothing$ is one level-1 block; within it, the *lettable* use $v$ ("Unoccupied - Unused") and the uses under works $W = F_\varnothing \setminus \{v\}$.

**Definition 2.2** (owners and occupiers). Each building $i$ of a district has a use $x_i$ and the family $h_i$ it last held occupied (the owner's fit-out, unseen by the census). Its neighbours $N_i$ are its $k$ nearest buildings. Given an affinity $a : U \times U \to \mathbb R$ and a demand $d(t)$, constant within each year, its rules fire in continuous time:

1. *Succession at level $\ell$*, at rate $\lambda_\ell$: the next use is drawn from $B_\ell(x_i)$ by the logit $p(y) \propto \pi(y)\, e^{\beta_n A_i(y)}$, where $A_i(y) = |N_i|^{-1}\sum_{j\in N_i} a(x_j, y)$.
2. *Vacating*, if $x_i$ is occupied, at rate $\nu e^{-\beta_d d(t)}$: the next use is drawn from $F_\varnothing$ by the same logit.
3. *Letting*, if $x_i \in F_\varnothing$, at rate $\eta e^{\beta_d d(t)}$: with probability $\rho$ the next use is drawn from the occupied uses of the family $h_i$, else from all occupied uses, each by the same logit.
4. *The incentive*, if $x_i = v$ and $i$ is enrolled, at the extra rate $\kappa\eta e^{\beta_d d(t)}$, drawn as a letting.

A redraw that returns the use held is not a change. The affinity is $a(u, y) = \ln \mathrm{CLQ}_{F(u)\to F(y)}$, the co-location quotient of $y$'s family about $u$'s (Leslie and Kronenfeld 2011, *Geogr. Anal.* 43, 306; `urban::affinity::affinity`), estimated from one census and clamped to $[-2, 2]$ (`family_affinity`).

*Rationale.* Succession is the nested chain's mechanism: an occupier is replaced, and the replacement is drawn at the level at which the building is let again. Vacating and letting are the owners' and occupiers' second mechanism, a building left with no successor and later let; §4.3 shows the record needs it. The logit is the owner's bounded choice (McFadden 1974; `nomos::agents::logit`): uses that fit the neighbours are favoured by $\beta_n$.

### 2.2 The Baseline

**Proposition 2.2** (the baseline is the nested chain). If $\nu = \eta = 0$, $\beta_n = 0$ and no building is enrolled, then each building's generator is
$$G = \sum_{\ell < L} \lambda_\ell (\Pi_\ell - I), \qquad \Pi_\ell(x, y) = [B_\ell(x) = B_\ell(y)]\, \pi(y)/\pi(B_\ell(x)),$$
the buildings are independent, and the law of a panel of censuses is the product of the chain's transition probabilities $P_\tau(x, y)$ of physis decision 0020.

*Proof.* With $\beta_n = 0$ the logit within $B_\ell(x)$ is $\pi$ restricted and renormalised, which is $\Pi_\ell(x, \cdot)$; rules 2–4 have rate zero. A building's rules read no other building, so the generator of the configuration is the sum of independent single-building generators, and its semigroup their tensor product. ∎

The check holds `generator` to physis's `NestedChain::generator` entry by entry (difference $7 \times 10^{-18}$), and the agent-based model's one-year law from every use, by Gillespie's simulation (20,000 runs a use), to the chain's kernel (largest $|z|$ 2.5).

### 2.3 Three Lumps

**Proposition 2.3** (strong lumpability). Fix the neighbours' uses and the demand. The chain of a building is strongly lumpable onto $\{O, V, W\}$ (occupied, lettable vacancy, works) in the sense of Kemeny and Snell (1960, *Finite Markov Chains*, Theorem 6.3.2): for lumps $A \neq B$, $\sum_{y \in B} G(x, y)$ is the same for every $x \in A$.

*Proof.* Level-0 succession draws from all of $U$ by a law that does not depend on $x$, so its rate into each lump is $\lambda_0 p(B)$ for every $x$. Level-1 succession stays within $x$'s family, hence within its lump if $x$ is occupied; within $F_\varnothing$ it moves between $V$ and $W$ by a law that does not depend on which vacant use is held. Vacating has the same rate from every occupied use and lands in $F_\varnothing$ by a law independent of $x$. Letting has the same rate from every vacant use and lands in $O$ whatever $h_i$ is. The incentive acts on $V$ alone. ∎

The lumped generator (`lumped_generator`) is the planner's model (§5). The check holds the lumps' rates equal across their uses to $10^{-16}$, the modulation on and the channel of vacancy open.

### 2.4 Schedule and Reproducibility

The model is a `nomos::agents` jump model run by Gillespie's direct method (Gillespie 1977), year by year with the year's demand, each year from a stream seeded by the run's seed mixed with the year (`year_seed`, splitmix64's finaliser). One seed gives one history to the bit, and an ensemble's statistics do not depend on the threads that run it (checked). Every firing appends a `Change` aligned with nomos's event log, so a run's history is exact, and `to_record` writes it in urban-science's record (ADR 0067), every event of provenance *simulated*: a simulated history is never read as an observed one.

## 3. The Model's Spectrum

**Proposition 3.1.** Let $f$ be a function of the use, and the chain stationary at $\pi$ (the baseline). Its autocovariance at lag $\tau$ is $C_f(\tau) = \sum_{m \ge 1} c_m e^{-\mu_m \tau}$ with $c_m = \langle f, Q_m f \rangle_\pi$, $Q_m = \Pi_m - \Pi_{m-1}$ and $\mu_m = \lambda_0 + \dots + \lambda_{m-1}$; sampled once a year, its spectrum is $S_f(\omega) = \sum_{m \ge 1} c_m (1 - r_m^2)/(1 - 2 r_m \cos\omega + r_m^2)$, $r_m = e^{-\mu_m}$.

*Proof.* $P_\tau = \sum_m e^{-\mu_m\tau} Q_m$ (decision 0020), $Q_m$ orthogonal projections in $L^2(\pi)$, so $C_f(\tau) = \langle f, P_\tau f\rangle_\pi - \langle f, Q_0 f\rangle_\pi = \sum_{m\ge 1} e^{-\mu_m\tau}\langle f, Q_m f\rangle_\pi$; each term is the autocovariance of an AR(1) of coefficient $r_m$ at integer lags, whose spectrum is the summand. ∎

The vacant family's indicator is a function of the level-1 block, so it lies in one mode ($m = 1$, memory $\lambda_0$); the lettable use's indicator lies in two. The case study draws the closed form beside the record's pooled periodogram; the computation is the study's, since it names no building (*Remark.* It belongs in physis's `nested_chain.hpp` when a second user needs it).

## 4. Calibration by Indirect Inference

### 4.1 The Estimator

The auxiliary model is the one `events` fits to a record, the nested chain, with the record's event rates beside it (the design's §10). Its statistics $\psi$ (`auxiliary`) are ten, over one-year pairs of censuses unless said:

| | statistic | identifies mainly |
|---|---|---|
| $\psi_0, \psi_1$ | $\ln\hat\lambda_0, \ln\hat\lambda_1$, the nested chain fitted to every transition (physis's `fit_rates`) | $\lambda_1$; $\lambda_0$ with $\nu,\eta$ |
| $\psi_2$ | falling vacant: occupied, then vacant | $\nu$ |
| $\psi_3$ | let: vacant, then occupied | $\eta$ |
| $\psi_4$ | falling into works: occupied, then under works | the logit within $F_\varnothing$ |
| $\psi_5$ | returning: vacancies between occupied censuses ending in the family before | $\rho$ |
| $\psi_6$ | the neighbours' pull: a new family's mean affinity to the neighbours, less its mean under $\pi$ | $\beta_n$ |
| $\psi_7, \psi_8$ | slopes of each year's falling vacant and let on $d(y)$ | $\beta_d$ |
| $\psi_9$ | moving between occupied families | $\lambda_0$ |

Simulated panels reproduce the record's structure: each building enters at its first census in its observed use and is observed at its censuses. With $\theta = (\ln\lambda_0, \ln\lambda_1, \ln\nu, \ln\eta, \operatorname{logit}\rho, \beta_n, \beta_d)$ and $\bar\psi(\theta)$ the mean over $S$ panels of common random numbers,
$$\hat\theta = \arg\min_\theta J(\theta), \qquad J(\theta) = \sum_k w_k\big(\psi_k^{\rm obs} - \bar\psi_k(\theta)\big)^2,$$
$w_k$ the inverse variances at the start (Gouriéroux, Monfort and Renault 1993; Smith 1993). The start reads the chain's one-year probabilities backwards (`starting_parameters`). The minimum is found by Levenberg and Marquardt's damped Gauss–Newton steps on the binding function, its Jacobian by central differences on the common random numbers, then polished by Nelder and Mead (`calibrate`). A simplex alone stopped short in seven dimensions on the piecewise constant $J$ (§4.2).

**Proposition 4.1** (Gouriéroux, Monfort and Renault 1993, Proposition 4). Under the model, as the record grows with $S$ fixed, $\hat\theta$ is consistent and asymptotically normal, with
$$\operatorname{Var}\hat\theta = (1 + S^{-1})(G^\top W G)^{-1} G^\top W \Sigma W G (G^\top W G)^{-1},$$
$G = \partial\bar\psi/\partial\theta$ and $\Sigma$ the covariance of $\psi$. Hansen's statistic $\frac{S'}{1+S'}g^\top\Sigma^{-1}g$, $g = \psi^{\rm obs} - \bar\psi(\hat\theta)$, is asymptotically $\chi^2$ on $10 - p$ degrees of freedom with the optimal weight. The module estimates $\Sigma$ from $S' = \max(2S, 10q)$ replicates, since the inverse of a covariance estimated from few replicates is unstable; with the diagonal weight the statistic is a diagnostic. Between nested specifications, the difference of the statistic at the two minima, with the larger model's $\Sigma$, is the distance-metric test.

### 4.2 The Estimator on Records of Known Parameters

The check simulates four records of 300 buildings over 21 years at known $\theta$ and calibrates each on its own seeds ($S = 12$). Every estimate lies within four of its standard errors of the truth (largest $|z|$ 2.0), and each parameter's mean over the four within three standard errors of that mean. $\beta_d$ is estimated high on average (0.62 ± 0.06 against 0.5): a finite-sample bias of the demand statistics, which the case study states. In a probe of the estimator, a simplex alone stalled on one of four such records at a $J$ near three times $J$ at the truth (27.7 against 10.2). That is why the binding function is inverted first.

### 4.3 What the Record Demands

On the Hoddle Grid's 1,535 buildings (33,083 censuses; the case study), three findings shaped the model.

1. **The chain alone cannot hold vacancy.** At its fitted rates the chain alone implies 0.011 of occupied buildings falling vacant in a year and 0.079 of vacancies let. The record has 0.033 and 0.220: vacancy is a short transit between occupiers, about three times as frequent as the family memory implies. The chain alone's $J$ is 11,809. Rules 2 and 3 are the mechanism this calls for.
2. **A vacancy remembers the use before it.** Of the vacancies seen to begin and end between occupied censuses, 54% end in the family that preceded them. A letting drawn afresh returns there 16% of the time, so the letting carries the family held with probability $\rho$ (rule 3). It is estimated at 0.51 ± 0.03.
3. **Neighbours shape the choice of use.** The new family's affinity to the neighbours exceeds its expectation by 0.070 in the record, against 0.024 without the neighbours' term. $\beta_n = 2.08 \pm 0.18$, and dropping it costs 129 on one degree of freedom (distance metric).

The full model fits each statistic within 2.8 simulated standard deviations ($J = 8.3$), but Hansen's statistic (62 on 3 degrees of freedom) rejects it. The largest misfit is falling into works, which the model over-predicts (0.0074 against 0.0060).

*Caveat.* Three properties of the record are outside the model, and the case study shows each:
- **CLUE surveys a property at least every two years.** A use is carried forward between surveys, so a vacancy is rarely seen let after one year (0.06) and often after two (0.48).
- **Vacancy doubled from 2020.** The demand index does not carry the jump ($\beta_d = -0.06 \pm 0.03$), and the pooled rates average the two regimes.
- **Vacancy is contagious.** Falling vacant rises from 0.029 to 0.059 as the share of vacant neighbours grows, and vacancy stays clustered. The neighbours' term acts in the owner's choice of use, not in the rate of vacating, so it reproduces neither. A neighbours' term in the vacating rate is the next modulation to test.

## 5. Interventions on Vacant Space as a Weakly Coupled POMDP

### 5.1 One Building

**Definition 5.1** (the planner's model of a building). The POMDP $(S, A, O, T, Z, R, \gamma)$ on $S = \{O, V, W\}$, a step a year (`building_model`):
- $A$ = {nothing, inspect, incentive}, with costs 0, 0.1 and 1 (incentive-years);
- $T_a = e^{Q_a}$, $Q_a$ the lumped generator of Proposition 2.3 with the building enrolled ($a$ = incentive) or not. The neighbours and the demand are held at the last census: this is the decoupling, declared;
- $O = S \cup \{\bot\}$: the census surveys the building in a year with chance $\rho_s$ and reads its lump, else reads $\bot$. An inspection reads it surely. $\rho_s = 1/2$ is declared from CLUE's survey of every commercial property at least every two years;
- $R(s, a) = -\int_0^1 e^{tQ_a}(s, V)\,dt$, the expected years of lettable vacancy within the year (Simpson's rule, 64 intervals). Works are not the planner's target.

### 5.2 The City

The buildings' sum, joined only by the year's budget $\sum_i c(a_i) \le B$, is a weakly coupled POMDP (the design's Proposition 1.4). For every $\lambda \ge 0$,
$$V^*(b) \le \sum_i V_i^\lambda(b_i) + \frac{\lambda B}{1-\gamma}$$
(Adelman and Mersereau 2008; Meuleau et al. 1998), with $V_i^\lambda$ building $i$'s value with reward $R - \lambda c$. tyche's `lagrangian_bound` minimises it over $\lambda$ with each $V_i^\lambda$ replaced by HSVI's upper bound, and `coupled_actions` ranks the buildings' actions within the budget. Buildings whose lumped models agree to $10^{-12}$ share one model. The check holds the bound above the product's value on three buildings, the product solved by HSVI (brute force, 27 states), and the ranked policy within its budget and below the bound.

### 5.3 Rules and Evaluations

On the Hoddle Grid (the case study, §6.3), the ranked policy reaches 97–98% of the Lagrangian bound's ceiling in the planner's model at budgets of 5 to 40 incentive-years a year, against 74–81% for the simple rules (every vacancy, by neighbours) and about 61% for the longest vacant first. In the agent-based model it averts 11.5 to 29.1 building-years more than every vacancy over 2025–2034 (paired, 1.4 to 4.0 standard errors).

The simple rules each enrol buildings last seen vacant and lettable, while the budget allows:
- every vacancy, in a random order;
- the longest vacant first;
- the most vacant neighbourhood first.

A policy is evaluated in two worlds (`evaluate_decoupled`, `evaluate_abm`):
- **The planner's own decoupled model.** The policy runs for $H$ years and then nothing, valued exactly beyond $H$. This value is a lower bound on $V^*$ in expectation, and the Lagrangian bound lies above it.
- **The agent-based model.** Here the neighbours' feedback is kept, and the years of lettable vacancy are measured from each run's changes, paired against doing nothing on common seeds.

*Caveat.* The values, costs, $\kappa = 1$ and $\rho_s$ are assumptions; the data give the dynamics, not the City's values or an incentive's effect. The planner's model holds the neighbours fixed, so spillovers between buildings appear only in the agent-based model's evaluation.

## 6. What Is Not Claimed

- The model is not a forecast. The record shows a change of regime in 2020 that the model does not hold (§4.3).
- No rents, leases or owners are observed: an owner's choice is a logit fitted to aggregates, not observed decisions.
- The case study is a study of methods on a real place, not a planning proposal.

## 7. References

- Adelman, D. and Mersereau, A. J. (2008), Relaxations of weakly coupled stochastic dynamic programs, *Oper. Res.* 56, 712–727.
- Gillespie, D. T. (1977), Exact stochastic simulation of coupled chemical reactions, *J. Phys. Chem.* 81, 2340–2361.
- Gouriéroux, C., Monfort, A. and Renault, E. (1993), Indirect inference, *J. Appl. Econometrics* 8, S85–S118.
- Hobolth, A. and Jensen, J. L. (2005), Statistical inference in evolutionary models of DNA sequences via the EM algorithm, *Stat. Appl. Genet. Mol. Biol.* 4, 18 (decision 0020's endpoint conditioning).
- Kemeny, J. G. and Snell, J. L. (1960), *Finite Markov Chains*, Van Nostrand.
- Leslie, T. F. and Kronenfeld, B. J. (2011), The colocation quotient: a new measure of spatial association between categorical subsets of points, *Geogr. Anal.* 43, 306–326.
- Marquardt, D. W. (1963), An algorithm for least-squares estimation of nonlinear parameters, *J. SIAM* 11, 431–441.
- McFadden, D. (1974), Conditional logit analysis of qualitative choice behavior, in P. Zarembka (ed.), *Frontiers in Econometrics*, Academic Press, 105–142.
- McKinnon, K. I. M. (1998), Convergence of the Nelder–Mead simplex method to a nonstationary point, *SIAM J. Optim.* 9, 148–158.
- Meuleau, N. et al. (1998), Solving very large weakly coupled Markov decision processes, *AAAI-98*, 165–172.
- Smith, A. A. (1993), Estimating nonlinear time-series models using simulated vector autoregressions, *J. Appl. Econometrics* 8, S63–S84.
- Smith, T. and Simmons, R. (2004), Heuristic search value iteration for POMDPs, *UAI-04*, 520–527.
