---
title: Agents and Decisions for Urban Planning. Agent-based models, decisions under partial observation, and four models of a city
status: seed
date: 2026-11-04
related: [tyche (a new module, decision), nomos (agents), calx (random, optimisation), physis (statistical_mechanics), urban-science (events, landscape, flow), architecture, the data library, the Urban Science Workbench]
---

# Agents and Decisions for Urban Planning

*This note designs two general methods and the urban models built on them, before most of them are built.*

- **Part I** states the formalism: decisions under partial observation (POMDPs and their variants), and agent-based models.
- **Part II** places them in the infrastructure. It also sets their interfaces and checks, and the order in which they are built.
- **Part III** designs four urban planning models on Melbourne's data: maintenance and heritage, housing and residential choice, pedestrians and streets, and land use and vacancy.

The user asked for both methods "for modern urban planning", reusable and kept apart from the applications. They chose, on 2026-11-04:
- agent-based models in nomos's planned `agents`, and POMDPs in a new tyche module `decision`;
- all four urban problems;
- all four families of POMDP: exact and point-based, online Monte Carlo, factored with mixed observability, and multi-agent and constrained.

They also asked that review papers be searched. The references in §13 were checked against Crossref, the publishers or DBLP. An entry marked *(unverified)* was not checked.

## 0. The Position in Brief

### 0.1 Eight Claims

1. **Two general methods, in the foundations.**
   - A decision under partial observation names no city, so it belongs in tyche. It is statistical decision theory, Bayes's rule carrying a belief that is acted on (MSC 62C, 90C40).
   - An agent-based model names no society, so it belongs in nomos, beside games, matching and choice, whose agents it animates.
   - urban-science holds only what names a street, a building or a household.
2. **One tuple, many variants.** Every variant is the POMDP tuple restricted, extended or multiplied (§1.6). The MDP sees its state; the MOMDP sees part of it; a constrained POMDP adds costs with budgets. A Dec-POMDP multiplies actions and observations over agents; a stochastic game multiplies the rewards too. A Bayes-adaptive POMDP puts the model's parameters into the state. Each solver states the class it solves and its guarantee.
3. **An answer is an interval.** Every offline solver returns a lower bound, the value of a plan it can carry out, and, where it can, an upper bound (HSVI's closes on $V^*(b_0)$ to a declared ε). A planning answer is a policy with an interval for its value, never a bare number. This is the infrastructure's rule that numbers carry their accuracy, applied to decisions.
4. **An agent-based model is a Markov process with local generators.** A schedule says how the local rules compose:
   - synchronously, as a product of kernels;
   - in a random sequence, as a composition;
   - asynchronously in continuous time, as a sum of generators.

   The last is a jump process, so a small model's law is physis's master equation exactly and can be checked against it. Synchronous updating creates artefacts (Huberman and Glance 1993); a model states its schedule.
5. **Reproducible by identity.** Every draw is keyed by (seed, rule, agent, step) on calx's counter-based streams. A run then does not depend on the number of threads, nor on the order in which agents are visited where the schedule is order-free.
6. **A simulated history is a record.** An agent-based model writes its events in the schema of urban-science's record (ADR 0067), with the provenance *simulated*. Every analysis of a record applies to simulations as to observations: the complex, its homology, the memories by level, the rates of events. Calibration compares the two by indirect inference, the nested chain of physis decision 0020 its auxiliary model (§6.4).
7. **Planning acts on the record.** A planner's actions are events of the record (interventions). A POMDP's hidden state is the record's unobserved part: the uses between censuses, the condition between inspections. Its observations are the record's observations, and its belief is the forward message of the record's contraction (physis decision 0020).
8. **A city's assets are weakly coupled.** Thousands of surfaces, trees and buildings each evolve on their own and are joined only by a budget. A Lagrangian relaxation separates them into one small POMDP each (Meuleau et al. 1998; Adelman and Mersereau 2008). It gives an upper bound on the city's optimum and a policy that ranks the assets. This is what makes maintenance planning at the scale of a city tractable.

### 0.2 What Is Not Claimed

- No POMDP or MDP for heritage conservation, housing policy or urban land use was found in the literature that the survey would cite (§13.4). The models of Part III are proposals, and their value is to be shown.
- The data hold almost no history of inspections, rents of dwellings or owners' decisions (§12). What the models cannot take from data they declare as assumptions, with their provenance *assumed*.

# Part I. The Formalism

## 1. Decisions Under Partial Observation

### 1.1 The Tuple and the Belief

**Definition 1.1 (POMDP).** A *partially observable Markov decision process* is a tuple $(S, A, O, T, Z, R, \gamma, b_0)$ whose parts are:
- finite sets $S, A, O$ of states, actions and observations;
- a transition kernel $T(s' \mid s, a)$;
- an observation kernel $Z(o \mid s', a)$ of the state entered;
- a reward $R(s, a)$;
- a discount $\gamma \in [0, 1]$, with $\gamma = 1$ only over a finite horizon;
- an initial law $b_0 \in \Delta(S)$.

A *history* is $h_t = (a_0, o_1, \dots, a_{t-1}, o_t)$. A *policy* maps histories to (laws on) actions, and its value is $V^\pi(b_0) = \mathbb E^\pi[\sum_t \gamma^t R(s_t, a_t)]$ (Kaelbling, Littman and Cassandra 1998).

**Theorem 1.2 (the belief is sufficient; Åström 1965; Smallwood and Sondik 1973).** The law $b_t = P(s_t \mid h_t) \in \Delta(S)$ is a sufficient statistic of the history for optimal control. It is carried forward by Bayes's rule,
$$\tau(b, a, o)(s') = \frac{Z(o \mid s', a) \sum_s T(s' \mid s, a)\, b(s)}{P(o \mid b, a)},$$
so that a POMDP is an MDP on $\Delta(S)$, the *belief MDP*, with reward $\rho(b, a) = \sum_s b(s) R(s, a)$.

*Remark.* In the language of Event Spacetime (ADR 0067), $\tau$ is one step of the forward contraction of physis's `event_chains`: an observation is a diagonal piece, and a transition a cell's propagator. A belief is the message the contraction carries.

### 1.2 Values

**Theorem 1.3 (PWLC; Smallwood and Sondik 1973).** For every finite horizon $H$, $V_H(b) = \max_{\alpha \in \Gamma_H} \alpha \cdot b$ for a finite set $\Gamma_H \subset \mathbb R^S$. Each $\alpha$ is the value of a conditional plan, labelled by the plan's first action.

The Bellman operator carries $\Gamma$ to $\Gamma' = \bigcup_a \bigoplus_o \Gamma^{a,o}$, where
$$\Gamma^{a,o} = \{ R(\cdot, a)/|O| + \gamma \textstyle\sum_{s'} T(s' \mid \cdot, a) Z(o \mid s', a)\, \alpha(s') : \alpha \in \Gamma \}$$
and $\oplus$ is the cross sum. It is a $\gamma$-contraction in the sup norm. So $\|V_H - V^*\|_\infty \le \gamma^H \max|R| / (1 - \gamma)$.

**Complexity.** These limits justify the design's insistence on bounds:
- finite-horizon POMDPs are PSPACE-complete (Papadimitriou and Tsitsiklis 1987);
- the infinite-horizon questions are undecidable (Madani, Hanks and Condon 2003).

### 1.3 Bounds From the States

All three are due to Hauskrecht (2000):
- **QMDP.** The MDP's action values as vectors (Littman, Cassandra and Kaelbling 1995). An upper bound: a decision maker who saw the state could do no worse.
- **The fast informed bound (FIB).** The fixed point of $\alpha_a(s) = R(s, a) + \gamma \sum_o \max_{a'} \sum_{s'} T(s' \mid s, a) Z(o \mid s', a)\, \alpha_{a'}(s')$. Tighter than QMDP, still above $V^*$.
- **Blind policies.** Each repeats one action whatever is seen: $\alpha_a = R(\cdot, a) + \gamma T_a \alpha_a$. A lower bound.

So $\text{QMDP} \ge \text{FIB} \ge V^* \ge \text{blind}$, pointwise on $\Delta(S)$.

### 1.4 Solvers Offline

| Solver | What it computes | Guarantee |
|---|---|---|
| Exact value iteration by incremental pruning (Cassandra, Littman and Zhang 1997), Lark's filter (White 1991) | $\Gamma_H$ | exact; $V^*$ within $\gamma^H \max\|R\|/(1-\gamma)$ |
| PBVI (Pineau, Gordon and Thrun 2003) | backups at a growing set $B$ | a lower bound; error $\le (R_{\max} - R_{\min})\,\delta_B / (1-\gamma)^2$ on reachable beliefs |
| Perseus (Spaan and Vlassis 2005) | randomised backups on a sampled $B$ | a lower bound, non-decreasing on $B$ |
| HSVI (Smith and Simmons 2004, 2005) | a lower and an upper bound refined along trials | $V^*(b_0) \in [\underline V(b_0), \overline V(b_0)]$ with width $\le \varepsilon$, after finitely many trials |
| SARSOP (Kurniawati, Hsu and Lee 2008) | HSVI's bounds on the optimally reachable beliefs | the same, faster (planned) |

Every vector made by a point backup of plans' values is itself a plan's value (Pineau, Gordon and Thrun 2003). So every lower bound above can be carried out as a policy.

### 1.5 Solvers Online

For state spaces too large to sweep, such as a district's cells:
- **The generative model.** A simulator $G(s, a) \to (s', o, r)$ replaces the tables.
- **POMCP** (Silver and Veness 2010). UCB1 over the tree of histories, with the belief as particles.
- **DESPOT** (Somani et al. 2013; Ye et al. 2017). Planning over $K$ sampled scenarios, regularised, with a bound sensitive to the size of its output.
- **POMCPOW** (Sunberg and Kochenderfer 2018). Weighted particles where observations are continuous. They prove that POMCP with double progressive widening converges to a suboptimal policy there.

The survey of online methods is Ross et al. (2008).

### 1.6 The Variants as One Family

Each variant is the tuple with a part restricted, extended or multiplied:

| Variant | Change to the tuple | Its structure | Solvers here |
|---|---|---|---|
| MDP | $Z$ reveals $s$ | values on $S$ | value iteration |
| MOMDP (Ong et al. 2010; Araya-López et al. 2010) | $S = X \times Y$, $X$ seen | beliefs on $Y$ given $x$ | the point-based solvers on $\{x\} \times \Delta(Y)$ |
| Factored POMDP | $S = \prod_i S_i$, $T$ a product of local kernels | a dynamic Bayesian network | online solvers on the generative model; products of marginals |
| Weakly coupled POMDPs (Meuleau et al. 1998; Adelman and Mersereau 2008; Whittle 1988) | $N$ factored components joined only by a constraint on their joint action each step | a Lagrangian separates them | one small POMDP per component and per multiplier |
| Constrained POMDP (Altman 1999; Isom, Meyn and Braatz 2008; Kim et al. 2011; Poupart et al. 2015) | costs $C_k(s, a)$ with budgets $\mathbb E[\sum \gamma^t C_k] \le c_k$ | Lagrangian duality; optimal policies may randomise | dual ascent over a point-based solver; CC-POMCP (Lee et al. 2018) online |
| ρ-POMDP (Araya et al. 2010) | a reward $\rho(b)$ of the belief | if $\rho$ is convex, so is the value | the value of inspection |
| Bayes-adaptive POMDP (Ross et al. 2011) | $T, Z$ unknown, their counts in the state | learning while planning | planned |
| Dec-POMDP (Bernstein et al. 2002; Oliehoek and Amato 2016) | $A = \prod_i A_i$, $O = \prod_i O_i$, each agent acting on its own observations | NEXP-complete over a finite horizon | exhaustive search for small horizons; JESP (Nair et al. 2003) |
| POSG (Hansen, Bernstein and Zilberstein 2004) | a reward for each agent | equilibria of policies | in nomos, on its games (planned) |
| I-POMDP (Gmytrasiewicz and Doshi 2005) | beliefs over the others' models | nested beliefs | planned, in nomos |

**Proposition 1.4 (Lagrangian bound for weakly coupled POMDPs).** Let components $i = 1, \dots, N$ have independent transitions and observations. Let them be joined only by $\sum_i c_i(a_i) \le B$ at every step. For every $\lambda \ge 0$,
$$V^*(b) \le \sum_i V_i^\lambda(b_i) + \frac{\lambda B}{1 - \gamma},$$
where $V_i^\lambda$ is the value of component $i$ alone with reward $R_i - \lambda c_i$, and $b = \otimes_i b_i$.
*Proof sketch:* relax the constraint into the objective; the relaxed problem separates. This is Adelman and Mersereau's (2008) result, carried from MDPs to beliefs with independent observations.

The policy that ranks components by their marginal gain in $V_i^\lambda$ at the minimising $\lambda$ is feasible after rounding. It is the asset-ranking policy of §10.

### 1.7 The Multi-Agent Boundary

A cooperative Dec-POMDP has one reward. Its solution is one joint policy, an optimisation, so it is tyche's. A POSG has a reward per agent. Its solution is an equilibrium, so it is nomos's, which holds the games. The two modules meet in nomos `agents`: an agent can plan with tyche's solvers against its belief about the others, and the population's play is checked against nomos's equilibria. Bernstein et al.'s (2002) NEXP-completeness bears on the agent-based models: their agents should be best responders or boundedly rational, not jointly optimal.

## 2. Agent-Based Models

### 2.1 Definition

**Definition 2.1 (ABM).** An *agent-based model* is a tuple $(I, X, E, N, \mathcal R, \sigma)$ whose parts are:
- a population $I_t$ of agents, which may change by births and deaths;
- a state space $X$ for each agent;
- an environment state $e \in E$, on a space (a graph, a lattice, the plane, a cell complex);
- neighbourhoods $N_i(\xi) \subset I_t$ read from the space;
- rules $r \in \mathcal R$, each a kernel $K_r(x_i' \mid x_i, x_{N_i}, e)$ that changes one agent, or one agent and its environment, from what it reads;
- a schedule $\sigma$.

The configuration $\xi_t = ((x_i)_{i \in I_t}, e_t)$ is then a Markov process.

### 2.2 Schedules

1. **Synchronous.** Every agent reads $\xi_t$ and all write $\xi_{t+1}$: a product kernel. Order-free. Its artefacts (oscillations, false synchrony) are known (Huberman and Glance 1993).
2. **Random sequential.** In each tick the agents are visited in a random order, each reading what the earlier ones wrote: a composition of single-agent kernels.
3. **Asynchronous in continuous time.** Rule $r$ fires at agent $i$ at rate $\lambda_r(\xi)$: the generator $L = \sum_i \sum_r \lambda_r(\xi)\, (\xi^{i,r} - \xi)$. This is a jump process (Gillespie 1977).

**Proposition 2.2.** Under schedule 3 with finitely many configurations, the law of $\xi_t$ solves the master equation $\dot p = p L$ (physis `master_equation`). For a population of $n$ agents whose rates depend only on counts, the scaled counts converge to the solution of the mean-field ODE as $n \to \infty$ (Kurtz 1970).

Both are checks of the engine: exact for small $n$, asymptotic for large.

### 2.3 Decisions of Agents

An agent's rule may be any of these:
- **A rule** (a threshold, a lookup).
- **A discrete choice.** $P(j) \propto \exp(\beta u_j)$ over options $j$ (McFadden 1974), nomos's econometrics estimating $u$.
- **A learning rule.**
- **A plan.** The agent solves a POMDP about what it cannot see (tyche `decision`), for example a household not knowing a dwelling's true condition, or an owner not knowing the demand for a use.

ODD+D (Müller et al. 2013) says how each is to be described.

### 2.4 Description, Records, Inference

- **ODD.** Every model fills ODD 2020 (Grimm et al. 2006, 2010, 2020), with its decision section from ODD+D. The engine keeps the description as a structure beside the model, not as prose apart from it.
- **The event log.** A run writes each rule's firing as an event: its time, its agent, what it touched, its kind. A domain maps the log to its record (§0.1, claim 6).
- **Ensembles.** Replicates by stream identity, with the observables' means and errors.
- **Calibration.** Three methods:
  - pattern-oriented modelling (Grimm et al. 2005);
  - approximate Bayesian computation and synthetic likelihood (Hartig et al. 2011; van der Vaart et al. 2015);
  - indirect inference (Gouriéroux, Monfort and Renault 1993; Smith 1993), in which the same auxiliary model is fitted to the observed record and to simulated ones, and the model's parameters are chosen to make the two fits agree.

  Sensitivity analysis follows ten Broeke, van Voorn and Ligtenberg (2016). The general procedures (ABC, indirect inference) are tyche's (`bayes`, `estimation`, planned); the engine gives them its simulator.

# Part II. The Engineering

## 3. Placement

| Unit | Holds | Depends on |
|---|---|---|
| calx `random` | `categorical.hpp`: one categorical draw by inversion; Walker's alias table (Vose's construction), for many draws from one law | — |
| tyche `decision` (new; tyche decision 0009) | the tuple and beliefs; the MDP; the bounds; exact value iteration; the point-based solvers; simulation; then online solvers, structure (MOMDP, factored, weakly coupled), constraints, Dec-POMDPs | calx `linalg`, `optimisation`, `random` |
| nomos `agents` (planned there since its decision 0003; nomos decision 0004) | populations, spaces, rules, schedules, streams by identity, event logs, observables and ensembles, the ODD description; agents that choose and agents that plan | tyche `decision`, calx `random`, morphe (graphs) |
| urban-science, four new modules (ADR 0071) | `maintenance`, `housing`, `movement` (planned since ADR 0060), `occupancy` | nomos `agents`, tyche `decision`, and its own `events`, `landscape`, `flow` |

The scope of tyche's CLAUDE.md gains *statistical decision theory: decisions under uncertainty, sequential and partially observed*. The domain adds nomos to the foundations declared in its `module.json`.

## 4. Interfaces

### 4.1 tyche `decision`

- **`pomdp.hpp`.**
  - `Pomdp` (dense $T$, $Z$, $R$) and `validate`;
  - `predict`, `observation_law` and `update` (Bayes's rule), and `expected_reward`;
  - `solve_mdp`, which is ε-optimal (Puterman 1994, Theorem 6.3.1), and `solve_mdp_horizon`.
- **`alpha.hpp`.**
  - `AlphaVector` and `ValueFunction`;
  - the three bounds `qmdp`, `fast_informed_bound` and `blind_lower_bound`;
  - `prune` (Lark's filter on calx's simplex) and `prune_dominated`;
  - `exact_backup` and `exact_value_iteration`;
  - `point_backup`.
- **`point_based.hpp`.**
  - `lookahead`, the exact search of the belief tree;
  - `pbvi` and `perseus`;
  - `UpperBound` (the sawtooth) and `hsvi`.
- **`simulate.hpp`.** `evaluate`: the Monte Carlo value of a policy of beliefs, with its standard error and truncation bound.
- **Planned.**
  - `online.hpp`: the generative-model concept, particle beliefs, POMCP and DESPOT.
  - `structured.hpp`: MOMDPs, factored models, and weakly coupled models with their Lagrangian.
  - `constrained.hpp`: constrained POMDPs by dual ascent, and CC-POMCP.
  - `multiagent.hpp`: Dec-POMDPs, by exhaustive policy trees and JESP.

### 4.2 nomos `agents`

- **`model.hpp`.** An agent model as a class: agents' states, rules with rates or probabilities, a space, and a schedule (synchronous, random sequential, Gillespie).
- **`space.hpp`.** Graphs (morphe), lattices with Moore and von Neumann neighbourhoods, and the plane with a grid index for radius queries.
- **`run.hpp`.** Stepping, the event log, observables, and ensembles of replicates by stream identity.
- **`describe.hpp`.** The ODD+D description as a structure.
- **`choice.hpp`.** Logit choice, threshold rules, and planning agents (through tyche `decision`).

## 5. Checks

| Check | Reference |
|---|---|
| Exact value iteration equals the exact lookahead on the belief tree, every horizon to 5 (the tiger, discounted and not; random POMDPs) | brute force (no shared code) |
| Every vector pruning keeps has a witness; pruning changes no value | the definition |
| The tiger at discount 0.75: $V^*(½) = 1.933439$ | pomdp.org's `tiger.aaai.POMDP`, its solution as the R package `pomdp` documents it (secondary) |
| QMDP ≥ FIB ≥ $V^*$ ≥ blind at 300 beliefs | Hauskrecht 2000 |
| PBVI's, Perseus's and HSVI's lower bounds below $V^*$; HSVI's upper bound above; its interval closing on $V^*(b_0)$ | the exact iterate and its error bound |
| A simulated policy earns its value within four standard errors | Monte Carlo |
| Dec-Tiger, horizons 2–6: $V^* = -4.000000, 5.190812, 4.802755, 7.026451, 10.381625$; the broadcast channel | Oliehoek et al. 2013, Table 2 |
| RockSample[4,4] about 18.0 and [7,8] about 21.27, as bounds | Smith and Simmons 2004; Kurniawati et al. 2008 |
| An ABM under Gillespie's schedule: its law against the master equation, for a small voter or SIR model | physis `master_equation` (exact) |
| Its mean field against the ODE, as $n$ grows | Kurtz 1970 |
| Schelling's model: segregation rising with the tolerance threshold; an ensemble's statistics independent of the number of threads | Schelling 1971 |

The values from published tables (Dec-Tiger, the broadcast channel) are exact optima, so they are fixtures. RockSample's, Tag's and Hallway's are simulated rewards of computed policies, so they are tested as bounds, not as equalities. The GPL tools (APPL/SARSOP, the MADP toolbox, AI-Toolbox) may run as external oracles in fixture scripts, never linked. POMDPs.jl and pomdp_py are MIT.

## 6. Stages

| Stage | Content | State |
|---|---|---|
| D1 | tyche `decision`: the tuple, beliefs, MDP, bounds, exact value iteration, PBVI, Perseus, HSVI, simulation; calx's categorical draws | written, checked (26 + 7) |
| D2 | online: generative models, particle beliefs, POMCP (with a leaf value); DESPOT and POMCPOW planned | POMCP written, checked (3) |
| D3 | structure and constraints: MOMDPs, components joined by a budget with the Lagrangian bound (grouped, weighted), constrained POMDPs by duality with finite-state controllers | written, checked (8 + 9) |
| D4 | several decision makers: Dec-POMDPs (exhaustive for small horizons, exact best responses, JESP), with Dec-Tiger's published optima | written, checked (11) |
| A1 | nomos `agents`: the engine and its checks (master equation, Kurtz, the voter martingale, Schelling) | written, checked (14) |
| P1 | tyche `processes`: Markov chains observed at discrete times (the panel likelihood), for the deterioration of surfaces | written, checked (5) |
| A2 | calibration: indirect inference and ABC in tyche, the engine as their simulator | later |
| U1 | `maintenance` on the Hoddle Grid's 9,630 rated road surfaces | written |
| U2–U4 | `housing`, `movement`, `occupancy` | next, each on its data |

# Part III. Four Models of a City

The models are urban-science's. Each is a question, its state, its agents or decision maker, its data and what the data cannot give. The datasets are those the survey of 2026-11-04 verified (CC BY 4.0 unless stated). The "Hoddle box" is the data library's box about the central grid.

## 7. Maintenance and Heritage (`maintenance`)

**The question.** Which streets, footpaths, furniture and trees to inspect, repair or renew each year under a budget, and which places to protect, when their condition is known only from the last inspection.

**The model.** A weakly coupled constrained POMDP (§1.6):
- **Component.** One asset. Its hidden condition is a class 1–5 (CoM's rating scale: 1 failed, 5 good).
- **Deterioration.** A Markov chain whose rate depends on age and design life, as in Madanat and Ben-Akiva (1994) and Papakonstantinou and Shinozuka (2014).
- **Actions.** Do nothing, inspect (a noisy observation of the class), maintain, renew.
- **Coupling.** A budget each year.
- **Heritage.** Protection adds a constraint, not a reward: a protected place may not be demolished, and its renewal costs more. Its value of information is a ρ-POMDP reward (Memarzadeh and Pozzi 2016).

**Data.**
- `road-segments-with-surface-type`: 9,630 rated surfaces in the box, each with its condition rating 1–5, evaluation date, installation date and design life.
- Street furniture: 3,688 rated assets, with the same fields.
- The urban forest's trees: 3,332 in the box, with a coarse useful-life class and planting dates.
- `building-permits`: repairs, demolitions and heritage permits, with estimated costs.
- Protection: the Victorian Heritage Register, already in the registry, and Vicmap's Heritage Overlay, 1,201 polygons in Melbourne, gazetted 1999–2025.

**What the data cannot give.**
- No asset has a history of inspections; each has only its latest rating.
- The deterioration law, the noise of an inspection, and costs other than the permits' estimates are assumptions, declared.
- The manifest should fetch the surfaces at every refresh, so that a panel of inspections accumulates.

## 8. Housing and Residential Choice (`housing`)

**The question.** How households sort among the dwellings of the centre by rent, type and place, and what changes in supply (new towers, student apartments) or in assistance do to rents and to who lives where.

**The model.** An agent-based model, after the reviews of Huang et al. (2014) and the microsimulation of Waddell (2002):
- **Agents.** Households, each with an income class, a size and a tenure.
- **Dwellings.** CLUE's residential dwellings by building, typed House/Townhouse, Residential Apartments or Student Apartments.
- **Choice.** Logit over vacant dwellings of a utility in rent, type, distance and neighbourhood (§2.3).
- **Rents.** Formed by excess demand by type and area.
- **Sorting.** Schelling-type, as the limiting check (Schelling 1971; Benenson, Omer and Hatna 2002).
- **Planning.** The planner's decision, which supply to permit where under uncertain demand, is an MOMDP whose seen part is the stock and whose hidden part is demand.

**Data.**
- CLUE `residential-dwellings`: the stock, 2002–2024, by building (30,526 dwellings in 2002, 114,150 in 2024).
- The Development Activity Monitor (registered): additions by bedrooms.
- Homes Victoria's quarterly *moving annual rents by suburb*: medians of new lettings. These are the calibration targets.
- ABS Census 2021 G37 and G40 at SA1: tenure and rent, 350 SA1s in Melbourne City, as cross-sections.

**What the data cannot give.**
- No rent or move of a single dwelling, and no vacancy between censuses.
- Preferences are fitted to aggregates, never observed.

## 9. Pedestrians and Streets (`movement`)

**The question.** How walkers distribute over the network of footpaths and crossings, and how street regimes (closures, events, signal phasing) move them, when only some sensors count them.

**The model.**
- **Agents.** Walkers, each with an origin, a destination and a route chosen by a logit over path costs on the pedestrian network: length, grade, crossings and signals (Hoogendoorn and Bovy 2004). Their arrivals are calibrated to the counters by hour and direction.
- **At the level of the street.** The domain's existing kinematic model (`flow`) stays the macroscopic reference. Social force (Helbing and Molnár 1995; Chen et al. 2018) is a later, finer scale.
- **The regime decision.** A POMDP: the hidden state is the crowd on unsensed links, the observations the counts at sensors, the actions a regime per street (open, closed, phased), the reward walkers' time and safety. It is solved online (§1.5) on the agent-based model as its generative model.

**Data.**
- The pedestrian counts (registered): 74 sensors in the box, by direction.
- `pedestrian-network`: 85,326 links (2019).
- `footpath-steepness`.
- `victorian-traffic-signals`: 96 sites in the box, and their configuration sheets.
- The AIRS transport activity counts on William Street, La Trobe Street and Swanston Street.
- The event permits (registered).

**What the data cannot give.**
- No origin–destination data and no signal timing logs.
- No record of street closures after 2018, and the counts' gap from November 2022 to October 2024.
- Route choice and the timing of each regime are assumptions.

## 10. Land Use and Vacancy (`occupancy`)

**The question.** How the uses of the district's buildings change and why they fall vacant, and which interventions (incentives for vacant space, permits for change of use, zoning) reduce vacancy, when the use is seen only at each census.

**The model.** It joins Event Spacetime's memory to agents and to a decision.
- **Agents.** Owners and occupiers on the district's complex, the buildings and their frontages (ADR 0067).
- **The baseline law.** The fitted nested chain (physis decision 0020: a use half kept for 2.9 years, its family for 7.3).
- **What modulates it.** A building's neighbours' uses (`affinity`'s co-location) and demand.
- **The planner's decision.** A weakly coupled POMDP over buildings:
  - the hidden state is each building's use between censuses;
  - the observations are the censuses;
  - the actions are interventions that change the redraw rates of vacant space, under a budget.
- **Calibration.** Indirect inference: the auxiliary model is the nested chain together with the event rates, fitted to the observed record and to simulated ones (§2.4).

**Data.**
- CLUE buildings and establishments, already in the registry.
- `floor-space-by-use-by-block`: 603 blocks × 40 uses, 2002–2024, with five vacancy columns. Small cells are suppressed, so column sums are lower bounds.
- `building-permits`: 4,129 changes of use, with costs and dates, as text to be geocoded.
- Vicmap Planning's zones: the uses each block admits.

**What the data cannot give.**
- No commercial rents, leases or owners.
- Owners' decision rules are assumptions, chosen among those the calibration does not reject.

## 11. What the Urban Models Share

- **Records.** All four write their simulated events into urban-science's record schema, so that the workbench's Event Spacetime can draw a simulated history beside the observed one. A simulated record is never drawn as an observed one: its provenance says so.
- **Places.** Places are named only through the data library's naming policy (`data/products/naming.json`, ADR 0069).
- **The workbench.** A plug-in, *Planning*, would draw policies and their intervals on the map. It is not designed here; the methods come first, as the user asked.

## 12. Gaps in the Data

These are the declared assumptions the models must carry (from the survey):
- **Maintenance.** No inspection histories, no condition of private or heritage fabric, no costs beyond permits' estimates.
- **Housing.** No dwelling-level rents or moves, vacancy only every five years, no open record of short stays.
- **Pedestrians.** No origin–destination data, no signal timings, no closures after 2018.
- **Land use.** No rents, leases or owners; suppressed small cells.

A registry fix the survey found: the City of Melbourne now gives its pedestrian counts the licence CC BY 4.0. The two pedestrian manifests still say "licence null".

## 13. References

### 13.1 POMDPs and Their Solvers

- Åström, K. J. (1965), Optimal control of Markov processes with incomplete state information, *J. Math. Anal. Appl.* 10, 174–205, doi:10.1016/0022-247X(65)90154-X.
- Smallwood, R. D. and Sondik, E. J. (1973), The optimal control of partially observable Markov processes over a finite horizon, *Oper. Res.* 21, 1071–1088, doi:10.1287/opre.21.5.1071.
- Sondik, E. J. (1978), The optimal control of partially observable Markov processes over the infinite horizon: discounted costs, *Oper. Res.* 26, 282–304, doi:10.1287/opre.26.2.282.
- Papadimitriou, C. H. and Tsitsiklis, J. N. (1987), The complexity of Markov decision processes, *Math. Oper. Res.* 12, 441–450, doi:10.1287/moor.12.3.441.
- Madani, O., Hanks, S. and Condon, A. (2003), On the undecidability of probabilistic planning and related stochastic optimization problems, *Artif. Intell.* 147, 5–34, doi:10.1016/S0004-3702(02)00378-8.
- Kaelbling, L. P., Littman, M. L. and Cassandra, A. R. (1998), Planning and acting in partially observable stochastic domains, *Artif. Intell.* 101, 99–134, doi:10.1016/S0004-3702(98)00023-X.
- Cassandra, A. R., Kaelbling, L. P. and Littman, M. L. (1994), Acting optimally in partially observable stochastic domains, *AAAI-94*, 1023–1028 (the tiger problem).
- Cassandra, A. R., Littman, M. L. and Zhang, N. L. (1997), Incremental pruning: a simple, fast, exact method for partially observable Markov decision processes, *UAI-97*, 54–61.
- White, C. C. (1991), A survey of solution techniques for the partially observed Markov decision process, *Ann. Oper. Res.* 32, 215–230 *(unverified)*.
- Littman, M. L., Cassandra, A. R. and Kaelbling, L. P. (1995), Learning policies for partially observable environments: scaling up, *ICML-95*, 362–370, doi:10.1016/B978-1-55860-377-6.50052-9.
- Hauskrecht, M. (2000), Value-function approximations for partially observable Markov decision processes, *JAIR* 13, 33–94, doi:10.1613/jair.678.
- Pineau, J., Gordon, G. and Thrun, S. (2003), Point-based value iteration: an anytime algorithm for POMDPs, *IJCAI-03*, 1025–1031.
- Spaan, M. T. J. and Vlassis, N. (2005), Perseus: randomized point-based value iteration for POMDPs, *JAIR* 24, 195–220, doi:10.1613/jair.1659.
- Smith, T. and Simmons, R. (2004), Heuristic search value iteration for POMDPs, *UAI-04*, 520–527.
- Smith, T. and Simmons, R. (2005), Point-based POMDP algorithms: improved analysis and implementation, *UAI-05*, 542–549.
- Kurniawati, H., Hsu, D. and Lee, W. S. (2008), SARSOP: efficient point-based POMDP planning by approximating optimally reachable belief spaces, *RSS IV*, 65–72, doi:10.15607/RSS.2008.IV.009.
- Shani, G., Pineau, J. and Kaplow, R. (2013), A survey of point-based POMDP solvers, *JAAMAS* 27, 1–51, doi:10.1007/s10458-012-9200-2.
- Ross, S., Pineau, J., Paquet, S. and Chaib-draa, B. (2008), Online planning algorithms for POMDPs, *JAIR* 32, 663–704, doi:10.1613/jair.2567.
- Silver, D. and Veness, J. (2010), Monte-Carlo planning in large POMDPs, *NIPS 23*, 2164–2172.
- Somani, A., Ye, N., Hsu, D. and Lee, W. S. (2013), DESPOT: online POMDP planning with regularization, *NIPS 26*.
- Ye, N., Somani, A., Hsu, D. and Lee, W. S. (2017), DESPOT: online POMDP planning with regularization, *JAIR* 58, 231–266, doi:10.1613/jair.5328.
- Sunberg, Z. and Kochenderfer, M. (2018), Online algorithms for POMDPs with continuous state, action, and observation spaces, *ICAPS* 28, 259–263, doi:10.1609/icaps.v28i1.13882.
- Ong, S. C. W., Png, S. W., Hsu, D. and Lee, W. S. (2010), Planning under uncertainty for robotic tasks with mixed observability, *IJRR* 29, 1053–1068, doi:10.1177/0278364910369861.
- Araya-López, M., Thomas, V., Buffet, O. and Charpillet, F. (2010), A closer look at MOMDPs, *ICTAI 2010*, 197–204, doi:10.1109/ICTAI.2010.101.
- Araya, M., Buffet, O., Thomas, V. and Charpillet, F. (2010), A POMDP extension with belief-dependent rewards, *NIPS 23*, 64–72.
- Ross, S., Pineau, J., Chaib-draa, B. and Kreitmann, P. (2011), A Bayesian approach for learning and planning in partially observable Markov decision processes, *JMLR* 12, 1729–1770.
- Lauri, M., Hsu, D. and Pajarinen, J. (2023), Partially observable Markov decision processes in robotics: a survey, *IEEE T-RO* 39, 21–40, doi:10.1109/TRO.2022.3200138.
- Kurniawati, H. (2022), Partially observable Markov decision processes and robotics, *Annu. Rev. Control Robot. Auton. Syst.* 5, 253–277, doi:10.1146/annurev-control-042920-092451.
- Puterman, M. L. (1994), *Markov Decision Processes: Discrete Stochastic Dynamic Programming*, Wiley.
- Kochenderfer, M. J., Wheeler, T. A. and Wray, K. H. (2022), *Algorithms for Decision Making*, MIT Press.

### 13.2 Several Decision Makers, Constraints and Weak Coupling

- Bernstein, D. S., Givan, R., Immerman, N. and Zilberstein, S. (2002), The complexity of decentralized control of Markov decision processes, *Math. Oper. Res.* 27, 819–840, doi:10.1287/moor.27.4.819.297.
- Oliehoek, F. A. and Amato, C. (2016), *A Concise Introduction to Decentralized POMDPs*, Springer, doi:10.1007/978-3-319-28929-8.
- Nair, R., Tambe, M., Yokoo, M., Pynadath, D. and Marsella, S. (2003), Taming decentralized POMDPs, *IJCAI-03*, 705–710.
- Seuken, S. and Zilberstein, S. (2007), Memory-bounded dynamic programming for DEC-POMDPs, *IJCAI-07*, 2009–2015.
- Hansen, E. A., Bernstein, D. S. and Zilberstein, S. (2004), Dynamic programming for partially observable stochastic games, *AAAI-04*, 709–715.
- Ooi, J. M. and Wornell, G. W. (1996), Decentralized control of a multiple access broadcast channel: performance bounds, *35th IEEE CDC*, 293–298, doi:10.1109/CDC.1996.574318.
- Oliehoek, F. A., Spaan, M. T. J., Amato, C. and Whiteson, S. (2013), Incremental clustering and expansion for faster optimal planning in Dec-POMDPs, *JAIR* 46, 449–509, doi:10.1613/jair.3804.
- Gmytrasiewicz, P. J. and Doshi, P. (2005), A framework for sequential planning in multi-agent settings, *JAIR* 24, 49–79, doi:10.1613/jair.1579.
- Altman, E. (1999), *Constrained Markov Decision Processes*, Chapman & Hall/CRC, doi:10.1201/9781315140223.
- Isom, J. D., Meyn, S. P. and Braatz, R. D. (2008), Piecewise linear dynamic programming for constrained POMDPs, *AAAI-08*, 291–296.
- Undurti, A. and How, J. P. (2010), An online algorithm for constrained POMDPs, *ICRA 2010*, 3966–3973, doi:10.1109/ROBOT.2010.5509743.
- Kim, D., Lee, J., Kim, K.-E. and Poupart, P. (2011), Point-based value iteration for constrained POMDPs, *IJCAI-11*, 1968–1973.
- Poupart, P., Malhotra, A., Pei, P., Kim, K.-E., Goh, B. and Bowling, M. (2015), Approximate linear programming for constrained partially observable Markov decision processes, *AAAI* 29, 3342–3348, doi:10.1609/aaai.v29i1.9655.
- Lee, J., Kim, G.-H., Poupart, P. and Kim, K.-E. (2018), Monte-Carlo tree search for constrained POMDPs, *NeurIPS 31*.
- Meuleau, N., Hauskrecht, M., Kim, K.-E., Peshkin, L., Kaelbling, L. P., Dean, T. and Boutilier, C. (1998), Solving very large weakly coupled Markov decision processes, *AAAI-98*, 165–172 *(unverified)*.
- Adelman, D. and Mersereau, A. J. (2008), Relaxations of weakly coupled stochastic dynamic programs, *Oper. Res.* 56, 712–727 *(unverified)*.
- Whittle, P. (1988), Restless bandits: activity allocation in a changing world, *J. Appl. Probab.* 25A, 287–298 *(unverified)*.

### 13.3 Agent-Based Models

- Heppenstall, A., Crooks, A., Malleson, N., Manley, E., Ge, J. and Batty, M. (2021), Future developments in geographical agent-based models: challenges and opportunities, *Geogr. Anal.* 53, 76–91, doi:10.1111/gean.12267.
- Heppenstall, A. J., Crooks, A. T., See, L. M. and Batty, M., eds. (2012), *Agent-Based Models of Geographical Systems*, Springer, doi:10.1007/978-90-481-8927-4.
- Parker, D. C., Manson, S. M., Janssen, M. A., Hoffmann, M. J. and Deadman, P. (2003), Multi-agent systems for the simulation of land-use and land-cover change: a review, *Ann. AAG* 93, 314–337, doi:10.1111/1467-8306.9302004.
- Matthews, R. B., Gilbert, N. G., Roach, A., Polhill, J. G. and Gotts, N. M. (2007), Agent-based land-use models: a review of applications, *Landsc. Ecol.* 22, 1447–1459, doi:10.1007/s10980-007-9135-1.
- Filatova, T. (2015), Empirical agent-based land market, *CEUS* 54, 397–413, doi:10.1016/j.compenvurbsys.2014.06.007.
- Huang, Q., Parker, D. C., Filatova, T. and Sun, S. (2014), A review of urban residential choice models using agent-based modeling, *Environ. Plan. B* 41, 661–689, doi:10.1068/b120043p.
- Benenson, I., Omer, I. and Hatna, E. (2002), Entity-based modeling of urban residential dynamics: the case of Yaffo, Tel Aviv, *Environ. Plan. B* 29, 491–512, doi:10.1068/b1287.
- Waddell, P. (2002), UrbanSim: modeling urban development for land use, transportation, and environmental planning, *JAPA* 68, 297–314, doi:10.1080/01944360208976274.
- Schelling, T. C. (1971), Dynamic models of segregation, *J. Math. Sociol.* 1, 143–186, doi:10.1080/0022250X.1971.9989794.
- Helbing, D. and Molnár, P. (1995), Social force model for pedestrian dynamics, *Phys. Rev. E* 51, 4282–4286, doi:10.1103/PhysRevE.51.4282.
- Helbing, D., Farkas, I. and Vicsek, T. (2000), Simulating dynamical features of escape panic, *Nature* 407, 487–490, doi:10.1038/35035023.
- Chen, X., Treiber, M., Kanagaraj, V. and Li, H. (2018), Social force models for pedestrian traffic – state of the art, *Transp. Rev.* 38, 625–653, doi:10.1080/01441647.2017.1396265.
- Duives, D. C., Daamen, W. and Hoogendoorn, S. P. (2013), State-of-the-art crowd motion simulation models, *Transp. Res. C* 37, 193–209, doi:10.1016/j.trc.2013.02.005.
- Hoogendoorn, S. P. and Bovy, P. H. L. (2004), Pedestrian route-choice and activity scheduling theory and models, *Transp. Res. B* 38, 169–190, doi:10.1016/S0191-2615(03)00007-9.
- Grimm, V. et al. (2006), A standard protocol for describing individual-based and agent-based models, *Ecol. Model.* 198, 115–126, doi:10.1016/j.ecolmodel.2006.04.023.
- Grimm, V. et al. (2010), The ODD protocol: a review and first update, *Ecol. Model.* 221, 2760–2768, doi:10.1016/j.ecolmodel.2010.08.019.
- Grimm, V. et al. (2020), The ODD protocol for describing agent-based and other simulation models: a second update, *JASSS* 23(2), 7, doi:10.18564/jasss.4259.
- Müller, B. et al. (2013), Describing human decisions in agent-based models: ODD + D, *Environ. Model. Softw.* 48, 37–48, doi:10.1016/j.envsoft.2013.06.003.
- Grimm, V. et al. (2005), Pattern-oriented modeling of agent-based complex systems: lessons from ecology, *Science* 310, 987–991, doi:10.1126/science.1116681.
- Hartig, F., Calabrese, J. M., Reineking, B., Wiegand, T. and Huth, A. (2011), Statistical inference for stochastic simulation models: theory and application, *Ecol. Lett.* 14, 816–827, doi:10.1111/j.1461-0248.2011.01640.x.
- van der Vaart, E., Beaumont, M. A., Johnston, A. S. A. and Sibly, R. M. (2015), Calibration and evaluation of individual-based models using approximate Bayesian computation, *Ecol. Model.* 312, 182–190, doi:10.1016/j.ecolmodel.2015.05.020.
- Thiele, J. C., Kurth, W. and Grimm, V. (2014), Facilitating parameter estimation and sensitivity analysis of agent-based models: a cookbook using NetLogo and R, *JASSS* 17(3), 11, doi:10.18564/jasss.2503.
- ten Broeke, G., van Voorn, G. and Ligtenberg, A. (2016), Which sensitivity analysis method should I use for my agent-based model?, *JASSS* 19(1), 5, doi:10.18564/jasss.2857.
- Huberman, B. A. and Glance, N. S. (1993), Evolutionary games and computer simulations, *PNAS* 90, 7716–7718 *(unverified)*.
- Gillespie, D. T. (1977), Exact stochastic simulation of coupled chemical reactions, *J. Phys. Chem.* 81, 2340–2361 *(unverified)*.
- Kurtz, T. G. (1970), Solutions of ordinary differential equations as limits of pure jump Markov processes, *J. Appl. Probab.* 7, 49–58 *(unverified)*.
- McFadden, D. (1974), Conditional logit analysis of qualitative choice behavior, in P. Zarembka (ed.), *Frontiers in Econometrics*, Academic Press, 105–142 *(unverified)*.
- Gouriéroux, C., Monfort, A. and Renault, E. (1993), Indirect inference, *J. Appl. Econometrics* 8, S85–S118 *(unverified)*.
- Smith, A. A. (1993), Estimating nonlinear time-series models using simulated vector autoregressions, *J. Appl. Econometrics* 8, S63–S84 *(unverified)*.

### 13.4 Infrastructure, Urban and Heritage Planning

- Madanat, S. (1993), Optimal infrastructure management decisions under uncertainty, *Transp. Res. C* 1, 77–88, doi:10.1016/0968-090X(93)90021-7.
- Madanat, S. (1993), Incorporating inspection decisions in pavement management, *Transp. Res. B* 27, 425–438, doi:10.1016/0191-2615(93)90015-3.
- Madanat, S. and Ben-Akiva, M. (1994), Optimal inspection and repair policies for infrastructure facilities, *Transp. Sci.* 28, 55–62, doi:10.1287/trsc.28.1.55.
- Ellis, H., Jiang, M. and Corotis, R. B. (1995), Inspection, maintenance, and repair with partial observability, *J. Infrastruct. Syst.* 1, 92–99, doi:10.1061/(ASCE)1076-0342(1995)1:2(92) (DOI read from the publisher's URL).
- Smilowitz, K. and Madanat, S. (2000), Optimal inspection and maintenance policies for infrastructure networks, *CACIE* 15, 5–13, doi:10.1111/0885-9507.00166.
- Papakonstantinou, K. G. and Shinozuka, M. (2014), Planning structural inspection and maintenance policies via dynamic programming and Markov processes, *RESS* 130: Part I, 202–213, doi:10.1016/j.ress.2014.04.005; Part II, 214–224, doi:10.1016/j.ress.2014.04.006.
- Papakonstantinou, K. G., Andriotis, C. P. and Shinozuka, M. (2018), POMDP and MOMDP solutions for structural life-cycle cost minimization under partial and mixed observability, *Struct. Infrastruct. Eng.* 14, 869–882, doi:10.1080/15732479.2018.1439973.
- Andriotis, C. P. and Papakonstantinou, K. G. (2019), Managing engineering systems with large state and action spaces through deep reinforcement learning, *RESS* 191, 106483, doi:10.1016/j.ress.2019.04.036.
- Andriotis, C. P. and Papakonstantinou, K. G. (2021), Deep reinforcement learning driven inspection and maintenance planning under incomplete information and constraints, *RESS* 212, 107551, doi:10.1016/j.ress.2021.107551.
- Morato, P. G., Papakonstantinou, K. G., Andriotis, C. P., Nielsen, J. S. and Rigo, P. (2022), Optimal inspection and maintenance planning for deteriorating structural components through dynamic Bayesian networks and Markov decision processes, *Struct. Saf.* 94, 102140, doi:10.1016/j.strusafe.2021.102140.
- Memarzadeh, M. and Pozzi, M. (2016), Value of information in sequential decision making: component inspection, permanent monitoring and system-level scheduling, *RESS* 154, 137–151, doi:10.1016/j.ress.2016.05.014.
- Haasnoot, M., Kwakkel, J. H., Walker, W. E. and ter Maat, J. (2013), Dynamic adaptive policy pathways: a method for crafting robust decisions for a deeply uncertain world, *Glob. Environ. Change* 23, 485–498, doi:10.1016/j.gloenvcha.2012.12.006.
