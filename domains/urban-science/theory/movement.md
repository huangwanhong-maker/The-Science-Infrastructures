---
title: Walkers as Agents. The logit over walks as a Markov chain, its agents, their calibration and a street's regime as a decision
status: draft
date: 2026-11-05
module: domains/urban-science/movement (ADR 0071, ADR 0075); checks: movement/tests/check_movement.cpp (24)
---

# Walkers as Agents

*This note states the model of the module `movement`: walkers as agents on a district's walking network, the law by which they choose their way, how they walk, how they are calibrated to the counters, and a street's regime as a decision under partial observation. It proves what the checks hold (`movement/tests/check_movement.cpp`).*

The module's place is set by ADR 0071 and its design (development/ideas/2026-11-04-agents-and-decisions-for-urban-planning.md, §9):
- the walkers are agents of nomos's engine;
- their route choice is the logit over walks of the district-economy flow study (research/experiments/2026-11-05-a-district-economy-as-a-flow), so that the mean of the micro model is that macro model;
- the kinematic street of `flow` stays the macroscopic reference for crowds;
- a street's regime is a POMDP solved online by tyche's POMCP on the agents.

## 1. The Standing of the Model

The walkers are a *representation* in the sense of Street Oceanography (theory/street-oceanography.md, Withholding 2.1): a mechanism offered for exclusion, never asserted to be how people walk.

- **What is data.** The network (OpenStreetMap through the city model), the counts each way (the City of Melbourne's counters), the event permits.
- **What is fitted.** The trips' amplitudes, hour by hour, to the counters each way (§5).
- **What is assumed.**
  - The route law's sensitivity θ is the flow study's, calibrated over a grid.
  - The widths come by class of street.
  - The signals' timing and the speed–density law are declared.
  - The festival's walkable width is declared, and so are the weights of the regime decision's reward.

  Each is named where it enters.

## 2. The Network

**Definition 2.1 (walkway).** A walkway is a plane graph of junctions $J$ and segments $E$. Each segment $e$ has a length $\ell_e$, a walkable width $w_e$, a parent length $\ell^\circ_e \ge \ell_e$ (the length of the segment as built of which it is a piece), a rise $r_e$, and a kind (street, laneway, mall). Its *states* are the directed segments $S = \{2e, 2e+1\}$: $2e$ runs from the tail $u_e$ to the head $v_e$, $2e+1$ back.

**Definition 2.2 (crossings and phases).** A junction $x$ where at least three carriageways (segments of kind *street*) meet is a road intersection. Passing through $x$ from $a$ into $b$, a walker crosses $c_{ab}$ carriageways in expectation:
- $c_{ab} = 1$ if it goes straight on (its heading turns by less than $45^\circ$);
- $c_{ab} = \tfrac12$ if it turns, since on its side of the street it crosses one of the two streets half the time;
- $c_{ab} = 0$ at other junctions, or on turning back.

A signal's two phases run along the junction's axis and across it. A walker leaving along $b$ waits for the phase of $b$'s direction.

## 3. The Route Law

### 3.1 The Logit over Walks

Let $\theta > 0$, $w_\text{ref} > 0$, $\kappa \ge 0$ (the weight of width), $\gamma \ge 0$ (of grade), $\beta_x, \beta_s \ge 0$ (of crossings and of the expected wait at a signal of cycle $C$ and green share $g$). Write $\varphi$ for Tobler's factor of speed on a grade, $\varphi(s) = e^{-3.5(|s+0.05|-0.05)}$ (Tobler 1993).

**Definition 3.1 (transitions).** For a state $b$ open to entry,
$$\omega_b = e^{-\theta \ell_b} \Big(\frac{w_b}{w_\text{ref}}\Big)^{\kappa \ell_b/\ell^\circ_b} e^{-\theta\gamma \ell_b (1/\varphi(g_b) - 1)},$$
and $M_{ab} = \omega_b\, \tau_{ab}$ for $b$ leaving the junction $x = \text{head}(a)$, with
$$\tau_{ab} = \exp\!\big(-\beta_x c_{ab} - \beta_s c_{ab} [x \text{ signalled}]\, (1-g)^2 C/2\big),$$
excluding $b = \bar a$ (the reverse of $a$) unless no other state leaves $x$ open. A closed state has no transition into it, and its transitions out are kept.

*Remark.* With $\kappa = 1$ and $\gamma = \beta_x = \beta_s = 0$ this is the flow study's law (its README §2.2). The check holds the module's $M$ equal to the flow study's, entry by entry, to $2.8\times10^{-17}$.

*Rationale.* With $\kappa = 1$, two parallel segments of the network as built divide the walkers in proportion to their widths, as a fluid's flux divides by conductance. The factor $(1-g)^2C/2$ is the mean wait at a two-phase signal for a walker arriving at random: it arrives in the red with probability $1-g$ and then waits on average half the red, $(1-g)C/2$.

**Definition 3.2 (target; walk sums).** A target $A \subseteq S$ is a set of states, absorbing: those entering a junction $d$, or the states of a set of segments. Its *walk sum* $z_A \colon S \to [0, \infty)$ is the least non-negative solution of
$$z_A(a) = 1 \ (a \in A), \qquad z_A(a) = \sum_b M_{ab}\, z_A(b) \ (a \notin A).$$
If $\rho(M) < 1$ it is unique: $z_A = \sum_{k \ge 0} (M_{TT})^k M_{TA} \mathbf 1$ on the transient states $T = S \setminus A$. It is the weight of every walk from $a$ that ends on first entering $A$.

**Definition 3.3 (the walkers' chain).** A walker bound from junction $o$ to $A$ leaves $o$ along $b$ with probability $p^0(b) = \omega_b z_A(b)/W_{oA}$, where $W_{oA} = \sum_{b \text{ out of } o} \omega_b z_A(b)$. From $a \notin A$ it takes $b$ with probability
$$P_{ab} = M_{ab}\, z_A(b) / z_A(a),$$
and it stops on entering $A$. $P$ is stochastic on the states with $z_A > 0$, by Definition 3.2. It is Doob's $h$-transform of the walk by $h = z_A$. A walker draws $b$ by nomos's logit choice with utilities $\ln(M_{ab} z_A(b))$.

**Proposition 3.4 (the walkers' mean is the logit flow).** Let $\rho(M) < 1$ and $W_{oA} > 0$. The mean number of visits of the walkers' chain to a state $c$ is
$$f_c(o \to A) = \pi_{oA}(c)\, z_A(c) / W_{oA}, \qquad \pi_{oA}(c) = \omega_c [c \text{ out of } o] + \sum_{a \notin A} \pi_{oA}(a)\, M_{ac},$$
the flow study's visits. The flow study computes them as $\pi_{od} = s_o^\top Z^{(d)}$ and $\sigma_d = Z^{(d)} \mathbf 1_{A_d}$, with $Z^{(d)}$ obtained by Woodbury.

*Proof.* The chain's mean visits are $v = (p^0)^\top \sum_k P_{TT}^k$ on $T$. Write $D = \operatorname{diag}(z_A)$, so that $P_{TT} = D^{-1} M_{TT} D$ and $P_{TT}^k = D^{-1} M_{TT}^k D$. Then
$$v(c) = \sum_k \sum_a \frac{\omega_a z_A(a)}{W_{oA}} \frac{(M^k_{TT})_{ac}\, z_A(c)}{z_A(a)} = \frac{z_A(c)}{W_{oA}} \sum_k (s_o^\top M_{TT}^k)_c,$$
where $s_o(a) = \omega_a[a \text{ out of } o]$. The sum is $\pi_{oA}(c)$, the least solution of its equation. For $c \in A$, the visits are the entries from $T$ and from $o$, each with $z_A(c) = 1$. ∎

The check holds Proposition 3.4 on a grid district by the chain's own fundamental matrix, to $1.1\times10^{-15}$. The module's Gauss–Seidel solves match the flow study's dense formulas to $2.2\times10^{-16}$ for $z$ and $2.2\times10^{-15}$ for $f$, over 156 trips with a dead end among the destinations.

**Corollary 3.5 (micro and macro).** Trips of a stream depart as a Poisson process of rate $\lambda$ per hour, their pairs $(o, A)$ drawn from a law $p$. The walkers' entries into $c$, counted by their hour of departure, then have mean $\lambda \sum_{(o,A)} p(o,A) f_c(o \to A)$ in every hour. That is the logit flow. Their variance is $\lambda \sum p\, \mathbb E[V_c^2]$, a compound Poisson variance, with $\mathbb E[V_c^2] = \sum_i p^0_i N_{ic}(2N_{cc} - 1)$ from the chain's fundamental matrix $N = (I - P_{TT})^{-1}$.

*Proof.* Poisson thinning makes the pairs' departures independent Poisson processes. Each trip's visits are independent of the others'. The second moment of the visits of an absorbing chain is Kemeny and Snell's (1960, *Finite Markov Chains*, Theorem 3.3.5). ∎

This is exact in expectation, whatever the walkers' speeds, queues and signals: they change *when* a walker enters a state, never *which* states it enters, since the route law reads no density. Counted by the hour of entry, the walkers differ from the quasi-static flow by the time trips take, the dynamic lag the flow study estimated at 2.5° of the daily line.

### 3.2 Solving

$z_A$ and $\pi$ are solved by Gauss–Seidel sweeps. These converge when $\rho(M) < 1$, which is checked when the routes are made. They need no inverse, so they serve networks of any size. A district's targets are solved in parallel.

**Each state to its own tolerance.** The sweeps stop when no state changes by more than $10^{-13}$ of its own value. A walk sum falls as $e^{-\theta\ell}$, so over a district it spans many orders of magnitude: $W_{oA} = 2.7\times10^{-13}$ for one of the case study's trips, to a dead end at the box's edge. The walkers' choices read the ratios $z_A(b)/z_A(a)$, and the identity $\sum_{c\in A} f_c = 1$ depends on the smallest values as much as on the largest.

A first version stopped when the largest change fell below $10^{-14}$ of the largest value. On that trip it left the visits near the target converged only to $5\times10^{-4}$, and its expected flow differed from the flow study's dense inverse by $1.8\times10^{-4}$ of the largest flow. The per-state rule cannot stop early: a state not yet reached changes by all of itself when it is reached. The check holds a ladder 1.6 km long, with $W = 2.9\times10^{-19}$, to the dense formulas state by state to $1.1\times10^{-14}$, and its trip enters its target exactly once to $2\times10^{-16}$.

## 4. The Agents

### 4.1 State, Rules and Schedule

**Definition 4.1 (walker).** A walker carries:
- its identity (its order of departure);
- its stream and pair;
- its current state $c$ and its place $s \in [0, \ell_c]$ along it;
- its free speed $v_0 \sim \mathcal N(1.34, 0.26^2)$ m/s, clipped to $[0.5, 2.5]$ (Weidmann 1993, §4.1);
- its status: walking, asking to enter a segment, waiting at a signal, or dwelling.

**The step** (Δt = 2 s), synchronous (nomos's `Schedule::synchronous`). Every walker reads the configuration at the step's start and writes only its own part.

1. **Walking.** The speed is $v = \max(v_\text{creep},\ v_0\, \varphi(g_c)\, (1 - \min(1, \rho_\text{ahead}/\rho_j)))$, with $\rho_j = 5.4\,\text{m}^{-2}$ (Weidmann's jam density, as `flow`'s kinematic street) and $v_\text{creep} = 0.05$ m/s. Here $\rho_\text{ahead}$ is the density of the others, both ways, in the next $10$ m of its segment (an assumption). A walker stands, $v = 0$, when those going its way ahead are at $0.95 \rho_j$.

   Walkers dwelling stand aside, at a stall or a shop window. They count in the density that slows the others, but not among the walkers going their way, and not in the room a segment offers at its mouth. A crowd that dwells thus thins a street's stream without stopping it. Counted as walkers ahead, the festival's visitors locked the district: 715,299 person-hours on a day whose ordinary walkers spend 63,571.

   Walkers waiting at a segment's end, at red or for room, are counted in the densities as a standing queue behind the end, at $\rho_q = 2.5\,\text{m}^{-2}$ (about 0.4 m² each; an assumption after Fruin's levels of service for queues). Counted all at one point, they read as an infinite density. The first run of the week counted them so, and the district locked by mid-afternoon: 115,214 walkers were held, because the queues at red on every corner blocked the mouths of the streets leaving it.
2. **At the end of a segment.** The walker either:
   - arrives;
   - begins its dwell, if its stream dwells (an exponential time of given mean, standing at a uniform place on the target segment), after which it walks on to the segment's end and back to its origin by the same law; or
   - draws its next state by Definition 3.3.

   If the step crosses a carriageway at a signalled junction, the walker waits for its phase.
3. **Admission.** A segment admits its *supply* (Godunov's receiving flow; Daganzo 1994). That is $q_\text{max} = v_0 \rho_j w/4$ walkers a second while the walkers going its way in the first $10$ m past its mouth are no denser than $\rho_j/2$, else $v_0 \rho (1 - \rho/\rho_j)\, w$, and none at $0.95\rho_j$ (or when everyone there, dwellers too, stands at $\rho_j$). Those coming the other way, among them a queue at red on the corner, slow a walker but do not take its room. Requests are ranked by the instant asked, then by identity. Those not admitted wait at the mouth.
4. **Tallies of crowding.** Fruin's walkway levels are tallied for walkers walking or dwelling, by the density of the others within 5 m either way. A walker waiting in a queue stands at the queue's assumed density, not at one the model reads, so it is not tallied.

**Proposition 4.2 (reproducibility).** A run depends on the seed alone, not on the number of threads.

*Proof.* Phase 1 writes each walker from the configuration at the step's start and from its own stream, keyed by (seed, "walk", identity, step). Its result does not depend on the order of evaluation. Phase 3 is a function of the requests ranked by a total order. The tallies are accumulated in the population's order, which is the order of departure. ∎

Draws are keyed by identity, not by the walker's index in the population. Two runs that differ elsewhere (a regime) therefore give a walker the same speed and the same draws while its history agrees. These are common random numbers for comparing regimes. (nomos's `tick` keys draws by index; see §8.)

**Proposition 4.3 (conservation).** At every step, walkers departed = arrived + failed + present. At every junction $x$, summed over a run in which every walker has finished,
$$\sum_{c \to x} n_c - \text{ended}(x) + \text{begun}(x) = \sum_{x \to c} n_c,$$
where $n_c$ is the entries into $c$.

*Proof.* Each walker is created once, removed once, and present between. Its path is a sequence of states joined head to tail. Every passage through $x$ adds one entry into a state entering $x$ and one out of $x$, except at the walk's beginning and end. ∎

### 4.2 The Kinematics as a Kinematic Wave

**Remark 4.4.** The density each walker reads is a local Greenshields law, $q(\rho) = v_0 \rho (1 - \rho/\rho_j)$. Admission by supply is Godunov's flux for the Lighthill–Whitham–Richards model at a cell boundary. A queue behind a bottleneck therefore stands in space at an over-critical density, and its head walks off at the bottleneck's capacity.

A segment's single density, read by every walker on it, gave the wrong answer. A segment fed at capacity filled past the critical density and discharged only at $q(\rho) < q_\text{max}$, a capacity drop the kinematic wave does not have (2,591 an hour through a lane of capacity 6,512). Read locally, the lane holds its critical density and passes 6,409 (98.4%), while the queue stands upstream at up to 4.18 m⁻². The check holds it within 5% and never above.

**Proposition 4.5 (signal wait).** Walkers arriving at random at a signal whose phase is green for the share $g$ of a cycle $C$ wait $(1-g)^2C/2$ on average: 15 s for $C = 120$ s, $g = \tfrac12$.

The check measures 14.93 ± 0.30 s over 4,066 walkers.

## 5. Calibration to the Counters

A counter on one footpath reads a share $\sigma$ of its segment's walkers: ½ on a street of two footpaths, 1 on a laneway (the flow study's assumption). Of each way it reads $(1 \pm \rho)$ of that share, its handedness, fitted by the flow study between −0.24 and +0.14.

The streams are grouped into patterns, each a fixed law of pairs times an hourly amplitude. By Corollary 3.5 the expected readings are linear in the amplitudes:
$$y_{kw}(h) \approx \sum_g a_g(h) \sum_{s \in g} r_s(h)\, \sigma_k (1 + \rho_{kw})\, u_{skw},$$
where $u_{skw}$ is the stream's unit passes on the counter's state for way $w$, from the exact flow. The amplitudes $a_g(h) \ge 0$ are tyche's NNLS hour by hour (`calibration.hpp`).

A stream that dwells counts its walk back in the hour it departs. This is the quasi-static assumption, and a caller that needs the lag composes its own columns (the case study's festival does).

The check recovers three streams' amplitudes from exact counts to $4.4 \times 10^{-16}$. The case study reuses the flow study's own calibration (gross R² 0.7214, net R² 0.5494, reproduced exactly) and drives the walkers with it.

## 6. A Street's Regime as a Decision

### 6.1 The Regimes

A street $R$ (a run of segments with an entry junction) is held in one of three regimes:
- **open**;
- **one-way** (its states running back toward the entry are closed to everyone);
- **closed to walkers passing through** (its states closed to class 0, open to class 1: the street's own walkers, those whose trips start or end on it, a crowd bound for it among them; the caller puts their streams in class 1, and the case study splits each stream's pairs so).

A regime is a pair of route laws, one per class, applied from each walker's next junction on.

### 6.2 The POMDP

- **State:** the configuration of the walkers, and a hidden scale $\kappa$ of the crowd the day brings.
- **Actions:** the three regimes, each held for an epoch.
- **Observations:** the binned readings of the sensed counters in the epoch.
- **Reward:** $-(\text{person-hours present}) - \lambda\,(\text{person-hours at Fruin's level F, } \rho \ge 2\,\text{m}^{-2})$, with $\lambda$ an assumption. This is walkers' time and their safety.

It is solved online by POMCP (Silver and Veness 2010): the generative model is the walkers themselves, a planner's copy running an epoch under the action with its draws reseeded. The planner's copies run lightened, with running totals only and a coarser step. Their beliefs are particles refilled by rejection on the binned observations. Silver and Veness's convergence holds as the simulations grow; with few, the answer is an estimate, and the case study reports its value against the fixed regimes by simulation.

### 6.3 An Exactly Solvable Case

**Definition 6.1 (corridor).** Two ways lead from one place to another, one sensed and one not. The demand $k \in \{\text{normal}, \text{surge}\}$ is hidden and Markov (it stays with probability $p$). In an epoch:
- $k \to k'$;
- $\text{Poisson}(\lambda_{k'})$ groups arrive, each choosing the unsensed way with probability $p_a$ set by the regime (nomos's logit);
- each of the $n$ groups on it leaves with probability $q$;
- it holds at most $N$, and a group beyond is turned away.

The observation is the sensed way's count, capped at $O-1$. The reward is $R((k,n), a) = -\text{crush}\cdot n^2 - \text{detour}_a\, \lambda_k$.

**Proposition 6.2.** The corridor is a finite POMDP with
$$T((k,n), a, (k',n')) = P(k' \mid k) \sum_j \text{Bin}(j; n, 1-q)\, P\big(\min(N, j + U) = n'\big), \qquad U \sim \text{Poisson}(\lambda_{k'} p_a),$$
$$Z(o \mid (k', n'), a) = P\big(\min(O-1, S) = o\big), \qquad S \sim \text{Poisson}(\lambda_{k'}(1-p_a)).$$
Its agents are a generative model of the same law.

*Proof.* The arrivals choosing each way are independent Poisson numbers by the thinning theorem. So the sensed count is independent of the unsensed arrivals given $k'$, and $Z$ depends only on the state entered and the action, as a POMDP's must. The survivors are binomial, independent of the arrivals. ∎

The check holds the agents' one-step law of $(s', o)$ to the tables, every state and action, within 4.1 standardised deviations over 20,000 draws each. POMCP on the agents earns −14.32 ± 0.42 against tyche's HSVI interval, which closes at −14.2609.

## 7. A Run as a Record

A run writes its events in urban-science's record (events/record.hpp; ADR 0067), each with provenance *simulated*:
- **Couplings:** the street and junction spaces, a junction named by the streets that meet there.
- **Primitives:** the cross-sections.
- **Events:**
  - a street's regime over its span;
  - crowds: minutes in which ten walkers or more, on average, stood at Fruin's level E (1.08 m⁻²);
  - queues;
  - each counter's simulated hourly reading.

Every analysis of a record (its complex, its homology, its rates) then applies to a simulation as to an observation. The workbench draws a simulated record apart from an observed one, by its provenance.

## 8. Gaps

- **A general method held here.** The logit over walks among directed segments (route.hpp: the walk sums, the walkers' chain, the expected flow) names no city. It belongs in a foundation, beside physis's `first_passage` or nomos's `choice`. It is written here, as the flow study wrote it in its own code, because this work may not edit the foundations.
- **nomos `agents`.** `tick` keys each agent's draws by its index, not its identity. A population with births and deaths, which walkers are, needs identity-keyed draws (nomos's own open question 3). The module keys its draws by identity inside its rule and ignores the stream `tick` offers.
- **tyche's POMCP on a simulator.** Each simulation draws afresh, so a node's actions are compared through independent noise. For the walkers that noise is larger than the regimes' differences. Common random numbers across a node's actions (the scenarios of PEGASUS and DESPOT) would serve. Neither can POMCP be handed a belief from outside, so a lost belief costs a new planner.
- **The record's vocabulary.** events/record.hpp lists provenance as recorded, inferred, modelled or synthetic. ADR 0071 and this module write *simulated*, which the comment should add.
- **Dwellers in the stream.** A festival's visitors dwelling counted as walkers ahead and as room taken at the mouths, and they locked the festival day. They now stand aside.
- **A counter at a segment's very end.** Such a counter missed the walkers who started from that end, because a pass at the tail itself was never counted. Two of the district's counters lost one way entirely (Chinatown–Little Bourke Street, and Elizabeth Street–Lonsdale Street). A counting line is now kept half a metre inside its segment.
- **Data.**
  - The data library holds no record of the signals, of the footpaths' steepness, or of the pedestrian network's own links. The signals are assumed at every junction of three carriageways, with fixed offsets; the grade is zero. Manifests for the City's `footpath-steepness` and `pedestrian-network`, and Victoria's traffic signals, would replace these assumptions.
  - Counterflow is not slowed more than a crowd going one way. One-way running therefore buys nothing in this model beyond keeping walkers off a street.

## References

- Akamatsu, T. (1996), Cyclic flows, Markov process and stochastic traffic assignment, *Transp. Res. B* 30, 369–386.
- Bell, M. G. H. (1995), Alternatives to Dial's logit assignment algorithm, *Transp. Res. B* 29, 287–295.
- Daganzo, C. F. (1994), The cell transmission model, *Transp. Res. B* 28, 269–287.
- Fosgerau, M., Frejinger, E. and Karlström, A. (2013), A link based network route choice model with unrestricted choice set, *Transp. Res. B* 56, 70–80.
- Fruin, J. J. (1971), *Pedestrian Planning and Design*, Metropolitan Association of Urban Designers and Environmental Planners, New York.
- Hoogendoorn, S. P. and Bovy, P. H. L. (2004), Pedestrian route-choice and activity scheduling theory and models, *Transp. Res. B* 38, 169–190.
- Kemeny, J. G. and Snell, J. L. (1960), *Finite Markov Chains*, Van Nostrand.
- Lighthill, M. J. and Whitham, G. B. (1955), On kinematic waves II, *Proc. R. Soc. A* 229, 317–345.
- Silver, D. and Veness, J. (2010), Monte-Carlo planning in large POMDPs, *NIPS 23*, 2164–2172.
- Tobler, W. (1993), Three presentations on geographical analysis and modeling, NCGIA Technical Report 93-1.
- Weidmann, U. (1993), *Transporttechnik der Fussgänger*, ETH Zürich, IVT Schriftenreihe 90.
