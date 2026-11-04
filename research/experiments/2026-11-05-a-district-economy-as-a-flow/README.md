# A District's Economy as a Flow

*Status: draft (2026-11-05). Research.* Study (2) of the shared problem [A District's Economy](../2026-11-05-a-district-economy-the-problem/README.md): the people of a few blocks of Melbourne's Hoddle Grid as a conserved flow on their streets, with its circulation, and a search that uses the model's structure, on the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)).

This is a study of methods on a real place, not a planning proposal. The shared problem's economics is a sketch whose parameters are assumptions, this study's model of the people is its own and provisional, and nothing here says what any frontage of Melbourne should hold.

## 1. The Question and Why

The user asked for a plan of a district's economy that uses a fluid model of its people, with their circulation, and an optimisation that uses the structured model: in space, in time and in frequency, with the people's flow and the arrangement of the buildings and their relations. The shared problem fixes everything but the people: the district, its 61 ground-floor frontages and eleven candidate links, the decisions (each frontage's kind and hours, each link made or not) and the value, the week's gross operating surplus. A study supplies the people passing each frontage at each hour of the mean week, P_f(h; plan). Study (1) reads them from the counters alone; study (3) from an action; this one from a fluid.

Two questions follow:

1. **What does a fluid model of the people say?** Where they are driven, where they eddy round blocks, how a passage changes their circulation, and what their week looks like in frequency.
2. **What does the model's structure buy a search?** Study (1)'s searches treat the value as a black box. Here the flow is linear in its sources and the value is a sum over the week's hours, and both can be used. The comparison is at equal budgets, with study (1)'s algorithms, seeds and budgets.

## 2. The Model of the People

The model is the study's own (`flow.hpp`), provisional, and written where the domains have no model of human flow: urban science's `movement` module is planned. Its numbers are calibrated in §3.

### 2.1 The Network and Its Widths

- **The network.** It is the shared problem's walking network as the urban-science domain builds it (`economy::Problem::with_links`), taken with every candidate link made: 369 junctions and 423 segments, eleven of them the links. A link not made is that network with its segment removed. The junctions a link would make on the streets it meets stay there as points along them.
- **The largest component.** 404 of the 412 segments as built form one component. The other five components are pieces of lanes and arcades the city model leaves unjoined; they carry nobody.
- **Widths.** Each segment has a walkable width, by its class in Hoddle's grid of 1837: streets of 99 ft (30.2 m) and little streets of 33 ft (10.1 m). The widths are assumptions:

| Class | Walkable width (m) |
|---|---|
| Swanston Street (closed to private cars since 1992) | 12 |
| the other main streets (two footpaths of 4.5 m) | 9 |
| little streets, other streets (two footpaths of 2 m) | 4 |
| the Bourke Street Mall | 12 |
| other pedestrian ways and arcades | 5 |
| laneways (the city model's default for a service way) | 3.5 |
| a lane carried through; a passage | 3.5; 5 |

### 2.2 The Walkers as Streams

People are many streams, each a trip from an origin junction to a destination junction. Walkers who share a destination share a potential: they are a class in the sense of Hughes (2002, *Transp. Res. B* 36, 507). Each stream is conserved, so that at every junction but its two ends what flows in flows out.

- **The choice of way.** A stream's walkers choose their way by a logit over every walk of the network (Bell 1995, *Transp. Res. B* 29, 287; Akamatsu 1996, *Transp. Res. B* 30, 369). A walk's states are the directed segments, so that a walker turns back only at a dead end: the link-based form of Fosgerau, Frejinger and Karlström (2013, *Transp. Res. B* 56, 70). Entering a directed segment e has the weight

  ω_e = exp(−θ ℓ_e) (w_e / w_ref)^(ℓ_e / ℓ_parent),

  so that a walk's log-weight is −θ times its length plus the log of the widths it takes. Between two parallel segments of the network as built, the walkers divide in proportion to width, as a fluid's flux divides by conductance. ℓ_parent is the length of the segment as built that a piece belongs to, so that a segment cut by a link's junction keeps its weight. w_ref = 12 m.
- **The walk sums.** With M the transitions between directed segments and Z = (I − M)⁻¹ = Σ_k M^k, the sum over every walk converges when M's spectral radius is below one. At the chosen θ = 0.02 per metre it is 0.821, with every link made too. Somewhere between θ = 0.005 and 0.01 per metre the sums start to diverge: the radius is 1.111 at 0.005 and 0.983 at 0.01.
- **A trip's flow.** A trip from o that ends on first entering d visits the directed segment c, on average, π_od(c) σ_d(c) / W_od times, with π_od = s_oᵀ Z⁽ᵈ⁾, σ_d = Z⁽ᵈ⁾ 1_{A_d} and W_od = s_oᵀ Z⁽ᵈ⁾ 1_{A_d}. Here s_o weights the segments leaving o, A_d is the set of segments entering d, and Z⁽ᵈ⁾ counts the walks with A_d absorbing. By Woodbury, Z⁽ᵈ⁾ = Z − Z[:, A_d] Z[A_d, A_d]⁻¹ (Z − I)[A_d, :], so that σ_d = Z[:, A_d] Z[A_d, A_d]⁻¹ 1. The flow is linear in the trips: streams superpose.

*Rationale.* A Darcy fluid in a network of pipes is driven by Kirchhoff's law through the conductances, and by Thomson's principle it dissipates least for its sources. Its circulation is therefore exactly zero (checked: shares 1.4 × 10⁻²⁷). Walkers are not such a fluid: they take near-shortest ways. The logit keeps what makes a fluid tractable (conservation and linearity) and adds the walkers' choice of way. The departure of their flow from a Darcy fluid's is then the circulation of §2.6.

### 2.3 The Sources and Sinks

Each hour of the mean week the trips are a sum of patterns. Each pattern is a fixed spatial law of origins and destinations times an amplitude x(h) ≥ 0, which is fitted (§3).

- **Entrances.** These are the network's 24 dead ends on the edge of its box, where its streets leave it. Each has a gross flow g_i(h), both ways, from the nearest of the City's counters outside the box on the same street, within 400 m (8 entrances). Failing that, a street or mall takes a counter within 120 m, scaled by width (6). Failing that, the median of the counters on streets, scaled by width (10). Counters at crossings, which count people crossing a street, are left out.
- **Buildings.** The 61 parcels and the 349 buildings of the context standing in the box carry their floor (footprint × storeys) on the two ends of their nearest segment. Work (2.68 M m² of offices and shops) is kept apart from homes (0.23 M m² of dwellings).
- **Attractors.** These are CLUE's 2024 establishments of the trading kinds in the box outside the 61 parcels: 485 retail, 399 food and drink, 636 services and 61 culture. Each is open by its kind's typical hours, with weight w_k (retail, food and culture 1, services ½; an assumption). The plan's frontages are attractors too, open by the plan's hours.
- **The patterns.**
  - *arrive* (work; home): from entrance i to the buildings, g_i(h) × the buildings' shares of floor;
  - *leave*: the reverse;
  - *through*: from entrance i to entrance j, g_i(h) g_j(h) / Σ g;
  - *local*: from building a to building b, in proportion to m_a m_b exp(−d_ab / 200 m) (λ = 200 m, an assumption);
  - *visit*: from the entrances and the buildings (half and half) to each open attractor, and back.
- **How the plan enters.** An open frontage of kind k draws w_k x_V(h)/W visits an hour, where W = 1,301.5 is the weight of every attractor as built. The passers at a frontage are therefore

  P_f(h) = B_f(h; links) + c(h) Σ_g w_{k_g} [g open at h] R_fg(links),

  affine in the plan's openings. B and R depend on the links through Z.

### 2.4 What a Frontage and a Counter See

A frontage and a counter each stand on one footpath. They see a share σ of their segment's walkers: ½ on streets and malls, which have two sides, and 1 on laneways (an assumption). A frontage therefore sees what a counter in its place would read.

### 2.5 The Links, by Woodbury

A link made adds the transitions into its two directed segments: two rank-one columns. If it carries a dead end on, it also removes the turning back there (one more). Z for any set of links is therefore Z0 + Z0 U (I − Vᵀ Z0 U)⁻¹ Vᵀ Z0, from one inverse of order 846 and a solve of order at most 33. It equals the direct inverse to 8.4 × 10⁻¹⁶ of the largest entry. All 2,048 sets of links are computed once, in 108 s on four threads, with the problem's own walking distances for each.

### 2.6 The Circulation

The net flow F = f(u → v) − f(v → u) of all the streams lives on the network's complex. Its junctions are the 0-cells, its segments the 1-cells, and the blocks physis traces (`plane_faces`) the 2-cells: 37 blocks as built, and 12 bounded faces with dead ends that stay holes, as the domain leaves them. The economics domain's `circulation` (on physis's flow with conductances c_e = w_e / ℓ_e) splits it into three flows, orthogonal in the energy Σ F²/c:

- the **driven** flow up an attraction α, which carries all of F's divergence. It is the Kirchhoff flow of the same sources (checked: relative difference 7 × 10⁻¹⁴), so it is what a Darcy fluid would carry;
- the **circulation round the blocks** that are 2-cells: the eddies;
- the **circulation round the holes**: the blocks with dead-end lanes, which the domain leaves open. Physically these are eddies round blocks too; the decomposition counts them as harmonic.

The spending is drawn as flux through the frontages: the week's spending at each frontage, entering it from its street.

### 2.7 Time and Frequency

The week is the circle of 168 civil hours, solved quasi-statically hour by hour. A walk across the district takes minutes, so at the daily line ω τ ≈ 2π × 10 min / 24 h ≈ 0.04 rad: the phase lag a dynamical model would add is about 2.5°. The model is linear in its sources, so each line of the week, X̂(m) = (1/168) Σ_h x_h e^{−2πimh/168}, is solved as a phasor: the sum of the sources' lines times their responses. This is checked equal to the lines of the hourly passers (relative 4.6 × 10⁻¹⁵). The lines reported are W1 (m = 1, the week), D1 and D2 (m = 7 and 14, one and two cycles a day), and D1 ± W1 (m = 6, 8, the days' modulation by the week).

## 3. Calibration

The amplitudes are fitted hour by hour by tyche's non-negative least squares to the City's counters inside the box: 16 of them, with the crossing at Degraves Street left out.

- **Each way.** 15 of the 16 counters record each way (`directions` in the data library's `street-counts`), placed on their segment's two directions by the compass. Their mean count each way at each hour of the week, over the reference year with Victoria's public holidays left out, is read by the urban-science domain's current meter (`urban::flow::current_meter`). The fit is to each way.

  *Rationale.* A trip and its reverse cross the same segments. Gross counts alone therefore cannot tell arrivals from departures. Fitted to gross counts only, the amplitudes of arrive and leave were zero at every hour, and the net flow vanished.
- **Handedness.** A counter on one footpath sees (1 + ρ_c)/2 of one way's walkers and (1 − ρ_c)/2 of the other's, since walkers keep to a side (Street Oceanography's handedness, `urban::flow::handedness`). ρ_c is fitted for each counter on a street by least squares over the week, alternating with the amplitudes in three rounds. The fitted values lie between −0.24 and +0.14.
- **Visits.** The counters cannot tell visits from the other trips: fitted freely, the visits' amplitude is zero at every hour. Their amplitude is therefore set each hour so that visits are ν = 0.1 of the hour's trips (an assumption), and the other modes are fitted to what is left. This costs the fit a little (§6.1).
- **θ and the modes.** θ was chosen over the grid {0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08} per metre by leave-one-counter-out R² of the gross. The homes' modes were then left out, because the fit is better without them.
- **The boundary.** As a check, the modelled gross on each entrance's segment is compared with the counter that feeds it.

## 4. The Structured Search

The value of a plan is a sum over the week's hours of the margin captured in each hour, less costs that belong each to one frontage or to one link:

V(plan) = Σ_h M_h(open set at h, kinds, P(h)) − Σ_f cost_f(k_f, hours_f) − Σ_links cost.

M_h depends only on the frontages open in hour h, their kinds and the passers then. The flow makes the passers affine in the openings, so a frontage's opening at h changes M_h and nothing else. `structured.hpp` uses this as follows.

- **The value hour by hour.** It equals the problem's `evaluate` to 1.0 × 10⁻¹⁴. Work is counted in margins of an hour: one evaluation of a plan is 168 of them.
- **A frontage's exact best response.** Hold every other frontage fixed. Then the value of frontage f with kind k and weekday hours [s, e), weekend hours [s′, e′), is Σ_{t∈[s,e)} Σ_{weekdays} gain_k(24d + t) + Σ_{t∈[s′,e′)} Σ_{weekend} gain_k(24d + t) − cost_f, where gain_k(h) = M_h(f open as k) − M_h(f closed). This is exact, not a linearisation. The best of the day's 300 intervals is found by enumeration, once for each day-type and each kind. A best response therefore costs 6 × 168 margins, six evaluations, and the value changes by exactly what it predicts (checked to 1.2 × 10⁻¹⁴).
- **BR.** Block-coordinate ascent: the frontages are swept in a random order (seeded), each taking its best response, until no frontage alone can do better. Then comes a descent over the links. Each link is flipped in a random order, the frontages are swept again from where they were, and the flip is kept if the value rises. This repeats until no flip helps or the budget is spent.
- **BR + (1+1) EA.** The same inner solve, with heuresis's (1+1) EA over the eleven links (`OnePlusLambda`, standard bit mutation conditioned on a change), warm-started from the incumbent's frontages. Its budget is counted in the inner solve's evaluations (`Evaluation::cost`).
- **GA, marginal kernels.** Study (1)'s GA with biased kernels. Each frontage's marginal value is computed once, at the plan as built: the best value of each of its six kinds with its best hours, 366 evaluations charged to the budget. The sampler takes, for each frontage with probability ½, a kind drawn with probability ∝ exp(V_f(k)/τ) (τ = A$1,000 a week, an assumption) and that kind's best hours. Otherwise it takes the problem's uniform kernel. The mutation does the same for each frontage with probability 1/61, after the problem's own mutation.

## 5. The Comparison

- **The black box, as study (1) runs it.** This is the same code, with its options, seeds 1–5 and budgets of 400 and 20,000 evaluations, on this study's flows:
  - the GA on the mixed encoding (μ = λ = 100, or 20 at 400; two elites; crossover 0.9; tournaments of two; the problem's kernels);
  - separable CMA-ES on the canonical relaxation (σ = 0.3, a random mean in [0.1, 0.9], restarted when it stops);
  - random search;
  - the GA started from the plan as built (its first generation the plan as built mutated five times over). Study (1) added this variant after this study's main run, so it was run afterwards on the same model (`--extra`), into `runs-extra.csv`.

  Each evaluation calls the problem's `evaluate`.
- **The structure.** BR, BR + (1+1) EA and GA with marginal kernels, at the same budgets and seeds. Work is counted in the same evaluations: margins of an hour ÷ 168, the table's cost included. BR and BR + (1+1) EA start from the plan as built.
- **The plans.** The best plan's value is read again by the problem's `evaluate`. The best plan and the best black-box plan are drawn beside the plan as built.

Study (1)'s MAP-Elites and NSGA-II search other objectives and are not run here.

*Caveat.* The values here and study (1)'s are on different flows: the plan as built is worth A$93,372 a week here and A$9,385 there. The comparison of searches holds within each study, not across them.

## 6. Results

### 6.1 Calibration

**The route choice.** Over the grid of θ, the leave-one-counter-out R² of the gross rises from −0.07 at 0.01 per metre to 0.51 at 0.015 and 0.02, then lies between 0.48 and 0.50 from 0.03 to 0.08. At 0.005 the walk sums diverge. With θ = 0.02 and the homes' modes left out (0.529 without them, 0.506 with), the fit is:

| Measure | Fitted | Each counter left out |
|---|---|---|
| R² of the gross, people an hour, 16 counters × 168 hours | 0.721 | 0.529 |
| R² in ln(1 + count) | 0.829 | 0.725 |
| R² of the net current, way₀ − way₁, 15 counters read each way | 0.549 | 0.088 |

- **By counter.** The median counter's R² is 0.50 (0.38 left out). The fit is best on Elizabeth Street at 197 (R² 0.98) and 124 (0.94), and on Swanston Street at the Town Hall (0.95). It is worst where the model's patterns do not reach: 231 Bourke Street on Royal Lane (0.03), Flinders Street at ACMI (0.04, over-predicted) and Flinders Lane at Degraves Street's south footpath (0.06, under-predicted). The model's day is smoother than the counters'. It keeps one crest at midday but not the commuters' sharp shoulders at 08:00 and 17:00 (figure 71).
- **What matters.** Without the through trips the left-out R² falls to 0.416, and without the local trips to 0.548. With every width alike (6 m) it is −0.68, so a segment's conducting as its width is essential to the fit. The set visits cost a little: without them it is 0.554.
- **The streams.** The through trips are positive at every hour (mean amplitude 0.34 of the entrances' flow) and the local trips in 142 hours. Arrivals at work are positive in 64 hours, departures in 5. The visits, set at ν = 0.1, are 0.089 of the trips over the week.
- **The boundary.** The modelled gross on the entrances' segments, over the gross their counters give, has a median of 0.95 (quartiles 0.85 to 1.10).
- **Handedness.** On every one of the 10 footpaths whose side the City names, the fitted ρ_c has the sign of left-hand traffic. Northbound walkers dominate on the west side, southbound on the east, eastbound on the north side and westbound on the south (|ρ| from 0.04 to 0.24). With no handedness this would happen with probability 2⁻¹⁰ ≈ 0.001 (a sign test, the footpaths taken as independent, which those of one street are not). It agrees with Street Oceanography's finding of 29 of 41 footpaths keeping left (P = 0.0058).

### 6.2 People Flow

- **The flow as built.** Summed over the segments of the network, the walkers number 268k an hour on a Wednesday at 08:00, 739k at noon, 677k at 17:00 and 807k at 14:00 on a Saturday (walker-segments an hour).
- **The busiest segments.** These are on Elizabeth Street, Collins Street and Little Collins Street, near 0.6 million walkers a week each.
- **The streams.** Over the week, the through trips carry 46% of the walking, the local trips 40%, the visits 9% and arrivals at work 5%. Departures carry 0.3% (figure 72).
- **What the plans change.** The passers barely move with the plan. Over the week the 61 frontages see 5.43 million passings as built, 5.35 million under the best plan (−1.5%) and 5.33 million under the best black-box plan. Their noon flows fall from 739k to 731k. A plan changes only where the visits go (9% of the trips), and the best plans open fewer frontages (35 against 49). The links reroute walkers: each made alone carries 8,400 to 88,000 walkers a week (§6.7).

### 6.3 Circulation

The net flow is the morning's tide. At Wednesday 08:00, Σ|net| is 29% of Σ gross, and its energy is 96.5% driven, 1.6% round the blocks and 2.0% round the holes (figure 73).

- **Where the eddies are.** The driven flow runs inward from the entrances to the floor, at most 1,991 people an hour on a segment. The largest eddies are 261 people an hour on a segment round the blocks and 163 round the holes. They lie along Swanston Street and Russell Street, where the walkers keep to the wide street more than a Darcy fluid would. The largest circulation round one block is −187 people an hour (clockwise) round a block of 1,442 m².
- **Through the week.** Weighted by energy over the 168 hours, the net flow is 96.6% driven, 1.5% round blocks and 1.9% round holes. The 38 hours that hold more than 1% of the largest hour's energy are the weekday mornings, the evening ebb and the weekend afternoons.
- **The quiet hours.** In the other hours the trips are nearly symmetric and the net flow nearly vanishes: at noon Σ|net| is below 0.1% of Σ gross, with an energy of 244 against 2.3 × 10⁸ at 08:00. What is left is mostly circulation. It comes from walks that pass back through their origin: an origin is not absorbing, so a trip and its reverse are not mirror images. This is a property of the logit over walks, not an eddy of the street, and the shares of those hours mean nothing.
- **How the links change it.** The best plan's three links make two more blocks and one more hole (39 and 13, against 37 and 12). At 08:00 the share round the holes falls from 1.97% to 1.64%, the share round the blocks rises from 1.56% to 1.64%, and the driven share from 96.47% to 96.72%. Carrying a lane through closes part of a hole into blocks, so part of the circulation that went round the hole now goes round blocks, as the topology requires. Made alone with the plan as built, no link moves the driven share at 17:00 by more than 0.3 points.
- **The spending as flux.** The best plan's frontages take in 1.88 million A$ a week, 1.70 million of it at 31 retail frontages. By the street it enters from, it is 857k on Swanston Street, 634k on Collins Street, 262k on Bourke Street and 132k on Russell Street. The lanes' and little streets' frontages are left vacant (figure 76).

*Remark.* A walkers' flow is close to a Darcy fluid's in its net. Its gross is not: the gross is mostly two-way traffic, which cancels in the net. The gross is what the frontages see. Circulation is a statement about the net flow, so it bears on the economy only through the tide, and here it is small.

### 6.4 Frequency

- **The tide on the streets.** As built, the daily line D1 of the gross flow has its crest between 14:50 and 16:05 on 90% of the segments (median 15:43). It is earlier on Flinders Lane and the southern streets and later on the lanes to the north (figure 74).
- **At the counters.** The model's D1 amplitudes follow the counters' (correlation 0.74, the model's a median 0.80 of the observed). Its crests do not: they lie between 15:11 and 15:56, where the counters' lie between 14:02 and 17:59 (correlation 0.24). The model has one tide for the whole district. Its patterns' amplitudes vary by the hour, but each pattern's spatial law is fixed, so it cannot give each street its own hour of high water.
- **How a plan changes the lines.** The passers' lines barely change (D1 amplitude of the passers at the frontages 35,302 as built, 34,761 under the best plan). The spending's do. Its D1 amplitude rises from 6,567 to 15,224 A$ an hour. The spread, the mean line's share of the power, rises from 0.419 to 0.489: from 70.3 to 82.1 effective hours a week. The daily lines' share falls from 0.570 to 0.490, and the week's other lines rise from 0.012 to 0.021. The best plan opens later into the evening (evening share 0.152 to 0.199; figure 77).

### 6.5 The Searches Compared

The best value found, in A$ a week (mean ± standard deviation over seeds 1–5; the plan as built is worth 93,372):

| Search | 400 evaluations | 20,000 evaluations | Seconds a run at 20,000 |
|---|---|---|---|
| random search | −88,981 ± 9,111 | −51,356 ± 6,378 | 2.8 (4 threads) |
| separable CMA-ES | −15,882 ± 28,565 | 393,140 ± 21,574 | 244 (4 threads) |
| GA | 389 ± 14,902 | 363,283 ± 19,136 | 2.7 (4 threads) |
| GA from the plan as built | 149,818 ± 11,274 | 350,567 ± 8,617 | 2.7 (4 threads) |
| GA, marginal kernels | 215,122 ± 29,657 | 448,738 ± 3,663 | 3.0 (4 threads) |
| **BR** | **448,443 ± 12,991** | **480,589 ± 13,741** (stopped at 11,744–13,579) | 1.5 (1 thread) |
| BR + (1+1) EA | 445,999 ± 10,854 | 473,530 ± 13,792 | 2.2 (1 thread) |

- **What the structure bought.** One sweep of exact best responses from the plan as built costs 367 evaluations. In every seed it reaches more (425,834 to 456,901) than any black-box search reaches with 20,000. The best black-box mean, CMA-ES's 393,140, is passed at 367 evaluations by BR and BR + EA in every seed. The GA with marginal kernels passes it at 1,866 to 2,666 evaluations. After 20,000 evaluations the best black-box mean is still 18% below BR's.
- **Convergence.** BR stops by itself, at a plan that no frontage alone and no single link can improve, after 11,744 to 13,579 evaluations and 23 link sets. Its five seeds reach four local optima: seeds 3 and 5 the same plan (487,996), the others 457,294 to 490,616. The (1+1) EA over the links, given the same inner solve, spends the whole budget and does a little worse (473,530). The descent over eleven bits is cheap, and the EA's random flips spend inner solves on sets of links already known to be worse.
- **The marginal kernels.** At 20,000 evaluations they lift study (1)'s GA by 24% and cut its spread across seeds by a factor of five (σ 3,663 against 19,136). At 400 they lift it from 389 to 215,122.
- **The start.** Starting the GA from the plan as built helps it early (149,818 at 400) but not late (350,567 at 20,000, below the uniform GA). BR's advantage is therefore not its start: one sweep from the same start is worth three times what the GA from it finds in 400 evaluations.
- **Time.** At equal evaluations BR is also the fastest by the clock. Separable CMA-ES's 244 s a run is heuresis's, not the model's (§12.2).

### 6.6 The Plans

| Plan | Value | Spending | Open | Kinds (vacant / retail / food / services / lobby / culture) | Links | Mix (nats) | Access (m) | Spread |
|---|---|---|---|---|---|---|---|---|
| as built | 93,372 | 735,678 | 49 | 12 / 14 / 18 / 13 / 4 / 0 | none | 0.883 | 149.5 | 0.419 |
| best black box (CMA-ES, seed 3) | 419,898 | 1,710,320 | 36 | 25 / 26 / 3 / 1 / 6 / 0 | 4 | 0.633 | 216.4 | 0.510 |
| best (BR, seed 1) | 490,616 | 1,884,982 | 35 | 26 / 31 / 4 / 0 / 0 / 0 | 3 | 0.081 | 227.7 | 0.489 |

- **The best plan.** It changes 45 of the 61 frontages. Retail takes 31 of them, most open on weekdays from 08:00 to 21:00 (14 frontages) or to 23:00 (9), and food and drink 4. The lanes' frontages and the little streets' are left vacant, and the services and lobbies go. It makes three links: Ramsden Place carried through to Collins Street, a lane to Regent Place, and Baptist Place to Collins Street.
- **Its accounts.** Its margin is 799,203 A$ a week, against 413,844 as built. Staffing falls from 283,640 to 252,070, since fewer frontages are open. Refits cost 31,277 and the links 3,894.
- **A monoculture, and why.** The value alone drives the district to retail. CLUE's co-location of retail with retail (CLQ 2.67, ln 0.98) is the strongest affinity after culture's. It lifts a retail frontage's chance to stop up to the clamp of 2 when retail stands within reach. Competition on one segment costs less than that catalysis gives, so the best plan puts 74 pairs of one kind on one segment (18 as built). Its mix falls to 0.08 nats and its walk to the nearest open daily need rises from 149 m to 228 m. A planner would not search value alone. The shared problem offers spread and mix as objectives for that reason.

### 6.7 What a Link Does

Each link made alone, with the plan as built, reroutes every stream at once (figure 78):

- **The walkers it carries.** These range from 8,373 a week (a lane to Athenaeum Place) to 88,277 (Ramsden Place carried through to Collins Street, 91 m) and 51,220 (Baptist Place to Collins Street, 36 m).
- **The frontages it passes by.** Every link moves passings between frontages, and ten of the eleven take more from the frontages' segments than they bring: Ramsden Place +107,555 and −131,025, Baptist Place +96,706 and −149,919. A link draws walkers off the frontages' streets into lanes that hold no frontage. The passage from Russell Place to Royal Lane is the exception (+33,573, −16,541).
- **The value with the plan as built.** Only one link pays for itself: the lane to Regent Place (+263 A$ a week, costing 721). The others lose between 91 and 3,749 a week.
- **From the best plan.** Flipping any one link and re-optimising the frontages lowers the value, by 431 to 9,186 A$ a week. The best plan is therefore a local optimum in its links as well as in its frontages. Its three links are worth their cost there and not as built. Removing Baptist Place carried through would cost 9,186 a week, against its own cost of 901: with retail clustered along Collins Street, the walkers it brings past those frontages are worth more than they were.

### 6.8 The Arrangement Used

| Plan | Value | Without catalysis (same flows) | Pairs of one kind on one segment | Catalysing pairs within reach | Party walls between trading frontages (one kind) |
|---|---|---|---|---|---|
| as built | 93,372 | 28,319 | 18 | 114 | 31 (10) |
| best black box | 419,898 | 104,157 | 50 | 90 | 21 (21) |
| best | 490,616 | 137,952 | 74 | 122 | 26 (25) |

- **Catalysis.** It is 70% of the value as built (65,053 of 93,372) and 72% of the best plan's (352,664 of 490,616). The best plan uses the arrangement more than the plan as built does. Party walls join frontages of one kind in 25 of its 26 trading pairs (10 of 31 as built). Its catalysing pairs within 127.6 m are mostly retail with retail (figure 79).
- **Competition.** The search accepts competition on a segment where catalysis outweighs it: 74 same-kind pairs share a segment, four times as many as built.

## 7. Checks

`check_flow` holds the model and the search to direct solves, conservation, closed forms and the problem's own value. All 11 pass:

| Check | Result |
|---|---|
| spectral radius of the transitions, every link made | 0.821 < 1 |
| Woodbury's update against the direct inverse, 6 sets of links | 8.4 × 10⁻¹⁶ of the largest entry |
| no walk sum negative | 26,752 entries exactly zero: segments no walk joins |
| a trip conserved at every junction, 8 trips | 1.6 × 10⁻¹⁵ |
| the absorbing destination by Woodbury against a direct solve | 1.3 × 10⁻¹⁵ |
| arrive, leave and visits against the sums of their trips | 4.2 × 10⁻¹⁷ |
| a Kirchhoff flow has no circulation (Thomson) | shares 1.4 × 10⁻²⁷ |
| the driven part of a walkers' flow is the Kirchhoff flow of its divergence | relative 7.1 × 10⁻¹⁴ |
| the passers' lines as phasors against the lines of the hourly passers | relative 4.6 × 10⁻¹⁵ |
| the value hour by hour against the problem's `evaluate`, 20 plans | relative 1.0 × 10⁻¹⁴ |
| a best response changes the value by exactly what it predicts, 8 plans | relative 1.2 × 10⁻¹⁴ |

## 8. Code Version and Running

| Library | Commit |
|---|---|
| The Science Infrastructures (root) | `512477b` |
| heuresis | `5730e70` |
| calx | `90aa710` |
| morphe | `14b3fb8` |
| physis | `e2b98d8` |
| tyche | `b4034b8` |
| nomos | `5382beb` |
| domains/urban-science, domains/architecture, domains/economics | as at `512477b` (unchanged since `0f4169c`) |
| the shared problem (`economy.hpp`, `results/problem.json`, `results/sensors.json`) | as written on 2026-11-05, not yet committed |

The code was built with g++ 11.5 and CMake 3.28 in WSL, in Release, without warnings from the study's code. The commands run from the repository's root:

```
cmake -S research/experiments/2026-11-05-a-district-economy-as-a-flow -B ~/.cache/science-build/exp-economy-flow -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-economy-flow -j4
P=research/experiments/2026-11-05-a-district-economy-the-problem/results
~/.cache/science-build/exp-economy-flow/check_flow $P/problem-inputs.json $P/problem.json $P/sensors.json
~/.cache/science-build/exp-economy-flow/as_a_flow $P/problem-inputs.json $P/problem.json $P/sensors.json research/experiments/2026-11-05-a-district-economy-as-a-flow/results --threads 4 --seeds 5
.venv/Scripts/python.exe research/experiments/2026-11-05-a-district-economy-as-a-flow/figures.py
```

The shared problem's `problem-inputs.json` and `district-geometry.json` are made by its `prepare.py` and `calibrate`. `check_flow` takes 10 s. `as_a_flow` takes 1,466 s on four threads of the development machine, shared with other work. Of that, the responses of all 2,048 sets of links take 108 s and separable CMA-ES 1,218 s. The same command with `--extra` runs only the GA from the plan as built (about 3 minutes), and `--calibrate` stops after the calibration.

## 9. Parameters

| Parameter | Value | Source |
|---|---|---|
| θ, the route choice's sensitivity to length | 0.02 per m | calibrated over a grid (§3) |
| walkable widths; w_ref | §2.1; 12 m | assumptions, by Hoddle's classes of street |
| σ, a footpath's share of its segment | ½ on streets and malls, 1 on laneways | assumption |
| an entrance's counter | same street within 400 m; else (streets, malls) within 120 m, scaled by width; else the median of the street counters × width / 9 m | assumption |
| λ, the local trips' distance | 200 m | assumption |
| w_k, the visits an open establishment draws | retail, food and drink, culture 1; services ½; vacant and lobby 0 | assumption |
| ν, the visits' share of the trips | 0.1 | assumption (the counters cannot tell it, §3) |
| the visits' origins | half the entrances by their flow, half the buildings by floor | assumption |
| ρ_c, a footpath's handedness | fitted, −0.24 to +0.14 | §3 |
| τ, the marginal kernels' temperature | A$1,000 a week | assumption |
| the economics | the shared problem's §4.3 | the problem |

## 10. Data and Provenance

| Data | Source | Licence | Retrieved |
|---|---|---|---|
| buildings, ways, the walking network | OpenStreetMap through the data library's city model, via the shared problem | ODbL 1.0, © OpenStreetMap contributors | 2026-10-31 |
| hourly pedestrian counts, each way; the counters' places | City of Melbourne, Pedestrian Counting System (the data library's `street-counts`, `street-sensors`) | the City's open data, CC BY 4.0 (the counts' metadata state none) | 2026-10-30 |
| establishments; floors | City of Melbourne, CLUE 2024 (via the shared problem's inputs) | CC BY 4.0 | 2026-11-04 |
| public holidays | Nager.Date (the data library's `street-holidays-victoria`) | Nager.Date's API terms | 2026-10-30 |

`results/flow-maps.json` holds the segments' geometry and is not committed (ODbL): `as_a_flow` makes it again. Everything else in `results/` is committed: `calibration.json`, `plans.json`, `summary.json`, `runs.csv`, `anytime.csv`, `runs-extra.csv`, `anytime-extra.csv`, `run.log` and `run-extra.log`.

## 11. Figures

`figures.py` draws ten figures from `results/` and the shared problem's geometry. `gallery.md` holds them with their captions.

| Figure | What it shows |
|---|---|
| `search-70-the-district-as-a-network` | the network by walkable width, floor, attractors, counters, entrances, the candidate links; the calibration over θ |
| `search-71-counters-and-flow` | the 16 counters' mean weeks against the fitted flow, and the flow with each counter left out |
| `search-72-people-on-the-streets` | the gross flow on every segment at four hours, and each stream's share of the walking |
| `search-73-circulation` | the morning's net flow: driven, round the blocks, round the holes; each part's energy through the week |
| `search-74-the-tide-of-the-district` | the daily line's amplitude and crest on every segment; against the counters |
| `search-75-what-the-structure-bought` | the searches at equal budgets: anytime curves and each seed's best |
| `search-76-the-plans-and-their-spending` | the plan as built, the best black-box plan and the best plan, with the spending as flux into the frontages |
| `search-77-the-week-of-the-plans` | the week's people, spending and frontages open; the spending's lines |
| `search-78-what-a-link-does` | each link's walkers, the passings it moves, and its value |
| `search-79-the-arrangement-used` | catalysis, competition and party walls in the best plan |

## 12. Caveats

### 12.1 What a Conserved Flow Leaves Out

- **People choose destinations as well as ways.** Here the trips' destinations are fixed patterns whose amplitudes are fitted. A plan draws visits in proportion to the weight of what it opens, from the entrances and the buildings in fixed shares. It changes no other trip: a shop draws no commuter off their way, and no visitor chooses one district over another.
- **People are not independent.** The streams superpose, so there is no crowding, no queue at a crossing and no slowing on a full footpath. Urban science's kinematic wave (`urban::flow::StreetDay`, on physis's `continuum`) holds those, and this study does not use them.
- **Trips have no chains.** A trip goes from one place to another and ends; a shopper's tour of several shops is several trips. Turning back is allowed only at dead ends.
- **A loop flow in quiet hours.** An origin is not absorbing, so a trip may pass back through its origin, and a trip and its reverse are not mirror images. Symmetric demand therefore leaves a small net flow that is pure circulation (§6.3). It is negligible beside the tide, but it means the circulation of the quiet hours is the model's, not the street's.
- **Gross counts do not fix direction.** The net flow rests on single footpaths' counts each way, read through a handedness that is fitted. Its leave-one-counter-out R² is 0.09. The fit puts most of the directional asymmetry into arrivals (positive in 64 hours) and almost none into departures (5 hours), so the buildings fill more over the week than they drain. The tide's flood is better fixed than its ebb.
- **Visits are set, not measured.** Their share ν = 0.1 is an assumption, and without it the fit is slightly better (§6.1). The plan's power to draw people rests on it.
- **The economics is the shared problem's.** Passers stop with fixed chances, and a visitor who came for a shop is counted as a passer who may stop (§4.1 of the problem).

### 12.2 Gaps and Bugs Found

- **urban-science `movement` (planned).** No domain holds a model of human flow on a network: trips, their choice of way, their calibration to counters. This study wrote one. The logit over walks among directed segments, its absorbing destinations by Woodbury and its rank-one link updates are general and belong in a foundation. Their natural place is beside physis's `statistical_mechanics` (first passage) or calx's linear algebra; their calibration to counters belongs in urban science.
- **heuresis `strategies`: separable CMA-ES is not O(n) per sample.** `CmaState::decompose` builds `inverse_root` as a dense n × n product even when separable (O(n³) per decomposition). `colour` and `whiten` multiply by the identity B in O(n²) per sample. At n = 621 a sample costs about 12 ms, against 0.14 ms for the GA's, both on four threads. Study (1)'s CMA-ES at 400 evaluations took 5 s for the same reason.
- **The shared problem.** A plan that maximises the value alone becomes a monoculture (§6.6). This follows from its same-kind affinities (retail with retail, CLQ 2.67, ln 0.98) and its clamp of catalysis at 2. It is not a bug, but a study that searches value alone should say so.
- **No other bug** was found in the domains or foundations used: the economics domain's `circulation`, physis's `flow_decomposition` and `plane_faces`, calx's dense solve, tyche's NNLS, urban science's `current_meter` and heuresis's GA, (1+1) EA, CMA-ES and random search all behaved as documented.
