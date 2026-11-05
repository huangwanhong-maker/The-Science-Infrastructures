---
title: Households among the Centre's Dwellings. A logit market, agents that form its rents, Schelling's limit, and a planner who cannot see demand
status: draft
date: 2026-11-05
module: domains/urban-science/housing
related: [nomos agents (schedules, draws by identity, logit choice, Schelling), tyche decision (MOMDPs, HSVI), tyche regression, ADR 0071 §8, ADR 0075]
---

# Households among the Centre's Dwellings

*This note states the model of `urban-science/housing` with its hypotheses and results, and says what the data can and cannot give it. The case study that runs it on the City of Melbourne is [research/experiments/2026-11-05-households-of-the-centre](../../../research/experiments/2026-11-05-households-of-the-centre/README.md).*

## 1. The Market

### 1.1 Definitions

**Definition 1.1 (segments and classes).** A *segment* j is a type of dwelling t(j) ∈ {house or townhouse, apartment, student apartment} (CLUE's classification) in an area a(j) (a CLUE small area). It has a stock S_j, a weekly rent r_j, a value δ_j and a distance d_j to the centre's anchor. A *class* c of households has a share w_c of the M households who consider the centre, a weekly income y_c and a taste γ_{c,t} for each type, −∞ where it may not live in the type (student apartments are for students; the model admits only one-person households to them).

**Definition 1.2 (utility).** With r̄ the rent of the outside option (living elsewhere in the metropolis),

$$V_{cj} = \delta_j + \gamma_{c,t(j)} - \alpha_c \ln(r_j/\bar r) - \theta d_j + \varepsilon_{cj}, \qquad \alpha_c = \alpha\,(\bar y / y_c)^{\pi},$$

ε independent Gumbel, the outside option's V = 0. Code: `Market::utility`; α = `alpha`, π = `income_power`, θ = `distance_cost`, r̄ = `reference_rent`.

So P_cj = e^{V_cj} / (1 + Σ_k e^{V_ck}) (McFadden 1974) and the demand is D_j = M Σ_c w_c P_cj (`demand`).

### 1.2 Results

**Proposition 1.3 (clearing).** If M exceeds the occupied stock of the segments each class may enter, there are unique rents with D_j(r) = (1 − v) S_j for every segment.

*Proof.* In z = ln r the demands are gross substitutes: ∂D_j/∂z_j = −Σ_c α_c M w_c P_cj (1 − P_cj) < 0 and ∂D_j/∂z_k = Σ_c α_c M w_c P_cj P_ck > 0. Demand is continuous, tends to 0 as z_j → ∞ and to the class-eligible M as z_j → −∞, so a solution exists; uniqueness under gross substitutability is Arrow, Block and Hurwicz's (1959). ∎

`clear` takes diagonal Newton steps on ln D_j − ln O_j with the own semi-elasticity ᾱ_j = Σ_c α_c D_cj (1 − P_cj)/D_j; with one class it is Berry's contraction read in rents.

**Corollary 1.4 (one class).** ln(O_j/(M − Σ_k O_k)) = δ_j − α ln(r_j/r̄) − θ d_j. The check holds the solver to it to 10⁻¹³.

**Proposition 1.5 (inversion; Berry 1994).** Given rents and occupied stock, the map δ ↦ δ + ln O − ln D(δ) is a contraction, and its fixed point is the δ at which those rents clear (`invert`).

**Proposition 1.6 (surplus; Small and Rosen 1981).** S = Σ_c M w_c ln(1 + Σ_j e^{V_cj}) / α_c, in units of ln r, satisfies Roy's identity ∂S/∂ ln r_j = −D_j (`surplus_log_rent`; checked by central differences to 3 × 10⁻¹⁰). In dollars, each class's term is multiplied by its mean rent paid (`surplus_dollars`, an approximation stated as such).

## 2. Agents

### 2.1 The Model

**Definition 2.1 (the city and its households).** Buildings, each in a zone (a CLUE block, in one area) with a capacity of each type; households, each with an income class, a tenure (owner or renter) and a size (one, or two or more), resident in one dwelling or searching. A segment's rent is one for all its dwellings.

**Definition 2.2 (a tick).** A quarter, by nomos's random sequential schedule (each household visited once, in an order drawn from (seed, "ordr", tick), each drawing from its own stream (seed, "agnt", household, tick)):

1. a resident leaves the centre at its tenure's exit hazard, whatever it offers;
2. a resident considers a move at its tenure's hazard, or when its block's share of its own class falls below f (Schelling's rule, if asked), or when its rent exceeds a share of its income (renters, if asked);
3. a searcher (a newcomer, a mover, a household evicted by a dwelling removed) chooses among the dwellings on offer, its own and the outside option, by the logit of u = δ_seg + γ − α_k ln(r/r̄) − θ d_b + η · like + ε;
4. each segment's log rent moves by κ ln((D_j + 1)/(V_j + 1)), capped at ±`cap`: D_j the searchers' *notional* demand, their probabilities over all the dwellings posted at the tick's start, and V_j those posted.

A dwelling vacated is offered from the next tick (`letting_lag`).

**Proposition 2.3 (options are exact aggregates).** In a zone, the dwellings of one type differ only by e^{−θ d_b}. The logit over the option (zone, type) with utility u + ln Σ_b n_b e^{−θ d_b}, its building then drawn with probability ∝ n_b e^{−θ d_b}, gives every dwelling the probability it has in the logit over dwellings one by one.

*Proof.* The logit's probability of a set is the sum of its members', Σ_b n_b e^{u − θ d_b} / Z = e^{u + ln W}/Z, and the conditional of b within the set is n_b e^{−θ d_b}/W. ∎

This makes a choice cost the number of zones on offer, not of buildings. The checks hold the probabilities equal to 10⁻¹⁵ and the draws within four standard errors.

**Remark 2.4 (two steps when no neighbour counts).** With η = 0 and the logit's scale 1, an option's utility is its segment's plus the log of its weight, so a searcher chooses a segment (or staying, or the outside option) by the logit of δ + γ − α ln(r/r̄) + ln W_seg and then an option in proportion to its weight: the same law, at the cost of the number of segments (`segment_choice`; checked to 10⁻¹⁵). Newcomers who look and do not enter then cost almost nothing, which a centre that is a small part of the households' choices requires (P₀ ≈ 1: in the case study, one newcomer in twenty settles). `compact` drops the households no longer active between ticks, renumbering the others (and so their streams) deterministically.

**Remark 2.5 (why the demand is notional).** If D_j counted only the probabilities of searchers who still found dwellings, it could not exceed what was posted: when searchers outnumber vacancies, the late ones find nothing and add nothing, and rents could not rise. A posted-vacancy market's excess demand is the demand for what was posted. A first version of the module counted the depleting list, and the rents of a full city stayed where the burn-in left them.

**Remark 2.6 (the natural vacancy is an outcome).** With a letting lag of one tick, the vacancy at a tick's end is the dwellings vacated in it plus those not let: turnover times the time to let, Rosen and Smith's (1983) natural vacancy, not a parameter.

### 2.2 Schelling's Limit

**Proposition 2.7.** Take a graph of buildings of one dwelling each, two classes, α = 0, every δ, γ, θ zero, a step utility η [like ≥ f], the best option taken with ties drawn uniformly over dwellings (in the options' order) and staying kept on a tie, no outside option, no hazard, and searching when like < f. Then from the same placement and seed the households' configuration equals that of `nomos::agents::Schelling` after every tick.

*Proof.* Both visit the agents in nomos's order. A household content at its cell draws nothing and stays, as Schelling's agent returns. A discontented one finds the maximal options: the vacant cells where it would be content (utility η), or, if there are none, no option better than staying (all 0), so it stays. Among the maxima, ordered by cell, it draws `below(count)` from its stream, as Schelling's agent draws among its `good` cells, which are the same cells in the same order. Its like share counts the occupied neighbours of the candidate, itself excluded, as Schelling's `like` does. ∎

The check holds it on a 20 × 20 torus for 30 ticks at f = 1/3 and 1/2, and the case study on a 40 × 40 torus for 40.

## 3. What the Aggregates Reveal

With one class, ln O_jq = ln M_q − ln(1 + Σ_k e^{V_kq}) + δ_j − α ln(r_jq/r̄_q) + ξ_jq, so

$$\ln r_{jq} = a_q + b_j - \tfrac{1}{\alpha}\ln S_{jq} + e_{jq}$$

with quarter and segment effects (`inverse_demand`, on tyche's least squares). The slope is −1/α; α(a_q − ln r̄_q) is the log of the demand index up to a constant, when the centre is a small part of the households' choices (P₀ ≈ 1).

**The identifying assumption.** A segment's stock is predetermined: towers are approved and built over years, so the stock is not correlated with the quarter's demand shock to that segment beyond what its effect holds. Developers build where demand grew, so this is a model, declared; the moving annual medians overlap, so the errors are clustered by segment.

## 4. The Planner

A MOMDP (Ong et al. 2010) on (x, y): x seen, the stock in units above a reference; y hidden, the demand in the same units, a Markov chain; a the units permitted, completed by the next period (x′ = min(x + a, X − 1)); the observation, the band of the clearing rent at (x′, y′) read with an error to each neighbouring band; the reward, the welfare at (x, y) less the cost of the units (`planner`). Beliefs live on {x} × Δ(Y); tyche's `to_pomdp` gives the POMDP whose HSVI interval bounds V*. The checks hold the MOMDP's beliefs to its POMDP's (10⁻¹⁶), a two-state case with demand known to its closed form V(0) = max(w₀/(1 − γ), w₀ − c + γ w₁/(1 − γ)), and HSVI's interval on a 2 × 2 case to exact value iteration.

## 5. What the Data Cannot Give

- **No dwelling's rent or move.** Rents are medians of new lettings by suburb group and type (Homes Victoria); moves are the Census's persons at a different address a year before, at one date.
- **No vacancy between censuses**, and the one Census (10 August 2021) was taken in a lockdown with the students away: its 27 % of private dwellings unoccupied is not a natural vacancy.
- **No rents of student apartments**, which let by the bed outside the bond data: the model starts them at the group's one-bedroom flats.
- **No incomes by tenure or by type of dwelling** at SA1 in the General Community Profile: classes, tenure and size are drawn independently (tenure by structure excepted).
- **Preferences are fitted to aggregates**, never observed: α to the rents and the stock, π to the Census's sorting, δ to the flows at 2002's rents, the owners' taste for houses to the Census's tenure. θ, the letting lag, the acceptance of one newcomer in five, the planner's costs and bands are assumptions.
- **The groups cut across the boundary.** Flemington-Kensington, South Yarra and Port Melbourne lie mostly outside the municipality, so the inverse demand uses only the six groups inside it.

## 6. References

- Arrow, K. J., Block, H. D. and Hurwicz, L. (1959), On the stability of the competitive equilibrium, II, *Econometrica* 27, 82–109.
- Arrow, K. J. and Hurwicz, L. (1958), On the stability of the competitive equilibrium, I, *Econometrica* 26, 522–552.
- Bayer, P., Ferreira, F. and McMillan, R. (2007), A unified framework for measuring preferences for schools and neighborhoods, *J. Polit. Econ.* 115, 588–638 *(unverified)*.
- Berry, S. T. (1994), Estimating discrete-choice models of product differentiation, *RAND J. Econ.* 25, 242–262.
- Duncan, O. D. and Duncan, B. (1955), A methodological analysis of segregation indexes, *Am. Sociol. Rev.* 20, 210–217.
- Huang, Q., Parker, D. C., Filatova, T. and Sun, S. (2014), A review of urban residential choice models using agent-based modeling, *Environ. Plan. B* 41, 661–689.
- McFadden, D. (1974), Conditional logit analysis of qualitative choice behavior, in P. Zarembka (ed.), *Frontiers in Econometrics*, 105–142.
- Ong, S. C. W., Png, S. W., Hsu, D. and Lee, W. S. (2010), Planning under uncertainty for robotic tasks with mixed observability, *IJRR* 29, 1053–1068.
- Rosen, K. T. and Smith, L. B. (1983), The price-adjustment process for rental housing and the natural vacancy rate, *Am. Econ. Rev.* 73, 779–786.
- Schelling, T. C. (1971), Dynamic models of segregation, *J. Math. Sociol.* 1, 143–186.
- Small, K. A. and Rosen, H. S. (1981), Applied welfare economics with discrete choice models, *Econometrica* 49, 105–130.
- Smith, T. and Simmons, R. (2004), Heuristic search value iteration for POMDPs, *UAI-04*, 520–527.
- Waddell, P. (2002), UrbanSim: modeling urban development for land use, transportation, and environmental planning, *JAPA* 68, 297–314.
