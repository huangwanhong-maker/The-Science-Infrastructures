# Street Oceanography Made Exact

*Status: draft, 2026-10-30.*

This note serves learning and research. It takes the user's working papers, W. Huang, *Street Oceanography: An Oceanographic Epistemology of Socio-Cultural Dynamical Systems* ("the paper", [SO]) and *Towards Street Oceanography: Ginza Street as a Case Study* ([G]), and states what the urban-science domain computes of them. Their quantities become definitions, their claims results with proofs or sources, and their candidate mechanisms representations whose standing is marked.

The paper proposes a field and withholds dynamical laws ([SO] §2). This note keeps that standing. It speaks in three registers, and every result the workbench shows carries its register:

- **Definition.** A quantity the paper carries (its Table 1, "carried" or "carried in part"), made computable from records.
- **Theorem.** A statement proved here or cited, which holds of any street whatever its dynamics.
- **Representation.** A model in the sense of [SO] §24, offered for sufficiency, for a consequence, or for exclusion, and never for prediction ([SO] Claim 24.1). The paper withholds that a crowd is a continuum (Withholding 2.1). A kinematic wave of pedestrians is therefore a representation, labelled as one wherever it is drawn.

The code:

- **The foundations** hold the general methods:
  - calx `signal`: tapers, multitaper and cross spectra, the bispectrum (calx decision 0013);
  - tyche `time_series`: spectral inference, harmonic analysis, Lomb–Scargle, red noise and surrogates; tyche `spatial`: co-location quotients and Clark–Evans (tyche decision 0005);
  - physis `continuum`: kinematic waves on networks and transport on graphs; physis `statistical_mechanics`: master equations, entropy production, first passage, mixing; physis `chemistry`: reaction–diffusion on graphs (physis decision 0016).
- **The domain** holds the street:
  - `street/`: the street system, its terrain, its hours and the visits of kinds of people;
  - `rhythm/`: calendars, records and their spectra;
  - `response/`: forcing and susceptibility;
  - `affinity/`: co-location, catalysis and spacing;
  - `flow/`: kinematic waves, budgets, transit and circulation;
  - `classification/`: the plane of segments.
- **The Urban Building Workbench** draws it, in its plug-in *Street Oceanography*.

§10 lists the amendments this note proposes to the paper.

## 1. The Object

### 1.1 The Street System

**Definition 1.1** ([SO] Definition 5.1). A *street system* is a quadruple (S, H, Φ, R):
- **The segment S.** A length of street [0, L] between intersections or marked changes of terrain, with arc length s in metres; or, for several segments, a graph of them.
- **The holdings H.** Each holding h has a position s_h, a set of uses u_h drawn from a finite set U, and a declared opening set O_h ⊂ W. Here W = ℝ/168ℤ is the circle of the hours of the civil week, Monday 00:00 at 0.
- **The flow Φ.** A random measure on time × the space Q of carried properties q (language, purpose, time available, prior acquaintance, …). A count is its zeroth moment: the number of units through a cordon in an interval, whatever they carry.
- **The record R.** The dated sequence of past states: tenancies, uses, rules and descriptions.

**Remark 1.2.** The paper adopts the segment as its unit provisionally ([SO] Claim 5.2). Here every quantity defined on holdings is also reported at the segment, and its within-segment variance is kept, so that the adequacy of the unit stays a question the data can answer.

### 1.2 Internal and External Terms

**Definition 1.3** ([SO] Claim 5.3). A *declaration* of a study divides its quantities into internal ones, whose evolution the study holds responsible to the others, and external ones (forcings), taken as given. The workbench records the declaration of every analysis. A forcing is either:
- periodic on a calendar (§3);
- dated (holidays, a rule's start and end, a lockdown); or
- irregular (weather).

## 2. Terrain

### 2.1 The Terrain Field and Its Contractions

**Definition 2.1** ([SO] eqs. 2–3). The *terrain* is a field T: S × W → ℝ^k. Its components are physical (width, cover), economic (price band), institutional (open or not, a rule in force) and symbolic (the use class and its standing). Each is measured in its own unit. An arrival of type a reads the terrain through a weighting w_a ∈ ℝ^k, giving the contraction

  τ_a(s, θ) = ⟨w_a, T(s, θ)⟩.

The workbench stores no contraction without its weighting, and reports the weighting with every scalar it draws.

**Theorem 2.2** (no canonical weighting; [SO] Claim 9.1 made exact). Let D = (ℝ_{>0})^k act on ℝ^k by changes of the components' units, (gT)_i = g_i T_i. For w ∈ ℝ^k write T ≼_w T′ when ⟨w, T⟩ ≤ ⟨w, T′⟩.

(a) No nonzero weighting is invariant: if ⟨w, gT⟩ = ⟨w, T⟩ for all g ∈ D and all T, then w = 0.

(b) The comparisons on which every non-negative weighting agrees form the dominance order: T ≼_w T′ for every w ∈ ℝ^k_{≥0} if and only if T_i ≤ T′_i for every i. This order is D-invariant.

(c) If the signs of the components' effects are not fixed (a price band attracts one arrival and repels another), the comparisons on which every w ∈ ℝ^k agrees are equalities only.

*Proof.*
- (a) Take T = e_i and g = 2 in the i-th place, 1 elsewhere. Then 2w_i = w_i, so w_i = 0.
- (b) If T′ − T ≥ 0 componentwise, then ⟨w, T′ − T⟩ ≥ 0 for every w ≥ 0. If instead (T′ − T)_i < 0 for some i, then w = e_i reverses the comparison. Positive scalings preserve the sign of each component, so the order is D-invariant.
- (c) For T ≠ T′ choose i with (T′ − T)_i ≠ 0, and w = ±e_i with either sign. ∎

**Corollary 2.3.** An index that ranks streets on one scale has chosen a weighting, and its ranking is a property of that choice ([SO] §9.1). The weighting-free content of a terrain survey is its dominance order. The workbench draws that order as a Hasse diagram, beside the rankings each named weighting produces. Part (c) is stronger than the paper's claim: where arrivals differ in the sign of a component's effect, nothing at all is comparable without a weighting.

### 2.2 The Length over Which Terrain Changes, and the Compression Ratio

The paper's compression ratio C = L_walk / L_grad ([SO] eq. 4) leaves L_grad undefined ("the distance over which the terrain changes appreciably"). For the categorical components, the uses, the following definition is adopted.

**Definition 2.4.** For holdings with uses, let P_same(r) be the probability that two holdings at separation r (in a bin about r) share a use. Let P_∞ = Σ_u p_u² be the probability for two drawn independently from the street's composition p. The *use correlation* is

  ρ(r) = (P_same(r) − P_∞) / (1 − P_∞).

The *terrain length* L_grad is ∫₀^{r₀} ρ(r) dr, where r₀ is the first zero of ρ (Corrsin's integral scale; Tennekes and Lumley 1972, §8.1). The compression ratio is C = L_walk / L_grad.

**Remark 2.5.**
- L_walk is the distance an arrival covers in a visit. No source used here measures it, and the workbench takes it as a declared parameter (400 m by default).
- The plane of §8 is invariant under a common rescaling of L_walk, so the classification does not depend on its value while it is common to the streets. It matters only where L_walk differs between streets, and then it must be measured.
- A finely interleaved street has ρ falling within a frontage or two, so L_grad is small and C large. A street of uniform use has ρ ≡ 0 because P_∞ = 1, so L_grad = 0. In that case C is undefined, and the workbench reports such a street as having a single regime.

### 2.3 The Departure from the Fully Mixed Arrangement

The paper offers, in place of an energy, "a measure of the departure of a segment's present arrangement from the arrangement it would occupy if fully mixed … computed from the composition alone", and leaves its construction open ([SO] §27.3, §36.2). It is constructed here.

**Definition 2.6.** Let n_{xu} count the holdings (or establishments) of use u at place x, p_{xu} = n_{xu}/N, with marginals p_x and p_u. The *departure from mixing* is

  M = Σ_{x,u} p_{xu} ln (p_{xu} / (p_x p_u)) = H(U) − H(U | X) = Σ_x p_x D(p(·|x) ‖ p(·)),

the mutual information between place and use (Theil and Finizza 1971; Mora and Ruiz-Castillo 2011). It vanishes exactly when every place holds the uses in the street's proportions, and it is at most H(U), which it reaches when every place holds one use only.

**Theorem 2.7** (decomposition across scales). If the places are grouped (holdings into segments, segments into districts), then exactly

  M = M_between + Σ_g p_g M_within(g),

where M_between is the departure computed on the groups' totals. *Proof.* This is the chain rule of mutual information, I(U; X) = I(U; G) + I(U; X | G), since G is a function of X (Cover and Thomas 2006, Theorem 2.5.2). ∎

The paper asks at what scale a street's heterogeneity lives ([SO] Claim 5.2, §19.4). Theorem 2.7 answers by an exact partition of M among the scales.

**Theorem 2.8** (exchange erodes the departure). Let a unit at place x move to place x′ with probability K(x′ | x), independently of its use: an exchange blind to kind. The new joint law p′(x′, u) = Σ_x p(x, u) K(x′ | x) has M(p′) ≤ M(p).

*Proof.* Since K does not depend on u, the variables U → X → X′ form a Markov chain. The data-processing inequality gives I(U; X′) ≤ I(U; X) (Cover and Thomas 2006, Theorem 2.8.1). ∎

**Corollary 2.9** (the form of the persistence criterion). The paper asks what allows the difference between two adjacent regimes to survive their exchange of arrivals ([SO] §9.2–9.3, §36.2). By Theorem 2.8, no process blind to kind can do so. A gradient of composition persists only where some process depends on kind and produces departure at least as fast as exchange destroys it: selection at an entry, retention that differs by use, a reaction among uses (§5), or the sorting of arrivals by the contraction of Definition 2.1. In continuous time:

  dM/dt = −(destruction by exchange) + (production by kind-dependent processes), the first term ≥ 0.

The paper calls the criterion "stated in form and not in content". With this corollary it has a content and a way to fail: a persistent gradient with no identifiable kind-dependent process would refute the account of §2.3.

**Representation 2.10** (pinning on a discrete street). Let each holding's state be bistable under a kind-dependent local process, and let holdings exchange with their neighbours at strength D. The front between two regimes is then pinned for D below a bound and moves above it (Keener 1987; Fáth 1998). For the piecewise-linear caricature du_n/dt = D(u_{n+1} − 2u_n + u_{n−1}) − u_n + H(u_n − a), the bound is in closed form (physis `chemistry`, `mckean_pinning`, decision 0016). A stationary front exists exactly when

  a ∈ (½(1 − 1/√(1 + 4D)), ½(1 + 1/√(1 + 4D))).

Its profile is u_n = μⁿ/(1 + μ) on one side, where μ = 2D/(2D + 1 + √(1 + 4D)). The interval shrinks to the single point a = ½ as D → ∞, the continuum limit. That a street is discrete, made of holdings, is what lets a front stay put over a range of imbalances: a continuum front moves unless the two regimes are exactly balanced. This is offered as a candidate mechanism, with its condition that can fail.

### 2.4 The Effective Terrain in Time

**Definition 2.11** ([SO] Claim 22.1, §14.1). The *open field* of holding h is 1_{O_h}(θ) on W. The street's *opening forcing* is

  F(θ) = Σ_h 1_{O_h}(θ),

the number of holdings open at hour θ of the week, and the *open fraction at s* is

  ω(s, θ) = Σ_h k(s − s_h) 1_{O_h}(θ) / Σ_h k(s − s_h),

with a kernel k of declared width. Declared hours are a declaration by the party best placed to make it ([SO] Claim 22.1). They are registered without instruments, and they are the first work of the campaign ([G] §13.2).

**Definition 2.12** (the spectrum of the opening forcing). F is a step function on the circle W. Its Fourier coefficients F̂(m), at frequencies m/168 h⁻¹, are its lines; each use class has its own. The lines at multiples of 7 (that is, of 1/24 h⁻¹) are the daily cycle. The others carry the week's modulation of the days: Saturday hours, Sunday closing, late-night Thursday.

**Definition 2.13** (the asymmetry of ascent and descent; [SO] §22.4). For each weekday, let o_h and c_h be the opening and closing times of the holdings open that day. The workbench reports three statistics:
- their dispersions σ_o and σ_c, with the asymmetry A = (σ_c − σ_o)/(σ_c + σ_o);
- Kendall's τ between the ranks of opening and of closing (τ > 0: the first to open close first; τ < 0: the first to open close last);
- the exact permutation null of τ, the holdings exchangeable.

The paper reads a regime switch from an ascent and a descent that retrace unequally. A ≠ 0, or τ significantly different from 0, is that inequality in the declared hours.

### 2.5 The Visits of Kinds of People

**Representation 2.14** ([SO] §8, §9; the user's question of which people go to which store when). Arrivals of type a are present on the street with a weekly profile N_a(θ) (office workers on weekdays at lunch, visitors at weekends). An arrival at s chooses among the holdings open at θ by the contraction of their uses with w_a, discounted by distance:

  V_a(h, θ) = N_a(θ) · 1_{O_h}(θ) exp(⟨w_a, e_{u_h}⟩ − β |s_a − s_h|) / Σ_{h′} (same).

This is a multinomial logit (McFadden 1974). V_a(·, θ) is the composition of the flow at each holding over the week. Its Fourier lines on W are the holding's clientele rhythm: a store visited by one type at one hour has lines spread across many harmonics, while one visited evenly has few. The representation shows how the observer-relative terrain ([SO] §7.3) separates kinds of people in time as well as in space. It predicts nothing, and its weights are declared, not fitted.

## 3. Rhythm: The Calendar Torus

### 3.1 Calendars as Angles

**Definition 3.1.** A *calendar* is an angle θ_c(t) ∈ ℝ/ℤ advancing at a frequency ω_c. The workbench uses three, taken over civil time in the zone of the street:
- **The civil week**, ω = 1/(7 d). Its harmonics include the civil day (7 cycles a week). The civil clock includes daylight saving time, the institutional reading of [SO] §7.4.
- **The synodic month**, ω = 1/29.530589 d. Its mean, with the mean new moon of 2000 January 6, 14:20 TT as epoch (Meeus 1998, *Astronomical Algorithms*, eq. 49.1 with k = 0).
- **The tropical year**, ω = 1/365.24219 d (Meeus 1998, ch. 27). The Gregorian year, 365.2425 d, differs from it by 2.3 × 10⁻⁹ d⁻¹, so for every record they are one constituent (Remark 3.5).

Dated forcings (holidays, rules) are not angles. They enter as indicators.

**Definition 3.2.** A forcing periodic on the calendars is F(t) = Φ(θ(t)) with Φ a function on the torus T^r and θ(t) = ωt + θ₀. Its *Fourier support* Σ ⊂ ℤ^r is the set of k with Φ̂(k) ≠ 0. Its lines lie at the lattice image {k·ω : k ∈ Σ}. Φ is *separable* when Φ(θ) = Σ_c φ_c(θ_c); then Σ lies on the axes.

**Remark 3.3** (the day and the week are commensurate). The day is exactly a seventh of the week, so the line traced by θ(t) on the (day, week) torus is closed. The torus of independent calendars is T³ (week, synodic month, year), not T⁴. "Sums and differences" of the daily and weekly lines are harmonics of the week, m/7 d⁻¹, whatever produced them. A test that reads them as interaction between the day and the week reads a property of the week alone.

The code names a week harmonic m = 7q + r (|r| ≤ 3, unique since 7 is odd) by days and weeks, D<q>±W<|r|>. Its order is |q| + |r|, the number of steps of a day or a week: W8 is D1+W1, of order 2.

### 3.2 Resolution and Aliasing

**Proposition 3.4** (the Rayleigh criterion; Godin 1972, *The Analysis of Tides*, §2.3; Foreman 1977). A record of length T separates two lines f₁, f₂ only when |f₁ − f₂| T ≥ R, with R = 1 by convention.

**Remark 3.5** (the record-length constraint, refined). [SO] §14.4 states that a component of period P is unavailable to a record shorter than P. That is the case f₂ = 0, separating a line from the mean. The binding constraint is the separation from the nearest line. Three consequences:
- The synodic month's sidebands of the daily cycle (1 ± 1/29.53 d⁻¹) need T ≥ 29.5 d.
- The annual sidebands of the weekly lines need T ≥ 1 yr.
- The tropical and the Gregorian year, 2.3 × 10⁻⁹ d⁻¹ apart, need T ≥ 1.2 × 10⁶ yr, so no record separates them.

The workbench computes, for a record and a list of constituents, which pairs it can separate.

**Proposition 3.6** (aliasing; [SO] §22.5). A line at f sampled at the rate f_s appears at |f − f_s · round(f/f_s)|. Hourly counts alias every component faster than 12 d⁻¹, so the separating observation of [SO] Table 5 (a finer interval) bears only on structure above that.

### 3.3 What a Response Can Contain

**Proposition 3.7** (the line support of responses). Let F = Φ(θ(t)) with Fourier support Σ, and let y be the response of a time-invariant system.
- (a) If the response is linear (y = h ∗ F), its lines lie in {k·ω : k ∈ Σ}.
- (b) If it is a Volterra series of order ≤ d, its lines lie in {k·ω : k ∈ Σ^{(≤d)}} ∪ {0}, where Σ^{(≤d)} is the set of sums of at most d elements of Σ ∪ (−Σ).
- (c) A line at a frequency outside {k·ω : k ∈ ⟨Σ⟩}, the subgroup Σ generates, is produced by no time-invariant response to F.

*Proof.* The Fourier transform carries a product of trigonometric polynomials to the convolution of their coefficients. A Volterra kernel of order j therefore maps lines k₁…k_j to k₁ + … + k_j, and a linear filter multiplies each line by its gain at that frequency. A time-invariant functional of finite memory is a limit of Volterra series (Boyd and Chua 1985, IEEE Trans. Circuits Syst. 32, 1150), and limits keep the support inside the closure of ⟨Σ⟩·ω. ∎

**Remark 3.8** (separability; an amendment to [SO] Discriminator 29.1). The paper's middle case reads lines "at sums, differences and multiples of the drivings" as disproportionate response. By Proposition 3.7(a), a *linear* response to a non-separable forcing already has them. Institutional forcings are non-separable by construction:
- hours that differ by weekday are non-separable on the week, though by Remark 3.3 that stays within the week's harmonics;
- hours that change with the season are non-separable on (week, year). An example is Ginza's vehicle-free arrangement, ending at 18:00 from April to September and at 17:00 from October to March ([G] Observation 8.1).

The discriminator's middle case therefore holds only relative to a forcing known to be separable, or after the linear response to the declared forcing has been removed. The workbench removes it: the declared hours give Φ, and the lines tested are those of the residual.

**Proposition 3.9** (the lattice is dense; the order of explanation). If ω has two or more rationally independent components, {k·ω : k ∈ ℤ^r} is dense in ℝ (Kronecker). Every frequency is then within any resolution δ = R/T of some combination of high enough order. "A component derivable from no driving" ([SO] Discriminator 29.1, third case) must therefore be qualified by an order d: no combination with |k|₁ ≤ d lies within δ. The fraction χ(d, T) of a band [0, f_max] within δ of such a combination is the probability that a random frequency would be "explained" by chance, the discriminator's false-explanation rate. The workbench computes χ exactly as a union of intervals. In the paper's constructed example (drivings of 7 and 29.5 d, 512 days hourly) χ grows with d, and the evidential weight of the third case falls with it. (This is the question of the idea note *The Tidal Spectrum as a Lattice*, here with the calendars as the lattice's generators.)

**Remark 3.9a** (the week is free). A forcing that is any function on the week's circle has every line m/168 h⁻¹ in its support. When the week is a declared driving, all its lines are therefore forced, the day's fourth and fifth harmonics as much as its first. The order that measures an explanation is then the *cross order* |k_M| + |k_Y|: the steps across the incommensurate calendars, with the week's own lines at cross order 0. The chance rate of the week's lines alone is already 2δ·7 per day of bandwidth. For a year of hourly counts read with a multitaper bandwidth of 0.01 d⁻¹, that is about 15% of the band. A street's "own" line has to fall in the rest.

**Remark 3.10** (phase coupling and its limit). Bicoherence detects quadratic phase coupling among components with random phases (Kim and Powers 1979). The lines of a deterministic periodic forcing keep their phases in every segment of the record, so three lines at f₁, f₂ and f₁ + f₂ have b² ≈ 1 whatever produced them. Bicoherence thus discriminates the response to *irregular* forcing (weather on a weekend afternoon), not to the calendar. A symmetric system has a vanishing bispectrum, so a nonzero bispectrum of the irregular part is the "asymmetry together with a nonlinearity" that [SO] §13.3 requires for converting zero-mean driving into standing activity.

### 3.4 The Standing Null of Integration

**Definition 3.11** ([SO] §28.2; Hasselmann 1976). The *integration null* is the AR(1) process with the record's lag-one autocorrelation and variance. Its density (2σ²/f_s)(1 − ρ²)/(1 + ρ² − 2ρ cos 2πf/f_s) is red. A line is reported only when it exceeds the null's level for its degrees of freedom (multitaper, χ²_ν). Claims of memory are tested against red-noise and phase-randomised surrogates (Theiler et al. 1992).

### 3.5 The Institutional Reading, Measured

**Representation 3.12** (daylight saving time as a natural experiment; [SO] §7.4, [G] Observation 8.1). At a transition the civil clock moves by an hour against the sun. Let φ_d be the UTC time of the daily maximum of activity, estimated by a harmonic fit (1 and 2 cycles a day) over the same weekdays in windows before and after the transition. The *civil share* is

  ψ = −Δφ / 1 h at the October transition, +Δφ / 1 h in April (Melbourne).

ψ = 1 means the street keeps the civil clock; ψ = 0 means it keeps the sun. The paper's fourth departure says the forcing reaches a street through an institutional reading, which predicts ψ ≈ 1 for streets driven by work and timetables. A street of leisure or of daylight could read the sun more. The test needs only hourly counts, which Melbourne has kept since 2009.

### 3.6 The Tide of a Street

The paper takes from oceanography its form, not its equations ([SO] §2). One of that form's oldest instruments applies to a street without borrowing any law: the tide gauge and its harmonic prediction.

**Definition 3.13** (the tide; `rhythm/tide.hpp`). Let y be a record of hourly counts, each placed at the middle of its hour in civil time τ, measured in hours from the week's origin. Its *tide* is the least-squares fit

  ŷ(τ) = Z₀ + Σ_c A_c cos(2π f_c τ − φ_c)

at the street constituents c that the record resolves (the week's harmonics, the year and its sidebands of the day; Rayleigh's criterion). It is fitted without a trend, as a port's tide is predicted about a fixed mean level (Foreman 1977; Pugh and Woodworth 2014, ch. 4). The *surge* is y − ŷ. It holds what the calendars do not drive: holidays, weather, events.

**Definition 3.14** (readings of a tide). Take the week's constituents alone (k_M = k_Y = 0), and the mean of the curves of Monday to Friday, and of Saturday and Sunday:
- *high water* and *low water*, their levels and civil hours. The hour of high water plays the part of a port's establishment;
- the *range*, high water − low water;
- the *form*, F = Σ_{q=1} A_c / Σ_{q=2} A_c over the constituents D<q>±W<r>. The lattice's day index q is the species of tidal analysis (Doodson's first index). Courtier's (1938) classes: F < 0.25 semidiurnal; 0.25 ≤ F < 1.5 mixed, mainly semidiurnal; 1.5 ≤ F ≤ 3 mixed, mainly diurnal; F > 3 diurnal.

*Remark 3.15.* The prediction is exact for a record lying on the lattice. The checks hold a constructed street to 10⁻⁸ people an hour, its high water to the brute force's minute, and its form to 10⁻¹³. A day's dip appears as the surge, less the share of it the least squares spread to the constituents: 2.3% of the dip, which reaches other days in proportion. Courtier's thresholds are the sea's, borrowed as a convention. On Melbourne's 100 outdoor sensors (October 2024 to October 2026), 72 streets are diurnal, 26 mixed and mainly diurnal, and 2 mixed and mainly semidiurnal. 50 reach their weekday high water between 16:30 and 18:30, and 41 between 11:30 and 14:30.

## 4. Forcing and Response

**Definition 4.1** ([SO] §13.1). The *response function* of a segment to irregular and dated forcings is estimated on the log record:

  ln(1 + y_t) = b(θ_week(t)) + Σ_j a_j cos/sin(2πj θ_year(t)) + Σ_f Σ_ℓ β_{fℓ} F_f(t − ℓ) + ε_t.

- b is a free value for each hour of the week.
- The forcings F_f are precipitation (mm/h, lags 0–3 h), the temperature's departure from its hour-of-year mean, and indicators of public holidays and of declared rules.
- The *susceptibility* to a forcing is Σ_ℓ β_{fℓ}, its proportional response.

Least squares with standard errors clustered by day (tyche `regression`), since the residuals are correlated within a day.

**Definition 4.2** (selectivity; [SO] §14.3, Figure 2). The paper's plane uses a "selectivity of response" it does not define. Group the fitted lines of the record into the drivings D: the day (lines at multiples of 1 d⁻¹), the week (the other multiples of 1/7 d⁻¹), the year and the synodic month. Let π_D be the share of the forced variance carried by D, Σ_D π_D = 1. The *selectivity* is

  S = (|D| Σ_D π_D² − 1) / (|D| − 1) ∈ [0, 1].

S = 0 when the drivings share the response equally, and S = 1 when one carries it all. It is the normalised inverse participation ratio. The paper's tension between amplification and breadth (§14.3) is S against the forced fraction of variance, which the workbench reports beside it.

**Definition 4.3** (the comparison across a repeating alteration; [SO] §26.3, [G] §11). Let a rule alternate a segment between configurations A and B on a schedule. Pair each A-day with the nearest B-day of the same season and weekday class. The effect is the mean of the paired differences of ln(1 + y). Its significance comes from randomisation inference, flipping signs within pairs. The workbench applies the design to the configurations a record supports, and states the anticipation that [SO] Claim 25.1 attaches to it.

## 5. Affinity and Reaction

**Definition 5.1** ([SO] §11.1). The *affinity* of use j for use i at separation r is estimated by the co-location quotient CLQ_{i→j} (Leslie and Kronenfeld 2011). Its significance comes from permuting uses, either freely or within strata of terrain, which is the null of independent response to terrain ([SO] §28.2). A_ij := ln CLQ_{i→j}. CLQ is asymmetric, as the paper requires of A ([SO] §11.1: "not symmetric").

**Representation 5.2** (catalysis in space and time: the marriage store and the ring store). A visitor arrives at a bridal store B at position s_B, with the need of a wedding dress. While at B, a second need, rings, arises at the rate κ. B is not consumed, which makes it a catalyst ([SO] §11.3, third mechanism). At a jeweller J at s_J the ring need is met at the rate μ.
- **Movement.** The visitor moves along the street as a random walk between frontages (rate D/Δs² to each neighbour) and leaves the street at the rate ℓ, or at its ends.
- **The process.** The visitor's state is (position, need) ∈ S × {dress, ring}, a Markov jump process with absorbing states *bought* and *left*.
- **The yield.** Y(d) is the probability of absorption at *bought* when d = |s_J − s_B|; physis's absorption probabilities give it exactly.
- **The temporal coupling.** The stores are open on O_B and O_J. The generator then changes hour by hour, and the yield averages over the hours of arrival at B. J closing before B cuts the pathway even when the two are adjacent.

Three consequences are computed exactly:
- the *coupling length* ξ_c, the e-folding distance of Y(d) beyond the stores' neighbourhood.
  - For a wandering visitor it is the decay of a killed walk's hitting probability: ξ_c = −Δs / ln λ, where λ + 1/λ = 2 + ℓΔs²/D. It tends to √(D/ℓ) on a fine street: about 42 m for a visitor wandering at D = 1 m²/s who leaves after half an hour on average.
  - A search speed v toward the ring store (the store seen, signposted or known: the symbolic terrain of §2.1) lengthens it toward v/ℓ. Legibility is part of affinity;
- the *coupling window*, Y against the offset of J's hours from B's;
- *creation without capture* ([SO] Claim 11.1): B's contribution to J's takings is Y times the value of a ring sale. An allocation of frontage by own takings does not register it, so removing B lowers J's takings by Y while costing B nothing it counted.

The pair of stores then has a function neither has alone, the ring sale. That function exists only within ξ_c and within the coupling window. This is the user's "marriage store and diamond store coupled in space–time", stated exactly.

**Representation 5.3** (spacing from local reinforcement and longer-range suppression; [SO] §11.3, fourth mechanism; Turing 1952). On the chain of holdings, an activator–inhibitor pair becomes Turing-unstable when d_v/d_u exceeds the critical ratio. Its fastest mode sets a spacing 2π/k. The prediction can fail: the observed spacing of clusters of a use is the dominant wavelength of its indicator along s, from a Lomb–Scargle periodogram of the uses' positions. The two lengths are compared, and the condition (suppression reaching further than reinforcement) is reported with them.

## 6. Flows, Retention and Budgets

**Theorem 6.1** (person-time is an exact budget; [SO] Claim 27.1). Let units i arrive at a_i and depart at d_i, and let N(t) be the number present. For every window [0, T],

  ∫₀^T N(t) dt = Σ_i |[a_i, d_i] ∩ [0, T]|.

*Proof.* N(t) = Σ_i 1[a_i ≤ t < d_i]; integrate over [0, T] and exchange the sum and the integral. ∎

Dividing by T gives Little's law L = λW in the limit (Little 1961; Kim and Whitt 2013 for finite windows), with boundary terms from the units present at 0 and at T. The paper's exact budget is this identity. From cordon counts, N(t) = N(0) + A(t) − D(t), and the budget closes to within the miscount of the cordons.

**Definition 6.2** (transit times; [SO] §10). Movement among the places of a segment, with exits as absorbing states, is a Markov jump process. The transit time is phase-type, P(τ > t) = α e^{Tt} 1 (Neuts 1981). Its tail, its moments and the probability of return to a place are computed exactly. The *duration threshold* of an exchange that needs τ* is P(τ ≥ τ*). No rate of arrival compensates if it vanishes ([SO] §10.3).

**Theorem 6.3** (circulation is the non-exact part of the rates). Let k_ij be the rates of a movement process on a graph whose edges carry rates both ways, and let a_ij = ln(k_ij / k_ji), a 1-cochain.
- (a) The process satisfies detailed balance if and only if a is exact (a = dV), that is, every cycle affinity A(C) = ∮_C a vanishes (Kolmogorov 1936).
- (b) Its entropy production σ = Σ_C J(C) A(C) ≥ 0 (Schnakenberg 1976) therefore depends on the rates only through the cycle affinities, which are the periods of a's non-exact part.
- (c) If the graph is the 1-skeleton of a 2-complex (streets bounding blocks), the Hodge decomposition a = dV + δΨ + h splits the affinities into those round the blocks (δΨ) and those round the holes (h). dV carries none.

*Proof.* (a) is Kolmogorov's criterion. (b) is Schnakenberg's network theorem. (c) holds because ∮_C dV = 0 for every cycle, and the periods of a closed form over a basis of cycles determine its cohomology class. ∎

**Remark 6.4.** This is the prediction of *Street Spin Foam* (the user's paper): circulation is sourced by the topological part of the valuation field alone, so that a region of large local disagreement and zero winding carries zero net circulation. Here it is restated on the infrastructure's Hodge decomposition (physis `discretisation`), with valuation read as log rate ratios. A model whose movement follows the gradient of a potential cannot circulate, whatever the potential.

**Representation 6.5** (the crowd as a kinematic wave; Withholding 2.1 marks its standing). A segment is a link with a fundamental diagram:
- free speed 1.34 m/s and jam density 5.4 m⁻² (Weidmann 1993), times the walkable width;
- Godunov's scheme (physis `continuum`);
- a rule that opens a carriageway on a schedule widens the link while it holds ([SO] §7.2; [G] §11: Ginza's Chūō-dōri, about 27 m, Saturdays, Sundays and holidays from 12:00 to 18:00, or 17:00 from October to March).

The representation shows the conversion of the same arrivals into different densities, queues and travel times under the two configurations. It is offered for exclusion: for instance, that the alternation could matter by capacity alone.

The workbench applies it to a street of its own, Swanston Street from Flinders Street to La Trobe Street. The street's blocks end at the cross streets where the city model's named ways cross its axis. A signal stands at each street, and none at a lane or a mall. Walkers enter at the mean rate counted northward at the street's start. A window of the day is kept second by second, holding the flow through every interface and the walkers' speed F/ρ there (zero at a red crossing). Nobody turns into a cross street, so the model's flow further up overstates the counts there. The difference measures what the cross streets take.

### 6.1 Currents Past a Counting Line

**Definition 6.6** (the current meter; `flow/currents.hpp`). A sensor counts d₁ walkers an hour toward its first direction and d₂ toward its second. Their means at each civil hour of the week give the transport d₁ + d₂, the net current u = d₁ − d₂, and the *residual current*, the mean of u over the week. They also give the hours of the mean weekday at which u is greatest (the *flood*, toward the first direction) and least (the *ebb*). For whole weeks of record, the means at the 168 hours are the projection onto the mean and the week's 84 harmonics.

**Proposition 6.7** (the measured edges bound the entropy production from below). Take the walkers as independent, each a Markov jump process among places. N of them jump from i to j at the mean rate n_ij = N p_i k_ij. Their entropy production in k_B per unit time is, at each instant (Schnakenberg 1976, eq. 7.6; Seifert 2012),

  σ = Σ_{edges {i, j}} (n_ij − n_ji) ln(n_ij / n_ji).

Every term is non-negative, and vanishes exactly when the edge balances. A counting line measures one edge's two rates, so the terms of the measured edges sum to a lower bound on σ.

*Proof.* N independent copies multiply Schnakenberg's single-walker σ by N, and p_i k_ij by N is n_ij. Each term has the form (x − y) ln(x/y) ≥ 0, since x − y and ln x − ln y have one sign. Dropping the unmeasured edges' terms lowers the sum. ∎

*Remark 6.8.* σ measures one-wayness, not crowding: a street with as many walking each way adds nothing, whatever their number. On Melbourne's 92 counting lines the bound peaks at about 43,000 k_B an hour on a weekday at 08:30, and falls to about 60 at 03:30. The checks hold the edges' terms equal to physis's σ from a rate matrix on a driven ring, to 10⁻¹³.

**Definition 6.9** (handedness). A sensor named for its side of a street ("Town Hall (West)") counts along one footpath. If walkers chose their side as left-hand traffic chooses its lane, the residual current would run east on north footpaths, west on south, north on west and south on east. With no handedness, each footpath's residual runs either way with probability ½. The number agreeing with the left-hand rule is then binomial(n, ½), a sign test. On Melbourne's 41 such footpaths, 29 agree (P(X ≥ 29) = 0.0058). The footpaths of one street share their walkers and are not independent, so the test overstates its evidence by an unknown amount. Still, the street's mirror symmetry, which exchanges its two sides, is broken in the direction of the road's.

## 7. Deposition and Memory

**Definition 7.1** ([SO] §15). The *record of tenancy* is the sequence, by census year, of the uses at each holding. CLUE gives it for Melbourne from 2002, with its holdings placed. For a segment, the workbench computes three quantities:
- *turnover*, the share of holdings whose use changed between consecutive censuses;
- *persistence*, the survival of a use at a holding over k censuses;
- the series M(year) of Definition 2.6.

A claim that a street retains its character while its contents change ([SO] §1, after Molotch et al. 2000) is then a claim about M(year) and the use composition against turnover. Under the integration null (Definition 3.11), a slow drift of M carries no memory.

## 8. The Classification (the Field's First Test)

**Definition 8.1** ([SO] Claim 29.2, Figure 2). Each segment is a point (C, S) of the plane: compression ratio (Definition 2.4) and selectivity of response (Definition 4.2). The field's first test is whether the points form groups. The workbench applies Clark and Evans's ratio R of mean nearest-neighbour distance to its expectation under complete spatial randomness on the plane's rectangle, after scaling both axes to [0, 1], with Donnelly's edge correction.
- R < 1, significantly: the segments cluster, and the vocabulary has sorted its object.
- R ≈ 1: they fill the plane without structure, which is the paper's first condition of failure ([SO] §29.3).

**Caveat 8.2.** Melbourne's sensors were placed for counting pedestrians, not by the principle of selection of [SO] §26.2. Their plane tests the vocabulary on an unchosen sample. It does not test the field's design, and it includes no case chosen to be uninteresting ([SO] Claim 26.1); Macaulay Road comes closest. Every report of the plane says so.

## 9. The Obligations of a Report

The paper sets three obligations on every report ([SO] §17.1, Requirement 28.1, Requirement 33.2):
- the instrument and its blind spots;
- the null account the finding competes with, and the separating observation;
- the uses the result makes available to a party with authority.

The workbench keeps the first two in plain words. Each result card says whether the result was measured, proved or modelled, how it was obtained, what it was tested against, and its limits. The papers' apparatus stays in this note, at the user's request.

## 10. Amendments to the Papers

1. **Claim 9.1 is a theorem, and stronger.** No nonzero weighting is invariant under changes of unit. The weighting-free content of a terrain survey is its dominance order. Where arrivals differ in the sign of a component's effect, nothing is comparable without a weighting (Theorem 2.2).
2. **L_grad needs a definition.** The integral scale of the use correlation is proposed (Definition 2.4). The classification is invariant to L_walk while L_walk is common to the streets (Remark 2.5).
3. **The departure from mixing is constructed.** It is the mutual information between place and use, exactly decomposable across scales (Definition 2.6, Theorem 2.7). This closes the open item of [SO] §36.2.
4. **The persistence criterion gets its form.** Exchange blind to kind never increases the departure (Theorem 2.8). A persistent gradient therefore needs a kind-dependent process, and the criterion is a balance of production against exchange (Corollary 2.9). On a discrete street, fronts can be pinned (Representation 2.10).
5. **The day and the week are one calendar.** Their combinations are harmonics of the week. The torus of independent calendars is T³ (Remark 3.3).
6. **The record-length constraint is about separation.** Separating neighbouring lines needs T ≥ 1/|Δf|, which can far exceed either period. The tropical and the civil year are never separable (Remark 3.5).
7. **Discriminator 29.1, middle case.** It needs a third limitation beside record length and the enumeration of drivings: the forcing must be separable on the calendar torus, or its declared non-separable structure must be removed first. Institutional forcings, such as hours by weekday and by season, are non-separable by construction (Proposition 3.7, Remark 3.8).
8. **Discriminator 29.1, third case.** It needs an order of explanation: with two incommensurate drivings every frequency is a combination of high enough order. The discriminator's chance rate χ(d, T) must be reported with it (Proposition 3.9).
9. **Bicoherence discriminates only the irregular part.** The lines of a calendar are phase-locked whatever their origin (Remark 3.10).
10. **Selectivity of response is defined** as the normalised inverse participation ratio over the drivings (Definition 4.2).
11. **The exact budget of person-time is Little's identity**, with its boundary terms (Theorem 6.1).
12. **Circulation is topological.** The irreversibility of a movement process depends only on the non-exact part of its log rate ratios. This joins *Street Oceanography* §13 with *Street Spin Foam*'s prediction on one decomposition (Theorem 6.3).

## 11. What Remains Open

- **The instrument of [SO] §20** (carried properties at entry and at departure, with a comparison population) is not built, and nothing here substitutes for it. Every representation of composition (2.14, 5.2) is declared, not measured.
- **The terrain's symbolic component** is registered here only through use classes, from ANZSIC codes and OSM tags. Display and description ([SO] §21) are not read.
- **Transit times** are not measured: Melbourne's sensors count crossings at a line, not entries and exits of a segment. Definition 6.2 is computed on representations.
- **The response function moves under intervention** ([SO] Claim 25.1). Every susceptibility reported is a statement about what a street did, under the arrangement in force.
- **Whether the plane of §8 separates on a designed sample** (with the uninteresting case) is the question the field stands on. Melbourne's sensors can only begin it.

## References

- Boyd, S. and Chua, L. O. (1985). Fading memory and the problem of approximating nonlinear operators with Volterra series. *IEEE Trans. Circuits Syst.* 32, 1150–1161.
- Clark, P. J. and Evans, F. C. (1954). Distance to nearest neighbor as a measure of spatial relationships in populations. *Ecology* 35, 445–453.
- Cover, T. M. and Thomas, J. A. (2006). *Elements of Information Theory*, 2nd ed. Wiley.
- Fáth, G. (1998). Propagation failure of traveling waves in a discrete bistable medium. *Physica D* 116, 176–190.
- Foreman, M. G. G. (1977). *Manual for Tidal Heights Analysis and Prediction*. Pacific Marine Science Report 77-10, Institute of Ocean Sciences, Victoria, B.C.
- Godin, G. (1972). *The Analysis of Tides*. University of Toronto Press.
- Hasselmann, K. (1976). Stochastic climate models. Part I: Theory. *Tellus* 28, 473–485.
- Huang, W. (2026). *Street Oceanography: An Oceanographic Epistemology of Socio-Cultural Dynamical Systems*. Working paper. [SO]
- Huang, W. (2026). *Towards Street Oceanography: Ginza Street as a Case Study*. Working paper. [G]
- Huang, W. (2026). *Street Spin Foam: A Preliminary Formulation*. Working paper.
- Keener, J. P. (1987). Propagation and its failure in coupled systems of discrete excitable cells. *SIAM J. Appl. Math.* 47, 556–572.
- Kim, S.-H. and Whitt, W. (2013). Statistical analysis with Little's law. *Oper. Res.* 61, 1030–1045.
- Kim, Y. C. and Powers, E. J. (1979). Digital bispectral analysis and its applications to nonlinear wave interactions. *IEEE Trans. Plasma Sci.* 7, 120–131.
- Leslie, T. F. and Kronenfeld, B. J. (2011). The colocation quotient: a new measure of spatial association between categorical subsets of points. *Geogr. Anal.* 43, 306–326.
- Little, J. D. C. (1961). A proof for the queuing formula L = λW. *Oper. Res.* 9, 383–387.
- McFadden, D. (1974). Conditional logit analysis of qualitative choice behavior. In P. Zarembka (ed.), *Frontiers in Econometrics*, 105–142. Academic Press.
- Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed. Willmann-Bell.
- Molotch, H., Freudenburg, W. and Paulsen, K. E. (2000). History repeats itself, but how? *Am. Sociol. Rev.* 65, 791–823.
- Mora, R. and Ruiz-Castillo, J. (2011). Entropy-based segregation indices. *Sociol. Methodol.* 41, 159–194.
- Neuts, M. F. (1981). *Matrix-Geometric Solutions in Stochastic Models*. Johns Hopkins University Press.
- Schnakenberg, J. (1976). Network theory of microscopic and macroscopic behavior of master equation systems. *Rev. Mod. Phys.* 48, 571–585.
- Tennekes, H. and Lumley, J. L. (1972). *A First Course in Turbulence*. MIT Press.
- Theil, H. and Finizza, A. J. (1971). A note on the measurement of racial integration of schools by means of informational concepts. *J. Math. Sociol.* 1, 187–193.
- Theiler, J., Eubank, S., Longtin, A., Galdrikian, B. and Farmer, J. D. (1992). Testing for nonlinearity in time series: the method of surrogate data. *Physica D* 58, 77–94.
- Turing, A. M. (1952). The chemical basis of morphogenesis. *Phil. Trans. R. Soc. B* 237, 37–72.
- Weidmann, U. (1993). *Transporttechnik der Fussgänger*. Schriftenreihe des IVT 90, ETH Zürich.
