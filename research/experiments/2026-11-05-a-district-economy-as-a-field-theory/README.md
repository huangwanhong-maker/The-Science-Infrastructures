# A District's Economy as a Field Theory

*Status: draft (2026-11-05). Research.* The third of three case studies of one shared problem, [a district's economy](../2026-11-05-a-district-economy-the-problem/README.md): the ground floors of the blocks between Swanston and Russell Streets, Flinders Lane and Bourke Street, in Melbourne's Hoddle Grid. Here the people are a field on the district's space–time complex. Their law is the stationary point of an action S, and their fluctuations are its Euclidean measure. The plan is a second field coupled to the first, and the searches use the action's structure.

This is a study of methods on a real place, not a planning proposal. Its economics is the shared problem's sketch, with every economic parameter an assumption (the problem's §4.3). Its people are a diffusive fluid calibrated to seventeen counters, and its plans say nothing about what any frontage should hold.

## 1. The Question and Why

The user asked for a plan of the district's economy that "utilizes a fluid field theory model with action S", with space, time and frequency, people flow, and the arrangement of the buildings and their relations all considered.

The shared problem fixes everything but the people (problem §1). A study supplies only the passers P_f(h; plan) at each frontage f and hour h of the mean week, and the problem's `evaluate()` gives the value, the week's lines of the spending, the mix and access. The three studies differ only in how much structure their model of the people has:

1. without a structured model ([study 1](../2026-11-05-a-district-economy-without-a-model/README.md));
2. a conserved flow on the network, with its circulation (study 2);
3. a field theory with an action, here.

The questions are these:

1. **The field.** What does a field theory of the people, calibrated to the counters, say of the passers at the 61 frontages, which hold no counter, and how sure can one be of it?
2. **The modes.** What do the action's modes (Hodge eigenmodes in space, Fourier lines in time), its conserved quantities and its fluctuations show of the district?
3. **The search.** What does the action's structure buy a search, at equal budgets with study (1)'s black-box searches?

*Rationale.* A field theory is the most structured of the three models. It names its degrees of freedom, its symmetries and its conservation laws, and its Gaussian part is solved exactly. If that structure is worth anything to planning, it should show in the calibration's honesty about what it does not know, and in searches that can use the locality of the action rather than treat the value as a black box.

## 2. Theory

### 2.1 The Base: the District's Space–Time Complex

- **Space.** X is the urban-science domain's district, a regular complex (`urban::landscape::district`, on the shared problem's network): 356 vertices (0-cells), 399 segments (1-cells) and F = 37 blocks (2-cells).
  - **The lattice.** Its segments are up to 224 m long, against a screening length of 65 m, so the study cuts every segment into equal pieces of at most 15 m. This is a subdivision, which changes no homology. The lattice has V = 1,029 vertices and E = 1,072 pieces, and each block's boundary runs through its segments' pieces (§8.3).
  - **Links.** A plan's links are new 1-cells on this lattice (§2.4).
  - **Its homology.** b₀ = 6: the main network and five lane pieces that no segment joins to it. b₁ = 12: the cycles that no block fills (the faces the domain leaves as holes). b₂ = 0. Euler: V − E + F = b₀ − b₁ + b₂ = −6, on the domain's cells and on the lattice alike.
  - **Orientation.** The cellular coboundaries d₀ : C⁰ → C¹ and d₁ : C¹ → C² (morphe's `chain_complex`) satisfy d₁d₀ = 0.
- **Time.** Z₁₆₈, the cycle of the week's hours: the mean week is periodic. Its modes are e^{iω_m t}, with ω_m = 2πm/168 rad h⁻¹ and m = 0 … 84.
- **Space–time.** 𝒳 = X × C₁₆₈, with spatial cells (σ, t) and temporal edges (v, t → t + 1).
- **Hodge stars.** These are diagonal:
  - ⋆₀(v) = ½ Σ_{e∋v} ℓ_e, the dual length of a vertex (m);
  - ⋆₁(e) = D c_e / ℓ_e, a conductance, with c_e by the Hoddle Grid's hierarchy: 1 on its 99-foot main streets, c_little on its 33-foot little streets, c_lane on the lanes. Leave-one-site-out chose c_little = c_lane = 1 (§3.1).
  - ⋆₂(B) = 1/area(B) on the blocks.
- **The Hodge Laplacians.** Δ_p = d_{p−1}d_{p−1}* + d_p*d_p, with d_p* = ⋆_p⁻¹ d_pᵀ ⋆_{p+1} (physis's `hodge_laplacian`).
  - Δ₀ = ⋆₀⁻¹ d₀ᵀ ⋆₁ d₀ acts on densities. Its eigenpairs Δ₀u_k = λ_k u_k are ⋆₀-orthonormal, u_kᵀ⋆₀u_l = δ_kl, and its kernel is spanned by the components' indicators (b₀ zero modes).
  - Δ₁ acts on currents. Its spectrum is Δ₀'s nonzero spectrum (gradients, V − b₀ = 1,023 modes on the lattice), Δ₂'s nonzero spectrum (curls round the blocks, F − b₂ = 37 modes) and b₁ = 12 harmonic modes (round the holes): 1,023 + 37 + 12 = 1,072, by Hodge's theorem for a finite complex (Eckmann 1944).

### 2.2 The Action

**The fields.**
- φ_v(t): people passing a point on one footpath an hour, the reading of a counter and the passers of a frontage. It is proportional to the line density of walkers, ρ = 2φ/v₀ people per metre on both footpaths, at walking speed v₀ (both footpaths alike, an assumption). The number of people in v's dual cell is (2/v₀)⋆₀(v)φ_v.
- J_e(t): the net current on segment e (people an hour, tail to head).
- r_v(t): the exchange of v with what is not the street: buildings, trams and stations, and the city beyond the area's edges.
- ψ_v(t): the Lagrange multiplier of continuity, Martin–Siggia–Rose's response field.

**The constraint.** Continuity is ⋆₀ ∂_t φ = d₀ᵀ J + r (in units in which ⋆₀φ counts people × v₀/2). It is the statement that the space–time current 𝒥 = (⋆₀φ on temporal edges, J on spatial edges), a 1-chain on 𝒳, has as its boundary the exchange: ∂𝒥 = r.

**The action**, on the week's circle:

```
S[φ, J, r, ψ] = ∫_{S¹} dt { Σ_e (J_e − 𝒥_e[φ; U])² / (4T ⋆₁(e) φ̄_e)          transport
                          + Σ_v (r_v − s_v + μ ⋆₀(v) φ_v)² / (4T Γ_v)       exchange, with the sources s
                          + i Σ_v ψ_v (⋆₀ ∂_t φ − d₀ᵀ J − r)_v }            continuity

𝒥 = −⋆₁ (d₀φ + φ̄ d₀U),   U_v(t) = −b Σ_f p_f(t) W_f(v)                     Fick's current and the plan's drift
```

**What each term is.**
- **Transport.** The current's deviation from its constitutive law, in the mobility σ_e = 2⋆₁φ̄_e of independent walkers: the action of macroscopic fluctuation theory (Bertini, De Sole, Gabrielli, Jona-Lasinio and Landim 2015), which on a graph is the large-deviation rate of the walkers' currents.
  - **The Benamou–Brenier form.** Write the transport term with θ_e the logarithmic mean of φ at e's ends (Maas 2011), so that θ_e d₀(ln φ) = d₀φ exactly, and leave U out. It is then
    - Σ_e J_e²/(4T⋆₁θ_e): Benamou and Brenier's kinetic term, ρ|v|²/2 with ρv = J;
    - plus Σ_e ⋆₁θ_e (d₀ ln φ)_e²/(4T): the Fisher information, Bohm's quantum potential in Madelung's fluid;
    - plus Σ_e J_e (d₀ ln φ)_e/(2T) = Σ_v ln φ_v (⋆₀∂_tφ − r)_v/(2T) = (1/2T)[d/dt Σ_v ⋆₀(φ ln φ − φ) − Σ_v r_v ln φ_v].
  - **On the closed week.** The total derivative integrates to zero, so the transport term is the Benamou–Brenier kinetic term plus the Fisher information, less an exchange term: the action of the Schrödinger problem, entropic optimal transport (Léonard 2014; Chen, Georgiou and Pavon 2016).
- **Exchange.** The linear-noise (Gaussian) form of walkers born at the sources at rate s_v and leaving the street at rate μ, with noise Γ_v.
- **The plan's potential.** U lowers the potential at every open frontage by b times p_f(t), the chance that a passer stops there: the shared problem's, by kind, band of the day and catalysis. W_f are f's Whitney weights on its segment. Walkers drift down U, towards open frontages. b = 5 is an assumption (§5).
- **The sources** s_v(t). They are not observed: they are calibrated to the counters (§3.1).
- **Crowding**, the interaction. Two forms are possible:
  - exclusion, whose mobility saturates as σ = 2⋆₁φ̄(1 − ρ/ρ_max);
  - the fall of walking speed with density (Weidmann 1993).

  Both are measured in §2.6, and both are negligible in this district.

**What is assumed.**
- The walkers are independent, apart from crowding.
- The noise is Gaussian (the linear-noise approximation about the mean density φ̄).
- Both footpaths carry alike.
- People leave the street at one rate μ everywhere.
- The plan acts on the people only through U: through the attraction of its open frontages, and through the network its links make.

### 2.3 Symmetries and Conservation Laws

- **The number of people (Noether).** S is invariant under ψ ↦ ψ + c 1_K, a constant shift of the response field on a component K of X. It changes by i c ∫ Σ_{v∈K} (⋆₀∂_tφ − d₀ᵀJ − r)_v, and Σ_{v∈K}(d₀ᵀJ)_v = (d₀1_K)ᵀJ = 0, since constants on a component are closed.
  - **Its current** is the space–time current 𝒥.
  - **Its charge** is the number of people on the component's streets, N_K = (2/v₀)Σ_{v∈K} ⋆₀φ_v. On every path and in every fluctuation it obeys

    dN_K/dt = (2/v₀) Σ_{v∈K} r_v,  and on the stationary path  dN_K/dt = (2/v₀) Σ_{v∈K} s_v − μ N_K.

  - **The charges are Δ₀'s kernel**, H⁰(X; ℝ): b₀ = 6 conserved numbers, one per component. In the modes, N̂(m) = (2/v₀)Σ ŝ(m)/(iω_m + μ). The run holds this at every line to 6 × 10⁻¹⁶.
- **Time translation.** The action without its sources and its plan is invariant under t ↦ t + 1 on Z₁₆₈, so the Fourier modes diagonalise it. The week's rhythm is the breaking of this symmetry by the sources and the plan's hours.
- **The currents' gauge.** Continuity sees only d₀ᵀJ. The currents' co-exact part (circulation round the blocks, 37 dimensions) and harmonic part (round the holes, b₁ = 12) are invisible to the density and to every counter.
  - The stationary current J = −⋆₁(d₀φ + φ̄d₀U) is a weighted gradient: Fick's current carries no circulation.
  - Its circulation and harmonic parts are fluctuations alone. This is a structural difference from study (2), whose conserved flow does circulate (§7.4).
- **No spatial symmetry.** The weighted complex has no automorphism to speak of, so the modes are not labelled by momenta. They are labelled by their rates Dλ_k, and their wavelengths 2π/√λ_k run from the district's width down to its shortest segments.

### 2.4 The Stationary Path and Its Modes

The Euler–Lagrange equations of S are satisfied by ψ = 0, J = 𝒥, r = s − μ⋆₀φ and continuity. On that path S = 0, its minimum on the real section. The path is the hydrodynamic equation, a diffusion with drift and exchange:

```
⋆₀ ∂_t φ = −L φ − μ ⋆₀ φ + s + L_φ̄ U,     L = d₀ᵀ ⋆₁ d₀ = ⋆₀ Δ₀,   L_φ̄ = d₀ᵀ diag(⋆₁ φ̄) d₀.
```

Without the plan's drift, it is diagonal in the modes. For sources band-limited to the week's 84 lines, this is the exact continuous-time solution, read at the hours:

```
φ̂(m) = Σ_k u_k  t_k(m)  u_kᵀ ŝ(m),     t_k(m) = 1 / (iω_m + D λ_k + μ).
```

- **The screening length.** ξ = √(D/μ). Leave-one-site-out chose ξ = 65 m.
- **The time constant.** The slowest nonzero mode relaxes at a rate Dλ_k + μ ≥ μ = 6 h⁻¹. So the district answers its sources within ten minutes: at the daily line the filter's modulus is 0.999, at six cycles a day 0.97 and at twelve 0.89 (search-82).
- **The plan's drift.** It is a small perturbation, b p ≤ 0.13, and is taken in linear response. The response of the passers to the attraction of g at hour h is

  R_fg(h) = b [O G₀ L_φ̄(h) W]_fg,  G₀ = Σ_k u_k u_kᵀ/(Dλ_k + μ),

  quasi-static (its error is of order (ω/μ)², 0.2% at the daily line). It is computed once on the network as built: a link changes it only near the link, a second-order effect.
- **The method of solution.** The path is solved exactly in the modes, by physis's Hodge Laplacian and the study's mode sums. physis's mechanics, which steps second-order discrete Euler–Lagrange equations forward from initial data, does not fit a first-order problem on a closed time. That is a gap (§8.2).

### 2.5 The Euclidean Measure: the Gaussian Theory

Integrating out J and r, which enter quadratically, leaves the Onsager–Machlup action of φ:

```
S_OM[φ] = (1/4T) ∫ ‖ ⋆₀∂_tφ + Lφ + μ⋆₀φ − s − L_φ̄U ‖²_{Γ_tot⁻¹},   Γ_tot = d₀ᵀ diag(2⋆₁φ̄) d₀ + diag(2Γ).
```

Its Gaussian measure e^{−S_OM} has the stationary path as its mean.

**The prior.** The study takes the noise white per unit length, Γ_tot → σ²(m) ⋆₀ at each line: the sources' residual is a white field with a spectrum in time. In the modes the measure is then diagonal:

```
⟨c_k(m) c̄_l(m')⟩_c = δ_kl δ_mm' σ²(m) / |iω_m + Dλ_k + μ|²,     c_k(m) = u_kᵀ ⋆₀ φ̂(m).
```

This is a Whittle–Matérn field on the network, the stochastic-PDE prior of Lindgren, Rue and Lindström (2011), whose precision is the action.

**The exact Gaussian expectations, as mode sums.**
- **The prior covariance** between any two points of the network (a counter, a frontage) is K_xy(m) = Σ_k B_xk |t_k(m)|² B_yk, with B = O U the modes read by Whitney's interpolation.
- **The posterior given the counters** is the Gaussian conditioning of this measure (§3.1). The lines decouple, because the measure is stationary in time and every counter reports every hour of its mean week.
- **Connected correlations are responses** (fluctuation–dissipation). In the Gaussian theory ∂⟨φ_x⟩/∂s_y = Σ_k u_k(x) t_k u_k(y) = G_xy, the propagator, which is also the field–source connected correlation ⟨φ_x ; ξ_y⟩/σ²⋆₀(y). Two consequences:
  - **A link** adds an edge of conductance κ between its ends. Its effect on the mean is −κ G b bᵀφ to first order (the connected correlation of φ with the new edge's energy ½κ(bᵀφ)²), and −κ G b (bᵀφ)/(1 + κ bᵀGb) exactly (Sherman–Morrison) for an edge between existing vertices (§7.6).
  - **The value's sensitivity** to each frontage's attraction is the adjoint Ā(h) = R(h)ᵀ A(h) (§2.7).
- **The free energy of the value.** For a plan held fixed, the value V = Σ A_f(h) P_f(h) − costs is linear in the passers. Under a Gaussian law of P it is Gaussian, so

  F_T = −T ln ⟨e^{−V/T}⟩ = E V − Var V/(2T)

  exactly: the free energy of the value as a source of the people's field, a planner's certainty equivalent at risk tolerance T. A free-energy difference ranks two plans by their expected difference less half its variance over T (§7.7).

### 2.6 The Interaction, and What Lies beyond the Gaussian

- **Crowding.** The busiest footpath-hour of the district is the weekend high water at the Town Hall (West), 3,841 passers an hour (problem §5.1). At v₀ = 1.34 m/s on a footpath 4 m wide (an assumption) its density is ρ = 3,841/(4,824 × 4) = 0.20 people/m²: Fruin's level of service A.
  - **The speed.** Weidmann's fundamental diagram, v/v₀ = 1 − exp(−1.913 (1/ρ − 1/5.4)), lowers the walking speed there by 1 × 10⁻⁴.
  - **The mobility.** Exclusion's mobility factor 1 − ρ/ρ_max is 0.963, a 3.7% change in the noise and the drift's response at that one place and hour. It leaves Fick's mean current unchanged (in the symmetric exclusion process the hydrodynamic equation stays linear).
  - **The verdict.** The people's theory is free to the order that matters for the economy. Its stationary path needs no Newton iteration, and its fluctuations need no Monte Carlo.
- **Positivity.** The ideal gas's entropy ρ ln ρ is not Gaussian where φ is small (night). The Gaussian field can dip below zero there, and the study holds passers at zero. In random plans 40 to 80 of the 10,248 frontage-hours are held (check_field).
- **Where the interaction is.** It is in the plan. The value's catalysis and competition couple frontages pairwise and more, and they make the plan a strongly interacting field (§2.7).

### 2.7 The Plan as a Field Coupled to the People

- **The plan field.** σ_f = (k_f, weekday hours, weekend hours) on each of the 61 frontages, and a bit for each link.
- **The value as its energy.**

  V(σ) = Σ_h Σ_f A_f(h; σ) [P_f(h; σ)]₊ − Σ_f C_f(σ_f) − Σ_links,  A_f(h) = [open] p_f c_f s_k m_k,  P = P⁰ + R(p − p⁰).

  A_f(h) is the problem's value per passer: the chance to stop with catalysis (p), competition's share (c), spend and margin. P⁰ is the field's stationary path as built, carried to the plan's network.
- **Its couplings.** A frontage is coupled only to:
  - its catalysis neighbours within 3ξ_cat = 127.6 m (the problem's relation *reach*, 335 pairs);
  - its competitors of one kind on its segment (*one segment*);
  - and, through R, everyone at the same hour, weakly.
- **Locality in time.** Every term lives at one hour: the value is a sum over hours, and a frontage's hours enter only through which hours it is open.
- **The exact local conditionals.** Fix all but frontage f. For each kind k, f's value is a constant plus the sums of an hourly margin m_k(h) over its weekday and weekend intervals, less their staffing. So all of f's 5 × 600 × 600 + 1 states (kinds × weekday (s, ℓ) × weekend (s, ℓ), and vacant) are valued exactly from 5 × 168 margins and a prefix sum.
  - **Each margin.** Its exact change is ΔV_h = Σ_g (A'_g [P_g + Σ_{g'∈N(f)} R_gg' Δp_g']₊ − A_g[P_g]₊): n|N(f)| products. Where no passer is held at zero, this is Σ_{g∈N(f)} Ā_g Δp_g plus local terms, with Ā = RᵀA the adjoint of the people's response.
  - **The heat bath.** The plan's Gibbs measure e^{V/T} has these conditionals. The study draws f's state from them exactly (a heat bath on the encoding's states as they are, a vacant frontage's 600² hour genes all alike), anneals T geometrically from A$ 3,000 to A$ 3 over 40 sweeps, and then sets T = 0 (iterated conditional modes) until a sweep changes nothing.
  - **The links.** Each link is drawn from its own two-state conditional, by an exact evaluation of the plan with and without it.

### 2.8 What Is Exact and What Is Approximate

| Exact (to rounding, held by check_field) | Approximate, and how far |
|---|---|
| the stationary path of the free theory in the modes; its conservation law at every line | the noise is Gaussian about the mean density (linear noise); positivity is imposed by holding passers at zero |
| the Gaussian conditioning on the counters, line by line; the posterior covariance at the frontages | the prior of the sources' residual is white per unit length, with a spectrum in time fitted line by line (an assumption of form) |
| the free energy F_T = EV − VarV/2T of a fixed plan under the Gaussian law | the plan's attraction in linear response, quasi-static (error ~(ω/μ)²) and on the network as built (second order in link × attraction) |
| each frontage's local conditional, and the value after every heat-bath update (equal to the problem's `evaluate()`) | crowding neglected (10⁻⁴ in speed, 4% in mobility at the busiest footpath-hour) |
| Sherman–Morrison's resummation of a new edge between existing vertices | a link drawn as one edge between its ends' Whitney points, against the domain's rebuilt complex (§7.6) |

## 3. Method

### 3.1 Calibration to the Counters

- **The observations.** The 17 counters on the district's network (those in its box, problem §5.1) observe φ at their feet on their nearest segment, by Whitney's interpolation. A counter of a crossing stands on the way it crosses to. Their data are their mean weeks, as Fourier lines ŷ_c(m), m = 0 … 84.
- **The sources.** They are a sum of source groups plus a residual. The groups are:
  - the area's edges (degree-one vertices at the box, weighted by the width class of their way);
  - the buildings' volume (each parcel's and neighbour's footprint × height, on its nearest segment);
  - the main streets', the little streets' and the lanes' length.

  Each group has its own complex line β_γ(m). The residual is white per unit length with variance σ²(m)⋆₀.
- **Line by line, universal kriging.** For each line:
  - generalised least squares gives β̂(m);
  - Gaussian conditioning on the counters gives the residual's posterior;
  - the posterior mean of the sources follows in closed form, ŝ(m) = Σ_γ a_γ β̂_γ + ⋆₀ ⊙ U(t̄ ⊙ Bᵀ Σ₁⁻¹ r);
  - REML gives σ̂²(m).

  The nugget τ² is relative, τ² × (level²/median level²) × K_cc, so that the fit weighs errors as the score does.
- **The 38 counters around the area.** They enter one variant of the edge group, feeding each of the area's edges by the mean week of its nearest outside counter. Leave-one-site-out preferred the edges weighted by width alone, so in the chosen calibration they do not enter.
- **Parameters chosen by leave-one-site-out.** Counters within 30 m are one site (the three at Flinders Lane and Degraves Street). A grid of 768 points was scored by the normalised error mean_c Σ_h(ŷ − y)²/Σ_h y²:
  - ξ ∈ {30, …, 420} m;
  - four ways' hierarchies (c_little, c_lane);
  - three nuggets, absolute or relative;
  - edges by width or from the outside counters;
  - with or without a group for the tram streets.
- **The baseline.** Each counter's week taken from the nearest counter of another site, which is roughly what a study without a model of space would do.

### 3.2 The Economy

- **The passers.** For a plan, P_f(h) = [P⁰_f(h; links) + Σ_g R_fg(h)(p_g(h) − p⁰_g(h))]₊.
  - P⁰(links) is the calibrated sources carried by vertex positions to the plan's network (the problem's `with_links`), and its stationary path read at the frontages' feet.
  - Every one of the 2,048 link sets is computed once and cached.
- **The value.** The passers go to the shared problem's `evaluate()`, which gives the value and every other measure. The study's own sum of the value (`plan.hpp`) is held equal to it.

### 3.3 The Searches

**At equal budgets of 20,000 evaluations, seeds 1–5.** One evaluation is one call of the black box's objective: the field's passers for a plan and the problem's `evaluate()`, 0.83 ms on one thread of the development machine with the caches warm.

**Study (1)'s black-box searches, with its settings** (its `without_a_model.cpp`):
- the GA on the mixed encoding (μ = λ = 100, two elites, crossover 0.9, tournaments of two, the problem's kernels);
- the GA from the plan as built (its first population the plan as built mutated five times);
- separable CMA-ES on the canonical relaxation (σ₀ = 0.3, bounds by transformation, restarts from random means);
- random search.

**Two more black boxes:**
- simulated annealing from the plan as built with the problem's mutation (heuresis's modified Lam schedule from T₀ = A$ 2,000);
- the categorical UMDA on the encoding as a product of categorical variables (λ = 100, μ = 30, with borders: heuresis's distributions, the information-geometric family that the plan's Gibbs measure suggests).

**The structured search.** The heat bath on the plan field (§2.7), from the plan as built and from a random plan.
- **Its accounting.** Its local work is converted at the measured cost of one evaluation, and each full rebuild of its state (one after each sweep, and one for each link tried) counts as one.
- **Its budget.** It needs about a tenth of the budget and stops, so it runs under budget.

## 4. Code Version

| Library | Commit |
|---|---|
| The Science Infrastructures (root) | `512477b` |
| calx | `90aa710` |
| morphe | `14b3fb8` |
| physis | `e2b98d8` |
| tyche | `b4034b8` |
| heuresis | `5730e70` |
| domains/urban-science, domains/architecture | unchanged since `485f82c` |
| the shared problem (`economy.hpp`, its `results/`) | as written on 2026-11-05, not yet committed at `512477b` |

Built with g++ 11.5 and CMake 3.28 in WSL, in Release, without warnings from this study's code. The run takes 13 minutes on four threads of the development machine, shared with the other two studies, and the links' stage 3 more. Most of it is the calibration's grid (2.7 minutes) and separable CMA-ES, whose five seeds run side by side (5 minutes).

## 5. Parameters

| Parameter | Value | Source |
|---|---|---|
| μ, the rate of leaving the street | 6 an hour (a stay of ten minutes) | assumption |
| ξ = √(D/μ), the screening length | 65 m; D = 25,350 m²/h (7.0 m²/s) | fitted: leave one site out, over 30–420 m |
| c_little, c_lane | 1, 1 | fitted: leave one site out, over {1, ⅓} × {1, ⅓, ⅙, ⅒} |
| τ², the nugget | 1 × the prior variance × the counter's squared level relative to the median | fitted: leave one site out, over {0.1, 0.3, 1}, absolute or relative |
| the source groups | the area's edges by width, the buildings' volume, the main streets', little streets' and lanes' length, the tram streets' length | chosen by leave one site out among four sets (edges by width or from the outside counters, with or without the tram streets); the tram streets (Swanston, Elizabeth, Collins, Bourke, Flinders) read from Public Transport Victoria's network map, an assumption |
| σ²(m), the residual's spectrum | REML at each line | fitted |
| the lattice's pieces | at most 15 m (ξ/4) | assumption, chosen to resolve ξ |
| v₀, walking speed | 1.34 m/s | Weidmann (1993) |
| b, the attraction | 5: an open frontage lowers the walkers' potential by 5 p_f(h) (5% more walkers at equilibrium where p = 1%) | assumption; the value is reported at b = 0 and 20 (§7.5) |
| the footpath's width, for crowding | 4 m | assumption (§2.6 only) |
| the fluctuations | the district's common deviation, its covariance between the hours of a week over 52 weeks; the local field's relative variance T_ζ (the median over the counters), its pooled autocorrelation, and the action's static correlation in space | the counters' residuals over the reference year (`counters.py`) |
| the heat bath | 40 sweeps from T = A$ 3,000 to A$ 3 (geometric), then T = 0 until a sweep changes nothing (at most 20) | assumption |
| the black boxes | as study (1) runs them; SA from T₀ = A$ 2,000, modified Lam; UMDA λ = 100, μ = 30 | study (1); assumption |

Every economic parameter is the shared problem's (its §4.3).

## 6. Data and Licences

| Data | Source | Licence | Retrieved |
|---|---|---|---|
| the district, its frontages, relations, links, the counters' mean weeks | the shared problem's `results/` (problem.json, sensors.json; its inputs from the data library's city model, CLUE and the street-sensor products) | ODbL 1.0 (OpenStreetMap) for the geometry; CC BY 4.0 (City of Melbourne) for the counts and CLUE | as the problem states (2026-10-30 to 2026-11-04) |
| the counters' hourly counts and their directions, over the reference year | City of Melbourne, Pedestrian Counting System, the data library's `street-counts` (read by `counters.py`) | the City's open data, CC BY 4.0, on the problem's understanding | 2026-10-30 |

`results/counters.json`, `calibration.json`, `field.json` and `search.json` are derived from the City's counts and the problem's committed results and are committed. `results/field-geometry.json` holds the network's positions and polylines (extracts of OpenStreetMap, ODbL): it is not committed (`.gitignore`), and `field_theory` writes it again.

## 7. Results

All numbers are from the run of 2026-11-05 (`field_theory … all`, then `… links`; results/).

### 7.1 The Field and Its Calibration (Spatial)

**The lattice.** The domain's district, with every segment cut into pieces of at most 15 m: 1,029 vertices, 1,072 pieces and 37 blocks, with the same homology (b₀ = 6, b₁ = 12, b₂ = 0). The domain's own cells, with segments up to 224 m long, do not resolve a screening length of 45–65 m (§7.6, §10).

**The calibration chosen by leave-one-site-out** (768 grid points):

| | ξ | c_little, c_lane | τ² | groups | error left out (normalised) | median R² of a week left out | RMS log error of weekly totals |
|---|---|---|---|---|---|---|---|
| the field | 65 m (D = 7.0 m²/s) | 1, 1 | 1, relative | the area's edges by width, the buildings' volume, the main, little and lane streets' lengths, the tram streets' | **0.449** | **0.17** | **0.49** (a factor 1.6) |
| the groups' regression alone | — | — | — | the same | 0.438 | 0.18 | — |
| the nearest counter of another site | — | — | — | — | 3.38 | −0.25 | 0.92 (a factor 2.5) |

- **Against the baseline.** The field predicts a held-out counter's week with an error 7.5 times smaller than the nearest other counter's week does.
  - Seven of the 17 counters are predicted at R² > 0.6 with their site left out: 124 Elizabeth Street 0.95, Little Collins Street at Swanston 0.90, 231 Bourke Street 0.87, 197 Elizabeth Street 0.85, Melbourne Central at Elizabeth 0.84, Chinatown at Swanston 0.72, Bourke Street Mall 0.64.
  - Its failures are the city's own structure. The two Swanston Street counters are predicted at 0.42 and 0.41 of their level: Swanston is the grid's pedestrian spine, which no source group singles out. The crossing at Degraves Street is predicted at 2.4 times its level: it counts people crossing a little street, not walking along it.
- **In sample.** The posterior interpolates the counters, median R² 0.96 (the nugget keeps it from interpolating exactly).
- **The kriged residual.** It adds nothing at a held-out site, where the groups' regression is as good, since the residual's correlation length (ξ) is short beside the 100–700 m between counters. It matters near the counters, where most frontages stand.
- **The 38 counters around the area** fed the edges in one variant. Leave-one-site-out preferred the edges weighted by width, so they do not enter the chosen calibration (§3.1).

**The field as built** (search-80). Weekly passers at the 61 frontages run from 51,700 to 160,100 (median 113,500). The plan as built is worth **A$ 159,412 a week** on the field's passers, with spending A$ 868,830, spread 0.413 (69 effective hours), mix 0.883 nats and access 151 m.

**What the counters leave unknown** (search-85). The posterior standard deviation of a frontage's weekly passers is a median 27% of its mean (16% to 103%). It is least beside the counters and greatest on the lanes and little streets between them. Propagated to the value (§7.7), it is A$ 41,119 a week for the plan as built.

### 7.2 The Modes (Frequency, in Space and Time)

- **Spatial frequencies** (search-81).
  - Δ₀'s rates Dλ_k run from 0.13 an hour (k = 6, the first above the six zero modes: a wavelength 2π/√λ of 2.8 km, a tilt across the district) to 15,000 an hour on the shortest pieces.
  - Only 38 nonzero modes (wavelengths above about 600 m) are slower than μ = 6 an hour. Every other mode relaxes faster than people leave, which is why ξ is the field's correlation length.
  - Δ₁'s 1,072 eigenvalues split exactly as Hodge's theorem says: 1,023 of them are Δ₀'s nonzero ones (gradients), 37 are Δ₂'s (curls round the blocks), and 12 are zero (harmonic, round the holes).
- **Temporal frequencies** (search-82).
  - The field's power is 57.9% in the mean line, 36.5% in the daily line D1, 2.6% in the day's harmonics D2–D12, and 3.0% in the week's other lines: the modulation of the days, weekday against weekend.
  - By spatial mode, 71.7% of the power is in the kernel (the numbers of people on each component) and 28.3% in the rest, concentrated in modes 8–131 (wavelengths of 2.0 km down to 180 m).
- **The action's filter in time.** |t_k(m)|/t_k(0) = (Dλ + μ)/|iω + Dλ + μ| is 0.999 at D1 for the slowest modes, 0.97 at six cycles a day and 0.89 at twelve. The district answers its sources with a lag of at most ten minutes, so the week's rhythm is its sources', and the action shapes space, not time.

### 7.3 People Flow

- **The conserved number** (search-83). The people on the district's streets, Noether's charge N(t), run from 14 at night to 6,253 at the Saturday crest, 2,528 on average. They obey dN/dt = Σs − μN at every line to 1.5 × 10⁻¹³. N(t) and the static limit Σs/μ are indistinguishable at the hourly scale (the ten-minute lag again).
- **The currents.** Fick's current is a weighted gradient: it runs down the density, from the busy streets into the lanes and out at the area's edges. On Saturday at 13:00 the largest is 562 people an hour.
- **Against the counters' directions.** All 17 counters of the district record the counts each way, so the net flow on each footpath over the reference year is known. Its correlation with Fick's current on the same segment is 0.01 on average (from −0.67 to +0.81), and its size is about three times Fick's (median RMS 92 against 29 walkers an hour).
  - The measured net flows are the commuters' tides. On Swanston Street at Flinders Lane the mean weekday's net flow is 880 walkers an hour one way at 08:00, and 185 the other way at 21:00.
  - A diffusive action has no such current. Its stationary current is the gradient of the density, and the circulating parts (co-exact, harmonic) that a directed flow needs are pure fluctuation in it (§2.3).
  - The field theory calibrated here knows where the people are, not where they are going. That is study (2)'s subject.

### 7.4 The Plan's Attraction, and What It Is Worth

- **The passers it draws.** In the heat bath's best plan, the attraction of open frontages draws 3.3% more passers past the 61 frontages over the week (search-88), from the streets around.
- **Its sensitivity.** The attraction b is an assumption, so the study reports the best plan's value at three strengths:

| b | the best plan's value (A$ a week) | the heat bath run afresh at that b (seed 1) | its kinds that differ from the best plan's |
|---|---|---|---|
| 0 | 644,336 | 633,486 | 9 of 61 |
| 5 (the study's) | 691,994 | 691,914 | — |
| 20 | 834,967 | 832,904 | 19 of 61 |

  - At b = 5 the attraction is worth A$ 47,700 a week (7%) to the best plan.
  - The plan as built is unchanged by b, by construction (the attraction enters as a change from it).
  - The heat bath is a local search: at b = 0 its single run found a plan A$ 10,850 below the b = 5 plan valued at b = 0.

### 7.5 The Arrangement of the Buildings and Their Relations

**The relations the plans use** (search-87):

| | catalysis pairs within reach, + and − | Σ a·w of the + pairs | competitors on one segment | party walls joining one kind |
|---|---|---|---|---|
| as built | 115 and 115 | 23.1 | 18 | 10 of 54 |
| the heat bath's best | 297 and 0 | 138.0 | 89 | 44 of 54 |
| the best black box (UMDA) | 249 and 3 | 117.7 | 82 | 40 of 54 |

- **The heat bath's best plan** is 48 retail frontages, 9 of food and drink and 4 vacant, with two links made: Rainbow Alley to Bourke Street, and a lane to Athenaeum Place. It is open a median 08:00–22:00 on weekdays and 09:00–23:00 at the weekend.
  - **It is a monoculture.** Under the shared problem's value, retail's catalysis of retail (CLQ 2.67, an affinity near the clamp) outweighs everything. The plan fills the reach graph with retail, accepts 89 pairs of competitors on one segment (each costs a share p/2 ≈ 0.5%), and keeps food and drink where food catalyses food and retail (Rainbow Alley and the Bourke Street corner).
  - **Every other measure pays for it.** The mix falls from 0.883 to 0.093 nats and access worsens from 151 m to 209 m. The spread rises from 0.413 to 0.510 (86 effective hours), because the longer hours fill the evening (evening share 0.15 to 0.23).
  - **The shared problem's value, searched alone, rewards a monoculture.** That is a property of the problem's economics (its catalysis and its clamp), not of this study's people, and a reason the problem offers the mix as a second objective.

### 7.6 The Links: Linear Response against the Network Rebuilt

**The model's link.** Each link is a new 1-cell of a lane's conductance between its ends' Whitney points on the lattice. Its exact effect is Woodbury's (Sherman–Morrison for one link): a k × k solve at each line, with no network rebuilt. The study's own resummation agrees with it to within 2% at every link; what remains is the passers held at zero.

**Three estimates of each link's effect.** The first-order linear response is the connected correlation of the field with the new edge's energy (FDT, §2.5). The domain's rebuilt complex is the problem's `with_links`, subdivided alike, with the sources carried by Whitney's weights. For each link alone, from the plan as built (search-89):

| | the change of the 61 frontages' weekly passers (of 6.86 million) | frontage by frontage, against the lattice's exact | the links ranked, against the lattice's exact (Spearman) |
|---|---|---|---|
| the lattice with the link (exact, the model's) | −115,900 (Baptist Place to Collins Street) to +15,200 (Rainbow Alley to Bourke Street) | — | — |
| first order (FDT), from the propagator alone | 1.6 to 6.4 times the exact | r = 1.000 at every link | 0.95 |
| the domain's rebuilt complex | 1.6 to 17 times the exact, two of opposite sign | median r = 0.88 | 0.50 |

- **The first order is proportional to the exact.** For one rank-one change the two differ, line by line, by the scalar 1/(1 + κ bᵀG_m b), which is 0.16 to 0.6 at the static line. So the propagator alone ranks the links almost exactly (Spearman 0.95), frontage by frontage in proportion, without rebuilding anything: the fluctuation–dissipation theorem at work.
- **The rebuilt complex drains more.** There the link is a street of its own length, with stars and the leaving μ along it. It is a sink that the massless Whitney edge lacks, so the rebuilt complex removes more people near both ends. The pattern is the same (median r = 0.88), the size is not.
  - The two are different representations of an arcade. Neither gives an arcade visitors of its own.
  - The searches used the model's link. With the rebuilt complex's, links would look worse still.
- **What the links are worth.**
  - Each link alone, from the plan as built, costs more than it brings: the value changes by −A$ 91 to −A$ 9,444 a week, its amortised cost included.
  - In the heat bath's best plan two links are made, Rainbow Alley to Bourke Street and a lane to Athenaeum Place. Unmaking them would lose A$ 2,800 and A$ 325 a week.
  - Their flows are nearly neutral (−A$ 296 and +A$ 397). They pay through catalysis: they bring retail frontages within one another's reach along the network (the problem's walking distances).
- **On the shared problem's caveat.** That a link is cheap and its worth depends almost entirely on the model of flow (problem §11) is borne out: here the flows make links slightly harmful, and the catalysis distances decide.

### 7.7 Free Energies: the Plans Ranked under the Euclidean Measure

**The week's fluctuations** (search-85).
- **The common deviation.** The district's common relative deviation from its mean week has a standard deviation of 12.3% at Wednesday 13:00 and 14.8% at Saturday 13:00, over the year's 52 weeks.
- **The local field.** After the common deviation is removed, the local relative field has variance T_ζ = 0.0171 (a standard deviation of 13%). Its autocorrelation is 0.62 at one hour, 0.39 at three and 0.43 at 24.
- **The action's spatial correlation.** At m = 0 it predicts the counters' local residual correlations with r = 0.38. The agreement is carried by the three counters at one site: the action predicts almost no correlation beyond 100 m, where the observed correlations scatter about −0.04. The common deviation's own definition, a weighted mean of the counters, pushes their residuals slightly negative.

**The plans' values, uncertainties and free energies** (A$ a week). The posterior s.d. is given the counters; the week's s.d. is common plus local:

| plan | E V | posterior s.d. | week's s.d. | F_T at T = 10⁴ | T = 3 × 10⁴ | T = 10⁵ | better than as built, given the counters: ΔE / posterior s.d. of the difference |
|---|---|---|---|---|---|---|---|
| as built | 159,412 | 41,119 | 42,306 | −14,620 | 101,401 | 142,009 | — |
| the heat bath's best | 691,994 | 105,033 | 110,514 | −470,268 | 304,573 | 575,767 | z = 8.12 |
| the best black box (UMDA) | 625,988 | 95,081 | 103,579 | −362,461 | 296,505 | 527,143 | z = 8.15 |
| the best SA | 554,239 | 85,814 | 90,731 | −225,573 | 294,302 | 476,258 | z = 7.98 |
| the best CMA-ES | 469,539 | 80,179 | 82,493 | −192,152 | 248,975 | 403,370 | z = 6.65 |
| the best GA | 390,114 | 66,971 | 68,592 | −69,391 | 236,946 | 344,164 | z = 7.65 |

- **Better than as built.** Every plan found is better than the plan as built given the counters, with probability Φ(z) > 1 − 2 × 10⁻¹¹.
- **How large the variances are.** The common deviation dominates the week's variance (41,043 of 42,306 for the plan as built). It is multiplicative, so a plan's s.d. grows with its value: 16% of E V for the best plan, 27% for the plan as built.
- **Each plan's own free energy** F_T = E V − Var V/2T is the certainty equivalent of the plan's holder, at risk tolerance T. It ranks the plan as built above the best plan when T < A$ 18,556 a week. At that crossover the best plan's extra variance, over twice its extra expected value, equals T; below it a planner prefers the smaller, surer value.
- **The free energy of a difference** −T ln E e^{−(V_A − V_B)/T} = ΔE − Var(V_A − V_B)/2T compares two plans under the one uncertain world they share. For the best plan against the plan as built it is +81,661 at T = 10⁴, +382,275 at 3 × 10⁴ and +487,490 at 10⁵: positive at every T shown, so the paired comparison never prefers the plan as built.
- **The two are different questions.** The first asks which plan one would rather hold; the second, whether one plan is better than another in the same week. The Gaussian theory answers both exactly.

### 7.8 The Search: What the Action's Structure Bought

At equal budgets of 20,000 evaluations over seeds 1–5 (search-86):

| search | best value, mean ± s.d. (A$ a week) | evaluations | wall time a run |
|---|---|---|---|
| **the heat bath, from as built** | **691,866 ± 104** | **3,461** (equivalents) | 3.1 s |
| **the heat bath, from a random plan** | **690,743 ± 2,627** | **3,416** (equivalents) | 3.1 s |
| UMDA | 578,153 ± 39,658 | 20,000 | 5.6 s |
| simulated annealing, from as built | 522,332 ± 30,842 | 20,000 | 19.6 s |
| CMA-ES (separable, restarts) | 441,839 ± 19,510 | 20,000 | 303 s |
| GA from as built | 371,826 ± 7,675 | 20,000 | 5.8 s |
| GA | 365,574 ± 26,916 | 20,000 | 6.2 s |
| random search | −43,051 ± 4,185 | 20,000 | 6.1 s |

- **What the structure bought.** Every heat-bath run passes the best value of all 30 black-box runs (UMDA's A$ 625,988) after 761 to 1,243 evaluations' worth of work, 4–6% of the budget.
  - It ends A$ 114,000 a week above the black boxes' best mean, A$ 66,000 above their single best run, and with a spread over seeds of A$ 104 rather than A$ 40,000, at a sixth of the budget.
  - From a random plan it reaches the same value (690,743 ± 2,627) as from the plan as built: the starting point hardly matters.
- **Where the gain comes from.** It is the action's locality, not its physics.
  - Each frontage's conditional law is computed exactly: every one of its 1.8 million states, from 5 × 168 hourly margins.
  - A black box sees one plan per evaluation, and its kernels move a frontage's kind or one hour gene at a time.
  - The heat bath's interval scan sets a frontage's kind and its weekday and weekend hours together, at their best given the neighbours. That is a move the mixed encoding's kernels make only by chance.
- **The physics' share.** The field's linear response R is what makes the attraction part of the local conditional exact: Ā = RᵀA, the value's response to each frontage's coupling. It is worth 7% of the best plan's value (§7.4).
- **The comparison with study (1).** The black-box set, settings and seeds are study (1)'s, but run on this study's passers. Its values are not comparable with study (1)'s, whose people are another model. What carries over is the ranking of methods on one value.
- **The accounting.** The heat bath's evaluations are equivalents: its local work's CPU time over the measured 0.83 ms of one black-box evaluation, with each full rebuild counted as one. The machine was shared with the other two studies, which inflates its wall time and so its equivalents. The equivalents are therefore an upper bound.

## 8. Checks, and What the Study Wrote Itself

### 8.1 Checks

`check_field` holds the study's code to closed forms, direct solves, physis's own propagator and the shared problem's `evaluate()`: 15 checks, all passing.

| Check | Largest error |
|---|---|
| the week's inverse transform undoes the forward one | 1.6 × 10⁻¹⁵ |
| the mode sum solves (iω⋆₀ + L + μ⋆₀) φ = s at six lines | 7.1 × 10⁻¹³ |
| Noether: (iω + μ) Σ⋆₀φ̂ = Σŝ at every line | 7.2 × 10⁻¹³ |
| with unit stars on the domain's cells, the mode sum is `physis::field_theory::free_propagator` | 3.6 × 10⁻¹⁵ |
| Δ₀ has one zero mode per component (6), constant on each | 2.3 × 10⁻¹³ |
| Δ₁'s harmonic dimension is E − V + b₀ − F = 12 on the subdivided lattice; Δ₂ has none | exact |
| kriging with no nugget interpolates the counters | 6.8 × 10⁻¹⁰ |
| the sources' posterior mean, carried by the stationary path, is the kriging predictor | 6.3 × 10⁻¹⁵ |
| the attraction's response R(h) equals a direct solve of (L + μ⋆₀)δφ = b L_φ̄ W δp | 1.7 × 10⁻¹⁴ |
| Sherman–Morrison for a new edge equals a direct solve | 7.5 × 10⁻¹⁴ |
| Woodbury for two of the plan's links equals a direct solve with their Whitney edges | 2.3 × 10⁻¹³ |
| the plan field's value equals `evaluate()` on the field's passers (the plan as built, 12 random plans) | 5.6 × 10⁻¹⁵ |
| after heat-bath updates (at T = 500 and T = 0), the incremental value equals a rebuild | 5.4 × 10⁻¹⁴ |
| … and equals `evaluate()` | 8.4 × 10⁻¹⁴ |
| an update at T = 0 never lowers the value | exact |

The run checks the heat bath's best plan against `evaluate()` again (A$ 691,993.5 both).

### 8.2 Gaps for the Foundations, and What This Study Wrote

The study writes no general method that a foundation already holds. It uses physis's `hodge_laplacian` for the modes and Hodge decomposition, morphe's chain complexes, calx's Cholesky, solves and dense algebra, heuresis's GA, CMA-ES, SA, UMDA and random search, and the problem's `evaluate()`. It wrote these itself, each marked provisional in its header:

- **A scalar field theory with Hodge stars** (`field.hpp`): the propagator (iω⋆₀ + L + μ⋆₀)⁻¹ as a mode sum, with Whitney's interpolation at points. physis's `field_theory::ScalarTheory` has the unweighted graph Laplacian in lattice units (½|dφ|² + m²φ²/2), which cannot carry a metric network's lengths and conductances. *For physis `field_theory`*: a scalar theory with diagonal stars, its propagator and its heat bath.
- **Gaussian conditioning of complex lines, with a trend** (`field.hpp`, `analysis.hpp`). This is universal kriging line by line: generalised least squares for the groups, the residual's posterior, REML for its scale, and the posterior covariance at new points. tyche's `spatial` conditions real fields one observation at a time, without a trend. *For tyche `spatial`*: conditioning with a regression trend (universal kriging), on complex (Fourier) data.
- **The subdivision of a cell complex** (`field.hpp`, `graph`). Every segment is cut into pieces, with the chain complex (blocks included) carried along and validated. The domain's `subdivide` rewrites one edge at a time through morphe's rewriting. *For morphe `topology`*: a barycentric-style subdivision of 1-cells that returns the subdivided chain complex.
- **A solver for stationary paths on a closed time.** It was not needed here, because the free theory is diagonal in the modes. physis's mechanics steps second-order discrete Euler–Lagrange equations forward from initial data, and has no periodic boundary-value solver for a first-order (Onsager–Machlup) action. *For physis `mechanics`*: Newton's method on calx's jets for the stationary path of a discrete action on a closed time; an interacting theory of the people would need it.
- **A heat bath with exact conditionals** (`plan.hpp`, `searches.hpp`). heuresis's annealing draws proposals from symmetric kernels and accepts them by Metropolis' rule. A problem that can compute its exact conditional laws (here, from the action's locality) has no seam to supply them. *For heuresis `annealing`*: a Gibbs (heat-bath) sweep whose conditionals the problem supplies, as an ask–tell algorithm.
- **The value's incremental evaluation** (`plan.hpp`). It is a second implementation of the problem's value, held equal to `evaluate()`. *For the shared problem*: an incremental evaluation of one frontage's change would let every study run local searches without writing its own.

### 8.3 Bugs and Oddities Found

- **The domain's cells under-resolve a field.** The district's segments are up to 224 m long; the field's screening length is 45–65 m. On the domain's own cells, a link's rebuilt complex, which splits one segment, changed the frontages' passers by up to 13%. The cause was not the link but the field's sag between corners, which the coarse lattice hid behind its linear interpolation. The study therefore subdivides to 15 m and represents links on the subdivided lattice (§2.1, §7.6). Any study of a field on the domain's district needs the same care.
- **The district is disconnected.** Its network has six components: the main one and five lane pieces that no segment joins to it. This may be worth a look in the designing study's preparation or the city model.
- **Separable CMA-ES spends about 300 s of its own** on 20,000 evaluations in 621 dimensions, against 6 s for the GA. Study (1) runs it too.
- **The shared problem's value, searched alone, rewards a retail monoculture** (§7.5).

## 9. Figures

`figures.py` draws ten figures from results/ (and the uncommitted geometry) in the house style, 1600 × 1000, into assets/gallery/.

| Figure | What it shows |
|---|---|
| `search-80-the-people-as-a-field` | the field's weekly mean on the network, the counters, the frontages' feet; leave-one-site-out against the nearest counter; the calibration's numbers |
| `search-81-the-districts-modes` | Δ₀'s first six modes above its kernel; Δ₁'s spectrum split into gradients, curls and harmonic currents |
| `search-82-space-time-and-frequency` | the field's power by spatial mode and line; the action's filter in time; the stationary path at four frontages; its lines |
| `search-83-people-flow` | Noether's charge N(t) and the inflow; Fick's current at Saturday 13:00; the counters' directions against it |
| `search-84-calibration` | the 17 counters' mean weeks, the field's left-out predictions, the nearest counter's; the grid |
| `search-85-the-euclidean-measure` | the common deviation's s.d. by hour; the residuals' correlation against the action's; the posterior s.d. at the frontages |
| `search-86-searches-at-equal-budgets` | anytime curves and final values of the heat bath and the six black boxes, five seeds each |
| `search-87-the-plan-and-its-relations` | the kinds as built and in the best plan; the catalysis, competition and party walls the plan uses |
| `search-88-the-week-of-the-plan` | spending, frontages open and passers drawn by hour; the spending's lines |
| `search-89-links-and-free-energies` | each link's effect by first-order response, exactly on the lattice and on the domain's rebuilt complex; free energies against T |

## 10. Caveats

- **This is not a planning proposal.** The economics is the shared problem's sketch (its §11), and its best plan, a retail monoculture, is a property of that sketch.
- **The people are a diffusive fluid.**
  - The action knows densities and their gradients, not purposes. It reproduces where people are (with a held-out error of a factor 1.6 in weekly totals), not where they are going: its current is uncorrelated with the counters' directions (§7.3).
  - Swanston Street, the grid's spine, is under-predicted when left out. Study (2)'s conserved flows are the natural complement.
- **Seventeen counters on 13 km of streets.** The posterior s.d. of a frontage's weekly passers is a median 27% of its mean. The plans' values carry that uncertainty (§7.7), and the plans are compared under it, not as if the passers were known.
- **The lattice.** The field is a lattice theory on the district's segments subdivided to 15 m. A finer subdivision, or the continuum metric graph, would change ξ's fit and the passers slightly. The domain's cells alone would not do (§8.3).
- **The plan's attraction.** b = 5 is an assumption; the best plan's value at b = 0 and 20 is reported (§7.4). The attraction is in linear response, quasi-static, and computed on the network as built.
- **The links.** A link is a massless edge on the lattice. The domain's rebuilt complex, in which the link is a street with its own leaving, drains more people (§7.6). Neither gives an arcade visitors of its own.
- **The fluctuations.**
  - They are calibrated to the counters' residuals: a common deviation and a local field with the action's spatial correlation.
  - The action's correlation length is far shorter than the counters' spacing, so the local field's spatial structure is barely tested.
  - The noise model, white per unit length with a spectrum in time, is an assumption of form.
- **The heat bath is local.** Annealed from A$ 3,000 to A$ 3 over 40 sweeps, it reached A$ 691,866 ± 104 from the plan as built over five seeds. At b = 0 one run fell A$ 10,850 short of a better plan. No search here is shown to be optimal.
- **The equal budget.** The heat bath's evaluations are equivalents of CPU time on a machine shared with two other studies, so they are an upper bound.

## 11. Running It

From the repository's root, after the shared problem's results exist (its `prepare.py` and `calibrate`):

```
.venv/Scripts/python.exe research/experiments/2026-11-05-a-district-economy-as-a-field-theory/counters.py
cmake -S research/experiments/2026-11-05-a-district-economy-as-a-field-theory -B ~/.cache/science-build/exp-economy-field -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-economy-field -j4
R=research/experiments/2026-11-05-a-district-economy-the-problem/results; X=research/experiments/2026-11-05-a-district-economy-as-a-field-theory
~/.cache/science-build/exp-economy-field/check_field $R/problem-inputs.json $R/problem.json $R/sensors.json $X/results/counters.json
~/.cache/science-build/exp-economy-field/field_theory $R/problem-inputs.json $R/problem.json $R/sensors.json $X/results/counters.json $X/results
~/.cache/science-build/exp-economy-field/field_theory $R/problem-inputs.json $R/problem.json $R/sensors.json $X/results/counters.json $X/results links
.venv/Scripts/python.exe $X/figures.py
```

`field_theory` takes a stage as its last argument:
- `calibrate`, the field alone;
- `field`, the field and its modes;
- `quick`, one seed at a small budget, a test of the pipeline;
- `links`, the links' comparison at a finished run's best plan, written back into search.json;
- `all`, the default.

The full run takes about 13 minutes, and `links` about 3.

## 12. References

- Benamou, J.-D. and Brenier, Y. (2000), A computational fluid mechanics solution to the Monge–Kantorovich mass transfer problem, *Numer. Math.* 84, 375–393.
- Bertini, L., De Sole, A., Gabrielli, D., Jona-Lasinio, G. and Landim, C. (2015), Macroscopic fluctuation theory, *Rev. Mod. Phys.* 87, 593–636.
- Chen, Y., Georgiou, T. T. and Pavon, M. (2016), On the relation between optimal transport and Schrödinger bridges: a stochastic control viewpoint, *J. Optim. Theory Appl.* 169, 671–691.
- Eckmann, B. (1944), Harmonische Funktionen und Randwertaufgaben in einem Komplex, *Comment. Math. Helv.* 17, 240–255.
- Fruin, J. J. (1971), *Pedestrian Planning and Design*, Metropolitan Association of Urban Designers and Environmental Planners.
- Léonard, C. (2014), A survey of the Schrödinger problem and some of its connections with optimal transport, *Discrete Contin. Dyn. Syst. A* 34, 1533–1574.
- Lindgren, F., Rue, H. and Lindström, J. (2011), An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach, *J. R. Stat. Soc. B* 73, 423–498.
- Maas, J. (2011), Gradient flows of the entropy for finite Markov chains, *J. Funct. Anal.* 261, 2250–2292.
- Martin, P. C., Siggia, E. D. and Rose, H. A. (1973), Statistical dynamics of classical systems, *Phys. Rev. A* 8, 423–437.
- Onsager, L. and Machlup, S. (1953), Fluctuations and irreversible processes, *Phys. Rev.* 91, 1505–1512.
- Weidmann, U. (1993), *Transporttechnik der Fussgänger*, Schriftenreihe des IVT 90, ETH Zürich.
- Woodbury, M. A. (1950), *Inverting Modified Matrices*, Memorandum Report 42, Statistical Research Group, Princeton University.
