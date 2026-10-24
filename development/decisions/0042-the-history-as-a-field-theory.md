# 0042. The History as a Field Theory, on the Foundations

- **Status:** Accepted, provisionally (the estimators of the process stay in the domain until the inference foundation exists); amends 0041
- **Date:** 2026-10-24
- **Decided by:** the user asked three things. First, whether the spins and amplitudes of a history's spin foam have a meaning and can be calculated. Second, whether one can pass back from the history's spacetime to a field theory's formalism, its integral, its perturbative expansion and its action S, perhaps by fitting a model. Third, that it be made "more accurate and physically and logically sound". When they saw the domain growing its own transfer matrix and sampler, they set a principle: what the theoretical libraries already hold is reused, not written again, because "reusable, integrative science foundation" is the purpose of the repository. Claude designed the computation, the field theory and the moves into the foundations, and wrote this record. The foundations' parts are physis [decision 0007](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0007-exact-z2-contraction-and-finite-heat-bath.md) and morphe's cycles, Q-analysis and canonical forms.

## 1. Context

ADR 0041 derived every label of the social spin foam from the character expansion of a ℤ₂ Gibbs weight on the refined complex K̂. Its numbers had three weaknesses:

- **The couplings were moment estimates**, each made as if its cells were independent. Two lenses of one tie share a rung, and the triangles at a congress share their rungs pairwise, so the same evidence was counted several times.
- **The spins, Wilson loops and posterior were sampled**, by a heat bath of the domain's own with a cluster move added so that strong balance would not freeze.
- **The field theory went one way only.** The history was read as a gauge field on a given complex. The passage back, from the complex to a field theory whose expansion generates it, was not made.

As the work began, the domain also gathered methods that name no history: a transfer matrix over ℤ₂, a sampler, dimensions of SU(2) invariants, cycles over GF(2), Atkin's Q-analysis, canonical forms of labelled graphs. The user's principle says where these belong.

## 2. Decision

### 2.1 The Spin Foam Computed Exactly

Z, the means of every rung and cell in the prior and the posterior, two-point functions, the laws at cuts, and the spin networks along which the foams cross a cut are computed exactly by physis's contraction of ℤ₂ gauge theories (theory note §6.6). The domain supplies the theory (links the rungs, 2-cells the lenses and triangles with their couplings, the attested signs fixed in the posterior) and the order of time. In that order the prior's states are confined to the ℤ/2 1-cycles of the network of the moment (Proposition 6.3), so the cost is set by the number of independent circles of open ties.

Every single canon is in reach: the largest state has 131 072 characters, and a build with the fit takes 26 s at most. The six canons joined are beyond reach and are sampled, by physis's heat bath of a finite group on the whole of K̂. That theory is the domain's model in temporal gauge, times 2^{#temporal} (Proposition 6.5), which the checks confirm on every canon in reach to 10⁻¹².

### 2.2 The Couplings Fitted

The couplings are the maximum a posteriori of the exact likelihood of the attested signs, log Z_post − log Z_prior. Each coupling's rate of change has a Jeffreys prior. The gradient and Hessian are exact, from calx's jets carried through the contraction; Newton's method runs from the moment estimates, and Laplace's approximation gives the errors. With independent cells the maximum is the moment estimate itself (Proposition 6.4). On the canons the correction is large where cells share rungs: balance falls from 0.83 to 0.64 in Europe 1785–1828 and from 0.95 to 0.60 in the alliances of 1871–1914. The application shows each coupling beside its moment estimate, with its interval.

### 2.3 The History as a Feynman Graph

A new module, `history::field`, reads the history as one Feynman graph of a group field theory over ℤ₂ (theory note §8):

- a segment of a worldline is a field quantum of rank n, n being the ties alive across it;
- an occasion is an interaction vertex, whose bubble, the link of the vertex, pairs the strands at closing, opening and passing corners;
- the bubbles are typed by morphe's canonical forms, with their symmetry factors.

The Feynman rules (δ's along propagators and at passing and opening corners, Gibbs kernels at closing corners, balance at the vertices) are eliminated over GF(2) for every record. Proposition 8.1 states that the amplitude of the observed graph is 2^{F + E} Z, each kernel taking its lens's holonomy on K̂, and the domain checks it on every canon. Nearly every occasion is a bubble type of its own. A history's action therefore has nearly as many interaction terms as occasions, and their couplings cannot be identified from one history.

### 2.4 The History in Time

Time-ordered, a history is one term of a Dyson series, and the same module estimates the laws of its occasions and strands as a marked point process:

- **occasions:** a rate λ_k per k-subset of continuing beings, with new beings created at the vertex;
- **strands:** the probabilities q_open and q_close that a tie begins or ends at a meeting;
- **binding:** ρ = μ₁/μ₀, the factor by which a live tie raises its pair's rate of meeting.

The theory's martingales, Dynkin's formula as the stochastic Schwinger–Dyson equations, test the fitted rates era by era. Constant rates are rejected for four of the six canons. The caveat is stated with the result: one canon cannot tell a history's tempo from its record's coverage.

### 2.5 Reuse of the Foundations

The domain keeps what is the record's, and the general methods live in the foundations:

| Method | Now in |
|---|---|
| exact contraction of ℤ₂ gauge theories, with derivatives | physis `gauge/contraction.hpp` |
| couplings by 2-cell, frozen links, the heat bath of finite groups | physis `gauge` (`field.hpp`, `updates.hpp`) |
| dimensions of SU(2) invariants | physis `symmetry/recoupling.hpp` |
| cycles over ℤ/2 (the space of foams) | morphe `topology/cycles.hpp` |
| Atkin's Q-analysis | morphe `topology/q_analysis.hpp` |
| canonical forms and automorphism groups (bubble types) | morphe `combinatorics/canonical.hpp` |
| binomial coefficients, linear solves, jets, binned errors | morphe `arithmetic`, calx `linalg` and `differentiation`, physis `field_theory` |

What stays in the domain is what is the record's:

- which cells there are, with which couplings, and which links the record fixes;
- the order of time and the cuts at dates;
- the loop around a tie's life;
- the elimination of the history's Feynman rules;
- the estimators of the process, until the inference foundation exists (README §4).

The domain's `module.json` declares physis among its foundations.

### 2.6 The Application

The History Complex gains a Field view. It shows the bubbles drawn as interaction vertices, the action term by term, the amplitude and the exact log-likelihood, the rates of the process with their martingale checks, and the χ² by era. The Moment view lists the spin networks at its cut, with their exact weights. Every answer says whether it is exact. A build's fitted couplings are passed to the later requests, and builds are kept in a cache, so nothing is fitted twice.

## 3. Alternatives Considered

- **Keep sampling.** The contraction reaches every single canon in seconds, and exact numbers can be checked against brute force and closed forms; sampled ones only within errors.
- **Fit a coupling for each bubble type.** Nearly every type occurs once, so such couplings are not identifiable from one history. The rates are grouped by the number of continuing beings instead.
- **A hierarchical prior over the colours' couplings.** It would pool evidence between colours. It belongs to the inference foundation, which is planned. Until then a colour with no evidence stays at the prior's mode, β = 0.
- **Keep the general methods in the domain.** Rejected by the user's principle: another 2-complex, another relation or another graph would have had to copy them.

## 4. Consequences

- **The domain** gains `foam/transfer.hpp`, an adapter of the record to physis's contraction together with the fit, and the module `field`. `foam/model.hpp` samples by physis and takes its space of foams from morphe. Its checks number 351, among them the model against physis's theory on all of K̂ on every canon in reach.
- **The foundations:** physis decision 0007; morphe's cycles over ℤ/2, Q-analysis and canonical forms.
- **The application's engine** compiles against physis as well.
- **CLAUDE.md** records the principle of reuse in §2.1, and §9 the state.
- **ADR 0041** is amended: the couplings are fitted exactly, the ensembles are exact where in reach, and the domain uses physis.
