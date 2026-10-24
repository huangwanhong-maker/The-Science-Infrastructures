# Hidden Markov Models and Their Semirings

*Status: draft. Serves learning and use.* Hidden Markov models as sums over the paths of a trellis; the three problems of an HMM as one recurrence read in two semirings; scaling; Baum–Welch as expectation maximisation and why it never goes down; the CpG-island model; and profile HMMs of a family of sequences. The code is the module [`hmm`](../hmm/): [`model.hpp`](../hmm/include/genomics/hmm/model.hpp), [`cpg.hpp`](../hmm/include/genomics/hmm/cpg.hpp), [`profile.hpp`](../hmm/include/genomics/hmm/profile.hpp). The pair HMM of two sequences is the align module's ([alignment.md](alignment.md)). The design is the idea note [The Genome Infrastructure](../../../development/ideas/2026-10-24-the-genome-infrastructure.md), §4.2.

## 1. Hidden Markov Models

**Definition 1.1 (HMM).** A hidden Markov model on a finite alphabet Σ (|Σ| = K) is a tuple (Q, π, A, E, ω) of:

- a finite set Q of S states;
- an initial distribution π on Q (the transitions a₀ₖ out of a silent begin state);
- a transition matrix A = (a_kl) with non-negative entries;
- an emission matrix E = (e_k(b)), each row a distribution on Σ;
- optionally an end vector ω = (a_k0), with Σ_l a_kl + a_k0 = 1 for every k; without one, Σ_l a_kl = 1.

A *path* is a word π = π₁ ⋯ π_T over Q, and the joint probability of a path and a sequence x = x₁ ⋯ x_T is

P(x, π) = a_{0π₁} e_{π₁}(x₁) · Π_{t=2}^{T} a_{π_{t−1}π_t} e_{π_t}(x_t) · a_{π_T 0},

the last factor present only with an end (Durbin et al. 1998, §3.2). In the code (`Model`) the fields are `initial`, `transition`, `emission` and `end`.

*Remark.* With an end state, Σ_T Σ_{x ∈ Σ^T} P(x) = 1 when the end is reached with probability one: the model is a distribution over all finite sequences, and the length is part of what it explains. Without one, it is a distribution over Σ^T for each T (Rabiner's convention). Both are needed: a profile HMM (§6) models whole sequences with their lengths; a CpG model (§5) reads a chromosome of given length.

**Definition 1.2 (missing symbol).** A code c ≥ K is *missing*: e_k(c) = 1 for every k. Then P(x) for an x with a missing position is the sum of P(x′) over the K completions x′, by distributivity, which is how N in DNA is read.

**The three problems** (Rabiner 1989, §III):

1. the *evaluation* of P(x) = Σ_π P(x, π);
2. the *decoding* of π* = argmax_π P(x, π), or of the marginals P(π_t = k | x);
3. the *learning* of (π, A, E, ω) from sequences.

## 2. One Sum in Two Semirings

**Proposition 2.1 (the matrix form).** Let D(b) = diag(e_1(b), …, e_S(b)). Then

P(x) = πᵀ D(x₁) A D(x₂) A ⋯ A D(x_T) ω,

with ω = (1, …, 1)ᵀ when the end is not modelled.

*Proof.* Expanding the product of matrices sums, over every choice of the indices π₁, …, π_T, the product of the chosen entries, which is P(x, π). ∎

The forward variables are the left partial products, α_t = πᵀ D(x₁) A ⋯ A D(x_t), so α_t(k) = P(x₁ ⋯ x_t, π_t = k), and the recursion α_t(l) = e_l(x_t) Σ_k α_{t−1}(k) a_kl is the product evaluated from the left (Rabiner, eqs. 19–21). The backward variables are the right partial products, β_t(k) = P(x_{t+1} ⋯ x_T, end | π_t = k), evaluated from the right (eqs. 23–24). That the two orders give one number is the associativity of the matrix product; the check `check_genomics_hmm_model` holds them to each other and to the enumeration of all S^T paths.

**Definition 2.2 (semirings of log weights).** On ℝ̄ = ℝ ∪ {−∞}:

- the *log semiring* (ℝ̄, ⊕, +) with a ⊕ b = log(e^a + e^b), zero −∞, one 0;
- the *tropical semiring* (ℝ̄, max, +), zero −∞, one 0.

**Proposition 2.3.** log is an isomorphism of semirings from (ℝ≥0, +, ×) onto the log semiring; and for h > 0, the operations a ⊕_h b = h log(e^{a/h} + e^{b/h}) converge to max(a, b) as h → 0, uniformly: max(a, b) ≤ a ⊕_h b ≤ max(a, b) + h log 2.

*Proof.* The first is the functional equation of the logarithm. For the second, e^{max/h} ≤ e^{a/h} + e^{b/h} ≤ 2 e^{max/h}; take h log. ∎

So the log-probability of the best path is the tropicalisation of the likelihood: the same product of matrices in the semiring (max, +) applied to the logarithms of the entries (Maslov's dequantisation; Pachter and Sturmfels 2005, §2.2, who read the likelihood of an HMM as a polynomial in its parameters and Viterbi as its tropical counterpart). `path_sum<R>` is that one recurrence, generic in R: with `LogSum` it returns log P(x), with `Tropical` the Viterbi score. **Viterbi's algorithm** (Viterbi 1967) is the tropical forward recursion δ_t(l) = log e_l(x_t) + max_k (δ_{t−1}(k) + log a_kl) together with the arguments of the maxima, from which the best path is read backwards (Rabiner, eqs. 32–35); the code breaks ties toward the smallest k, so the path is a function of the model and x.

*Remark (posterior decoding).* The most probable path and the path of most probable states differ: π̂_t = argmax_k P(π_t = k | x) may not even be a path of positive probability. They answer different losses (Durbin et al. §3.2): Viterbi minimises the expected 0–1 loss on whole paths, posterior decoding the expected number of wrong positions.

## 3. Scaling

A product of T probabilities underflows a double for T in the thousands (the check's 20 000 symbols have P(x) ≈ e^{−29461}). Two cures exist, logarithms and scaling; the code uses logarithms for Viterbi (sums are exact and fast) and scaling for forward and backward (sums of probabilities need no exponentials).

**Proposition 3.1 (Rabiner 1989, §V.A).** Define α̂₁ = α₁/c₁ and α̂_t = α̃_t / c_t with α̃_t(l) = e_l(x_t) Σ_k α̂_{t−1}(k) a_kl and c_t = Σ_l α̃_t(l). Then α_t = α̂_t Π_{s≤t} c_s, and

log P(x) = Σ_t log c_t + log Σ_k α̂_T(k) a_k0.

*Proof.* By induction, α̃_t = α_t / Π_{s<t} c_s, since the recursion is linear; dividing by c_t gives the claim for α̂_t, and P(x) = Σ_k α_T(k) a_k0. ∎

Each α̂_t sums to one, so nothing underflows that matters. The backward variables are scaled by constants d_t of their own, so that `backward()` computes log P(x) independently of `forward()`, and the check compares them.

**Proposition 3.2 (posteriors from any scaling).** γ_t(k) = P(π_t = k | x) = α_t(k) β_t(k) / P(x) (Durbin et al. §3.2), and γ_t is proportional to α̂_t ⊙ β̂_t for any positive scalings of α_t and β_t. Likewise ξ_t(k, l) = P(π_t = k, π_{t+1} = l | x) is proportional to α̂_t(k) a_kl e_l(x_{t+1}) β̂_{t+1}(l).

*Proof.* Σ_k γ_t(k) = 1, so γ_t is determined by any vector proportional to it; and α_t(k) β_t(k) = P(x, π_t = k). The same for ξ_t, which sums to one over (k, l). ∎

This is why `posterior()` normalises the product at each t instead of tracking the scales.

## 4. Baum–Welch as Expectation Maximisation

Let θ be the parameters and x⁽¹⁾, …, x⁽ⁿ⁾ independent sequences. The paths are the hidden data.

**Definition 4.1 (expected counts).** Under θ,

- I_k = Σ_j γ₁⁽ʲ⁾(k), the expected number of paths starting in k;
- A_kl = Σ_j Σ_t ξ_t⁽ʲ⁾(k, l), the expected number of transitions k → l;
- E_k(b) = Σ_j Σ_{t : x_t⁽ʲ⁾ = b} γ_t⁽ʲ⁾(k), the expected number of emissions of b from k;
- F_k = Σ_j γ_T⁽ʲ⁾(k), the expected number of endings from k.

**Proposition 4.2 (the M-step; Rabiner eqs. 40a–c and 109–110; Durbin et al. §3.3).** The maximiser over θ′ of Q(θ′ | θ) = E_θ[log P(x, π | θ′) | x] is

π′_k = I_k / n, a′_kl = A_kl / (Σ_l′ A_kl′ + F_k), a′_k0 = F_k / (Σ_l′ A_kl′ + F_k), e′_k(b) = E_k(b) / Σ_b′ E_k(b′).

*Proof.* log P(x, π | θ′) is a sum of terms log π′_k, log a′_kl, log e′_k(b), each counted as often as the path uses it, so Q is Σ I_k log π′_k + Σ A_kl log a′_kl + Σ F_k log a′_k0 + Σ E_k(b) log e′_k(b). Its parts separate by the rows of the distributions, and on a simplex Σ c_i log p_i is maximised at p_i = c_i/Σ c (Gibbs' inequality, or a Lagrange multiplier). ∎

**Theorem 4.3 (EM never decreases the likelihood; Baum et al. 1970, Dempster, Laird and Rubin 1977).** If θ′ maximises Q(· | θ), then P(x | θ′) ≥ P(x | θ).

*Proof.* For any θ′, log P(x | θ′) = Q(θ′ | θ) − H(θ′ | θ), where H(θ′ | θ) = E_θ[log P(π | x, θ′) | x]. By Gibbs' inequality H(θ′ | θ) ≤ H(θ | θ). Hence log P(x | θ′) − log P(x | θ) ≥ Q(θ′ | θ) − Q(θ | θ) ≥ 0. ∎

**Proposition 4.4 (pseudocounts are a prior).** Adding r to each expected count of a row maximises Q(θ′ | θ) + Σ r log θ′, the expected complete log posterior under a Dirichlet prior with parameters r + 1. Theorem 4.3 then holds for the log posterior log P(x | θ) + Σ r log θ, not for the likelihood.

*Proof.* The prior adds r log p_i to each term of the row, and the maximiser of Σ (c_i + r) log p_i is (c_i + r)/Σ (c + r). The proof of 4.3 goes through with log P(x | θ) replaced by the log posterior, since the prior does not depend on π. ∎

The code reports both quantities (`Training::log_likelihood`, `Training::log_posterior`) and the checks verify the right one is monotone in each case: with pseudocounts the likelihood itself may fall slightly, as the prior pulls the parameters off the maximum of the likelihood. A parameter the model sets to zero has zero expected count and receives no pseudocount, so it stays zero: a left-to-right model, or a profile, keeps its topology.

*Remark (what EM does not promise).* It converges to a stationary point, which may be a local maximum, and the states may come out permuted (label switching). The check therefore starts from a model halfway between the truth and uniform, and verifies a property that holds at the global maximum and is testable: the fitted likelihood is at least that of the true parameters (on 60 000 symbols, −78 708.4 against −78 715.7 nats). The parameters then lie within five standard errors √(p(1 − p)/n_k) of the truth, n_k the expected visits of state k.

**Remark (determinism).** The expected counts are sums over sequences; the code cuts the sequences into shards of a fixed size and adds the shards' counts in order, so that the floating-point result does not depend on the number of threads (`baum_welch` with a `ThreadPool`: bit for bit on one thread and four).

## 5. CpG Islands

**Definition 5.1 (the two chains; Durbin et al. §3.1).** Two first-order Markov chains on {A, C, G, T}: a⁺ estimated from putative CpG islands of human DNA, a⁻ from the sequence between them. Their tables (p. 50) are `CPG_PLUS` and `CPG_MINUS`. The log-odds of a sequence is

S(x) = log₂ P(x | +)/P(x | −) = Σ_i β_{x_{i−1} x_i}, β_st = log₂ (a⁺_st / a⁻_st),

in bits (the book's table on p. 51). Its largest entry is β_CG = 1.812: in the − chain C is followed by G with probability 0.078, because methylated C in CG deaminates to T, while islands are unmethylated (Bird 1987).

**Definition 5.2 (the eight-state HMM; §3.2).** States s⁺ and s⁻ for each base s, each emitting its own base; within a label the transitions are those of its chain scaled by 1 − p (in) or 1 − q (out), and a(s⁺ → t⁻) = p a⁻_st, a(s⁻ → t⁺) = q a⁺_st.

**Proposition 5.3.** Under this model the label sequence is a Markov chain on {+, −} with P(+ → −) = p and P(− → +) = q, whatever the bases; so island lengths are geometric with mean 1/p, gaps with mean 1/q, and the stationary probability of + is q/(p + q). Given the labels, the bases are the chain of each label, the base after a switch drawn from the new label's row.

*Proof.* The total probability out of s⁺ into the − states is Σ_t p a⁻_st = p, and into + is (1 − p) Σ_t a⁺_st = 1 − p; the same for −. These do not depend on s, so the label is lumpable (Kemeny and Snell 1960, §6.3). ∎

The book leaves p and q to the modeller; the defaults (1 kb islands, 100 kb between) are the human genome's order of magnitude. Decoding by Viterbi or by posterior (`cpg_islands`) returns the maximal runs of + states. The check draws a sequence of 600 kb from the model itself (so the truth is the path's labels) and decodes it with base-level sensitivity 0.987 and specificity 0.988.

## 6. Profile HMMs

**Definition 6.1 (the profile; Krogh et al. 1994; Durbin et al. Figure 5.2).** For m match states: M₁, …, M_m emitting by e_{M_k}, insert states I₀, …, I_m emitting by e_{I_k}, silent deletes D₁, …, D_m; M₀ is the begin and M_{m+1} the end. Out of node k there are nine transitions, from M_k, I_k and D_k to M_{k+1}, I_k and D_{k+1}, those to D₀ or D_{m+1} absent.

A profile is a hidden Markov model whose states are laid out along the columns of a family's alignment: the match states are the conserved columns, inserts and deletes the gaps.

**Construction 6.2 (Durbin et al. §5.3).** From an alignment of n rows:

1. a column is a *match column* when at most half its rows are gaps (the 50% rule);
2. each row determines a path: a residue in the k-th match column is M_k, a gap there is D_k, a residue in an insert column after the k-th match column is I_k, a gap there is nothing;
3. the transitions and emissions along the n paths are counted, pseudocounts added, and each row of counts normalised.

*Remark.* Step 2 is where the alignment's information enters: the paths are no longer hidden, so the estimate is the complete-data maximum likelihood of §4 (with the prior of the pseudocounts), not EM. The check works Durbin et al.'s globin example (Figure 5.3) by hand: of its ten columns the 50% rule keeps eight; GLB1_GLYDI's two residues in columns 4–5 are I₃I₃, so out of I₃ the counts are (IM, II, ID) = (1, 1, 0), (2, 2, 1)/5 with Laplace's rule.

**Definition 6.3 (log-odds scores; §5.4).** Against the null model R of i.i.d. residues with frequencies q, the Viterbi score is log₂ max_π P(x, π | M) / P(x | R) and the forward score log₂ Σ_π P(x, π | M) / P(x | R). Dividing each emission by q makes every factor of P(x | R) cancel against one of the path, so the recurrences of §2 run on the odds e/q:

V^M_j(i) = log (e_{M_j}(x_i)/q_{x_i}) + max { V^M_{j−1}(i − 1) + log a_{M_{j−1}M_j}, V^I_{j−1}(i − 1) + log a_{I_{j−1}M_j}, V^D_{j−1}(i − 1) + log a_{D_{j−1}M_j} },

and likewise V^I_j(i) from row i − 1 at node j, V^D_j(i) from row i at node j − 1 (no residue).

**Proposition 6.4 (the order of evaluation).** The rows i = 0, …, L in order, and within a row the nodes j = 0, …, m in order, is a topological order of the recurrence.

*Proof.* M_j(i) depends on row i − 1; I_j(i) on row i − 1; D_j(i) on node j − 1 of row i. No dependency points to a later row, or within a row to a later node. ∎

**Modes.** Global: the whole sequence against the whole model. Glocal: the model whole, the sequence's ends free, by entering M₀ at any row and summing (or maximising) the end over all rows. Local: entry into any M_k with probability 2/(m(m + 1)) and exit from any M_k, the uniform distribution over the m(m + 1)/2 fragments [k, k′] of the model (Eddy 2008). The check enumerates every path of a four-state profile through short sequences, in each mode, and holds both scores to the maximum and the sum.

*Caveat.* The free flanks have no length distribution here, so a free-mode forward score sums over every placement with weight one; HMMER's flanking states and null model make the score a true likelihood ratio, lower by up to about log₂ L bits.

**Scaling the forward of a profile.** The rows of the profile's lattice play the part of the positions t of §3: every value of row i is computed from row i − 1 or from earlier nodes of row i, so dividing a row by its sum c_i rescales everything after it by the same factor. The free begin state, worth 1 at every row, enters row i as 1/C_{i−1}, C_{i−1} = Π_{s<i} c_s; since it is part of the row's sum, C_i ≥ 1 and the term never overflows. This forward (`profile_forward`) runs 25 times faster than the same sum in the log semiring (`profile_forward_log`), to which the check holds it on 20 kb to 10⁻¹³.

**The alignment of a sequence to a profile** is the Viterbi path: its match states give the A2M row of the sequence in the family's match columns, its inserts what lies between. On a simulated family (`simulate::evolve_family`), the alignment places 94% of held-out members' residues in the match states of their true columns; the rest fall where the true alignment is one of several equally good ones, an indel in a run of equal bases.

## 7. Where the Code Is

| Concept | Code |
|---|---|
| Definition 1.1, the joint probability | `Model`, `path_log_probability` |
| Proposition 2.1, Definition 2.2 | `path_sum<LogSum>`, `path_sum<Tropical>` |
| Viterbi | `viterbi` (back-pointers of one byte for S ≤ 256) |
| Propositions 3.1–3.2 | `forward`, `backward`, `log_likelihood`, `posterior`, `posterior_decode` |
| Theorem 4.3, Proposition 4.4 | `baum_welch`, `Training` |
| Definitions 5.1–5.2 | `CPG_PLUS`, `CPG_MINUS`, `log_odds_table`, `log_odds`, `cpg_model`, `cpg_islands` |
| Construction 6.2 | `build_profile` |
| Definition 6.3, modes | `profile_viterbi`, `profile_forward`, `profile_forward_log`, `align_to_profile`, `profile_path_score` |

## References

- Baum, L. E., Petrie, T., Soules, G. and Weiss, N. (1970), *Annals of Mathematical Statistics* 41, 164–171.
- Bird, A. P. (1987), CpG islands as gene markers in the vertebrate nucleus, *Trends in Genetics* 3, 342–347.
- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977), *Journal of the Royal Statistical Society B* 39, 1–38.
- Durbin, R., Eddy, S. R., Krogh, A. and Mitchison, G. (1998), *Biological Sequence Analysis*, Cambridge University Press, ch. 3 and 5.
- Eddy, S. R. (1998), Profile hidden Markov models, *Bioinformatics* 14, 755–763; (2008), *PLoS Computational Biology* 4, e1000069.
- Kemeny, J. G. and Snell, J. L. (1960), *Finite Markov Chains*, Van Nostrand.
- Krogh, A., Brown, M., Mian, I. S., Sjölander, K. and Haussler, D. (1994), *Journal of Molecular Biology* 235, 1501–1531.
- Pachter, L. and Sturmfels, B. (eds.) (2005), *Algebraic Statistics for Computational Biology*, Cambridge University Press.
- Rabiner, L. R. (1989), *Proceedings of the IEEE* 77, 257–286.
- Viterbi, A. J. (1967), Error bounds for convolutional codes and an asymptotically optimum decoding algorithm, *IEEE Transactions on Information Theory* 13, 260–269.
