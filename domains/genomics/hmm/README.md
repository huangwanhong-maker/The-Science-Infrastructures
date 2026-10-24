# hmm: Hidden Markov Models, CpG Islands and Profiles

*Status: draft (2026-10-24).* Hidden Markov models on a finite alphabet and their three problems (the probability of a sequence, its most probable path, the estimation of the parameters), and the two models of biological sequence analysis built on them here: the CpG-island model and profile HMMs of a family of sequences. Serves use and learning; the theory is [theory/hmm.md](../theory/hmm.md).

## 1. Concepts

- **One recurrence, two semirings.** The probability of a sequence is a sum over the paths of the trellis of products of weights; the probability of the best path is the same expression with the sum replaced by a maximum. `path_sum<R>` runs the recurrence in the log semiring (log-sum-exp, +) or the tropical semiring (max, +); the fast paths (Viterbi with back-pointers, forward and backward scaled as Rabiner's §V.A) are held to it and to the enumeration of every path.
- **EM never goes down.** Baum–Welch is expectation maximisation; with pseudocounts it is EM for the maximum a posteriori estimate under Dirichlet priors, so the quantity that never decreases is the log posterior, which it reports beside the likelihood. A zero parameter stays zero, so a model's topology (left to right, a profile's) is kept.
- **Determinism.** Sampling draws from `calx::random::Stream` (purpose `SAMPLE`); Baum–Welch's E-step runs on a `ThreadPool` in shards of a fixed number of sequences, merged in order, so its result is the same on any number of threads, bit for bit.
- **Missing data.** A code at or above the number of symbols (N in DNA) is emitted with probability 1 by every state: the sum over what it might have been.

*Remark.* The generic machinery is statistics, not genomics: it names no object in nature, and moves to the planned inference foundation when that exists (CLAUDE.md §8, open decision 3). It is written for any number of states and symbols, its only genomic dependency being the `Code` type, so that the move is mechanical.

## 2. Content

| Header | Content |
|---|---|
| [`model.hpp`](include/genomics/hmm/model.hpp) | `Model` (S states, K symbols, initial, transition, emission, optional end); `path_sum<Tropical \| LogSum>`; `viterbi`, `path_log_probability`; `forward`, `backward`, `log_likelihood` (scaled); `posterior`, `posterior_decode`; `baum_welch` with pseudocounts, partial re-estimation and a thread pool; `sample`, `sample_to_end` |
| [`cpg.hpp`](include/genomics/hmm/cpg.hpp) | Durbin et al.'s + and − chains (p. 50), their log-odds table β, `log_odds` of a window and in windows, `estimate_chain`; the eight-state `cpg_model`; `cpg_islands` by Viterbi or posterior decoding, each island with its log-odds, GC and CpG o/e |
| [`profile.hpp`](include/genomics/hmm/profile.hpp) | `Profile` (match, insert, delete states; nine transitions a node); `build_profile` from an alignment by the 50% rule with Laplace or background pseudocounts and optional weights; `profile_viterbi`, `profile_forward` (scaled), `profile_forward_log` (the reference) in global, glocal and local modes, in bits; `align_to_profile` with the A2M row and inserts; `profile_path_score` |

## 3. Dependencies

`genomics::sequence` (codes); calx (`random` for sampling, `scheduler` for Baum–Welch in parallel). The checks also use `genomics::simulate` (families with true alignments).

## 4. Checks

`check_genomics_hmm_model`, `check_genomics_hmm_cpg`, `check_genomics_hmm_profile`:

| Check | What it verifies |
|---|---|
| model | Viterbi, forward, backward and posteriors against the enumeration of every path of 3-state models (with and without an end, with structural zeros, with a missing symbol); scaling against the log semiring on 20 000 symbols; Baum–Welch: the likelihood (no pseudocounts) and the log posterior (pseudocounts) never decrease, the optimum at least as likely as the truth, parameters within 5 standard errors, zeros kept, one thread equal to four; the end probability recovered; sample frequencies |
| cpg | The chains' rows; the log-odds equal to the book's table β (p. 51) to 0.006 bits; the eight-state model valid with the sojourns asked; a chain re-estimated from its samples; islands drawn from the model found in 600 kb (base-level Sn 0.987, Sp 0.988; 32 of 34 islands) by Viterbi and by posterior decoding |
| profile | Durbin et al.'s globin alignment (Figure 5.3): eight match columns, emissions and transitions by Laplace's rule worked by hand; Viterbi and forward (scaled and log) against every path of a small profile in three modes; alignments rescoring to their scores; a simulated family: held-out members above shuffled sequences by more than 40 bits, 94% of residues aligned to their true match states; the scaled forward equal to the log one on 20 kb |

## 5. Performance

On one core of the i7-11370H: Viterbi of the eight-state CpG model over 5 Mb in 0.52 s and its forward in 0.39 s; a profile of 302 match states against 100 kb at 140 million cells a second (Viterbi) and 80 million (scaled forward), against 3 million for the log-semiring forward.

## 6. References

- Durbin, R., Eddy, S. R., Krogh, A. and Mitchison, G. (1998), *Biological Sequence Analysis*, Cambridge University Press, ch. 3 (Markov chains and HMMs, the CpG islands; the chains' tables on p. 50) and ch. 5 (profile HMMs).
- Rabiner, L. R. (1989), A tutorial on hidden Markov models and selected applications in speech recognition, *Proceedings of the IEEE* 77, 257–286.
- Baum, L. E., Petrie, T., Soules, G. and Weiss, N. (1970), A maximization technique occurring in the statistical analysis of probabilistic functions of Markov chains, *Annals of Mathematical Statistics* 41, 164–171.
- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977), Maximum likelihood from incomplete data via the EM algorithm, *Journal of the Royal Statistical Society B* 39, 1–38.
- Krogh, A., Brown, M., Mian, I. S., Sjölander, K. and Haussler, D. (1994), Hidden Markov models in computational biology: applications to protein modeling, *Journal of Molecular Biology* 235, 1501–1531.
- Eddy, S. R. (1998), Profile hidden Markov models, *Bioinformatics* 14, 755–763; Eddy, S. R. (2008), A probabilistic model of local sequence alignment that simplifies statistical significance estimation, *PLoS Computational Biology* 4, e1000069.
- Pachter, L. and Sturmfels, B. (eds.) (2005), *Algebraic Statistics for Computational Biology*, Cambridge University Press, ch. 1–2.

## 7. Open Questions

1. **The move to the inference foundation** (CLAUDE.md §8, open decision 3): `model.hpp` goes whole; `cpg.hpp` and `profile.hpp` stay here.
2. **Profile priors.** Dirichlet mixtures (Sjölander et al. 1996) and sequence weights (Henikoff and Henikoff 1994) would improve profiles from few or redundant sequences.
3. **Statistics of profile scores.** E-values for profile searches (Eddy 2008's Gumbel and exponential tails) are not yet computed; nor are the length distributions of HMMER's flanking states, so free-mode forward scores exceed HMMER's by up to log₂ L bits.
4. **Vector lanes.** The profile recurrences are scalar; Farrar's striped layout over `calx::cpu::Lanes` (as the align module uses for Smith–Waterman) would multiply Viterbi's speed.
