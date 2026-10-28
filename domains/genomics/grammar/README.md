# grammar: Group-Equivariant Grammars of Nucleic Acids

*Status: draft (2026-10-25).*

The user's draft *Towards a Group Context-Sensitive Grammar for Genome Analysis* (W. Huang, 2026), made exact for nucleic acids. The module covers:
- the symmetry of the bases and of the strands, as groups acting on letters and on words;
- the alphabets that breaking a group derives;
- grammars of structure, and their symmetry;
- the genetic code as a function on a group;
- the channels into which a group splits a genome's statistics.

It serves learning, use and research. The theory, with the proofs and the amendments to the draft, is [theory/group-grammars.md](../theory/group-grammars.md). The general machinery is morphe's: `algebra/action.hpp`, `characters.hpp` and `abelian_fourier.hpp`; and `formal/grammar.hpp`, `earley.hpp`, `inside_outside.hpp`, `equivariant.hpp` and `breaking.hpp` ([morphe decision 0006](https://github.com/huangwanhong-maker/morphe/blob/main/development/decisions/0006-representations-and-grammars-with-symmetry.md)). This module writes none of it again.

## 1. Concepts

- **The bases are a torsor of the Klein group.** Watson–Crick complementarity κ (A ↔ T, C ↔ G) and transition τ (A ↔ G, C ↔ T) generate V₄ = C₂ × C₂, which acts on the four bases freely and transitively.
  - **Subgroups.** The cosets of its three subgroups of order 2 are IUPAC's two-letter codes: R/Y, W/S and M/K.
  - **Characters.** Its three nontrivial characters are the classical binary classifications: purine–pyrimidine, weak–strong and amino–keto.
  - **Breaking.** Choosing a base point identifies the bases with the group. It fixes which class each character calls +1, but not the partition. In the domain's two-bit codes the action is exclusive or.
- **Breakings form chains.** The orbit maps of a descending chain of subgroups give the alphabets N → R/Y (or W/S, or M/K) → A, C, G, T. Grammars over a coarse alphabet are images of grammars over the bases, and parsing first at the coarse level masks the fine parse exactly.
- **The strands are a group on words.** Complement (by letters) and reversal (by anti-automorphism, a sign) generate the strand group {1, κ, ρ, κρ}; with the Klein group they generate a group of order 8.
  - **Canonical k-mers** are counted under these groups by Burnside's lemma.
  - **Measured symmetry.** Which elements a genome's k-mer statistics keep is measured: for *Escherichia coli*, only the reverse complement (Chargaff's second rule).
  - **Channels.** The isotypic decomposition of the frequency vector under the group's characters splits the statistics into orthogonal channels, one per character.
- **A grammar's symmetry is readable from its rules.** The hairpin grammar emits stem and loop from the outside in, so that it is closed under reversal. With Watson–Crick pairs it is invariant under the whole strand group; the G·T wobble pair breaks the complement and keeps reversal. morphe's equivariant parsers run on these grammars.
- **The genetic code is a function on (C₂)⁶.** Codons are a torsor of V₄³.
  - **The stabiliser** is the group of elements that keep every amino acid. It is trivial for the standard code, and is {1, τ₃}, the third-position transition, for the vertebrate mitochondrial code.
  - **The spectrum.** The code's power spectrum is the exact Walsh–Hadamard transform of its invariance profile.

## 2. Content

| Header | Content |
|---|---|
| [`nucleotides.hpp`](include/genomics/grammar/nucleotides.hpp) | `klein()`: the group, its action on the bases and its three involutions; `classifications`: the three subgroups of order 2 with their cosets as IUPAC letters; `character_value` at a base point; `quotient`: the orbit map of a subgroup, one step of a breaking; codons as six bits, `codon_bits`, `codon_bases`, `at_position` |
| [`strands.hpp`](include/genomics/grammar/strands.hpp) | `base_and_strand_group()` (order 8), `strand_group()`, `reverse_complement_group()` as morphe actions with signs; `kmer_frequencies`, `act_on_kmer`; `kmer_invariance`, the total variation between a genome's k-mer statistics and their image under each element |
| [`channels.hpp`](include/genomics/grammar/channels.hpp) | `kmer_channels`: the isotypic decomposition of the k-mer frequency vector by morphe's exact character table, each channel's energy and fraction |
| [`structures.hpp`](include/genomics/grammar/structures.hpp) | `HairpinSpec`, `allowed_pairs`; `hairpin_grammar` and `contains_hairpin_grammar` (ε-free) as morphe grammars; `is_hairpin`, `contains_hairpin` by brute force, for checks |
| [`code_symmetry.hpp`](include/genomics/grammar/code_symmetry.hpp) | `code_on_group`, `invariance_profile`, `stabiliser`, `failures`; `code_power_spectrum` by morphe's Walsh–Hadamard transform; `ignores_positions` |

## 3. Dependencies

- `genomics::sequence`: the codes, IUPAC masks, NCBI's genetic codes, and the GenBank reader in the checks.
- morphe:
  - `algebra`: group actions with signs, character tables, the Walsh–Hadamard transform;
  - `formal`: grammars, Earley's parser, inside–outside, invariant grammars, breakings.
- calx's `random` and `verification`, in the checks.

## 4. Checks

`check_genomics_grammar_grammar` holds the module to the following, with 26 checks.

| Part | What it verifies |
|---|---|
| The bases | The action free and transitive; κ, τ and σ = κτ as stated; the cosets of the subgroups are exactly the bases, the six two-letter codes and N, and B, D, H, V none; the characters multiply (χ_RY χ_WS = χ_MK); a new base point changes signs, not partitions; the breaking chain N → R/Y → bases |
| Words | Canonical k-mers under the reverse complement and under the strand group, by Burnside's lemma, against their closed forms up to k = 8 and against enumeration up to k = 6 |
| Hairpins | The grammars with and without wobble pairs generate the hairpins (every word up to length 8); the Watson–Crick grammar invariant under the strand group; with wobble, invariant under reversal and not the complement; the equivariant Earley parser (complement) and span sharing (strand group) agree with Earley's on 400 words; the flanked grammar against search |
| The code | The standard code's stabiliser trivial; τ₃ keeping 60 codons of 64 and failing at ATA/ATG and TGA/TGG; the vertebrate mitochondrial code's stabiliser {1, τ₃}; the power spectrum the Walsh–Hadamard transform of the profile, exactly, for both codes; the mitochondrial spectrum vanishing on the characters odd under τ₃; the standard code's spectrum mostly in characters of the first two bases (73.4%) |
| *Escherichia coli* | Its 86 tRNA genes carry a Watson–Crick hairpin (stem ≥ 5, loop 3–9) in 98.8% of cases, against 40.7% of random windows; its 4-mer statistics keep the reverse complement alone (TV 0.0054, the others ≥ 0.167); the eight channels' energies sum to ‖f‖², and the reverse-complement-odd channels hold 5.0 × 10⁻⁵ of it against 4.6 × 10⁻² for the complement-odd; graduated breaking on ten tRNAs, R/Y or W/S first, exact. Skipped without the data cache. |

## 5. Open Questions

1. **Groups acting on segments.** The draft's ℤ_n (tandem shifts) and S_k (exon permutations) act on positions and segments, not letters. A grammar of segmented words with S_k acting on its segments is the natural first non-abelian case; morphe has S_n's representations.
2. **Larger groups of k-mer statistics.** The hyperoctahedral group acting on k-mers contains the order-8 group. Which of its subgroups a genome's statistics keep is the draft's Q2 at the next scale.
3. **The codes' spectra across NCBI's 25 tables.** Their stabilisers and spectra compared, as an instance of the draft's Q4.

## 6. References

- Cornish-Bowden, A. (1985). Nomenclature for incompletely specified bases in nucleic acid sequences: recommendations 1984. *Nucleic Acids Research* 13, 3021–3030.
- Huang, W. (2026). Towards a group context-sensitive grammar for genome analysis. Draft.
- Rudner, R., Karkas, J. D. and Chargaff, E. (1968). *PNAS* 60, 921–922.
- Serre, J.-P. (1977). *Linear Representations of Finite Groups*. Springer.
- Terras, A. (1999). *Fourier Analysis on Finite Groups and Applications*. Cambridge.
