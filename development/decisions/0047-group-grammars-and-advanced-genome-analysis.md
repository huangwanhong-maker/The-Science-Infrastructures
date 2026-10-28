# 0047. Group Grammars for Genomes, and a Survey of Advanced Genome Analysis

- **Status:** Accepted, provisionally; amends 0043
- **Date:** 2026-10-25
- **Decided by:**
  - **The request.** The user gave Claude their draft, *Towards a Group Context-Sensitive Grammar for Genome Analysis* (W. Huang, 2026). They asked it to "survey advanced genome analysis algorithms and implement them, also include ours (paper attached, implement it first)".
  - **The guidance.** They added that the draft's formalism is preliminary, and that where it cannot be reproduced the implementation should follow the idea's direction correctly.
  - **What Claude did.** It wrote the precise form of the draft's mathematics and its proofs, the morphe headers of actions, grammars and their parsers, and the genomics module `grammar`. It chose the survey's algorithms and wrote this record and the survey. Agents wrote the rest to the contracts given them, each module checked:
    - representations and cyclotomic numbers;
    - the r-index and the PBWT;
    - signals, classification, seeds and whole-genome comparison;
    - haplotypes and selection scans;
    - RNA structure;
    - the Workbench's integration.

## 1. Context

The draft makes a finite group G the primitive of a grammar of genomes. Concrete letters arise from G's group algebra by symmetry-breaking maps, parsing takes rules up to G, and the Wedderburn blocks of ℂ[G] are proposed as independent parsing channels. Its open questions:
- Q1: grammar induction by an equivariant inside–outside;
- Q2: discovering G from data;
- Q4: one G with different breakings across species;
- Q5: non-abelian groups.

Some of its statements hold as written; others need restating to be true; a few do not hold. Following the idea's direction rather than the letter is what the user asked for.

The genomics and evolution domains (ADR 0043) already held the classical algorithms of alignment, indexing, mapping, assembly, variants, annotation, expression, trees and populations. Several families of advanced methods were missing:
- RNA structure;
- spectral analysis;
- taxonomic classification;
- whole-genome comparison;
- haplotype models and selection scans;
- compressed indexes of repetitive collections and panels;
- seeds robust to indels.

## 2. Decision

1. **The precise form of the draft.** The note `domains/genomics/theory/group-grammars.md` states the framework as definitions, proves its results, and lists eleven amendments to the draft, each with its reason. In brief:
   - Burnside's bound runs the other way;
   - matrix units are not idempotents;
   - ℂ[G] as an algebra does not determine G (as a Hopf algebra it does);
   - purine–pyrimidine is a character of the Klein group, not derived from ℤ₂ alone;
   - the per-position orbit scan decides a larger language;
   - the worst-case cost is not divided by |G|;
   - a G-invariant language lies in one Wedderburn block, so the blocks parallelise equivariant linear maps, not recognition;
   - the breakings of a torsor differ by right translations;
   - G_DNA's factors act on positions or segments, not on letters;
   - the complexity class is unchanged;
   - the two phases are a grammar and a letter-to-letter map.
2. **The general mathematics in morphe** (morphe decision 0006):
   - group actions with signs;
   - cyclotomic numbers;
   - exact character tables of any finite group (Dixon–Schneider);
   - group algebras and their Wedderburn blocks;
   - Young's seminormal form and the Fourier transform on S_n;
   - exact abelian and Walsh–Hadamard transforms;
   - grammars (Earley's parser, inside–outside over semirings, monotone search);
   - grammars with symmetry: rule orbits, the equivariant Earley parser with coset binding, shared spans, orbit-canonical search, breakings and exact coarse-to-fine masks;
   - the r-index and the positional BWT;
   - increasing subsequences.
   `formal` now depends on `algebra`.
3. **The genomics module `grammar`**: the draft made exact for nucleic acids.
   - **The bases.** The Klein group of complementarity and transition acts simply transitively on A, C, G, T. IUPAC's two-letter codes are its cosets, and the purine–pyrimidine, weak–strong and amino–keto classifications are its characters.
   - **Breakings and their chains.**
   - **Strands.** The strand group on words, with canonical k-mers counted by Burnside's lemma.
   - **The symmetries a genome's statistics keep, measured:** Chargaff's second rule as the group {1, κρ}, with its isotypic channels.
   - **Hairpin grammars,** whose symmetry the wobble pair breaks.
   - **The genetic code** as a function on (C₂)⁶, with its stabiliser and exact spectrum.
4. **The survey**, `domains/genomics/theory/the-landscape-of-genome-analysis.md`, maps the field's problems to their algorithms and sources, and says which are held, added here, or planned. The added modules:
   - **genomics `rna`**: Nussinov's pairing, stochastic grammars of secondary structure trained by inside–outside with the orbit-tied (equivariant) EM, pair probabilities and maximum expected accuracy structures.
   - **genomics `signal`**: spectral analysis on calx's FFT, the period-3 signal of coding DNA.
   - **genomics `classify`**: Kraken's lowest-common-ancestor classification, on a taxonomy read from GenBank lineages.
   - **genomics `kmers`**: strobemers. Syncmers were already in `sequence`; their density is now checked exactly.
   - **genomics `comparative`**: whole-genome comparison by maximal unique matches and their chains.
   - **evolution `haplotypes`**: Li and Stephens's copying model with phasing and imputation, and EHH and iHS.
   - **evolution `population`**: a forward simulator of sweeps.
5. **The Genome Workbench**: a Grammar view and Toolbox tools that make these usable and visible, through engine operations on the domains.

## 3. Alternatives Considered

- **Implementing the draft literally.** Rejected where it fails: its scan rule accepts words outside the language (the checks show it), and its complexity and channel claims do not hold. The user asked for the idea's direction, implemented correctly.
- **A domain of its own for grammars of nucleic acids.** Rejected: they are genomics, and the general part is morphe's.
- **Free-energy RNA folding** (Zuker, McCaskill with nearest-neighbour parameters). Deferred: its parameters are data that would have to be fetched, not typed; the grammar-based folding needs none.
- **Deep equivariant encoders** (the draft's §8). Deferred to the Python environment. The exact parts are morphe's characters and group algebras: isotypic decomposition, and Schur's block-diagonal equivariant maps.

## 4. Consequences

- CLAUDE.md §§1 and 9, the README, the domains' and morphe's READMEs and the survey record the new modules; the theory notes give their proofs.
- The domains write no general method of their own: morphe holds the groups, grammars and indexes, and calx the FFT.
- The draft's open questions now stand as follows; the note's §9 has the detail.
  - Q1 (equivariant EM) is answered as far as convergence goes, and run on real tRNAs.
  - Q2 is answered by measurement for genomic statistics.
  - Q4 has a small exact instance: the genetic codes as functions on one group.
  - Q5 has its tools, the non-abelian representations.
- **Open.**
  - Free-energy folding.
  - Long-read assembly.
  - Graph genomes.
  - Structural variants.
  - The draft's Q3 (treewidth) and Q6.
