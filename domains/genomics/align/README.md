# align

*Status: draft. Module of the genomics domain.* Pairwise alignment: the biological scoring of sequences (substitution matrices, affine gaps, Karlin–Altschul statistics) on morphe's generic alignment over semirings, which defines the answers, and the algorithms that compute them fast, each a schedule of the same sum checked against the generic one. Serves use (alignment for mapping, assembly, variants and multiple alignment) and learning (the theory note [alignment.md](../theory/alignment.md)).

## 1. Concepts

- **An alignment is a path; a scoring is a weight in a semiring; an algorithm is a schedule.** The alignments of a query and a target are the paths of their edit graph, a grid of (m + 1)(n + 1) vertices with diagonal, vertical and horizontal edges. A scoring weighs every path in a semiring, and an alignment problem is the sum over the paths a mode allows: the optimal score in (max, +), the edit distance in (min, +), the probability of the pair HMM in (+, ×), the number of optimal alignments in the tropical semiring with multiplicities, and every optimal alignment for every parameter at once in the polytope algebra (theory §2, §6).
- **The definition is morphe's.** One generic, slow, obviously correct dynamic programme over the edit graph, in any semiring (`morphe/formal/alignment.hpp`), names no object in nature and lives in the foundation. `reference.hpp` gives it the domain's scorings; every algorithm of this module is checked against it.
- **The algorithms are schedules.** Gotoh's recurrences row by row, Myers and Miller's divide and conquer in linear space, the band, the X-drop's pruned frontier, and the kernels on vector lanes (bit vectors, stripes, wavefronts) compute the same tropical sum in different orders and representations (theory §4).
- **The four modes** say where paths begin and end: global, local, semiglobal (the query end to end, the target's ends free) and overlap (free end gaps on both). A gap of length L scores −(open + L · extend); scores are maximised, penalties minimised (`types.hpp`).

## 2. Content by Header

| Header | Content | Status |
|---|---|---|
| `types.hpp` | Modes, matrices, scorings and penalties, CIGARs, alignments; `rescore`, `penalty`, `edit_distance`, `identity`, `render`, `explicit_matches` | draft |
| `scoring.hpp` | The NCBI protein matrices BLOSUM45, 50, 62, 80, 90 and PAM30, 70, 250, transcribed exactly in the code order (`matrix("BLOSUM62")`, `matrix_names()`), and BLAST's built-in half-bit BLOSUM80 (`"BLOSUM80-BLAST"`, §7); DNA presets of BLAST, BWA-MEM and minimap2; `blastp()`; Robinson and Robinson's amino-acid frequencies | draft |
| `reference.hpp` | The domain's scoring as morphe's three-state affine scheme in (max, +): `reference_score`, `reference_alignment` (generic traceback), `reference_count` (optimum and the exact number of optimal alignments), `reference_alignments` (all of them) | draft |
| `gotoh.hpp` | `gotoh(query, target, scoring, mode)`: Gotoh's recurrences with traceback (one byte of directions per vertex), deterministic ties (left-aligned gaps); `gotoh_score`, `gotoh_end` and `local_end` in O(n) memory; `gotoh_field`, the whole H, I, D field for drawing | draft |
| `hirschberg.hpp` | `hirschberg(query, target, scoring, mode)`: Myers and Miller's linear-space alignment with affine gaps, global, and the other modes by locating the ends first | draft |
| `banded.hpp` | `banded(query, target, scoring, mode, low, high)`: the band of diagonals j − i ∈ [low, high]; `banded_is_optimal`, the certificate for global alignment; `banded_global`, doubling bands until certified | draft |
| `xdrop.hpp` | `ungapped_xdrop`; `extend_right`, `extend_left` and `xdrop`: gapped X-drop extension from an anchor, with minimap2's Z-drop (`Drop {x, z}`) | draft |
| `statistics.hpp` | `score_distribution`, `karlin_lambda`, `karlin_H`, `karlin_K` (closed forms and the full series), `karlin_K_blast` (BLAST's run-time truncation), `karlin_altschul`; `bit_score`, `evalue`, `pvalue`; NCBI BLAST's gapped table (`gapped_parameters`, `gapped_table`) | draft |
| `parametric.hpp` | `parametric(query, target, dna)`: the alignment polygon of (mismatches, gaps) by polytope propagation, each vertex with its counts, a witness CIGAR and its region in the (α, β) plane; `optimal`, `optimal_alignment` at a parameter | draft |
| `pairhmm.hpp` | `PairHmm` (DNA, or from a matrix's target frequencies); `viterbi`, `forward`, `backward` in log space, `posterior`, `mea_alignment` | draft |
| `bitparallel.hpp` | Myers' bit-vector edit distance (unit costs, N matches nothing): blocks of 64 rows, Ukkonen's band in blocks with proved rules, a doubling threshold; every semiglobal end within k (`myers_search`); global and semiglobal distances and alignments with their CIGAR from the kept band (`myers_distance`, `myers_align`); `MyersStats` | draft |
| `striped.hpp` | Farrar's striped Smith–Waterman with affine gaps on calx lanes: the query profile, 8-bit biased lanes then 16 and 32 on overflow, lazy F with E (exactly Gotoh's scores), the end by Gotoh's rule, the start from the reversed prefixes, the path by `banded_global` (`StripedProfile`, `striped_end`, `striped_local`); `StripedStats` | draft |
| `inter.hpp` | One query against many targets, one target to a lane (SWIPE): local scores in 16-bit lanes, saturated lanes rescored by the striped kernel (`inter_scores`); `InterStats` | draft |
| `wavefront.hpp` | The wavefront algorithm with gap-affine penalties: global and ends-free with traceback (`wavefront`), the penalty alone in O(s) memory (`wavefront_penalty`), BiWFA in O(s) memory with breakpoints inside gaps (`biwfa`), WF-adaptive (off by default), the penalty form of a scoring and back (Eizenga and Paten 2022: `penalty_form`, `score_from_penalty`, `wavefront_global`); `WavefrontStats`, and `WavefrontRecord` for drawing the fronts | draft |
| `codon.hpp` | `codon_alignment(cds_x, cds_y, code, scoring)`: two coding sequences aligned through their proteins under any NCBI genetic code, the gaps carried back to whole codons (PAL2NAL's method); final stops and the partial codons of mitochondrial genes left out and counted | draft |
| `poa.hpp` | Partial-order alignment: `PoaGraph` with aligned columns ordered by Kahn's algorithm, global and local affine alignment of a sequence to the DAG, fusion, the heaviest-bundle consensus, the induced multiple alignment | draft |

## 3. Conventions

- **Query and target.** The query indexes the rows (i = 0 .. m), the target the columns (j = 0 .. n); CIGAR I is a query base against a gap, D a target base against a gap.
- **Ties.** Every algorithm returns one alignment, the same on every machine. Gotoh's rule (`gotoh.hpp`): the end is the first optimal vertex in row-major order; at a vertex, start before pair before I before D; in a gap, opening before extending. Gaps are therefore as far left as the score allows. Other schedules return an optimal path that may differ among co-optimal ones; the checks hold every returned path to its claimed score by `rescore`.
- **Local alignment of score 0** is the empty alignment at (0, 0).
- **Codes** index the matrix directly (`sequence/alphabet.hpp`): DNA A, C, G, T, N = 0 .. 4; protein in the NCBI order ARNDCQEGHILKMFPSTWYVBZX*.

## 4. Dependencies

`genomics::sequence`; `morphe::formal` (the generic alignment, the semirings, the polytope algebra) and `morphe::arithmetic` (exact counts); `calx::solvers` (Brent's method for λ); `calx::cpu` (the lanes of the fast kernels). The checks use `genomics::simulate`, `calx::random` and `calx::verification`.

## 5. Checks

`node domains/genomics/tests/run-tests.mjs`, or one program at a time from the build (`check_genomics_align_<topic>`).

| Check | What it verifies |
|---|---|
| gotoh | The NCBI matrices symmetric with their files' entries and extremes; Robinson and Robinson's frequencies a distribution; Durbin et al.'s HEAGAWGHEE and PAWHEAE (BLOSUM50, linear gap −8): global 1 (HEAGAWGHE-E / --P-AW-HEAE), local 28 (AWGHE / AW-HE), overlap 25 (GAWGHEE / PAW-HEA); Gotoh equal to morphe's generic sum on 400 random DNA and protein pairs in each mode, every path allowed by its mode and rescored to its claim, the linear-memory score and end the traceback's; the field for drawing; homologous pairs from the simulator; the k co-optimal alignments of A^k against A^(k−1) and the leftmost gap; the speed on 10 kb × 10 kb |
| hirschberg | Myers–Miller's score equal to Gotoh's on 800 random pairs per mode under five scorings (linear gaps, large opens, free extension, BLOSUM62), paths rescored; gaps through the middle row; a homologous pair of 4 kb |
| banded | The full band is Gotoh path for path; the band of an optimal path's diagonals reaches the optimum; narrower bands never more, paths inside the band; the certificate never wrong; doubling bands across indels of 41–50 bases; bands without a path refused |
| xdrop | With X = ∞ the extensions are the best anchored alignments, by brute force over prefixes and suffixes; bounded and anchored for finite X, rescored; the ungapped extension the best prefix sums and growing with X; a planted core recovered from its middle; the Z-drop stopping at a divergent block that the X-drop alone crosses |
| statistics | λ, K, H of the eight matrices under Robinson and Robinson against NCBI BLAST's ungapped values (BLOSUM62: 0.3176, 0.134, 0.4012; BLOSUM80 in BLAST's half-bit form, the FTP file's at λ ≈ ln 2 / 3); BLAST's truncated K above the full series; the series for K against its closed forms; target frequencies a distribution; E = m n 2^{−S′}; BLAST's gapped table (BLOSUM62 11/1: 0.267, 0.041, 0.14); simulation of 3000 random protein pairs against the extreme-value law |
| parametric | At random (α, β) the polygon's optimum equals the best enumerated alignment, in the selected vertex's region; witnesses with their vertex's counts; regions covering the plane; each vertex Gotoh's unique optimum at an interior direction; protein; polygon sizes within the lattice bound |
| bitparallel | Textbook distances (kitten/sitting 3, vintner/writers 5, survey/surgery); global and semiglobal distances and ends against the O(mn) recurrence for lengths 1–5000 across every block boundary; thresholds refused below the distance; every CIGAR's edit distance its claim; the band against the unbanded run on 30 kb (5% of the cells); a read located in 200 kb |
| striped | Scores and ends equal to Gotoh's local optimum on 306 DNA and protein pairs under six scorings (blastn, BWA-MEM, minimap2, blastp, BLOSUM80, linear); full alignments rescored; the 8 → 16 → 32-bit passes driven by scores of 200, 1916 and 38506; lazy F on long insertions |
| inter | Batch scores equal to Gotoh's for DNA and protein; 16-bit saturation rescored |
| wavefront | The WFA2-lib example's penalty 24 by WFA, BiWFA and Gotoh; penalties equal to Gotoh's optimum on 382 pairs under seven penalty sets, global, and 378 ends-free; CIGARs holding their penalty; score-only equal to full; BiWFA optimal everywhere, and on 30 kb at 1% in 253 KB against the WFA's 16.5 MB; 8100 sub-problems over every begin and end state against a three-state reference; unit penalties give Myers' distance; the Eizenga–Paten identity on 180 paths; WF-adaptive valid; the record's fronts rising to the end |
| codon | A gene against itself; planted codon indels found whole at their place; every codon once and in order and every gap on a codon boundary on 200 related pairs; synonymous changes leaving the proteins ungapped; the vertebrate mitochondrial code (TGA tryptophan, AGA a stop) and its partial stops |
| poa | On a chain the alignment is Gotoh's, global and local under four scorings (800 cases); on graphs of 2–4 sequences, brute force over every source-to-sink path (600); edges forward in the order; the consensus of 20 long reads (8% errors) of 2 kb within 4 edits of the truth (2 found); MSA rows the sequences with gaps; sum-of-pairs recall on a simulated family |
| pairhmm | Forward = backward; P(x, y) the generic sum in (+, ×) and the sum over enumerated alignments; Viterbi the generic sum in (max, +) of log weights; posteriors in [0, 1], rows and columns summing to at most 1, totalling the expected aligned pairs; MEA at least Viterbi's expected pairs; 1 kb sequences in log space; BLOSUM62's target frequencies |

The checks of morphe's generic alignment itself (every mode and gap model against the enumeration of all paths, Delannoy numbers, the semiring and homomorphism laws, the alignment polygon against the hull of all vectors) are morphe's (`foundations/morphe/formal/tests/`).

## 6. Performance

Measured on the development machine (i7-11370H, one thread, GCC 11.5 `-O3 -march=native`, 10 kb × 10 kb homologous DNA, BWA-MEM's scoring), while other builds shared the machine: Gotoh with traceback 141 Mcells/s (0.71 s), score only 347 Mcells/s, local with traceback 110 Mcells/s at a load of about 3; 70–100, 200–250 and 60–90 Mcells/s at a load of 5–10. A branch-free variant of the inner loop showed no reliable gain and was not kept. The kernels on lanes and the wavefronts, on one thread ([theory](../theory/fast-alignment.md), §7): Smith–Waterman of 1 kb against 100 kb, scalar Gotoh 0.23 GCUPS, striped 16.4 GCUPS in 8-bit lanes (×70) and 6.6 GCUPS when the score needs 16 bits, Myers' search with k = 150 18.5 GCUPS; 320 targets of 10 kb one to a lane 8.5 GCUPS. Global alignment of 10 kb pairs (penalties 4, 6, 2): Gotoh 0.5 s at every divergence; the WFA 0.8 ms at 1%, 50 ms at 5%, 191 ms at 10% (with its path; 0.8, 21 and 60 ms for the penalty alone), BiWFA 3.3, 96 and 293 ms in 0.11, 0.43 and 0.77 MB against the WFA's 1.7, 47.5 and 151.5 MB. On eight threads in the Genome Atlas, 16,000 reads of 150 bases against a query of 120: Gotoh 1.2 GCUPS, striped 26, one target to a lane 31, every score equal.

## 7. Data

The protein matrices are NCBI's files, public domain, retrieved with curl from `https://ftp.ncbi.nih.gov/blast/matrices/` on 2026-10-24 and transcribed into `scoring.hpp` by script; their SHA-256:

| File | SHA-256 |
|---|---|
| BLOSUM45 | `19a077d5060a94a6212ad609df303a240463873a69694e9f576c5762fbffdb09` |
| BLOSUM50 | `a152e30b8afcab42281ff8604bc51358be98253155fde62ae20707805a0d3d67` |
| BLOSUM62 | `85510d3846ee6d5f4778e425cf8daf6e0dbb889b306f2d13434e1254780efb40` |
| BLOSUM80 | `2e597176b8092808ba199cfc9d80dc8d657b5a3e5015dfb81d557285819ce187` |
| BLOSUM90 | `4de6ce225ed8ee84e041c4f63c2daf2f43aa86cdd1b120868074dc8495c7214f` |
| PAM30 | `8fc5688d1798e2326eab3719b9ee9299952a20c9434e4c5528bf4532e4586e63` |
| PAM70 | `d180e22d80f25ada412d6eb63070fe7227c9547c6673727291fd9b2d452946e5` |
| PAM250 | `936cc5ef23e0262c5fdaf74d5d81f58848eb74cca0d5bef0f396d27b1689704e` |

BLAST's own matrices (`src/util/tables/sm_<matrix>.c` of the NCBI C++ Toolkit, retrieved 2026-10-24) equal the FTP files on the twenty standard residues, and differ in the entries of B, Z and X, for every matrix but BLOSUM80: BLAST's BLOSUM80 is at a scale of ln 2 / 2 ("Entries for the BLOSUM80 matrix at a scale of ln(2)/2.0"), the FTP file at ln 2 / 3 (made from `blosum80_3.iij`), and they differ in 324 of the 400 standard entries. BLAST's statistics for "BLOSUM80" belong to its own matrix, which `scoring.hpp` therefore also holds, as `BLOSUM80-BLAST`, transcribed by script from `sm_blosum80.c` (SHA-256 `5e2f92198415f2f7cfcd3a2482a10d89d4b47172c2d2f16949fb8488ec046297`) without its column J; `statistics.hpp` files BLAST's BLOSUM80 rows under that name.

Robinson and Robinson's frequencies and the gapped parameters are NCBI BLAST's, from `blast_stat.c` of the NCBI C++ Toolkit (`https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/blast/core/blast_stat.c`, retrieved 2026-10-24, SHA-256 `d3f2839626931d136f0de8c27aad8731f15b8ad93ac1df25be12d5cb527300a3`), the gapped rows transcribed by script.

## 8. References

- Altschul, S. F. and Erickson, B. W. (1986), Optimal sequence alignment using affine gap costs, *Bull. Math. Biol.* 48, 603–616.
- Altschul, S. F. et al. (1990), Basic local alignment search tool, *J. Mol. Biol.* 215, 403–410; (1997), Gapped BLAST and PSI-BLAST, *Nucleic Acids Res.* 25, 3389–3402.
- Altschul, S. F. and Gish, W. (1996), Local alignment statistics, *Methods Enzymol.* 266, 460–480.
- Chao, K.-M., Pearson, W. R. and Miller, W. (1992), Aligning two sequences within a specified diagonal band, *CABIOS* 8, 481–487.
- Do, C. B., Mahabhashyam, M. S. P., Brudno, M. and Batzoglou, S. (2005), ProbCons, *Genome Res.* 15, 330–340.
- Durbin, R., Eddy, S. R., Krogh, A. and Mitchison, G. (1998), *Biological Sequence Analysis*, Cambridge University Press, ch. 2 and 4.
- Gotoh, O. (1982), An improved algorithm for matching biological sequences, *J. Mol. Biol.* 162, 705–708.
- Gusfield, D., Balasubramanian, K. and Naor, D. (1994), Parametric optimization of sequence alignment, *Algorithmica* 12, 312–326.
- Henikoff, S. and Henikoff, J. G. (1992), Amino acid substitution matrices from protein blocks, *PNAS* 89, 10915–10919.
- Holmes, I. and Durbin, R. (1998), Dynamic programming alignment accuracy, *J. Comput. Biol.* 5, 493–504.
- Karlin, S. and Altschul, S. F. (1990), Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes, *PNAS* 87, 2264–2268.
- Li, H. (2018), Minimap2: pairwise alignment for nucleotide sequences, *Bioinformatics* 34, 3094–3100.
- Myers, E. W. and Miller, W. (1988), Optimal alignments in linear space, *CABIOS* 4, 11–17.
- Pachter, L. and Sturmfels, B. (2004), Parametric inference for biological sequence analysis, *PNAS* 101, 16138–16143.
- Robinson, A. B. and Robinson, L. R. (1991), Distribution of glutamine and asparagine residues and their near neighbors in peptides and proteins, *PNAS* 88, 8880–8884.
- Smith, T. F. and Waterman, M. S. (1981), Identification of common molecular subsequences, *J. Mol. Biol.* 147, 195–197.
- Waterman, M. S., Eggert, M. and Lander, E. (1992), Parametric sequence comparisons, *PNAS* 89, 6090–6093.
- Zhang, Z., Schwartz, S., Wagner, L. and Miller, W. (2000), A greedy algorithm for aligning DNA sequences, *J. Comput. Biol.* 7, 203–214.

## 9. Open Questions

1. **The local score of Durbin et al.'s example.** The book's Fig. 2.6 gives 28 (AWGHE / AW-HE), which the code reproduces; the page numbers of Figs. 2.5–2.8 are cited by figure, not by page, until the book is at hand.
2. **Gapped statistics.** λ and K of gapped alignment are taken from BLAST's simulated table; estimating them for any scoring by simulation (the islands method of Altschul et al. 2001) would free the module from the table.
3. **B, Z and X.** NCBI's FTP files and BLAST's built-in matrices disagree on the ambiguity codes; the module follows the FTP files, as asked, and BLAST only for BLOSUM80's statistics. Which to prefer for ambiguous residues is open.
4. **Edge corrections.** E-values here are the asymptotic K m n e^{−λS}; BLAST's effective lengths (α, β of the table, or Spouge's finite-size correction with C, α_v, σ) are not yet applied.
5. **Parametric alignment beyond global.** The plane projection is exact only for global alignment; local and semiglobal parametric alignment need the polytope of (matches, mismatches, gaps) in three dimensions, which the polytope algebra of morphe does not yet have.
6. **Affine parametric alignment.** Three parameters (mismatch, open, extend) also need a 3-polytope.
