# Karyotypes: Their Language, Their Algebra and Their Copy Number

*Status: draft (2026-11-01).* The theory of the module [`karyotype`](../karyotype/): the description of chromosomes by ISCN as a formal language, the chromosomes of a cell as words over oriented segments of the normal ones, structural rearrangements as operations on such words, copy number as the multiset of their segments, and the descriptions of microarray and sequencing results as the images of copy-number profiles and of graphs of junctions. It serves learning, and says for every rule where it comes from.

## 1. Sources

ISCN 2024 (Hastings, Moore and Chia, eds.) and ISCN 2020 (McGowan-Jordan, Hastings and Moore, eds.) are not open. The rules used here come from open sources that quote or summarise them, each named where it is used:

- **[H25]** Hastings et al. (2025), *ISCN 2024: summary of revisions and new nomenclature* (the undisclosed sex U; der for de novo or malsegregated rearrangements, rec for recombinants; spacing; gain before loss across clones; the karyotype and microarray formats);
- **[MM21]** Miller and Madan (2021), *ISCN 2020 compared to ISCN 2016*, E.C.A. Newsletter 47 (the order of breakpoints from pter to qter for rearrangements within one chromosome; the order of inserted segments; duplications by the order of bands; dmat, dpat, dinh; the tilde; ploidies written <3n>; the sex chromosomes first in arr);
- **[W19]** the Karyotype Ontology of Warrender, Moorman and Lord (2019), whose code quotes about 180 examples of ISCN 2009 with their pages, their explanations and, for 27 of them, their detailed systems;
- **[O14]** Ordulu et al. (2014), the nomenclature of sequencing results that ISCN adopted from 2016;
- **[C25]** Conlin et al. (2025), with ISCN 2020's microarray and detailed forms;
- **[A19]** CytoGPS (Abrams et al. 2019) and its loss-gain-fusion model (Abrams et al. 2015), whose computed results ship with RCytoGPS 1.2.14;
- **[G22]** GenQA's assessment of ISCN 2020 (Hastings, Tabiner and Deans 2022);
- **[B04]** Bakhshi et al. (2004), quoting ISCN 1995's ranges of ploidy;
- the symbol tables of the Coriell Institute and the Wisconsin State Laboratory of Hygiene (secondary sources, used only for the meaning of symbols: // separates the donor and recipient lines of a chimera).

*Caveat.* Where only secondary sources were found, the rule is marked [secondary]; where the sources were silent and a choice had to be made, it is marked [convention] and the choice is the one the examples of [W19] and [A19] satisfy.

## 2. The Language

### 2.1 The Alphabet and the Grammar

**Definition 2.1.** An ISCN description is a word over bytes. Its grammar ([`grammar.hpp`](../karyotype/include/genomics/karyotype/grammar.hpp)) is context-free; morphe's BNF reader compiles it to a grammar over the atoms of the byte classes it distinguishes and parses it by Earley's algorithm, returning the derivation as a named tree, or the longest viable prefix of the input and the terminals that may follow it.

**Proposition 2.2.** *The language is not regular, and its natural grammar is unambiguous on the structure the interpretation reads.* The first because derivatives nest parentheses and breakpoint groups are matched to chromosome lists only by counting (which the grammar leaves to the rules of sense, §2.3); the second because each band, breakpoint group and abnormality begins with a symbol that determines its kind (an arm letter, an abbreviation, a sign), so that the parse trees of a text differ at most in hidden helper nodes. *Proof sketch.* The language contains the words der(1)t(1;3)(p32;q21)…, whose parentheses nest; intersecting with a regular language gives a non-regular one. For unambiguity: every concatenation in the grammar is of items whose first letters are disjoint (p, q, pter, qter and cen begin breakpoints; abbreviations begin structural abnormalities; + and − begin numerical ones; a digit begins a count), and Earley's parser therefore has at most one way to complete each visible item. ∎

### 2.2 What the Grammar Covers

Counts (46, 45~48, 44-45, ploidies <3n>); sex chromosomes (XX, XY, X, XXY…; c constitutional, ? uncertain; U undisclosed [H25]); numerical abnormalities; the structural ones add, del, der, dic, dup, fis, fra, hsr, i, idic, ider, ins, inv, neo, qdp, r, rec, rob, t, tas, trc, trp with the prefixes psu, inv, dir, dic, trc, der; breakpoints by arm, region, band and sub-band, pter, qter, cen, p10 and q10, uncertain (q2?, q?) and intervals (q21~q24; ISCN 2020's q11.2-13.1 [MM21]); multipliers x2; c and the inheritances dn, mat, pat, inh, dmat, dpat, dinh [MM21]; markers, rings of unknown origin, double minutes, inc; alternatives with or; cell counts [20] and composite karyotypes [cp10]; cell lines separated by / and, for chimeras, // [secondary]; idem, sl, sdl and sdl*n*; the detailed system with → or ->, :: and :, ? for material of unknown origin, rings opened and closed by ::, nucleotide positions and nontemplated bases [O14]; arr with build, bands, coordinates with or without commas [MM21], underscore or dash, copies x1, x2~3, proportions [0.5], hmz and inheritance, whole chromosomes (X,1-22)x2 [C25], the extended form with the normal flanks [C25]; seq with structural abnormalities in either system, regions, rearrangement sequences 1q23.3(+)(164,789,896)::8q12.2(+)(61,760,650) [O14], complexes (2,3,8)cx; ish, nuc ish and rsa recognised as text.

### 2.3 Rules of Sense

A text in the language may describe nothing. The rules that make it describe something are checked after parsing, each fault with its span: a translocation has one group of breakpoints per chromosome, each of one breakpoint or, for exchanged segments, two (ISCN 2009 p. 79 [W19]); a deletion one or two; inversions, duplications, triplications and quadruplications two; an insertion between two chromosomes one on the recipient and two on the donor, within one chromosome three; isochromosomes, additions, fissions one (fissions at p10 or q10); only derivatives (der, ider, rec) are followed by the abnormalities that made them; chromosomes are 1–22, X and Y; idem, sl and sdl refer to clones written before.

### 2.4 Editions

**Rule 2.4** [MM21, on ISCN 2020 §4.3.1.1, §9.2.5 and §9.2.9]. In ISCN 2020 and later the breakpoints of a rearrangement within one chromosome are written from pter to qter; the bands of an inserted segment with the one nearer the recipient's pter first; a duplication is direct when its bands are written in the order from pter to qter and inverted otherwise. ISCN 2016 and earlier wrote them from the centromere outward. On the long arm the two orders agree; on the short arm they are opposite, so dup(1)(p31p34) is direct in 2016 and inverted in 2020, and ins(5;2)(p14;q22q32) of 2016 is ins(5;2)(p14;q32q22) of 2020. The reader takes the edition as an argument and holds every value in 2020's order.

## 3. Bands and Points

**Definition 3.1.** For each chromosome c and arm a ∈ {p, q}, the bands of the ideogram at its finest level form a sequence B_{c,a} along the arm, each named by a word over the digits (region, band, sub-bands). A *band name* w names the set [w] of the finest bands whose words have w as a prefix.

**Proposition 3.2.** *For every name w, [w] is a run of consecutive bands; the names ordered by the prefix order form a forest (one tree per arm) whose leaves are the finest bands, and the coarser levels of NCBI's ideogram are antichains of it whose runs have the same coordinates.* The first is the numbering convention of ISCN (regions and bands numbered outward from the centromere, sub-bands subdividing their band), and is checked on the data rather than assumed: `check_bands` verifies for every prefix of every band of GRCh38 that its class is a run, and that every band of the 400 and 550 levels is the run of 850-level bands its name begins, with its own coordinates. ∎

**Definition 3.3.** A *point* of chromosome c is a position x of the boundary between two bases together with an interval I = [lo, hi) of positions it may hold: a breakpoint written as the name w lies in [w]'s extent; pter, qter and the centromere (cen, and the p10 and q10 of whole-arm rearrangements, ISCN 2009 p. 80 [W19]) are exact; a nucleotide position of a sequencing result is exact. The representative x is the middle of I; the breakpoints of one abnormality that share an interval are spread through it in the order written. Two breakpoints written at the same band of the same chromosome in different groups of one abnormality are one point (the p22 of der(1)del(1)(p34p22)ins(1;17)(p22;q11q25), ISCN 2009 p. 64 [W19]: "the proximal one is listed as the point of insertion"); breakpoints on homologues (7 and 7*) are distinct.

## 4. The Algebra of Rearrangements

### 4.1 Words

**Definition 4.1.** Let Σ be the set of *segments* (c, x, y) of reference chromosomes between points x ≠ y, written x→y, with the involution (c, x, y)⁻¹ = (c, y, x); add the letters ? (material of unknown origin), hsr and the nontemplated words over {A, C, G, T}. A chromosome of the cell is a word over Σ (linear) or a word up to rotation (a ring). The reversal of a word w = s₁⋯s_n is w⁻¹ = s_n⁻¹⋯s₁⁻¹; w and w⁻¹ are the same chromosome read from its two ends.

**Definition 4.2.** For a breakpoint x of chromosome c with centromere κ_c, the *centric piece* prox(c, x) is the word ending at x that holds κ_c (pter→x if x lies on q, qter→x if on p), and the *acentric piece* dist(c, x) the word starting at x that does not (x→qter, or x→pter). At the centromere, p10 and q10 name the arm that is the piece: prox(c, p10) = pter→p10, prox(c, q10) = qter→q10, and dist alike from the centromere outward.

### 4.2 The Operations

The structural abnormalities act as follows (each a statement of ISCN's examples quoted in [W19], pages of ISCN 2009):

| Abnormality | Word(s) | Replaces |
|---|---|---|
| t(c₁;…;c_n)(x₁;…;x_n) | der(c_k) = prox(c_k, x_k) · dist(c_{k−1}, x_{k−1}), indices mod n (pp. 77–79: the segment distal to x₁ goes to c₂, …, that of c_n to c₁) | one of each c_k |
| t(a;b)(p10;q10), whole arms | prox(a, x)·dist(b, y) and dist-of-the-other-arms (p. 80) | a, b |
| t(…)(x₁x₂;…) of segments | c_k with [x₁, x₂] replaced by c_{k−1}'s segment in the order written (p. 79) | each c_k |
| del(c)(x), del(c)(x y) | pter→x (or x→qter); pter→x · y→qter (p. 61) | c |
| dup, trp, qdp(c)(x y) | pter→x · (x→y)ⁿ⁺¹ · y→qter, the copies alternately reversed when inverted, merged (p. 69, 75, 83) | c |
| inv(c)(x y) | pter→x · y→x · y→qter (p. 72) | c |
| ins(a;b)(x;y₁y₂) | pter→x · y₁→y₂ · x→qter on a; b without [y₁, y₂] (p. 71) | a, b |
| ins(c)(x y₁ y₂) | the segment moved within c (p. 71) | c |
| add(c)(x) | pter→x · ? (or ? · x→qter) (p. 60) | c |
| hsr(c)(x) | pter→x · hsr · x→qter (p. 70) | c |
| i(c)(q10), idic(c)(x) | w·w⁻¹ for w = prox(c, x) (pp. 68, 72) | c |
| dic(a;b)(x;y), rob, der(a;b) | prox(a, x)·prox(b, y)⁻¹ (pp. 67, 81) | a, b |
| trc(a;b;c)(x;y₁y₂;z) | prox(a, x) · y₁→y₂ · prox(c, z)⁻¹ (p. 82) | a, b, c |
| r(c)(x y), r(a;b)(…) | the ring (x→y), or the ring of the segments written (p. 75) | the chromosomes |
| fis(c)(p10), (q10) | pter→p10:, qter→q10: (p. 69) | c |
| der(c)…, der(a;b)… | the abnormalities applied in order to c's word, those of a second named chromosome bringing its centric piece (pp. 60–65) | c, or a and b |
| ider(c)(q10)… | the q side of the derivative doubled (p. 64) | c |
| rec(c)dup(cq)inv(c)(x y) | qter→y · x→qter, the recombinant of a pericentric inversion (rec of [H25]; [C25]) | c |
| tas(a;b) | the chromosomes kept, associated end to end without a proven fusion ([MM21] on §9.2.16) | none in number |

**Remark 4.3.** Inside a derivative an abnormality acts on the current word: t(c;d)(x;y) cuts the word at x and keeps the side holding the derivative's centromere, replacing the other by d's piece; del removes the factor between its breakpoints; inv reverses it; dup repeats it; ins inserts at its point (reversed when the segment holding it is) or removes its segment where the derivative is its donor. der(1)t(1;3)(p32;q21)t(3;7)(q28;q11.2) is thus 7qter→7q11.2::3q28→3q21::1p32→1qter: the second translocation cuts the 3 the first attached (ISCN 2009 p. 63).

### 4.3 Orientation

**Rule 4.4** [O14, "Strand direction"; [convention] for whole arms]. A derivative is read so that the segment holding its centromere strictly runs from pter toward qter; where no segment holds a centromere strictly (the junction lies at it), from a pter if the word has one, else from the chromosome named first. This reproduces every detailed form of [W19] (der(2)t(2;5)(p12;q31) is 5qter→5q31::2p12→2qter; t(1;3)(p10;p10)'s second derivative 1qter→1q10::3q10→3qter) and of [A19] (der(2;13)(2qter→2q10::13q10→13qter)).

### 4.4 Balance

**Definition 4.5.** The *content* of a word is the multiset of the elementary segments it covers: refine all segments at all points of the cell; each segment x→y covers the elementary segments between x and y once, in either direction. The content of a cell is the sum over its chromosomes.

**Proposition 4.6.** *A translocation, an inversion, an insertion between two chromosomes and a whole-arm translocation keep the content of the chromosomes they replace; a deletion removes [x, y]; a duplication adds it once, a triplication twice, a quadruplication thrice.* *Proof.* For t(c₁;…;c_n): the derivatives are prox(c_k, x_k)·dist(c_{k−1}, x_{k−1}); summing over k, each chromosome c_j contributes prox(c_j, x_j) (in der(c_j)) and dist(c_j, x_j) (in der(c_{j+1})), whose contents are complementary in c_j because the two pieces share the point x_j (Definition 3.3: one point, not two). For whole arms the four arms are distributed between the two derivatives. For inv, the factor reversed has the same content. For ins, the segment leaves b and enters a with the same points. For dup, the word pter→x·(x→y)²·y→qter covers [x, y] twice and the rest once; trp and qdp likewise. ∎

*Remark.* The proposition is why the points of a rearrangement must be shared: had the two derivatives of t(9;22) been built with independently placed breaks in 9q34, the interpretation would report a spurious loss and gain within the band. `check_interpret` verifies balance for every example of ISCN 2009 made only of translocations, inversions and insertions.

## 5. The Cell and Its Copy Number

### 5.1 The Normal Complement

**Rule 5.1.** The ploidy k is the one written (<3n>, [MM21]) or the level whose range 23k ± 11 holds the count (near-triploid 58–80, near-tetraploid 81–103, ISCN 1995 as quoted by [B04]; [secondary]). The cell starts with k copies of each autosome and the sex chromosomes written.

**Rule 5.2** (ISCN 2009 pp. 56–57, 60, 72, 78 [W19]). The sex chromosomes written are the normal ones present: an abnormal sex chromosome is written among the abnormalities (46,X,i(X)(q10); 46,Y,t(X;8)(p22.3;q24.1); 46,t(X;Y)(q22;q11.2)) and replaces none of them; the −X and −Y of an acquired loss are already reflected (45,X,−Y), unless the sex chromosomes are marked constitutional (46,XXYc,−X: Klinefelter cells that lost an X hold XY) or the line inherits them from the clone it refers to.

**Rule 5.3** (ISCN 2009 pp. 60, 80–82 [W19]). A structural abnormality replaces the normal chromosomes it is made from (the table of §4.2); with + it is additional; with − before a numerical change a normal chromosome is missing, and before a structural one in a subclone the abnormal chromosome of the clone referred to is missing while the chromosomes it replaced stay missing (−dic(1;20) in a sideline). The −10 of 45,XY,−10,der(10)t(10;17)(q22;p12) is the loss of the homologue the derivative did not replace ("the der(10) replaces a normal chromosome 10; the homologous chromosome 10 is lost").

**Rule 5.4.** idem refers to the stemline (the first clone), sl to the stemline, sdl*n* to the *n*-th clone after it ([A19]: "idem in this clone refers to the first clone"); double minutes are not counted among the chromosomes [convention]. The number of chromosomes made is compared with the count written; a disagreement is reported (a note for composite and incomplete karyotypes, whose counts are ranges over cells).

**Rule 5.5.** The reference is the normal diploid of the constitutional sex (male when a Y is written anywhere) times the ploidy: autosomes k, X and Y k and 0 for a female, k − ⌊k/2⌋ and ⌊k/2⌋ for a male. Under U the sex chromosomes have no reference.

### 5.2 Copy Number

**Definition 5.6.** The *coverage* of chromosome c in a cell is the function on positions counting the segments of c, over all the cell's chromosomes, that cover each position; it is a step function on the representative points. Its maximal runs at a value other than the reference are the *changes*: each with the bands at its ends as written (pter and qter as the terminal bands, the centromere as the band beside it), its representative extent, its *inner* extent from the far end of its first point's interval to the near end of its last's (certainly affected), and its *outer* extent from the near end of the first to the far end of the last (possibly affected).

**Proposition 5.7.** *If no breakpoint's interval meets a band B, the copy number on B is the same for every placement of the breakpoints within their intervals consistent with the order written.* *Proof.* The coverage changes only at points; a placement moves each point within its interval, which does not meet B, so B lies on the same side of every point in every placement, and its coverage is the same. ∎

The copy number of a band is reported as the set of values of the runs whose outer extents meet it: a single value under the hypothesis of Proposition 5.7, an interval otherwise.

*Caveat.* When breakpoints on different homologues share a band, their relative order is not given by the description; they are placed in the order written, and the band's interval of values is that of this placement. The inner and outer extents of each change are exact whatever the placement.

### 5.3 Loss, Gain and Fusion

**Definition 5.8** (the LGF model of [A19]). For each abnormality, its *loss* is the set of bands where its chromosomes cover fewer copies than the normal chromosomes it replaces, its *gain* where more, over the outer extents; its *fusion* the bands of the points of its junctions (the consecutive segments of its words and the closing of rings), a centromeric junction named by its arm (2q10). The model of a cell line is the union over its abnormalities.

*Remark.* The model is of events, not of net copy number: in 44~45,XX,+2,der(2;13)(q10;q10),… the band 2p25.3 is both lost (by the derivative) and gained (by +2), as CytoGPS computes it; the net copy number of §5.2 says 2p has two copies. `check_cytogps` holds the model to CytoGPS's results band for band.

## 6. Describing Results

### 6.1 Microarrays

A copy-number segment [s, e) of chromosome c at copies k is written c, the names of the bands holding its first and last bases (one name when they are the same), (s + 1_e), xk (ISCN 2020 ch. 14 as [C25] shows it: arr[GRCh38] 1q21.1q21.2(147,105,904_147,917,509)x1); a mosaic range x2~3, a proportion [0.3], hmz for a run of homozygosity; a whole chromosome (13)x3, chromosomes at the same copies grouped as (1-20,22)x3 [G22]; the sex chromosomes first [MM21]; a normal result arr(X,1-22)x2 or arr(X,Y)x1,(1-22)x2 [C25].

### 6.2 Sequencing: The Graph of Junctions

**Definition 6.1.** Given junctions (adjacencies) of a genome, each joining one side of a break on one chromosome to one side of a break on another, cut every chromosome at every break into *pieces*. The *junction graph* has the two ends of every piece as vertices, an edge between the two ends of each piece, an edge for each junction, and, at a break neither side of which is joined, the reference edge between the two pieces' facing ends.

**Proposition 6.2.** *Every vertex has degree at most two, so the graph is a disjoint union of paths and cycles; the paths whose ends are telomeres or unjoined ends, read along their piece edges, are the linear chromosomes the junctions make, and the cycles are rings.* *Proof.* Each end lies on exactly one piece edge, and on at most one junction or reference edge (an end joined twice is reported and the later junction kept). A graph of maximum degree two is a union of paths and cycles. ∎

A path holding a centromere strictly is a derivative, named by it and read by Rule 4.4; acentric paths are reported as lost or unplaced pieces. Two derivatives exchanging the distal pieces of two chromosomes are a translocation t(a;b)(x;y); one derivative of one chromosome skipping a piece forward is a deletion, returning backward a tandem duplication, with one piece reversed an inversion; one derivative with another chromosome's piece and no partner an unbalanced der(a)t(a;b)(x;y); anything else a complex rearrangement (a,b,…)cx written derivative by derivative [O14]. Each junction is also written as a rearrangement sequence, with its strand and its base: a left piece is read forward to its last base, a right piece from its first [O14].

*Remark.* The circle closes: a sequencing description written from junctions, read back by the grammar and interpreted, gives the same copy number at the same coordinates (`check_describe`).

## 7. Where the Code Implements It

| Section | Code |
|---|---|
| §2 | `grammar.hpp` (the BNF), `iscn.hpp` (`parse`, the rules of sense, `to_edition_2020`), `format.hpp` (the canonical form) |
| §3 | `bands.hpp` (`Ideogram::find`, `band_at`, `range_name`), `interpret.hpp` (`Interpreter::resolve`, `resolve_all`) |
| §4 | `interpret.hpp` (`prox`, `dist`, `realise`, `apply`, `orient`, `merged`) |
| §5 | `interpret.hpp` (`interpret`, `coverage`, the changes, band copies, events and LGF), `long_form`, `detailed` |
| §6 | `describe.hpp` (`describe_arr`, `describe_seq`) |
