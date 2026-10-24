# Fast Alignment

*Status: draft. Serves learning and use.* The fast kernels of alignment as schedules of one sum: bit vectors that hold a column of the unit-cost table as its differences and carry a match down the column by an addition; stripes and lanes that compute the local recurrences many cells at a time; wavefronts that compute the tropical sum in the order of its values; and the edit graph of a partial order, which aligns a sequence to many at once. The definitions are those of [Alignment as a Sum over Paths](alignment.md) (the edit graph, §1; semirings, §2; affine gaps, §3); the code is the module [`align`](../align/): `bitparallel.hpp`, `striped.hpp`, `inter.hpp`, `wavefront.hpp`, `poa.hpp`; the design is the idea note [The Genome Infrastructure](../../../development/ideas/2026-10-24-the-genome-infrastructure.md), §1.1.

## 1. Schedules of One Sum

The value of an alignment is a sum over the paths of the edit graph G(x, y) in a semiring (alignment.md, Definition 2.2); with affine gaps, over the paths of the product of G with the three-state automaton (alignment.md, Proposition 3.2). A *schedule* of the sum chooses which vertices to visit, in what order, in what representation, and what to keep. Gotoh's recurrences row by row are one schedule; each kernel here is another, for one semiring, and each is checked against Gotoh's (and so, through `check_gotoh`, against morphe's generic definition).

| Kernel | Semiring | Order of the vertices | Representation | Work |
|---|---|---|---|---|
| Myers (`bitparallel.hpp`) | (min, +), unit costs | column by column, 64 rows a word | the differences of a column, 2 bits a vertex | O(⌈k/64⌉ n) for the threshold k |
| Farrar (`striped.hpp`) | (max, +), affine, local | column by column, W stripes of the query | lanes of 8, 16 or 32 bits | O(⌈m/W⌉ n) vector operations |
| SWIPE (`inter.hpp`) | (max, +), affine, local | W problems at once, the same cell in each lane | lanes of 16 bits | O(m n) vector operations for W targets |
| WFA (`wavefront.hpp`) | (min, +), affine penalties | by increasing penalty, along the diagonals | the furthest offset of each diagonal | O((n + m) s) at worst, O(n + s²) typical |
| POA (`poa.hpp`) | (max, +), affine | a topological order of the DAG × the sequence | rows of the nodes | O((V + E) n) |

*Remark.* The semiring decides what a schedule may skip. The tropical semirings are selective, a ⊕ b ∈ {a, b}, so a vertex whose value cannot reach the optimum may be dropped without changing it: the band of the bit vectors and the wavefronts both rest on this, and neither exists in the counting or probability semirings, where every path contributes.

## 2. Bit-Parallel Columns

Let C be the unit-cost table of x (the query, length m, rows i) and y (the target, length n, columns j): C(i, 0) = i, C(0, j) = j in global mode and 0 in semiglobal mode, and C(i, j) = min(C(i − 1, j − 1) + [x_i ≠ y_j], C(i − 1, j) + 1, C(i, j − 1) + 1), where x_i = y_j counts only for codes below σ (N matches nothing).

**Proposition 2.1 (Ukkonen 1985).** For all i, j ≥ 1:

- C(i, j) − C(i − 1, j) ∈ {−1, 0, 1} and C(i, j) − C(i, j − 1) ∈ {−1, 0, 1};
- C(i, j) − C(i − 1, j − 1) ∈ {0, 1}.

*Proof.* C(i, j) ≤ C(i − 1, j) + 1 by the I step. Conversely, delete x_i from an optimal alignment of x₁…x_i with y₁…y_j: if x_i stood against a gap, its column goes (−1); if it stood against y_k, that column becomes a gap column (at most +1). So C(i − 1, j) ≤ C(i, j) + 1. The horizontal case is symmetric. For the diagonal, C(i, j) ≤ C(i − 1, j − 1) + 1 by the pair step; and deleting both x_i and y_j from an optimal alignment of the prefixes never raises its cost: if they are paired, their column goes; otherwise one of them, say x_i, stands in the last column against a gap (−1), and y_j either stands against a gap (−1) or against some x_l, whose column becomes a gap column (at most +1). ∎

So a column is its first value and its vertical differences. Write them as two bit vectors: Pv with bit i − 1 set when C(i, j) − C(i − 1, j) = +1, and Mv when it is −1; Eq(c) has bit i − 1 set when x_i = c (the *match vectors* of the query, one per symbol, built once).

**Proposition 2.2 (the column step; Myers 1999, in the formulation of Hyyrö 2001).** Let Pv, Mv be the vertical differences of column j − 1 and Eq = Eq(y_j). With Xv = Eq ∨ Mv and

Xh = (((Eq ∧ Pv) + Pv) ⊕ Pv) ∨ Eq,

the horizontal differences of column j are Ph = Mv ∨ ¬(Xh ∨ Pv) (+1) and Mh = Pv ∧ Xh (−1), and its vertical differences are, after shifting Ph and Mh up by one row (bringing in the difference of row 0), Pv′ = Mh ∨ ¬(Xv ∨ Ph) and Mv′ = Ph ∧ Xv.

*The reasoning, which Myers proves.* At one cell, with the differences Δv entering from the left column and Δh entering from the row above, the diagonal difference is 0 exactly when x_i = y_j, or Δv = −1, or Δh = −1; then Δh_out = Δd − Δv and Δv_out = Δd − Δh. Everything is local except Δh, which comes from the row above in the same column: a chain down the column. Xh marks the rows whose Δh is −1: a row is such when it matches, or when the row above had Δh = −1 and this row had Δv = +1 (a run of +1 differences carries a match downwards). The addition (Eq ∧ Pv) + Pv carries a 1 from each matching row with Pv set through the following run of Pv bits, and the exclusive or with Pv leaves the bits the carry passed: one machine addition resolves the whole chain of 64 rows. This is why a column of 64 cells is a dozen word operations.

**Proposition 2.3 (blocks; Myers 1999).** Cut the query into blocks of 64 rows. A block computes its column from its own Pv, Mv and Eq and the horizontal difference hin ∈ {−1, 0, +1} entering at its top row: hin = +1 shifts a 1 into Ph, hin = −1 a 1 into Mh and also sets the lowest bit of Eq (a −1 entering from above makes the top row's diagonal difference 0, as a match does; Hyyrö 2003). The block passes on the horizontal difference of its last row. So a column of a pattern of any length is ⌈m/64⌉ such steps, chained.

`myers_block` is this step; `MyersPattern` keeps the match vectors symbol-major, one word per block, with a row of zeros for the symbols that match nothing. The bits above the last row of the last block hold garbage that is never read: the addition carries upwards only, so no bit depends on a higher one.

## 3. The Band and the Traceback of the Bit Vectors

For a threshold k only the cells that can lie on a path of cost at most k are needed. Let h(i, j) = |(m − i) − (n − j)| in global mode, the least cost of any path from (i, j) to (m, n) (each step changes the diagonal by at most one, and a step that changes it costs one), and h = 0 in semiglobal mode, where every cell of the last row is an end.

**Definition 3.1 (relevant).** A cell is *relevant* when C(i, j) + h(i, j) ≤ k.

**Lemma 3.2.** The optimal predecessor of a relevant cell is relevant.

*Proof.* h is consistent: h(p) ≤ cost(p → q) + h(q) for every edge p → q, because an edge changes the diagonal by at most one and an edge that changes it costs one. If p is the optimal predecessor of q, C(p) + cost(p → q) = C(q), so C(p) + h(p) ≤ C(q) + h(q) ≤ k. ∎

**Proposition 3.3 (the band moves by one row).** If every relevant cell of column j − 1 lies in rows at most R, every relevant cell of column j lies in rows at most R + 1. In global mode, if every relevant cell of column j − 1 lies in rows at least T > 0 and row 0 of column j is not relevant, every relevant cell of column j lies in rows at least T.

*Proof.* Let (i, j) be relevant with i > R + 1. By Proposition 2.1, C(i − 1, j − 1) ≤ C(i, j), and h(i − 1, j − 1) = h(i, j) since the two cells share a diagonal; so (i − 1, j − 1) is relevant, with i − 1 > R: a contradiction. At the top, follow optimal predecessors from a relevant (i, j) with i < T: by Lemma 3.2 they are relevant; one that lies in column j − 1 has a row at most i < T, which is excluded, so they stay in column j and descend to row 0, which is not relevant: a contradiction. ∎

*Remark.* Whether row 0 is relevant in column j depends on j + |m − n + j| ≤ k, a sum that never decreases with j: once row 0 leaves the band it never returns, which is what the top of the band needs.

**Proposition 3.4 (cells outside the band).** Compute the blocks of the band only, taking the cells outside it at upper bounds: below the band, vertical differences of +1 from its last row; above it, a horizontal difference of +1. Then every computed value is at least the true one, and every relevant cell's computed value is exact.

*Proof.* The recurrence is a minimum of sums, monotone in its inputs, so inputs at least the true values give outputs at least the true values, by induction in the order of computation; and the upper bounds are upper bounds by Proposition 2.1. A relevant cell lies in the band, and so do its chain of optimal predecessors (Lemma 3.2, Proposition 3.3); by induction along that chain its computed value is at most the predecessor's computed (exact) value plus the step's cost, which is the true value. ∎

So `myers_pass` keeps, for each column, a range [first, last] of blocks with the score of each block's last row, and:

- adds the block below when the last row of the band, in the previous column, had C + h ≤ k (by Proposition 3.3 nothing deeper can be relevant), initialising it from the upper bound;
- drops a block at the bottom, and in global mode at the top (when row 0 is out), when every cell of it has C + h > k, by the bound C(i) ≥ score − (rows below i) each column and by the exact values every eighth;
- keeps row 0 as part of the band in semiglobal mode, so that an emptied band can grow again from it.

Every cell within k is exact (Proposition 3.4), so the ends of semiglobal matches within k and the global distance, when at most k, are exact; and Ukkonen's doubling of k from 64 (Ukkonen 1985) finds the distance d in O(⌈d/64⌉ n) word operations, since the last threshold is less than 2d. This is the band of Edlib (Šošić and Šikić 2017), whose conditions were tuned by experiment; here each is a consequence of Propositions 3.3 and 3.4.

**The traceback.** Kept for every column, the band's (Pv, Mv) and the scores of the blocks' last rows give any cell's value: C(i, j) is the score of its block's last row less the vertical differences below row i, two popcounts. From (m, n), or from the leftmost end of least distance, the path is read back by choosing a predecessor whose value plus the step's cost is the cell's, a pair first, then I, then D; by Proposition 3.4 such a predecessor is exact and relevant, so the walk never leaves the band. Memory: 20 bytes per block of the band per column.

## 4. Stripes and Lanes

### 4.1 The Striped Order

Farrar (2007) computes the local recurrences (gotoh.hpp; alignment.md §3) a column of the target at a time, on vectors of W lanes. The query is cut into W stripes of S = ⌈m/W⌉ positions; lane l of segment s holds the query position s + l S. Then:

- the diagonal dependency H(i − 1, j − 1) of segment s is segment s − 1 of the previous column in the same lane, and for s = 0 the previous column's last segment shifted up one lane;
- the horizontal dependency E(i, j − 1) is the same segment of the previous column;
- the vertical dependency F(i − 1, j) is the previous segment of the same column, in the same lane, except across the stripe boundaries, where it passes from the last segment of lane l to the first of lane l + 1.

The query profile (Rognes and Seeberg 2000) holds, for each symbol c of the target, the S vectors of s(x_{s + lS}, c), so that the inner loop loads its scores. Only the vertical dependency across lanes is out of order.

**Proposition 4.1 (lazy F; Farrar 2007).** Compute a column with the F entering each lane's first segment set to 0. Then carry F from each lane's last segment into the next lane's first, segment by segment, raising H and the E of the next column where F exceeds them, and stop at the first segment where F ≤ H − (o + e) in every lane. The column is then exact.

*Proof.* The local recurrence floors every value at 0, so 0 is a lower bound of the true F entering a stripe, and the first pass computes lower bounds of every H and F, exact wherever the optimal F-chain does not cross a stripe boundary. The carry computes the missing chain: F(i) = max(H(i − 1) − (o + e), F(i − 1) − e) along the column in its true order. When the carried F satisfies F ≤ H(i) − (o + e) in every lane at segment i, it cannot raise H(i), and the next candidate F − e is below H(i) − (o + e), which the first pass already counted in F(i + 1) with the final H(i); so no later cell can change. The E of the next column is raised with every raised H, which keeps Gotoh's transition from F through H into E (Farrar's update; SSW drops it and forbids that transition). ∎

The loop seldom runs more than once around the lanes: a vertical gap long enough to cross a stripe must pay for itself. Cheap gaps and long insertions make it run longer, and the check drives it through millions of corrections to see them right.

### 4.2 Saturation

Narrow lanes compute more cells per instruction: 64 lanes of 8 bits on AVX-512. SSW's device (Zhao et al. 2013) makes unsigned bytes serve a signed recurrence: the profile is biased by β = max(1, −min s), so that it is never negative, and a cell is H + s computed as (H + (s + β)) − β, each operation saturating; the saturation at 0 of the subtraction is the local floor itself.

**Proposition 4.2.** If the best value computed in 8-bit lanes is below 255 − β, every value is exact.

*Proof.* The only operation that can lose information is the saturating addition, which saturates exactly when H + s + β ≥ 255; its result is then 255 − β after the subtraction, a value of the column, so the best would be at least 255 − β. The subtractions of gap costs saturate at 0, where the local recurrence's maximum with 0 puts them anyway, and the maxima are exact. ∎

When the best reaches 255 − β the pass is redone in 16-bit signed lanes (floored at 0 by a maximum), and past 32767 − max s in 32-bit lanes; the stride of the lanes is the build's (`calx::cpu::Lanes`), so the same code is SSE2, AVX2 or AVX-512.

### 4.3 The End, the Start and the Path

The kernel keeps the best value and the first cell holding it in row-major order (the smallest query end, then the smallest target end), gotoh.hpp's rule, from the columns whose maximum reaches the best.

**Proposition 4.3.** Let (i, j) be that cell and b its value. The best local alignment of the reversed prefixes (x₁…x_i)ᴿ and (y₁…y_j)ᴿ has score b, and every alignment of score b there ends, read forwards, at (i, j). The global optimum of the region between its start and (i, j) is b.

*Proof.* A local alignment of the prefixes is a local alignment of x and y, so it scores at most b, and the optimal one ending at (i, j) scores b. An alignment of score b of the prefixes that ends, forwards, at (i′, j′) with i′ ≤ i and j′ ≤ j makes H(i′, j′) ≥ b, so (i′, j′) is a cell of value b; were it ≠ (i, j) it would come first in row-major order. For the region: its global alignments are local alignments of x and y (at most b), and the local one found is a global alignment of the region (b). ∎

So `striped_local` runs the kernel on the reversed prefixes for the start, and a1's `banded_global` on the region for the path: an alignment in O(⌈m/W⌉ n) plus the region's cost.

### 4.4 One Target per Lane

`inter.hpp` puts W different problems in the lanes, the query against W targets, one cell (i, j) of each in a vector (Rognes 2011, SWIPE; Wozniak 1997): every dependency stays in its lane and there is no correction. The scores, s(x_i, y⁽ˡ⁾_j) with a different target base in each lane, are gathered once per column into a vector for each query symbol (σW scalar writes against the m vector rows that use them). Targets are batched by length, and past its end a target scores −16384 a base, which can only lower H. A lane whose best reaches the 16-bit limit is rescored by the striped kernel.

## 5. Wavefronts: the Tropical Sum in the Order of Its Values

With penalties (match 0, mismatch x, a gap of length L costing o + L e), the optimal alignment is a shortest path in the product graph of alignment.md §3. Dijkstra's algorithm would visit the vertices in increasing order of their distance; the wavefront algorithm (Marco-Sola et al. 2021) does so, and represents the visited set by its frontier on each diagonal.

### 5.1 Furthest Points and Dominance

Let k = h − v be a diagonal (h the target offset, v = h − k the query offset). Write F(v, h) for the least penalty of a path from the start to (v, h) and B(v, h) for the least penalty of a path from (v, h) to the end.

**Proposition 5.1 (diagonal monotonicity).** B(v, h) ≥ B(v + 1, h + 1), and F(v, h) ≤ F(v + 1, h + 1), whenever both points lie in the matrix.

*Proof.* Let P be a path from (v, h) to the end; build P′ from (v + 1, h + 1), no dearer, that agrees with P from some point on. If P begins with a pair, P′ is the rest of P. Otherwise suppose P begins with a run D^L (L ≥ 1) of target bases against gaps, after which it first advances in the query: either by an I step from (v, h + L) to (v + 1, h + L), or by a pair from (v, h + L) to (v + 1, h + L + 1). In the first case P′ takes D^{L−1} along row v + 1 to (v + 1, h + L) and follows P from there: it drops one column of P's D run and one of the I run that follows, each saving e, or a whole o + e when the run vanishes; the remainder of P's I run, now begun afresh, pays the opening P paid once already. In the second case P′ takes D^L to (v + 1, h + L + 1) and follows P from there: it drops the pair, and a D run of P that follows now extends P′'s instead of opening. A run of I first is symmetric. The second inequality is the first on the reversed sequences. ∎

So on each diagonal the points reached at penalty at most s form an initial segment, and its furthest point stands for it: from the furthest point the rest of the alignment is no dearer than from any nearer one. The wavefront M̃(s, k) of a penalty s is that furthest offset, and the algorithm need keep nothing else.

### 5.2 The Recurrences of the Wavefronts

With Ĩ, D̃ the furthest points reached by a path ending in a gap of each kind (Marco-Sola et al. 2021, the gap-affine recurrences, with I a query base against a gap, from diagonal k + 1, and D a target base, from k − 1):

- Ĩ(s, k) = max(M̃(s − o − e, k + 1), Ĩ(s − e, k + 1)),
- D̃(s, k) = max(M̃(s − o − e, k − 1), D̃(s − e, k − 1)) + 1,
- M̃(s, k) = extend(max(M̃(s − x, k) + 1, Ĩ(s, k), D̃(s, k))),

where extend follows the diagonal while the bases match (matches cost 0: the zero-weight edges of the graph, collapsed). Each candidate is dropped when it leaves the matrix. The first penalty s at which M̃(s, m − n) = m (with v = n) is the optimum (Marco-Sola et al. 2021, from the dominance of Proposition 5.1). The extension compares eight bytes at a time: the trailing zeros of the exclusive or of two words, divided by eight, is the first difference; codes of N are replaced by bytes that differ between the two sequences, and sentinels after the ends stop it.

The work: a wavefront of penalty s spans at most about 2s/e diagonals, so the offsets computed number O(s²/e), and the extension walks each diagonal once in all, O(n + m) per diagonal at worst. In practice, on similar sequences, the work is O(n + s²) against Gotoh's O(nm) (§7). The traceback keeps every wavefront, O(s²/e) offsets; the score alone keeps the last max(x, o + e) + 1, O(s).

### 5.3 Penalties and Scores

**Proposition 5.2 (Eizenga and Paten 2022).** Let a scoring have match α, mismatch −β, and gaps −(q + L r). For every global alignment of lengths m and n,

S = (α (m + n) − P) / 2, with P the penalty under x = 2(α + β), o = 2q, e = 2r + α.

*Proof.* An alignment with a matches, b mismatches, g gaps and L gapped bases has 2a + 2b + L = m + n. Then S = α a − β b − q g − r L = α(m + n − L)/2 − (α + β) b − q g − r L, and 2S = α(m + n) − [2(α + β) b + 2q g + (2r + α) L]. ∎

The two orders are the same on the global alignments of one pair, so the WFA on the penalties (reduced by their greatest common divisor) computes Gotoh's global optimum under the scoring (`penalty_form`, `wavefront_global`). With N in the sequences the identity needs the matrix to score N as a mismatch.

### 5.4 The Ends of a Sub-Problem

BiWFA (§5.5) splits an alignment at a point inside a gap, and its halves must then begin or end inside that gap. Two kinds of end serve:

- *open*: a gap of one kind is open across the end. At the beginning, the path's first run of that kind pays no opening; at the end, its last run of that kind is refunded its opening. Neither forbids anything: the path may also not continue the gap.
- *forced*: the path's last step must be a step of the gap.

**Proposition 5.3.** The cost to the end B is non-increasing along a diagonal (Proposition 5.1) when the end is plain or open, but not when it is forced.

*Proof.* The construction of Proposition 5.1 keeps the tail of P, so it keeps P's last run and its refund; it may shorten P's first run, where an open beginning's discount applies equally. For a forced end it can consume the forced step: with x = A, y = A and the last step forced to be D, the only path is I then D, of penalty 2(o + e), from (0, 0); but from (1, 1), the furthest point of diagonal 0 at penalty 0, the end is the empty path, which does not end with D. The WFA keeps (1, 1) and loses (0, 0), and never finds the path. ∎

A forced end is therefore peeled: a path that must end with a D step is a path to (n, m − 1) with a D gap open at its end, followed by that step at o + e; when the path before it ended with a run of D, the refund makes the step cost e, an extension. So every sub-problem has plain or open ends, toward which the furthest points dominate, and the WFA solves it exactly (`wf_solve`; the check compares all nine combinations of ends with a three-state recurrence on 8100 sub-problems).

### 5.5 Two Directions Meeting

BiWFA (Marco-Sola et al. 2023) runs the WFA forwards from the start and backwards, on the reversed sequences, from the end (whose open gap is the reverse direction's open beginning), in turn, keeping only the last max(x, o + e) + 1 wavefronts of each. Forward diagonal k is reverse diagonal (m − n) − k, and a forward offset h_f and a reverse offset h_r on it *overlap* when h_f + h_r ≥ m.

**Proposition 5.4 (soundness).** An overlap of the forward wavefront of penalty a with the reverse one of penalty b, in M, gives a path of penalty at most a + b; in I or D, at most a + b − o; and the forward point is on such a path, with the cost from it to the end at most b (or b − o).

*Proof.* The forward point X_f and the reverse point X_r (in forward coordinates) lie on one diagonal with X_r ≤ X_f. By Proposition 5.1 (and 5.3, the ends being plain or open), the cost from X_f to the end is at most the cost from X_r, which is at most b. In a gap, the forward path ends with a run of the gap and the reverse path begins with one, each having paid its opening (or none, if the run is the open gap of its end); joined at X_f inside one gap they pay it once, hence − o. ∎

**Proposition 5.5 (completeness).** Let s* be the optimum and p = max(x, o + e). Comparing each new wavefront with the other direction's last p + 1, the least sum found is s* by the time the two penalties sum to s* + o + p.

*Proof (sketch).* Take an optimal path and, when the forward penalty is a_f, the last point of the path whose prefix costs at most a_f; its prefix costs c > a_f − p, since one step costs at most p. Its suffix costs s* − c, or s* + o − c in a gap (counted with its own opening), and the reverse wavefront of that penalty reaches the point (Proposition 5.1). The two wavefronts lie within p of the two current penalties once these sum to about s* + p (s* + o + p in a gap), and the pair was compared when the later of the two was computed. ∎

The point found splits the problem into two of known penalties, aligned recursively; halves of penalty at most 256 by the WFA with its traceback. Each level keeps O(s) offsets, so the memory is O(s) against the WFA's O(s²): 0.24 MB against 14.5 MB for two 30 kb sequences at 1% (§7). The time is about two to three times the WFA's.

### 5.6 The Adaptive Heuristic

WF-adaptive (Marco-Sola et al. 2021, the heuristic; WFA2-lib's defaults) drops, after each wavefront of at least 10 diagonals, those whose remaining distance max(n − v, m − h) exceeds the least by more than 50. It is a heuristic: the alignment is valid and its penalty an upper bound, which on similar sequences is the optimum (the check finds it on all its pairs at 1%). It is off by default.

## 6. Partial Orders

### 6.1 The Graph and Its Columns

A partial-order graph (Lee, Grasso and Sharlow 2002) is a DAG of bases whose edges join the consecutive bases of the sequences that pass through it. Nodes that alignments placed in one column with different bases are *aligned*, and a column is a set of nodes with distinct bases closed under this relation.

**Proposition 6.1.** Adding a sequence along an alignment to the graph keeps the column graph (the quotient of the DAG by the columns) acyclic.

*Proof.* The alignment's path visits nodes in increasing order of the graph (a pair or a deletion moves along an edge; an insertion stays). Each base of the sequence becomes a node in the column of the node it is paired with (fused, or a new member) or a new column of its own, so the sequence's nodes lie in columns of increasing order, with the new columns between their neighbours; the new edges join consecutive ones forwards. A linear order of the columns therefore extends to one of the new columns, and the quotient stays acyclic. ∎

`PoaGraph::sort` keeps such an order by Kahn's algorithm on the columns, the ready column of the smallest node first, its nodes adjacent: the order of the graph and the columns of its multiple alignment at once.

### 6.2 The Edit Graph of a DAG

The edit graph of a DAG D and a word x has the vertices (u, j) for the nodes u and positions j, with a pair edge from (p, j − 1) to (u, j) for each edge p → u of D, an I edge from (u, j − 1) and a D edge from (p, j). It is the product of D with the path of x, itself a DAG, so the sum over its paths is computed in a topological order of D × {0, …, n} (alignment.md, Proposition 2.3), and with Gotoh's three states for affine gaps. A virtual start S precedes the sources. Global: from S to a sink of D with x consumed; local: anywhere, with the floor 0. The check compares it with the best over every source-to-sink path of small graphs (brute force) and, on a chain, with Gotoh.

### 6.3 The Heaviest Bundle and the Multiple Alignment

Lee (2003) reads a consensus from the graph as a heavy path: each node keeps its heaviest incoming edge (ties to the predecessor of higher score), its score the edge's weight plus that predecessor's; the path back from the best node is the consensus, here continued forward along the heaviest outgoing edges to a sink. It is greedy, not the heaviest path, and that is its virtue: an edge of weight 1 (one read's error) never wins against the sequence's own path, whose edges weigh nearly the number of reads. On 20 reads of a 2 kb sequence at 8% error the consensus is within 2 edits of the sequence (`check_poa`).

The multiple alignment the graph holds has one row for each sequence, its bases in the columns of its nodes and gaps elsewhere; by Proposition 6.1 each row's columns increase. On a simulated family of 8 (600 bp; per branch 2% substitutions and 0.5% indels) it recovers 98.6% of the true pairs; at 4% and 1% per branch, 92%.

## 7. The Kernels Measured

On an Intel i7-11370H (AVX-512, one thread), from `align/benchmarks/bench_kernels.cpp`, 2026-10-24. Cell updates per second (GCUPS) for a 1 kb query against a 100 kb target:

| Kernel | Time | GCUPS |
|---|---|---|
| Gotoh, scalar (`gotoh_end`) | 0.43 s | 0.23 |
| striped, 8-bit lanes (score below 255) | 6.1 ms | 16.4 |
| striped, 8 bits then 16 (score 1208) | 15 ms | 6.6 |
| inter-sequence, 16-bit, 320 targets of 10 kb | 0.38 s | 8.5 |
| Myers, every end within 150 edits | 5.4 ms | 18.5 (32% of the cells in the band) |

Global alignment of 10 kb pairs, penalties (4, 6, 2):

| Divergence | Penalty | Gotoh O(nm) | WFA with path | WFA score | BiWFA | Memory, WFA / BiWFA |
|---|---|---|---|---|---|---|
| 1% | 752 | 0.50 s | 0.8 ms | 0.8 ms | 3.3 ms | 1.7 MB / 0.11 MB |
| 5% | 3978 | 0.48 s | 50 ms | 21 ms | 96 ms | 47.5 MB / 0.43 MB |
| 10% | 7104 | 0.49 s | 191 ms | 60 ms | 293 ms | 151.5 MB / 0.77 MB |

*Remark.* The wavefronts' advantage is the ratio of mn to the number of offsets they visit, which falls as the penalty grows: 600 times at 1%, 10 at 5%, 2.6 at 10% with the traceback. Past a few percent of divergence on long sequences, banded Gotoh or the striped kernel on the region compete.

## 8. Where the Code Is

| Concept | Code |
|---|---|
| Match vectors, the column step, blocks | `MyersPattern`, `detail::myers_block` (`bitparallel.hpp`) |
| The band, the doubling, the traceback | `detail::myers_pass`, `myers_distance`, `myers_search`, `myers_align` |
| The striped order, lazy F, saturation | `detail::striped_pass`, `StripedProfile` (`striped.hpp`) |
| The start from the reversed prefixes, the path | `striped_local` |
| One target per lane | `inter_scores` (`inter.hpp`) |
| Wavefronts, extension, traceback | `detail::WfEngine`, `detail::wf_run`, `detail::wf_backtrace`, `wavefront`, `wavefront_penalty` (`wavefront.hpp`) |
| Penalties and scores | `penalty_form`, `score_from_penalty`, `wavefront_global` |
| Ends of sub-problems, BiWFA | `detail::wf_solve`, `detail::wf_breakpoint`, `detail::wf_bidirectional`, `biwfa` |
| The record for drawing | `WavefrontRecord`, `MyersStats`, `StripedStats` |
| Partial orders | `PoaGraph`: `align`, `add`, `consensus`, `msa` (`poa.hpp`) |

## References

- Eizenga, J. M. and Paten, B. (2022). Improving the time and space complexity of the WFA algorithm and generalizing its scoring. bioRxiv 2022.01.12.476087.
- Farrar, M. (2007). Striped Smith–Waterman speeds database searches six times over other SIMD implementations. *Bioinformatics* 23, 156–161.
- Gotoh, O. (1982). An improved algorithm for matching biological sequences. *J. Mol. Biol.* 162, 705–708.
- Hyyrö, H. (2001). Explaining and extending the bit-parallel approximate string matching algorithm of Myers. Tech. Rep. A-2001-10, University of Tampere.
- Hyyrö, H. (2003). A bit-vector algorithm for computing Levenshtein and Damerau edit distances. *Nordic J. Computing* 10, 29–39.
- Lee, C., Grasso, C. and Sharlow, M. F. (2002). Multiple sequence alignment using partial order graphs. *Bioinformatics* 18, 452–464.
- Lee, C. (2003). Generating consensus sequences from partial order multiple sequence alignment graphs. *Bioinformatics* 19, 999–1008.
- Marco-Sola, S., Moure, J. C., Moreto, M. and Espinosa, A. (2021). Fast gap-affine pairwise alignment using the wavefront algorithm. *Bioinformatics* 37, 456–463.
- Marco-Sola, S., Eizenga, J. M., Guarracino, A., Paten, B., Garrison, E. and Moreto, M. (2023). Optimal gap-affine alignment in O(s) space. *Bioinformatics* 39, btad074.
- Myers, G. (1999). A fast bit-vector algorithm for approximate string matching based on dynamic programming. *J. ACM* 46, 395–415.
- Rognes, T. and Seeberg, E. (2000). Six-fold speed-up of Smith–Waterman sequence database searches using parallel processing on common microprocessors. *Bioinformatics* 16, 699–706.
- Rognes, T. (2011). Faster Smith–Waterman database searches with inter-sequence SIMD parallelisation. *BMC Bioinformatics* 12, 221.
- Šošić, M. and Šikić, M. (2017). Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. *Bioinformatics* 33, 1394–1395.
- Ukkonen, E. (1985). Algorithms for approximate string matching. *Information and Control* 64, 100–118.
- Wozniak, A. (1997). Using video-oriented instructions to speed up sequence comparison. *CABIOS* 13, 145–150.
- Zhao, M., Lee, W.-P., Garrison, E. P. and Marth, G. T. (2013). SSW library: an SIMD Smith–Waterman C/C++ library for use in genomic applications. *PLoS ONE* 8, e82138.
