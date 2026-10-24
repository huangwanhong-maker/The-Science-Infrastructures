// Progressive multiple alignment along a guide tree, and its refinement: the three stages of
// MUSCLE (Edgar 2004, Nucleic Acids Research 32, 1792–1797; BMC Bioinformatics 5, 113) on the
// profile alignment of CLUSTAL W (Thompson, Higgins and Gibson 1994, Nucleic Acids Research 22,
// 4673–4680), after Feng and Doolittle (1987, Journal of Molecular Evolution 25, 351–360) and
// Higgins and Sharp (1988, Gene 73, 237–244).
//
// 1. A tree of k-mer distances by UPGMA (guide_tree.hpp); the sequences are aligned at its
//    nodes, deepest first, each node's alignment the profile alignment of its children's, and
//    "once a gap, always a gap": a node never changes what its children aligned.
// 2. The tree again, by UPGMA of the corrected distances of the first alignment (Kimura for
//    protein, Jukes–Cantor for DNA), and the progressive alignment along it when its clusters
//    differ; the second alignment is kept when its SP score is not lower. (MUSCLE re-aligns only
//    the subtrees whose order changed; here the whole tree is re-aligned, since the weights
//    change everywhere with the branch lengths.)
// 3. Refinement by tree-dependent restricted partitioning: each edge of the tree, in order of
//    decreasing distance from the root, splits the rows in two; the two sub-alignments (their
//    columns of gaps only removed) are re-aligned as profiles, and the result is kept if it
//    raises the SP score (scores.hpp). Only the pairs across the split can change, since a
//    sub-alignment's pairs project to themselves, so the change is computed on them alone. A few
//    rounds, until one changes nothing.
//
// Profiles. The profile of an alignment holds, in each column, the weighted frequency f(a) of
// each code a below the matrix's size and the weighted occupancy o = Σ_a f(a) (the fraction of
// rows with a residue there, 1 − the gap fraction), the weights being CLUSTAL W's from the tree
// (or 1). Two profiles x (m columns) and y (n columns) are aligned by Gotoh's recurrences
// (align/gotoh.hpp) with sum-of-pairs costs, each normalised per pair of rows:
//
//   a column pair     s(i, j) = Σ_a Σ_b f_x,i(a) f_y,j(b) S(a, b)
//                     (the expected score of the residue pairs across the two columns; x is kept
//                     sparse, y as the vectors S f_y,j, so a cell costs the residues of x's column)
//   a column of x against a new gap in y, inserted after y's column j:
//                     extension e · o_x,i        (each residue of x pays one extension per row of y)
//                     opening   (g/2) · o_x,i · o_y,j at its first column and
//                               (g/2) · o_x,i · o_y,j+1 at its last
//   and symmetrically for a column of y against a gap in x, with o = 1 beyond either end. A gap
//   opens for the pairs whose rows had residues on both sides of where it is inserted; the
//   fractions o_y,j and o_y,j+1 count those that do (half of the opening at each end, the
//   "gap open" and "gap close" of MUSCLE), so that a gap is cheap where the other profile already
//   has gaps, and two single sequences are aligned exactly as Gotoh aligns them (a gap of length
//   L costing g + L e). With terminal_open false, a gap at either end of a profile pays no
//   opening. This is an approximation of the exact sum-of-pairs alignment of two alignments
//   with affine gaps, which is NP-complete (Kececioglu and Starrett 2004, Proceedings of RECOMB); the
//   objective the refinement keeps improvements of is the exact SP score.
//
// Ties are broken as Gotoh's (a pair, then a column of x, then of y; opening before extending),
// so the result is a function of the input. Parallel: the k-mer and corrected distances by rows,
// the nodes of the tree as the tasks of a calx::scheduler::TaskGraph (siblings at once), and the
// SP score by rows; every task writes its own slot, so the alignment is the same on any number of
// threads.
//
// Complexity: the distances O(n² L); each node O(m n r) for profiles of m and n columns and r
// residues per column of x (at most the alphabet's size), with (m + 1)(n + 1) bytes of
// traceback; n − 1 nodes, so O(n L²) for n sequences of length L. Refinement: per edge one
// profile alignment and the SP across the split, O(L² + n² L).

#pragma once

#include <genomics/align/scoring.hpp>
#include <genomics/align/types.hpp>
#include <genomics/msa/guide_tree.hpp>
#include <genomics/msa/msa.hpp>
#include <genomics/msa/scores.hpp>

#include <calx/scheduler/task_graph.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace genomics::msa {

// ---- profiles ------------------------------------------------------------------------------------
struct Profile {
  std::size_t columns = 0;
  int sigma = 0;                     // codes 0 .. sigma − 1 (the matrix's size)
  std::vector<double> occupancy;     // the weighted fraction of rows with a residue
  std::vector<std::uint32_t> begin;  // columns + 1: column c's residues are [begin[c], begin[c + 1])
  std::vector<Code> code;            // the codes present in each column
  std::vector<double> freq;          // their weighted frequencies
};

// The profile of an alignment under row weights (all 1 if `weights` is empty).
inline Profile make_profile(const Msa& m, const std::vector<double>& weights, int sigma) {
  require_rectangular(m, "genomics::msa::make_profile");
  if (!weights.empty() && weights.size() != m.size()) throw std::invalid_argument("genomics::msa::make_profile: one weight per row");
  Profile p;
  p.columns = m.columns();
  p.sigma = sigma;
  std::vector<double> w = weights.empty() ? std::vector<double>(m.size(), 1.0) : weights;
  double total = std::accumulate(w.begin(), w.end(), 0.0);
  if (!(total > 0)) { std::fill(w.begin(), w.end(), 1.0); total = static_cast<double>(m.size()); }
  for (double& x : w) x /= total;
  std::vector<double> f(static_cast<std::size_t>(sigma));
  p.begin.reserve(p.columns + 1);
  p.begin.push_back(0);
  p.occupancy.reserve(p.columns);
  for (std::size_t c = 0; c < p.columns; c++) {
    std::fill(f.begin(), f.end(), 0.0);
    double occ = 0;
    for (std::size_t r = 0; r < m.size(); r++) {
      const Code x = m.rows[r][c];
      if (x == GAP) continue;
      if (x >= sigma) throw std::invalid_argument("genomics::msa::make_profile: a code outside the matrix");
      f[x] += w[r];
      occ += w[r];
    }
    for (int a = 0; a < sigma; a++)
      if (f[static_cast<std::size_t>(a)] > 0) {
        p.code.push_back(static_cast<Code>(a));
        p.freq.push_back(f[static_cast<std::size_t>(a)]);
      }
    p.occupancy.push_back(occ);
    p.begin.push_back(static_cast<std::uint32_t>(p.code.size()));
  }
  return p;
}

// ---- profile–profile alignment ------------------------------------------------------------------------------------
struct ProfileAlignment {
  double score = 0;
  std::string path;  // 'M' a column of each; 'I' a column of x against a gap; 'D' a column of y against a gap
};

inline ProfileAlignment align_profiles(const Profile& x, const Profile& y, const align::Scoring& s, bool terminal_open = false) {
  const int sigma = s.matrix.size;
  if (x.sigma != sigma || y.sigma != sigma) throw std::invalid_argument("genomics::msa::align_profiles: profiles and matrix of different alphabets");
  if (s.gap_open < 0 || s.gap_extend < 0) throw std::invalid_argument("genomics::msa::align_profiles: negative gap costs");
  const std::size_t m = x.columns, n = y.columns, W = n + 1;
  const std::size_t S = static_cast<std::size_t>(sigma);
  // (S f_y)(a) for each column j of y, stored by a and then j, so that a row of pair scores is a
  // sum of contiguous vectors, one for each residue of x's column.
  std::vector<double> sf(S * W, 0.0);
  for (std::size_t j = 1; j <= n; j++)
    for (std::uint32_t k = y.begin[j - 1]; k < y.begin[j]; k++) {
      const double f = y.freq[k];
      const std::size_t b = y.code[k];
      for (std::size_t a = 0; a < S; a++) sf[a * W + j] += f * s.matrix.score[a * S + b];
    }
  const double half = 0.5 * s.gap_open, e = s.gap_extend, NEG = -std::numeric_limits<double>::infinity();
  // Occupancies by 1-based column, 1 beyond either end; the factors of the openings and closings
  // of gaps in x (I steps, by y's column) with the ends zeroed when they are free.
  std::vector<double> oy(n + 2, 1.0), open_i(W), close_i(W);
  for (std::size_t j = 1; j <= n; j++) oy[j] = y.occupancy[j - 1];
  for (std::size_t j = 0; j <= n; j++) {
    const bool end = (j == 0 || j == n) && !terminal_open;
    open_i[j] = end ? 0.0 : half * oy[j];
    close_i[j] = end ? 0.0 : half * oy[j + 1];
  }
  const auto ox = [&](std::size_t c) { return c == 0 || c > m ? 1.0 : x.occupancy[c - 1]; };
  std::vector<double> Hp(W, NEG), Hc(W, NEG), Ip(W, NEG), Ic(W, NEG), pair(W, 0.0);
  std::vector<std::uint8_t> dir((m + 1) * W, 0);
  enum : std::uint8_t { FROM_PAIR = 0, FROM_I = 1, FROM_D = 2, I_EXTENDS = 4, D_EXTENDS = 8 };
  for (std::size_t i = 0; i <= m; i++) {
    const double oxi = ox(i);
    const bool end_row = (i == 0 || i == m) && !terminal_open;  // a gap in y here touches an end of x
    const double open_d = end_row ? 0.0 : half * oxi, close_d = end_row ? 0.0 : half * ox(i + 1);
    std::uint8_t* t = dir.data() + i * W;
    if (i == 0) {
      // The first row: columns of y against gaps only.
      Hc[0] = 0;
      Ic[0] = NEG;
      t[0] = FROM_PAIR;
      double hl = 0, dl = NEG;
      for (std::size_t j = 1; j <= n; j++) {
        const double a = hl - open_d * oy[j], b = dl;
        const bool dext = b > a;
        const double dv = (dext ? b : a) - e * oy[j];
        Hc[j] = dv - close_d * oy[j];
        Ic[j] = NEG;
        t[j] = static_cast<std::uint8_t>(FROM_D | (dext ? D_EXTENDS : 0));
        hl = Hc[j];
        dl = dv;
      }
      std::swap(Hp, Hc);
      std::swap(Ip, Ic);
      continue;
    }
    // The pair scores of the row: Σ over x's residues of f · (S f_y)(a), contiguous in j.
    std::fill(pair.begin(), pair.end(), 0.0);
    for (std::uint32_t k = x.begin[i - 1]; k < x.begin[i]; k++) {
      const double f = x.freq[k];
      const double* v = sf.data() + static_cast<std::size_t>(x.code[k]) * W;
      for (std::size_t j = 1; j <= n; j++) pair[j] += f * v[j];
    }
    // I(i, j) for the whole row: it needs only the row above.
    const double ext_i = e * oxi;
    for (std::size_t j = 0; j <= n; j++) {
      const double a = Hp[j] - oxi * open_i[j], b = Ip[j];
      const bool iext = b > a;
      Ic[j] = (iext ? b : a) - ext_i;
      t[j] = iext ? I_EXTENDS : 0;
    }
    // H and D along the row.
    Hc[0] = Ic[0] - oxi * close_i[0];
    t[0] |= FROM_I;
    double hl = Hc[0], dl = NEG;
    for (std::size_t j = 1; j <= n; j++) {
      const double a = hl - open_d * oy[j], b = dl;
      const bool dext = b > a;
      const double dv = (dext ? b : a) - e * oy[j];
      double h = Hp[j - 1] + pair[j];
      std::uint8_t src = FROM_PAIR;
      const double ci = Ic[j] - oxi * close_i[j];
      if (ci > h) { h = ci; src = FROM_I; }
      const double cd = dv - close_d * oy[j];
      if (cd > h) { h = cd; src = FROM_D; }
      Hc[j] = h;
      t[j] = static_cast<std::uint8_t>(t[j] | src | (dext ? D_EXTENDS : 0));
      hl = h;
      dl = dv;
    }
    std::swap(Hp, Hc);
    std::swap(Ip, Ic);
  }
  ProfileAlignment out;
  out.score = Hp[n];
  std::size_t i = m, j = n;
  enum { IN_H, IN_I, IN_D } state = IN_H;
  std::string rev;
  rev.reserve(m + n);
  while (i > 0 || j > 0) {
    const std::uint8_t d = dir[i * W + j];
    if (state == IN_H) {
      const std::uint8_t src = d & 3;
      if (src == FROM_PAIR) { rev += 'M'; i--; j--; }
      else state = src == FROM_I ? IN_I : IN_D;
      continue;
    }
    if (state == IN_I) {
      rev += 'I';
      state = (d & I_EXTENDS) ? IN_I : IN_H;
      i--;
    } else {
      rev += 'D';
      state = (d & D_EXTENDS) ? IN_D : IN_H;
      j--;
    }
  }
  out.path.assign(rev.rbegin(), rev.rend());
  return out;
}

// The alignment of x's rows and y's rows that a path of align_profiles makes: x's rows first.
inline Msa merge(const Msa& x, const Msa& y, const std::string& path) {
  Msa out;
  out.names = x.names;
  out.names.insert(out.names.end(), y.names.begin(), y.names.end());
  out.rows.reserve(x.size() + y.size());
  for (const auto& r : x.rows) {
    Codes row(path.size(), GAP);
    std::size_t c = 0;
    for (std::size_t k = 0; k < path.size(); k++)
      if (path[k] != 'D') row[k] = r[c++];
    out.rows.push_back(std::move(row));
  }
  for (const auto& r : y.rows) {
    Codes row(path.size(), GAP);
    std::size_t c = 0;
    for (std::size_t k = 0; k < path.size(); k++)
      if (path[k] != 'I') row[k] = r[c++];
    out.rows.push_back(std::move(row));
  }
  return out;
}

namespace detail {

inline std::vector<double> subset(const std::vector<double>& w, const std::vector<std::size_t>& ids) {
  std::vector<double> out;
  out.reserve(ids.size());
  for (std::size_t i : ids) out.push_back(w.empty() ? 1.0 : w[i]);
  return out;
}

// The alignment `m`, whose row k is the input's ids[k], with its rows in the input's order.
inline Msa in_input_order(const Msa& m, const std::vector<std::size_t>& ids) {
  Msa out;
  out.names.resize(m.size());
  out.rows.resize(m.size());
  for (std::size_t k = 0; k < ids.size(); k++) {
    out.names[ids[k]] = m.names[k];
    out.rows[ids[k]] = m.rows[k];
  }
  return out;
}

// The clusters of a dendrogram, each a sorted list of leaves, sorted.
inline std::vector<std::vector<std::size_t>> clusters(const Dendrogram& t) {
  std::vector<std::vector<std::size_t>> out;
  for (std::size_t k = 0; k < t.merges.size(); k++) {
    auto c = t.leaves_under(t.leaves + k);
    std::sort(c.begin(), c.end());
    out.push_back(std::move(c));
  }
  std::sort(out.begin(), out.end());
  return out;
}

inline void validate_input(const std::vector<Codes>& seqs, const std::vector<std::string>& names, const align::Scoring& s, const char* who) {
  if (!names.empty() && names.size() != seqs.size()) throw std::invalid_argument(std::string(who) + ": one name per sequence");
  for (const auto& x : seqs)
    for (Code c : x)
      if (c >= s.matrix.size) throw std::invalid_argument(std::string(who) + ": a code outside the matrix (a gap in an unaligned sequence?)");
}

inline std::vector<std::string> names_or_default(const std::vector<std::string>& names, std::size_t n) {
  if (!names.empty()) return names;
  std::vector<std::string> out;
  for (std::size_t i = 0; i < n; i++) out.push_back("s" + std::to_string(i + 1));
  return out;
}

}  // namespace detail

// The progressive alignment of the sequences along a tree whose leaves are they, under row
// weights (empty: all 1). Rows in the input's order.
inline Msa align_along(const std::vector<Codes>& seqs, const std::vector<std::string>& names_in, const Dendrogram& tree, const std::vector<double>& weights,
                       const align::Scoring& s, bool terminal_open, calx::scheduler::ThreadPool& pool) {
  detail::validate_input(seqs, names_in, s, "genomics::msa::align_along");
  const std::size_t n = seqs.size();
  const auto names = detail::names_or_default(names_in, n);
  if (tree.leaves != n) throw std::invalid_argument("genomics::msa::align_along: the tree's leaves are not the sequences");
  if (n == 0) return {};
  struct Group {
    Msa msa;
    std::vector<std::size_t> ids;
  };
  std::vector<Group> node(tree.nodes());
  for (std::size_t i = 0; i < n; i++) node[i] = {single(seqs[i], names[i]), {i}};
  if (n == 1) return node[0].msa;
  calx::scheduler::TaskGraph graph;
  std::vector<std::size_t> task_of(tree.nodes(), SIZE_MAX);
  for (std::size_t k = 0; k < tree.merges.size(); k++) {
    const Merge mg = tree.merges[k];
    std::vector<std::size_t> after;
    if (task_of[mg.a] != SIZE_MAX) after.push_back(task_of[mg.a]);
    if (task_of[mg.b] != SIZE_MAX) after.push_back(task_of[mg.b]);
    task_of[n + k] = graph.add(
        "node " + std::to_string(n + k),
        [&, mg, k](std::size_t) {
          Group& X = node[mg.a];
          Group& Y = node[mg.b];
          const Profile px = make_profile(X.msa, detail::subset(weights, X.ids), s.matrix.size);
          const Profile py = make_profile(Y.msa, detail::subset(weights, Y.ids), s.matrix.size);
          const ProfileAlignment pa = align_profiles(px, py, s, terminal_open);
          Group& Z = node[n + k];
          Z.msa = merge(X.msa, Y.msa, pa.path);
          Z.ids = X.ids;
          Z.ids.insert(Z.ids.end(), Y.ids.begin(), Y.ids.end());
          X = Group{};
          Y = Group{};
        },
        std::move(after), static_cast<double>(mg.size));
  }
  graph.run(pool);
  const Group& root = node[tree.root()];
  return detail::in_input_order(root.msa, root.ids);
}

// ---- refinement ------------------------------------------------------------------------------------
struct Refinement {
  std::size_t rounds = 0, trials = 0, accepted = 0;
  std::int64_t gain = 0;  // the rise of the SP score
};

// Tree-dependent restricted partitioning: up to `rounds` passes over the edges of the tree.
inline Msa refine(Msa msa, const Dendrogram& tree, const std::vector<double>& weights, const align::Scoring& s, bool terminal_open, int rounds,
                  calx::scheduler::ThreadPool& pool, Refinement* report = nullptr) {
  require_rectangular(msa, "genomics::msa::refine");
  const std::size_t n = msa.size();
  Refinement rep;
  if (tree.leaves != n) throw std::invalid_argument("genomics::msa::refine: the tree's leaves are not the rows");
  if (n < 2) {
    if (report) *report = rep;
    return msa;
  }
  // The edges by decreasing distance from the root: increasing height, then node. The root's
  // second child splits as its first does, and is left out.
  std::vector<std::size_t> edges;
  const std::size_t skip = tree.merges.back().b;
  for (std::size_t v = 0; v + 1 < tree.nodes(); v++)
    if (v != skip) edges.push_back(v);
  std::stable_sort(edges.begin(), edges.end(), [&](std::size_t a, std::size_t b) { return tree.height(a) < tree.height(b); });
  std::vector<bool> inside(n);
  for (int round = 0; round < rounds; round++) {
    rep.rounds++;
    bool changed = false;
    for (std::size_t v : edges) {
      auto left = tree.leaves_under(v);
      std::sort(left.begin(), left.end());
      std::fill(inside.begin(), inside.end(), false);
      for (std::size_t r : left) inside[r] = true;
      std::vector<std::size_t> right;
      for (std::size_t r = 0; r < n; r++)
        if (!inside[r]) right.push_back(r);
      const Msa X = sub_alignment(msa, left), Y = sub_alignment(msa, right);
      const Profile px = make_profile(X, detail::subset(weights, left), s.matrix.size);
      const Profile py = make_profile(Y, detail::subset(weights, right), s.matrix.size);
      const ProfileAlignment pa = align_profiles(px, py, s, terminal_open);
      std::vector<std::size_t> ids = left;
      ids.insert(ids.end(), right.begin(), right.end());
      Msa candidate = detail::in_input_order(merge(X, Y, pa.path), ids);
      rep.trials++;
      const std::int64_t before = sp_across(msa, left, right, s, terminal_open, pool);
      const std::int64_t after = sp_across(candidate, left, right, s, terminal_open, pool);
      if (after > before) {
        msa = std::move(candidate);
        rep.accepted++;
        rep.gain += after - before;
        changed = true;
      }
    }
    if (!changed) break;
  }
  if (report) *report = rep;
  return msa;
}

// ---- the three stages ------------------------------------------------------------------------------------
struct ProgressiveOptions {
  align::Scoring scoring;          // the residue scores and gap costs (default_scoring)
  bool terminal_open = false;      // whether a gap at an end pays its opening
  bool weighted = true;            // CLUSTAL W's weights in the profiles
  int k = 0;                       // the k-mer length; 0: 4 (DNA), 4 on Dayhoff(6) (protein)
  bool second_stage = true;        // re-estimate the tree from the first alignment
  int refine_rounds = 2;           // passes of tree-dependent restricted partitioning
};

// The default scorings, chosen by the mean Q over simulated families (the module's README, §6:
// twelve DNA families of 8 and 32 members, twelve protein families): DNA +2/−1 with open 4,
// extend 1; protein BLOSUM62 with open 6, extend 1 (a gap of length L costing open + L · extend,
// align/types.hpp). Other scorings near these give nearly the same accuracy.
inline align::Scoring default_scoring(Alphabet a) {
  if (a == Alphabet::dna) return align::scoring(align::Matrix::dna(2, -1, -1), 4, 1);
  return align::scoring(align::matrix("BLOSUM62"), 6, 1);
}
inline ProgressiveOptions default_options(Alphabet a) {
  ProgressiveOptions o;
  o.scoring = default_scoring(a);
  return o;
}

struct ProgressiveResult {
  Msa msa;
  Dendrogram tree;              // the tree the kept alignment followed
  std::vector<double> weights;  // its CLUSTAL W weights
  std::int64_t sp_first = 0, sp_second = 0, sp = 0;  // after stage 1, stage 2 (the kept one's), and at the end
  bool second_kept = false;
  Refinement refinement;
};

inline ProgressiveResult progressive(const std::vector<Codes>& seqs, const std::vector<std::string>& names, Alphabet alphabet, const ProgressiveOptions& o,
                                     calx::scheduler::ThreadPool& pool) {
  detail::validate_input(seqs, names, o.scoring, "genomics::msa::progressive");
  ProgressiveResult res;
  const std::size_t n = seqs.size();
  if (n == 0) return res;
  const Compression comp = alphabet == Alphabet::dna ? dna_letters() : dayhoff6();
  const int k = o.k > 0 ? o.k : 4;
  const DistanceMatrix d1 = kmer_distances(seqs, k, comp, pool);
  res.tree = upgma(d1);
  res.weights = o.weighted ? clustal_weights(res.tree) : std::vector<double>(n, 1.0);
  res.msa = align_along(seqs, names, res.tree, res.weights, o.scoring, o.terminal_open, pool);
  res.sp_first = res.sp_second = sp_score(res.msa, o.scoring, o.terminal_open, pool);
  if (o.second_stage && n > 2) {
    const Dendrogram t2 = upgma(corrected_distances(res.msa, alphabet, pool));
    if (detail::clusters(t2) != detail::clusters(res.tree)) {
      const auto w2 = o.weighted ? clustal_weights(t2) : std::vector<double>(n, 1.0);
      Msa a2 = align_along(seqs, names, t2, w2, o.scoring, o.terminal_open, pool);
      const std::int64_t sp2 = sp_score(a2, o.scoring, o.terminal_open, pool);
      if (sp2 >= res.sp_first) {
        res.msa = std::move(a2);
        res.tree = t2;
        res.weights = w2;
        res.sp_second = sp2;
        res.second_kept = true;
      }
    }
  }
  res.msa = refine(std::move(res.msa), res.tree, res.weights, o.scoring, o.terminal_open, o.refine_rounds, pool, &res.refinement);
  res.sp = res.sp_second + res.refinement.gain;
  return res;
}

}  // namespace genomics::msa
