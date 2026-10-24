// Profile hidden Markov models: a family of sequences as a chain of match, insert and delete
// states built from its multiple alignment, and the scores and alignments of a sequence
// against it (Krogh, Brown, Mian, Sjölander and Haussler 1994, J. Mol. Biol. 235, 1501–1531;
// Durbin, Eddy, Krogh and Mitchison 1998, Biological Sequence Analysis, ch. 5; Eddy 1998,
// Bioinformatics 14, 755–763).
//
// - The architecture (Durbin et al. Figure 5.2). Nodes k = 0 … m: match states M_1 … M_m
//   emitting with e_{M_k}, insert states I_0 … I_m emitting with e_{I_k}, silent delete states
//   D_1 … D_m; M_0 is the begin state and M_{m+1} the end. Out of node k, nine transitions:
//   M_k, I_k and D_k each to M_{k+1}, I_k and D_{k+1} (written MM, MI, MD, IM, II, ID, DM, DI, DD),
//   where D_0 and D_{m+1} do not exist.
// - Construction from an alignment (Durbin et al. §5.3). Rows of codes with sequence::GAP; a
//   column is a match column when at most half its rows (by weight) are gaps, the 50% rule;
//   each row's state path is then read off (a residue in a match column is M_k, a gap there
//   D_k, a residue in an insert column I_k, a gap there nothing), and the transitions and
//   emissions along the paths are counted. Pseudocounts: Laplace's rule (one per symbol and
//   per transition, Durbin et al. §5.3), or A·q_b proportional to the background (§5.6, the
//   simplest Dirichlet prior); Dirichlet mixtures (Sjölander et al. 1996) are not yet here.
//   Inserts emit the background by default, so that their log-odds is zero (§5.3).
// - Scores are log-odds against the null model of i.i.d. residues with frequencies q_b
//   (Durbin et al. §5.4), in bits: the Viterbi score log₂ max_π P(x, π | M)/P(x | R) and the
//   forward score log₂ Σ_π P(x, π | M)/P(x | R), by the recurrences of §5.4 for V^M_j(i),
//   V^I_j(i), V^D_j(i) (and their sums), run row by row over the residues i = 0 … L and, within
//   a row, node by node (a delete state depends on the node before it in the same row). The
//   forward sum runs in probability space with each row rescaled (Rabiner 1989, §V.A), checked
//   against the same sum in the log semiring.
// - Modes, by which ends are free:
//     global   the whole sequence against the whole model (§5.4);
//     glocal   the whole model against a part of the sequence: the residues before and after
//              are flanks that emit the background and cost nothing (log-odds 0), as the free
//              end gaps of semi-global alignment; HMMER2's "ls" mode;
//     local    a fragment of the model against a part of the sequence: entry into any M_k with
//              probability 2/(m(m + 1)) and exit from any M_k with probability 1, the uniform
//              fragment distribution of HMMER3 (Eddy 2008, PLoS Comput. Biol. 4, e1000069),
//              without I_0 and I_m (their place is the flanks').
//   *Caveat.* The flanks carry no length distribution, so a forward score in the free modes
//   sums over every placement with weight one, and exceeds HMMER's (whose N, C and null-model
//   loops make it a probability ratio) by up to about log₂ L bits.
// - The alignment of a sequence to the profile: the Viterbi path, with the residue of each match
//   state (an A2M row of the family's match columns) and the residues inserted after each node.
//
// Complexity: O(L m) time for a sequence of length L against m match states (nine semiring
// operations a cell: on one core of an i7-11370H, about 140 million cells a second for Viterbi,
// 80 million for the scaled forward, 3 million for the log-semiring forward, measured on a
// profile of 302 match states against 100 kb); O(m) memory for scores, O(L m)
// for an alignment (three values a cell, from which the traceback recomputes its choices).
// Building: O(n W) for n rows of width W.
//
// Limits: the scores are exact for the model; the model is only as good as its alignment and
// its prior, and sequences are not weighted unless weights are given (Henikoff and Henikoff
// 1994's position-based weights, or Krogh and Mitchison's maximum entropy weights, are left to
// the caller).

#pragma once

#include <genomics/hmm/model.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace genomics::hmm {

enum ProfileTransition : std::uint8_t { MM, MI, MD, IM, II, ID, DM, DI, DD };

struct Profile {
  std::size_t length = 0;   // m
  std::size_t symbols = 0;  // K
  std::vector<double> match;       // e_{M_k}(b), (m + 1) × K, row 0 unused
  std::vector<double> insert;      // e_{I_k}(b), (m + 1) × K
  std::vector<double> background;  // q_b, K
  std::vector<std::array<double, 9>> transition;  // out of node k, k = 0 … m
  std::vector<std::size_t> columns;  // the alignment columns of M_1 … M_m

  double e_match(std::size_t k, std::size_t b) const { return match[k * symbols + b]; }
  double e_insert(std::size_t k, std::size_t b) const { return insert[k * symbols + b]; }
  double t(std::size_t k, ProfileTransition which) const { return transition[k][which]; }
};

struct ProfileOptions {
  double gap_fraction = 0.5;          // a column with more gaps than this fraction of the rows is an insert column
  double transition_pseudocount = 1;  // added to each transition that exists (Laplace)
  double emission_pseudocount = 1;    // Laplace: added to each symbol; with background_prior, A = this × K
  bool background_prior = false;      // pseudocounts A q_b instead of Laplace's
  bool insert_background = true;      // inserts emit q
  std::vector<double> background;     // q; empty: the alignment's residue frequencies (with one pseudocount each)
  std::vector<double> weights;        // one per row; empty: all one
};

enum class ProfileMode { global, glocal, local };

namespace detail {

// Whether a transition out of node k exists in a profile of length m.
constexpr bool transition_exists(std::size_t k, std::size_t m, int which) {
  if (k == 0 && (which == DM || which == DI || which == DD)) return false;  // no D_0
  if (k == m && (which == MD || which == ID || which == DD)) return false;  // no D_{m+1}
  return true;
}

}  // namespace detail

// The profile of an alignment, its rows of codes (sequence::GAP for gaps; a code ≥ K other than
// GAP occupies its column but adds no emission count).
inline Profile build_profile(const std::vector<Codes>& rows, std::size_t symbols, const ProfileOptions& o = {}) {
  if (rows.empty() || symbols == 0) throw std::invalid_argument("genomics::hmm::build_profile: an alignment of at least one row, and symbols");
  const std::size_t n = rows.size(), W = rows[0].size(), K = symbols;
  for (const auto& r : rows)
    if (r.size() != W) throw std::invalid_argument("genomics::hmm::build_profile: the rows of an alignment have one width");
  if (!o.weights.empty() && o.weights.size() != n) throw std::invalid_argument("genomics::hmm::build_profile: one weight per row");
  const auto weight = [&](std::size_t r) { return o.weights.empty() ? 1.0 : o.weights[r]; };
  double total = 0;
  for (std::size_t r = 0; r < n; r++) total += weight(r);

  Profile p;
  p.symbols = K;
  std::vector<bool> is_match(W, false);
  for (std::size_t c = 0; c < W; c++) {
    double gaps = 0;
    for (std::size_t r = 0; r < n; r++) gaps += rows[r][c] == sequence::GAP ? weight(r) : 0;
    is_match[c] = gaps <= o.gap_fraction * total;
    if (is_match[c]) p.columns.push_back(c);
  }
  const std::size_t m = p.length = p.columns.size();
  if (m == 0) throw std::invalid_argument("genomics::hmm::build_profile: the alignment has no match column");

  if (!o.background.empty()) {
    if (o.background.size() != K) throw std::invalid_argument("genomics::hmm::build_profile: one background frequency per symbol");
    p.background = o.background;
  } else {
    p.background.assign(K, 1.0);
    for (std::size_t r = 0; r < n; r++)
      for (Code c : rows[r])
        if (c < K) p.background[c] += weight(r);
  }
  {
    double s = 0;
    for (double q : p.background) s += q;
    for (double& q : p.background) q /= s;
  }

  std::vector<double> mc((m + 1) * K, 0), ic((m + 1) * K, 0);
  std::vector<std::array<double, 9>> tc(m + 1);
  for (auto& a : tc) a.fill(0);
  enum Kind { M_, I_, D_ };
  for (std::size_t r = 0; r < n; r++) {
    const double w = weight(r);
    Kind kind = M_;
    std::size_t node = 0, k = 0;
    const auto move = [&](Kind to) {
      // From (kind, node) to (to, node') with node' = node + 1 for M and D, node for I.
      tc[node][3 * kind + to] += w;
      node = to == I_ ? node : node + 1;
      kind = to;
    };
    for (std::size_t c = 0; c < W; c++) {
      const Code x = rows[r][c];
      const bool residue = x != sequence::GAP;
      if (is_match[c]) {
        k++;
        move(residue ? M_ : D_);
        if (residue && x < K) mc[k * K + x] += w;
      } else if (residue) {
        move(I_);
        if (x < K) ic[k * K + x] += w;
      }
    }
    tc[node][3 * kind + M_] += w;  // into M_{m+1}, the end
  }
  p.transition.resize(m + 1);
  for (std::size_t k = 0; k <= m; k++)
    for (int from = 0; from < 3; from++) {
      double s = 0;
      int allowed = 0;
      for (int to = 0; to < 3; to++)
        if (detail::transition_exists(k, m, 3 * from + to)) { s += tc[k][3 * from + to] + o.transition_pseudocount; allowed++; }
      for (int to = 0; to < 3; to++) {
        const int which = 3 * from + to;
        if (!detail::transition_exists(k, m, which)) p.transition[k][which] = 0;
        else p.transition[k][which] = s > 0 ? (tc[k][which] + o.transition_pseudocount) / s : 1.0 / allowed;
      }
    }
  const auto emissions = [&](std::vector<double>& out, const std::vector<double>& counts, std::size_t first) {
    out.assign((m + 1) * K, 0);
    for (std::size_t k = first; k <= m; k++) {
      double s = 0;
      for (std::size_t b = 0; b < K; b++) s += counts[k * K + b];
      const double A = o.emission_pseudocount * static_cast<double>(K);  // the total pseudocount, either way
      for (std::size_t b = 0; b < K; b++) {
        const double prior = o.background_prior ? A * p.background[b] : o.emission_pseudocount;
        out[k * K + b] = s + A > 0 ? (counts[k * K + b] + prior) / (s + A) : p.background[b];
      }
    }
  };
  emissions(p.match, mc, 1);
  if (o.insert_background) {
    p.insert.assign((m + 1) * K, 0);
    for (std::size_t k = 0; k <= m; k++)
      for (std::size_t b = 0; b < K; b++) p.insert[k * K + b] = p.background[b];
  } else {
    emissions(p.insert, ic, 0);
  }
  return p;
}

namespace detail {

struct ProfileTables {
  std::size_t m = 0, K = 0;
  std::vector<double> lm, li;  // log-odds (nats), (m + 1) × (K + 1), column K the missing symbol (0)
  std::vector<std::array<double, 9>> lt;
  double entry = 0;  // log 2/(m(m + 1))
};

inline ProfileTables tables(const Profile& p) {
  ProfileTables t;
  t.m = p.length;
  t.K = p.symbols;
  const std::size_t K1 = t.K + 1;
  t.lm.assign((t.m + 1) * K1, NEG_INF);
  t.li.assign((t.m + 1) * K1, NEG_INF);
  for (std::size_t k = 0; k <= t.m; k++) {
    for (std::size_t b = 0; b < t.K; b++) {
      if (k >= 1) t.lm[k * K1 + b] = safe_log(p.e_match(k, b)) - std::log(p.background[b]);
      t.li[k * K1 + b] = safe_log(p.e_insert(k, b)) - std::log(p.background[b]);
    }
    t.lm[k * K1 + t.K] = k >= 1 ? 0 : NEG_INF;
    t.li[k * K1 + t.K] = 0;
  }
  t.lt.resize(t.m + 1);
  for (std::size_t k = 0; k <= t.m; k++)
    for (int w = 0; w < 9; w++) t.lt[k][w] = safe_log(p.transition[k][w]);
  t.entry = std::log(2.0 / (static_cast<double>(t.m) * static_cast<double>(t.m + 1)));
  return t;
}

// The sum over the paths of the profile in R, in nats. With `store`, the values of M, I and D
// of every row are kept (3 (L + 1)(m + 1) of them), for a traceback; otherwise two rows.
template <class R>
double profile_sum(const ProfileTables& tb, const Codes& x, ProfileMode mode, std::vector<double>* store) {
  const std::size_t m = tb.m, L = x.size(), W = m + 1, K1 = tb.K + 1;
  const bool local = mode == ProfileMode::local, free_ends = mode != ProfileMode::global;
  const std::size_t rows = store ? L + 1 : 2;
  std::vector<double> scratch;
  std::vector<double>& buf = store ? *store : scratch;
  buf.assign(3 * rows * W, R::zero());
  const auto row = [&](std::size_t i, int which) { return &buf[((store ? i : (i & 1)) * 3 + static_cast<std::size_t>(which)) * W]; };
  double total = R::zero();
  for (std::size_t i = 0; i <= L; i++) {
    double* M = row(i, 0);
    double* I = row(i, 1);
    double* D = row(i, 2);
    if (!store && i >= 2) {
      std::fill(M, M + W, R::zero());
      std::fill(I, I + W, R::zero());
      std::fill(D, D + W, R::zero());
    }
    const double* Mp = i ? row(i - 1, 0) : nullptr;
    const double* Ip = i ? row(i - 1, 1) : nullptr;
    const double* Dp = i ? row(i - 1, 2) : nullptr;
    const double begin = (i == 0 || free_ends) ? R::one() : R::zero();
    M[0] = local ? R::zero() : begin;
    const std::size_t sym = i ? std::min<std::size_t>(x[i - 1], tb.K) : 0;
    for (std::size_t k = 0; k <= m; k++) {
      if (i >= 1 && k >= 1) {
        const auto& t = tb.lt[k - 1];
        double s = R::plus(R::plus(R::times(Mp[k - 1], t[MM]), R::times(Ip[k - 1], t[IM])), R::times(Dp[k - 1], t[DM]));
        if (local) s = R::plus(s, tb.entry);  // from the begin state, after a free flank
        M[k] = R::times(s, tb.lm[k * K1 + sym]);
      }
      if (i >= 1 && !(local && (k == 0 || k == m))) {
        const auto& t = tb.lt[k];
        const double s = R::plus(R::plus(R::times(Mp[k], t[MI]), R::times(Ip[k], t[II])), R::times(Dp[k], t[DI]));
        I[k] = R::times(s, tb.li[k * K1 + sym]);
      }
      if (k >= 1) {
        const auto& t = tb.lt[k - 1];
        D[k] = R::plus(R::plus(R::times(M[k - 1], t[MD]), R::times(I[k - 1], t[ID])), R::times(D[k - 1], t[DD]));
      }
    }
    if (!free_ends && i != L) continue;
    double e = R::zero();
    if (local)
      for (std::size_t k = 1; k <= m; k++) e = R::plus(e, M[k]);
    else {
      const auto& t = tb.lt[m];
      e = R::plus(R::plus(R::times(M[m], t[MM]), R::times(I[m], t[IM])), R::times(D[m], t[DM]));
    }
    total = R::plus(total, e);
  }
  return total;
}

// The forward sum in probability space, each row of the dynamic programme rescaled to sum to
// one (Rabiner 1989, §V.A, applied to the rows i of the profile's lattice), in nats. A row
// computed from the row before is in that row's units C_{i−1}; the begin state, worth 1 at every
// row in the free modes, enters as 1/C_{i−1}, and since it belongs to the row's sum, C_i ≥ 1 and
// that term never overflows; the end's contributions are added in log space. It equals
// profile_sum<LogSum> to rounding (the checks), at the cost of multiplications instead of
// exponentials and logarithms.
inline double profile_forward_scaled(const ProfileTables& tb, const Codes& x, ProfileMode mode) {
  const std::size_t m = tb.m, L = x.size(), W = m + 1, K1 = tb.K + 1;
  const bool local = mode == ProfileMode::local, free_ends = mode != ProfileMode::global;
  std::vector<double> om(tb.lm.size()), oi(tb.li.size());
  for (std::size_t i = 0; i < om.size(); i++) om[i] = std::exp(tb.lm[i]);
  for (std::size_t i = 0; i < oi.size(); i++) oi[i] = std::exp(tb.li[i]);
  std::vector<std::array<double, 9>> t(m + 1);
  for (std::size_t k = 0; k <= m; k++)
    for (int w = 0; w < 9; w++) t[k][w] = std::exp(tb.lt[k][w]);
  const double entry = std::exp(tb.entry);
  std::vector<double> Mp(W, 0), Ip(W, 0), Dp(W, 0), M(W), I(W), D(W);
  double logC = 0, total = NEG_INF;
  for (std::size_t i = 0; i <= L; i++) {
    const double unit = std::exp(-logC);  // the begin state, 1, in the units of the row before
    const double begin = (i == 0 || free_ends) ? unit : 0;
    const std::size_t sym = i ? std::min<std::size_t>(x[i - 1], tb.K) : 0;
    M[0] = local ? 0 : begin;
    I[0] = 0;
    D[0] = 0;
    for (std::size_t k = 0; k <= m; k++) {
      if (k >= 1) {
        const auto& a = t[k - 1];
        M[k] = i >= 1 ? om[k * K1 + sym] * (Mp[k - 1] * a[MM] + Ip[k - 1] * a[IM] + Dp[k - 1] * a[DM] + (local ? unit * entry : 0.0)) : 0;
        D[k] = M[k - 1] * a[MD] + I[k - 1] * a[ID] + D[k - 1] * a[DD];
      }
      if (i >= 1 && !(local && (k == 0 || k == m))) {
        const auto& a = t[k];
        I[k] = oi[k * K1 + sym] * (Mp[k] * a[MI] + Ip[k] * a[II] + Dp[k] * a[DI]);
      } else {
        I[k] = 0;
      }
    }
    if (free_ends || i == L) {
      double e = 0;
      if (local)
        for (std::size_t k = 1; k <= m; k++) e += M[k];
      else
        e = M[m] * t[m][MM] + I[m] * t[m][IM] + D[m] * t[m][DM];
      if (e > 0) total = log_add(total, std::log(e) + logC);
    }
    // The row's sum, the begin state's value included (M_0 holds it except in local mode, where
    // the begin enters every match state directly).
    double s = local ? unit : 0;
    for (std::size_t k = 0; k <= m; k++) s += M[k] + I[k] + D[k];
    if (!(s > 0)) break;  // no path reaches this row (global mode): nothing further can end
    const double inv = 1 / s;
    for (std::size_t k = 0; k <= m; k++) {
      M[k] *= inv;
      I[k] *= inv;
      D[k] *= inv;
    }
    logC += std::log(s);
    std::swap(M, Mp);
    std::swap(I, Ip);
    std::swap(D, Dp);
  }
  return total;
}

}  // namespace detail

// The Viterbi log-odds score of x against the profile, in bits (−∞ if no path emits x).
inline double profile_viterbi(const Profile& p, const Codes& x, ProfileMode mode = ProfileMode::glocal) {
  return detail::profile_sum<Tropical>(detail::tables(p), x, mode, nullptr) / std::numbers::ln2;
}

// The forward log-odds score, in bits, by the scaled recursion.
inline double profile_forward(const Profile& p, const Codes& x, ProfileMode mode = ProfileMode::glocal) {
  return detail::profile_forward_scaled(detail::tables(p), x, mode) / std::numbers::ln2;
}
// The same in the log semiring: the reference the scaled recursion is checked against.
inline double profile_forward_log(const Profile& p, const Codes& x, ProfileMode mode = ProfileMode::glocal) {
  return detail::profile_sum<LogSum>(detail::tables(p), x, mode, nullptr) / std::numbers::ln2;
}

struct ProfileStep {
  char state = 'M';                  // 'M', 'I' or 'D'
  std::size_t node = 0;              // k
  std::size_t position = SIZE_MAX;   // the residue emitted (M, I), SIZE_MAX for D
};

struct ProfileAlignment {
  double score = NEG_INF;           // the Viterbi log-odds, bits
  std::size_t begin = 0, end = 0;   // the residues aligned, [begin, end); the rest are flanks
  std::vector<ProfileStep> steps;   // in order
  Codes match_row;                  // m codes: the residue of M_k, GAP where deleted or outside a fragment
  std::vector<std::size_t> inserted;  // residues inserted after node k, k = 0 … m
};

// The Viterbi alignment of x to the profile.
inline ProfileAlignment align_to_profile(const Profile& p, const Codes& x, ProfileMode mode = ProfileMode::glocal) {
  const detail::ProfileTables tb = detail::tables(p);
  std::vector<double> buf;
  const double best = detail::profile_sum<Tropical>(tb, x, mode, &buf);
  ProfileAlignment out;
  out.score = best / std::numbers::ln2;
  out.match_row.assign(p.length, sequence::GAP);
  out.inserted.assign(p.length + 1, 0);
  if (best == NEG_INF) return out;
  const std::size_t m = tb.m, L = x.size(), W = m + 1;
  const bool local = mode == ProfileMode::local, free_ends = mode != ProfileMode::global;
  const auto V =[&](std::size_t i, int which, std::size_t k) { return buf[(i * 3 + static_cast<std::size_t>(which)) * W + k]; };
  // The first argument of the maximum, as std::max chooses: the first of equal values.
  const auto first_max = [](const double* c, int n) {
    int arg = 0;
    double v = c[0];
    for (int j = 1; j < n; j++)
      if (v < c[j]) { v = c[j]; arg = j; }
    return arg;
  };
  // The end: the row and the state it was reached from.
  std::size_t i = L;
  int state = 0;  // 0 M, 1 I, 2 D
  std::size_t k = m;
  {
    double bestv = NEG_INF;
    for (std::size_t r = free_ends ? 0 : L; r <= L; r++) {
      double c[3];
      std::size_t kk = m;
      int st = 0;
      double v;
      if (local) {
        v = NEG_INF;
        for (std::size_t j = 1; j <= m; j++)
          if (v < V(r, 0, j)) { v = V(r, 0, j); kk = j; }
        st = 0;
      } else {
        const auto& t = tb.lt[m];
        c[0] = V(r, 0, m) + t[MM];
        c[1] = V(r, 1, m) + t[IM];
        c[2] = V(r, 2, m) + t[DM];
        st = first_max(c, 3);
        v = c[st];
      }
      if (bestv < v) { bestv = v; i = r; k = kk; state = st; }
    }
  }
  out.end = i;
  // Back from the end, recomputing each choice as the recurrence made it. M_0, D_0 (and in
  // local mode I_0) hold the semiring's zero where they are not the begin state, so they are
  // never chosen wrongly; reaching M_0, or an entry in local mode, is the begin.
  std::vector<ProfileStep> steps;
  for (;;) {
    if (state == 0) {
      if (k == 0) break;
      steps.push_back({'M', k, i - 1});
      out.match_row[k - 1] = x[i - 1];
      const auto& t = tb.lt[k - 1];
      const double c[4] = {V(i - 1, 0, k - 1) + t[MM], V(i - 1, 1, k - 1) + t[IM], V(i - 1, 2, k - 1) + t[DM], local ? tb.entry : NEG_INF};
      const int j = first_max(c, 4);
      i--;
      if (j == 3) break;
      state = j;
      k--;
      if (state == 0 && k == 0) break;
    } else if (state == 1) {
      steps.push_back({'I', k, i - 1});
      out.inserted[k]++;
      const auto& t = tb.lt[k];
      double c[3] = {V(i - 1, 0, k) + t[MI], V(i - 1, 1, k) + t[II], V(i - 1, 2, k) + t[DI]};
      state = first_max(c, 3);
      i--;
      if (state == 0 && k == 0) break;
    } else {
      steps.push_back({'D', k, SIZE_MAX});
      const auto& t = tb.lt[k - 1];
      double c[3] = {V(i, 0, k - 1) + t[MD], V(i, 1, k - 1) + t[ID], V(i, 2, k - 1) + t[DD]};
      state = first_max(c, 3);
      k--;
      if (state == 0 && k == 0) break;
    }
  }
  out.begin = i;
  out.steps.assign(steps.rbegin(), steps.rend());
  return out;
}

// The log-odds score of an alignment's path, in bits, recomputed from the profile: what
// align_to_profile claims, held to the path it returns.
inline double profile_path_score(const Profile& p, const Codes& x, const ProfileAlignment& a, ProfileMode mode = ProfileMode::glocal) {
  const detail::ProfileTables tb = detail::tables(p);
  const std::size_t K1 = tb.K + 1, m = tb.m;
  const bool local = mode == ProfileMode::local;
  if (a.steps.empty()) return NEG_INF;
  double s = local ? tb.entry : 0;
  int kind = 0;  // the state before: M_0 (begin)
  std::size_t node = 0;
  for (std::size_t j = 0; j < a.steps.size(); j++) {
    const ProfileStep& st = a.steps[j];
    const int to = st.state == 'M' ? 0 : st.state == 'I' ? 1 : 2;
    if (!(local && j == 0)) s += tb.lt[node][3 * kind + to];
    const std::size_t sym = st.position == SIZE_MAX ? 0 : std::min<std::size_t>(x[st.position], tb.K);
    if (to == 0) s += tb.lm[st.node * K1 + sym];
    if (to == 1) s += tb.li[st.node * K1 + sym];
    kind = to;
    node = st.node;
  }
  if (!local) s += tb.lt[m][3 * kind + 0];
  return s / std::numbers::ln2;
}

}  // namespace genomics::hmm
