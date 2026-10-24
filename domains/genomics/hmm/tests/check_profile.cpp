// Checks of profile HMMs: the construction from Durbin et al.'s alignment of seven globins
// (1998, Figure 5.3) gives its eight match columns and the probabilities that Laplace's rule
// gives its counts, worked by hand from the alignment; on a small profile every state path of
// short sequences is enumerated, and the Viterbi and forward scores of the three modes are held
// to the maximum and the sum over them; an alignment's path scores what Viterbi claims and
// consumes its residues in order; and a profile built from a simulated family's true alignment
// scores held-out members far above the same sequences shuffled, and aligns their residues to
// the match states they truly belong to.

#include <genomics/hmm/profile.hpp>
#include <genomics/simulate/family.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::hmm;

namespace {

std::string fixed(double v, int digits = 1) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*f", digits, v);
  return b;
}

Codes protein_row(const std::string& s) {
  Codes out;
  for (char c : s) out.push_back(c == '-' ? sequence::GAP : sequence::protein_code(c));
  return out;
}
Codes dna_row(const std::string& s) {
  Codes out;
  for (char c : s) out.push_back(c == '-' ? sequence::GAP : sequence::dna_code(c));
  return out;
}

// Every path of the profile through x, by recursion on (state, node, residues consumed), with
// weights computed from the profile's probabilities directly: the definition the dynamic
// programme is held to.
struct Enumeration {
  double best = NEG_INF, total = NEG_INF;
  std::size_t paths = 0;
};

struct Walker {
  const Profile& p;
  const Codes& x;
  ProfileMode mode;
  Enumeration out;
  double lm(std::size_t k, Code b) const { return b < p.symbols ? std::log(p.e_match(k, b) / p.background[b]) : 0; }
  double li(std::size_t k, Code b) const { return b < p.symbols ? std::log(p.e_insert(k, b) / p.background[b]) : 0; }
  double lt(std::size_t k, int w) const { return p.transition[k][w] > 0 ? std::log(p.transition[k][w]) : NEG_INF; }
  void record(double w) {
    if (w == NEG_INF) return;
    out.best = std::max(out.best, w);
    out.total = log_add(out.total, w);
    out.paths++;
  }
  // kind: 0 M, 1 I, 2 D.
  void walk(int kind, std::size_t node, std::size_t i, double w) {
    if (w == NEG_INF) return;
    const std::size_t m = p.length, L = x.size();
    const bool local = mode == ProfileMode::local;
    if (!local && node == m && (mode != ProfileMode::global || i == L)) record(w + lt(m, 3 * kind + 0));
    if (local && kind == 0 && node >= 1) record(w);
    if (node + 1 <= m && i < L) walk(0, node + 1, i + 1, w + lt(node, 3 * kind + 0) + lm(node + 1, x[i]));
    if (i < L && !(local && (node == 0 || node == m))) walk(1, node, i + 1, w + lt(node, 3 * kind + 1) + li(node, x[i]));
    if (node + 1 <= m) walk(2, node + 1, i, w + lt(node, 3 * kind + 2));
  }
  Enumeration run() {
    const std::size_t m = p.length, L = x.size();
    if (mode == ProfileMode::global) walk(0, 0, 0, 0);
    else if (mode == ProfileMode::glocal)
      for (std::size_t i = 0; i <= L; i++) walk(0, 0, i, 0);
    else {
      const double entry = std::log(2.0 / (static_cast<double>(m) * static_cast<double>(m + 1)));
      for (std::size_t i = 0; i < L; i++)
        for (std::size_t k = 1; k <= m; k++) walk(0, k, i + 1, entry + lm(k, x[i]));
    }
    out.best /= std::log(2.0);
    out.total /= std::log(2.0);
    return out;
  }
};

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- Durbin et al.'s globin alignment (Figure 5.3) --------------------------------------------------
  {
    // Ten columns of seven globins, Durbin et al. 1998, Figure 5.3 (also the FIU CAP5510 lecture 8,
    // 2011, users.cs.fiu.edu/~giri/teach/Bioinf/S11/Lec8.pdf, retrieved 2026-10-24). The book stars
    // columns 1–3 and 6–10 as match columns.
    const std::vector<Codes> rows = {protein_row("VGA--HAGEY"), protein_row("V----NVDEV"), protein_row("VEA--DVAGH"), protein_row("VKG------D"),
                                     protein_row("VYS--TYETS"), protein_row("FNA--NIPKH"), protein_row("IAGADNGAGV")};
    const Profile p = build_profile(rows, 20);
    s.require("the 50% rule makes columns 1–3 and 6–10 the match columns: eight match states, as the book's stars", p.length == 8 && p.columns == std::vector<std::size_t>{0, 1, 2, 5, 6, 7, 8, 9});
    const Code V = sequence::protein_code('V'), F = sequence::protein_code('F'), I = sequence::protein_code('I'), N = sequence::protein_code('N'),
               H = sequence::protein_code('H'), W = sequence::protein_code('W');
    // Worked by hand from the alignment, with Laplace's rule (one pseudocount per symbol and per transition):
    // M1 sees V five times, F and I once: (5 + 1)/(7 + 20) = 6/27, 2/27, and 1/27 for the others.
    double d = std::fabs(p.e_match(1, V) - 6.0 / 27) + std::fabs(p.e_match(1, F) - 2.0 / 27) + std::fabs(p.e_match(1, I) - 2.0 / 27) + std::fabs(p.e_match(1, W) - 1.0 / 27);
    // M4 (column 6) sees N three times, H, D and T once, in six residues: 4/26 and 2/26.
    d += std::fabs(p.e_match(4, N) - 4.0 / 26) + std::fabs(p.e_match(4, H) - 2.0 / 26);
    s.check("match emissions with Laplace's rule: M1 has V 6/27, F and I 2/27, others 1/27; M4 has N 4/26, H 2/26", d, 1e-15);
    // Out of node 3: M3 → M4 four times (HBA, MYG, GLB5, LGB2), M3 → I3 once (GLB1), M3 → D4 once (GLB3):
    // (5, 2, 2)/9; I3 → I3 once and I3 → M4 once: (2, 2, 1)/5 for IM, II, ID; D3 → M4 once (HBB): (2, 1, 1)/4.
    double t = std::fabs(p.t(3, MM) - 5.0 / 9) + std::fabs(p.t(3, MI) - 2.0 / 9) + std::fabs(p.t(3, MD) - 2.0 / 9);
    t += std::fabs(p.t(3, IM) - 2.0 / 5) + std::fabs(p.t(3, II) - 2.0 / 5) + std::fabs(p.t(3, ID) - 1.0 / 5);
    t += std::fabs(p.t(3, DM) - 2.0 / 4) + std::fabs(p.t(3, DI) - 1.0 / 4) + std::fabs(p.t(3, DD) - 1.0 / 4);
    // Begin: seven rows to M1 of nine pseudocounted (8, 1, 1)/10; out of M1: six to M2, one (HBB) to D2: (7, 1, 2)/10;
    // node 8 into the end: (MM, MI) = (8, 1)/9, the delete transitions absent.
    t += std::fabs(p.t(0, MM) - 0.8) + std::fabs(p.t(0, MD) - 0.1) + std::fabs(p.t(1, MM) - 0.7) + std::fabs(p.t(1, MD) - 0.2);
    t += std::fabs(p.t(8, MM) - 8.0 / 9) + std::fabs(p.t(8, MI) - 1.0 / 9) + p.t(8, MD) + p.t(0, DM);
    s.check("transitions with Laplace's rule: node 3 (5, 2, 2)/9, (2, 2, 1)/5, (2, 1, 1)/4; begin (8, 1, 1)/10; M1 (7, 1, 2)/10; into the end (8, 1)/9", t, 1e-15);
    const Codes hba = protein_row("VGAHAGEY"), glb1 = protein_row("IAGADNGAGV");
    const ProfileAlignment a = align_to_profile(p, glb1, ProfileMode::global);
    std::size_t inserted = 0;
    for (auto n : a.inserted) inserted += n;
    s.require("GLB1_GLYDI aligns globally with its two residues of columns 4–5 as inserts after node 3", a.inserted[3] == 2 && inserted == 2 && a.match_row == protein_row("IAGNGAGV"));
    s.require("a member of the alignment scores higher than its reversal", profile_viterbi(p, hba, ProfileMode::global) > profile_viterbi(p, Codes(hba.rbegin(), hba.rend()), ProfileMode::global));
  }

  // ---- every path of a small profile ---------------------------------------------------------------------------
  {
    const std::vector<Codes> rows = {dna_row("AC-GT"), dna_row("AC-G-"), dna_row("T-TGA"), dna_row("AC--A")};
    ProfileOptions o;
    o.insert_background = false;  // inserts with emissions of their own, so that their log-odds is not zero
    o.background = {0.3, 0.2, 0.2, 0.3};
    const Profile p = build_profile(rows, 4, o);
    calx::random::Stream r(5, 0x70726f66, 0, 0);
    double worst_v = 0, worst_f = 0;
    std::size_t paths = 0, cases = 0;
    bool aligned = true;
    for (int rep = 0; rep < 24; rep++) {
      Codes x(static_cast<std::size_t>(rep % 7));
      for (auto& c : x) c = static_cast<Code>(r.below(4));
      if (rep % 5 == 4 && !x.empty()) x[0] = sequence::N;  // missing
      for (ProfileMode mode : {ProfileMode::global, ProfileMode::glocal, ProfileMode::local}) {
        if (mode == ProfileMode::local && x.empty()) continue;
        Walker w{p, x, mode, {}};
        const Enumeration e = w.run();
        cases++;
        paths += e.paths;
        const double v = profile_viterbi(p, x, mode), f = profile_forward(p, x, mode), fl = profile_forward_log(p, x, mode);
        worst_v = std::max(worst_v, std::fabs(v - e.best));
        worst_f = std::max({worst_f, std::fabs(f - e.total), std::fabs(fl - e.total)});
        const ProfileAlignment a = align_to_profile(p, x, mode);
        bool ok = std::fabs(a.score - v) < 1e-12 && std::fabs(profile_path_score(p, x, a, mode) - v) < 1e-9;
        std::size_t next = a.begin;
        for (const auto& st : a.steps)
          if (st.state != 'D') ok = ok && st.position == next++;
        ok = ok && next == a.end && (mode != ProfileMode::global || (a.begin == 0 && a.end == x.size()));
        aligned = aligned && ok;
      }
    }
    s.check("Viterbi scores of global, glocal and local modes equal the maximum over every path (" + std::to_string(cases) + " cases, " + std::to_string(paths) + " paths)", worst_v, 1e-9, "bits");
    s.check("forward scores, scaled and in the log semiring, equal the log of the sum over every path", worst_f, 1e-9, "bits");
    s.require("each alignment's path rescores to its Viterbi score, consuming its residues in order (and all of them in global mode)", aligned);
  }

  // ---- a simulated family ---------------------------------------------------------------------------------------
  {
    const Codes root = simulate::random_sequence(150, 0.5, 70);
    const simulate::Family fam = simulate::evolve_family(root, 32, 0.06, 0.01, 71);
    // Every fourth member held out, so that the held-out ones are spread over the tree.
    std::vector<Codes> training;
    std::vector<std::size_t> held_out;
    for (std::size_t j = 0; j < fam.alignment.size(); j++) {
      if (j % 4 == 3) held_out.push_back(j);
      else training.push_back(fam.alignment[j]);
    }
    ProfileOptions o;
    o.background = {0.25, 0.25, 0.25, 0.25};
    const Profile p = build_profile(training, 4, o);
    calx::random::Stream r(72, 0x73687566, 0, 0);
    double min_member = INFINITY, max_shuffled = -INFINITY, min_local = INFINITY, max_local_shuffled = -INFINITY;
    double right = 0, total = 0;
    for (std::size_t j : held_out) {
      const Codes& x = fam.sequences[j];
      Codes y = x;
      for (std::size_t i = y.size(); i > 1; i--) std::swap(y[i - 1], y[r.below(i)]);  // Fisher–Yates
      min_member = std::min(min_member, profile_viterbi(p, x));
      max_shuffled = std::max(max_shuffled, profile_viterbi(p, y));
      min_local = std::min(min_local, profile_forward(p, x, ProfileMode::local));
      max_local_shuffled = std::max(max_local_shuffled, profile_forward(p, y, ProfileMode::local));
      // The true place of each residue in a match column, against the alignment's.
      std::vector<std::size_t> index(fam.alignment[j].size(), SIZE_MAX);
      for (std::size_t c = 0, n = 0; c < index.size(); c++)
        if (fam.alignment[j][c] != sequence::GAP) index[c] = n++;
      const ProfileAlignment a = align_to_profile(p, x);
      std::vector<std::size_t> placed(p.length + 1, SIZE_MAX);
      for (const auto& st : a.steps)
        if (st.state == 'M') placed[st.node] = st.position;
      for (std::size_t k = 1; k <= p.length; k++) {
        const std::size_t truth = index[p.columns[k - 1]];
        if (truth == SIZE_MAX) continue;
        total++;
        right += placed[k] == truth;
      }
    }
    s.require("glocal Viterbi: every held-out member of a simulated family scores above every shuffled one, by more than 40 bits",
              min_member - max_shuffled > 40, "members from " + fixed(min_member) + " bits, shuffled up to " + fixed(max_shuffled) + " bits");
    s.require("local forward: members above shuffled sequences by more than 40 bits", min_local - max_local_shuffled > 40,
              "members from " + fixed(min_local) + " bits, shuffled up to " + fixed(max_local_shuffled) + " bits");
    // Where an indel falls in a run of equal bases, or inserted bases happen to match, the true
    // alignment is one of several equally good ones: agreement is high, not complete.
    // The scaled forward against the log semiring on long sequences, where scaling matters.
    double worst = 0;
    const Codes far = simulate::random_sequence(20000, 0.5, 73);
    Codes with_member = far;
    with_member.insert(with_member.begin() + 7000, fam.sequences[held_out[0]].begin(), fam.sequences[held_out[0]].end());
    for (const Codes* y : std::vector<const Codes*>{&far, &with_member, &fam.sequences[held_out[1]]})
      for (ProfileMode mode : {ProfileMode::global, ProfileMode::glocal, ProfileMode::local}) {
        const double a = profile_forward(p, *y, mode), b = profile_forward_log(p, *y, mode);
        worst = std::max(worst, std::fabs(a - b) / std::max(1.0, std::fabs(b)));
      }
    s.check("the scaled forward equals the log-semiring forward on 20 kb, with and without a member inside, in every mode (relative)", worst, 1e-10);
    s.require("the alignment to the profile places at least 90% of held-out members' residues at their true match states", right / total >= 0.90, fixed(100 * right / total, 1) + "% of " + std::to_string(static_cast<int>(total)));
  }

  // ---- arguments ----------------------------------------------------------------------------------------------------
  {
    bool threw = false;
    try { build_profile({dna_row("AC"), dna_row("A")}, 4); } catch (const std::invalid_argument&) { threw = true; }
    s.require("build_profile rejects rows of different widths", threw);
  }
  return s.done();
}
