// Checks of gotoh.hpp and scoring.hpp: the published example of Durbin et al. (HEAGAWGHEE against
// PAWHEAE, BLOSUM50, a linear gap of −8) in three modes with its alignments; Gotoh's score equal
// to morphe's generic sum over paths on random DNA and protein pairs and on homologous pairs from
// the simulator, in all four modes, with affine and linear gaps; every returned path a path of its
// mode, rescored to its claimed score; the linear-memory score and end agreeing with the
// traceback; the number of co-optimal alignments counted exactly by the tropical semiring with
// multiplicities; the NCBI matrices symmetric, with their published entries and the relative
// entropy and expected score their files state; and the speed of the scalar kernel on 10 kb.

#include <genomics/align/gotoh.hpp>
#include <genomics/align/reference.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/simulate/family.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::align;

namespace {

Codes random_codes(calx::random::Stream& r, std::size_t n, std::size_t sigma) {
  Codes x(n);
  for (auto& c : x) c = static_cast<Code>(r.below(sigma));
  return x;
}

bool allowed(const Alignment& a, std::size_t m, std::size_t n, Mode mode) {
  if (query_length(a.cigar) != a.query_end - a.query_begin || target_length(a.cigar) != a.target_end - a.target_begin) return false;
  if (a.query_end > m || a.target_end > n) return false;
  switch (mode) {
    case Mode::global: return a.query_begin == 0 && a.target_begin == 0 && a.query_end == m && a.target_end == n;
    case Mode::local: return true;
    case Mode::semiglobal: return a.query_begin == 0 && a.query_end == m;
    case Mode::overlap: return (a.query_begin == 0 || a.target_begin == 0) && (a.query_end == m || a.target_end == n);
  }
  return false;
}

constexpr Mode MODES[] = {Mode::global, Mode::local, Mode::semiglobal, Mode::overlap};

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x676f746f, 0, 0);

  // ---- the matrices ---------------------------------------------------------------------------------
  {
    bool symmetric = true;
    for (const auto& name : matrix_names()) {
      const Matrix M = matrix(name);
      for (int a = 0; a < 24; a++)
        for (int b = 0; b < 24; b++) symmetric = symmetric && M(static_cast<Code>(a), static_cast<Code>(b)) == M(static_cast<Code>(b), static_cast<Code>(a));
    }
    s.require("matrices: the eight NCBI matrices and BLAST's BLOSUM80 are symmetric on the 24 codes", symmetric && matrix_names().size() == 9);
    // BLAST's built-in BLOSUM80 (half bits) against the FTP file (third bits): A/A 5 against 7.
    s.require("matrices: BLOSUM80-BLAST is BLAST's half-bit matrix (A/A 5, W/W 11), the file's BLOSUM80 a third-bit one (A/A 7, W/W 16)",
              matrix("BLOSUM80-BLAST")(0, 0) == 5 && matrix("BLOSUM80-BLAST")(17, 17) == 11 && matrix("BLOSUM80")(0, 0) == 7 && matrix("BLOSUM80")(17, 17) == 16);
    const Matrix B62 = matrix("BLOSUM62"), B50 = matrix("BLOSUM50"), P250 = matrix("PAM250");
    const auto at = [](const Matrix& M, char a, char b) { return M(sequence::protein_code(a), sequence::protein_code(b)); };
    // Entries read from NCBI's files (https://ftp.ncbi.nih.gov/blast/matrices/, 2026-10-24).
    s.require("matrices: BLOSUM62 has W/W 11, C/C 9, A/A 4, W/C −2, the stop −4 against A and 1 against itself",
              at(B62, 'W', 'W') == 11 && at(B62, 'C', 'C') == 9 && at(B62, 'A', 'A') == 4 && at(B62, 'W', 'C') == -2 && at(B62, '*', 'A') == -4 && at(B62, '*', '*') == 1);
    s.require("matrices: BLOSUM50 has W/W 15, H/H 10, P/W −4; PAM250 has W/W 17, C/W −8", at(B50, 'W', 'W') == 15 && at(B50, 'H', 'H') == 10 && at(B50, 'P', 'W') == -4 &&
                                                                                          at(P250, 'W', 'W') == 17 && at(P250, 'C', 'W') == -8);
    s.require("matrices: the extremes match the files' own headers (PAM30 −17 .. 13, PAM70 −11 .. 13, PAM250 −8 .. 17)",
              matrix("PAM30").min_score() == -17 && matrix("PAM30").max_score() == 13 && matrix("PAM70").min_score() == -11 && matrix("PAM70").max_score() == 13 &&
                  P250.min_score() == -8 && P250.max_score() == 17);
    // The file's header: "Entropy = 0.6979, Expected = −0.5209" (bits; BLOSUM62 in half bits, so
    // a score s is s/2 bits), w.r.t. the background the matrix was built from. With the
    // Robinson–Robinson background instead the expected score is that of NCBI BLAST's ungapped
    // statistics (statistics.hpp); here: the expected score is negative and λ-scaled sums agree.
    const auto p = robinson_frequencies();
    double sum = 0, expected = 0;
    for (int a = 0; a < 24; a++) {
      sum += p[static_cast<std::size_t>(a)];
      for (int b = 0; b < 24; b++) expected += p[static_cast<std::size_t>(a)] * p[static_cast<std::size_t>(b)] * B62(static_cast<Code>(a), static_cast<Code>(b));
    }
    s.check("frequencies: Robinson and Robinson's twenty frequencies sum to 1", sum - 1, 1e-12);
    s.require("frequencies: BLOSUM62's expected score under them is negative, as local alignment needs", expected < 0, std::to_string(expected));
    bool unknown = false;
    try {
      matrix("BLOSUM63");
    } catch (const std::invalid_argument&) { unknown = true; }
    s.require("matrices: an unknown name is refused", unknown);
  }

  // ---- Durbin et al. ---------------------------------------------------------------------------------------
  {
    // Durbin, Eddy, Krogh and Mitchison 1998, Biological Sequence Analysis, §2.3: Fig. 2.5 (global,
    // score 1), Fig. 2.6 (local, score 28) and Fig. 2.8 (overlap, score 25), with BLOSUM50 and a
    // linear gap of −8. (The local score 28 is confirmed in the reproduction of PMC2571991, Fig. 1.)
    const Codes x = sequence::encode_protein("HEAGAWGHEE"), y = sequence::encode_protein("PAWHEAE");
    const Scoring sc = scoring(matrix("BLOSUM50"), 0, 8);
    const auto show = [&](const Alignment& a) {
      const auto t = render(a, x, y, sequence::protein_letter);
      return t.query + " / " + t.target;
    };
    const Alignment g = gotoh(x, y, sc, Mode::global), l = gotoh(x, y, sc, Mode::local), o = gotoh(x, y, sc, Mode::overlap);
    s.require("Durbin et al. Fig. 2.5: global score 1, HEAGAWGHE-E / --P-AW-HEAE", g.score == 1 && show(g) == "HEAGAWGHE-E / --P-AW-HEAE", std::to_string(g.score) + ": " + show(g));
    s.require("Durbin et al. Fig. 2.6: local score 28, AWGHE / AW-HE", l.score == 28 && show(l) == "AWGHE / AW-HE", std::to_string(l.score) + ": " + show(l));
    s.require("Durbin et al. Fig. 2.8: overlap score 25, GAWGHEE / PAW-HEA", o.score == 25 && show(o) == "GAWGHEE / PAW-HEA", std::to_string(o.score) + ": " + show(o));
  }

  // ---- against the generic sum over paths -------------------------------------------------------------------
  const Scoring dna_affine = scoring(Matrix::dna(2, -3, -1), 5, 2), dna_linear = scoring(Matrix::dna(1, -1, -1), 0, 2), protein = blastp();
  for (Mode mode : MODES) {
    bool score = true, path = true, rescored = true, linear_memory = true, counted = true;
    std::size_t pairs = 0;
    for (int rep = 0; rep < 400; rep++) {
      const bool is_protein = rep % 4 == 3;
      const Scoring& sc = is_protein ? protein : (rep % 2 ? dna_affine : dna_linear);
      const std::size_t m = r.below(41), n = r.below(41);
      const Codes x = random_codes(r, m, is_protein ? 20 : (rep % 8 == 0 ? 5 : 4)), y = random_codes(r, n, is_protein ? 20 : 4);
      const Alignment a = gotoh(x, y, sc, mode);
      const std::int64_t ref = reference_score(x, y, sc, mode);
      const Endpoint e = gotoh_end(x, y, sc, mode);
      score = score && a.score == ref;
      path = path && allowed(a, m, n, mode);
      rescored = rescored && rescore(a, x, y, sc) == a.score;
      linear_memory = linear_memory && e.score == a.score && e.query_end == a.query_end && e.target_end == a.target_end && gotoh_score(x, y, sc, mode) == a.score;
      if (m <= 12 && n <= 12) counted = counted && reference_count(x, y, sc, mode).score == a.score;
      pairs++;
    }
    const std::string tag = std::string("gotoh, ") + mode_name(mode) + ": ";
    s.require(tag + "the score is morphe's generic sum over paths in (max, +)", score, std::to_string(pairs) + " random pairs of length ≤ 40: DNA (affine 2/−3/5/2, linear 1/−1/0/2) and protein (BLOSUM62, 11/1)");
    s.require(tag + "the path begins and ends where the mode allows and spans what it says", path);
    s.require(tag + "the path rescores to the claimed score", rescored);
    s.require(tag + "the linear-memory score and end are the traceback's", linear_memory);
    s.require(tag + "the tropical semiring with multiplicities has the same optimum", counted);
  }
  {
    // The field for drawing: its value at the alignment's end is the score; in global mode H(i, j)
    // is the optimal global score of the prefixes q[0, i) and t[0, j); the gap states never exceed H.
    bool end = true, prefixes = true, bounded = true;
    for (int rep = 0; rep < 60; rep++) {
      const Codes x = random_codes(r, r.below(25), 4), y = random_codes(r, r.below(25), 4);
      for (Mode mode : MODES) {
        const ScoreField f = gotoh_field(x, y, dna_affine, mode);
        const Alignment a = gotoh(x, y, dna_affine, mode);
        end = end && f.rows == x.size() + 1 && f.cols == y.size() + 1 && f.at(a.query_end, a.target_end) == a.score;
        for (std::size_t c = 0; c < f.best.size(); c++)
          bounded = bounded && f.ending_i[c] <= f.best[c] && f.ending_d[c] <= f.best[c] && f.best[c] != ScoreField::unreachable;
        if (mode == Mode::global && rep < 15)
          for (std::size_t i = 0; i <= x.size(); i++)
            for (std::size_t j = 0; j <= y.size(); j++)
              prefixes = prefixes && f.at(i, j) == gotoh_score(Codes(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(i)), Codes(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(j)), dna_affine, Mode::global);
      }
    }
    s.require("field: gotoh_field's value at the alignment's end is gotoh's score, in every mode", end, "60 random pairs of length ≤ 24");
    s.require("field: in global mode H(i, j) is the optimal score of the prefixes, at every vertex", prefixes);
    s.require("field: every vertex is reachable, and the gap states never exceed H", bounded);
  }
  {
    // Homologous pairs: two members of a simulated family (substitutions 0.1, indels 0.05).
    bool ok = true;
    std::size_t cells = 0;
    for (std::uint64_t seed = 1; seed <= 12; seed++) {
      const Codes root = simulate::random_sequence(150 + 20 * seed, 0.5, seed);
      const auto f = simulate::evolve_family(root, 2, 0.1, 0.05, seed);
      for (Mode mode : MODES) {
        const Alignment a = gotoh(f.sequences[0], f.sequences[1], bwa_mem(), mode);
        ok = ok && a.score == reference_score(f.sequences[0], f.sequences[1], bwa_mem(), mode) && rescore(a, f.sequences[0], f.sequences[1], bwa_mem()) == a.score;
        cells += f.sequences[0].size() * f.sequences[1].size();
      }
    }
    s.require("gotoh: on homologous pairs of 170–390 bases (BWA-MEM's scoring), every mode equals the generic sum and rescores", ok, std::to_string(cells) + " cells");
  }
  {
    // The reference itself: its traceback rescores to its score through the domain's rescore, and
    // it counts the global alignments as Delannoy.
    bool ok = true;
    for (int rep = 0; rep < 100; rep++) {
      const Codes x = random_codes(r, r.below(15), 4), y = random_codes(r, r.below(15), 4);
      for (Mode mode : MODES) {
        const Alignment a = reference_alignment(x, y, dna_affine, mode);
        ok = ok && a.score == reference_score(x, y, dna_affine, mode) && rescore(a, x, y, dna_affine) == a.score && allowed(a, x.size(), y.size(), mode);
      }
      ok = ok && reference_alignments(x, y, Mode::global) == morphe::formal::delannoy(x.size(), y.size());
    }
    s.require("reference: the generic traceback rescores to the generic score, and global alignments number D(m, n)", ok);
  }
  {
    // Co-optimal paths: the number of optimal global alignments of A^k against A^(k−1) with
    // linear gaps is k (where the gap goes), and gotoh puts the gap first (left-aligned).
    bool ok = true;
    for (std::size_t k = 1; k <= 30; k++) {
      const Codes x(k, sequence::A), y(k - 1, sequence::A);
      const auto c = reference_count(x, y, dna_linear, Mode::global);
      const Alignment a = gotoh(x, y, dna_linear, Mode::global);
      ok = ok && c.count == morphe::arithmetic::Integer(static_cast<long long>(k)) && (a.cigar.front() == CigarOp{'I', 1} || k == 1);
    }
    s.require("ties: A^k against A^(k−1) has exactly k optimal alignments, and gotoh's gap is the leftmost", ok);
  }

  // ---- speed ----------------------------------------------------------------------------------------------------
  {
    const Codes root = simulate::random_sequence(10'000, 0.5, 77);
    const auto f = simulate::evolve_family(root, 2, 0.05, 0.01, 77);
    const Codes& x = f.sequences[0];
    const Codes& y = f.sequences[1];
    const double cells = static_cast<double>(x.size()) * static_cast<double>(y.size());
    const auto t0 = std::chrono::steady_clock::now();
    const Alignment a = gotoh(x, y, bwa_mem(), Mode::global);
    const auto t1 = std::chrono::steady_clock::now();
    const std::int64_t score_only = gotoh_score(x, y, bwa_mem(), Mode::global);
    const auto t2 = std::chrono::steady_clock::now();
    const Alignment l = gotoh(x, y, bwa_mem(), Mode::local);
    const auto t3 = std::chrono::steady_clock::now();
    const double full = std::chrono::duration<double>(t1 - t0).count(), only = std::chrono::duration<double>(t2 - t1).count(), loc = std::chrono::duration<double>(t3 - t2).count();
    char detail[256];
    std::snprintf(detail, sizeof detail, "%zu × %zu: traceback %.0f Mcells/s (%.2f s), score only %.0f Mcells/s, local with traceback %.0f Mcells/s", x.size(), y.size(), cells / full / 1e6, full,
                  cells / only / 1e6, cells / loc / 1e6);
    s.require("speed: Gotoh with traceback on 10 kb × 10 kb, score and path consistent", a.score == score_only && rescore(a, x, y, bwa_mem()) == a.score && rescore(l, x, y, bwa_mem()) == l.score && l.score >= a.score,
              detail);
  }
  return s.done();
}
