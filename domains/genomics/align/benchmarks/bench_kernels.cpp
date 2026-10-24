// Measurements of the fast kernels of alignment (not a check: its numbers depend on the
// machine). Each kernel is timed against the scalar schedule of the same sum, on simulated
// sequences from fixed seeds, and its answer compared with that schedule's:
//   1. cell updates per second (GCUPS) of the striped Smith–Waterman kernel against Gotoh's
//      scalar recurrences, a 1 kb query against a 100 kb target; the inter-sequence kernel on
//      a batch; Myers' bit vectors on the same search, as unit-cost edit distance;
//   2. the wavefront algorithm (with its traceback, by its score alone, and BiWFA) against
//      Gotoh's O(nm) recurrences on 10 kb pairs at 1%, 5% and 10% divergence;
//   3. the memory of BiWFA's wavefronts against the WFA's on two 30 kb sequences at 1%.
//
// Build from the repository root, in WSL or with any g++ ≥ 11, on one line:
//   g++ -std=c++20 -O3 -march=native -I domains/genomics/align/include -I domains/genomics/sequence/include
//       -I domains/genomics/simulate/include -I foundations/calx/compute/backends/cpu/include
//       -I foundations/calx/numerics/random/include domains/genomics/align/benchmarks/bench_kernels.cpp -o bench_kernels

#include <genomics/align/bitparallel.hpp>
#include <genomics/align/gotoh.hpp>
#include <genomics/align/inter.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/align/striped.hpp>
#include <genomics/align/wavefront.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/cpu/target.hpp>
#include <calx/random/stream.hpp>

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::align;

namespace {

// Seconds of the fastest of `reps` runs of f.
template <class F>
double seconds(int reps, F&& f) {
  double best = 1e300;
  for (int k = 0; k < reps; k++) {
    const auto t0 = std::chrono::steady_clock::now();
    f();
    const auto t1 = std::chrono::steady_clock::now();
    best = std::min(best, std::chrono::duration<double>(t1 - t0).count());
  }
  return best;
}

// A copy with substitutions, insertions and deletions of 1 to 3 bases, a third each, at `rate`.
Codes diverged(const Codes& x, double rate, std::uint64_t seed) {
  calx::random::Stream r(seed, 0x42454e43, 0, 0);
  Codes y;
  for (std::size_t i = 0; i < x.size(); i++) {
    const double u = r.uniform();
    if (u < rate / 3) { y.push_back(simulate::substitute(r, x[i])); continue; }
    if (u < 2 * rate / 3) { i += r.below(3); continue; }
    if (u < rate)
      for (std::size_t k = 0, len = 1 + r.below(3); k < len; k++) y.push_back(static_cast<Code>(r.below(4)));
    y.push_back(x[i]);
  }
  return y;
}

}  // namespace

int main() {
  std::printf("host: %s\n\n", calx::cpu::describe().c_str());

  // ---- 1. GCUPS ----------------------------------------------------------------------------------------
  {
    const Scoring sc = minimap2_map_ont();
    const Codes target = simulate::random_sequence(100'000, 0.5, 1);
    Codes query(target.begin() + 40'000, target.begin() + 41'000);
    query = diverged(query, 0.1, 2);
    query.resize(1000, sequence::A);
    const double cells = static_cast<double>(query.size()) * static_cast<double>(target.size());
    Endpoint scalar_end, striped_e;
    const double t_scalar = seconds(3, [&] { scalar_end = gotoh_end(query, target, sc, Mode::local); });
    const StripedProfile profile(query, sc);
    StripedStats st;
    const double t_striped = seconds(5, [&] { st = {}; striped_e = profile.end(target, &st); });
    std::printf("Smith-Waterman, 1 kb query x 100 kb target (score %lld, %s)\n", static_cast<long long>(striped_e.score),
                striped_e.score == scalar_end.score && striped_e.target_end == scalar_end.target_end ? "same end as Gotoh" : "DIFFERENT from Gotoh");
    std::printf("  Gotoh, scalar             %8.3f s  %7.3f GCUPS\n", t_scalar, cells / t_scalar / 1e9);
    std::printf("  striped, %2d-bit lanes      %8.4f s  %7.3f GCUPS  (x%.0f; the 8-bit pass first, %llu lazy-F corrections)\n", st.bits, t_striped,
                cells / t_striped / 1e9, t_scalar / t_striped, static_cast<unsigned long long>(st.lazy_f));
    // Against an unrelated target the score stays in a byte: the 8-bit pass alone.
    const Codes other = simulate::random_sequence(100'000, 0.5, 4);
    StripedStats st8;
    Endpoint e8;
    const double t8 = seconds(5, [&] { st8 = {}; e8 = profile.end(other, &st8); });
    std::printf("  striped, unrelated target %8.4f s  %7.3f GCUPS  (%d-bit lanes, score %lld)\n", t8, cells / t8 / 1e9, st8.bits, static_cast<long long>(e8.score));
    // The inter-sequence kernel: the same query against 320 targets of 10 kb.
    std::vector<Codes> batch;
    for (std::uint32_t k = 0; k < 320; k++) batch.push_back(simulate::random_sequence(10'000, 0.5, 3, k));
    InterStats is;
    std::vector<std::int64_t> scores;
    const double t_inter = seconds(3, [&] { is = {}; scores = inter_scores(query, batch, sc, &is); });
    const double batch_cells = static_cast<double>(query.size()) * 320.0 * 10'000.0;
    std::printf("  inter-sequence, 16-bit    %8.4f s  %7.3f GCUPS  (320 targets of 10 kb, one per lane)\n", t_inter, batch_cells / t_inter / 1e9);
    // Myers' bit vectors on the same search, unit costs, every end within 150 edits.
    const MyersPattern mp(query);
    std::size_t hits = 0;
    MyersStats ms;
    const double t_myers = seconds(5, [&] { ms = {}; hits = myers_search(mp, target, 150, &ms).size(); });
    std::printf("  Myers, edit distance <=150 %7.4f s  %7.3f GCUPS  (%zu ends; %.1f%% of the cells in the band)\n\n", t_myers, cells / t_myers / 1e9, hits,
                100.0 * static_cast<double>(ms.cells()) / cells);
  }

  // ---- 2. WFA against the O(nm) recurrences ----------------------------------------------------------------
  {
    const Penalties p{4, 6, 2};
    Scoring sc;
    sc.matrix = Matrix::dna(0, -4, -4);
    sc.gap_open = 6;
    sc.gap_extend = 2;
    std::printf("Global alignment of 10 kb pairs, penalties (4, 6, 2)\n");
    std::printf("  divergence  penalty  Gotoh O(nm)   WFA+path    WFA score   BiWFA      peak memory WFA / BiWFA\n");
    for (double rate : {0.01, 0.05, 0.10}) {
      const Codes q = simulate::random_sequence(10'000, 0.5, 11);
      const Codes t = diverged(q, rate, 12);
      std::int64_t dp = 0, w = 0, ws = 0, wb = 0;
      const double t_dp = seconds(1, [&] { dp = -gotoh_score(q, t, sc, Mode::global); });
      WavefrontStats full, bi;
      const double t_w = seconds(3, [&] { full = {}; w = -wavefront(q, t, p, Mode::global, {}, &full)->score; });
      const double t_ws = seconds(3, [&] { ws = wavefront_penalty(q, t, p); });
      const double t_b = seconds(3, [&] { bi = {}; wb = -biwfa(q, t, p, &bi).score; });
      std::printf("  %5.0f%%     %7lld  %9.3f s  %8.4f s  %8.4f s  %8.4f s  %7.1f MB / %.2f MB%s\n", rate * 100, static_cast<long long>(dp), t_dp, t_w, t_ws, t_b,
                  static_cast<double>(full.peak_bytes) / 1e6, static_cast<double>(bi.peak_bytes) / 1e6, dp == w && w == ws && ws == wb ? "" : "  MISMATCH");
    }
    std::printf("\n");
  }

  // ---- 3. BiWFA's memory ------------------------------------------------------------------------------------
  {
    const Penalties p{4, 6, 2};
    const Codes q = simulate::random_sequence(30'000, 0.5, 21);
    const Codes t = diverged(q, 0.01, 22);
    WavefrontStats full, bi;
    std::int64_t w = 0, b = 0;
    const double t_w = seconds(1, [&] { full = {}; w = -wavefront(q, t, p, Mode::global, {}, &full)->score; });
    const double t_b = seconds(1, [&] { bi = {}; b = -biwfa(q, t, p, &bi).score; });
    std::printf("Two 30 kb sequences at 1%%: penalty %lld (BiWFA %lld)\n", static_cast<long long>(w), static_cast<long long>(b));
    std::printf("  WFA   %8.4f s, wavefronts at their peak %9.3f MB\n", t_w, static_cast<double>(full.peak_bytes) / 1e6);
    std::printf("  BiWFA %8.4f s, wavefronts at their peak %9.3f MB\n", t_b, static_cast<double>(bi.peak_bytes) / 1e6);
  }
  return 0;
}
