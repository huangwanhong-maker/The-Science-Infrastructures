// Checks of the all-against-all distance matrix: 40 genomes of a star phylogeny (a root of 200 kb
// and leaves diverged from it by substitutions at rates up to 4 per cent), sketched on the pool;
// the matrix symmetric with a zero diagonal, each entry the Mash distance of its two sketches, the
// same to the bit on one and eight threads and for blocks of 1, 7 and 32 genomes, with one task a
// block on or above the diagonal in the trace; and the distances of the star additive, D(i, j)
// ≈ D(i, root) + D(j, root), as distances along a tree are.

#include <genomics/sequence/alphabet.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/sketch/distance.hpp>
#include <genomics/sketch/minhash.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;
using sequence::Codes;

namespace {

constexpr std::uint32_t CHECK = 0x44495354;

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool one(1), eight(8);

  // The root is genome 0; leaf i is the root with each base substituted under a Poisson process of
  // rate d_i, so that the Mash distance from the root estimates d_i and between leaves d_i + d_j.
  const std::size_t n = 40, L = 200000;
  const Codes root = simulate::random_sequence(L, 0.5, 31);
  std::vector<std::vector<Codes>> genomes{{root}};
  std::vector<double> rate{0};
  for (std::size_t i = 1; i < n; i++) {
    const double d = 0.002 + 0.038 * static_cast<double>(i - 1) / static_cast<double>(n - 2);
    calx::random::Stream r(31, CHECK, static_cast<std::uint32_t>(i), 0);
    Codes x = root;
    for (auto& c : x)
      if (r.uniform() < -std::expm1(-d)) c = simulate::substitute(r, c);
    genomes.push_back({x});
    rate.push_back(d);
  }
  const auto t0 = std::chrono::steady_clock::now();
  const auto sketches = sketch::minhash_all(eight, genomes, 21, 2000);
  const double sketching = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  const auto alone = sketch::minhash_all(one, genomes, 21, 2000);
  bool pooled = true;
  for (std::size_t i = 0; i < n; i++) pooled = pooled && sketches[i].hashes == alone[i].hashes && sketches[i].kmers == alone[i].kmers;
  char took[64];
  std::snprintf(took, sizeof took, "%zu genomes in %.3f s on eight threads", n, sketching);
  s.require("sketches made on eight threads equal sketches made on one", pooled, took);

  const auto ref = sketch::mash_distances(one, sketches, 32);
  bool symmetric = true, diagonal = true, entries = true;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      symmetric = symmetric && ref(i, j) == ref(j, i);
      if (i == j) diagonal = diagonal && ref(i, i) == 0;
      entries = entries && ref(i, j) == sketch::compare(sketches[i], sketches[j]).distance;
    }
  s.require("the matrix is symmetric", symmetric);
  s.require("its diagonal is zero", diagonal);
  s.require("each entry is the Mash distance of its two sketches", entries);

  bool same = true, tasks = true;
  std::string detail;
  for (std::size_t block : {1, 7, 32}) {
    for (auto* pool : {&one, &eight}) {
      const auto d = sketch::mash_distances(*pool, sketches, block);
      same = same && d.values == ref.values;
      const std::size_t b = (n + block - 1) / block;
      tasks = tasks && d.tasks == b * (b + 1) / 2 && d.trace.records.size() == d.tasks && d.trace.threads == pool->size();
      if (pool == &eight) {
        char buf[96];
        std::snprintf(buf, sizeof buf, "%sblock %zu: %zu tasks, utilisation %.2f", detail.empty() ? "" : "; ", block, d.tasks, d.trace.utilisation());
        detail += buf;
      }
    }
  }
  s.require("the matrix is the same to the bit on one and eight threads, for blocks of 1, 7 and 32", same);
  s.require("the trace holds one task for each block on or above the diagonal", tasks, detail);

  // Additivity on the star: D(i, j) ≈ D(i, 0) + D(j, 0), within the sketch's error at D(i, j).
  double worst = 0;
  for (std::size_t i = 1; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) worst = std::max(worst, std::fabs(ref(i, j) - (ref(i, 0) + ref(j, 0))));
  s.check("the distances of a star phylogeny are additive: the largest |D(i, j) − D(i, root) − D(j, root)|", worst, 0.02);
  double rates = 0;
  for (std::size_t i = 1; i < n; i++) rates = std::max(rates, std::fabs(ref(i, 0) - rate[i]));
  s.check("the distance of each leaf from the root estimates its rate d_i (largest error)", rates, 0.01);

  return s.done();
}
