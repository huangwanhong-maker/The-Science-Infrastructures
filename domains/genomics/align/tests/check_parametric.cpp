// Checks of parametric.hpp: on short random DNA pairs, at random (α, β) the optimal score of the
// parametric polygon is the brute-force optimum over every global alignment, and the vertex it
// selects has that parameter in its region; every witness is a global alignment with its vertex's
// counts; at an interior direction of each region the domain's own Gotoh (match 0, mismatch −A,
// gap −B per column) returns exactly that vertex; the regions cover the parameter plane; protein
// pairs; and the size of the polygon of homologous pairs of a few hundred bases.

#include <genomics/align/gotoh.hpp>
#include <genomics/align/parametric.hpp>
#include <genomics/simulate/family.hpp>
#include <genomics/simulate/genome.hpp>

#include <morphe/formal/alignment.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <array>
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

// Mismatches and gap columns of a global alignment, from its CIGAR.
std::pair<std::int64_t, std::int64_t> counts(const Alignment& a, const Codes& q, const Codes& t, Code limit) {
  std::int64_t x = 0, g = 0;
  std::size_t i = a.query_begin, j = a.target_begin;
  for (const auto& o : a.cigar)
    for (std::uint32_t k = 0; k < o.length; k++) {
      if (o.op == 'M') {
        x += !same(q[i], t[j], limit);
        i++;
        j++;
      } else {
        g++;
        (o.op == 'I' ? i : j)++;
      }
    }
  return {x, g};
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x70617261, 0, 0);
  {
    bool brute = true, region = true, witness = true, cover = true;
    std::size_t vertices = 0;
    for (int rep = 0; rep < 60; rep++) {
      const Codes q = random_codes(r, r.below(8), 4), t = random_codes(r, r.below(8), 4);
      const Parametric P = parametric(q, t);
      vertices += P.vertices.size();
      for (const auto& v : P.vertices) {
        const auto [x, g] = counts(v.alignment, q, t, 4);
        witness = witness && x == v.mismatches && g == v.gaps && v.alignment.query_end == q.size() && v.alignment.target_end == t.size() &&
                  2 * (v.matches + v.mismatches) + v.gaps == static_cast<std::int64_t>(q.size() + t.size());
      }
      // Every global alignment's (matches, mismatches, gaps), by enumeration.
      std::vector<std::array<std::int64_t, 3>> all;
      morphe::formal::enumerate_paths(q.size(), t.size(), morphe::formal::Mode::global, [&](const morphe::formal::EditPath& path) {
        std::int64_t mt = 0, mm = 0, gp = 0;
        std::size_t i = 0, j = 0;
        for (char op : path.ops) {
          if (op == 'M') { (q[i] == t[j] ? mt : mm)++; i++; j++; }
          else { gp++; (op == 'I' ? i : j)++; }
        }
        all.push_back({mt, mm, gp});
      });
      for (int k = 0; k < 20; k++) {
        const double alpha = 3 * r.uniform(), beta = 3 * r.uniform();
        double best = -1e300;
        for (const auto& a : all) best = std::max(best, a[0] - alpha * a[1] - beta * a[2]);
        const std::size_t o = P.optimal(alpha, beta);
        brute = brute && std::fabs(P.score(o, alpha, beta) - best) < 1e-9;
        region = region && P.vertices[o].optimal_at(alpha, beta);
        bool covered = false;
        for (const auto& v : P.vertices) covered = covered || v.optimal_at(alpha, beta);
        cover = cover && covered;
      }
    }
    s.require("parametric: at random (α, β) ∈ [0, 3]² the polygon's optimum is the best of every enumerated global alignment", brute, "60 random DNA pairs of length ≤ 7, " + std::to_string(vertices) + " vertices");
    s.require("parametric: the selected vertex has the parameter in its region", region);
    s.require("parametric: every witness is a global alignment with its vertex's mismatches and gaps", witness);
    s.require("parametric: the regions cover the parameters", cover);
  }
  {
    // Each vertex is the unique optimum of Gotoh at an interior direction of its region: minimising
    // A·mismatches + B·gaps is the scoring match 0, mismatch −A, gap −B per column.
    bool unique = true;
    std::size_t tried = 0;
    for (int rep = 0; rep < 30; rep++) {
      const Codes root = random_codes(r, 40 + r.below(40), 4);
      const auto f = simulate::evolve_family(root, 2, 0.15, 0.05, 300 + static_cast<std::uint64_t>(rep));
      const Codes& q = f.sequences[0];
      const Codes& t = f.sequences[1];
      const Parametric P = parametric(q, t);
      for (const auto& v : P.vertices) {
        if (v.kind != morphe::formal::NormalCone::Kind::sector) continue;
        const std::int64_t A = v.first.x + v.second.x, B = v.first.y + v.second.y;  // inside the region's cone of c
        if (A < 0 || B < 0) continue;  // Gotoh needs nonnegative penalties
        Scoring sc;
        sc.matrix = Matrix::dna(0, static_cast<int>(-A), static_cast<int>(-A));
        sc.gap_open = 0;
        sc.gap_extend = static_cast<int>(B);
        const Alignment a = gotoh(q, t, sc, Mode::global);
        const auto [x, g] = counts(a, q, t, 4);
        unique = unique && x == v.mismatches && g == v.gaps && a.score == -(A * v.mismatches + B * v.gaps);
        tried++;
      }
    }
    s.require("parametric: at an interior direction of each region Gotoh returns that vertex, and only it", unique && tried > 0, std::to_string(tried) + " vertices of 30 homologous pairs of 40–80 bases");
  }
  {
    // Protein: B, Z, X and * never match.
    const Codes q = sequence::encode_protein("HEAGAWGHEE"), t = sequence::encode_protein("PAWHEAE");
    const Parametric P = parametric(q, t, false);
    bool ok = !P.vertices.empty();
    for (const auto& v : P.vertices) {
      const auto [x, g] = counts(v.alignment, q, t, 20);
      ok = ok && x == v.mismatches && g == v.gaps;
    }
    // At α = β = 0 the score is the number of matches, so the optimum is a longest common
    // subsequence: AWHEE, of length 5 (P is not in the first, and no A follows its last H).
    std::vector<std::vector<int>> L(q.size() + 1, std::vector<int>(t.size() + 1, 0));
    for (std::size_t i = 1; i <= q.size(); i++)
      for (std::size_t j = 1; j <= t.size(); j++) L[i][j] = q[i - 1] == t[j - 1] ? L[i - 1][j - 1] + 1 : std::max(L[i - 1][j], L[i][j - 1]);
    s.require("parametric: protein HEAGAWGHEE and PAWHEAE, witnesses consistent and the optimum at (0, 0) a longest common subsequence (AWHEE, 5)",
              ok && P.vertices[P.optimal(0, 0)].matches == 5 && L[q.size()][t.size()] == 5, std::to_string(P.vertices.size()) + " vertices");
  }
  {
    std::string sizes;
    bool ok = true;
    for (std::size_t n : {100, 200, 400}) {
      const Codes root = simulate::random_sequence(n, 0.5, n);
      const auto f = simulate::evolve_family(root, 2, 0.1, 0.03, n);
      const auto t0 = std::chrono::steady_clock::now();
      const Parametric P = parametric(f.sequences[0], f.sequences[1]);
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      char buf[96];
      std::snprintf(buf, sizeof buf, "%s%zu vertices at n = %zu (%.2f s)", sizes.empty() ? "" : ", ", P.vertices.size(), n, secs);
      sizes += buf;
      // A convex lattice polygon in an N × N box has at most about 12 (4π²)^{−1/3} N^{2/3} ≈ 3.52 N^{2/3}
      // vertices (Acketa and Žunić 1995); here N = m + n, and the margin is 4 N^{2/3}.
      ok = ok && static_cast<double>(P.vertices.size()) <= 4 * std::cbrt(std::pow(static_cast<double>(f.sequences[0].size() + f.sequences[1].size()), 2));
    }
    s.require("parametric: homologous pairs of 100–400 bases have polygons within the lattice bound, 4 (m + n)^{2/3} vertices", ok, sizes);
  }
  return s.done();
}
