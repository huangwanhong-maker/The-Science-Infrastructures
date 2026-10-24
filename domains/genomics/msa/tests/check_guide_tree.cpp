// Checks of guide trees: the k-mer profile against a count by brute force; MUSCLE's fractional
// common k-mer count on a hand-worked pair (F = 2/4); Dayhoff's six classes making residues of
// one group indistinguishable; the distance matrix the same on 1 and 8 threads; the corrected
// distances against their closed forms; UPGMA on the worked example of 5S rRNA distances
// (Erdmann and Wolters 1986, as the UPGMA article of Wikipedia works it: merges at 8.5, 11, 14,
// 16.5); the nearest-neighbour chain against the greedy algorithm on random matrices and on
// points of the plane (the same merges, sizes and distances); ultrametric matrices giving back
// their trees; CLUSTAL W's weights on a tree worked by hand.

#include <genomics/msa/guide_tree.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace genomics::msa;
namespace seq = genomics::sequence;

namespace {

bool same_dendrogram(const Dendrogram& a, const Dendrogram& b, double tolerance) {
  if (a.leaves != b.leaves || a.merges.size() != b.merges.size()) return false;
  for (std::size_t k = 0; k < a.merges.size(); k++) {
    const Merge &x = a.merges[k], &y = b.merges[k];
    if (x.a != y.a || x.b != y.b || x.size != y.size || std::fabs(x.distance - y.distance) > tolerance * std::max(1.0, std::fabs(y.distance))) return false;
  }
  return true;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x47545245, 0, 0);

  // ---- k-mers ------------------------------------------------------------------------------------
  {
    bool counted = true;
    for (int t = 0; t < 100 && counted; t++) {
      Codes x(r.below(200));
      for (auto& c : x) c = r.uniform() < 0.02 ? seq::N : static_cast<Code>(r.below(4));
      const int k = 1 + static_cast<int>(r.below(8));
      std::map<std::string, std::uint32_t> brute;
      std::size_t total = 0;
      for (std::size_t i = 0; i + static_cast<std::size_t>(k) <= x.size(); i++) {
        bool ok = true;
        std::string w;
        for (int j = 0; j < k; j++) {
          ok = ok && x[i + static_cast<std::size_t>(j)] < 4;
          w += seq::dna_letter(x[i + static_cast<std::size_t>(j)]);
        }
        if (ok) { brute[w]++; total++; }
      }
      const KmerProfile p = kmer_profile(x, k, dna_letters());
      counted = p.total == total && p.kmers.size() == brute.size();
      std::size_t i = 0;
      for (const auto& [w, n] : brute) counted = counted && p.counts[i++] == n;  // both in lexicographic order
    }
    s.require("the k-mer profile of 100 random sequences with Ns equals a count by brute force (k = 1 to 8)", counted);

    // ACGTACGT has ACG, CGT twice, GTA, TAC; ACGTTT has ACG, CGT, GTT, TTT: shared 2 of the shorter's 4.
    const auto x = kmer_profile(seq::encode_dna("ACGTACGT"), 3, dna_letters()), y = kmer_profile(seq::encode_dna("ACGTTT"), 3, dna_letters());
    s.check("the fractional common 3-mer count of ACGTACGT and ACGTTT is 2/4", kmer_similarity(x, y) - 0.5, 1e-15);
    s.check("a sequence shares all its k-mers with itself", kmer_similarity(x, x) - 1.0, 0);
    // S↔T and D↔E lie within Dayhoff groups.
    const Codes p1 = seq::encode_protein("MSDKLLVTWAGHIRE"), p2 = seq::encode_protein("MTEKLLVSWAGHIRD");
    s.require("residues of one Dayhoff group are one class: F = 1 on Dayhoff(6), below 1 on twenty letters",
              kmer_similarity(kmer_profile(p1, 3, dayhoff6()), kmer_profile(p2, 3, dayhoff6())) == 1.0 &&
                  kmer_similarity(kmer_profile(p1, 3, protein_letters()), kmer_profile(p2, 3, protein_letters())) < 1.0);

    std::vector<Codes> seqs;
    for (int i = 0; i < 37; i++) {
      Codes z(100 + r.below(200));
      for (auto& c : z) c = static_cast<Code>(r.below(4));
      seqs.push_back(z);
    }
    calx::scheduler::ThreadPool one(1), eight(8);
    const auto d1 = kmer_distances(seqs, 4, dna_letters(), one), d8 = kmer_distances(seqs, 4, dna_letters(), eight);
    bool symmetric = true;
    for (std::size_t i = 0; i < d1.n; i++)
      for (std::size_t j = 0; j < d1.n; j++) symmetric = symmetric && d1(i, j) == d1(j, i) && d1(i, j) >= 0 && d1(i, j) <= 1 && (i != j || d1(i, j) == 0);
    s.require("the k-mer distances of 37 sequences are symmetric, in [0, 1], 0 on the diagonal, and the same on 1 and 8 threads", symmetric && d1.d == d8.d);
  }

  // ---- corrected distances ------------------------------------------------------------------------------------
  {
    s.check("Jukes–Cantor at p = 0.3: −(3/4) ln(1 − 0.4) = 0.3831192", jukes_cantor(0.3) - 0.3831192, 1e-7);
    s.check("Kimura's protein distance at p = 0.3: −ln(1 − 0.3 − 0.018) = 0.3827256", kimura_protein(0.3) - 0.3827256, 1e-7);
    s.require("both vanish at p = 0 and are held finite beyond saturation", jukes_cantor(0) == 0 && kimura_protein(0) == 0 && std::isfinite(jukes_cantor(0.9)) &&
                                                                              std::isfinite(kimura_protein(1.0)));
    Msa m;
    m.names = {"a", "b"};
    m.rows = {seq::encode_dna("ACGTACGTAC"), seq::encode_dna("ACGTTCGTA-")};
    s.check("the p-distance counts columns where both rows have a base: 1 of 9", p_distance(m.rows[0], m.rows[1], 4) - 1.0 / 9, 1e-15);
  }

  // ---- UPGMA ------------------------------------------------------------------------------------
  {
    // 5S rRNA distances of a Bacillus subtilis, b B. stearothermophilus, c Lactobacillus viridescens,
    // d Acholeplasma modicum, e Micrococcus luteus.
    DistanceMatrix d(5);
    const double v[5][5] = {{0, 17, 21, 31, 23}, {17, 0, 30, 34, 21}, {21, 30, 0, 28, 39}, {31, 34, 28, 0, 43}, {23, 21, 39, 43, 0}};
    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++) d.at(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) = v[i][j];
    const Dendrogram t = upgma(d);
    const bool ok = t.merges.size() == 4 && t.merges[0] == Merge{0, 1, 17, 2} && t.merges[1] == Merge{4, 5, 22, 3} && t.merges[2] == Merge{2, 3, 28, 2} &&
                    t.merges[3] == Merge{6, 7, 33, 5};
    s.require("UPGMA of the 5S rRNA distances merges (a, b) at height 8.5, then e at 11, (c, d) at 14, and the root at 16.5", ok);
    s.require("its Newick string", t.newick({"a", "b", "c", "d", "e"}) == "((e:11,(a:8.5,b:8.5):2.5):5.5,(c:14,d:14):2.5);", t.newick({"a", "b", "c", "d", "e"}));
    s.require("the greedy algorithm gives the same dendrogram", same_dendrogram(upgma_naive(d), t, 0));

    bool chain = true, plane = true;
    for (int trial = 0; trial < 200; trial++) {
      const std::size_t n = 2 + r.below(70);
      DistanceMatrix m(n), e(n);
      std::vector<double> px(n), py(n);
      for (std::size_t i = 0; i < n; i++) { px[i] = r.uniform(); py[i] = r.uniform(); }
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = i + 1; j < n; j++) {
          m.at(i, j) = m.at(j, i) = r.uniform();
          e.at(i, j) = e.at(j, i) = std::hypot(px[i] - px[j], py[i] - py[j]);
        }
      chain = chain && same_dendrogram(upgma(m), upgma_naive(m), 1e-12);
      plane = plane && same_dendrogram(upgma(e), upgma_naive(e), 1e-12);
    }
    s.require("the nearest-neighbour chain equals the greedy algorithm on 200 random matrices of 2 to 71 points", chain);
    s.require("and on 200 sets of random points of the plane", plane);

    // An ultrametric from a random tree: UPGMA gives its clusters back, at their heights.
    bool ultra = true;
    for (int trial = 0; trial < 50; trial++) {
      const std::size_t n = 3 + r.below(40);
      // Random merges at increasing heights.
      std::vector<std::vector<std::size_t>> clusters;
      for (std::size_t i = 0; i < n; i++) clusters.push_back({i});
      DistanceMatrix m(n);
      double h = 0;
      std::vector<std::vector<std::size_t>> made;
      while (clusters.size() > 1) {
        h += 0.1 + r.uniform();
        const std::size_t a = r.below(clusters.size());
        std::size_t b = r.below(clusters.size() - 1);
        if (b >= a) b++;
        for (std::size_t x : clusters[a])
          for (std::size_t y : clusters[b]) m.at(x, y) = m.at(y, x) = 2 * h;
        auto merged = clusters[a];
        merged.insert(merged.end(), clusters[b].begin(), clusters[b].end());
        std::sort(merged.begin(), merged.end());
        made.push_back(merged);
        clusters.erase(clusters.begin() + static_cast<std::ptrdiff_t>(std::max(a, b)));
        clusters.erase(clusters.begin() + static_cast<std::ptrdiff_t>(std::min(a, b)));
        clusters.push_back(merged);
      }
      const Dendrogram t = upgma(m);
      for (std::size_t k = 0; k < t.merges.size(); k++) {
        auto c = t.leaves_under(n + k);
        std::sort(c.begin(), c.end());
        ultra = ultra && c == made[k];
      }
    }
    s.require("UPGMA of 50 ultrametrics gives back the trees they were made from, merge by merge", ultra);
  }

  // ---- CLUSTAL W weights ------------------------------------------------------------------------------------
  {
    // ((A:1, B:1):2, C:3): A and B carry 1 + 2/2, C carries 3; normalised to mean 1: 6/7, 6/7, 9/7.
    DistanceMatrix d(3);
    d.at(0, 1) = d.at(1, 0) = 2;
    d.at(0, 2) = d.at(2, 0) = d.at(1, 2) = d.at(2, 1) = 6;
    const auto w = clustal_weights(upgma(d));
    s.check("CLUSTAL W's weights on ((A:1, B:1):2, C:3) are 6/7, 6/7, 9/7", std::fabs(w[0] - 6.0 / 7) + std::fabs(w[1] - 6.0 / 7) + std::fabs(w[2] - 9.0 / 7), 1e-14);
    DistanceMatrix z(4);
    const auto u = clustal_weights(upgma(z));
    s.require("a tree of zero branches gives every sequence weight 1", u == std::vector<double>(4, 1.0));
  }
  return s.done();
}
