// Checks of the period-3 signal on a real genome: when the data cache has Escherichia coli K-12
// MG1655 (NC_000913.3), windows of 120 bases every 3 bases are divided into those lying wholly in
// the annotated CDSs (pseudogenes left out) and those touching no annotated gene of any kind, and
// the relative strength P = S(W/3)/S̄ is held to what it was measured to be (2026-10-24): coding
// windows stronger on average by more than twice, and above intergenic ones in more than 0.83 of
// pairs (the area under the ROC curve, measured 0.857); intergenic windows near 1 − Σ p_b², the
// value for random sequence of their composition (measured 0.744 against 0.741). Longer windows
// separate the two more (W = 351: a mean ratio of 4.6, AUC 0.972), reported.

#include <genomics/sequence/genbank.hpp>
#include <genomics/signal/spectrum.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace genomics;

namespace {

std::string fixed(double v, int digits = 3) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*f", digits, v);
  return b;
}

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

// The probability that a coding window outscores an intergenic one, ties counting one half
// (the Mann–Whitney statistic over the product of the two sizes).
double auc(const std::vector<double>& in, const std::vector<double>& out) {
  std::vector<std::pair<double, int>> all;
  all.reserve(in.size() + out.size());
  for (double v : in) all.push_back({v, 1});
  for (double v : out) all.push_back({v, 0});
  std::sort(all.begin(), all.end());
  double ranks = 0;
  for (std::size_t i = 0; i < all.size();) {
    std::size_t j = i;
    while (j < all.size() && all[j].first == all[i].first) j++;
    const double r = static_cast<double>(i + j + 1) / 2;
    for (std::size_t t = i; t < j; t++)
      if (all[t].second) ranks += r;
    i = j;
  }
  const double a = static_cast<double>(in.size()), b = static_cast<double>(out.size());
  return (ranks - a * (a + 1) / 2) / (a * b);
}

double mean(const std::vector<double>& v) {
  double s = 0;
  for (double x : v) s += x;
  return v.empty() ? 0 : s / static_cast<double>(v.size());
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const char* env = std::getenv("SCIENCE_CACHE");
  const std::filesystem::path cache = env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
  const std::filesystem::path file = cache / "ncbi-genomes" / "NC_000913.3.gb";
  if (!std::filesystem::exists(file)) {
    s.skip("Escherichia coli K-12 (data/cache/ncbi-genomes/NC_000913.3.gb) is not in the data cache: its period-3 signal is not measured");
    return s.done();
  }
  const auto recs = sequence::read_genbank(read_file(file));
  const sequence::Codes genome = sequence::encode_dna(recs.at(0).sequence);
  const std::size_t n = genome.size();
  // Bases in a CDS (not a pseudogene), and bases in any gene (gene, CDS, rRNA, tRNA, ncRNA, ...).
  std::vector<std::uint32_t> coding(n + 1, 0), genic(n + 1, 0);
  {
    std::vector<char> c(n, 0), g(n, 0);
    for (const auto& f : recs[0].features) {
      const bool cds = f.key == "CDS";
      if (!(cds || f.key == "gene" || f.key == "rRNA" || f.key == "tRNA" || f.key == "ncRNA" || f.key == "tmRNA" || f.key == "misc_RNA")) continue;
      for (const auto& p : f.location.parts)
        for (std::size_t i = p.begin; i < p.end && i < n; i++) {
          if (cds && !f.has("pseudo")) c[i] = 1;
          g[i] = 1;
        }
    }
    for (std::size_t i = 0; i < n; i++) {
      coding[i + 1] = coding[i] + static_cast<std::uint32_t>(c[i]);
      genic[i + 1] = genic[i] + static_cast<std::uint32_t>(g[i]);
    }
  }
  for (std::size_t W : {120, 351}) {
    const auto profile = signal::period3_profile(genome, W, 3);
    std::vector<double> in, out;
    std::array<double, 4> composition{};
    for (const auto& w : profile) {
      const std::size_t c = coding[w.begin + W] - coding[w.begin], g = genic[w.begin + W] - genic[w.begin];
      if (c == W) {
        in.push_back(w.strength());
      } else if (g == 0) {
        out.push_back(w.strength());
        for (std::size_t i = w.begin; i < w.begin + W; i++)
          if (genome[i] < 4) composition[genome[i]] += 1;
      }
    }
    double total = 0, sum_sq = 0;
    for (double v : composition) total += v;
    for (double v : composition) sum_sq += (v / total) * (v / total);
    const double mi = mean(in), mo = mean(out), a = auc(in, out);
    const std::string head = "E. coli K-12, windows of " + std::to_string(W) + " every 3 bases (" + std::to_string(in.size()) + " in CDSs, " + std::to_string(out.size()) + " intergenic)";
    const std::string detail = "mean P " + fixed(mi) + " in CDSs, " + fixed(mo) + " intergenic, ratio " + fixed(mi / mo) + ", AUC " + fixed(a, 4) + "; 1 − Σ p² of the intergenic windows " + fixed(1 - sum_sq);
    if (W == 120) {
      s.require(head + ": the period-3 strength is higher in CDSs, on average more than twice", mi > 2 * mo, detail);
      s.require(head + ": a CDS window outscores an intergenic one in more than 0.83 of pairs", a > 0.83, detail);
      s.check(head + ": intergenic windows near 1 − Σ p_b², the value for random sequence of their composition", mo - (1 - sum_sq), 0.03);
    } else {
      s.require(head + ": the separation grows with the window (reported; ratio above 4, AUC above 0.96)", mi > 4 * mo && a > 0.96, detail);
    }
  }
  return s.done();
}
