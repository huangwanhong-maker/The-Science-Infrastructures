// Checks of sequence-to-expression learning: the k-mer baseline and the convolutional network
// recover a planted regulatory code from simulated sequences (truth known), the network's filters
// are matched to the planted motifs; MEME's motif format, the motif matcher and the grouped split;
// and, when the data cache holds them, a stated subset of the Random Promoter DREAM Challenge's test
// set (data/registry/dream-promoters.json) learnt and evaluated on held-out families, the filters
// matched to JASPAR 2024's fungal profiles (jaspar-fungi.json).

#include <genomics/regulation/learn.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <optional>
#include <set>
#include <sstream>
#include <string>

using namespace genomics;
using namespace genomics::regulation;

namespace {

std::optional<std::string> read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  if (!f) return std::nullopt;
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

std::string fixed(double v, int d = 3) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*f", d, v);
  return b;
}

std::vector<double> matrix_of(const std::string& word) {
  std::vector<double> p(4 * word.size(), 0.0);
  for (std::size_t i = 0; i < word.size(); i++) p[4 * i + sequence::dna_code(word[i])] = 1;
  return p;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool pool(4);

  // ---- MEME's format and the matcher ------------------------------------------------------------------
  const std::string meme =
      "MEME version 4\n\nALPHABET= ACGT\n\nMOTIF M1 GATA\nletter-probability matrix: alength= 4 w= 4 nsites= 10 E= 0\n"
      " 0 0 1 0\n 1 0 0 0\n 0 0 0 1\n 1 0 0 0\n\nMOTIF M2 CACGTG\nletter-probability matrix: alength= 4 w= 6 nsites= 10 E= 0\n"
      " 0 1 0 0\n 1 0 0 0\n 0 1 0 0\n 0 0 1 0\n 0 0 0 1\n 0 0 1 0\nURL x\n";
  const auto known = read_meme_motifs(meme);
  s.require("MEME's format: two motifs, their ids, names and widths", known.size() == 2 && known[0].id == "M1" && known[1].name == "CACGTG" && known[1].width == 6);
  {
    const auto self = match_motif(matrix_of("TTCACGTGAA"), 10, known, 5, 50);
    s.require("a query holding CACGTG matches M2 with similarity 1 at offset 2 (p " + fixed(self.p_value) + ")", self.motif == 1 && std::fabs(self.similarity - 1) < 1e-12 && self.offset == 2);
    const auto rc = match_motif(matrix_of("TATC"), 4, known, 4, 0);
    s.require("TATC matches GATA on the reverse strand", rc.motif == 0 && rc.reverse && std::fabs(rc.similarity - 1) < 1e-12);
  }

  // ---- a planted code -------------------------------------------------------------------------------------
  {
    // Random 60-mers; activity = 1.5 · [TTACCCG present] − 1.0 · [GATAAG present] + 0.5 · (GC fraction − ½) · 4
    // + noise of sd 0.2. Families of one sequence each.
    const std::string up = "TTACCCG", down = "GATAAG";
    calx::random::Stream r(31, 0x434f4445, 0, 0);
    Measured m;
    for (int i = 0; i < 6000; i++) {
      std::string x;
      for (int t = 0; t < 60; t++) x += "ACGT"[r.below(4)];
      const bool a = r.bernoulli(0.4), b = r.bernoulli(0.4);
      if (a) x.replace(5 + r.below(20), up.size(), up);
      if (b) x.replace(32 + r.below(20), down.size(), down);
      double gc = 0;
      for (char c : x) gc += c == 'G' || c == 'C';
      m.sequences.push_back(sequence::encode_dna(x));
      m.values.push_back(1.5 * a - 1.0 * b + 2.0 * (gc / 60 - 0.5) + r.normal(0, 0.2));
    }
    const auto parts = grouped_split(m, {{0, 60}}, 0.1, 0.2, 3);
    const auto tr = subset(m, parts.train), va = subset(m, parts.validation), te = subset(m, parts.test);
    KmerBaselineOptions bo;
    bo.k_max = 3;
    bo.gapped_l = 6;
    bo.gapped_k = 4;
    const auto base = fit_kmer_baseline(tr, bo);
    const auto bfit = agreement(predict(base, te.sequences), te.values);
    s.require("the k-mer baseline (k ≤ 3 and gapped 4-mers in 6, elastic net, λ by 5-fold CV) on the planted code: held-out r² " + fixed(bfit.r2) + " (" +
                  std::to_string(base.nonzero) + " of " + std::to_string(base.names.size()) + " coefficients non-zero), above 0.6",
              bfit.r2 > 0.6);
    CnnOptions co;
    co.filters = 8;
    co.width = 9;
    co.pool = 4;
    co.hidden = 16;
    co.epochs = 12;
    co.rate = 3e-3;
    const auto t0 = std::chrono::steady_clock::now();
    const auto cnn = train_cnn(tr, va, co, &pool);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const auto cfit = agreement(predict(cnn, te.sequences, &pool), te.values);
    s.require("the CNN (8 filters of 9, pooling, dense 16) on the planted code: held-out r² " + fixed(cfit.r2) + " in " + fixed(seconds, 1) + " s, above 0.8 (the noise of sd 0.2 caps it near 0.95)",
              cfit.r2 > 0.8);
    // The filters, read from their activating windows, match the planted words.
    std::vector<KnownMotif> planted(2);
    planted[0] = {"P1", up, up.size(), matrix_of(up)};
    planted[1] = {"P2", down, down.size(), matrix_of(down)};
    const auto motifs = filter_motifs(cnn, te.sequences);
    std::set<std::size_t> found;
    for (const auto& f : motifs) {
      if (f.sites < 20) continue;
      const auto [a, b] = informative_core(f.pwm);
      if (b - a < 5) continue;
      std::vector<double> core(f.pwm.probability.begin() + static_cast<std::ptrdiff_t>(4 * a), f.pwm.probability.begin() + static_cast<std::ptrdiff_t>(4 * b));
      const auto mm = match_motif(core, b - a, planted, 5, 0);
      if (mm.similarity > 0.8) found.insert(mm.motif);
    }
    s.require("both planted words are among the filters' motifs (column correlation above 0.8)", found.size() == 2, std::to_string(found.size()) + " found");
  }

  // ---- the DREAM Challenge's measured promoters ----------------------------------------------------------
  const char* env = std::getenv("SCIENCE_CACHE");
  const std::filesystem::path cache = env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
  const auto text = read_file(cache / "dream-promoters" / "filtered_test_data_with_MAUDE_expression.txt");
  if (!text) {
    s.skip("the DREAM Challenge's promoters are not in the data cache (data/cache/dream-promoters): no model is learnt from measured activity");
    return s.done();
  }
  const auto all = read_measured(*text);
  s.require("71,103 measured promoters of 110 bases read", all.problems.empty() && all.size() == 71103 && all.sequences[0].size() == 110);
  // A stated subset: families keyed by the variable region's bases 17–36 and 77–96; 12% of the
  // sequences to train, 2% to validate, 5% to test (about 8,500, 1,400 and 3,600), the rest unused.
  const auto parts = grouped_split(all, {{17, 20}, {77, 20}}, 0.02, 0.05, 7);
  std::vector<std::size_t> tr_rows(parts.train.begin(), parts.train.begin() + static_cast<std::ptrdiff_t>(std::min<std::size_t>(parts.train.size(), 8500)));
  const auto tr = subset(all, tr_rows), va = subset(all, parts.validation), te = subset(all, parts.test);
  KmerBaselineOptions bo;
  bo.from = 17;
  bo.length = 80;
  bo.k_max = 4;
  bo.gapped_k = 0;
  bo.folds = 3;
  bo.lambda_count = 12;
  const auto t1 = std::chrono::steady_clock::now();
  const auto base = fit_kmer_baseline(tr, bo);
  const double base_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
  const auto bfit = agreement(predict(base, te.sequences), te.values);
  CnnOptions co;
  co.epochs = 8;
  const auto t2 = std::chrono::steady_clock::now();
  const auto cnn = train_cnn(tr, va, co, &pool);
  const double cnn_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t2).count();
  const auto cfit = agreement(predict(cnn, te.sequences, &pool), te.values);
  s.require("DREAM promoters, " + std::to_string(tr.size()) + " to train and " + std::to_string(te.size()) + " held-out families' sequences: the k-mer baseline (k ≤ 4, λ by 3-fold CV over 12 values, " + fixed(base_seconds, 1) +
                " s) Pearson r " + fixed(bfit.pearson) + ", Spearman ρ " + fixed(bfit.spearman) + "; above 0.3",
            bfit.pearson > 0.3);
  s.require("the CNN (32 filters of 15, " + std::to_string(cnn.training.curve.size()) + " epochs, " + fixed(cnn_seconds, 1) + " s) Pearson r " + fixed(cfit.pearson) + ", Spearman ρ " +
                fixed(cfit.spearman) + "; above 0.3",
            cfit.pearson > 0.3);
  const auto jaspar = read_file(cache / "jaspar-fungi" / "JASPAR2024_CORE_fungi_non-redundant_pfms_meme.txt");
  if (!jaspar) {
    s.skip("JASPAR's fungal profiles are not in the data cache (data/cache/jaspar-fungi): the filters are not matched");
    return s.done();
  }
  const auto profiles = read_meme_motifs(*jaspar);
  s.require("JASPAR 2024's 178 fungal profiles read", profiles.size() == 178);
  const auto motifs = filter_motifs(cnn, te.sequences);
  std::string named;
  std::size_t matched = 0;
  for (const auto& f : motifs) {
    if (f.sites < 20) continue;
    const auto [a, b] = informative_core(f.pwm);
    if (b - a < 5) continue;
    std::vector<double> core(f.pwm.probability.begin() + static_cast<std::ptrdiff_t>(4 * a), f.pwm.probability.begin() + static_cast<std::ptrdiff_t>(4 * b));
    const auto mm = match_motif(core, b - a, profiles, 5, 100, f.filter + 1);
    if (mm.p_value <= 0.01 && mm.similarity > 0.7) {
      matched++;
      if (named.size() < 160) named += (named.empty() ? "" : ", ") + profiles[mm.motif].name + " (" + fixed(mm.similarity, 2) + ")";
    }
  }
  s.require("filters matching a JASPAR fungal profile (column correlation above 0.7, empirical p ≤ 0.01): " + std::to_string(matched) + " — " + named, matched >= 1);
  return s.done();
}
