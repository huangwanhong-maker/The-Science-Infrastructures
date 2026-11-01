// Checks of the σ⁷⁰ promoter model: its EM recovers a planted two-box model with its spacer law
// from simulated windows (truth known); learnt from Escherichia coli's transcription start sites
// (Ettwiller et al. 2016, data/registry/ecoli-tss-cappable.json) it finds Harley and Reynolds's
// consensus boxes and the 17-base spacer, separates start-site windows from random genomic ones,
// and ranks the sites' measured strength; the scan's coordinates on both strands.

#include <genomics/regulation/sigma70.hpp>

#include <tyche/descriptive/rank.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <optional>
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

int hamming(const std::string& a, const std::string& b) {
  int d = 0;
  for (std::size_t i = 0; i < a.size(); i++) d += a[i] != b[i];
  return d;
}

// The area under the ROC curve of scores of positives against negatives (Mann–Whitney's U / n₁n₀).
double auc(const std::vector<double>& pos, const std::vector<double>& neg) {
  std::vector<double> all(pos);
  all.insert(all.end(), neg.begin(), neg.end());
  const auto r = tyche::descriptive::ranks(all);
  double sum = 0;
  for (std::size_t i = 0; i < pos.size(); i++) sum += r[i];
  const double n1 = static_cast<double>(pos.size()), n0 = static_cast<double>(neg.size());
  return (sum - n1 * (n1 + 1) / 2) / (n1 * n0);
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool pool(4);

  // ---- recovery of a planted model --------------------------------------------------------------------
  {
    // Boxes drawn from matrices that favour CTGACA and GATAAT at 0.7 a column, spacers 15–19 with
    // weights (1, 2, 6, 2, 1)/12, −10 offsets −13 … −11, 60% of windows with a promoter.
    const std::string b35 = "CTGACA", b10 = "GATAAT";
    const double spacer_law[5] = {1.0 / 12, 2.0 / 12, 6.0 / 12, 2.0 / 12, 1.0 / 12};
    calx::random::Stream r(11, 0x504c4e54, 0, 0);
    std::vector<Codes> windows;
    const auto draw_box = [&](const std::string& box, Codes& x, std::size_t at) {
      for (std::size_t k = 0; k < box.size(); k++) {
        if (r.uniform() < 0.7) x[at + k] = sequence::dna_code(box[k]);
        else x[at + k] = static_cast<Code>(r.below(4));
      }
    };
    for (int i = 0; i < 1500; i++) {
      Codes x(70);
      for (auto& c : x) c = static_cast<Code>(r.below(4));
      if (r.uniform() < 0.6) {
        const int off = -13 + static_cast<int>(r.below(3));
        double u = r.uniform();
        std::size_t sp = 15;
        for (int k = 0; k < 5 && u >= spacer_law[k]; k++) { u -= spacer_law[k]; sp++; }
        const std::size_t a = static_cast<std::size_t>(60 + off);
        draw_box(b10, x, a);
        draw_box(b35, x, a - sp - 6);
      }
      windows.push_back(std::move(x));
    }
    const auto m = learn_sigma70(windows, {}, &pool);
    bool monotone = true;
    for (std::size_t k = 1; k < m.trace.size(); k++) monotone = monotone && m.trace[k] >= m.trace[k - 1] - 1e-9;
    s.require("planted boxes recovered: −10 " + consensus(m.minus10) + " (GATAAT), −35 " + consensus(m.minus35) + " (CTGACA)", consensus(m.minus10) == b10 && consensus(m.minus35) == b35);
    double tv = 0;
    for (int k = 0; k < 5; k++) tv += 0.5 * std::fabs(m.pi_s[static_cast<std::size_t>(k)] - spacer_law[k]);
    s.check("the spacer law recovered: total variation from (1, 2, 6, 2, 1)/12", tv, 0.08);
    s.check("the fraction of windows with a promoter recovered (0.6)", m.gamma - 0.6, 0.08);
    s.require("EM's log posterior never falls (" + std::to_string(m.iterations) + " iterations)", monotone);
  }

  // ---- scanning coordinates -----------------------------------------------------------------------------
  {
    Sigma70 m;
    m.minus10 = annotation::pwm_from_sites({sequence::encode_dna("TATAAT")}, 0.01);
    m.minus35 = annotation::pwm_from_sites({sequence::encode_dna("TTGACA")}, 0.01);
    m.pi_r.assign(m.offsets(), 1.0 / static_cast<double>(m.offsets()));
    m.pi_r[static_cast<std::size_t>(-12 - m.r_min)] = 0.5;
    m.pi_s = {0.1, 0.1, 0.6, 0.1, 0.1};
    calx::random::Stream r(4, 0x53434e4e, 0, 0);
    std::string text;
    for (int i = 0; i < 300; i++) text += "ACGT"[r.below(4)];
    // A promoter at 100 on the forward strand (spacer 17): TTGACA at 100, TATAAT at 123.
    text.replace(100, 6, "TTGACA");
    text.replace(123, 6, "TATAAT");
    const auto hits = scan_promoters(m, sequence::encode_dna(text), 15);
    bool fwd = false;
    for (const auto& h : hits) fwd = fwd || (!h.reverse && h.minus35 == 100 && h.minus10 == 123 && h.spacer == 17 && h.start == 135);
    // The same promoter on the reverse strand: its reverse complement placed at 150.
    std::string rc = sequence::reverse_complement(std::string_view(text.substr(100, 29)));
    std::string t2 = text;
    t2.replace(100, 29, std::string(29, 'A'));
    t2.replace(150, 29, rc);
    const auto hits2 = scan_promoters(m, sequence::encode_dna(t2), 15);
    bool rev = false;
    for (const auto& h : hits2) rev = rev || (h.reverse && h.minus10 == 150 && h.minus35 == 173 && h.spacer == 17 && h.start == 143);
    s.require("a planted promoter is found on the forward strand: −35 at 100, −10 at 123, spacer 17, +1 at 135", fwd);
    s.require("and on the reverse strand, in forward coordinates: −10 at 150, −35 at 173, +1 at 143", rev);
  }

  // ---- Escherichia coli's start sites --------------------------------------------------------------------
  const char* env = std::getenv("SCIENCE_CACHE");
  const std::filesystem::path cache = env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
  const auto gtf = read_file(cache / "ecoli-tss-cappable" / "ettwiller-2016-tss.gtf");
  const auto fasta = read_file(cache / "ecoli-tss-cappable" / "U00096.2.fasta");
  if (!gtf || !fasta) {
    s.skip("Ettwiller et al.'s start sites are not in the data cache (data/cache/ecoli-tss-cappable): the model is not learnt from E. coli");
    return s.done();
  }
  std::vector<sequence::Problem> problems;
  auto sites = read_start_sites(*gtf, &problems);
  const auto records = sequence::read_fasta(*fasta);
  s.require("16,359 start sites read, and the genome U00096.2 (4,639,675 bases)", problems.empty() && sites.size() == 16359 && records.size() == 1 && records[0].sequence.size() == 4639675,
            std::to_string(sites.size()) + " sites");
  const Codes genome = sequence::encode_dna(records[0].sequence);
  // The strongest sites learn the model; every other one of them is held out.
  std::sort(sites.begin(), sites.end(), [](const StartSite& a, const StartSite& b) { return a.strength > b.strength; });
  std::vector<Codes> train, test;
  std::vector<double> test_strength;
  for (std::size_t i = 0; i < 4000 && i < sites.size(); i++) {
    auto w = site_window(genome, sites[i]);
    if (w.empty()) continue;
    if (i % 2 == 0) train.push_back(std::move(w));
    else {
      test.push_back(std::move(w));
      test_strength.push_back(sites[i].strength);
    }
  }
  const auto t0 = std::chrono::steady_clock::now();
  const auto m = learn_sigma70(train, {}, &pool);
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  const std::string c10 = consensus(m.minus10), c35 = consensus(m.minus35);
  char spacer_text[64];
  std::snprintf(spacer_text, sizeof spacer_text, "%.2f %.2f %.2f %.2f %.2f", m.pi_s[0], m.pi_s[1], m.pi_s[2], m.pi_s[3], m.pi_s[4]);
  s.require("learnt from E. coli's 2000 strongest start sites (every other of the top 4000) in " + std::to_string(seconds).substr(0, 4) + " s: −10 box " + c10 + ", within one base of TATAAT (Harley and Reynolds 1987)",
            hamming(c10, "TATAAT") <= 1);
  // The −35 box is the weaker (Harley and Reynolds 1987: its TTG the most conserved part), and a
  // box one base off the hexamer with a spacer one longer describes nearly the same promoters; the
  // check reads the spacer from where the box's TTG lies.
  const auto ttg = c35.find("TTG");
  const long canonical_spacer = static_cast<long>(m.spacer_mode()) - (ttg == std::string::npos ? 99 : static_cast<long>(ttg));
  s.require("−35 box " + c35 + ": it holds TTG, the conserved half of TTGACA, at offset " + std::to_string(static_cast<long>(ttg)), ttg != std::string::npos && ttg <= 1);
  s.require("the commonest spacer, " + std::to_string(m.spacer_mode()) + " from this box (law over 15 … 19: " + std::string(spacer_text) + "), is 17 from the canonical hexamer; the −10 box most often at " +
                std::to_string(m.offset_mode()) + ", a promoter in " + std::to_string(m.gamma).substr(0, 4) + " of the windows",
            canonical_spacer == 17);
  // Held out: start-site windows against random windows of the genome.
  std::vector<double> pos, neg;
  for (const auto& w : test) pos.push_back(promoter_score(m, w));
  calx::random::Stream r(21, 0x4e454757, 0, 0);
  for (int i = 0; i < 4000; i++) {
    StartSite fake;
    fake.position = 100 + r.below(genome.size() - 200);
    fake.reverse = r.bernoulli(0.5);
    const auto w = site_window(genome, fake);
    if (!w.empty()) neg.push_back(promoter_score(m, w));
  }
  const double a = auc(pos, neg);
  s.require("held-out start sites against random windows of the genome: area under the ROC curve " + std::to_string(a).substr(0, 5) + ", above 0.8", a > 0.8);
  const double rho = tyche::descriptive::spearman(pos, test_strength);
  s.require("held-out sites' promoter scores and their measured strength (RRS): Spearman ρ = " + std::to_string(rho).substr(0, 6) +
                " — positive, and weak: the boxes are one of many determinants of a start site's output",
            rho > 0);
  return s.done();
}
