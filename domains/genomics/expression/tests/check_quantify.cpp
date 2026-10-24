// Checks of quantification: effective lengths against their closed forms (a point mass at f:
// l − f + 1; any law: Σ_{f ≤ l} p(f)(l − f + 1) / P(F ≤ l)); the EM's fixed point against the
// maximum of the likelihood found by a fine search on two transcripts sharing a class (and the
// closed form n₀/(n₀ + n₁) at equal lengths); the log-likelihood never decreasing from round to
// round; on a simulated transcriptome of genes with isoforms sharing exons and paralogs,
// abundances recovered from 400 000 simulated reads (correlation, median and mean relative error
// against the true counts, reported), better than splitting each class's reads equally, overall
// and on genes of several isoforms; accurate where a transcript has reads of its own and on the
// genes' totals (an isoform whose every k-mer lies in other isoforms is not identifiable from
// classes, theory note §5.3); TPM summing to 10⁶, and the genes' TPM against the truth; the
// bootstrap the same on 1 and 8 threads, its ±2 deviations covering the truth on identifiable
// transcripts about 95% of the time, and its spread wider where reads are shared.

#include <genomics/expression/pseudoalign.hpp>
#include <genomics/expression/quantify.hpp>
#include <genomics/expression/transcriptome.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics::expression;

namespace {

double pearson(const std::vector<double>& a, const std::vector<double>& b) {
  const double n = static_cast<double>(a.size());
  double ma = 0, mb = 0;
  for (std::size_t i = 0; i < a.size(); i++) { ma += a[i] / n; mb += b[i] / n; }
  double sab = 0, saa = 0, sbb = 0;
  for (std::size_t i = 0; i < a.size(); i++) {
    sab += (a[i] - ma) * (b[i] - mb);
    saa += (a[i] - ma) * (a[i] - ma);
    sbb += (b[i] - mb) * (b[i] - mb);
  }
  return sab / std::sqrt(saa * sbb);
}

// The mean of |estimate − truth| / truth over the transcripts with truth ≥ 10 (and in `only` if given).
double mean_relative_error(const std::vector<double>& est, const std::vector<double>& truth, const std::vector<bool>* only = nullptr) {
  double s = 0;
  std::size_t n = 0;
  for (std::size_t t = 0; t < truth.size(); t++)
    if (truth[t] >= 10 && (!only || (*only)[t])) { s += std::fabs(est[t] - truth[t]) / truth[t]; n++; }
  return n ? s / static_cast<double>(n) : 0;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool one(1), eight(8);

  // ---- effective lengths ------------------------------------------------------------------------------------
  {
    FragmentLengths point;
    point.p.assign(201, 0.0);
    point.p[200] = 1.0;
    const auto e = effective_lengths({150, 200, 201, 1000}, point);
    s.require("a point mass at 200 gives l − 199 where a fragment fits, and l where none does", e[0] == 150 && e[1] == 1 && e[2] == 2 && e[3] == 801);
    const auto d = FragmentLengths::normal(180, 40);
    bool sums = true;
    for (std::size_t l : {50u, 120u, 180u, 260u, 900u}) {
      double num = 0, mass = 0;
      for (std::size_t f = 1; f <= std::min<std::size_t>(l, d.max()); f++) {
        num += d.p[f] * static_cast<double>(l - f + 1);
        mass += d.p[f];
      }
      const double want = num / mass >= 1 ? num / mass : static_cast<double>(l);
      sums = sums && std::fabs(effective_lengths({l}, d)[0] - want) < 1e-9 * want;
    }
    s.require("for a normal law of mean 180 and sd 40, l − μ_l + 1 equals Σ_{f ≤ l} p(f)(l − f + 1) / P(F ≤ l) at five lengths", sums);
  }

  // ---- the EM's fixed point ------------------------------------------------------------------------------------
  {
    ClassTable classes;
    const std::uint32_t a[] = {0}, b[] = {1}, ab[] = {0, 1};
    classes.insert(a);
    classes.insert(b);
    classes.insert(ab);
    const std::vector<std::uint32_t> n = {120, 40, 300};
    EmOptions tight;
    tight.change = 1e-12;
    tight.change_limit = 0;
    tight.max_rounds = 100000;
    const auto eq = quantify(classes, n, {1000, 1000}, tight, one);
    s.check("at equal lengths the EM's estimate of transcript 0 is N n0/(n0 + n1) = 345", eq.counts[0] - 460.0 * 120 / 160, 1e-6, "fragments");
    // Unequal lengths: the maximum of log L over α₀ by golden-section search.
    const std::vector<double> len = {800, 2000};
    const auto f = [&](double x) {
      return 120 * std::log(x / len[0]) + 40 * std::log((1 - x) / len[1]) + 300 * std::log(x / len[0] + (1 - x) / len[1]);
    };
    double lo = 1e-9, hi = 1 - 1e-9;
    const double g = (std::sqrt(5.0) - 1) / 2;
    for (int i = 0; i < 200; i++) {
      const double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
      if (f(x1) < f(x2)) lo = x1; else hi = x2;
    }
    const auto un = quantify(classes, n, len, tight, one);
    s.check("at lengths 800 and 2000 the EM's fixed point is the maximum of log L (golden-section search): |α̂₀ − α*₀|", un.counts[0] / 460 - 0.5 * (lo + hi), 1e-7);
  }

  // ---- a simulated sample ------------------------------------------------------------------------------------
  {
    TranscriptomeSpec spec;
    spec.genes = 400;
    spec.paralogs = 0.15;
    spec.paralog_divergence = 0.03;
    const auto sim = simulate_transcriptome(spec, 101);
    const Transcriptome& tx = sim.transcriptome;
    const auto rho = simulate_abundances(tx.size(), 1.5, 0.1, 103);
    const auto d = FragmentLengths::normal(200, 30);
    ReadModel model;
    model.length = 75;
    model.error = 0.002;
    const std::size_t N = 400000;
    const auto frags = simulate_fragments(tx, rho, d, model, N, 107);
    const Index index(tx, 31);
    const Pseudoalignment p = pseudoalign(index, frags, eight);
    const auto eff = effective_lengths(tx.lengths(), d);
    EmOptions o;
    o.trace = true;
    const Abundances em = quantify(p, eff, o, eight);
    // The truth: the fragments drawn from each transcript, among those pseudoaligned, so that both
    // sides count the same fragments.
    std::vector<double> truth(tx.size(), 0.0);
    for (std::size_t i = 0; i < frags.size(); i++)
      if (p.fragment_class[i] != ClassTable::NONE) truth[frags[i].transcript] += 1;
    const auto split = equal_split(p.classes, p.counts, tx.size());
    // Transcripts of genes with several isoforms.
    std::vector<std::size_t> per_gene(tx.genes.size(), 0);
    for (auto g : tx.gene_of) per_gene[g]++;
    std::vector<bool> isoform(tx.size());
    for (std::size_t t = 0; t < tx.size(); t++) isoform[t] = per_gene[tx.gene_of[t]] > 1;
    const double r_em = pearson(em.counts, truth), r_split = pearson(split, truth);
    const double e_em = mean_relative_error(em.counts, truth), e_split = mean_relative_error(split, truth);
    const double i_em = mean_relative_error(em.counts, truth, &isoform), i_split = mean_relative_error(split, truth, &isoform);
    // The share of each transcript's reads in classes of its own.
    std::vector<double> own(tx.size(), 0.0), total(tx.size(), 0.0);
    for (std::uint32_t c = 0; c < p.counts.size(); c++)
      for (std::uint32_t t : p.classes.members(c)) {
        total[t] += p.counts[c];
        if (p.classes.members(c).size() == 1) own[t] += p.counts[c];
      }
    std::vector<bool> identifiable(tx.size());
    for (std::size_t t = 0; t < tx.size(); t++) identifiable[t] = total[t] > 0 && own[t] >= 0.1 * total[t];
    const double u_em = mean_relative_error(em.counts, truth, &identifiable);
    std::vector<double> rel;
    for (std::size_t t = 0; t < tx.size(); t++)
      if (truth[t] >= 10) rel.push_back(std::fabs(em.counts[t] - truth[t]) / truth[t]);
    std::sort(rel.begin(), rel.end());
    const double median_rel = rel[rel.size() / 2];
    // Genes: the sums over their isoforms.
    std::vector<double> gene_truth(tx.genes.size(), 0.0), gene_em(tx.genes.size(), 0.0), gene_split(tx.genes.size(), 0.0);
    for (std::size_t t = 0; t < tx.size(); t++) {
      gene_truth[tx.gene_of[t]] += truth[t];
      gene_em[tx.gene_of[t]] += em.counts[t];
    }
    const double g_em = mean_relative_error(gene_em, gene_truth);
    std::printf("%zu transcripts of %zu genes, %zu fragments, %zu pseudoaligned, %zu classes; EM: %d rounds%s\n", tx.size(), tx.genes.size(), N, p.aligned, p.classes.size(),
                em.rounds, em.converged ? ", converged" : "");
    std::printf("  correlation with the true counts: EM %.4f, equal split %.4f\n", r_em, r_split);
    std::printf("  mean relative error (true count >= 10): EM %.4f (median %.4f), equal split %.4f; on genes of several isoforms: EM %.4f, equal split %.4f\n", e_em,
                median_rel, e_split, i_em, i_split);
    std::printf("  on transcripts with 10%% or more of their reads their own: EM %.4f; on gene totals: EM %.4f\n", u_em, g_em);
    s.check("the EM's counts correlate with the true counts: 1 − r", 1 - r_em, 0.01);
    s.require("the EM resolves shared reads better than splitting them equally: lower mean relative error, overall and on genes of several isoforms",
              e_em < e_split && i_em < i_split);
    s.check("the EM's median relative error over transcripts of 10 or more fragments", median_rel, 0.10);
    s.check("its mean relative error on the transcripts with 10% or more of their reads in classes of their own", u_em, 0.06);
    s.check("and on the genes' totals", g_em, 0.01);
    bool monotone = true;
    for (std::size_t k = 1; k < em.log_likelihood.size(); k++) monotone = monotone && em.log_likelihood[k] >= em.log_likelihood[k - 1] - 1e-9 * std::fabs(em.log_likelihood[k - 1]);
    s.require("the log-likelihood never decreases from one round of the EM to the next", monotone && em.log_likelihood.size() == static_cast<std::size_t>(em.rounds),
              std::to_string(em.log_likelihood.front()) + " to " + std::to_string(em.log_likelihood.back()));
    double tpm_sum = 0;
    for (double x : em.tpm) tpm_sum += x;
    s.check("TPM sums to 10^6", tpm_sum / 1e6 - 1, 1e-12);
    // TPM of genes, estimated against simulated.
    std::vector<double> gt(tx.genes.size(), 0.0), ge(tx.genes.size(), 0.0);
    for (std::size_t t = 0; t < tx.size(); t++) {
      gt[tx.gene_of[t]] += 1e6 * rho[t];
      ge[tx.gene_of[t]] += em.tpm[t];
    }
    for (auto* v : {&gt, &ge})
      for (double& x : *v) x = std::log1p(x);
    s.check("log(1 + TPM) of the genes, estimated against simulated: 1 − r", 1 - pearson(ge, gt), 0.01);
    const Abundances em1 = quantify(p, eff, o, one);
    s.require("the EM gives the same counts on 1 and 8 threads", em1.counts == em.counts && em1.rounds == em.rounds);

    // The bootstrap.
    EmOptions bo;
    const Bootstrap b1 = bootstrap(p.classes, p.counts, eff, 30, 109, bo, one), b8 = bootstrap(p.classes, p.counts, eff, 30, 109, bo, eight);
    s.require("the bootstrap (30 replicates) is the same on 1 and 8 threads", b1.mean == b8.mean && b1.sd == b8.sd);
    std::size_t covered = 0, counted = 0;
    for (std::size_t t = 0; t < tx.size(); t++) {
      if (truth[t] < 10 || !identifiable[t]) continue;
      counted++;
      covered += std::fabs(em.counts[t] - truth[t]) <= 2 * b1.sd[t];
    }
    const double coverage = static_cast<double>(covered) / static_cast<double>(counted);
    s.require("on transcripts of 10% or more reads their own, the estimate lies within 2 bootstrap deviations of the truth for 90% of them or more (95% expected)",
              coverage >= 0.90, std::to_string(covered) + " of " + std::to_string(counted));
    double cv_unique = 0, cv_shared = 0;
    std::size_t nu = 0, ns = 0;
    for (std::size_t t = 0; t < tx.size(); t++) {
      if (em.counts[t] < 100) continue;
      const double cv = b1.sd[t] / em.counts[t], f = own[t] / total[t];
      if (f > 0.9) { cv_unique += cv; nu++; }
      else if (f < 0.2) { cv_shared += cv; ns++; }
    }
    cv_unique /= static_cast<double>(std::max<std::size_t>(nu, 1));
    cv_shared /= static_cast<double>(std::max<std::size_t>(ns, 1));
    s.require("the bootstrap's relative spread is wider for transcripts whose reads are mostly shared than for those whose reads are their own",
              nu > 0 && ns > 0 && cv_shared > cv_unique, "mean CV " + std::to_string(cv_shared) + " (" + std::to_string(ns) + " shared) against " + std::to_string(cv_unique) + " (" +
                                                        std::to_string(nu) + " unique)");
  }
  return s.done();
}
