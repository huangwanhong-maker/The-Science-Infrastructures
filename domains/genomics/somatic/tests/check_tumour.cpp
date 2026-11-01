// A tumour and its normal end to end, against the simulator's truth. A 500 kb reference with
// microsatellites planted in it; a normal diploid (SNVs and indels); a tumour of purity 0.75 with a
// trunk and two nested subclones (CCF 1, 0.55, 0.3), each with its own synthetic signature, an
// allele-specific profile of six states (1+1, 2+1, 1+0, 2+0, 2+2, 3+1) and microsatellites slipped;
// 100× of the tumour and 40× of the normal, 150-base reads, mapped. Checked: the somatic calls'
// precision and recall (by expected allele fraction); the same calls on one thread and eight;
// purity, ploidy and the integer copy numbers by base; the clones' CCFs and the true tree among
// those the sum rule allows; the catalogue's exposures to the three signatures; MSI's sites and
// score; the burden.

#include <genomics/mapping/mapper.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/mutate.hpp>
#include <genomics/simulate/reads.hpp>
#include <genomics/simulate/tumour.hpp>
#include <genomics/somatic/catalogue.hpp>
#include <genomics/somatic/clonality.hpp>
#include <genomics/somatic/copy_number.hpp>
#include <genomics/somatic/msi.hpp>
#include <genomics/somatic/mutect.hpp>
#include <genomics/somatic/signatures.hpp>
#include <genomics/variant/caller.hpp>
#include <genomics/variant/evaluate.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <set>
#include <string>
#include <tuple>
#include <vector>

using namespace genomics;
using sequence::Codes;

namespace {

double since(std::chrono::steady_clock::time_point t0) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool one(1), pool(8);

  // ---- the reference, with microsatellites ------------------------------------------------------------
  const std::size_t L = 500000;
  Codes ref = simulate::random_sequence(L, 0.45, 101);
  {
    calx::random::Stream r(102, 0, 0, 0);
    for (std::size_t at = 2000; at + 200 < L; at += 2500) {
      const std::size_t unit = 1 + r.below(4), copies = unit == 1 ? 10 + r.below(8) : 6 + r.below(8);
      Codes motif(unit);
      do {
        for (auto& b : motif) b = static_cast<sequence::Code>(r.below(4));
      } while (unit > 1 && std::all_of(motif.begin(), motif.end(), [&](sequence::Code b) { return b == motif[0]; }));
      for (std::size_t k = 0; k < unit * copies; k++) ref[at + k] = motif[k % unit];
    }
  }
  const std::vector<Codes> contigs{ref};
  const auto sites = somatic::find_microsatellites(contigs, {10, 50, 5, 40, 5, 5});
  std::vector<simulate::RepeatSite> repeat_sites;
  for (const auto& m : sites) repeat_sites.push_back({m.contig, m.begin, m.unit, m.copies});

  // ---- the normal and the tumour ---------------------------------------------------------------------
  simulate::MutationRates rates;
  rates.snv = 1.5e-3;
  rates.indel = 1e-4;
  const simulate::Diploid normal = simulate::mutate_diploid(ref, rates, 103);
  const auto sig = simulate::synthetic_signatures();
  simulate::TumourSpec spec;
  spec.purity = 0.75;
  spec.clones = {{-1, 1.0, 900, 40, {1, 0, 0}}, {0, 0.55, 450, 15, {0, 1, 0}}, {1, 0.3, 350, 10, {0, 0, 1}}};
  spec.signatures = {sig[0].profile, sig[1].profile, sig[3].profile};
  spec.copy_number = {{0, 0, 90000, {1, 1}}, {0, 90000, 170000, {2, 1}}, {0, 170000, 250000, {1, 0}}, {0, 250000, 330000, {2, 0}}, {0, 330000, 410000, {2, 2}}, {0, 410000, L, {3, 1}}};
  spec.msi.sites = repeat_sites;
  spec.msi.rate = 0.3;
  const simulate::Tumour tumour = simulate::make_tumour(contigs, {normal}, spec, 104);
  std::printf("      truth: %zu somatic variants, %zu microsatellites (%zu slipped), ploidy %.3f\n", tumour.truth.size(), tumour.msi_sites, tumour.slipped.size(), tumour.ploidy());

  simulate::ShortReadModel model;
  auto t0 = std::chrono::steady_clock::now();
  const auto treads = simulate::tumour_reads(tumour, simulate::reads_for_depth(tumour, 100, model.length), model, 105);
  const auto nreads = simulate::short_reads(simulate::normal_source(tumour), simulate::reads_for_depth(tumour, 40, model.length), model, 106);
  const mapping::Mapper mapper(pool, contigs, mapping::MapperOptions::short_reads());
  const auto tmap = mapper.map_all(pool, treads), nmap = mapper.map_all(pool, nreads);
  std::printf("      %zu tumour and %zu normal reads simulated and mapped in %.1f s\n", treads.size(), nreads.size(), since(t0));
  const auto tviews = variant::read_views(treads), nviews = variant::read_views(nreads);
  const somatic::MappedReads T(tmap, tviews), N(nmap, nviews);

  // ---- somatic calls ---------------------------------------------------------------------------------
  t0 = std::chrono::steady_clock::now();
  const auto calls = somatic::call_somatic(pool, contigs, T, N);
  std::printf("      %zu candidates weighed, %zu calls in %.1f s; callable %zu bases\n", calls.sites_weighed, calls.calls.size(), since(t0), calls.callable);
  {
    const auto calls1 = somatic::call_somatic(one, contigs, T, N);
    bool same = calls1.calls.size() == calls.calls.size() && calls1.callable == calls.callable;
    for (std::size_t i = 0; same && i < calls.calls.size(); i++)
      same = calls1.calls[i].position == calls.calls[i].position && calls1.calls[i].alt == calls.calls[i].alt && calls1.calls[i].lod_tumour == calls.calls[i].lod_tumour;
    s.require("the same calls on one thread and on eight", same);
  }
  std::vector<variant::Call> as;
  for (const auto& c : calls.calls) as.push_back(somatic::as_call(c));
  std::vector<variant::TruthVariant> truth;
  for (const auto& v : tumour.truth) truth.push_back({v.contig, v.position, v.ref, v.alt, 1});
  const auto cmp = variant::compare(as, truth, contigs);
  std::printf("      SNVs: TP %zu FP %zu FN %zu (precision %.4f recall %.4f); indels: TP %zu FP %zu FN %zu\n", cmp.snv.tp, cmp.snv.fp, cmp.snv.fn, cmp.snv.precision(),
              cmp.snv.recall(), cmp.insertion.tp + cmp.deletion.tp, cmp.insertion.fp + cmp.deletion.fp, cmp.insertion.fn + cmp.deletion.fn);
  s.check("somatic SNVs: 1 − precision", 1 - cmp.snv.precision(), 0.03);
  s.check("somatic indels (microsatellite slips included): 1 − precision", 1 - static_cast<double>(cmp.insertion.tp + cmp.deletion.tp) / std::max<double>(1, cmp.insertion.tp + cmp.deletion.tp + cmp.insertion.fp + cmp.deletion.fp), 0.1);
  {
    // Recall by expected allele fraction.
    std::set<std::tuple<std::size_t, std::string, std::string>> found;
    for (const auto& c : calls.calls)
      if (c.passed()) found.insert({c.position, c.ref, c.alt});
    std::size_t hi = 0, hi_found = 0, lo = 0, lo_found = 0;
    for (const auto& v : tumour.truth) {
      if (!v.snv()) continue;
      const bool f = found.count({v.position, v.ref, v.alt}) > 0;
      if (v.vaf >= 0.15) { hi++; hi_found += f; } else { lo++; lo_found += f; }
    }
    std::printf("      SNV recall at expected VAF ≥ 0.15: %zu of %zu; below: %zu of %zu\n", hi_found, hi, lo_found, lo);
    s.check("somatic SNVs at expected VAF ≥ 0.15: 1 − recall", 1 - static_cast<double>(hi_found) / static_cast<double>(hi), 0.05);
    s.check("somatic SNVs below 0.15 (subclonal at one copy): 1 − recall", 1 - static_cast<double>(lo_found) / static_cast<double>(lo), 0.5);
  }

  // ---- copy number -----------------------------------------------------------------------------------
  t0 = std::chrono::steady_clock::now();
  const auto ncalls = variant::call_variants(pool, contigs, nmap, nviews);
  const auto profile = somatic::copy_number(pool, contigs, T, N, ncalls);
  std::printf("      copy number in %.1f s: %zu sites, %zu segments; ρ %.2f ψ %.2f (pass %d, goodness %.1f%%), ploidy %.3f\n", since(t0), profile.sites.size(), profile.segments.size(),
              profile.fit.purity, profile.fit.psi, profile.fit.pass, profile.fit.goodness, profile.ploidy);
  s.require("ASCAT finds a fit", profile.fit.found);
  s.check("purity within 0.05 of 0.75", profile.fit.purity - 0.75, 0.05);
  s.check("ploidy (the integers' mean) within 0.1 of the truth", profile.ploidy - tumour.ploidy(), 0.1);
  {
    std::size_t right = 0;
    for (const auto& g : profile.segments)
      for (std::size_t p = g.begin; p < g.end; p += 100) {
        const auto& x = tumour.segment_at(0, p);
        right += 100 * (g.major() == std::max(x.copies[0], x.copies[1]) && g.minor() == std::min(x.copies[0], x.copies[1]));
      }
    std::printf("      segments:");
    for (const auto& g : profile.segments) std::printf(" [%zu, %zu) %d+%d", g.begin, g.end, g.major(), g.minor());
    std::printf("\n");
    s.check("the integer copy numbers right over 95% of the bases (1 − fraction)", 1 - static_cast<double>(right) / static_cast<double>(L), 0.05);
  }

  // ---- clonality -------------------------------------------------------------------------------------
  auto ccf = somatic::cancer_cell_fractions(calls.calls, profile.segments, profile.fit.purity);
  const auto clusters = somatic::cluster_mutations(ccf, 6);
  std::printf("      clones:");
  for (const auto& c : clusters.clones) std::printf(" CCF %.3f (%zu)", c.ccf, c.mutations);
  std::printf("\n");
  s.require("three clones", clusters.clones.size() == 3, std::to_string(clusters.clones.size()));
  if (clusters.clones.size() == 3) {
    s.check("the trunk's CCF", clusters.clones[0].ccf - 1, 0.05);
    s.check("the first subclone's CCF (0.55)", clusters.clones[1].ccf - 0.55, 0.07);
    s.check("the second subclone's CCF (0.3)", clusters.clones[2].ccf - 0.3, 0.08);
    std::vector<double> c;
    for (const auto& x : clusters.clones) c.push_back(x.ccf);
    const auto trees = somatic::clone_trees(c);
    bool chain = false;
    for (const auto& t : trees.trees) chain = chain || (t.parent == std::vector<int>{-1, 0, 1});
    s.require("the true tree (a chain) is among those the sum rule allows", chain, std::to_string(trees.trees.size()) + " trees");
  }
  {
    std::size_t clonal_right = 0, n = 0;
    for (const auto& m : ccf) {
      const auto it = std::find_if(tumour.truth.begin(), tumour.truth.end(), [&](const simulate::SomaticVariant& v) { return v.position == m.position; });
      if (it == tumour.truth.end()) continue;
      n++;
      clonal_right += m.clonal == (it->clone == 0);
    }
    std::printf("      clonal or not, by the 95%% interval: %zu of %zu right\n", clonal_right, n);
    s.check("clonal or subclonal by McGranahan's rule (1 − fraction right)", 1 - static_cast<double>(clonal_right) / static_cast<double>(std::max<std::size_t>(n, 1)), 0.15);
    std::size_t m_right = 0, m_n = 0;
    for (const auto& m : ccf) {
      const auto it = std::find_if(tumour.truth.begin(), tumour.truth.end(), [&](const simulate::SomaticVariant& v) { return v.position == m.position; });
      if (it == tumour.truth.end() || it->clone != 0) continue;
      m_n++;
      m_right += m.multiplicity == it->multiplicity;
    }
    s.check("the trunk's multiplicities (1 − fraction right)", 1 - static_cast<double>(m_right) / static_cast<double>(std::max<std::size_t>(m_n, 1)), 0.1);
  }

  // ---- signatures ------------------------------------------------------------------------------------
  {
    const auto cat = somatic::catalogue(calls.calls, contigs);
    const auto fit = somatic::refit(cat.sbs96, spec.signatures);
    // The truth: the called true substitutions by the signature that made them.
    std::vector<double> want(3, 0);
    std::set<std::size_t> called;
    for (const auto& c : calls.calls)
      if (c.passed() && !c.indel) called.insert(c.position);
    for (const auto& v : tumour.truth)
      if (v.sbs96 >= 0 && called.count(v.position)) want[static_cast<std::size_t>(v.signature)]++;
    double err = 0, total = 0;
    for (std::size_t j = 0; j < 3; j++) { err += std::fabs(fit.exposures[j] - want[j]); total += want[j]; }
    std::printf("      exposures refitted %.0f %.0f %.0f, truth %.0f %.0f %.0f, cosine %.4f\n", fit.exposures[0], fit.exposures[1], fit.exposures[2], want[0], want[1], want[2], fit.cosine);
    s.check("the catalogue's exposures to the three signatures (Σ|ê − e| / Σe)", err / total, 0.12);
    s.check("the reconstruction explains the catalogue (1 − cosine)", 1 - fit.cosine, 0.05);
  }

  // ---- MSI ---------------------------------------------------------------------------------------------
  {
    const auto msi = somatic::msi(pool, sites, T, N);
    std::set<std::size_t> slipped;
    for (const auto& v : tumour.truth)
      if (v.microsatellite) slipped.insert(v.site);
    std::size_t tp = 0, fp = 0, fn = 0;
    for (const auto& x : msi.sites) {
      if (!x.tested) continue;
      const bool truth_slip = slipped.count(x.site) > 0;
      if (x.somatic && truth_slip) tp++;
      else if (x.somatic) fp++;
      else if (truth_slip) fn++;
    }
    std::printf("      MSI: %zu sites tested, %zu somatic, score %.1f%% (%s); TP %zu FP %zu FN %zu\n", msi.tested, msi.somatic, msi.score, msi.high ? "MSI-high" : "stable", tp, fp, fn);
    s.require("most sites tested", msi.tested >= sites.size() * 8 / 10, std::to_string(msi.tested) + " of " + std::to_string(sites.size()));
    s.check("MSI sites: 1 − precision", 1 - static_cast<double>(tp) / std::max<double>(1, tp + fp), 0.05);
    s.check("MSI sites: 1 − recall", 1 - static_cast<double>(tp) / std::max<double>(1, tp + fn), 0.2);
    s.require("MSI-high", msi.high);
  }

  // ---- burden ------------------------------------------------------------------------------------------
  {
    const auto b = somatic::mutational_burden(calls);
    std::printf("      burden: %.0f per megabase (%zu substitutions, %zu indels over %zu callable bases)\n", b.per_megabase, b.substitutions, b.indels, b.callable);
    s.require("the callable territory covers 95% of the reference", b.callable >= L * 95 / 100, std::to_string(b.callable));
  }
  return s.done();
}
