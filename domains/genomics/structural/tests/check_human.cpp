// Checks of the detection on real human sequence: two regions of GRCh38 (NC_000022.11 bases
// 30,000,001–35,000,000, bands 22q12.1–q12.3; NC_000021.9 bases 30,000,001–31,800,000, with two
// runs of N), from the data library's cache (ncbi-grch38-chromosomes;
// skipped when absent). Two simulated individuals (labelled so; their variants are drawn here, not
// observed): one with a balanced reciprocal translocation between the regions, a heterozygous
// deletion of 1.2 Mb, a copy-neutral loss of heterozygosity of 500 kb and smaller variants; one with
// a trisomy of the chromosome-21 region. Both read at 30× of 2 × 150 bases. The translocation is
// called as one balanced event at its two breakpoints; the large deletion is a deletion at copy
// number 1; the loss of heterozygosity is copy-neutral, and the deletion's own runs of homozygosity
// are hemizygous, so the two are told apart by depth and B-allele fractions together; the trisomy is
// an aneuploidy of three copies whose heterozygous sites have their minor allele at 1/3.

#include "scenario.hpp"

#include <genomics/karyotype/bands.hpp>
#include <genomics/structural/describe.hpp>

#include <calx/verification/check.hpp>

#include <filesystem>

int main() {
  using namespace scenario;
  calx::verification::Suite s;
  calx::scheduler::ThreadPool pool(8);
  const auto dir = std::filesystem::path(cache()) / "ncbi-grch38-chromosomes";
  if (!std::filesystem::exists(dir / "NC_000022.11.fa") || !std::filesystem::exists(dir / "NC_000021.9.fa")) {
    s.skip("GRCh38 chromosomes 21 and 22 are not in the data cache (python data/acquisition/fetch.py ncbi-grch38-chromosomes)");
    return s.done();
  }
  const auto region = [&](const char* file, std::size_t begin, std::size_t end) {
    const auto records = sequence::read_fasta(read_file((dir / file).string()));
    const Codes all = sequence::encode_dna(records.at(0).sequence);
    return Codes(all.begin() + static_cast<std::ptrdiff_t>(begin), all.begin() + static_cast<std::ptrdiff_t>(end));
  };
  const std::size_t off22 = 30'000'000, off21 = 30'000'000;
  const std::vector<Codes> contigs{region("NC_000022.11.fa", off22, 35'000'000), region("NC_000021.9.fa", off21, 31'800'000)};
  std::size_t n_bases = 0;
  for (const auto x : contigs[1]) n_bases += x > 3;
  s.require("the chromosome-21 region carries its two runs of N (30 bases)", n_bases == 30, std::to_string(n_bases));

  simulate::PairModel pm;
  pm.fragment.gc_bias = 0.5;
  pm.fragment.gc_optimum = 0.45;
  pm.fragment.gc_width = 0.08;
  // ---- the first individual ------------------------------------------------------------------------
  std::vector<simulate::StructuralVariant> svs;
  const auto add = [&](simulate::SvType t, std::uint32_t c, std::size_t a, std::size_t b, std::array<bool, 2> on) -> simulate::StructuralVariant& {
    simulate::StructuralVariant v;
    v.type = t;
    v.contig = c;
    v.begin = a;
    v.end = b;
    v.on = on;
    svs.push_back(v);
    return svs.back();
  };
  {
    auto& t = add(simulate::SvType::translocation, 0, 2'200'000, 2'200'000, {false, true});
    t.partner = 1;
    t.partner_position = 800'000;
  }
  add(simulate::SvType::deletion, 0, 500'000, 1'700'000, {true, false});
  add(simulate::SvType::loss_of_heterozygosity, 0, 2'400'000, 2'900'000, {true, false});
  add(simulate::SvType::deletion, 1, 300'000, 305'000, {true, true});
  add(simulate::SvType::tandem_duplication, 1, 1'200'000, 1'230'000, {false, true});
  add(simulate::SvType::inversion, 1, 1'500'000, 1'510'000, {true, false});
  const Run a = run(pool, contigs, svs, 30, pm, 51);
  report(a, "GRCh38 22q12 and 21q21 (simulated individual 1)");
  const auto& d = a.detection;
  // The translocation, as one event.
  bool tra = false;
  for (const auto& c : d.calls)
    if (c.kind == structural::SvKind::translocation && c.passed() && c.contig == 0 && c.partner == 1 && std::llabs(static_cast<long long>(c.begin) - 2'200'000) < 50 &&
        std::llabs(static_cast<long long>(c.partner_position) - 800'000) < 50)
      tra = true;
  s.require("the balanced translocation is one event at its two breakpoints, t(21;22)", tra);
  bool del = false;
  for (const auto& c : d.calls)
    if ((c.kind == structural::SvKind::deletion || c.kind == structural::SvKind::copy_loss) && c.contig == 0 && c.begin < 520'000 && c.end > 1'680'000 && c.copies > 0.8 && c.copies < 1.2) del = true;
  s.require("the 1.2 Mb deletion is a deletion of one copy", del);
  bool cnloh = false;
  std::size_t hemizygous = 0;  // bases of the deletion in hemizygous runs (a false heterozygous site may cut a run)
  for (const auto& x : d.allelic) {
    if (x.contig == 0 && x.state == structural::AllelicState::copy_neutral_loh && x.begin < 2'450'000 && x.end > 2'850'000) cnloh = true;
    if (x.contig == 0 && x.state == structural::AllelicState::hemizygous_deletion) {
      const std::size_t lo = std::max<std::size_t>(x.begin, 500'000), hi = std::min<std::size_t>(x.end, 1'700'000);
      if (hi > lo) hemizygous += hi - lo;
    }
  }
  s.require("the loss of heterozygosity is copy-neutral (BAF without heterozygous sites, depth at two copies)", cnloh);
  s.check("and the deletion's runs of homozygosity are hemizygous (depth at one copy) over its length: told apart", 1 - static_cast<double>(hemizygous) / 1.2e6, 0.05);
  s.check("copy-number concordance over mappable bins", 1 - a.copy.fraction(), 0.05);

  // ---- the first individual said in VCF, HGVS and ISCN ----------------------------------------------
  {
    const auto ideogram_file = std::filesystem::path(cache()) / "ncbi-ideogram" / "ideogram_9606_GCF_000001305.16_850_V1.tsv";
    if (!std::filesystem::exists(ideogram_file)) {
      s.skip("the GRCh38 ideogram is not in the data cache (ncbi-ideogram)");
    } else {
      const auto ideogram = karyotype::read_ideogram(read_file(ideogram_file.string()), 850);
      const std::vector<structural::Placement> at{{"NC_000022.11", karyotype::chromosome_index("22"), off22}, {"NC_000021.9", karyotype::chromosome_index("21"), off21}};
      const auto desc = structural::describe_sample(contigs, at, d.calls, d.segments, d.allelic, &ideogram);
      std::printf("# ISCN arr: %s\n# ISCN seq: %s\n#           %s\n", desc.arr.c_str(), desc.seq.short_form.c_str(), desc.seq.detailed.c_str());
      bool hgvs_deletion = false;
      for (std::size_t i = 0; i < d.calls.size(); i++) {
        const auto& c = d.calls[i];
        std::printf("#   %s %s | %s | %s\n", structural::sv_kind_name(c.kind), c.passed() ? "PASS" : c.filters.front().c_str(), desc.calls[i].hgvs.c_str(), desc.calls[i].iscn.c_str());
        const auto& h = desc.calls[i].hgvs;
        if (c.kind == structural::SvKind::deletion && c.contig == 1 && h.rfind("NC_000021.9:g.30300", 0) == 0 && h.size() > 3 && h.compare(h.size() - 3, 3, "del") == 0) hgvs_deletion = true;
      }
      s.require("the 5 kb deletion is NC_000021.9:g.30300…_30305…del in HGVS", hgvs_deletion);
      s.require("ISCN arr gives chromosome 22's one-copy loss in its bands and its copy-neutral run of homozygosity (hmz)",
                desc.arr.find("22q12") != std::string::npos && desc.arr.find("x1") != std::string::npos && desc.arr.find("hmz") != std::string::npos, desc.arr);
      std::string tra;
      for (std::size_t i = 0; i < d.calls.size(); i++)
        if (d.calls[i].kind == structural::SvKind::translocation) tra = desc.calls[i].iscn;
      s.require("ISCN seq gives the reciprocal translocation as t(21;22) with its bands", tra.find("t(21;22)(q") != std::string::npos, tra);
    }
  }

  // ---- the second individual: a trisomy ------------------------------------------------------------
  std::vector<simulate::StructuralVariant> tri(1);
  tri[0].type = simulate::SvType::aneuploidy;
  tri[0].contig = 1;
  tri[0].begin = 0;
  tri[0].end = contigs[1].size();
  tri[0].copies = 1;
  tri[0].on = {true, false};
  const Run b = run(pool, contigs, tri, 30, pm, 61);
  report(b, "GRCh38 (simulated individual 2, trisomy of the chromosome-21 region)");
  const auto& an = b.detection.aneuploidy;
  s.require("the trisomy is an aneuploidy of the chromosome-21 region at three copies, and of nothing else", an.size() == 1 && an[0].contig == 1 && an[0].copies == 3);
  bool imbalance = false;
  for (const auto& x : b.detection.allelic)
    if (x.contig == 1 && x.copies == 3 && x.minor == 1) imbalance = true;
  s.require("its heterozygous sites carry the minor allele once in three (2 + 1)", imbalance);
  return s.done();
}
