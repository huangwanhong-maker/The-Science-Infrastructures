// The toolbox's operations of phylogeny and populations (src/tools/data.js and alignment.js): a
// tree built from the user's sequences, aligned here or already aligned, by the same core as the
// Trees view's gene trees (shared_trees.hpp); a tree in Newick read and rooted; dN/dS of two coding
// sequences, by the same core as trees.dnds; and the samples of a VCF on their principal
// components, with the diversity of their sites.
//
// Every computation is the domains' (genomics::msa, evolution::phylogeny, evolution::selection,
// evolution::population): the engine reads the inputs, holds them within bounds and reports.

#pragma once

#include "genomes.hpp"
#include "inputs.hpp"
#include "shared.hpp"
#include "shared_trees.hpp"

#include <genomics/msa/progressive.hpp>
#include <genomics/variant/vcf.hpp>

#include <evolution/phylogeny/newick.hpp>
#include <evolution/population/diversity.hpp>
#include <evolution/population/genotypes.hpp>
#include <evolution/population/structure.hpp>
#include <evolution/population/vcf.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <string>
#include <vector>

namespace workbench {

namespace pop = evolution::population;

// A tree of DNA sequences: aligned progressively (MUSCLE's three stages) unless they are aligned
// already, then by UPGMA or neighbour joining on distances corrected by the model, or by maximum
// likelihood under it, with the bootstrap; written in Newick with its supports.
BENCH_OP("tool.tree.build") {
  const auto o = detail::tree_options(request);
  const bool aligned = flag(request, "aligned", false);
  detail::Bounds b;
  b.bytes = std::size_t{8} << 20;
  b.records = 300;
  detail::Reading reading;
  const auto records = all_records(request, "fasta", b, reading);
  if (records.size() < 3) throw Fault("a tree needs at least three sequences");
  std::size_t total = 0, longest = 0;
  for (const auto& r : records) {
    require_dna(r, "a tree here");
    total += r.sequence.size();
    longest = std::max(longest, r.sequence.size());
  }
  if (total > 600'000 || longest > 40'000) throw Fault("a tree here takes up to 600,000 bases, none longer than 40,000");
  const auto names = distinct_names(records);
  std::vector<std::string> rows;
  double align_seconds = 0;
  std::string method = "given";
  if (aligned) {
    for (const auto& r : records) {
      std::string row = r.sequence;
      for (char& c : row) c = c == '.' ? '-' : static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
      if (!rows.empty() && row.size() != rows[0].size())
        throw Fault("aligned rows have one length: " + r.id + " has " + std::to_string(row.size()) + " columns, " + records[0].id + " " + std::to_string(rows[0].size()));
      rows.push_back(std::move(row));
    }
  } else {
    namespace msa = genomics::msa;
    std::vector<seq::Codes> seqs;
    for (const auto& r : records) {
      std::string letters;
      for (char c : r.sequence)
        if (c != '-' && c != '.') letters += c;
      seqs.push_back(seq::encode_dna(letters));
    }
    const auto t = std::chrono::steady_clock::now();
    const auto res = msa::progressive(seqs, names, msa::Alphabet::dna, msa::default_options(msa::Alphabet::dna), worker.pool);
    align_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
    for (const auto& r : res.msa.rows) rows.push_back(seq::decode_dna(r));
    method = "progressive (MUSCLE's k-mer guide tree, re-estimated tree and refinement); sum of pairs " + std::to_string(res.sp);
  }
  auto built = detail::tree_of_alignment(names, rows, o, worker.pool);
  Json out = std::move(built.json);
  const Json columns = detail::alignment_columns(rows);
  for (const char* k : {"columns", "variable", "informative"}) out.set(k, columns[k]);
  phy::NewickWriting w;
  w.precision = 6;
  out.set("newick", phy::write_newick(built.tree, w));
  out.set("taxa", names.size());
  out.set("method", o.method);
  out.set("model", o.model);
  Json alignment = Json::object({{"by", method}, {"seconds", align_seconds}});
  if (!aligned && rows.size() * (rows.empty() ? 0 : rows[0].size()) <= 2'000'000) {
    alignment.set("names", names);
    alignment.set("rows", rows);
  }
  out.set("alignment", alignment);
  out.set("seconds", Json::object({{"align", align_seconds}, {"tree", built.seconds}}));
  out.set("reading", reading.json());
  return out;
}

// A tree read from Newick (Felsenstein's grammar; numeric internal labels read as supports), rooted
// as written or at the midpoint of its longest leaf-to-leaf path.
BENCH_OP("tool.tree.read") {
  const std::string text = str(request, "text");
  if (text.find_first_not_of(" \t\r\n") == std::string::npos) throw Fault("give a tree in Newick");
  if (text.size() > (std::size_t{8} << 20)) throw Fault("a tree here has at most 8 MB of Newick");
  const std::string root = str(request, "root", "midpoint");
  if (root != "as" && root != "midpoint") throw Fault("root the tree as written or at the midpoint");
  const auto read = phy::read_newick(text);
  Json problems = Json::array();
  for (std::size_t i = 0; i < read.problems.size() && i < 20; i++) problems.push(Json::object({{"line", read.problems[i].line}, {"message", read.problems[i].message}}));
  if (read.trees.empty())
    throw Fault("no tree could be read" + (read.problems.empty() ? std::string() : ": line " + std::to_string(read.problems[0].line) + ", " + read.problems[0].message));
  phy::Tree t = read.trees.front();
  if (root == "midpoint") {
    if (!t.has_lengths()) throw Fault("a midpoint root needs a length on every branch");
    t = phy::ladderize(phy::midpoint_root(t));
  }
  phy::NewickWriting w;
  w.precision = 8;
  return Json::object({{"leaves", t.leaf_count()}, {"length", t.total_length()}, {"binary", t.is_binary()}, {"rooted", t.rooted}, {"lengths", t.has_lengths()},
                       {"tree", detail::tree_json(t)}, {"newick", phy::write_newick(t, w)}, {"trees", read.trees.size()}, {"root", root},
                       {"problems", problems}, {"problem_count", read.problems.size()}});
}

// dN/dS of two coding sequences: aligned by their codons through their proteins, then synonymous
// and non-synonymous sites and differences by Nei and Gojobori (1986), under a genetic code.
BENCH_OP("tool.align.dnds") {
  const auto& code = genetic_code_of(request, 1);
  const auto letters = [](const std::string& s) {
    std::string out;
    for (char c : s)
      if (std::isalpha(static_cast<unsigned char>(c))) out += c;
    return out;
  };
  const std::string a = letters(str(request, "a")), b = letters(str(request, "b"));
  if (a.size() < 3 || b.size() < 3) throw Fault("give two coding sequences of at least one codon");
  if (a.size() > 200'000 || b.size() > 200'000) throw Fault("a coding sequence here has at most 200,000 bases");
  if (detail::looks_protein(a) || detail::looks_protein(b)) throw Fault("dN/dS is of coding sequences of DNA");
  const evolution::selection::CodonTable codons(code);
  Json out;
  try {
    out = detail::dnds_pair(seq::encode_dna(a), seq::encode_dna(b), code, codons, true);
  } catch (const std::invalid_argument& e) {
    throw Fault(detail::message_of(e));
  }
  out.set("table", code.id);
  out.set("code", std::string(code.name));
  out.set("lengths", Json::array({a.size(), b.size()}));
  return out;
}

// The samples of a VCF: its biallelic sites as dosages (evolution::population, from the genomics
// domain's reader), their principal components (Patterson, Price and Reich 2006) and the diversity
// of the sites: segregating sites, Watterson's θ, π and Tajima's D (1989).
BENCH_OP("tool.pop.vcf") {
  detail::Reading reading;
  const std::string text = vcf_text(request, std::size_t{512} << 20, reading);
  const auto vcf = genomics::variant::read_vcf(text, &reading.problems);
  if (vcf.header.samples.empty()) throw Fault(reading.source + " has no samples");
  pop::VcfConversion conversion;
  const auto g = pop::genotypes_from_vcf(vcf, &conversion);
  reading.records = vcf.records.size();
  if (g.sites == 0) throw Fault(reading.source + " has no biallelic sites");
  const auto counts = pop::allele_counts(g);
  const auto S = pop::segregating_sites(counts);
  const double pi = pop::nucleotide_diversity(counts), theta = pop::theta_watterson(counts), d = pop::tajima_d(counts);
  Json out = Json::object({{"samples", Json(std::vector<std::string>(vcf.header.samples.begin(), vcf.header.samples.begin() + static_cast<std::ptrdiff_t>(std::min<std::size_t>(5000, vcf.header.samples.size()))))},
                           {"sample_count", g.samples}, {"sites", g.sites}, {"chromosomes", g.chromosomes.size()},
                           {"diversity", Json::object({{"segregating", S}, {"pi", detail::finite_or_null(pi)}, {"theta", detail::finite_or_null(theta)}, {"tajima_d", detail::finite_or_null(d)}})}});
  if (g.samples >= 3 && g.samples <= 5000) {
    const auto pc = pop::pca(g, std::min<std::size_t>(4, g.samples), &worker.pool);
    std::vector<double> pc1(g.samples), pc2(g.samples, 0.0);
    for (std::size_t i = 0; i < g.samples; i++) {
      pc1[i] = pc.vectors[i * pc.components];
      if (pc.components > 1) pc2[i] = pc.vectors[i * pc.components + 1];
    }
    std::vector<double> explained;
    for (const double v : pc.values) explained.push_back(v / static_cast<double>(g.samples - 1));
    out.set("pca", Json::object({{"pc1", numbers(pc1, 5)}, {"pc2", numbers(pc2, 5)}, {"values", numbers(pc.values, 4)}, {"explained", numbers(explained, 4)},
                                 {"tracy_widom", numbers(pc.tracy_widom, 3)}, {"sites_used", pc.sites_used}}));
  }
  out.set("conversion", Json::object({{"records", conversion.records}, {"sites", conversion.sites}, {"not_biallelic", conversion.not_biallelic},
                                      {"other_ploidy", conversion.other_ploidy}, {"unknown_allele", conversion.unknown_allele}, {"absent", conversion.absent},
                                      {"phased_records", conversion.phased_records}, {"phased", conversion.phased}}));
  out.set("reading", reading.json());
  return out;
}

}  // namespace workbench
