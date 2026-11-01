// Measured abundances of a genome's proteins, and how well an index of codon usage ranks them.
//
// - PaxDb's abundance files (Wang, Herrmann, Simonovic, Szklarczyk and von Mering 2015, Proteomics
//   15, 3163–3168; Huang et al. 2023, Mol. Cell. Proteomics 22, 100640): a header of lines "#key:
//   value" (#name, #score, #coverage, #weights, …), then the columns gene_name, string_external_id,
//   abundance, tab-separated, the abundance in parts per million of the proteome's molecules. A
//   STRING identifier is "<taxon>.<locus>", which for E. coli K-12 (taxon 511145) is the Blattner
//   b-number that RefSeq's NC_000913.3 gives as /locus_tag.
// - Agreement. Each gene is matched by its locus tag (a locus with several CDSs, as infB's three
//   translation starts, is counted once, by its first CDS); genes without a measured abundance, or
//   with an abundance of zero, are left out; Spearman's ρ of the index with the abundance over the
//   rest, and Pearson's r of the index with log₁₀ abundance.
//
// Readers never throw on content.

#pragma once

#include <genomics/codon/genes.hpp>

#include <tyche/descriptive/rank.hpp>
#include <tyche/descriptive/summary.hpp>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace genomics::codon {

struct Abundance {
  std::string name, id, locus;  // locus: the identifier after the taxon
  double ppm = 0;
};

struct AbundanceTable {
  std::string name, coverage, score;
  std::vector<Abundance> rows;
  std::unordered_map<std::string, double> by_locus;
  std::vector<Problem> problems;
};

inline AbundanceTable read_paxdb(std::string_view text) {
  AbundanceTable t;
  sequence::detail::for_each_line(text, [&](std::string_view line, std::size_t number) {
    if (line.empty()) return;
    if (line[0] == '#') {
      const auto colon = line.find(':');
      if (colon == std::string_view::npos) return;
      const std::string key(sequence::detail::trim(line.substr(1, colon - 1)));
      const std::string value(sequence::detail::trim(line.substr(colon + 1)));
      if (key == "name") t.name = value;
      else if (key == "coverage") t.coverage = value;
      else if (key == "score") t.score = value;
      return;
    }
    std::vector<std::string_view> col;
    std::size_t from = 0;
    for (;;) {
      const auto tab = line.find('\t', from);
      col.push_back(line.substr(from, tab == std::string_view::npos ? std::string_view::npos : tab - from));
      if (tab == std::string_view::npos) break;
      from = tab + 1;
    }
    if (col.size() < 3) { t.problems.push_back({number, "a row with " + std::to_string(col.size()) + " columns, not three"}); return; }
    Abundance a;
    a.name = std::string(col[0]);
    a.id = std::string(col[1]);
    const auto dot = a.id.find('.');
    a.locus = dot == std::string::npos ? a.id : a.id.substr(dot + 1);
    const std::string v(sequence::detail::trim(col[2]));
    char* end = nullptr;
    a.ppm = std::strtod(v.c_str(), &end);
    if (v.empty() || *end) { t.problems.push_back({number, "the abundance '" + v + "' is not a number"}); return; }
    t.by_locus.emplace(a.locus, a.ppm);
    t.rows.push_back(std::move(a));
  });
  return t;
}

struct Agreement {
  std::size_t genes = 0, matched = 0;
  double spearman = 0, pearson_log = 0;
  std::vector<std::size_t> gene;        // the genes compared (indices into the gene list)
  std::vector<double> index, abundance;
};

inline Agreement agreement(const std::vector<CodingGene>& genes, const std::vector<double>& index, const AbundanceTable& table) {
  Agreement a;
  a.genes = genes.size();
  std::unordered_set<std::string> seen;
  for (std::size_t i = 0; i < genes.size() && i < index.size(); i++) {
    const auto& tag = genes[i].locus_tag;
    if (tag.empty() || !seen.insert(tag).second) continue;
    const auto it = table.by_locus.find(tag);
    if (it == table.by_locus.end() || !(it->second > 0)) continue;
    a.gene.push_back(i);
    a.index.push_back(index[i]);
    a.abundance.push_back(it->second);
  }
  a.matched = a.gene.size();
  if (a.matched >= 3) {
    a.spearman = tyche::descriptive::spearman(a.index, a.abundance);
    std::vector<double> lg(a.abundance.size());
    for (std::size_t k = 0; k < lg.size(); k++) lg[k] = std::log10(a.abundance[k]);
    a.pearson_log = tyche::descriptive::correlation(a.index, lg);
  }
  return a;
}

}  // namespace genomics::codon
