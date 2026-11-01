// Helpers shared by the codon module's checks: the data cache, files, fixtures.

#pragma once

#include <genomics/codon/genes.hpp>
#include <genomics/sequence/genbank.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace check {

inline std::filesystem::path cache() {
  const char* env = std::getenv("SCIENCE_CACHE");
  return env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
}

inline std::optional<std::string> read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  if (!f) return std::nullopt;
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

inline std::optional<genomics::sequence::GenBank> genome(const std::string& accession) {
  const auto text = read_file(cache() / "ncbi-genomes" / (accession + ".gb"));
  if (!text) return std::nullopt;
  auto records = genomics::sequence::read_genbank(*text);
  if (records.empty()) return std::nullopt;
  return records.front();
}

// The rows of a tab-separated fixture (lines beginning with '#' and the header skipped).
inline std::vector<std::vector<std::string>> fixture(const std::string& name, bool header = true) {
  std::vector<std::vector<std::string>> rows;
  const auto text = read_file(std::filesystem::path(GENOMICS_FIXTURES) / name);
  if (!text) return rows;
  std::istringstream in(*text);
  std::string line;
  bool first = true;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (first && header) { first = false; continue; }
    first = false;
    std::vector<std::string> cols;
    std::size_t from = 0;
    for (;;) {
      const auto tab = line.find('\t', from);
      cols.push_back(line.substr(from, tab == std::string::npos ? std::string::npos : tab - from));
      if (tab == std::string::npos) break;
      from = tab + 1;
    }
    rows.push_back(std::move(cols));
  }
  return rows;
}

// The key the fixtures use for a gene: locus tag, start, end, strand.
inline std::string key(const genomics::codon::CodingGene& g) {
  return g.locus_tag + "|" + std::to_string(g.start) + "|" + std::to_string(g.stop) + "|" + (g.reverse ? "-" : "+");
}

// The codon of letters at an index of NCBI's order.
inline std::string codon_letters(int index) {
  const auto b = genomics::sequence::codon_of(index);
  return {genomics::sequence::dna_letter(b[0]), genomics::sequence::dna_letter(b[1]), genomics::sequence::dna_letter(b[2])};
}

}  // namespace check
