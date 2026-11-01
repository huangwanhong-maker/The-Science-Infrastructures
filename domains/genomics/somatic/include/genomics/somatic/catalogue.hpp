// Mutational catalogues: the somatic substitutions of a tumour counted by class (SBS-96) and its
// doublets (DBS-78), from calls and the reference (sequence/contexts.hpp defines the classes);
// tables of signatures or catalogues read from text by the names of their classes.
//
// - SBS-96. Each passing substitution of one base, its class from the reference's bases before
//   and after it, pyrimidine-centred (Alexandrov et al. 2013, Nature 500, 415–421).
// - DBS-78. Two passing substitutions on adjacent bases of one contig are one doublet (as
//   SigProfilerMatrixGenerator merges them; Bergstrom et al. 2019, BMC Genomics 20, 685, §"Doublet-
//   base substitutions"), counted in DBS-78 and not in SBS-96; three or more adjacent substitutions
//   (a multi-base substitution) are counted in neither. *Caveat.* Two adjacent calls may lie on
//   different haplotypes; read-backed phasing would tell, and is not done.
// - Tables (read_profiles): COSMIC's and SigProfiler's text matrices, tab-separated, a header of
//   column names and one row per class named in its first column ("A[C>A]A" or "AC>CA"), any
//   row order. The reader never throws on content: rows it cannot place, duplicate rows and
//   columns of the wrong length become problems with their line numbers. COSMIC's own signature
//   tables are not fetched by the data library (their terms forbid redistribution; README §3):
//   a user holding them passes the path.
//
// Complexity: O(c) for c calls; O(rows · columns) for a table.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/contexts.hpp>
#include <genomics/somatic/mutect.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::somatic {

struct Catalogue {
  std::vector<double> sbs96 = std::vector<double>(sequence::SBS96, 0);
  std::vector<double> dbs78 = std::vector<double>(sequence::DBS78, 0);
  std::size_t substitutions = 0, doublets = 0, multiple = 0, unclassified = 0;
  // For each substitution counted in sbs96, its call and class; for each doublet, its first call and class.
  std::vector<std::pair<std::size_t, int>> sbs, dbs;
};

// The catalogue of the passing (or all) substitutions of `calls` against the reference.
inline Catalogue catalogue(const std::vector<SomaticCall>& calls, const std::vector<Codes>& contigs, bool pass_only = true) {
  Catalogue c;
  std::vector<std::size_t> snv;
  for (std::size_t i = 0; i < calls.size(); i++)
    if (!calls[i].indel && calls[i].ref.size() == 1 && calls[i].alt.size() == 1 && (!pass_only || calls[i].passed())) snv.push_back(i);
  std::sort(snv.begin(), snv.end(), [&](std::size_t a, std::size_t b) {
    return calls[a].contig != calls[b].contig ? calls[a].contig < calls[b].contig : calls[a].position < calls[b].position;
  });
  for (std::size_t a = 0; a < snv.size();) {
    std::size_t b = a + 1;
    while (b < snv.size() && calls[snv[b]].contig == calls[snv[a]].contig && calls[snv[b]].position == calls[snv[b - 1]].position + 1) b++;
    const std::size_t run = b - a;
    if (run == 1) {
      const SomaticCall& s = calls[snv[a]];
      const int k = s.contig < contigs.size() ? sequence::sbs96_at(contigs[s.contig], s.position, sequence::dna_code(s.alt[0])) : -1;
      if (k >= 0) {
        c.sbs96[static_cast<std::size_t>(k)]++;
        c.substitutions++;
        c.sbs.push_back({snv[a], k});
      } else {
        c.unclassified++;
      }
    } else if (run == 2) {
      const SomaticCall &x = calls[snv[a]], &y = calls[snv[a + 1]];
      const int k = sequence::dbs78_index(sequence::dna_code(x.ref[0]), sequence::dna_code(y.ref[0]), sequence::dna_code(x.alt[0]), sequence::dna_code(y.alt[0]));
      if (k >= 0) {
        c.dbs78[static_cast<std::size_t>(k)]++;
        c.doublets++;
        c.dbs.push_back({snv[a], k});
      } else {
        c.unclassified += 2;
      }
    } else {
      c.multiple += run;
    }
    a = b;
  }
  return c;
}

struct ProfileTable {
  std::vector<std::string> names;                // the columns after the first
  std::vector<std::vector<double>> profiles;     // one per column, over the classes in this file's order (SBS-96's or DBS-78's)
  std::size_t classes = 0;                       // 96 or 78
  std::vector<std::pair<std::size_t, std::string>> problems;  // line, what
};

namespace detail {

inline int class_of_name(std::string_view n) {
  if (n.size() == 7 && n[1] == '[' && n[3] == '>' && n[5] == ']') {
    using sequence::dna_code;
    const sequence::Code l = dna_code(n[0]), r = dna_code(n[2]), a = dna_code(n[4]), x = dna_code(n[6]);
    if (r != sequence::C && r != sequence::T) return -1;
    return sequence::sbs96_index(l, r, a, x);
  }
  if (n.size() == 5 && n[2] == '>') {
    for (std::size_t k = 0; k < sequence::DBS78; k++)
      if (sequence::DBS78_NAMES[k] == n) return 1000 + static_cast<int>(k);
  }
  return -1;
}

}  // namespace detail

// A table of profiles (signatures or catalogues) over SBS-96 or DBS-78 classes.
inline ProfileTable read_profiles(std::string_view text) {
  ProfileTable t;
  std::istringstream in{std::string(text)};
  std::string line;
  std::size_t number = 0;
  bool header = false;
  std::vector<char> seen;
  while (std::getline(in, line)) {
    number++;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    std::vector<std::string> f;
    std::string cell;
    std::istringstream row(line);
    while (std::getline(row, cell, '\t')) f.push_back(cell);
    if (!header) {
      header = true;
      if (f.size() < 2) { t.problems.push_back({number, "a header of at least two columns"}); return t; }
      t.names.assign(f.begin() + 1, f.end());
      continue;
    }
    const int k = detail::class_of_name(f[0]);
    if (k < 0) { t.problems.push_back({number, "an unknown class " + f[0]}); continue; }
    const std::size_t classes = k >= 1000 ? sequence::DBS78 : sequence::SBS96, idx = static_cast<std::size_t>(k >= 1000 ? k - 1000 : k);
    if (t.classes == 0) {
      t.classes = classes;
      t.profiles.assign(t.names.size(), std::vector<double>(classes, 0));
      seen.assign(classes, 0);
    }
    if (classes != t.classes) { t.problems.push_back({number, "a class of another catalogue: " + f[0]}); continue; }
    if (seen[idx]) { t.problems.push_back({number, "a class given twice: " + f[0]}); continue; }
    if (f.size() != t.names.size() + 1) { t.problems.push_back({number, "a row of " + std::to_string(f.size()) + " cells for " + std::to_string(t.names.size() + 1) + " columns"}); continue; }
    seen[idx] = 1;
    for (std::size_t j = 0; j < t.names.size(); j++) {
      char* end = nullptr;
      const double v = std::strtod(f[j + 1].c_str(), &end);
      if (end == f[j + 1].c_str() || !std::isfinite(v)) {
        t.problems.push_back({number, "not a number in column " + t.names[j]});
        continue;
      }
      t.profiles[j][idx] = v;
    }
  }
  for (std::size_t i = 0; i < seen.size(); i++)
    if (!seen[i]) t.problems.push_back({0, "a class missing: " + (t.classes == sequence::SBS96 ? sequence::sbs96_name(i) : sequence::dbs78_name(i))});
  return t;
}

}  // namespace genomics::somatic
