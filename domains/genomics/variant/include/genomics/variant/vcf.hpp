// VCF 4.2, the text form of variant calls (The Variant Call Format Specification, VCFv4.2,
// samtools/hts-specs; Danecek et al. 2011, Bioinformatics 27, 2156–2158), as records written and
// read back; calls become records in caller.hpp.
//
// - Header (§1.2–1.4): ##fileformat=VCFv4.2, then meta-information lines (##key=value, or
//   ##key=<ID=…,…> for INFO, FILTER, FORMAT and contig), then the header line #CHROM POS ID REF
//   ALT QUAL FILTER INFO, followed by FORMAT and the sample names when there are samples.
// - Records (§1.6): CHROM; POS, 1-based, of REF's first base (0-based inside); ID ('.'); REF; ALT,
//   alleles separated by commas or '.'; QUAL, Phred, or '.'; FILTER, PASS, codes separated by ';',
//   or '.'; INFO, key=value or flag entries separated by ';', or '.'; FORMAT, keys separated by
//   ':'; one column per sample of values in FORMAT's order. The caller (caller.hpp, vcf_record
//   and caller_header) writes INFO DP (and FS, the strand bias) and FORMAT GT:GQ:DP:AD, AD with
//   Number=R (one depth per allele, REF's first), as VCF 4.2 allows.
// - Genotypes (§1.6.2, GT): allele indices separated by '/' (unphased) or '|' (phased), '.' for a
//   missing allele: 0/1, 1|0, ./., 1/2, and a haploid 0. A genotype is phased here when every
//   separator is '|'; a haploid call, with none, is not.
// - Reading never throws on content: a line that is not a record (fewer than eight columns, a
//   POS that is not a positive integer, an empty REF, a QUAL that is not a number) is a problem
//   with its line number and is skipped; sample columns not matching the header, a GT that is not
//   one, a meta line before ##fileformat are problems and the record is kept.
//
// Complexity: linear in the text.

#pragma once

#include <genomics/sequence/fasta.hpp>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace genomics::variant {

struct VcfGenotype {
  std::vector<int> alleles;  // −1 for a missing allele
  bool phased = false;
  bool missing() const {
    for (int a : alleles)
      if (a >= 0) return false;
    return true;
  }
  friend bool operator==(const VcfGenotype&, const VcfGenotype&) = default;
};

struct VcfRecord {
  std::string chrom;
  std::size_t position = 0;  // 0-based
  std::string id = ".";
  std::string ref;
  std::vector<std::string> alts;  // empty for '.'
  double qual = std::numeric_limits<double>::quiet_NaN();  // NaN for '.'
  std::vector<std::string> filters;                        // empty for '.'
  std::vector<std::pair<std::string, std::string>> info;   // a flag has the value ""
  std::vector<std::string> format;
  std::vector<std::vector<std::string>> samples;           // values in FORMAT's order
  std::size_t line = 0;

  bool has_qual() const { return !std::isnan(qual); }
  bool passed() const { return filters.size() == 1 && filters[0] == "PASS"; }
  std::string_view info_value(std::string_view key) const {
    for (const auto& [k, v] : info)
      if (k == key) return v;
    return {};
  }
  // A sample's value for a FORMAT key, or "" when absent (VCF's trailing fields may be dropped).
  std::string_view sample_value(std::size_t sample, std::string_view key) const;
  VcfGenotype genotype(std::size_t sample) const;
};

struct VcfHeader {
  std::string fileformat = "VCFv4.2";
  std::vector<std::string> meta;                             // the ## lines after ##fileformat, without "##"
  std::vector<std::pair<std::string, std::size_t>> contigs;  // from ##contig=<ID=…,length=…>
  std::vector<std::string> samples;
};

struct Vcf {
  VcfHeader header;
  std::vector<VcfRecord> records;
};

namespace detail {

template <class F>
void vcf_lines(std::string_view text, F&& f) {
  std::size_t number = 0, i = 0;
  while (i < text.size()) {
    std::size_t j = text.find('\n', i);
    if (j == std::string_view::npos) j = text.size();
    std::string_view line = text.substr(i, j - i);
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    f(line, ++number);
    i = j + 1;
  }
}

inline std::vector<std::string_view> fields(std::string_view s, char sep) {
  std::vector<std::string_view> f;
  std::size_t i = 0;
  for (;;) {
    const std::size_t j = s.find(sep, i);
    f.push_back(s.substr(i, j == std::string_view::npos ? std::string_view::npos : j - i));
    if (j == std::string_view::npos) break;
    i = j + 1;
  }
  return f;
}

template <class T>
bool parse_number(std::string_view s, T& out) {
  if (s.empty()) return false;
  const auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), out);
  return ec == std::errc() && p == s.data() + s.size();
}

// A Phred value with at most two decimals, without trailing zeros.
inline std::string qual_string(double q) {
  if (std::isnan(q)) return ".";
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.2f", q);
  std::string s(buf);
  while (!s.empty() && s.back() == '0') s.pop_back();
  if (!s.empty() && s.back() == '.') s.pop_back();
  return s;
}

}  // namespace detail

// A GT value; `ok` false when it is not one.
inline VcfGenotype parse_genotype(std::string_view gt, bool* ok = nullptr) {
  VcfGenotype g;
  bool good = !gt.empty(), any_separator = false, all_phased = true;
  std::size_t i = 0;
  while (i < gt.size() && good) {
    if (gt[i] == '.') {
      g.alleles.push_back(-1);
      i++;
    } else if (gt[i] >= '0' && gt[i] <= '9') {
      int v = 0;
      while (i < gt.size() && gt[i] >= '0' && gt[i] <= '9') v = v * 10 + (gt[i++] - '0');
      g.alleles.push_back(v);
    } else {
      good = false;
      break;
    }
    if (i < gt.size()) {
      if (gt[i] == '/' || gt[i] == '|') {
        any_separator = true;
        all_phased = all_phased && gt[i] == '|';
        i++;
        if (i == gt.size()) good = false;
      } else {
        good = false;
      }
    }
  }
  g.phased = good && any_separator && all_phased;
  if (ok) *ok = good;
  return g;
}

inline std::string genotype_string(const VcfGenotype& g) {
  if (g.alleles.empty()) return ".";
  std::string s;
  for (std::size_t i = 0; i < g.alleles.size(); i++) {
    if (i) s += g.phased ? '|' : '/';
    s += g.alleles[i] < 0 ? "." : std::to_string(g.alleles[i]);
  }
  return s;
}

inline std::string_view VcfRecord::sample_value(std::size_t sample, std::string_view key) const {
  if (sample >= samples.size()) return {};
  for (std::size_t k = 0; k < format.size(); k++)
    if (format[k] == key) return k < samples[sample].size() ? std::string_view(samples[sample][k]) : std::string_view{};
  return {};
}

inline VcfGenotype VcfRecord::genotype(std::size_t sample) const {
  const auto gt = sample_value(sample, "GT");
  if (gt.empty()) return {};
  return parse_genotype(gt);
}

// ---- reading --------------------------------------------------------------------------------------
inline Vcf read_vcf(std::string_view text, std::vector<sequence::Problem>* problems = nullptr) {
  Vcf out;
  bool header_line = false, fileformat = false;
  const auto problem = [&](std::size_t line, std::string m) { if (problems) problems->push_back({line, std::move(m)}); };
  detail::vcf_lines(text, [&](std::string_view line, std::size_t number) {
    if (line.empty()) return;
    if (line.substr(0, 2) == "##") {
      const std::string_view body = line.substr(2);
      if (body.substr(0, 11) == "fileformat=") {
        out.header.fileformat = std::string(body.substr(11));
        fileformat = true;
        return;
      }
      if (!fileformat) problem(number, "a meta-information line before ##fileformat");
      out.header.meta.emplace_back(body);
      if (body.substr(0, 8) == "contig=<") {
        std::string id;
        std::size_t length = 0;
        const auto inner = body.substr(8, body.size() >= 9 ? body.size() - 9 : 0);
        for (const auto kv : detail::fields(inner, ',')) {
          if (kv.substr(0, 3) == "ID=") id = std::string(kv.substr(3));
          if (kv.substr(0, 7) == "length=") detail::parse_number(kv.substr(7), length);
        }
        if (!id.empty()) out.header.contigs.push_back({id, length});
      }
      return;
    }
    if (line.front() == '#') {
      const auto f = detail::fields(line, '\t');
      if (f.size() < 8 || f[0] != "#CHROM") problem(number, "a header line that is not #CHROM POS ID REF ALT QUAL FILTER INFO");
      for (std::size_t i = 9; i < f.size(); i++) out.header.samples.emplace_back(f[i]);
      header_line = true;
      return;
    }
    if (!header_line) problem(number, "a record before the #CHROM line (read as one)");
    const auto f = detail::fields(line, '\t');
    if (f.size() < 8) { problem(number, "a record of fewer than eight columns"); return; }
    VcfRecord r;
    r.line = number;
    r.chrom = std::string(f[0]);
    std::size_t pos = 0;
    if (!detail::parse_number(f[1], pos) || pos == 0) { problem(number, "POS is not a positive integer"); return; }
    r.position = pos - 1;
    r.id = std::string(f[2]);
    if (f[3].empty() || f[3] == ".") { problem(number, "an empty REF"); return; }
    r.ref = std::string(f[3]);
    if (f[4] != ".")
      for (const auto a : detail::fields(f[4], ',')) r.alts.emplace_back(a);
    if (f[5] != ".") {
      double q = 0;
      if (!detail::parse_number(f[5], q)) { problem(number, "QUAL is not a number"); return; }
      r.qual = q;
    }
    if (f[6] != ".")
      for (const auto x : detail::fields(f[6], ';')) r.filters.emplace_back(x);
    if (f[7] != ".")
      for (const auto kv : detail::fields(f[7], ';')) {
        const std::size_t eq = kv.find('=');
        if (eq == std::string_view::npos) r.info.push_back({std::string(kv), ""});
        else r.info.push_back({std::string(kv.substr(0, eq)), std::string(kv.substr(eq + 1))});
      }
    if (f.size() > 8) {
      for (const auto k : detail::fields(f[8], ':')) r.format.emplace_back(k);
      for (std::size_t i = 9; i < f.size(); i++) {
        std::vector<std::string> v;
        for (const auto x : detail::fields(f[i], ':')) v.emplace_back(x);
        r.samples.push_back(std::move(v));
      }
    }
    if (r.samples.size() != out.header.samples.size()) problem(number, "the record has " + std::to_string(r.samples.size()) + " sample columns, the header " + std::to_string(out.header.samples.size()));
    for (std::size_t i = 0; i < r.samples.size(); i++) {
      const auto gt = r.sample_value(i, "GT");
      bool ok = true;
      if (!gt.empty()) parse_genotype(gt, &ok);
      if (!ok) problem(number, "a GT that is not a genotype: " + std::string(gt));
    }
    out.records.push_back(std::move(r));
  });
  return out;
}

// ---- writing --------------------------------------------------------------------------------------
inline std::string vcf_line(const VcfRecord& r) {
  std::string s = r.chrom + '\t' + std::to_string(r.position + 1) + '\t' + (r.id.empty() ? "." : r.id) + '\t' + r.ref + '\t';
  if (r.alts.empty()) s += '.';
  for (std::size_t i = 0; i < r.alts.size(); i++) s += (i ? "," : "") + r.alts[i];
  s += '\t' + detail::qual_string(r.qual) + '\t';
  if (r.filters.empty()) s += '.';
  for (std::size_t i = 0; i < r.filters.size(); i++) s += (i ? ";" : "") + r.filters[i];
  s += '\t';
  if (r.info.empty()) s += '.';
  for (std::size_t i = 0; i < r.info.size(); i++) s += (i ? ";" : "") + r.info[i].first + (r.info[i].second.empty() ? "" : "=" + r.info[i].second);
  if (!r.format.empty()) {
    s += '\t';
    for (std::size_t i = 0; i < r.format.size(); i++) s += (i ? ":" : "") + r.format[i];
    for (const auto& v : r.samples) {
      s += '\t';
      for (std::size_t i = 0; i < v.size(); i++) s += (i ? ":" : "") + v[i];
    }
  }
  s += '\n';
  return s;
}

inline std::string write_vcf(const VcfHeader& h, const std::vector<VcfRecord>& records) {
  std::string out = "##fileformat=" + h.fileformat + "\n";
  for (const auto& m : h.meta) out += "##" + m + "\n";
  out += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
  if (!h.samples.empty()) {
    out += "\tFORMAT";
    for (const auto& s : h.samples) out += '\t' + s;
  }
  out += '\n';
  for (const auto& r : records) out += vcf_line(r);
  return out;
}

}  // namespace genomics::variant
