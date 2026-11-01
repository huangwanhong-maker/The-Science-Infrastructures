// The human ideogram: chromosomes, arms, regions, bands and sub-bands, their coordinates on the
// assembly and their Giemsa stains, read from NCBI's Genome Decoration Page or UCSC's cytoBand.
//
// - Names (ISCN 2024 §2; ISCN 2020 ch. 2). A band is named by its chromosome, its arm (p, the
//   short arm; q, the long) and a word of digits read outward from the centromere: region, band,
//   then sub-bands after a point, "1p36.33" being region 3, band 6, sub-bands 3 and 3 of the short
//   arm of chromosome 1. A band seen at a higher resolution is divided into sub-bands that extend
//   its name, so the hierarchy of the bands is the prefix order on the words of digits: 1p3 ⊇
//   1p36 ⊇ 1p36.3 ⊇ 1p36.33. Here the word is held without its point ("3633"); the point is
//   written after the second digit.
// - Files. NCBI's ideogram files (ftp.ncbi.nlm.nih.gov/pub/gdp, README of 2012) give one band a
//   line, tab-separated: chromosome, arm, band, the ISCN ideogram's start and stop, the band's
//   first and last base (1-based, closed), stain (gneg, gpos, acen, gvar, stalk) and density
//   (25, 50, 75, 100 of gpos). For GRCh38 (assembly unit GCF_000001305.16) NCBI gives the levels
//   of 400, 550 and 850 bands; ISCN's 300-band level is not published. UCSC's cytoBand table
//   (chrom, chromStart, chromEnd, name, gieStain; 0-based, half-open) of hg38 holds the same 862
//   bands as NCBI's 850 level with the same coordinates (compared on 2026-11-01).
// - Inside, coordinates are 0-based and half-open, as everywhere in the domain. The centromere is
//   the boundary between the last band of p and the first of q (between the two acen bands), and
//   the breakpoints p10 and q10 of whole-arm rearrangements are that point (ISCN 2024 §5).
// - Readers never throw on content: a line that cannot be read is a problem with its line; a
//   chromosome whose bands leave a gap or overlap is a problem too, and is kept.
//
// Complexity: reading is linear in the text; finding a band by name or by position is a binary
// search, O(log n) in the bands of the chromosome.

#pragma once

#include <genomics/sequence/fasta.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace genomics::karyotype {

using sequence::Problem;

// Chromosomes 1–22, X and Y as the indices 0–23.
inline constexpr int chromosome_count = 24;
inline constexpr int chr_x = 22, chr_y = 23;

inline int chromosome_index(std::string_view name) {
  if (name.size() > 3 && (name.substr(0, 3) == "chr" || name.substr(0, 3) == "Chr")) name.remove_prefix(3);
  if (name == "X") return chr_x;
  if (name == "Y") return chr_y;
  if (name.empty() || name.size() > 2) return -1;
  int v = 0;
  for (char c : name) {
    if (c < '0' || c > '9') return -1;
    v = v * 10 + (c - '0');
  }
  if (name[0] == '0') return -1;
  return v >= 1 && v <= 22 ? v - 1 : -1;
}
inline std::string chromosome_name(int index) {
  if (index == chr_x) return "X";
  if (index == chr_y) return "Y";
  if (index < 0 || index >= chromosome_count) return "?";
  return std::to_string(index + 1);
}
inline bool is_sex(int index) { return index == chr_x || index == chr_y; }

// Giemsa stains: negative, positive by density, centromeric, variable heterochromatin, stalks of
// the acrocentric short arms.
enum class Stain : std::uint8_t { gneg, gpos25, gpos50, gpos75, gpos100, acen, gvar, stalk };
inline std::string_view stain_name(Stain s) {
  static constexpr std::string_view names[] = {"gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar", "stalk"};
  return names[static_cast<int>(s)];
}
inline std::optional<Stain> stain_from(std::string_view s, std::string_view density = {}) {
  if (s == "gpos") {
    if (density == "25") return Stain::gpos25;
    if (density == "50") return Stain::gpos50;
    if (density == "75") return Stain::gpos75;
    if (density == "100" || density.empty()) return Stain::gpos100;
    return std::nullopt;
  }
  static constexpr std::string_view names[] = {"gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar", "stalk"};
  for (int i = 0; i < 8; i++)
    if (s == names[i]) return static_cast<Stain>(i);
  return std::nullopt;
}

// A band's name from its arm and word of digits: ('p', "3633") → "p36.33".
inline std::string band_label(char arm, std::string_view digits) {
  std::string out(1, arm);
  if (digits.size() <= 2) return out + std::string(digits);
  return out + std::string(digits.substr(0, 2)) + "." + std::string(digits.substr(2));
}
// The word of digits of a written band ("36.33" → "3633"), or nullopt when it is not one: digits,
// with at most one point, after the second digit.
inline std::optional<std::string> band_digits(std::string_view written) {
  std::string d;
  bool point = false;
  for (std::size_t i = 0; i < written.size(); i++) {
    const char c = written[i];
    if (c == '.') {
      if (point || i != 2) return std::nullopt;
      point = true;
    } else if (c >= '0' && c <= '9') {
      d += c;
    } else {
      return std::nullopt;
    }
  }
  if (d.empty() || (point && d.size() < 3)) return std::nullopt;
  return d;
}

struct Band {
  int chromosome = 0;
  char arm = 'p';
  std::string digits;                  // "3633"
  std::uint64_t start = 0, end = 0;    // on the assembly, 0-based, half-open
  double iscn_start = 0, iscn_end = 0; // on the ideogram of its level
  Stain stain = Stain::gneg;
  std::string name() const { return band_label(arm, digits); }                          // "p36.33"
  std::string full_name() const { return chromosome_name(chromosome) + name(); }       // "1p36.33"
  std::uint64_t length() const { return end - start; }
};

// A run of consecutive bands of one chromosome, [first, last] in the ideogram's order, and the
// interval of the assembly they cover.
struct BandRange {
  std::size_t first = 0, last = 0;
  std::uint64_t start = 0, end = 0;
};

struct Ideogram {
  std::string assembly = "GRCh38";
  int level = 850;
  std::vector<Band> bands;  // by chromosome (1–22, X, Y), then from pter to qter
  std::array<std::size_t, chromosome_count + 1> first{};  // the bands of chromosome c are [first[c], first[c + 1])
  std::array<std::uint64_t, chromosome_count> length{};      // the end of its last band
  std::array<std::uint64_t, chromosome_count> centromere{};  // the boundary between its arms
  std::array<double, chromosome_count> iscn_length{};        // its length on the ideogram

  bool has(int c) const { return c >= 0 && c < chromosome_count && first[c + 1] > first[c]; }
  std::size_t begin_of(int c) const { return first[c]; }
  std::size_t end_of(int c) const { return first[c + 1]; }

  // The bands of chromosome c on arm `arm` whose words begin with `digits` (all of them for an
  // empty word): a band, region or sub-band at any resolution, as one run. Nullopt when none.
  std::optional<BandRange> find(int c, char arm, std::string_view digits) const {
    if (!has(c)) return std::nullopt;
    std::optional<BandRange> out;
    for (std::size_t i = first[c]; i < first[c + 1]; i++) {
      const Band& b = bands[i];
      if (b.arm != arm || b.digits.compare(0, digits.size(), digits) != 0) continue;
      if (!out) out = BandRange{i, i, b.start, b.end};
      else { out->last = i; out->start = std::min(out->start, b.start); out->end = std::max(out->end, b.end); }
    }
    return out;
  }
  // The same by a written name: "p36.33", "q21", "p1", "q".
  std::optional<BandRange> find(int c, std::string_view name) const {
    if (name.empty() || (name[0] != 'p' && name[0] != 'q')) return std::nullopt;
    if (name.size() == 1) return find(c, name[0], "");
    const auto d = band_digits(name.substr(1));
    if (!d) return std::nullopt;
    return find(c, name[0], *d);
  }
  // The arm p or q of chromosome c.
  std::optional<BandRange> arm(int c, char a) const { return find(c, a, ""); }

  // The band holding position x of chromosome c (x < length; the last band for x = length).
  std::size_t band_at(int c, std::uint64_t x) const {
    if (!has(c)) throw std::invalid_argument("genomics::karyotype::Ideogram::band_at: no chromosome " + chromosome_name(c));
    const auto lo = bands.begin() + static_cast<std::ptrdiff_t>(first[c]);
    const auto hi = bands.begin() + static_cast<std::ptrdiff_t>(first[c + 1]);
    auto it = std::upper_bound(lo, hi, x, [](std::uint64_t v, const Band& b) { return v < b.end; });
    if (it == hi) --it;
    return static_cast<std::size_t>(it - bands.begin());
  }
  // The bands covering [start, end) of chromosome c (at least one band).
  BandRange bands_over(int c, std::uint64_t start, std::uint64_t end) const {
    const std::size_t a = band_at(c, start);
    const std::size_t b = band_at(c, end > start ? end - 1 : start);
    return {a, b, bands[a].start, bands[b].end};
  }
  // The ISCN name of the bands covering [start, end): one band ("q21.1"), or the first and the
  // last ("q21.1q21.2") as microarray descriptions write them (ISCN 2024 ch. 8).
  std::string range_name(int c, std::uint64_t start, std::uint64_t end) const {
    const BandRange r = bands_over(c, start, end);
    if (r.first == r.last) return bands[r.first].name();
    return bands[r.first].name() + bands[r.last].name();
  }
  // The smallest named band, at any resolution, containing [start, end) of one arm: the longest
  // common prefix of the words of the bands it covers; nullopt when it crosses the centromere.
  std::optional<std::string> smallest_band(int c, std::uint64_t start, std::uint64_t end) const {
    const BandRange r = bands_over(c, start, end);
    const Band& a = bands[r.first];
    const Band& b = bands[r.last];
    if (a.arm != b.arm) return std::nullopt;
    std::size_t k = 0;
    while (k < a.digits.size() && k < b.digits.size() && a.digits[k] == b.digits[k]) k++;
    // A prefix of the digits names a band only when every band of the run shares it, which a
    // common prefix of the first and the last guarantees (prefix classes are runs).
    return band_label(a.arm, a.digits.substr(0, k));
  }
  // The bands of chromosome c as a name → index lookup is find(); this is the terminal band.
  std::size_t pter_band(int c) const { return first[c]; }
  std::size_t qter_band(int c) const { return first[c + 1] - 1; }
};

namespace detail {

inline bool parse_u64(std::string_view s, std::uint64_t& out) {
  if (s.empty()) return false;
  const auto r = std::from_chars(s.data(), s.data() + s.size(), out);
  return r.ec == std::errc() && r.ptr == s.data() + s.size();
}
inline bool parse_double(std::string_view s, double& out) {
  if (s.empty()) return false;
  std::uint64_t whole = 0;
  const std::size_t dot = s.find('.');
  if (!parse_u64(s.substr(0, dot), whole)) return false;
  out = static_cast<double>(whole);
  if (dot != std::string_view::npos) {
    double scale = 0.1;
    for (char c : s.substr(dot + 1)) {
      if (c < '0' || c > '9') return false;
      out += (c - '0') * scale;
      scale /= 10;
    }
  }
  return true;
}
inline std::vector<std::string_view> split_tabs(std::string_view line) {
  std::vector<std::string_view> f;
  std::size_t b = 0;
  for (;;) {
    const std::size_t e = line.find('\t', b);
    f.push_back(line.substr(b, e == std::string_view::npos ? std::string_view::npos : e - b));
    if (e == std::string_view::npos) break;
    b = e + 1;
  }
  return f;
}

// Sorts the bands by chromosome and position, fills the indexes, and reports gaps and overlaps.
inline void finish(Ideogram& g, std::vector<Problem>* problems) {
  std::stable_sort(g.bands.begin(), g.bands.end(), [](const Band& a, const Band& b) {
    return a.chromosome != b.chromosome ? a.chromosome < b.chromosome : a.start < b.start;
  });
  g.first.fill(0);
  std::size_t i = 0;
  for (int c = 0; c < chromosome_count; c++) {
    g.first[c] = i;
    std::uint64_t at = 0;
    bool seen_q = false;
    g.centromere[c] = 0;
    double iscn = 0;
    while (i < g.bands.size() && g.bands[i].chromosome == c) {
      const Band& b = g.bands[i];
      if (b.start != at && problems)
        problems->push_back({0, "chromosome " + chromosome_name(c) + ": band " + b.name() + " begins at " + std::to_string(b.start) +
                                    (b.start > at ? ", leaving a gap after " : ", overlapping the band ending at ") + std::to_string(at)});
      if (b.arm == 'q' && !seen_q) { g.centromere[c] = b.start; seen_q = true; }
      if (b.arm == 'p' && seen_q && problems)
        problems->push_back({0, "chromosome " + chromosome_name(c) + ": band " + b.name() + " of the short arm lies after the long arm's first band"});
      at = b.end;
      iscn = std::max(iscn, b.iscn_end);
      i++;
    }
    g.length[c] = at;
    g.iscn_length[c] = iscn;
    if (!seen_q) g.centromere[c] = at;
  }
  g.first[chromosome_count] = i;
}

}  // namespace detail

// Reads one of NCBI's ideogram files (Genome Decoration Page): `level` is the band level its name
// gives (400, 550, 850).
inline Ideogram read_ideogram(std::string_view text, int level = 850, std::vector<Problem>* problems = nullptr,
                              std::string assembly = "GRCh38") {
  Ideogram g;
  g.level = level;
  g.assembly = std::move(assembly);
  std::size_t line_no = 0, b = 0;
  while (b < text.size()) {
    std::size_t e = text.find('\n', b);
    if (e == std::string_view::npos) e = text.size();
    std::string_view line = text.substr(b, e - b);
    b = e + 1;
    line_no++;
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    if (line.empty() || line[0] == '#') continue;
    const auto f = detail::split_tabs(line);
    const auto fail = [&](const std::string& why) {
      if (problems) problems->push_back({line_no, why});
    };
    if (f.size() < 8) { fail("fewer than eight columns"); continue; }
    Band band;
    band.chromosome = chromosome_index(f[0]);
    if (band.chromosome < 0) { fail("not a chromosome of the human ideogram: " + std::string(f[0])); continue; }
    if (f[1] != "p" && f[1] != "q") { fail("the arm is neither p nor q: " + std::string(f[1])); continue; }
    band.arm = f[1][0];
    const auto d = band_digits(f[2]);
    if (!d) { fail("not a band's name: " + std::string(f[2])); continue; }
    band.digits = *d;
    std::uint64_t first = 0, last = 0;
    if (!detail::parse_u64(f[5], first) || !detail::parse_u64(f[6], last) || last < first) {
      fail("the base-pair positions are not start ≤ stop");
      continue;
    }
    // NCBI's files differ in their starts: the first base of the band (1-based, as the README of
    // the Genome Decoration Page says), or the previous band's stop (the GRCh37 files of the 400
    // and 850 levels). Each start is read against the band before it.
    const bool continues = !g.bands.empty() && g.bands.back().chromosome == band.chromosome;
    if (!continues) band.start = first <= 1 ? 0 : first - 1;
    else if (first == g.bands.back().end) band.start = first;
    else band.start = first == 0 ? 0 : first - 1;
    band.end = last;
    if (!detail::parse_double(f[3], band.iscn_start) || !detail::parse_double(f[4], band.iscn_end)) {
      fail("the ISCN positions are not numbers");
      continue;
    }
    const auto st = stain_from(f[7], f.size() > 8 ? f[8] : std::string_view{});
    if (!st) { fail("not a stain: " + std::string(f[7])); continue; }
    band.stain = *st;
    g.bands.push_back(std::move(band));
  }
  detail::finish(g, problems);
  return g;
}

// Reads UCSC's cytoBand table (its level is that of the bands it holds; 850 for hg38). Rows of
// contigs other than chr1–22, chrX and chrY are skipped.
inline Ideogram read_cytoband(std::string_view text, int level = 850, std::vector<Problem>* problems = nullptr,
                              std::string assembly = "GRCh38") {
  Ideogram g;
  g.level = level;
  g.assembly = std::move(assembly);
  std::size_t line_no = 0, b = 0;
  while (b < text.size()) {
    std::size_t e = text.find('\n', b);
    if (e == std::string_view::npos) e = text.size();
    std::string_view line = text.substr(b, e - b);
    b = e + 1;
    line_no++;
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    if (line.empty() || line[0] == '#') continue;
    const auto f = detail::split_tabs(line);
    if (f.size() < 5) { if (problems) problems->push_back({line_no, "fewer than five columns"}); continue; }
    const int c = chromosome_index(f[0]);
    if (c < 0) continue;
    Band band;
    band.chromosome = c;
    if (f[3].empty() || (f[3][0] != 'p' && f[3][0] != 'q')) { if (problems) problems->push_back({line_no, "not a band's name: " + std::string(f[3])}); continue; }
    band.arm = f[3][0];
    const auto d = band_digits(f[3].substr(1));
    if (!d) { if (problems) problems->push_back({line_no, "not a band's name: " + std::string(f[3])}); continue; }
    band.digits = *d;
    if (!detail::parse_u64(f[1], band.start) || !detail::parse_u64(f[2], band.end) || band.end <= band.start) {
      if (problems) problems->push_back({line_no, "the positions are not start < end"});
      continue;
    }
    const auto st = stain_from(f[4]);
    if (!st) { if (problems) problems->push_back({line_no, "not a stain: " + std::string(f[4])}); continue; }
    band.stain = *st;
    band.iscn_start = static_cast<double>(band.start);
    band.iscn_end = static_cast<double>(band.end);
    g.bands.push_back(std::move(band));
  }
  detail::finish(g, problems);
  return g;
}

}  // namespace genomics::karyotype
