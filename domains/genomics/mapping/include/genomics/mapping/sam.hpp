// SAM, the text form of alignments (The SAM/BAM Format Specification, version 1.6, SAM/BAM
// Format Specification Working Group, 2024 revision; Li et al. 2009, Bioinformatics 25,
// 2078–2079), written from mappings and read back as records.
//
// - Header (§1.3): @HD with VN:1.6 and SO:unsorted; one @SQ per reference sequence with SN and
//   LN; one @PG with ID, PN and, when given, CL.
// - Records (§1.4), eleven mandatory fields separated by tabs: QNAME; FLAG (0x4 unmapped, 0x10
//   the reverse strand, 0x100 a secondary alignment, 0x800 a supplementary one); RNAME; POS, the
//   1-based leftmost reference base (0 when unmapped); MAPQ; CIGAR, with the soft clips (S) of
//   the read's unaligned ends, so that the operations consuming the query sum to the length of
//   SEQ; RNEXT, PNEXT and TLEN for single-end reads ('*', 0, 0); SEQ and QUAL on the strand
//   mapped (the reverse complement and the reversed qualities for FLAG 0x10, §1.4.10–11), QUAL
//   Phred + 33 or '*'. Optional fields (§1.5) TAG:TYPE:VALUE: NM:i (edit distance), AS:i
//   (alignment score), tp:A (P primary or supplementary, S secondary, as minimap2 writes it),
//   s1:i (the chaining score, minimap2's).
// - Reading never throws on content: a line that is not a record (fewer than eleven fields, a
//   number that is not one, a CIGAR that is not one) is a problem with its line number, and the
//   rest is read on; header lines are kept, and @SQ's names and lengths read.
//
// Complexity: linear in the text. Coordinates inside are 0-based; POS is 1-based in the text only.

#pragma once

#include <genomics/align/types.hpp>
#include <genomics/mapping/mapper.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/fasta.hpp>

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace genomics::mapping {

inline constexpr std::uint16_t SAM_UNMAPPED = 0x4, SAM_REVERSE = 0x10, SAM_SECONDARY = 0x100, SAM_SUPPLEMENTARY = 0x800;

struct SamTag {
  std::string name;  // two characters
  char type = 'Z';   // A, i, f, Z, H, B
  std::string value;
  friend bool operator==(const SamTag&, const SamTag&) = default;
};

struct SamRecord {
  std::string qname;
  std::uint16_t flag = 0;
  std::string rname = "*";
  std::int64_t position = -1;  // 0-based; −1 when POS is 0
  int mapq = 255;
  align::Cigar cigar;          // empty for '*'
  std::string rnext = "*";
  std::int64_t pnext = -1;     // 0-based; −1 when PNEXT is 0
  std::int64_t tlen = 0;
  std::string seq = "*", qual = "*";
  std::vector<SamTag> tags;
  std::size_t line = 0;        // where it was read

  bool unmapped() const { return flag & SAM_UNMAPPED; }
  bool reverse() const { return flag & SAM_REVERSE; }
  bool secondary() const { return flag & SAM_SECONDARY; }
  bool supplementary() const { return flag & SAM_SUPPLEMENTARY; }
  const SamTag* tag(std::string_view name) const {
    for (const auto& t : tags)
      if (t.name == name) return &t;
    return nullptr;
  }
  // An integer tag's value, or `otherwise`.
  std::int64_t tag_int(std::string_view name, std::int64_t otherwise = 0) const {
    const SamTag* t = tag(name);
    if (!t || t->type != 'i') return otherwise;
    std::int64_t v = otherwise;
    std::from_chars(t->value.data(), t->value.data() + t->value.size(), v);
    return v;
  }
};

struct Sam {
  std::vector<std::string> header;                              // the header lines, without their ending
  std::vector<std::pair<std::string, std::size_t>> references;  // @SQ: SN and LN, in order
  std::vector<SamRecord> records;
};

// ---- writing ---------------------------------------------------------------------------------------
inline std::string sam_header(const std::vector<std::string>& names, const std::vector<std::size_t>& lengths, std::string_view program = "genomics-mapping",
                              std::string_view command_line = "") {
  if (names.size() != lengths.size()) throw std::invalid_argument("genomics::mapping::sam_header: one length for each name");
  std::string h = "@HD\tVN:1.6\tSO:unsorted\n";
  for (std::size_t i = 0; i < names.size(); i++) h += "@SQ\tSN:" + names[i] + "\tLN:" + std::to_string(lengths[i]) + "\n";
  h += "@PG\tID:" + std::string(program) + "\tPN:" + std::string(program);
  if (!command_line.empty()) h += "\tCL:" + std::string(command_line);
  h += "\n";
  return h;
}

// The record of a mapping of a read (its bases and qualities as sequenced; quality may be empty).
inline SamRecord sam_record(const Mapping& m, std::string_view name, const Codes& read, std::string_view quality, const std::vector<std::string>& contig_names) {
  if (m.contig >= contig_names.size()) throw std::invalid_argument("genomics::mapping::sam_record: a contig without a name");
  if (read.size() != m.query_length) throw std::invalid_argument("genomics::mapping::sam_record: the read is not the one mapped");
  SamRecord r;
  r.qname = std::string(name);
  r.flag = static_cast<std::uint16_t>((m.reverse ? SAM_REVERSE : 0) | (m.primary ? 0 : SAM_SECONDARY) | (m.supplementary ? SAM_SUPPLEMENTARY : 0));
  r.rname = contig_names[m.contig];
  r.position = static_cast<std::int64_t>(m.target_begin);
  r.mapq = m.mapq;
  r.cigar = m.cigar();
  r.seq = sequence::decode_dna(m.reverse ? sequence::reverse_complement(read) : read);
  if (quality.size() == read.size()) {
    r.qual = std::string(quality);
    if (m.reverse) std::reverse(r.qual.begin(), r.qual.end());
  }
  r.tags.push_back({"NM", 'i', std::to_string(m.nm)});
  r.tags.push_back({"AS", 'i', std::to_string(m.alignment.score)});
  r.tags.push_back({"tp", 'A', m.primary ? "P" : "S"});
  r.tags.push_back({"s1", 'i', std::to_string(m.chain_score)});
  return r;
}

inline SamRecord sam_unmapped(std::string_view name, const Codes& read, std::string_view quality) {
  SamRecord r;
  r.qname = std::string(name);
  r.flag = SAM_UNMAPPED;
  r.mapq = 0;
  r.seq = sequence::decode_dna(read);
  if (quality.size() == read.size()) r.qual = std::string(quality);
  return r;
}

inline std::string sam_line(const SamRecord& r) {
  std::string s;
  s.reserve(r.seq.size() * 2 + 96);
  s += r.qname.empty() ? "*" : r.qname;
  s += '\t' + std::to_string(r.flag) + '\t' + r.rname + '\t' + std::to_string(r.position + 1) + '\t' + std::to_string(r.mapq) + '\t' + align::cigar_string(r.cigar);
  s += '\t' + r.rnext + '\t' + std::to_string(r.pnext + 1) + '\t' + std::to_string(r.tlen) + '\t' + r.seq + '\t' + r.qual;
  for (const auto& t : r.tags) s += '\t' + t.name + ':' + t.type + ':' + t.value;
  s += '\n';
  return s;
}

inline std::string write_sam(std::string_view header, const std::vector<SamRecord>& records) {
  std::string out(header);
  for (const auto& r : records) out += sam_line(r);
  return out;
}

// ---- reading ---------------------------------------------------------------------------------------
namespace detail {

// Calls f(line, number) for each line of the text, without its ending (\n or \r\n).
template <class F>
void lines(std::string_view text, F&& f) {
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

template <class T>
bool number(std::string_view s, T& out) {
  if (s.empty()) return false;
  const auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), out);
  return ec == std::errc() && p == s.data() + s.size();
}

inline std::vector<std::string_view> split(std::string_view line, char sep) {
  std::vector<std::string_view> f;
  std::size_t i = 0;
  for (;;) {
    const std::size_t j = line.find(sep, i);
    f.push_back(line.substr(i, j == std::string_view::npos ? std::string_view::npos : j - i));
    if (j == std::string_view::npos) break;
    i = j + 1;
  }
  return f;
}

}  // namespace detail

inline Sam read_sam(std::string_view text, std::vector<sequence::Problem>* problems = nullptr) {
  Sam out;
  const auto problem = [&](std::size_t line, std::string m) { if (problems) problems->push_back({line, std::move(m)}); };
  detail::lines(text, [&](std::string_view line, std::size_t number) {
    if (line.empty()) return;
    if (line.front() == '@') {
      out.header.emplace_back(line);
      if (line.substr(0, 3) == "@SQ") {
        std::string name;
        std::size_t length = 0;
        bool has_length = false;
        for (const auto f : detail::split(line, '\t')) {
          if (f.substr(0, 3) == "SN:") name = std::string(f.substr(3));
          if (f.substr(0, 3) == "LN:") has_length = detail::number(f.substr(3), length);
        }
        if (name.empty() || !has_length) problem(number, "an @SQ line without SN or LN");
        else out.references.push_back({name, length});
      }
      return;
    }
    const auto f = detail::split(line, '\t');
    if (f.size() < 11) { problem(number, "a record of fewer than eleven fields"); return; }
    SamRecord r;
    r.line = number;
    r.qname = std::string(f[0]);
    unsigned flag = 0;
    std::int64_t pos = 0, pnext = 0;
    if (!detail::number(f[1], flag) || flag > 0xFFFF) { problem(number, "FLAG is not an integer below 65536"); return; }
    r.flag = static_cast<std::uint16_t>(flag);
    r.rname = std::string(f[2]);
    if (!detail::number(f[3], pos) || pos < 0) { problem(number, "POS is not a nonnegative integer"); return; }
    r.position = pos - 1;
    if (!detail::number(f[4], r.mapq) || r.mapq < 0 || r.mapq > 255) { problem(number, "MAPQ is not an integer in [0, 255]"); return; }
    try {
      r.cigar = align::parse_cigar(f[5]);
    } catch (const std::invalid_argument&) {
      problem(number, "CIGAR is not one: " + std::string(f[5]));
      return;
    }
    r.rnext = std::string(f[6]);
    if (!detail::number(f[7], pnext) || pnext < 0) { problem(number, "PNEXT is not a nonnegative integer"); return; }
    r.pnext = pnext - 1;
    if (!detail::number(f[8], r.tlen)) { problem(number, "TLEN is not an integer"); return; }
    r.seq = std::string(f[9]);
    r.qual = std::string(f[10]);
    if (r.seq != "*" && !r.cigar.empty() && align::query_length(r.cigar) != r.seq.size())
      problem(number, "the CIGAR's query length differs from SEQ's (kept)");
    if (r.qual != "*" && r.seq != "*" && r.qual.size() != r.seq.size()) problem(number, "QUAL and SEQ differ in length (kept)");
    for (std::size_t i = 11; i < f.size(); i++) {
      const auto t = f[i];
      if (t.size() < 5 || t[2] != ':' || t[4] != ':') { problem(number, "an optional field not of the form TAG:TYPE:VALUE: " + std::string(t)); continue; }
      r.tags.push_back({std::string(t.substr(0, 2)), t[3], std::string(t.substr(5))});
    }
    out.records.push_back(std::move(r));
  });
  return out;
}

// The mapping a record describes, against the references named in order (contig i is names[i]);
// false when the record is unmapped or names an unknown reference.
inline bool mapping_of(const SamRecord& r, const std::vector<std::string>& names, Mapping& m) {
  if (r.unmapped() || r.position < 0) return false;
  const auto it = std::find(names.begin(), names.end(), r.rname);
  if (it == names.end()) return false;
  m = Mapping{};
  m.contig = static_cast<std::uint32_t>(it - names.begin());
  m.reverse = r.reverse();
  m.primary = !r.secondary();
  m.supplementary = r.supplementary();
  m.mapq = r.mapq;
  m.query_length = align::query_length(r.cigar);
  std::size_t left = 0, right = 0;
  align::Cigar core = r.cigar;
  if (!core.empty() && core.front().op == 'S') { left = core.front().length; core.erase(core.begin()); }
  if (!core.empty() && core.back().op == 'S') { right = core.back().length; core.pop_back(); }
  m.alignment.cigar = std::move(core);
  m.alignment.query_begin = left;
  m.alignment.query_end = m.query_length - right;
  m.alignment.target_begin = static_cast<std::size_t>(r.position);
  m.alignment.target_end = m.alignment.target_begin + align::target_length(m.alignment.cigar);
  m.alignment.score = r.tag_int("AS");
  m.nm = static_cast<std::size_t>(r.tag_int("NM"));
  m.chain_score = r.tag_int("s1");
  m.target_begin = m.alignment.target_begin;
  m.target_end = m.alignment.target_end;
  return true;
}

}  // namespace genomics::mapping
