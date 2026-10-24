// The features of an annotated genome as tracks an application can draw: from a GenBank
// record (genomics::sequence::read_genbank) and from GFF3.
//
// - Tracks from GenBank. Each feature but the source becomes a compact record: its kind (gene,
//   CDS, rRNA, tRNA, other RNAs, repeat region, mobile element, miscellaneous, other), its key,
//   a name (/gene, else /locus_tag, else /product, else the qualifier that names its kind:
//   /mobile_element_type, /rpt_family, /standard_name, /label, /note), its strand, its parts as
//   segments of the genome, and its extent. The location's parts (INSDC Feature Table, §3.4)
//   are put in the forward strand's order; on a circular genome a part that begins before the
//   one listed before it marks the crossing of the origin (join(2901..3000,1..203)), and the
//   feature's extent [start, stop) then wraps: [start, length) ∪ [0, stop). A CDS is translated
//   by sequence::translate_cds (its /codon_start and /transl_table, the standard code when it
//   has none, INSDC §7.4.5) and compared with its /translation. Features are given lanes,
//   kind by kind, by first fit in order of start, so that a lane's features never overlap.
// - GFF3 (the Sequence Ontology's specification, version 1.26, L. Stein, 18 August 2020): nine
//   tab-separated columns (seqid, source, type, start, end, score, strand, phase, attributes),
//   '.' for an undefined field; percent-encoding (RFC 3986) decoded in every column; attributes
//   tag=value;… with comma-separated values; coordinates 1-based and inclusive, read here as
//   [start − 1, end), with end beyond the landmark's length for a feature that crosses the
//   origin of a circular landmark (Is_circular=true); the directives ##gff-version,
//   ##sequence-region, ### and ##FASTA (or a line beginning with '>'), after which the rest is
//   FASTA. Lines that share an ID are one feature (a discontinuous CDS, say); Parent names
//   the features a feature is part of, and the part-of relation is resolved into a directed
//   acyclic graph of features with its roots (acyclic, not a forest: an exon may be part of
//   several transcripts), a cycle in it being a fault (spec, "Parent (part_of) Relationships").
//
// Readers never throw on content: every fault is a sequence::Problem with its line.
//
// Complexity: O(F log F) for F features (sorting and lanes, first fit over lanes in use); GFF3
// is read in O(text) with hashing of IDs.
//
// Limits: remote parts of a GenBank location (another record's) are left out, as the location
// reader leaves them; the part-of relation is not checked against the Sequence Ontology's
// (that needs the ontology), only for existence and cycles.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/fasta.hpp>
#include <genomics/sequence/genbank.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace genomics::annotation {

using sequence::Code;
using sequence::Codes;
using sequence::Problem;

// ---- tracks of a GenBank record -------------------------------------------------------------------
enum class Kind : std::uint8_t { gene, cds, rrna, trna, rna, repeat, mobile, misc, other };
inline constexpr std::size_t KINDS = 9;

inline const char* kind_name(Kind k) {
  static constexpr const char* names[KINDS] = {"gene", "CDS", "rRNA", "tRNA", "RNA", "repeat", "mobile element", "misc", "other"};
  return names[static_cast<std::size_t>(k)];
}

inline Kind kind_of(std::string_view key) {
  if (key == "gene") return Kind::gene;
  if (key == "CDS") return Kind::cds;
  if (key == "rRNA") return Kind::rrna;
  if (key == "tRNA") return Kind::trna;
  if (key == "ncRNA" || key == "tmRNA" || key == "misc_RNA" || key == "mRNA" || key == "precursor_RNA") return Kind::rna;
  if (key == "repeat_region") return Kind::repeat;
  if (key == "mobile_element" || key == "mobile_element_type") return Kind::mobile;
  if (key == "misc_feature") return Kind::misc;
  return Kind::other;
}

enum class Translation : std::uint8_t { none, agrees, differs, unchecked };  // unchecked: pseudo, an exception, remote parts

struct Segment {
  std::size_t begin = 0, end = 0;  // within [0, length)
};

struct TrackFeature {
  Kind kind = Kind::other;
  std::string key, name, locus_tag, gene, product;
  bool reverse = false, mixed_strands = false, partial = false, pseudo = false;
  std::vector<Segment> segments;     // in the genome's order
  std::size_t start = 0, stop = 0;   // the extent [start, stop), or [start, length) ∪ [0, stop) when it wraps
  bool wraps = false;                // crosses the origin
  std::size_t bases = 0;             // the length of the parts
  Translation translation = Translation::none;
  std::size_t feature = 0;           // its index among the record's features
  std::uint32_t lane = 0;            // within its kind
};

struct Tracks {
  std::string accession;
  std::size_t length = 0;
  bool circular = false;
  std::vector<TrackFeature> features;  // by start
  std::array<std::size_t, KINDS> counts{};
  std::array<std::uint32_t, KINDS> lanes{};  // lanes used by each kind
};

struct TrackOptions {
  bool include_source = false;
  bool check_translations = true;
  int default_table = 1;  // the INSDC default when a CDS has no /transl_table
};

namespace detail {

inline std::string name_of(const sequence::Feature& f) {
  for (const char* q : {"gene", "locus_tag", "product", "mobile_element_type", "rpt_family", "standard_name", "label"}) {
    const std::string v = f.qualifier(q);
    if (!v.empty()) return v;
  }
  std::string note = f.qualifier("note");
  if (note.size() > 48) note = note.substr(0, 45) + "...";
  return note.empty() ? f.key : note;
}

// First fit of features into lanes. A lane is (where it becomes free, where a wrapping
// feature's tail takes it again, or the length). Wrapping features come first, each in a lane of
// its own (any two contain the origin, so they overlap); the others follow in order of start.
inline std::uint32_t place(std::vector<std::pair<std::size_t, std::size_t>>& lanes, const TrackFeature& f, std::size_t length) {
  if (f.wraps) {
    lanes.push_back({f.stop, f.start});
    return static_cast<std::uint32_t>(lanes.size() - 1);
  }
  for (std::size_t l = 0; l < lanes.size(); l++)
    if (lanes[l].first <= f.start && f.stop <= lanes[l].second) { lanes[l].first = f.stop; return static_cast<std::uint32_t>(l); }
  lanes.push_back({f.stop, length});
  return static_cast<std::uint32_t>(lanes.size() - 1);
}

}  // namespace detail

inline Tracks feature_tracks(const sequence::GenBank& record, const TrackOptions& o = {}, std::vector<Problem>* problems = nullptr) {
  Tracks t;
  t.accession = record.version.empty() ? record.accession.empty() ? record.locus : record.accession : record.version;
  t.length = record.length ? record.length : record.sequence.size();
  t.circular = record.circular();
  const Codes genome = o.check_translations ? sequence::encode_dna(record.sequence) : Codes{};
  const auto problem = [&](std::size_t line, std::string m) { if (problems) problems->push_back({line, std::move(m)}); };
  for (std::size_t i = 0; i < record.features.size(); i++) {
    const sequence::Feature& f = record.features[i];
    if (f.key == "source" && !o.include_source) continue;
    const auto& parts = f.location.parts;
    if (parts.empty()) continue;
    TrackFeature tf;
    tf.kind = kind_of(f.key);
    tf.key = f.key;
    tf.name = detail::name_of(f);
    tf.locus_tag = f.qualifier("locus_tag");
    tf.gene = f.qualifier("gene");
    tf.product = f.qualifier("product");
    tf.feature = i;
    tf.pseudo = f.has("pseudo") || f.has("pseudogene");
    tf.reverse = f.location.reverse();
    for (const auto& p : parts) {
      tf.mixed_strands = tf.mixed_strands || p.reverse != parts.front().reverse;
      tf.partial = tf.partial || p.partial_lower || p.partial_upper;
      tf.bases += p.end - p.begin;
    }
    // The parts in the forward strand's order: as read, or reversed for the reverse strand.
    std::vector<sequence::Part> forward(parts.begin(), parts.end());
    if (tf.reverse && !tf.mixed_strands) std::reverse(forward.begin(), forward.end());
    std::size_t crossing = 0;  // the index of the first part after the origin
    if (t.circular && !tf.mixed_strands)
      for (std::size_t j = 1; j < forward.size(); j++)
        if (forward[j].begin < forward[j - 1].begin) { crossing = j; break; }
    tf.wraps = crossing > 0;
    for (const auto& p : forward)
      if (p.end > p.begin) tf.segments.push_back({p.begin, p.end});
    std::sort(tf.segments.begin(), tf.segments.end(), [](const Segment& a, const Segment& b) { return a.begin < b.begin; });
    if (tf.wraps) {
      tf.start = forward.front().begin;
      for (std::size_t j = crossing; j < forward.size(); j++) tf.stop = std::max(tf.stop, forward[j].end);
    } else {
      tf.start = f.location.lower();
      tf.stop = f.location.upper();
    }
    if (f.location.remote) problem(f.line, f.key + " " + tf.name + ": parts in another record are left out of its track");
    if (tf.kind == Kind::cds && f.has("translation")) {
      if (!o.check_translations || tf.pseudo || f.location.remote || f.has("exception") || f.has("ribosomal_slippage") || f.has("transl_except")) {
        tf.translation = Translation::unchecked;
      } else {
        const std::string table = f.qualifier("transl_table");
        const bool known = table.empty() || [&] {
          try { sequence::genetic_code(std::stoi(table)); return true; } catch (...) { return false; }
        }();
        if (!known) {
          tf.translation = Translation::unchecked;
          problem(f.line, "CDS " + tf.name + ": the translation table " + table + " is not one of NCBI's");
        } else {
          tf.translation = sequence::translate_cds(genome, f, o.default_table) == f.qualifier("translation") ? Translation::agrees : Translation::differs;
          if (tf.translation == Translation::differs) problem(f.line, "CDS " + tf.name + " (" + f.location.text + "): its translation differs from its /translation");
        }
      }
    }
    t.features.push_back(std::move(tf));
  }
  std::stable_sort(t.features.begin(), t.features.end(), [](const TrackFeature& a, const TrackFeature& b) { return a.start < b.start; });
  std::array<std::vector<std::pair<std::size_t, std::size_t>>, KINDS> lanes;
  // Wrapping features first, so that they take the ends of their lanes before anything else.
  for (int pass = 0; pass < 2; pass++)
    for (auto& f : t.features) {
      if (f.wraps != (pass == 0)) continue;
      const std::size_t k = static_cast<std::size_t>(f.kind);
      f.lane = detail::place(lanes[k], f, t.length);
    }
  for (const auto& f : t.features) t.counts[static_cast<std::size_t>(f.kind)]++;
  for (std::size_t k = 0; k < KINDS; k++) t.lanes[k] = static_cast<std::uint32_t>(lanes[k].size());
  return t;
}

// ---- GFF3 ---------------------------------------------------------------------------------------------------
struct GffFeature {
  std::string seqid, source, type;
  std::size_t begin = 0, end = 0;  // 0-based, half-open; end beyond the landmark when it crosses the origin
  double score = NAN;              // NaN for '.'
  char strand = '.';               // '+', '-', '.' or '?'
  int phase = -1;                  // 0, 1, 2, or −1 for '.'
  std::vector<std::pair<std::string, std::vector<std::string>>> attributes;  // decoded, in order
  std::string id, name;
  std::vector<std::string> parents;
  std::size_t line = 0;
  std::size_t entity = 0;  // the feature, among Gff::entities, that this line is part of

  const std::vector<std::string>* attribute(std::string_view tag) const {
    for (const auto& [k, v] : attributes)
      if (k == tag) return &v;
    return nullptr;
  }
};

// A feature: the lines that share an ID (or one line without one), with its place in the
// part-of graph (acyclic; an exon may be part of several transcripts).
struct GffEntity {
  std::string id, type, seqid;
  std::vector<std::size_t> lines;     // indices into Gff::features
  std::vector<std::size_t> parents;   // entity indices
  std::vector<std::size_t> children;  // entity indices
  std::size_t begin = 0, end = 0;     // the union of its lines' extents
};

struct SequenceRegion {
  std::string seqid;
  std::size_t begin = 0, end = 0;  // 0-based, half-open
  std::size_t line = 0;
};

struct Gff {
  std::string version;
  std::vector<SequenceRegion> regions;
  std::vector<GffFeature> features;  // one per line, in the file's order
  std::vector<GffEntity> entities;   // in the order of their first lines
  std::vector<std::size_t> roots;    // entities that are part of nothing
  std::vector<std::string> circular; // the seqids of landmarks with Is_circular=true
  std::vector<sequence::Record> fasta;
};

namespace detail {

inline int hex(char c) {
  if (c >= '0' && c <= '9') return c - '0';
  if (c >= 'a' && c <= 'f') return c - 'a' + 10;
  if (c >= 'A' && c <= 'F') return c - 'A' + 10;
  return -1;
}

}  // namespace detail

// RFC 3986 percent-decoding; `ok` becomes false at a '%' not followed by two hex digits (which is kept).
inline std::string percent_decode(std::string_view s, bool& ok) {
  std::string out;
  out.reserve(s.size());
  for (std::size_t i = 0; i < s.size(); i++) {
    if (s[i] == '%') {
      const int a = i + 2 < s.size() ? detail::hex(s[i + 1]) : -1, b = i + 2 < s.size() ? detail::hex(s[i + 2]) : -1;
      if (a >= 0 && b >= 0) {
        out += static_cast<char>(16 * a + b);
        i += 2;
        continue;
      }
      ok = false;
    }
    out += s[i];
  }
  return out;
}

inline Gff read_gff3(std::string_view text, std::vector<Problem>* problems = nullptr) {
  Gff g;
  const auto problem = [&](std::size_t line, std::string m) { if (problems) problems->push_back({line, std::move(m)}); };
  std::unordered_map<std::string, std::size_t> by_id;
  bool fasta = false, first = true;
  std::size_t fasta_line = 0;
  std::string fasta_text;
  sequence::detail::for_each_line(text, [&](std::string_view line, std::size_t number) {
    if (fasta) {
      fasta_text.append(line);
      fasta_text += '\n';
      return;
    }
    if (line.empty()) return;
    const bool was_first = first;
    first = false;
    if (was_first && line.substr(0, 13) != "##gff-version") problem(number, "the first line is not ##gff-version 3");
    if (line.front() == '>') {
      fasta = true;
      fasta_line = number - 1;
      fasta_text.append(line);
      fasta_text += '\n';
      return;
    }
    if (line.substr(0, 2) == "##") {
      std::string_view d = line.substr(2);
      if (d.substr(0, 11) == "gff-version") {
        g.version = std::string(sequence::detail::trim(d.substr(11)));
        if (g.version.empty() || g.version[0] != '3') problem(number, "##gff-version is not 3");
        if (!was_first) problem(number, "##gff-version is not the first line");
      } else if (d.substr(0, 15) == "sequence-region") {
        std::vector<std::string> w;
        std::string cur;
        for (char c : d.substr(15)) {
          if (c == ' ' || c == '\t') { if (!cur.empty()) { w.push_back(cur); cur.clear(); } }
          else cur += c;
        }
        if (!cur.empty()) w.push_back(cur);
        char* e1 = nullptr;
        char* e2 = nullptr;
        if (w.size() != 3) { problem(number, "##sequence-region needs a seqid, a start and an end"); return; }
        const unsigned long long b = std::strtoull(w[1].c_str(), &e1, 10), e = std::strtoull(w[2].c_str(), &e2, 10);
        if (*e1 || *e2 || b < 1 || e < b) { problem(number, "##sequence-region: the range " + w[1] + " " + w[2] + " is not one"); return; }
        bool ok = true;
        const std::string seqid = percent_decode(w[0], ok);
        for (const auto& r : g.regions)
          if (r.seqid == seqid) problem(number, "a second ##sequence-region for " + seqid);
        g.regions.push_back({seqid, static_cast<std::size_t>(b - 1), static_cast<std::size_t>(e), number});
      } else if (d == "FASTA") {
        fasta = true;
        fasta_line = number;
      }
      return;  // ### and other directives
    }
    if (line.front() == '#') return;
    // Nine columns, split on tabs only.
    std::array<std::string_view, 9> col;
    std::size_t n = 0, from = 0;
    for (;;) {
      const std::size_t tab = line.find('\t', from);
      if (n < 9) col[n] = line.substr(from, tab == std::string_view::npos ? std::string_view::npos : tab - from);
      n++;
      if (tab == std::string_view::npos) break;
      from = tab + 1;
    }
    if (n != 9) { problem(number, "a feature line has " + std::to_string(n) + " tab-separated columns, not nine"); return; }
    GffFeature f;
    f.line = number;
    bool encoded = true;
    f.seqid = percent_decode(col[0], encoded);
    f.source = percent_decode(col[1], encoded);
    f.type = percent_decode(col[2], encoded);
    char* e1 = nullptr;
    char* e2 = nullptr;
    const std::string s4(col[3]), s5(col[4]);
    const unsigned long long b = std::strtoull(s4.c_str(), &e1, 10), e = std::strtoull(s5.c_str(), &e2, 10);
    if (s4.empty() || s5.empty() || *e1 || *e2 || b < 1) { problem(number, "start and end must be positive integers: " + s4 + ", " + s5); return; }
    if (e < b) { problem(number, "start " + s4 + " is after end " + s5); return; }
    f.begin = static_cast<std::size_t>(b - 1);
    f.end = static_cast<std::size_t>(e);
    if (col[5] != ".") {
      const std::string sc(col[5]);
      char* es = nullptr;
      f.score = std::strtod(sc.c_str(), &es);
      if (sc.empty() || *es) { problem(number, "the score '" + sc + "' is not a number"); f.score = NAN; }
    }
    if (col[6].size() != 1 || std::string_view("+-.?").find(col[6][0]) == std::string_view::npos) problem(number, "the strand '" + std::string(col[6]) + "' is not +, -, . or ?");
    else f.strand = col[6][0];
    if (col[7] == ".") {
      if (f.type == "CDS" || f.type == "SO:0000316") problem(number, "a CDS without its phase");
    } else if (col[7].size() == 1 && col[7][0] >= '0' && col[7][0] <= '2') {
      f.phase = col[7][0] - '0';
    } else {
      problem(number, "the phase '" + std::string(col[7]) + "' is not 0, 1, 2 or .");
    }
    if (col[8] != "." && !col[8].empty()) {
      std::string_view a = col[8];
      while (!a.empty()) {
        const std::size_t sc = a.find(';');
        std::string_view pair = a.substr(0, sc);
        while (!pair.empty() && pair.front() == ' ') pair.remove_prefix(1);
        if (!pair.empty()) {
          const std::size_t eq = pair.find('=');
          if (eq == std::string_view::npos) problem(number, "the attribute '" + std::string(pair) + "' has no '='");
          else {
            std::string tag = percent_decode(pair.substr(0, eq), encoded);
            std::vector<std::string> values;
            std::string_view v = pair.substr(eq + 1);
            for (;;) {
              const std::size_t comma = v.find(',');
              values.push_back(percent_decode(v.substr(0, comma), encoded));
              if (comma == std::string_view::npos) break;
              v = v.substr(comma + 1);
            }
            if (tag == "ID") {
              if (values.size() != 1) problem(number, "an ID with more than one value");
              f.id = values[0];
            } else if (tag == "Name") {
              f.name = values[0];
            } else if (tag == "Parent") {
              f.parents.insert(f.parents.end(), values.begin(), values.end());
            } else if (tag == "Is_circular" && values[0] == "true") {
              if (std::find(g.circular.begin(), g.circular.end(), f.seqid) == g.circular.end()) g.circular.push_back(f.seqid);
            }
            f.attributes.emplace_back(std::move(tag), std::move(values));
          }
        }
        if (sc == std::string_view::npos) break;
        a = a.substr(sc + 1);
      }
    }
    if (!encoded) problem(number, "a '%' that is not followed by two hexadecimal digits");
    // One feature per ID: lines that share it are its parts.
    const std::size_t index = g.features.size();
    if (!f.id.empty()) {
      const auto it = by_id.find(f.id);
      if (it != by_id.end()) {
        GffEntity& ent = g.entities[it->second];
        if (ent.type != f.type || ent.seqid != f.seqid) problem(number, "the ID " + f.id + " is shared by lines of another type or landmark");
        ent.lines.push_back(index);
        ent.begin = std::min(ent.begin, f.begin);
        ent.end = std::max(ent.end, f.end);
        f.entity = it->second;
        g.features.push_back(std::move(f));
        return;
      }
      by_id.emplace(f.id, g.entities.size());
    }
    GffEntity ent;
    ent.id = f.id;
    ent.type = f.type;
    ent.seqid = f.seqid;
    ent.lines.push_back(index);
    ent.begin = f.begin;
    ent.end = f.end;
    f.entity = g.entities.size();
    g.entities.push_back(std::move(ent));
    g.features.push_back(std::move(f));
  });
  // The part-of relation.
  for (std::size_t i = 0; i < g.entities.size(); i++) {
    GffEntity& ent = g.entities[i];
    for (std::size_t line : ent.lines)
      for (const auto& p : g.features[line].parents) {
        const auto it = by_id.find(p);
        if (it == by_id.end()) {
          problem(g.features[line].line, "the Parent " + p + " is not the ID of any feature");
          continue;
        }
        if (std::find(ent.parents.begin(), ent.parents.end(), it->second) == ent.parents.end()) ent.parents.push_back(it->second);
      }
  }
  for (std::size_t i = 0; i < g.entities.size(); i++) {
    for (std::size_t p : g.entities[i].parents) g.entities[p].children.push_back(i);
    if (g.entities[i].parents.empty()) g.roots.push_back(i);
  }
  // Cycles, by depth-first search with three colours (iterative).
  {
    std::vector<std::uint8_t> colour(g.entities.size(), 0);
    for (std::size_t s = 0; s < g.entities.size(); s++) {
      if (colour[s]) continue;
      std::vector<std::pair<std::size_t, std::size_t>> stack{{s, 0}};
      colour[s] = 1;
      while (!stack.empty()) {
        auto& [v, next] = stack.back();
        if (next < g.entities[v].children.size()) {
          const std::size_t w = g.entities[v].children[next++];
          if (colour[w] == 1) problem(g.features[g.entities[w].lines.front()].line, "the Parent relation has a cycle through " + (g.entities[w].id.empty() ? g.entities[w].type : g.entities[w].id));
          else if (colour[w] == 0) { colour[w] = 1; stack.push_back({w, 0}); }
        } else {
          colour[v] = 2;
          stack.pop_back();
        }
      }
    }
  }
  // Bounds of the sequence regions (a circular landmark's features may run past its end).
  for (const auto& f : g.features)
    for (const auto& r : g.regions)
      if (r.seqid == f.seqid && (f.begin < r.begin || f.end > r.end) && std::find(g.circular.begin(), g.circular.end(), f.seqid) == g.circular.end())
        problem(f.line, f.type + " lies outside the ##sequence-region of " + f.seqid);
  if (fasta) {
    std::vector<Problem> fp;
    g.fasta = sequence::read_fasta(fasta_text, &fp);
    for (auto& p : fp) problem(p.line + fasta_line, p.message);
  }
  return g;
}

// The segments of a GFF3 feature on a landmark of `length`, split at the origin when it crosses it.
inline std::vector<Segment> segments(const Gff& g, const GffEntity& e, std::size_t length) {
  std::vector<Segment> out;
  for (std::size_t line : e.lines) {
    const GffFeature& f = g.features[line];
    if (length && f.end > length) {
      out.push_back({f.begin, length});
      out.push_back({0, f.end - length});
    } else {
      out.push_back({f.begin, f.end});
    }
  }
  std::sort(out.begin(), out.end(), [](const Segment& a, const Segment& b) { return a.begin < b.begin; });
  return out;
}

}  // namespace genomics::annotation
