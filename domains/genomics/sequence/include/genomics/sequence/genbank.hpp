// GenBank flat files: the record of a sequence with its features (NCBI, GenBank Release
// Notes, §3.4; the INSDC Feature Table Definition, version 11.x, §3.4 for locations).
//
// A record is read from its LOCUS line to its "//": the name, length, molecule and topology
// of the LOCUS line; DEFINITION, ACCESSION, VERSION, the ORGANISM and its lineage; the
// FEATURES, each a key (columns 6–20), a location and qualifiers (/name=value from column 22,
// continued over lines: /translation joined without spaces, other text with one); and the
// ORIGIN, the sequence in letters.
//
// A location is a small context-free language, read here by recursive descent:
//
//   location   ::= span | complement(location) | join(location, ...) | order(location, ...)
//                | accession:span                                      (remote; noted, not read)
//   span       ::= [<]b..[>]e | b | b^e | b.e                            (1-based, inclusive)
//
// It becomes a list of parts in the order they are read (complement reverses the order and
// the strands of what it encloses), each part 0-based and half-open on the forward strand,
// with partial ends kept as the lower and upper bounds of the part ('<' and '>' always
// qualify the lower and the upper coordinate). A span a^b, between two bases, has no bases;
// one b.e, somewhere within, is read as the whole range and noted as uncertain.
//
// extract() joins the parts of a location from a genome, reverse-complementing the parts on
// the reverse strand; a join across the origin of a circular genome (join(3981..5386,1..136)
// in phiX174) needs nothing more. translate_cds() applies /codon_start, /transl_table and
// /transl_except (selenocysteine at TGA, stops completed by polyadenylation), and agrees with
// RefSeq's /translation on all 4,522 CDSs of the fourteen genomes of
// data/registry/ncbi-genomes.json (the checks, 2026-10-24).

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/fasta.hpp>
#include <genomics/sequence/genetic_code.hpp>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace genomics::sequence {

struct Part {
  std::size_t begin = 0, end = 0;  // 0-based, half-open, forward strand
  bool reverse = false;
  bool partial_lower = false, partial_upper = false;
};

struct Location {
  std::vector<Part> parts;  // in the order they are read
  std::string kind;         // "span", "join" or "order" (the outermost combinator)
  bool remote = false;      // some part lies in another record, and was left out
  bool uncertain = false;   // a b.e or b^e span
  std::string text;         // as written
  std::size_t lower() const { std::size_t m = SIZE_MAX; for (const auto& p : parts) m = std::min(m, p.begin); return parts.empty() ? 0 : m; }
  std::size_t upper() const { std::size_t m = 0; for (const auto& p : parts) m = std::max(m, p.end); return m; }
  std::size_t length() const { std::size_t n = 0; for (const auto& p : parts) n += p.end - p.begin; return n; }
  bool reverse() const { return !parts.empty() && parts.front().reverse; }
};

struct Feature {
  std::string key;  // gene, CDS, rRNA, tRNA, misc_feature, ...
  Location location;
  std::vector<std::pair<std::string, std::string>> qualifiers;  // in order; a bare /pseudo has value ""
  std::size_t line = 0;

  // The first value of a qualifier, or "" when it is absent.
  std::string qualifier(std::string_view name) const {
    for (const auto& [k, v] : qualifiers)
      if (k == name) return v;
    return "";
  }
  bool has(std::string_view name) const {
    for (const auto& q : qualifiers)
      if (q.first == name) return true;
    return false;
  }
};

struct GenBank {
  std::string locus, molecule, topology, division, date;  // the LOCUS line
  std::size_t length = 0;                                  // as LOCUS gives it
  std::string definition, accession, version, organism;
  std::vector<std::string> lineage;
  std::vector<Feature> features;
  std::string sequence;  // ORIGIN, letters in lower case as written
  bool circular() const { return topology == "circular"; }
};

namespace detail {

class LocationParser {
 public:
  LocationParser(std::string_view s, std::size_t size) : s_(s), size_(size) {}
  bool parse(Location& out, std::string& error) {
    try {
      skip();
      out.parts = location(out);
      skip();
      if (i_ != s_.size()) throw std::string("unexpected '") + s_[i_] + "' at " + std::to_string(i_ + 1);
      return true;
    } catch (const std::string& e) {
      error = e;
      return false;
    }
  }

 private:
  std::vector<Part> location(Location& loc) {
    if (word("complement")) {
      expect('(');
      auto inner = location(loc);
      expect(')');
      std::vector<Part> out(inner.rbegin(), inner.rend());
      for (auto& p : out) p.reverse = !p.reverse;
      if (loc.kind.empty()) loc.kind = "span";
      return out;
    }
    for (const char* combinator : {"join", "order"}) {
      if (word(combinator)) {
        if (loc.kind.empty() || loc.kind == "span") loc.kind = combinator;
        expect('(');
        std::vector<Part> out;
        for (;;) {
          auto inner = location(loc);
          out.insert(out.end(), inner.begin(), inner.end());
          skip();
          if (peek() == ',') { i_++; continue; }
          break;
        }
        expect(')');
        return out;
      }
    }
    // A remote span, accession[.version]:span, is noted and left out.
    const std::size_t colon = s_.find(':', i_);
    const std::size_t stop = s_.find_first_of(",()", i_);
    if (colon != std::string_view::npos && (stop == std::string_view::npos || colon < stop)) {
      loc.remote = true;
      i_ = colon + 1;
      span(loc);
      return {};
    }
    if (loc.kind.empty()) loc.kind = "span";
    return {span(loc)};
  }

  Part span(Location& loc) {
    Part p;
    skip();
    if (peek() == '<') { p.partial_lower = true; i_++; }
    const std::size_t b = number();
    std::size_t e = b;
    if (peek() == '.' && i_ + 1 < s_.size() && s_[i_ + 1] == '.') {
      i_ += 2;
      if (peek() == '>') { p.partial_upper = true; i_++; }
      if (peek() == '<') { p.partial_lower = true; i_++; }  // seen in old records: 1..<20
      e = number();
    } else if (peek() == '>') {
      p.partial_upper = true;
      i_++;
    } else if (peek() == '^') {
      i_++;
      number();
      loc.uncertain = true;
      p.begin = p.end = b;  // between bases b and b + 1: no bases
      return p;
    } else if (peek() == '.') {
      i_++;
      e = number();
      loc.uncertain = true;
    }
    if (b < 1 || e < b) throw std::string("the span ") + std::to_string(b) + ".." + std::to_string(e) + " is not a range of bases";
    if (size_ && e > size_) throw std::string("the span ends at ") + std::to_string(e) + ", beyond the sequence's " + std::to_string(size_);
    p.begin = b - 1;
    p.end = e;
    return p;
  }

  std::size_t number() {
    skip();
    if (i_ >= s_.size() || !std::isdigit(static_cast<unsigned char>(s_[i_]))) throw std::string("expected a number at ") + std::to_string(i_ + 1);
    std::size_t v = 0;
    while (i_ < s_.size() && std::isdigit(static_cast<unsigned char>(s_[i_]))) v = v * 10 + static_cast<std::size_t>(s_[i_++] - '0');
    return v;
  }
  bool word(std::string_view w) {
    skip();
    if (s_.substr(i_, w.size()) == w) {
      std::size_t j = i_ + w.size();
      while (j < s_.size() && s_[j] == ' ') j++;
      if (j < s_.size() && s_[j] == '(') { i_ += w.size(); return true; }
    }
    return false;
  }
  void expect(char c) {
    skip();
    if (peek() != c) throw std::string("expected '") + c + "' at " + std::to_string(i_ + 1);
    i_++;
  }
  char peek() { skip(); return i_ < s_.size() ? s_[i_] : '\0'; }
  void skip() { while (i_ < s_.size() && (s_[i_] == ' ' || s_[i_] == '\t')) i_++; }

  std::string_view s_;
  std::size_t size_ = 0, i_ = 0;
};

inline bool starts_with(std::string_view s, std::string_view p) { return s.substr(0, p.size()) == p; }
inline std::string_view field(std::string_view line) { return line.size() > 12 ? trim(line.substr(12)) : std::string_view{}; }

}  // namespace detail

// Reads a location as written in a feature table; size, when known, bounds it.
inline bool parse_location(std::string_view text, Location& out, std::string& error, std::size_t size = 0) {
  out = Location{};
  out.text = std::string(text);
  return detail::LocationParser(text, size).parse(out, error);
}

// Every record of a GenBank text.
inline std::vector<GenBank> read_genbank(std::string_view text, std::vector<Problem>* problems = nullptr) {
  std::vector<GenBank> out;
  const auto problem = [&](std::size_t line, std::string m) { if (problems) problems->push_back({line, std::move(m)}); };
  enum class Section { none, header, organism, features, origin } section = Section::none;
  GenBank* g = nullptr;
  Feature* feature = nullptr;
  std::string location_text;
  std::size_t location_line = 0;
  std::string* value = nullptr;  // the qualifier being continued
  bool joined = false;           // continued without spaces (/translation)
  bool open_quote = false;
  std::string keyword;           // the last keyword of the header
  bool lineage_open = false;     // the lineage's last name continues on the next line

  const auto finish_location = [&] {
    if (!feature || location_text.empty()) return;
    std::string error;
    if (!parse_location(location_text, feature->location, error, g->length))
      problem(location_line, "feature " + feature->key + ": the location \"" + location_text + "\" cannot be read: " + error);
    location_text.clear();
  };
  const auto finish_value = [&] {
    if (value && value->size() >= 1 && value->front() == '"') {
      std::string v = value->substr(1);
      if (!v.empty() && v.back() == '"') v.pop_back();
      // Inside a quoted value, "" stands for one quotation mark.
      std::string u;
      for (std::size_t i = 0; i < v.size(); i++) {
        u += v[i];
        if (v[i] == '"' && i + 1 < v.size() && v[i + 1] == '"') i++;
      }
      *value = u;
    }
    value = nullptr;
    open_quote = false;
  };

  detail::for_each_line(text, [&](std::string_view line, std::size_t number) {
    if (detail::starts_with(line, "LOCUS")) {
      out.emplace_back();
      g = &out.back();
      feature = nullptr;
      value = nullptr;
      section = Section::header;
      // LOCUS  name  length bp  molecule  topology  division  date (whitespace-separated)
      std::vector<std::string> w;
      std::string cur;
      for (char c : line.substr(5)) {
        if (c == ' ' || c == '\t') { if (!cur.empty()) { w.push_back(cur); cur.clear(); } }
        else cur += c;
      }
      if (!cur.empty()) w.push_back(cur);
      if (!w.empty()) g->locus = w[0];
      for (std::size_t k = 1; k < w.size(); k++) {
        if (k + 1 < w.size() && (w[k + 1] == "bp" || w[k + 1] == "aa")) g->length = std::stoull(w[k]);
        else if (w[k] == "circular" || w[k] == "linear") g->topology = w[k];
        else if (w[k].find("DNA") != std::string::npos || w[k].find("RNA") != std::string::npos) g->molecule = w[k];
      }
      if (w.size() >= 2) g->date = w.back();
      if (w.size() >= 3) g->division = w[w.size() - 2];
      return;
    }
    if (!g) return;
    if (detail::starts_with(line, "//")) {
      finish_value();
      finish_location();
      if (g->length && g->sequence.size() != g->length)
        problem(number, g->locus + ": ORIGIN has " + std::to_string(g->sequence.size()) + " bases where LOCUS says " + std::to_string(g->length));
      g = nullptr;
      section = Section::none;
      return;
    }
    if (!line.empty() && line[0] != ' ') {
      // A keyword of the header, FEATURES or ORIGIN.
      finish_value();
      finish_location();
      feature = nullptr;
      if (detail::starts_with(line, "FEATURES")) { section = Section::features; return; }
      if (detail::starts_with(line, "ORIGIN")) { section = Section::origin; return; }
      section = Section::header;
      keyword = std::string(line.substr(0, line.find(' ')));
      if (detail::starts_with(line, "DEFINITION")) g->definition = std::string(detail::field(line));
      else if (detail::starts_with(line, "ACCESSION")) {
        const auto f = detail::field(line);
        g->accession = std::string(f.substr(0, f.find(' ')));
      } else if (detail::starts_with(line, "VERSION")) {
        const auto f = detail::field(line);
        g->version = std::string(f.substr(0, f.find(' ')));
      }
      return;
    }
    switch (section) {
      case Section::header: {
        if (detail::starts_with(line, "  ORGANISM")) { g->organism = std::string(detail::field(line)); section = Section::organism; lineage_open = false; return; }
        if (line.size() > 2 && line[2] != ' ') { keyword = std::string(detail::trim(line.substr(0, 12))); return; }  // a sub-keyword
        // A continuation of DEFINITION.
        if (keyword == "DEFINITION" && line.size() > 12 && line.substr(0, 12) == std::string(12, ' '))
          g->definition += " " + std::string(detail::field(line));
        return;
      }
      case Section::organism: {
        if (line.size() > 12 && line.substr(0, 12) == std::string(12, ' ')) {
          std::string_view f = detail::field(line);
          bool first = true;
          while (!f.empty()) {
            std::size_t sc = f.find(';');
            std::string_view item = detail::trim(f.substr(0, sc));
            const bool ends = sc != std::string_view::npos || (!item.empty() && item.back() == '.');
            if (!item.empty() && item.back() == '.') item.remove_suffix(1);
            // A name wrapped from the line before (phiX174's "Sinsheimervirus / phiX174") is one name.
            if (!item.empty() && first && lineage_open && !g->lineage.empty()) g->lineage.back() += " " + std::string(item);
            else if (!item.empty()) g->lineage.emplace_back(item);
            lineage_open = !ends;
            first = false;
            if (sc == std::string_view::npos) break;
            f = f.substr(sc + 1);
          }
        } else {
          section = Section::header;
        }
        return;
      }
      case Section::features: {
        if (line.size() < 22) return;
        if (line.substr(0, 5) == "     " && line[5] != ' ') {
          // A new feature: its key in columns 6–20, its location from column 22.
          finish_value();
          finish_location();
          g->features.emplace_back();
          feature = &g->features.back();
          feature->key = std::string(detail::trim(line.substr(5, 16)));
          feature->line = number;
          location_text = std::string(detail::trim(line.substr(21)));
          location_line = number;
          return;
        }
        if (!feature) return;
        std::string_view rest = line.substr(21);
        if (!open_quote && !rest.empty() && rest.front() == '/') {
          finish_value();
          finish_location();
          const std::size_t eq = rest.find('=');
          std::string name(rest.substr(1, eq == std::string_view::npos ? std::string_view::npos : eq - 1));
          std::string v = eq == std::string_view::npos ? "" : std::string(detail::trim(rest.substr(eq + 1)));
          feature->qualifiers.emplace_back(std::move(name), std::move(v));
          value = &feature->qualifiers.back().second;
          joined = feature->qualifiers.back().first == "translation";
          // A quoted value is open until its closing quotation mark (a doubled one does not close it).
          if (!value->empty() && value->front() == '"') {
            std::size_t quotes = 0;
            for (char c : value->substr(1)) quotes += c == '"';
            open_quote = quotes % 2 == 0;
          }
          return;
        }
        if (value) {
          const std::string_view piece = detail::trim(rest);
          if (!joined && !value->empty()) *value += ' ';
          *value += piece;
          std::size_t quotes = 0;
          for (char c : piece) quotes += c == '"';
          if (quotes % 2 == 1) open_quote = !open_quote;
          return;
        }
        // A location continued over lines.
        location_text += std::string(detail::trim(rest));
        return;
      }
      case Section::origin: {
        for (char c : line)
          if (std::isalpha(static_cast<unsigned char>(c))) g->sequence += c;
        return;
      }
      case Section::none: return;
    }
  });
  if (g) problem(0, g->locus + ": the record ends without its \"//\"");
  return out;
}

// The bases of a location, in the order it reads them.
inline Codes extract(const Codes& genome, const Location& location) {
  Codes out;
  out.reserve(location.length());
  for (const auto& p : location.parts) {
    if (p.end > genome.size() || p.begin > p.end) continue;
    if (!p.reverse) out.insert(out.end(), genome.begin() + static_cast<std::ptrdiff_t>(p.begin), genome.begin() + static_cast<std::ptrdiff_t>(p.end));
    else
      for (std::size_t i = p.end; i-- > p.begin;) out.push_back(complement(genome[i]));
  }
  return out;
}

// The protein of a CDS as GenBank writes its /translation: from /codon_start, by its
// /transl_table (or `default_table`), the first codon as M when it is a start and the CDS is
// not partial at its 5' end, every /transl_except applied, the final stop dropped.
inline std::string translate_cds(const Codes& genome, const Feature& cds, int default_table = 1) {
  const Codes x = extract(genome, cds.location);
  const std::string start = cds.qualifier("codon_start"), table = cds.qualifier("transl_table");
  const std::size_t frame = start.empty() ? 0 : static_cast<std::size_t>(std::max(0, std::stoi(start) - 1));
  const GeneticCode& code = genetic_code(table.empty() ? default_table : std::stoi(table));
  const bool five_prime_partial = !cds.location.parts.empty() &&
                                  (cds.location.reverse() ? cds.location.parts.front().partial_upper : cds.location.parts.front().partial_lower);
  std::string p = translate(x, code, frame, false);
  if (!p.empty() && !five_prime_partial && frame + 3 <= x.size() && is_start(code, x[frame], x[frame + 1], x[frame + 2])) p[0] = 'M';
  // /transl_except=(pos:location,aa:Sec): a codon read otherwise than the table says
  // (selenocysteine at a TGA, a stop completed by polyadenylation; INSDC Feature Table §7.3).
  for (const auto& [key, value] : cds.qualifiers) {
    if (key != "transl_except") continue;
    const auto pos = value.find("pos:"), comma = value.find(",aa:");
    if (pos == std::string::npos || comma == std::string::npos) continue;
    Location at;
    std::string error;
    if (!parse_location(value.substr(pos + 4, comma - pos - 4), at, error) || at.parts.empty()) continue;
    // The offset in the CDS of the codon's first base, read on the CDS's strand.
    const Part& first = at.parts.front();
    const std::size_t g = first.reverse ? first.end - 1 : first.begin;
    std::size_t offset = 0;
    bool found = false;
    for (const auto& q : cds.location.parts) {
      if (g >= q.begin && g < q.end) {
        offset += q.reverse ? q.end - 1 - g : g - q.begin;
        found = true;
        break;
      }
      offset += q.end - q.begin;
    }
    if (!found || offset < frame) continue;
    const std::size_t codon = (offset - frame) / 3;
    std::string aa = value.substr(comma + 4);
    while (!aa.empty() && (aa.back() == ')' || aa.back() == ' ')) aa.pop_back();
    static const std::pair<const char*, char> three[] = {
        {"Ala", 'A'}, {"Arg", 'R'}, {"Asn", 'N'}, {"Asp", 'D'}, {"Cys", 'C'}, {"Gln", 'Q'}, {"Glu", 'E'}, {"Gly", 'G'}, {"His", 'H'}, {"Ile", 'I'},
        {"Leu", 'L'}, {"Lys", 'K'}, {"Met", 'M'}, {"Phe", 'F'}, {"Pro", 'P'}, {"Ser", 'S'}, {"Thr", 'T'}, {"Trp", 'W'}, {"Tyr", 'Y'}, {"Val", 'V'},
        {"Sec", 'U'}, {"Pyl", 'O'}, {"Asx", 'B'}, {"Glx", 'Z'}, {"Xle", 'J'}, {"TERM", '*'}, {"OTHER", 'X'}};
    for (const auto& [name, letter] : three)
      if (aa == name) {
        if (codon < p.size()) p[codon] = letter;
        else if (codon == p.size()) p += letter;  // an incomplete final codon
      }
  }
  if (!p.empty() && p.back() == '*') p.pop_back();
  return p;
}

}  // namespace genomics::sequence
