// Checks of the HGVS grammar, its value and its writing (grammar.hpp, hgvs.hpp), on every example of
// the HGVS Nomenclature 21.1 (tests/fixtures/hgvs-examples.tsv, made by make_hgvs_examples.py from
// the standard's sources at a pinned commit):
//
// - every example the standard gives is in the language, with no fault of meaning, and has exactly
//   one derivation (counted by an independent counter of derivations over the compiled grammar, ε
//   rules included, saturating at 2);
// - every location the standard writes alone (c.88+1, g.(123_456)_(789_1000)) is read by its family's
//   location rule;
// - every counter-example is rejected by the grammar or reported by the checks of meaning, except
//   those whose fault takes the reference sequence (listed in the fixture with why);
// - written back, every example is the text it was read from, in the style it used (three- or
//   one-letter, Ter or *, a stated deleted sequence kept), and the canonical form reads back to the
//   same value;
// - messages: the furthest viable prefix and what may follow it, as morphe's Earley chart gives them;
// - random descriptions of every kind, written by the formatter, read back to themselves, unambiguous;
// - readings in plain words.

#include <genomics/nomenclature/hgvs.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace nm = genomics::nomenclature;
namespace mf = morphe::formal;

namespace {

// The number of derivations of a word from a symbol, saturated at 2: memoised on spans, ε-rules
// included; a derivation cycle (A ⇒+ A over one span) is reported.
struct Derivations {
  const mf::Grammar& g;
  std::vector<mf::Symbol> w;
  std::vector<std::vector<std::size_t>> by_lhs;
  std::vector<char> nullable;
  std::vector<std::vector<char>> first;
  std::unordered_map<std::uint64_t, int> memo, seqmemo;
  bool cycle = false;
  explicit Derivations(const mf::Grammar& gr) : g(gr) {
    by_lhs.resize(g.nonterminals.size());
    for (std::size_t r = 0; r < g.rules.size(); r++) by_lhs[g.index(g.rules[r].lhs[0])].push_back(r);
    nullable = g.nullable();
    first.assign(g.nonterminals.size(), std::vector<char>(g.terminals, 0));
    for (bool ch = true; ch;) {
      ch = false;
      for (const auto& r : g.rules) {
        const std::size_t a = g.index(r.lhs[0]);
        for (mf::Symbol s : r.rhs) {
          if (g.terminal(s)) { if (!first[a][s]) { first[a][s] = 1; ch = true; } break; }
          const std::size_t b = g.index(s);
          for (std::size_t t = 0; t < g.terminals; t++)
            if (first[b][t] && !first[a][t]) { first[a][t] = 1; ch = true; }
          if (!nullable[b]) break;
        }
      }
    }
  }
  static std::uint64_t key(std::uint64_t a, std::uint64_t b, std::uint64_t i, std::uint64_t j) { return (((a * 64 + b) * 4096 + i) * 4096) + j; }
  int count(mf::Symbol s, std::size_t i, std::size_t j) {
    if (g.terminal(s)) return (j == i + 1 && w[i] == s) ? 1 : 0;
    const std::size_t a = g.index(s);
    if (i == j && !nullable[a]) return 0;
    if (i < j && !first[a][w[i]]) return 0;
    const std::uint64_t k = key(a, 0, i, j);
    if (auto it = memo.find(k); it != memo.end()) {
      if (it->second < 0) { cycle = true; return 0; }
      return it->second;
    }
    memo[k] = -1;
    int total = 0;
    for (std::size_t r : by_lhs[a]) {
      total += seq(r, 0, i, j);
      if (total >= 2) { total = 2; break; }
    }
    memo[k] = total;
    return total;
  }
  int seq(std::size_t rule, std::size_t d, std::size_t i, std::size_t j) {
    const mf::Form& rhs = g.rules[rule].rhs;
    if (d == rhs.size()) return i == j ? 1 : 0;
    const std::uint64_t k = key(rule, d, i, j);
    if (auto it = seqmemo.find(k); it != seqmemo.end()) return it->second;
    int total = 0;
    for (std::size_t m = i; m <= j && total < 2; m++) {
      const int rest = seq(rule, d + 1, m, j);
      if (rest) total += count(rhs[d], i, m) * rest;
    }
    total = std::min(total, 2);
    seqmemo[k] = total;
    return total;
  }
};

int derivations(const mf::Bnf& bnf, const std::string& text, bool& cycle) {
  Derivations c(bnf.grammar());
  for (char ch : text) c.w.push_back(bnf.atom(static_cast<unsigned char>(ch)));
  const int n = c.count(bnf.grammar().start(), 0, c.w.size());
  cycle = c.cycle;
  return n;
}

struct Example {
  std::string text, status, page, reason;
};

std::vector<Example> examples() {
  std::ifstream in(std::string(GENOMICS_FIXTURES) + "/hgvs-examples.tsv");
  std::vector<Example> out;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::stringstream ss(line);
    Example e;
    std::string number;
    std::getline(ss, e.text, '\t');
    std::getline(ss, e.status, '\t');
    std::getline(ss, e.page, '\t');
    std::getline(ss, number, '\t');
    std::getline(ss, e.reason, '\t');
    out.push_back(e);
  }
  return out;
}

bool has_error(const nm::HgvsParse& p) {
  for (const auto& d : p.diagnostics)
    if (d.severity == nm::Diagnostic::Severity::error) return true;
  return false;
}

// ---- random descriptions -------------------------------------------------------------------------------
struct Gen {
  calx::random::Stream& r;
  std::size_t below(std::size_t n) { return static_cast<std::size_t>(r.below(n)); }
  std::string seq(bool rna, std::size_t n) {
    std::string s;
    for (std::size_t i = 0; i < n; i++) s += (rna ? "acgu" : "ACGT")[below(4)];
    return s;
  }
  nm::Position pos(nm::Coordinate t, std::int64_t at) {
    nm::Position p;
    p.number = at;
    if (t == nm::Coordinate::c || t == nm::Coordinate::r) {
      const auto k = below(6);
      if (k == 0) p.utr = '-';
      if (k == 1) p.utr = '*';
      if (k == 2) { p.offset_sign = '+'; p.offset = 1 + static_cast<std::int64_t>(below(200)); }
      if (k == 3) { p.offset_sign = '-'; p.offset = 1 + static_cast<std::int64_t>(below(200)); }
    }
    if (t == nm::Coordinate::n && below(3) == 0) { p.offset_sign = '+'; p.offset = 1 + static_cast<std::int64_t>(below(50)); }
    return p;
  }
  nm::Variant nucleotide(nm::Coordinate t) {
    const bool rna = t == nm::Coordinate::r;
    nm::Variant v;
    const std::int64_t a = 1 + static_cast<std::int64_t>(below(5000));
    v.location.start.first = v.location.start.last = pos(t, a);
    v.location.end = v.location.start;
    const auto range = [&] {
      v.location.range = true;
      v.location.end.first = v.location.end.last = pos(t, a + 1 + static_cast<std::int64_t>(below(30)));
      if (below(5) == 0) {  // an uncertain start
        v.location.start.uncertain = v.location.start.span = true;
        v.location.start.last = v.location.start.first;
        v.location.start.first.number = std::max<std::int64_t>(1, a - 10);
      }
    };
    switch (below(8)) {
      case 0: v.kind = nm::Kind::substitution; v.was = seq(rna, 1); v.now = seq(rna, 1); if (v.now == v.was) v.now = rna ? (v.was == "a" ? "c" : "a") : (v.was == "A" ? "C" : "A"); break;
      case 1: v.kind = nm::Kind::deletion; if (below(2)) range(); break;
      case 2: v.kind = nm::Kind::duplication; if (below(2)) range(); break;
      case 3: v.kind = nm::Kind::insertion; range(); v.inserted.push_back({nm::Piece::Kind::sequence, seq(rna, 1 + below(6)), {}, {}, false, {}, t}); break;
      case 4: v.kind = nm::Kind::delins; if (below(2)) range(); v.inserted.push_back({nm::Piece::Kind::sequence, seq(rna, 1 + below(6)), {}, {}, false, {}, t}); break;
      case 5: v.kind = nm::Kind::inversion; range(); break;
      case 6: v.kind = nm::Kind::repeat; v.units.push_back({seq(rna, 3), {nm::Copies::Kind::exact, 2 + static_cast<std::int64_t>(below(30)), 0, false}}); v.units.back().copies.high = v.units.back().copies.low; break;
      default: v.kind = nm::Kind::identity; if (below(2)) range(); break;
    }
    return v;
  }
  nm::Variant protein() {
    nm::Variant v;
    v.protein = true;
    static const std::string aa = "ARNDCQEGHILKMFPSTWYV";
    const auto res = [&](std::int64_t n) { return nm::Residue{aa[below(aa.size())], n}; };
    const std::int64_t a = 2 + static_cast<std::int64_t>(below(800));
    v.residues.start.first = v.residues.start.last = res(a);
    v.residues.end = v.residues.start;
    const auto range = [&](std::int64_t len) {
      v.residues.range = true;
      v.residues.end.first = v.residues.end.last = res(a + len);
    };
    switch (below(7)) {
      case 0: v.kind = nm::Kind::substitution; v.now = std::string(1, aa[below(aa.size())]); break;
      case 1: v.kind = nm::Kind::deletion; if (below(2)) range(1 + static_cast<std::int64_t>(below(5))); break;
      case 2: v.kind = nm::Kind::duplication; if (below(2)) range(1 + static_cast<std::int64_t>(below(5))); break;
      case 3: v.kind = nm::Kind::insertion; range(1); v.residues_inserted = std::string(1, aa[below(aa.size())]) + aa[below(aa.size())]; break;
      case 4: v.kind = nm::Kind::frameshift; v.fs_new = std::string(1, aa[below(aa.size())]); v.has_stop = true; v.stop_at = 2 + static_cast<std::int64_t>(below(60)); break;
      case 5: v.kind = nm::Kind::extension; v.residues.start.first.residue = v.residues.start.last.residue = v.residues.end.first.residue = v.residues.end.last.residue = '*'; v.now = "Q"; v.has_stop = true; v.stop_at = 1 + static_cast<std::int64_t>(below(40)); break;
      default: v.kind = nm::Kind::delins; range(1); v.residues_inserted = "WV"; break;
    }
    return v;
  }
};

}  // namespace

int main() {
  calx::verification::Suite s;
  const auto& bnf = nm::hgvs_grammar();
  s.require("the grammar compiles to morphe's grammar", bnf.grammar().rules.size() > 300,
            std::to_string(bnf.rule_names().size()) + " rules written, " + std::to_string(bnf.grammar().rules.size()) + " compiled, " + std::to_string(bnf.atoms()) + " atoms of bytes");

  const auto ex = examples();
  std::size_t counts[5] = {0, 0, 0, 0, 0};
  std::string rejected, ambiguous, faulted, not_caught, rewritten, position_rejected, reason_missing;
  std::size_t written_back = 0, valid = 0;
  for (const auto& e : ex) {
    if (e.status == "valid") {
      counts[0]++;
      valid++;
      const auto p = nm::parse_hgvs(e.text);
      if (!p.accepted) { rejected += e.text + " (" + p.error + "); "; continue; }
      bool cycle = false;
      if (derivations(bnf, e.text, cycle) != 1 || cycle) ambiguous += e.text + "; ";
      if (has_error(p)) faulted += e.text + "; ";
      // Written back in some style, it is the text; the canonical form reads back to the same value.
      bool same = false;
      for (int k = 0; k < 8 && !same; k++) same = nm::format(p.value, nm::Style{(k & 1) != 0, (k & 2) != 0, (k & 4) != 0}) == e.text;
      written_back += same;
      if (!same) rewritten += e.text + " -> " + p.canonical + "; ";
      const auto q = nm::parse_hgvs(p.canonical);
      if (!q.accepted || nm::format(q.value) != p.canonical) rewritten += "[canonical " + p.canonical + " does not read back]; ";
    } else if (e.status == "position") {
      counts[1]++;
      std::string t = e.text;
      if (const auto c = t.find(':'); c != std::string::npos) t = t.substr(c + 1);
      const char k = t[0] == 'm' || t[0] == 'o' ? 'g' : t[0];
      if (!bnf.parse(t.substr(2), std::string(1, k) + "_loc").accepted) position_rejected += e.text + "; ";
    } else if (e.status == "invalid") {
      counts[2]++;
      const auto p = nm::parse_hgvs(e.text);
      if (p.accepted && p.diagnostics.empty()) not_caught += e.text + "; ";
    } else if (e.status == "in-context") {
      counts[3]++;
      if (e.reason.empty()) reason_missing += e.text + "; ";
    } else {
      counts[4]++;
      if (e.reason.empty()) reason_missing += e.text + "; ";
    }
  }
  s.require("the fixture holds the standard's examples", ex.size() > 1000,
            std::to_string(counts[0]) + " descriptions, " + std::to_string(counts[1]) + " locations, " + std::to_string(counts[2]) + " counter-examples, " +
                std::to_string(counts[3]) + " counter-examples in context, " + std::to_string(counts[4]) + " excluded with reasons");
  s.require("every description the standard gives is in the language", rejected.empty(), rejected);
  s.require("each has exactly one derivation (counted independently of Earley's parser)", ambiguous.empty(), ambiguous);
  s.require("none breaks a condition of meaning", faulted.empty(), faulted);
  s.require("written back in its own style, each is the text it was read from; the canonical form reads back", rewritten.empty() && written_back == valid,
            std::to_string(written_back) + " of " + std::to_string(valid) + (rewritten.empty() ? "" : ": " + rewritten));
  s.require("every location the standard writes alone is read by its family's location rule", position_rejected.empty(), position_rejected);
  s.require("every counter-example is rejected by the grammar or reported by the checks of meaning", not_caught.empty(), not_caught);
  s.require("every counter-example that needs the reference, and every exclusion, says why", reason_missing.empty(), reason_missing);

  // ---- messages -----------------------------------------------------------------------------------------------
  {
    const auto p = nm::parse_hgvs("NM_004006.2:c.4375");
    const auto has = [&](const std::string& x) { return std::find(p.expected.begin(), p.expected.end(), x) != p.expected.end(); };
    s.require("NM_004006.2:c.4375: incomplete after the position, an edit or more of the position expected",
              !p.accepted && p.furthest == 18 && has("'>'") == false && has("'del'") && has("'dup'") && has("'ins'") && has("'inv'") && has("'_'") && has("[+\\-]"), p.error);
    const auto q = nm::parse_hgvs("NM_004006.2:c.4375C>");
    s.require("NM_004006.2:c.4375C>: a base expected after '>'", !q.accepted && q.furthest == 20 && q.error.find("[ACGTBDHKMNRSVWY]") != std::string::npos, q.error);
    const auto r = nm::parse_hgvs("NM_004006.2:c.4375dex");
    s.require("NM_004006.2:c.4375dex: at column 21, 'l' (of 'del') expected", !r.accepted && r.furthest == 20 && r.error.find("at column 21") == 0 && r.error.find("'l' (of 'del')") != std::string::npos, r.error);
    const auto t = nm::parse_hgvs("g.123insG");
    s.require("g.123insG reads, and the check of meaning names its fault", t.accepted && has_error(t) && t.diagnostics[0].message.find("two flanking positions") != std::string::npos);
    const auto u = nm::parse_hgvs("p.Trp24Cys ");
    s.require("a trailing space is not part of a description (general.md: spaces are not permitted)", !u.accepted && u.furthest == 10, u.error);
  }

  // ---- random descriptions -----------------------------------------------------------------------------------
  {
    calx::random::Stream r(20261101, 0x4e4f4d45, 1, 0);
    Gen gen{r};
    std::size_t ok = 0, unambiguous = 0, n = 0;
    std::string first_bad;
    const nm::Coordinate types[] = {nm::Coordinate::g, nm::Coordinate::m, nm::Coordinate::c, nm::Coordinate::n, nm::Coordinate::r, nm::Coordinate::p};
    for (int i = 0; i < 600; i++) {
      nm::Hgvs h;
      nm::Description d;
      d.has_reference = gen.below(3) != 0;
      d.reference.accession = "NM_00" + std::to_string(1000 + gen.below(9000));
      d.reference.version = std::to_string(1 + gen.below(9));
      d.variants.coordinate = types[gen.below(6)];
      nm::Alternative alt;
      nm::Group g;
      const auto variant = [&] { return d.variants.coordinate == nm::Coordinate::p ? gen.protein() : gen.nucleotide(d.variants.coordinate); };
      const auto shape = gen.below(4);
      if (shape == 0) {
        g.chain.items.push_back({{variant()}, gen.below(2) != 0, 0});
      } else if (shape == 1) {  // in cis
        g.bracketed = true;
        nm::Allele a;
        a.chains.push_back({{{{variant()}, false, 0}}});
        a.chains.push_back({{{{variant()}, false, 0}}});
        a.separators = ";";
        g.alleles.push_back(a);
      } else if (shape == 2) {  // in trans
        g.bracketed = true;
        for (int k = 0; k < 2; k++) {
          nm::Allele a;
          a.chains.push_back({{{{variant()}, gen.below(2) != 0, 0}}});
          g.alleles.push_back(a);
        }
      } else {  // predicted alternatives
        g.chain.items.push_back({{variant(), variant()}, true, 0});
      }
      alt.groups.push_back(g);
      if (gen.below(4) == 0) {
        nm::Group g2;
        g2.chain.items.push_back({{variant()}, false, 0});
        alt.groups.push_back(g2);
      }
      d.variants.alternatives.push_back(alt);
      h.alternatives.push_back(d);
      const std::string text = nm::format(h);
      const auto p = nm::parse_hgvs(text);
      n++;
      const bool same = p.accepted && nm::format(p.value) == text;
      ok += same;
      bool cycle = false;
      unambiguous += p.accepted && derivations(bnf, text, cycle) == 1 && !cycle;
      if (!same && first_bad.empty()) first_bad = text + (p.accepted ? " -> " + nm::format(p.value) : " : " + p.error);
    }
    s.require("600 random descriptions of every coordinate type and kind read back to themselves", ok == n, std::to_string(ok) + " of " + std::to_string(n) + (first_bad.empty() ? "" : "; " + first_bad));
    s.require("each with one derivation", unambiguous == n, std::to_string(unambiguous) + " of " + std::to_string(n));
  }

  // ---- the value, the canonical form, readings ---------------------------------------------------------------------
  {
    const auto p = nm::parse_hgvs("NG_012232.1(NM_004006.2):c.93+1G>T");
    const auto& d = p.value.alternatives.at(0);
    const auto& v = d.variants.alternatives.at(0).groups.at(0).chain.items.at(0).choices.at(0);
    s.require("NG_012232.1(NM_004006.2):c.93+1G>T: reference, selector, c., position 93 offset +1, G>T",
              d.reference.accession == "NG_012232" && d.reference.version == "1" && d.reference.selector() == "NM_004006.2" && d.variants.coordinate == nm::Coordinate::c &&
                  v.kind == nm::Kind::substitution && v.location.start.first.number == 93 && v.location.start.first.offset_sign == '+' && v.location.start.first.offset == 1 && v.was == "G" && v.now == "T");
    const auto q = nm::parse_hgvs("NP_003997.1:p.W24*");
    s.require("one-letter codes and * are read, the canonical form is three-letter with Ter", q.canonical == "NP_003997.1:p.Trp24Ter", q.canonical);
    const auto a = nm::parse_hgvs("NC_000023.11:g.33344590_33344592delGAT");
    s.require("a stated deleted sequence is dropped from the canonical form, and noted", a.canonical == "NC_000023.11:g.33344590_33344592del" && !a.diagnostics.empty(), a.canonical);
    const auto b = nm::parse_hgvs("NC_000003.12:g.(63912602_63912844)delN[15]");
    s.require("an uncertain deletion keeps its N[15]", b.canonical == "NC_000003.12:g.(63912602_63912844)delN[15]", b.canonical);
    const std::string r1 = nm::reading(nm::parse_hgvs("NM_000492.4:c.1521_1523del").value);
    s.require("reading: NM_000492.4:c.1521_1523del", r1 == "a deletion of 3 bases, c.1521 to c.1523 of NM_000492.4", r1);
    const std::string r2 = nm::reading(nm::parse_hgvs("NP_0123456.1:p.(Arg97ProfsTer23)").value);
    s.require("reading: a predicted frameshift", r2.find("a frameshift from Arg97, which becomes Pro, with a new stop at codon 23 of the new frame") == 0 && r2.find("(predicted)") != std::string::npos, r2);
    const std::string r3 = nm::reading(nm::parse_hgvs("NC_000023.10:g.(?_32238146)_(32984039_?)del").value);
    s.require("reading: a deletion with unknown ends", r3 == "a deletion of somewhere before g.32238146 to somewhere after g.32984039 of NC_000023.10", r3);
  }
  return s.done();
}
