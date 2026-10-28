// Checks of genomics/grammar: the Klein group of the bases and its subgroups, characters and
// breakings; the strand groups on words and their canonical k-mers; hairpin grammars, their
// languages and symmetry, with morphe's equivariant parsers on them; the genetic code's symmetry
// on (C₂)⁶; and, when the data library holds the NCBI genomes, Escherichia coli's tRNA genes
// against random windows, and the symmetry of its k-mer statistics (Chargaff's second rule).

#include <genomics/grammar/channels.hpp>
#include <genomics/grammar/code_symmetry.hpp>
#include <genomics/grammar/nucleotides.hpp>
#include <genomics/grammar/strands.hpp>
#include <genomics/grammar/structures.hpp>
#include <genomics/sequence/genbank.hpp>

#include <morphe/formal/earley.hpp>
#include <morphe/formal/breaking.hpp>
#include <morphe/formal/equivariant.hpp>
#include <morphe/formal/inside_outside.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <utility>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace gg = genomics::grammar;
namespace gs = genomics::sequence;
namespace ma = morphe::algebra;
namespace mf = morphe::formal;

namespace {

std::vector<std::vector<mf::Symbol>> all_words(std::size_t n) {
  std::vector<std::vector<mf::Symbol>> out;
  std::vector<mf::Symbol> w(n, 0);
  for (;;) {
    out.push_back(w);
    std::size_t i = 0;
    while (i < n && ++w[i] == 4) w[i++] = 0;
    if (i == n) break;
  }
  return out;
}

std::uint64_t pow4(std::size_t e) { return std::uint64_t{1} << (2 * e); }

std::filesystem::path cache() {
  const char* env = std::getenv("SCIENCE_CACHE");
  return env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the Klein group of the bases ------------------------------------------------------------------
  {
    const auto& k = gg::klein();
    bool regular = true;
    for (gs::Code x = 0; x < 4; x++) regular = regular && k.bases.orbit(x).size() == 4 && k.bases.stabiliser(x).size() == 1;
    s.require("the Klein group acts on the bases freely and transitively (a torsor)", regular);
    s.require("κ is Watson–Crick complementarity: A ↔ T, C ↔ G",
              k.bases.act(k.complement, gs::A) == gs::T && k.bases.act(k.complement, gs::C) == gs::G);
    s.require("τ is transition: A ↔ G, C ↔ T; and κτ = σ: A ↔ C, G ↔ T",
              k.bases.act(k.transition, gs::A) == gs::G && k.bases.act(k.transition, gs::C) == gs::T &&
                  k.group.mul(k.complement, k.transition) == k.transversion && k.bases.act(k.transversion, gs::A) == gs::C);
    const auto cl = gg::classifications();
    std::string letters;
    for (const auto& c : cl) letters += std::string{c.letters[0], '/', c.letters[1], ' '};
    s.require("the cosets of the three subgroups of order 2 are IUPAC's R/Y, W/S, M/K", letters == "R/Y W/S M/K ", letters);
    // Every coset of every subgroup, as a mask; the three-letter codes are none of them.
    std::set<int> cosets;
    const std::vector<std::vector<std::uint32_t>> subgroups{{0}, {0, k.transition}, {0, k.complement}, {0, k.transversion}, {0, 1, 2, 3}};
    for (const auto& h : subgroups)
      for (gs::Code x = 0; x < 4; x++) {
        int m = 0;
        for (auto g : h) m |= 1 << k.bases.act(g, x);
        cosets.insert(m);
      }
    bool three = true;
    for (char c : std::string("BDHV")) three = three && !cosets.count(gs::iupac_mask(c));
    s.require("the cosets are the four bases, the six two-letter codes and N; B, D, H, V are none", cosets.size() == 11 && three && cosets.count(gs::iupac_mask('N')));
    bool group = true, flips = true;
    for (gs::Code x = 0; x < 4; x++) {
      group = group && gg::character_value(k.transition, gs::A, x) * gg::character_value(k.complement, gs::A, x) == gg::character_value(k.transversion, gs::A, x);
      // Moving the base point from A to T = κ·A multiplies each character by its value at κ.
      flips = flips && gg::character_value(k.transition, gs::T, x) == -gg::character_value(k.transition, gs::A, x) &&
              gg::character_value(k.complement, gs::T, x) == gg::character_value(k.complement, gs::A, x);
    }
    s.require("the classifications are characters: χ_RY χ_WS = χ_MK", group);
    s.require("a new base point changes a character's sign by its value there, not its partition", flips);
    const auto ry = gg::quotient({0, k.transition}), n = gg::quotient({0, 1, 2, 3}), bases = gg::quotient({0});
    s.require("a breaking chain: N, then R/Y, then the bases", n.to == 1 && ry.to == 2 && bases.to == 4 && ry(gs::A) == ry(gs::G) && ry(gs::A) != ry(gs::C));
  }

  // ---- words: canonical k-mers by Burnside ----------------------------------------------------------
  {
    const auto& rc = gg::reverse_complement_group();
    const auto& st = gg::strand_group();
    bool closed = true, brute = true;
    for (std::size_t k = 1; k <= 8; k++) {
      const auto want_rc = k % 2 ? pow4(k) / 2 : (pow4(k) + pow4(k / 2)) / 2;
      const auto want_st = (pow4(k) + pow4((k + 1) / 2) + (k % 2 ? 0 : pow4(k / 2))) / 4;
      closed = closed && rc.words.word_orbits(k) == morphe::arithmetic::Integer(static_cast<unsigned long long>(want_rc)) &&
               st.words.word_orbits(k) == morphe::arithmetic::Integer(static_cast<unsigned long long>(want_st));
      if (k <= 6) {
        std::set<std::vector<mf::Symbol>> c1, c2;
        for (const auto& w : all_words(k)) { c1.insert(rc.words.canonical(w)); c2.insert(st.words.canonical(w)); }
        brute = brute && c1.size() == want_rc && c2.size() == want_st;
      }
    }
    s.require("canonical k-mers under the reverse complement: (4^k + 4^{k/2})/2 or 4^k/2, by Burnside", closed);
    s.require("and under the strand group, (4^k + 4^⌈k/2⌉ + [k even] 4^{k/2})/4; both against enumeration", brute);
  }

  // ---- hairpin grammars -------------------------------------------------------------------------------
  {
    const gg::HairpinSpec wc{2, 1, 3, false}, wobble{2, 1, 3, true};
    const auto g = gg::hairpin_grammar(wc), gw = gg::hairpin_grammar(wobble);
    const mf::Earley e(g), ew(gw);
    bool lang = true;
    std::size_t members = 0;
    for (std::size_t n = 1; n <= 8; n++)
      for (const auto& w : all_words(n)) {
        const bool a = e.accepts(w), b = ew.accepts(w);
        lang = lang && a == gg::is_hairpin(w, wc) && b == gg::is_hairpin(w, wobble);
        members += a;
      }
    s.require("the hairpin grammars generate the hairpins, with and without wobble pairs (every word up to 8)", lang, std::to_string(members) + " Watson–Crick hairpins");
    const auto& st = gg::strand_group();
    const mf::GrammarAction strands{&st.words, {}};
    // The reversal subgroup alone, and the complement alone (no reversal).
    const auto z2 = ma::cyclic(2);
    const auto reversal = ma::GroupAction::from_images(z2, 4, {{0, 1, 2, 3}, {0, 1, 2, 3}}, {1, -1});
    const auto complement = ma::GroupAction::from_images(z2, 4, {{0, 1, 2, 3}, {3, 2, 1, 0}});
    const mf::GrammarAction by_reversal{&reversal, {}}, by_complement{&complement, {}};
    s.require("the Watson–Crick hairpin grammar is invariant under the whole strand group", mf::invariant(g, strands));
    s.require("with wobble pairs it keeps reversal and loses the complement", mf::invariant(gw, by_reversal) && !mf::invariant(gw, by_complement));
    mf::EquivariantEarley eq(g, by_complement);
    const mf::Earley plain(g);
    bool same = true;
    std::size_t eq_items = 0, plain_items = 0, reused = 0, spans = 0;
    calx::random::Stream rng(2026, 0, 0, 0);
    for (int t = 0; t < 400; t++) {
      std::vector<mf::Symbol> w(4 + rng.below(9));
      for (auto& x : w) x = static_cast<mf::Symbol>(rng.below(4));
      // Half of them hairpins by construction: a random stem closed by its reverse complement.
      if (t % 2) {
        const std::size_t stem = 2 + rng.below(3), loop = 1 + rng.below(3);
        w.assign(2 * stem + loop, 0);
        for (std::size_t i = 0; i < stem + loop; i++) w[i] = static_cast<mf::Symbol>(rng.below(4));
        for (std::size_t i = 0; i < stem; i++) w[w.size() - 1 - i] = 3 - w[i];
      }
      const auto r = eq.parse(w);
      const auto p = plain.parse(w);
      same = same && r.accepts == p.accepts && r.accepts == gg::is_hairpin(w, wc);
      eq_items += r.items;
      plain_items += p.items;
      const auto sr = mf::span_recognition(g, strands, w);
      same = same && sr.accepts == p.accepts;
      reused += sr.reused;
      spans += sr.spans;
    }
    s.require("on hairpins, the equivariant parser (complement) and span sharing (strand group) agree with Earley's", same,
              std::to_string(eq_items) + " equivariant items against " + std::to_string(plain_items) + "; " + std::to_string(reused) + " of " +
                  std::to_string(spans) + " spans shared by strand symmetry");
    const auto gc = gg::contains_hairpin_grammar({3, 2, 4, false});
    const mf::Earley ec(gc);
    bool contains = true;
    for (int t = 0; t < 600; t++) {
      std::vector<mf::Symbol> w(6 + rng.below(9));
      for (auto& x : w) x = static_cast<mf::Symbol>(rng.below(4));
      contains = contains && ec.accepts(w) == gg::contains_hairpin(w, {3, 2, 4, false});
    }
    s.require("a hairpin anywhere: the flanked grammar against search", contains);
  }

  // ---- the genetic code on (C₂)⁶ ----------------------------------------------------------------------------
  {
    const auto& k = gg::klein();
    const std::uint32_t tau3 = gg::at_position(k.transition, 2);
    const auto standard = gg::code_on_group(gs::genetic_code(1));
    const auto mito = gg::code_on_group(gs::genetic_code(2));
    const auto ps = gg::invariance_profile(standard);
    const auto fail = gg::failures(standard, tau3);
    std::string f;
    for (const auto& x : fail) f += x + " ";
    s.require("the standard code is invariant under no element but the identity", gg::stabiliser(standard) == std::vector<std::uint32_t>{0});
    s.require("the third-position transition keeps 60 codons of 64, failing at ATA/ATG and TGA/TGG", ps[tau3] == 60 && f == "ATA/ATG TGA/TGG ", f);
    s.require("the vertebrate mitochondrial code is invariant under it: its stabiliser is {1, τ₃}", gg::stabiliser(mito) == std::vector<std::uint32_t>{0, tau3});
    // The spectrum: Wiener–Khinchin exactly, and invariance as vanishing on the odd characters.
    bool wk = true;
    for (const auto& code : {standard, mito}) {
      const auto p = gg::code_power_spectrum(code);
      const auto prof = gg::invariance_profile(code);
      const auto wh = ma::walsh_hadamard(std::vector<long long>(prof.begin(), prof.end()));
      for (std::size_t u = 0; u < 64; u++) wk = wk && wh[u] == p[u];
    }
    s.require("the code's power spectrum is the Walsh–Hadamard transform of its invariance profile (exact, both codes)", wk);
    const auto pm = gg::code_power_spectrum(mito);
    bool vanish = true;
    for (std::uint32_t u = 0; u < 64; u++)
      if (__builtin_popcount(u & tau3) % 2) vanish = vanish && pm[u] == 0;
    s.require("the mitochondrial code's spectrum vanishes on every character odd under τ₃ (invariance, dually)", vanish);
    const auto p1 = gg::code_power_spectrum(standard);
    long long all = 0, first_two = 0, first = 0;
    for (std::uint32_t u = 1; u < 64; u++) {
      all += p1[u];
      if (gg::ignores_positions(u, {2})) first_two += p1[u];
      if (gg::ignores_positions(u, {1, 2})) first += p1[u];
    }
    char detail[200];
    std::snprintf(detail, sizeof detail, "of the standard code's non-constant spectrum, %.1f%% lies in characters of the first two bases, %.1f%% in those of the first alone",
                  100.0 * static_cast<double>(first_two) / static_cast<double>(all), 100.0 * static_cast<double>(first) / static_cast<double>(all));
    s.require("the code's spectrum lies mostly in characters that ignore the third base", 2 * first_two > all, detail);
  }

  // ---- real genomes --------------------------------------------------------------------------------------
  {
    const auto path = cache() / "ncbi-genomes" / "NC_000913.3.gb";
    if (!std::filesystem::exists(path)) {
      s.skip("Escherichia coli is not in the data cache: its tRNAs and k-mer symmetry are not checked");
    } else {
      std::ifstream in(path, std::ios::binary);
      std::stringstream text;
      text << in.rdbuf();
      const auto records = gs::read_genbank(text.str());
      const auto& rec = records.at(0);
      const auto genome = gs::encode_dna(rec.sequence);
      const gg::HairpinSpec spec{5, 3, 9, false};
      const auto g = gg::contains_hairpin_grammar(spec);
      const mf::Earley e(g);
      std::size_t trnas = 0, with = 0, random_with = 0;
      calx::random::Stream rng(7, 0, 0, 0);
      for (const auto& f : rec.features) {
        if (f.key != "tRNA") continue;
        const auto x = gs::extract(genome, f.location);
        if (x.size() < 60 || x.size() > 100) continue;
        std::vector<mf::Symbol> w(x.begin(), x.end());
        if (std::any_of(w.begin(), w.end(), [](mf::Symbol b) { return b > 3; })) continue;
        trnas++;
        with += e.accepts(w);
        // A window of the same length at a random place.
        const std::size_t at = rng.below(genome.size() - w.size());
        std::vector<mf::Symbol> r(genome.begin() + static_cast<std::ptrdiff_t>(at), genome.begin() + static_cast<std::ptrdiff_t>(at + w.size()));
        if (std::any_of(r.begin(), r.end(), [](mf::Symbol b) { return b > 3; })) continue;
        random_with += e.accepts(r);
      }
      const double a = static_cast<double>(with) / static_cast<double>(trnas), b = static_cast<double>(random_with) / static_cast<double>(trnas);
      char detail[160];
      std::snprintf(detail, sizeof detail, "%zu tRNAs: %.1f%% contain a stem of 5 Watson–Crick pairs closing a loop of 3–9; random windows %.1f%%", trnas, 100 * a, 100 * b);
      s.require("Escherichia coli's tRNA genes carry hairpins far more often than random windows", trnas >= 80 && a > b + 0.25, detail);
      const auto inv = gg::kmer_invariance(genome, 4, gg::base_and_strand_group());
      const auto& w = gg::base_and_strand_group();
      double rc = 1, others = 1;
      std::string list;
      for (const auto& v : inv) {
        if (v.element == w.reverse_complement) rc = v.total_variation;
        else if (v.element != 0) others = std::min(others, v.total_variation);
        char line[96];
        std::snprintf(line, sizeof line, "%s %.4f; ", v.name.c_str(), v.total_variation);
        list += line;
      }
      s.require("its 4-mer statistics keep the reverse complement alone among the 8 symmetries (Chargaff's second rule)", rc < 0.05 && others > 4 * rc, list);

      // The Wedderburn channels of the 4-mer frequencies under the group of order 8.
      double norm = 0;
      const auto channels = gg::kmer_channels(genome, 4, w, &norm);
      double sum = 0, odd_rc = 0, odd_complement = 0;
      for (const auto& c : channels) {
        sum += c.energy;
        if (c.signs.empty()) continue;
        if (c.signs[w.reverse_complement] < 0) odd_rc += c.fraction;
        if (c.signs[w.complement] < 0) odd_complement += c.fraction;
      }
      s.check("the eight channels' energies sum to ‖f‖²", std::fabs(sum - norm) / norm, 1e-12);
      char cd[200];
      std::snprintf(cd, sizeof cd, "%zu channels; odd under the reverse complement %.2e of ‖f‖², odd under the complement %.2e",
                    channels.size(), odd_rc, odd_complement);
      s.require("the channels odd under the reverse complement hold almost nothing; those odd under the complement do not", channels.size() == 8 && odd_rc < 1e-3 && odd_complement > 20 * odd_rc, cd);

      // Graduated breaking on real tRNAs: parse over R/Y (or W/S) first, then over the bases.
      const auto& k = gg::klein();
      const auto fine = gg::contains_hairpin_grammar({4, 3, 8, false});
      const auto bfine = mf::BinaryGrammar::from(fine);
      const std::vector<morphe::arithmetic::Integer> unit(fine.rules.size(), morphe::arithmetic::Integer(1));
      bool exact = true;
      std::size_t spans_ry = 0, kept_ry = 0, kept_ws = 0, used = 0;
      for (const auto& f : rec.features) {
        if (f.key != "tRNA" || used >= 10) continue;
        const auto x = gs::extract(genome, f.location);
        std::vector<mf::Symbol> wd(x.begin(), x.end());
        if (x.size() > 100 || std::any_of(wd.begin(), wd.end(), [](mf::Symbol b) { return b > 3; })) continue;
        used++;
        const auto full = mf::inside<mf::Counting>(bfine, unit, wd);
        for (const auto& [sub, kept] : {std::pair{std::vector<std::uint32_t>{0, k.transition}, &kept_ry}, std::pair{std::vector<std::uint32_t>{0, k.complement}, &kept_ws}}) {
          mf::CoarseToFine st;
          const auto mask = mf::coarse_mask(fine, gg::quotient(sub), wd, &st);
          const auto pruned = mf::inside<mf::Counting>(bfine, unit, wd, &mask);
          exact = exact && pruned.at(0, 0, wd.size()) == full.at(0, 0, wd.size());
          *kept += st.kept;
          if (kept == &kept_ry) spans_ry += st.spans;
        }
      }
      char gd[200];
      std::snprintf(gd, sizeof gd, "%zu tRNAs; spans kept open: %.1f%% by R/Y, %.1f%% by W/S", used,
                    100.0 * static_cast<double>(kept_ry) / static_cast<double>(spans_ry), 100.0 * static_cast<double>(kept_ws) / static_cast<double>(spans_ry));
      s.require("graduated breaking: the coarse pass masks the fine one exactly (every parse counted)", exact && used == 10, gd);
    }
  }
  return s.done();
}
