// Checks of describe.hpp and spdi.hpp: a change's class of equivalent representations and the three
// canonical choices from it, held to brute force; the HGVS description by the standard's rules, held
// to the standard's and SPDI's examples and to Mutalyzer 3's description extractor.
//
// - Theorem 1 (theory/nomenclature.md §3), by brute force on random references with tandem repeats:
//   the minimal representations of a change's haplotype (every window of the reference whose
//   replacement makes it, no base shared at the two sides' ends) are one window, or, for a pure
//   deletion or insertion, the placements p of the interval [|h| − lcs, lcp] (deletions) or
//   [|r| − lcs, lcp] (insertions), lcp and lcs the longest common prefix and suffix of r and h.
// - Theorem 2: VCF's normal form (Tan et al., variant/normalize.hpp) is the interval's left end with
//   its padding base, HGVS's 3′ form its right end, SPDI's VOCA (as Holmes et al. state it, voca())
//   the hull; every description written makes the haplotype again (change_of).
// - The standard's examples: the 3′ rule in an A-stretch (DNA/deletion.md), CFTR's ATCTTTGGT, a
//   duplication not an insertion, an inversion, the insertion of a copy g.17_18ins5_16
//   (DNA/insertion.md), repeats; SPDI's Table 2 (ATGAGACT: 2:GAGA:GA, n.5_6del; ATGACT: 2:GA:GAGA,
//   n.3_4dup).
// - Mutalyzer 3's description_extract on 95 synthetic sequences (tests/fixtures/nomenclature-extractor.tsv):
//   where it describes the change as one variant, the description here is the same text.

#include <genomics/nomenclature/describe.hpp>
#include <genomics/nomenclature/spdi.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace nm = genomics::nomenclature;
using genomics::sequence::Codes;

namespace {

std::string letters(const Codes& x) { return genomics::sequence::decode_dna(x); }

struct Rep {
  std::size_t p, d;  // the window [p, p + d) of the reference
  std::string ins;
  bool operator<(const Rep& o) const { return std::tie(p, d, ins) < std::tie(o.p, o.d, o.ins); }
  bool operator==(const Rep& o) const { return std::tie(p, d, ins) == std::tie(o.p, o.d, o.ins); }
};

// Every minimal representation of h against r.
std::vector<Rep> minimal_representations(const Codes& r, const Codes& h) {
  std::vector<Rep> out;
  const std::size_t n = r.size(), m = h.size();
  for (std::size_t p = 0; p <= n; p++)
    for (std::size_t e = p; e <= n; e++) {
      const std::size_t suffix = n - e;
      if (m < p + suffix) continue;
      if (!std::equal(r.begin(), r.begin() + static_cast<std::ptrdiff_t>(p), h.begin())) continue;
      if (!std::equal(r.end() - static_cast<std::ptrdiff_t>(suffix), r.end(), h.end() - static_cast<std::ptrdiff_t>(suffix))) continue;
      const std::string D = letters(Codes(r.begin() + static_cast<std::ptrdiff_t>(p), r.begin() + static_cast<std::ptrdiff_t>(e)));
      const std::string I = letters(Codes(h.begin() + static_cast<std::ptrdiff_t>(p), h.end() - static_cast<std::ptrdiff_t>(suffix)));
      if (D == I) continue;
      if (!D.empty() && !I.empty() && (D.front() == I.front() || D.back() == I.back())) continue;
      out.push_back({p, e - p, I});
    }
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- Theorems 1 and 2 by brute force ---------------------------------------------------------------------
  {
    calx::random::Stream r(20261101, 0x44455343, 0, 0);
    std::size_t cases = 0, indels = 0, interval_ok = 0, closed_ok = 0, vcf_ok = 0, hgvs_ok = 0, spdi_ok = 0, voca_ok = 0, meaning_ok = 0, total_reps = 0;
    std::string first_bad;
    for (int k = 0; k < 600; k++) {
      // A reference: random flanks around a tandem repeat of a short unit.
      std::string unit, ref;
      for (std::size_t i = 0, u = 1 + r.below(4); i < u; i++) unit += "ACGT"[r.below(4)];
      for (std::size_t i = 0, f = 3 + r.below(6); i < f; i++) ref += "ACGT"[r.below(4)];
      for (std::size_t i = 0, c = 1 + r.below(5); i < c; i++) ref += unit;
      for (std::size_t i = 0, f = 3 + r.below(6); i < f; i++) ref += "ACGT"[r.below(4)];
      const Codes rc = genomics::sequence::encode_dna(ref);
      // A change: a deletion, an insertion (often of the unit), a substitution or a delins.
      nm::Change c;
      c.position = r.below(ref.size());
      const auto kind = r.below(4);
      if (kind == 0) c.deleted = ref.substr(c.position, std::min<std::size_t>(1 + r.below(4), ref.size() - c.position));
      else if (kind == 1) { c.inserted = r.below(2) ? unit : std::string(1, "ACGT"[r.below(4)]); if (r.below(3) == 0) c.inserted += c.inserted; }
      else if (kind == 2) { c.deleted = ref.substr(c.position, 1); c.inserted = std::string(1, "ACGT"[(genomics::sequence::dna_code(ref[c.position]) + 1 + r.below(3)) % 4]); }
      else { c.deleted = ref.substr(c.position, std::min<std::size_t>(2 + r.below(3), ref.size() - c.position)); for (std::size_t i = 0, L = 1 + r.below(4); i < L; i++) c.inserted += "ACGT"[r.below(4)]; }
      const Codes h = nm::apply_change(rc, c);
      if (h == rc) continue;
      cases++;
      const auto reps = minimal_representations(rc, h);
      total_reps += reps.size();
      const nm::Placements pl = nm::placements(rc, c);
      // Theorem 1: the placements.
      std::set<Rep> want(reps.begin(), reps.end()), have;
      if (pl.indel) {
        indels++;
        for (std::size_t p = pl.first; p <= pl.last; p++) {
          if (pl.deletion) have.insert({p, pl.length, ""});
          else have.insert({p, 0, letters(Codes(h.begin() + static_cast<std::ptrdiff_t>(p), h.begin() + static_cast<std::ptrdiff_t>(p + pl.length)))});
        }
        // The closed form by lcp and lcs.
        std::size_t lcp = 0, lcs = 0;
        while (lcp < rc.size() && lcp < h.size() && rc[lcp] == h[lcp]) lcp++;
        while (lcs < rc.size() && lcs < h.size() && rc[rc.size() - 1 - lcs] == h[h.size() - 1 - lcs]) lcs++;
        const std::size_t lo = pl.deletion ? h.size() - lcs : rc.size() - lcs, hi = std::min(lcp, pl.deletion ? h.size() : rc.size());
        closed_ok += lo == pl.first && hi == pl.last;
      } else {
        have.insert({pl.minimal.position, pl.minimal.deleted.size(), pl.minimal.inserted});
        closed_ok++;
      }
      const bool interval = want == have;
      interval_ok += interval;
      // Theorem 2: VCF's left end, HGVS's right end, SPDI's hull.
      const auto vcf = nm::vcf_form(rc, c);
      const nm::Change left = pl.at_first();
      std::string vref = left.deleted, valt = left.inserted;
      std::size_t vpos = left.position;
      if (vref.empty() || valt.empty()) {
        if (vpos > 0) { vref.insert(vref.begin(), ref[vpos - 1]); valt.insert(valt.begin(), ref[vpos - 1]); vpos--; }
        else { const char after = ref[left.deleted.size()]; vref += after; valt += after; }
      }
      const bool v_ok = vcf.position == vpos && vcf.ref == vref && vcf.alts.at(0) == valt;
      vcf_ok += v_ok;
      const nm::Described d = nm::describe(rc, c);
      const bool h_ok = d.change.position == pl.at_last().position || d.variant.kind == nm::Kind::delins;  // a delins at an end has no flank
      hgvs_ok += h_ok;
      const nm::Spdi sp = nm::spdi_form(rc, c);
      std::size_t hb = pl.indel ? pl.first : pl.minimal.position;
      std::size_t he = pl.indel ? pl.last + (pl.deletion ? pl.length : 0) : pl.minimal.position + pl.minimal.deleted.size();
      const bool s_ok = sp.position == hb && sp.deleted == ref.substr(hb, he - hb) && nm::apply_change(rc, nm::change_of(sp)) == h;
      spdi_ok += s_ok;
      const nm::Spdi vo = nm::voca(rc, nm::Spdi{"", c.position, c.deleted, c.inserted});
      voca_ok += vo == sp;
      // The description makes the haplotype again.
      const auto back = nm::change_of(rc, d.variant);
      const bool m_ok = back && nm::apply_change(rc, *back) == h;
      meaning_ok += m_ok;
      if ((!interval || !v_ok || !h_ok || !s_ok || !(vo == sp) || !m_ok) && first_bad.empty())
        first_bad = ref + " @" + std::to_string(c.position) + " " + c.deleted + ">" + c.inserted + " → " + d.text() + " SPDI " + sp.text() + " VOCA " + vo.text();
    }
    s.require("Theorem 1: the minimal representations are one window or an interval of placements, by brute force", interval_ok == cases,
              std::to_string(interval_ok) + " of " + std::to_string(cases) + " changes (" + std::to_string(indels) + " indels, " + std::to_string(total_reps) + " representations)" + (first_bad.empty() ? "" : "; " + first_bad));
    s.require("the interval's closed form [|h| − lcs, lcp] (deletions), [|r| − lcs, lcp] (insertions)", closed_ok == cases, std::to_string(closed_ok) + " of " + std::to_string(cases));
    s.require("Theorem 2: VCF's normal form (Tan et al.) is the left end, with its padding base", vcf_ok == cases, std::to_string(vcf_ok) + " of " + std::to_string(cases));
    s.require("HGVS's 3′ rule is the right end", hgvs_ok == cases, std::to_string(hgvs_ok) + " of " + std::to_string(cases));
    s.require("SPDI's contextual allele is the hull, and makes the haplotype", spdi_ok == cases, std::to_string(spdi_ok) + " of " + std::to_string(cases));
    s.require("VOCA as Holmes et al. state it (shrink, bud, blossom) gives the hull", voca_ok == cases, std::to_string(voca_ok) + " of " + std::to_string(cases));
    s.require("every description written makes the haplotype again", meaning_ok == cases, std::to_string(meaning_ok) + " of " + std::to_string(cases));
  }

  // ---- the standard's and SPDI's examples ----------------------------------------------------------------------------------
  {
    const auto desc = [](const std::string& ref, std::size_t p, const std::string& del, const std::string& ins, nm::DescribeOptions o = {}) {
      return nm::describe(genomics::sequence::encode_dna(ref), nm::Change{p, del, ins}, o).text();
    };
    // DNA/deletion.md: ATTGAAAAAAAATTAG, one A of the 8-A stretch c.5690..c.5697 deleted: the last one.
    s.require("the 3′ rule: an A deleted from an A-stretch is its last (DNA/deletion.md, c.5697del)", desc("ATTGAAAAAAAATTAG", 5, "A", "") == "12del", desc("ATTGAAAAAAAATTAG", 5, "A", ""));
    s.require("one A inserted there is a duplication of the last (c.5697dup)", desc("ATTGAAAAAAAATTAG", 7, "", "A") == "12dup");
    s.require("two As inserted are a repeated sequence over the stretch", desc("ATTGAAAAAAAATTAG", 7, "", "AA") == "5_12A[10]", desc("ATTGAAAAAAAATTAG", 7, "", "AA"));
    // CFTR's ATCTTTGGT: TCT, CTT and TTT deleted at c.1520_1522, c.1521_1523, c.1522_1524 (DNA/deletion.md).
    s.require("CFTR: CTT and TCT deleted from ATCTTTGGT are one class, written c.1521_1523del's place", desc("ATCTTTGGT", 1, "TCT", "") == desc("ATCTTTGGT", 2, "CTT", "") && desc("ATCTTTGGT", 2, "CTT", "") == "3_5del",
              desc("ATCTTTGGT", 1, "TCT", ""));
    s.require("and TTT deleted is another (c.1522_1524del)", desc("ATCTTTGGT", 3, "TTT", "") == "4_6del");
    s.require("an inserted copy of the sequence before it is a duplication (DNA/duplication.md, c.20dup)", desc("AGAAGTAGAGG", 6, "", "T") == "6dup");
    s.require("TCAG replaced by CTGA is an inversion (DNA/inversion.md)", desc("CATCAGCCT", 2, "TCAG", "CTGA") == "3_6inv");
    s.require("one base by one is a substitution, two or more a delins", desc("ACGT", 1, "C", "T") == "2C>T" && desc("ACGT", 1, "CG", "TA") == "2_3delinsTA");
    nm::DescribeOptions copies;
    copies.copy_length = 12;
    s.require("ATCGATCGATCGATCGAGGG with ATCGATCGATCG inserted after g.17: g.17_18ins5_16, the nearest copy (DNA/insertion.md)",
              desc("ATCGATCGATCGATCGAGGGTCCC", 17, "", "ATCGATCGATCG", copies) == "17_18ins5_16", desc("ATCGATCGATCGATCGAGGGTCCC", 17, "", "ATCGATCGATCG", copies));
    nm::DescribeOptions seven;
    seven.copy_length = 7;
    s.require("an inverted copy of the sequence before it is an insertion with inv (g.234_235ins123_234inv)", desc("TTTTACGGTCAGGG", 11, "", "TGACCGT", seven) == "11_12ins5_11inv",
              desc("TTTTACGGTCAGGG", 11, "", "TGACCGT", seven));
    // SPDI's Table 2 (Holmes et al. 2020).
    const Codes t2 = genomics::sequence::encode_dna("ATGAGACT"), t2i = genomics::sequence::encode_dna("ATGACT");
    const nm::Spdi a = nm::spdi_form(t2, {2, "GA", ""}), b = nm::spdi_form(t2, {4, "GA", ""});
    const nm::Spdi c = nm::spdi_form(t2i, {2, "", "GA"}), d = nm::spdi_form(t2i, {4, "", "GA"});
    s.require("SPDI Table 2: the first or the second GA of ATGAGACT deleted are 2:GAGA:GA, n.5_6del",
              a.position == 2 && a.deleted == "GAGA" && a.inserted == "GA" && a == b && desc("ATGAGACT", 2, "GA", "") == "5_6del", a.text());
    s.require("SPDI Table 2: GA inserted before or after the GA of ATGACT is 2:GA:GAGA, n.3_4dup",
              c.position == 2 && c.deleted == "GA" && c.inserted == "GAGA" && c == d && desc("ATGACT", 2, "", "GA") == "3_4dup", c.text());
    const auto parsed = nm::parse_spdi("NC_000001.11:12345:1:G");
    s.require("SPDI with a deletion length reads (Seq1:4:1:G)", parsed.ok && parsed.deletion_as_length && parsed.deleted_length == 1 && parsed.spdi.inserted == "G");
    s.require("a malformed SPDI is reported, not thrown", !nm::parse_spdi("NC_000001.11:x:A:G").ok && !nm::parse_spdi("a:1:A").ok);
  }

  // ---- Mutalyzer 3's description extractor -------------------------------------------------------------------------------
  {
    std::ifstream in(std::string(GENOMICS_FIXTURES) + "/nomenclature-extractor.tsv");
    std::string line, version;
    std::size_t single = 0, agree = 0, decomposed = 0;
    std::string disagree;
    while (std::getline(in, line)) {
      if (line.rfind("# retrieved", 0) == 0) version = line.substr(2);
      if (line.empty() || line[0] == '#') continue;
      std::stringstream ss(line);
      std::string ref, obs, theirs;
      std::getline(ss, ref, '\t');
      std::getline(ss, obs, '\t');
      std::getline(ss, theirs, '\t');
      if (!theirs.empty() && theirs.front() == '"') theirs = theirs.substr(1, theirs.size() - 2);
      if (theirs.empty() || theirs.front() == '[' || theirs.front() == '{') { decomposed++; continue; }
      // The change: the whole of the reference replaced by the observed, trimmed by describe.
      single++;
      const std::string ours = nm::describe(genomics::sequence::encode_dna(ref), nm::Change{0, ref, obs}).text();
      if (ours == theirs) agree++;
      else disagree += ref + ">" + obs + ": " + ours + " vs " + theirs + "; ";
    }
    if (single == 0) s.skip("the Mutalyzer extractor fixture is absent (tests/fixtures/nomenclature-extractor.tsv)");
    else
      s.require("Mutalyzer's description extractor: every change it describes as one variant is described the same", agree == single,
                std::to_string(agree) + " of " + std::to_string(single) + " (" + std::to_string(decomposed) + " it decomposes into several variants; " + version + ")" + (disagree.empty() ? "" : ": " + disagree));
  }
  return s.done();
}
