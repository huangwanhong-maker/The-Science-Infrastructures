// The predicted consequence of a genomic change for a transcript and its protein: the p. description
// by the HGVS Nomenclature's protein recommendations (version 21.1, protein/*.md), the predicted
// r. description, and the consequence terms of the Sequence Ontology with Ensembl VEP's order of
// severity (McLaren et al. 2016, Genome Biology 17, 122; Ensembl release 116, Constants.pm).
//
// - The protein is predicted by translating the variant transcript with the record's genetic code:
//   the reference mRNA (its exons joined), the change applied at its place along the transcript (the
//   3′-most placement in the transcript's direction, transcript.hpp), the observed protein read from
//   the start codon to the first stop in frame, into the 3′ UTR if need be. The description compares
//   the reference protein P (with its stop) and the observed Q, as protein/*.md asks ("protein variant
//   descriptions should be derived from comparing the variant protein sequence with the reference
//   protein sequence"):
//   · Q = P: synonymous, p.(Leu5=), the residue of the first codon changed (protein/substitution.md);
//   · the start codon changed to one that is not a start: p.(Met1?) ("the consequence … of a variant
//     affecting the translation initiation codon can not be predicted"; Mutalyzer 3 writes p.?);
//   · the first difference at P's stop: an extension, p.(Ter110GlnextTer17), the new stop counted
//     in the residues added after the one replacing the stop (Mutalyzer 3's count:
//     p.(Ter148TyrextTer21) for HBB c.444A>T; theory §6), or Ter? when there is none (protein/extension.md);
//   · a change of length not divisible by three before the stop: a frameshift from the first changed
//     residue, p.(Arg97ProfsTer23), the new stop counted from it as 1, or a nonsense substitution
//     when that residue is itself a stop, p.(Tyr4Ter) (protein/frameshift.md);
//   · otherwise in frame: the common prefix and suffix trimmed, a substitution, a deletion or an
//     insertion placed by the 3′ rule on the protein (the most C-terminal residue changed), a
//     duplication when the inserted residues repeat those before them, a delins otherwise
//     (protein/deletion.md, duplication.md, insertion.md, delins.md). When the change introduces a
//     stop the residues after it are not described (p.(Met3_His4insGlyTer)); how many reference
//     residues the new ones replace is the number of codons the DNA change touches, which is how the
//     standard's own examples count (p.(Pro578_Lys579delinsLeuTer) for c.1733_1735delinsTTT;
//     theory/nomenclature.md §6 notes that this needs the DNA, against the rule quoted above).
// - Not predicted: a change at a splice site or in its region (VEP's splice terms) or across an
//   exon's edge, p.? ("an effect … is expected, but … it is not possible to give a reliable
//   prediction", uncertain.md; Mutalyzer 3 gives no protein there); a deletion of the whole
//   transcript, p.0?. An intronic change away from the splice sites is p.(=), as Mutalyzer has it.
//   Changes in the UTRs and outside the transcript are p.(=); new upstream starts (p.Met1ext-5)
//   are not predicted.
// - Consequence terms follow VEP's predicates (Bio::EnsEMBL::Variation::Utils::VariationEffect and
//   BaseTranscriptVariationAllele::_intron_effects, release 116): splice donor and acceptor sites are
//   an intron's first and last two bases; the splice region is 1–3 bases of the exon or 3–8 of the
//   intron; the donor's 5th base, its region (+3 to +6) and the polypyrimidine tract (−3 to −17);
//   upstream and downstream within 5,000 bases. Terms are listed most severe first.
//
// The genetic code is the transcript's (/transl_table, NCBI's tables in sequence/genetic_code.hpp).

#pragma once

#include <genomics/nomenclature/transcript.hpp>
#include <genomics/sequence/genetic_code.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::nomenclature {

// ---- consequence terms ---------------------------------------------------------------------------------
struct ConsequenceTerm {
  std::string_view term, accession, impact;
  int rank;
};

// Ensembl VEP's consequences in their order of severity (Ensembl release 116,
// modules/Bio/EnsEMBL/Variation/Utils/Constants.pm, OVERLAP_CONSEQUENCES' 'rank'; the table of
// https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html, retrieved
// 2026-11-01). tests/fixtures/vep-consequences.tsv holds the same table as made from the source.
inline const std::vector<ConsequenceTerm>& consequence_terms() {
  static const std::vector<ConsequenceTerm> terms = {
      {"transcript_ablation", "SO:0001893", "HIGH", 1},
      {"splice_acceptor_variant", "SO:0001574", "HIGH", 2},
      {"splice_donor_variant", "SO:0001575", "HIGH", 3},
      {"stop_gained", "SO:0001587", "HIGH", 4},
      {"frameshift_variant", "SO:0001589", "HIGH", 5},
      {"stop_lost", "SO:0001578", "HIGH", 6},
      {"start_lost", "SO:0002012", "HIGH", 7},
      {"transcript_amplification", "SO:0001889", "HIGH", 8},
      {"feature_elongation", "SO:0001907", "HIGH", 9},
      {"feature_truncation", "SO:0001906", "HIGH", 10},
      {"inframe_insertion", "SO:0001821", "MODERATE", 11},
      {"inframe_deletion", "SO:0001822", "MODERATE", 12},
      {"missense_variant", "SO:0001583", "MODERATE", 13},
      {"protein_altering_variant", "SO:0001818", "MODERATE", 14},
      {"splice_donor_5th_base_variant", "SO:0001787", "LOW", 15},
      {"splice_region_variant", "SO:0001630", "LOW", 16},
      {"splice_donor_region_variant", "SO:0002170", "LOW", 17},
      {"splice_polypyrimidine_tract_variant", "SO:0002169", "LOW", 18},
      {"incomplete_terminal_codon_variant", "SO:0001626", "LOW", 19},
      {"start_retained_variant", "SO:0002019", "LOW", 20},
      {"stop_retained_variant", "SO:0001567", "LOW", 21},
      {"synonymous_variant", "SO:0001819", "LOW", 22},
      {"coding_sequence_variant", "SO:0001580", "MODIFIER", 23},
      {"mature_miRNA_variant", "SO:0001620", "MODIFIER", 24},
      {"5_prime_UTR_variant", "SO:0001623", "MODIFIER", 25},
      {"3_prime_UTR_variant", "SO:0001624", "MODIFIER", 26},
      {"non_coding_transcript_exon_variant", "SO:0001792", "MODIFIER", 27},
      {"intron_variant", "SO:0001627", "MODIFIER", 28},
      {"NMD_transcript_variant", "SO:0001621", "MODIFIER", 29},
      {"non_coding_transcript_variant", "SO:0001619", "MODIFIER", 30},
      {"coding_transcript_variant", "SO:0001968", "MODIFIER", 31},
      {"upstream_gene_variant", "SO:0001631", "MODIFIER", 32},
      {"downstream_gene_variant", "SO:0001632", "MODIFIER", 33},
      {"TFBS_ablation", "SO:0001895", "MODERATE", 34},  // MODERATE in Constants.pm, MODIFIER on the web page
      {"TFBS_amplification", "SO:0001892", "MODIFIER", 35},
      {"TF_binding_site_variant", "SO:0001782", "MODIFIER", 36},
      {"regulatory_region_ablation", "SO:0001894", "MODIFIER", 37},
      {"regulatory_region_amplification", "SO:0001891", "MODIFIER", 38},
      {"regulatory_region_variant", "SO:0001566", "MODIFIER", 39},
      {"intergenic_variant", "SO:0001628", "MODIFIER", 40},
      {"sequence_variant", "SO:0001060", "MODIFIER", 41},
  };
  return terms;
}

inline const ConsequenceTerm* consequence(std::string_view term) {
  for (const auto& t : consequence_terms())
    if (t.term == term) return &t;
  return nullptr;
}

inline void sort_by_severity(std::vector<std::string>& terms) {
  const auto rank = [](const std::string& t) { const ConsequenceTerm* c = consequence(t); return c ? c->rank : 1000; };
  std::sort(terms.begin(), terms.end(), [&](const std::string& a, const std::string& b) { return rank(a) < rank(b); });
  terms.erase(std::unique(terms.begin(), terms.end()), terms.end());
}

// ---- prediction ------------------------------------------------------------------------------------------------
struct ProteinPrediction {
  std::string transcript, protein;     // NM_000518.5, NP_000509.1
  Projection projection;               // the c. or n. description
  bool described = false;              // a p. description was made
  std::string whole;                   // "=", "?", "0?" when the description is a whole-protein statement
  Variant variant;                     // otherwise the p. variant (predicted)
  std::string reference, observed;     // the proteins, one-letter, '*' for the stop
  std::string rna;                     // r.(…) predicted, when the change is exonic or outside the transcript
  std::vector<std::string> consequences;  // SO terms, most severe first
  std::string note;                    // why the protein is not described, or how it was
  std::string text(const Style& st = {}) const {
    if (!described) return std::string();
    const std::string id = protein.empty() ? std::string() : protein + ":";
    return id + "p." + (whole.empty() ? "(" + format_variant(variant, st) + ")" : whole);
  }
  std::string impact() const {
    if (consequences.empty()) return "MODIFIER";
    const ConsequenceTerm* c = consequence(consequences.front());
    return c ? std::string(c->impact) : "MODIFIER";
  }
};

namespace detail {

inline std::string translate_from(const Codes& m, std::size_t from, const sequence::GeneticCode& code) {
  std::string q;
  for (std::size_t i = from; i + 3 <= m.size(); i += 3) {
    const char a = sequence::amino_acid(code, m[i], m[i + 1], m[i + 2]);
    q += a;
    if (a == '*') break;
  }
  return q;
}

inline Residue residue_at(const std::string& p, std::size_t i) { return Residue{i < p.size() ? p[i] : '*', static_cast<std::int64_t>(i + 1)}; }

inline Location<Residue> residue_range(const std::string& p, std::size_t a, std::size_t b) {  // 0-based inclusive
  Location<Residue> l;
  l.start.first = l.start.last = residue_at(p, a);
  l.end.first = l.end.last = residue_at(p, b);
  l.range = a != b;
  if (!l.range) l.end = l.start;
  return l;
}

// VEP's intron effects of a change [vs, ve] (1-based genomic, an insertion vs = ve + 1) on a
// transcript's introns (BaseTranscriptVariationAllele::_intron_effects, VariationEffect::_intron_overlap).
struct IntronEffects {
  bool donor = false, acceptor = false, fifth = false, donor_region = false, polypyrimidine = false, region = false, intronic = false;
};
inline bool overlap(std::int64_t s1, std::int64_t e1, std::int64_t s2, std::int64_t e2) { return e1 >= s2 && s1 <= e2; }

inline IntronEffects intron_effects(const Transcript& t, std::int64_t vs, std::int64_t ve) {
  IntronEffects x;
  const bool insertion = vs == ve + 1;
  for (std::size_t k = 0; k + 1 < t.exons.size(); k++) {
    const Exon &e = t.exons[k], &f = t.exons[k + 1];
    // The intron, 1-based inclusive on the forward strand.
    const std::int64_t is = static_cast<std::int64_t>(t.reverse ? f.end : e.end) + 1, ie = static_cast<std::int64_t>(t.reverse ? e.begin : f.begin);
    if (ie < is) continue;
    const std::int64_t lo = std::min(vs, ve), hi = std::max(vs, ve);
    if (!overlap(lo, hi, is - 3, ie + 3) && !(insertion && overlap(lo - 1, hi + 1, is - 3, ie + 3))) continue;
    if (overlap(vs, ve, is + 2, ie - 2) || (insertion && (vs == is + 2 || ve == ie - 2))) x.intronic = true;
    const bool start_site = overlap(vs, ve, is, is + 1), end_site = overlap(vs, ve, ie - 1, ie);
    if (t.reverse ? end_site : start_site) x.donor = true;
    if (t.reverse ? start_site : end_site) x.acceptor = true;
    if (t.reverse ? overlap(vs, ve, ie - 4, ie - 4) : overlap(vs, ve, is + 4, is + 4)) x.fifth = true;
    if (t.reverse ? overlap(vs, ve, ie - 5, ie - 2) : overlap(vs, ve, is + 2, is + 5)) x.donor_region = true;
    if (t.reverse ? overlap(vs, ve, is + 2, is + 16) : overlap(vs, ve, ie - 16, ie - 2)) x.polypyrimidine = true;
    if (!start_site && !end_site &&
        (overlap(vs, ve, is + 2, is + 7) || overlap(vs, ve, ie - 7, ie - 2) || overlap(vs, ve, is - 3, is - 1) || overlap(vs, ve, ie + 1, ie + 3) ||
         (insertion && (vs == is || ve == ie || vs == is + 2 || ve == ie - 2))))
      x.region = true;
  }
  return x;
}

// Where along the transcript a genomic change lies: [a, a + d) replaced (an insertion d = 0 before a),
// when every base it names (an insertion's two flanks) is exonic and consecutive along it.
inline bool along_change(const Transcript& t, const Change& c, std::size_t& a, std::size_t& d, std::string& inserted) {
  d = c.deleted.size();
  if (d > 0) {
    std::size_t lo = SIZE_MAX, hi = 0;
    for (std::size_t g = c.position; g < c.position + d; g++) {
      const std::size_t x = along(t.exons, t.reverse, g);
      if (x == SIZE_MAX) return false;
      lo = std::min(lo, x);
      hi = std::max(hi, x);
    }
    if (hi - lo + 1 != d) return false;
    a = lo;
  } else {
    if (c.position == 0) return false;
    const std::size_t x = along(t.exons, t.reverse, c.position - 1), y = along(t.exons, t.reverse, c.position);
    if (x == SIZE_MAX || y == SIZE_MAX) return false;
    if (t.reverse ? x != y + 1 : y != x + 1) return false;
    a = t.reverse ? x : y;
  }
  inserted = t.reverse ? sequence::reverse_complement(c.inserted) : c.inserted;
  return true;
}

}  // namespace detail

// The mature transcript's sequence (its exons joined, in its direction).
inline Codes transcript_sequence(const Codes& genome, const Transcript& t) {
  Codes m;
  m.reserve(t.length());
  for (const auto& e : t.exons) {
    if (!t.reverse) m.insert(m.end(), genome.begin() + static_cast<std::ptrdiff_t>(e.begin), genome.begin() + static_cast<std::ptrdiff_t>(e.end));
    else
      for (std::size_t i = e.end; i-- > e.begin;) m.push_back(sequence::complement(genome[i]));
  }
  return m;
}

inline ProteinPrediction predict(const Codes& genome, const Change& change, const Transcript& t) {
  ProteinPrediction out;
  out.transcript = t.id;
  out.protein = t.protein_id;
  out.projection = project(genome, change, t);
  const Change& c = out.projection.change;
  auto& terms = out.consequences;
  const auto vs = static_cast<std::int64_t>(c.position) + 1;
  const auto ve = static_cast<std::int64_t>(c.position + c.deleted.size());  // an insertion: ve = vs − 1
  const auto lo = static_cast<std::int64_t>(t.lowest()) + 1, hi = static_cast<std::int64_t>(t.highest());
  const auto finish = [&] { sort_by_severity(terms); return out; };
  // Outside the transcript.
  if (std::max(vs, ve) < lo || std::min(vs, ve) > hi) {
    const bool before = (std::max(vs, ve) < lo) != t.reverse;
    const std::int64_t distance = std::max(vs, ve) < lo ? lo - std::max(vs, ve) : std::min(vs, ve) - hi;
    terms.push_back(distance <= 5000 ? (before ? "upstream_gene_variant" : "downstream_gene_variant") : "intergenic_variant");
    if (t.coding) { out.described = true; out.whole = "(=)"; }
    out.rna = "r.(=)";
    out.note = "outside the transcript";
    return finish();
  }
  if (!c.deleted.empty() && vs <= lo && ve >= hi) {
    terms.push_back("transcript_ablation");
    if (t.coding) { out.described = true; out.whole = "0?"; }
    out.rna = "r.0?";
    out.note = "the whole transcript is deleted";
    return finish();
  }
  const detail::IntronEffects ie = detail::intron_effects(t, vs, ve);
  if (ie.donor) terms.push_back("splice_donor_variant");
  if (ie.acceptor) terms.push_back("splice_acceptor_variant");
  if (ie.fifth) terms.push_back("splice_donor_5th_base_variant");
  else if (ie.donor_region) terms.push_back("splice_donor_region_variant");
  if (ie.polypyrimidine) terms.push_back("splice_polypyrimidine_tract_variant");
  if (ie.region && !ie.donor && !ie.acceptor && !ie.fifth && !ie.donor_region) terms.push_back("splice_region_variant");
  if (ie.intronic) terms.push_back("intron_variant");
  std::size_t a = 0, d = 0;
  std::string inserted;
  const bool exonic = detail::along_change(t, c, a, d, inserted);
  if (!exonic) {
    const bool splice = out.projection.exonic || ie.donor || ie.acceptor || ie.fifth || ie.donor_region || ie.region;
    if (t.coding) out.described = true;
    if (splice) {
      out.whole = "?";
      out.rna = "r.spl?";
      out.note = "the change reaches a splice site or its region, or crosses an exon's edge: its effect on the protein cannot be predicted";
      if (out.projection.exonic) terms.push_back(t.coding ? "coding_sequence_variant" : "non_coding_transcript_exon_variant");
    } else {
      out.whole = "(=)";
      out.rna = "r.(=)";
      out.note = "intronic, away from the splice sites: no change of the protein is expected";
      if (!ie.intronic && terms.empty()) terms.push_back("intron_variant");
    }
    if (!t.coding) terms.push_back("non_coding_transcript_variant");
    return finish();
  }
  // The variant mRNA.
  const Codes m = transcript_sequence(genome, t);
  Codes mv(m.begin(), m.begin() + static_cast<std::ptrdiff_t>(a));
  for (char x : inserted) mv.push_back(sequence::dna_code(x));
  mv.insert(mv.end(), m.begin() + static_cast<std::ptrdiff_t>(a + d), m.end());
  {
    // r. predicted: the c. or n. description in RNA's letters.
    std::string r = format_variant(out.projection.variant);
    for (char& x : r) {
      if (x == 'T') x = 'u';
      else if (x >= 'A' && x <= 'Z' && std::string_view("ACGUNBDHKMRSVWY").find(x) != std::string_view::npos) x = static_cast<char>(x - 'A' + 'a');
    }
    out.rna = "r.(" + r + ")";
  }
  if (!t.coding) {
    terms.push_back("non_coding_transcript_exon_variant");
    out.note = "a non-coding transcript";
    return finish();
  }
  const sequence::GeneticCode& code = sequence::genetic_code(t.genetic_code);
  const std::size_t cb = t.cds_begin, ce = t.cds_end;
  const std::size_t end = a + d;
  if ((d > 0 && end <= cb) || (d == 0 && a <= cb)) {
    terms.push_back("5_prime_UTR_variant");
    out.described = true;
    out.whole = "(=)";
    out.note = "in the 5′ UTR";
    return finish();
  }
  if (a >= ce) {
    terms.push_back("3_prime_UTR_variant");
    out.described = true;
    out.whole = "(=)";
    out.note = "in the 3′ UTR";
    return finish();
  }
  // The proteins.
  std::string P = detail::translate_from(m, cb, code);
  if (P.size() > (ce - cb) / 3) P.resize((ce - cb) / 3);  // the CDS's own end (an incomplete stop is completed below)
  if (!P.empty() && cb + 3 <= m.size() && sequence::is_start(code, m[cb], m[cb + 1], m[cb + 2])) P[0] = 'M';
  if (P.empty() || P.back() != '*') P += '*';  // a stop completed by polyadenylation (mitochondrial genes)
  std::string Q = detail::translate_from(mv, cb, code);
  const bool q_stops = !Q.empty() && Q.back() == '*';
  if (!q_stops && t.exons.size() == 1 && mv.size() == cb + (ce - cb) + inserted.size() - d && (mv.size() - cb) % 3 != 0) {
    // A CDS-only transcript ending in an incomplete codon: the stop is completed as the reference's.
    Q += '*';
  }
  if (!Q.empty() && cb + 3 <= mv.size() && sequence::is_start(code, mv[cb], mv[cb + 1], mv[cb + 2])) Q[0] = 'M';
  out.reference = P;
  out.observed = Q;
  out.described = true;
  // The start codon.
  const bool touches_start = (d > 0 && a < cb + 3 && end > cb) || (d == 0 && a > cb && a < cb + 3);
  if (touches_start) {
    const bool still = cb + 3 <= mv.size() && sequence::is_start(code, mv[cb], mv[cb + 1], mv[cb + 2]);
    if (!still || inserted.size() != d) {
      terms.push_back("start_lost");
      out.whole = "(Met1?)";
      out.note = "the start codon is lost: the protein's consequence cannot be predicted";
      return finish();
    }
    terms.push_back("start_retained_variant");
  }
  const std::size_t n = P.size() - 1;  // P's stop
  std::size_t k = 0;
  while (k < P.size() && k < Q.size() && P[k] == Q[k]) k++;
  const std::int64_t delta = static_cast<std::int64_t>(inserted.size()) - static_cast<std::int64_t>(d);
  Variant& v = out.variant;
  v.protein = true;
  if (k >= P.size() && Q.size() == P.size()) {
    // The same protein.
    const std::size_t codon = (a - cb) / 3;
    v.kind = Kind::identity;
    v.residues = detail::residue_range(P, std::min(codon, n), std::min(codon, n));
    terms.push_back(codon >= n ? "stop_retained_variant" : "synonymous_variant");
    return finish();
  }
  if (k == n) {
    // The stop lost: an extension.
    v.kind = Kind::extension;
    v.residues = detail::residue_range(P, n, n);
    v.now = std::string(1, k < Q.size() ? Q[k] : 'X');
    v.has_stop = true;
    if (q_stops && Q.size() - 1 > n) v.stop_at = static_cast<std::int64_t>(Q.size() - 1 - n);
    else v.stop_unknown = true;
    terms.push_back("stop_lost");
    if (delta % 3 != 0) terms.push_back("frameshift_variant");
    return finish();
  }
  if (delta % 3 != 0) {
    if (k < Q.size() && Q[k] == '*') {
      v.kind = Kind::substitution;
      v.residues = detail::residue_range(P, k, k);
      v.now = "*";
      terms.push_back("stop_gained");
      terms.push_back("frameshift_variant");
      return finish();
    }
    v.kind = Kind::frameshift;
    v.residues = detail::residue_range(P, k, k);
    v.fs_new = std::string(1, k < Q.size() ? Q[k] : 'X');
    v.has_stop = true;
    if (q_stops) v.stop_at = static_cast<std::int64_t>(Q.size() - 1 - k + 1);
    else v.stop_unknown = true;
    terms.push_back("frameshift_variant");
    return finish();
  }
  // In frame: a premature stop is one before where P's stop moves by the change (delta / 3 residues).
  const std::int64_t moved_stop = static_cast<std::int64_t>(n) + delta / 3;
  if (q_stops && static_cast<std::int64_t>(Q.size()) - 1 < moved_stop) {
    terms.push_back("stop_gained");
    if (Q[k] == '*') {
      v.kind = Kind::substitution;
      v.residues = detail::residue_range(P, k, k);
      v.now = "*";
      return finish();
    }
    // The new residues up to and including the stop replace the reference residues of the codons
    // the change touches (the standard's examples; see the header).
    const std::size_t first_codon = (a - cb) / 3;
    std::int64_t last_codon;
    if (d > 0) last_codon = static_cast<std::int64_t>((end - 1 - cb) / 3);
    else last_codon = (a - cb) % 3 == 0 ? static_cast<std::int64_t>(first_codon) - 1 : static_cast<std::int64_t>(first_codon);
    const std::int64_t count = last_codon - static_cast<std::int64_t>(k) + 1;
    v.residues_inserted = Q.substr(k);
    if (count <= 0) {
      v.kind = Kind::insertion;
      v.residues = detail::residue_range(P, k - 1, k);
      v.residues.range = true;
    } else {
      v.kind = Kind::delins;
      v.residues = detail::residue_range(P, k, k + static_cast<std::size_t>(count) - 1);
    }
    if (delta > 0) terms.push_back("inframe_insertion");
    else if (delta < 0) terms.push_back("inframe_deletion");
    return finish();
  }
  // The common suffix (of the proteins up to their stops).
  std::size_t s = 0;
  while (s < P.size() - k && s < Q.size() - k && P[P.size() - 1 - s] == Q[Q.size() - 1 - s]) s++;
  const std::string D = P.substr(k, P.size() - s - k), N = Q.substr(k, Q.size() - s - k);
  if (D.size() == 1 && N.size() == 1) {
    v.kind = Kind::substitution;
    v.residues = detail::residue_range(P, k, k);
    v.now = N;
    terms.push_back(N == "*" ? "stop_gained" : "missense_variant");
    return finish();
  }
  if (N.empty()) {
    // A deletion, placed C-terminally: slide right while P[i] = P[i + L], before the stop.
    const std::size_t L = D.size();
    std::size_t i = k;
    while (i + L < n && P[i] == P[i + L]) i++;
    v.kind = Kind::deletion;
    v.residues = detail::residue_range(P, i, i + L - 1);
    terms.push_back("inframe_deletion");
    return finish();
  }
  if (D.empty()) {
    // An insertion between k − 1 and k, placed C-terminally (u rotating), a duplication when it repeats what is before it.
    std::string u = N;
    std::size_t i = k, lo_i = k;
    while (i < n && P[i] == u.front()) { std::rotate(u.begin(), u.begin() + 1, u.end()); i++; }
    std::string w = N;
    while (lo_i > 0 && P[lo_i - 1] == w.back()) { std::rotate(w.rbegin(), w.rbegin() + 1, w.rend()); lo_i--; }
    const std::size_t L = u.size();
    terms.push_back("inframe_insertion");
    if (i - lo_i >= L) {
      v.kind = Kind::duplication;
      v.residues = detail::residue_range(P, i - L, i - 1);
      return finish();
    }
    v.kind = Kind::insertion;
    v.residues = detail::residue_range(P, i - 1, i);
    v.residues.range = true;
    v.residues_inserted = u;
    return finish();
  }
  v.kind = Kind::delins;
  v.residues = detail::residue_range(P, k, k + D.size() - 1);
  v.residues_inserted = N;
  terms.push_back(N.size() > D.size() ? "inframe_insertion" : N.size() < D.size() ? "inframe_deletion" : "missense_variant");
  return finish();
}

}  // namespace genomics::nomenclature
