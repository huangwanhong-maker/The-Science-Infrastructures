// The genes an analysis of codon usage reads from an annotated genome: its protein-coding genes
// as codons, with their places on the genome, and its tRNA genes as anticodons with their amino
// acids; and the reference sets of genes named by their products.
//
// - Coding genes. Each CDS of a GenBank record (genomics::sequence::read_genbank) is extracted by
//   its location (INSDC Feature Table §3.4; across the origin of a circular genome, on either
//   strand) from its /codon_start, cut to whole codons, and its final stop dropped. A CDS is left
//   out, with the reason counted, when it is a pseudogene (/pseudo, /pseudogene), partial (a '<' or
//   '>' end), read with a frameshift (/ribosomal_slippage) or another /exception, has parts in
//   another record, is read otherwise than its table says (/transl_except: selenocysteine at a
//   TGA), has a length that is not a whole number of codons, or has a stop inside it: such codons
//   do not mean what the genetic code says they mean, and codon usage counts only what does. A
//   /transl_except that only completes the final stop (aa:TERM, a stop finished by
//   polyadenylation in vertebrate mitochondria) keeps the gene, its incomplete codon dropped. Its
//   genetic code is its /transl_table, or the caller's default (11 for bacteria, 2 for vertebrate
//   mitochondria; INSDC's default, when a record names none, is 1).
// - tRNA genes. Each tRNA feature gives its amino acid and anticodon, read 5′→3′ (positions
//   34, 35, 36, so that the first letter is the wobble base): from /anticodon=(pos:…,aa:Phe,seq:gaa)
//   (INSDC Feature Table §7.3.1), which mitochondrial RefSeq records carry; or, where it is absent,
//   from a /note of the form "tRNA-Ile(CAU)" and the /product "tRNA-Ile", as RefSeq writes
//   Escherichia coli K-12's (NC_000913.3); "tRNA-initiator Met(CAU)" marks an initiator, and
//   "tRNA-Sec" selenocysteine's tRNA, neither of which reads codons during elongation.
// - Reference sets. A reference set of highly expressed genes is chosen by product, here as a
//   regular expression over /product. The standard choice for bacteria is the ribosomal proteins
//   (Sharp and Li 1987 used ribosomal proteins with a few other abundant proteins; Sharp, Bailes,
//   Grocock, Peden and Sockett 2005, Nucleic Acids Research 33, 1141–1153, the ribosomal proteins
//   alone): products "30S ribosomal subunit protein …" and "50S ribosomal subunit protein …" in
//   RefSeq's E. coli K-12, 56 genes. Karlin and Mrázek's classes of chaperones (CH) and of
//   translation and transcription factors (TF) are their Table 2 (J. Bacteriol. 182, 5238–5250,
//   2000), given below by E. coli's gene names.
//
// Complexity: O(genome + features). Readers never throw on content; what is left out is counted.

#pragma once

#include <genomics/annotation/features.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/genbank.hpp>
#include <genomics/sequence/genetic_code.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <map>
#include <regex>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::codon {

using sequence::Code;
using sequence::Codes;
using sequence::GeneticCode;
using sequence::Problem;

// Counts of the 64 codons in NCBI's order (TCAG on each base; sequence::codon_index).
using CodonCounts = std::array<double, 64>;

struct CodingGene {
  std::string locus_tag, gene, product, protein_id;
  std::size_t start = 0, stop = 0;  // the extent on the forward strand, [start, stop), or [start, length) ∪ [0, stop) when it wraps
  bool reverse = false, wraps = false;
  int table = 11;                   // the NCBI genetic code it is read with
  Codes cds;                        // the coding bases from /codon_start, whole codons, the final stop dropped
  std::size_t feature = 0;          // its index among the record's features
  std::size_t codons() const { return cds.size() / 3; }
};

// Why CDSs were left out.
struct Skipped {
  std::size_t pseudo = 0, partial = 0, exception = 0, remote = 0, transl_except = 0, length = 0, internal_stop = 0, short_gene = 0;
  std::size_t total() const { return pseudo + partial + exception + remote + transl_except + length + internal_stop + short_gene; }
};

struct GeneOptions {
  int default_table = 11;         // when a CDS has no /transl_table
  std::size_t min_codons = 0;     // genes shorter than this (stop excluded) are left out
};

struct CodingGenes {
  std::string accession;
  std::size_t length = 0;
  bool circular = false;
  std::vector<CodingGene> genes;  // in the order of the record's features
  Skipped skipped;
  std::vector<Problem> problems;
};

namespace detail {

inline std::string lower(std::string s) {
  for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

}  // namespace detail

inline CodingGenes coding_genes(const sequence::GenBank& record, const GeneOptions& o = {}) {
  CodingGenes out;
  out.accession = record.version.empty() ? record.accession : record.version;
  out.length = record.length ? record.length : record.sequence.size();
  out.circular = record.circular();
  const Codes genome = sequence::encode_dna(record.sequence);
  // The extents of the features, as the annotation module draws them (wrapping across the origin).
  annotation::TrackOptions to;
  to.check_translations = false;
  const auto tracks = annotation::feature_tracks(record, to);
  std::vector<const annotation::TrackFeature*> track_of(record.features.size(), nullptr);
  for (const auto& t : tracks.features) track_of[t.feature] = &t;

  for (std::size_t i = 0; i < record.features.size(); i++) {
    const sequence::Feature& f = record.features[i];
    if (f.key != "CDS") continue;
    if (f.has("pseudo") || f.has("pseudogene")) { out.skipped.pseudo++; continue; }
    if (f.location.remote) { out.skipped.remote++; continue; }
    bool partial = false;
    for (const auto& p : f.location.parts) partial = partial || p.partial_lower || p.partial_upper;
    if (partial) { out.skipped.partial++; continue; }
    if (f.has("ribosomal_slippage") || f.has("exception") || f.has("artificial_location")) { out.skipped.exception++; continue; }
    // A /transl_except that only completes the final stop (aa:TERM: a stop finished by
    // polyadenylation, as in vertebrate mitochondria) leaves the codons as they are; any other
    // (selenocysteine at a TGA) recodes one, and the gene is left out.
    bool recoded = false, completed_stop = false;
    for (const auto& [key, value] : f.qualifiers)
      if (key == "transl_except") (value.find("aa:TERM") != std::string::npos ? completed_stop : recoded) = true;
    if (recoded) { out.skipped.transl_except++; continue; }
    const std::string start = f.qualifier("codon_start"), table = f.qualifier("transl_table");
    const int t = table.empty() ? o.default_table : std::atoi(table.c_str());
    const GeneticCode* code = nullptr;
    for (const auto& g : sequence::genetic_codes())
      if (g.id == t) code = &g;
    if (!code) {
      out.problems.push_back({f.line, "CDS " + f.qualifier("locus_tag") + ": NCBI has no translation table " + table});
      continue;
    }
    Codes x = sequence::extract(genome, f.location);
    const std::size_t frame = start.empty() ? 0 : static_cast<std::size_t>(std::max(0, std::atoi(start.c_str()) - 1));
    if (frame >= x.size()) { out.skipped.length++; continue; }
    x.erase(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(frame));
    if (completed_stop) x.resize(x.size() - x.size() % 3);  // the incomplete stop's one or two bases
    if (x.size() % 3 != 0) { out.skipped.length++; continue; }
    if (x.size() >= 3 && sequence::is_stop(*code, x[x.size() - 3], x[x.size() - 2], x[x.size() - 1])) x.resize(x.size() - 3);
    bool internal = false;
    for (std::size_t k = 0; k + 3 <= x.size(); k += 3) internal = internal || sequence::is_stop(*code, x[k], x[k + 1], x[k + 2]);
    if (internal) { out.skipped.internal_stop++; continue; }
    if (x.size() / 3 < o.min_codons) { out.skipped.short_gene++; continue; }
    CodingGene g;
    g.locus_tag = f.qualifier("locus_tag");
    g.gene = f.qualifier("gene");
    g.product = f.qualifier("product");
    g.protein_id = f.qualifier("protein_id");
    g.reverse = f.location.reverse();
    if (const auto* tr = track_of[i]) {
      g.start = tr->start;
      g.stop = tr->stop;
      g.wraps = tr->wraps;
    } else {
      g.start = f.location.lower();
      g.stop = f.location.upper();
    }
    g.table = t;
    g.cds = std::move(x);
    g.feature = i;
    out.genes.push_back(std::move(g));
  }
  return out;
}

// The counts of a gene's codons (codons with an ambiguous base are skipped).
inline CodonCounts codon_counts(const Codes& cds) {
  CodonCounts n{};
  for (std::size_t i = 0; i + 3 <= cds.size(); i += 3) {
    const int k = sequence::codon_index(cds[i], cds[i + 1], cds[i + 2]);
    if (k >= 0) n[static_cast<std::size_t>(k)] += 1;
  }
  return n;
}
inline CodonCounts codon_counts(const CodingGene& g) { return codon_counts(g.cds); }

// The summed counts of a set of genes.
inline CodonCounts pooled_counts(const std::vector<CodingGene>& genes, const std::vector<std::size_t>& which) {
  CodonCounts n{};
  for (std::size_t i : which) {
    const auto c = codon_counts(genes[i]);
    for (std::size_t k = 0; k < 64; k++) n[k] += c[k];
  }
  return n;
}
inline CodonCounts pooled_counts(const std::vector<CodingGene>& genes) {
  CodonCounts n{};
  for (const auto& g : genes) {
    const auto c = codon_counts(g);
    for (std::size_t k = 0; k < 64; k++) n[k] += c[k];
  }
  return n;
}

// ---- reference sets ---------------------------------------------------------------------------------
// The indices of the genes whose /product matches a regular expression (ECMAScript, case-insensitive).
inline std::vector<std::size_t> genes_by_product(const std::vector<CodingGene>& genes, const std::string& pattern) {
  const std::regex re(pattern, std::regex::ECMAScript | std::regex::icase);
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < genes.size(); i++)
    if (std::regex_search(genes[i].product, re)) out.push_back(i);
  return out;
}

// The ribosomal proteins, as RefSeq names them in bacteria ("30S ribosomal subunit protein S12",
// "50S ribosomal protein L2"); not their modifying enzymes or assembly factors (a product that goes
// on, as "… protein L16-arginine 3-hydroxylase", is excluded by the anchors). In NC_000913.3: 56.
inline constexpr const char* RIBOSOMAL_PROTEINS = R"(^(30S|50S) ribosomal (subunit )?protein [LS][0-9]+[A-Z]?$)";

// The indices of the genes with one of the given /gene names.
inline std::vector<std::size_t> genes_by_name(const std::vector<CodingGene>& genes, const std::vector<std::string>& names) {
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < genes.size(); i++)
    if (std::find(names.begin(), names.end(), genes[i].gene) != names.end()) out.push_back(i);
  return out;
}
// Karlin and Mrázek (2000), Table 2, in E. coli's gene names: the chaperones and degradation
// proteins (CH: GroEL, DnaK, GrpE, HtpG, ClpP, ClpB, ClpX, FtsH, HslU, trigger factor, a PPIase
// (here cyclophilin B, ppiB: the choice is ours, the table names none), thioredoxin reductase; E.
// coli has no ClpC) and the translation and transcription factors (TF: EF-Tu, EF-G, EF-Ts, IF-2,
// IF-3, RNA polymerase α, β, β′, NusA, NusB, NusG, ribosome recycling factor).
inline const std::vector<std::string>& ecoli_chaperones() {
  static const std::vector<std::string> names = {"groL", "dnaK", "grpE", "htpG", "clpP", "clpB", "clpX", "ftsH", "hslU", "tig", "ppiB", "trxB"};
  return names;
}
inline const std::vector<std::string>& ecoli_translation_factors() {
  static const std::vector<std::string> names = {"tufA", "tufB", "fusA", "tsf", "infB", "infC", "rpoA", "rpoB", "rpoC", "nusA", "nusB", "nusG", "frr"};
  return names;
}

// ---- tRNA genes -----------------------------------------------------------------------------------------
struct TrnaGene {
  std::string amino_acid;   // three letters as INSDC writes them: Ala, …, fMet (initiator), Sec, Pyl
  Codes anticodon;          // 5′→3′: positions 34 (the wobble base), 35, 36
  std::string locus_tag, gene;
  std::size_t start = 0, stop = 0;
  bool reverse = false;
  bool initiator = false;   // reads the start codon only
  std::size_t feature = 0;
};

struct TrnaGenes {
  std::vector<TrnaGene> genes;
  std::size_t unread = 0;   // tRNA features whose amino acid or anticodon could not be read
  std::vector<Problem> problems;
};

namespace detail {

inline Codes anticodon_codes(std::string_view s) {
  Codes c;
  for (char ch : s)
    if (std::isalpha(static_cast<unsigned char>(ch))) c.push_back(sequence::dna_code(ch));  // U reads as T
  return c;
}

}  // namespace detail

inline TrnaGenes trna_genes(const sequence::GenBank& record) {
  TrnaGenes out;
  for (std::size_t i = 0; i < record.features.size(); i++) {
    const sequence::Feature& f = record.features[i];
    if (f.key != "tRNA" || f.has("pseudo") || f.has("pseudogene")) continue;
    TrnaGene t;
    t.locus_tag = f.qualifier("locus_tag");
    t.gene = f.qualifier("gene");
    t.start = f.location.lower();
    t.stop = f.location.upper();
    t.reverse = f.location.reverse();
    t.feature = i;
    const std::string ac = f.qualifier("anticodon"), note = f.qualifier("note"), product = f.qualifier("product");
    std::string aa, seq;
    if (!ac.empty()) {
      // (pos:611..613,aa:Phe,seq:gaa)
      const auto a = ac.find("aa:"), s = ac.find("seq:");
      if (a != std::string::npos) aa = ac.substr(a + 3, ac.find_first_of(",)", a) - a - 3);
      if (s != std::string::npos) seq = ac.substr(s + 4, ac.find_first_of(",)", s) - s - 4);
    } else {
      // /note="tRNA-Ile(CAU)" or "tRNA-initiator Met(CAU)"
      static const std::regex re(R"(tRNA-(initiator )?([A-Za-z]{3,4})\(([ACGTUacgtu]{3})\))");
      std::smatch m;
      if (std::regex_search(note, m, re)) {
        t.initiator = m[1].matched;
        aa = m[2];
        seq = m[3];
      }
    }
    if (aa.empty() && product.size() > 5 && product.substr(0, 5) == "tRNA-") aa = product.substr(5);
    if (aa == "fMet" || note.find("initiator") != std::string::npos) t.initiator = true;
    if (aa == "fMet") aa = "Met";
    t.amino_acid = aa;
    t.anticodon = detail::anticodon_codes(seq);
    if (aa.empty() || t.anticodon.size() != 3 || sequence::count_ambiguous(t.anticodon)) {
      out.unread++;
      out.problems.push_back({f.line, "tRNA " + t.locus_tag + ": its amino acid or anticodon cannot be read"});
      continue;
    }
    out.genes.push_back(std::move(t));
  }
  return out;
}

// The one-letter code of a three-letter amino acid name, as the genetic code tables write it
// ('U' for Sec, 'O' for Pyl, 0 when unknown).
inline char amino_acid_letter(std::string_view three) {
  static const std::pair<const char*, char> table[] = {
      {"Ala", 'A'}, {"Arg", 'R'}, {"Asn", 'N'}, {"Asp", 'D'}, {"Cys", 'C'}, {"Gln", 'Q'}, {"Glu", 'E'}, {"Gly", 'G'}, {"His", 'H'}, {"Ile", 'I'}, {"Leu", 'L'},
      {"Lys", 'K'}, {"Met", 'M'}, {"Phe", 'F'}, {"Pro", 'P'}, {"Ser", 'S'}, {"Thr", 'T'}, {"Trp", 'W'}, {"Tyr", 'Y'}, {"Val", 'V'}, {"Sec", 'U'}, {"Pyl", 'O'}};
  for (const auto& [name, letter] : table)
    if (three == name) return letter;
  return 0;
}

}  // namespace genomics::codon
