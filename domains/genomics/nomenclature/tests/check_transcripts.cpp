// Checks of transcript.hpp and protein.hpp: the numbering of transcripts and the projection of
// genomic changes through them to proteins.
//
// - Numbering on a constructed transcript (introns of even and odd length, both UTRs): every genomic
//   base maps to a c. position and back; c.-1 and c.*1 flank the CDS, there is no c.0; the middle of an
//   intron changes from + to −, its central base of an odd intron numbered with + (numbering.md).
// - Strand symmetry: the reverse complement of the genome with the transcript mirrored gives every
//   base the same number, and every change projected on it the same c. description (the 3′ rule in the
//   transcript's direction).
// - Real genes (when the data library has fetched them): HBB (minus strand), CFTR (plus), BRCA1
//   (minus, with every RefSeq isoform), the mitochondrial genome (genetic code 2): every CDS's
//   translation is RefSeq's /translation; the c. numbering of exon edges, intron ends and the middle of
//   every intron is Ensembl VEP's (tests/fixtures/nomenclature-positions.tsv); real variants
//   (nomenclature-variants.tsv) are described in g., SPDI and VCF as Mutalyzer 3 and NCBI's Variation
//   Services describe them, projected to c. and p. as Mutalyzer does, with VEP's most severe
//   consequence.

#include <genomics/nomenclature/protein.hpp>
#include <genomics/nomenclature/spdi.hpp>
#include <genomics/sequence/genbank.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace nm = genomics::nomenclature;
namespace sq = genomics::sequence;
using sq::Codes;

namespace {

std::string read_file(const std::filesystem::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

std::vector<std::vector<std::string>> table(const std::string& name, std::string* header = nullptr) {
  std::ifstream in(std::string(GENOMICS_FIXTURES) + "/" + name);
  std::vector<std::vector<std::string>> rows;
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("# retrieved", 0) == 0 && header) *header = line.substr(2);
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> f;
    std::stringstream ss(line);
    std::string x;
    while (std::getline(ss, x, '\t')) f.push_back(x);
    while (f.size() < 14) f.emplace_back();
    rows.push_back(f);
  }
  return rows;
}

// The part of a description after its reference's "c." (or "p."): NC_000011.10(NM_000518.5):c.20A>T → 20A>T.
std::string after(const std::string& d, const std::string& type) {
  const auto at = d.find(":" + type + ".");
  return at == std::string::npos ? std::string() : d.substr(at + 3);
}

nm::Transcript mirror(const nm::Transcript& t, std::size_t n) {
  nm::Transcript m = t;
  m.reverse = !t.reverse;
  for (auto& e : m.exons) e = {n - e.end, n - e.begin};
  return m;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- VEP's table ------------------------------------------------------------------------------------------------
  {
    const auto rows = table("vep-consequences.tsv");
    bool same = rows.size() == nm::consequence_terms().size();
    for (std::size_t i = 0; same && i < rows.size(); i++) {
      const auto& t = nm::consequence_terms()[i];
      same = rows[i][0] == t.term && rows[i][1] == t.accession && rows[i][2] == t.impact && std::stoi(rows[i][3]) == t.rank;
    }
    s.require("the consequence terms, accessions, impacts and ranks are Ensembl release 116's Constants.pm", same, std::to_string(rows.size()) + " terms");
  }

  // ---- numbering on a constructed transcript ------------------------------------------------------------------
  {
    calx::random::Stream r(20261101, 0x5452414e, 0, 0);
    Codes genome(420);
    for (auto& b : genome) b = static_cast<sq::Code>(r.below(4));
    nm::Transcript t;
    t.id = "T1";
    t.exons = {{50, 80}, {100, 131}, {150, 200}, {260, 330}};  // introns of 20, 19 and 60 bases
    t.coding = true;
    t.cds_begin = 10;
    t.cds_end = 10 + 3 * 40;  // inside the third exon
    bool round = true;
    for (std::size_t g = 0; g < genome.size(); g++) {
      const nm::Mapped m = nm::map_genomic(t, g);
      const auto back = nm::map_transcript(t, m.position);
      round = round && back && *back == g;
    }
    s.require("every base of a constructed transcript and its flanks maps to its c. position and back", round);
    const auto at = [&](std::size_t g) { return nm::format_position(nm::map_genomic(t, g).position); };
    s.require("c.-1 is the base before the A of ATG, c.1 the A; c.*1 the base after the stop; no c.0",
              at(59) == "-1" && at(60) == "1" && at(278) == "120" && at(279) == "*1",
              at(59) + " " + at(60) + " " + at(278) + " " + at(279));
    s.require("an intron of 20: c.20+10 then c.21-10 at its middle", at(89) == "20+10" && at(90) == "21-10", at(89) + " " + at(90));
    s.require("an intron of 19: its central base c.51+10, from the upstream exon with + (numbering.md)", at(140) == "51+10" && at(141) == "52-9", at(140) + " " + at(141));
    s.require("upstream and downstream of the transcript the numbering goes on", at(0) == "-60" && at(419).front() == '*', at(0) + " " + at(419));
    bool wrong_anchor_refused = !nm::map_transcript(t, nm::parse_hgvs("c.19+1A>G").value.alternatives[0].variants.alternatives[0].groups[0].chain.items[0].choices[0].location.start.first);
    s.require("an offset not counted from an exon's edge is refused (c.19+1 where c.20 ends the exon)", wrong_anchor_refused);

    // Strand symmetry.
    const Codes rc = sq::reverse_complement(genome);
    const nm::Transcript tm = mirror(t, genome.size());
    bool same_numbers = true;
    for (std::size_t g = 0; g < genome.size(); g++)
      same_numbers = same_numbers && nm::format_position(nm::map_genomic(t, g).position) == nm::format_position(nm::map_genomic(tm, genome.size() - 1 - g).position);
    s.require("the reverse complement with the transcript mirrored numbers every base the same", same_numbers);
    std::size_t same = 0, n = 0;
    std::string first_bad;
    for (int k = 0; k < 400; k++) {
      nm::Change c;
      c.position = 40 + r.below(300);
      const auto kind = r.below(4);
      if (kind == 0) c.deleted = sq::decode_dna(Codes(genome.begin() + static_cast<std::ptrdiff_t>(c.position), genome.begin() + static_cast<std::ptrdiff_t>(c.position + 1 + r.below(3))));
      else if (kind == 1) for (std::size_t i = 0, L = 1 + r.below(3); i < L; i++) c.inserted += "ACGT"[r.below(4)];
      else { c.deleted = std::string(1, sq::dna_letter(genome[c.position])); c.inserted = std::string(1, "ACGT"[(genome[c.position] + 1 + r.below(3)) % 4]); }
      if (kind == 3) c.inserted = sq::decode_dna(Codes(genome.begin() + static_cast<std::ptrdiff_t>(c.position) - 2, genome.begin() + static_cast<std::ptrdiff_t>(c.position)));
      const nm::Change mc{genome.size() - c.position - c.deleted.size(), sq::reverse_complement(c.deleted), sq::reverse_complement(c.inserted)};
      const std::string a = nm::project(genome, c, t).text(), b = nm::project(rc, mc, tm).text();
      n++;
      same += a == b;
      if (a != b && first_bad.empty()) first_bad = a + " vs " + b;
    }
    s.require("400 changes projected on both strands give the same c. description", same == n, std::to_string(same) + " of " + std::to_string(n) + (first_bad.empty() ? "" : "; " + first_bad));
  }

  // ---- real genes ------------------------------------------------------------------------------------------------------
  const char* env = std::getenv("SCIENCE_CACHE");
  const std::filesystem::path cache = env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
  struct Record {
    std::string accession;
    std::int64_t offset = 0;
    sq::GenBank gb;
    Codes genome;
    std::vector<nm::Transcript> transcripts;
  };
  std::map<std::string, Record> records;
  const std::vector<std::tuple<std::string, std::string, std::int64_t>> files = {
      {"HBB", "ncbi-human-genes/NC_000011.10_5222001-5232000.gb", 5222000},
      {"CFTR", "ncbi-human-genes/NC_000007.14_117479001-117670000.gb", 117479000},
      {"BRCA1", "ncbi-human-genes/NC_000017.11_43043001-43172000.gb", 43043000},
      {"MT", "ncbi-genomes/NC_012920.1.gb", 0},
  };
  for (const auto& [gene, file, offset] : files) {
    const auto path = cache / file;
    if (!std::filesystem::exists(path)) continue;
    auto recs = sq::read_genbank(read_file(path));
    if (recs.empty()) continue;
    Record rec;
    rec.gb = std::move(recs.front());
    rec.accession = rec.gb.version;
    rec.offset = offset;
    rec.genome = sq::encode_dna(rec.gb.sequence);
    std::vector<sq::Problem> problems;
    rec.transcripts = nm::transcripts(rec.gb, &problems);
    records[gene] = std::move(rec);
  }
  if (records.empty()) {
    s.skip("the genes of the data library are not in the cache (python data/acquisition/fetch.py ncbi-human-genes; ncbi-genomes): the real transcripts are not checked");
    return s.done();
  }
  const auto find = [&](const std::string& gene, const std::string& id) -> const nm::Transcript* {
    for (const auto& t : records[gene].transcripts)
      if (t.id == id) return &t;
    return nullptr;
  };
  // Translations.
  for (auto& [gene, rec] : records) {
    std::size_t coding = 0, agree = 0;
    std::string bad;
    for (const auto& t : rec.transcripts) {
      if (!t.coding || t.translation.empty()) continue;
      coding++;
      const Codes m = nm::transcript_sequence(rec.genome, t);
      std::string p = sq::translate(Codes(m.begin() + static_cast<std::ptrdiff_t>(t.cds_begin), m.begin() + static_cast<std::ptrdiff_t>(std::min(t.cds_end, m.size()))), sq::genetic_code(t.genetic_code), 0, true);
      if (p == t.translation) agree++;
      else if (bad.empty()) bad = t.id;
    }
    s.require(gene + ": every CDS translated through its transcript model is RefSeq's /translation", coding > 0 && agree == coding,
              std::to_string(agree) + " of " + std::to_string(coding) + " transcripts" + (bad.empty() ? "" : ", first differing " + bad));
  }
  // Positions against VEP.
  {
    std::string header;
    const auto rows = table("nomenclature-positions.tsv", &header);
    std::size_t n = 0, agree = 0;
    std::string bad;
    for (const auto& f : rows) {
      const std::string gene = f[1] == "NM_000518.5" ? "HBB" : f[1] == "NM_000492.4" ? "CFTR" : "BRCA1";
      if (!records.count(gene)) continue;
      const nm::Transcript* t = find(gene, f[1]);
      if (!t) continue;
      const std::int64_t g = std::stoll(f[2]) - records[gene].offset - 1;
      n++;
      const std::string ours = nm::format_position(nm::map_genomic(*t, static_cast<std::size_t>(g)).position);
      if (ours == f[3]) agree++;
      else if (bad.size() < 400) bad += f[1] + " g." + f[2] + ": " + ours + " vs " + f[3] + "; ";
    }
    if (!n) s.skip("the position fixture is absent (tests/fixtures/nomenclature-positions.tsv)");
    else s.require("the c. numbering of exon edges, intron ends and every intron's middle is VEP's", agree == n, std::to_string(agree) + " of " + std::to_string(n) + " positions (" + header + ")" + (bad.empty() ? "" : ": " + bad));
  }
  // Variants against Mutalyzer, NCBI and VEP.
  {
    std::string header;
    const auto rows = table("nomenclature-variants.tsv", &header);
    std::size_t n = 0, g_ok = 0, spdi_ok = 0, vcf_ok = 0, c_ok = 0, c_n = 0, p_ok = 0, p_n = 0, vep_ok = 0, vep_n = 0;
    std::string bad;
    const auto note = [&](const std::string& x) { if (bad.size() < 1500) bad += x + "; "; };
    for (const auto& f : rows) {
      const std::string& gene = f[0];
      if (!records.count(gene) || f[6].empty()) continue;
      const Record& rec = records[gene];
      n++;
      // The genomic change, from Mutalyzer's genomic equivalent.
      const auto gp = nm::parse_hgvs(f[6]);
      if (!gp.accepted) { note("unreadable " + f[6]); continue; }
      const auto& v = gp.value.alternatives.at(0).variants.alternatives.at(0).groups.at(0).chain.items.at(0).choices.at(0);
      const auto change = nm::change_of(rec.genome, v, rec.offset + 1);
      if (!change) { note("no change of " + f[6]); continue; }
      const char type = gene == "MT" ? 'm' : 'g';
      nm::DescribeOptions o;
      o.first = rec.offset + 1;
      const std::string ours_g = rec.accession + ":" + type + "." + nm::describe(rec.genome, *change, o).text();
      if (ours_g == f[6]) g_ok++; else note(ours_g + " vs Mutalyzer " + f[6]);
      // SPDI and VCF, against NCBI (positions of the whole chromosome).
      nm::Spdi sp = nm::spdi_form(rec.genome, *change, rec.accession);
      const auto vf = nm::vcf_form(rec.genome, *change);
      sp.position += static_cast<std::size_t>(rec.offset);
      if (!f[7].empty()) {
        if (sp.text() == f[7]) spdi_ok++; else note(sp.text() + " vs NCBI " + f[7]);
        const std::string vcf = rec.accession + ":" + std::to_string(vf.position + 1 + static_cast<std::size_t>(rec.offset)) + ":" + vf.ref + ":" + vf.alts.at(0);
        if (vcf == f[8]) vcf_ok++; else note(vcf + " vs NCBI VCF " + f[8]);
      }
      // The transcript and the protein.
      const std::string transcript = [&] {
        const auto a = f[1].find('('), b = f[1].find(')');
        return a == std::string::npos ? std::string() : f[1].substr(a + 1, b - a - 1);
      }();
      const nm::Transcript* t = transcript.empty() ? nullptr : find(gene, transcript);
      if (gene == "MT") {
        // The gene named by VEP (MT-ATP6 → ATP6), its protein against VEP's.
        const std::string symbol = f[10].substr(0, f[10].find('|'));
        std::string bare = symbol.rfind("MT-", 0) == 0 ? symbol.substr(3) : symbol;
        if (bare.size() <= 3 && bare[0] == 'T') bare = "TRN" + bare.substr(1);  // HGNC's MT-TL1 is RefSeq's TRNL1
        for (const auto& x : rec.transcripts)
          if (x.gene == bare || x.id == bare) t = &x;
        if (!t) { note("no transcript for " + symbol); continue; }
        const nm::ProteinPrediction pp = nm::predict(rec.genome, *change, *t);
        const std::string vep_terms = f[10].substr(f[10].find('|') + 1);
        vep_n++;
        if (!pp.consequences.empty() && vep_terms.substr(0, vep_terms.find(',')) == pp.consequences.front()) vep_ok++;
        else note(f[1] + ": " + (pp.consequences.empty() ? "" : pp.consequences.front()) + " vs VEP " + vep_terms);
        if (!f[12].empty()) {
          p_n++;
          const std::string theirs = f[12].substr(f[12].find(":p.") + 3);
          const std::string ours = pp.whole.empty() ? nm::format_variant(pp.variant) : pp.whole;
          if (ours == theirs) p_ok++; else note(f[1] + ": p." + ours + " vs VEP p." + theirs);
        }
        continue;
      }
      if (!t) { note("no transcript " + transcript); continue; }
      const nm::ProteinPrediction pp = nm::predict(rec.genome, *change, *t);
      c_n++;
      const std::string ours_c = nm::format_variant(pp.projection.variant), theirs_c = after(f[3], "c");
      if (ours_c == theirs_c) c_ok++; else note(f[1] + ": c." + ours_c + " vs Mutalyzer c." + theirs_c);
      // The protein against Mutalyzer, and against VEP where VEP describes one (VEP writes no
      // parentheses). Their conventions differ where the standard leaves a choice (theory §6):
      // Mutalyzer writes p.(=) for a synonymous change where VEP and the standard's example name the
      // residue, p.(His3=); p.? for a lost start where the standard writes p.(Met1?); nothing where
      // splicing may change, where the standard writes p.?.
      {
        p_n++;
        const std::string ours_p = pp.text().empty() ? std::string() : after(pp.text(), "p");
        const std::string mut = after(f[4], "p"), vep = f[12].empty() ? std::string() : "(" + f[12].substr(f[12].find(":p.") + 3) + ")";
        const bool synonymous = ours_p.size() > 2 && ours_p.substr(ours_p.size() - 2) == "=)" && mut == "(=)";
        const bool by_mutalyzer = ours_p == mut || synonymous || (ours_p == "(Met1?)" && mut == "?") || (ours_p == "?" && mut.empty());
        const bool by_vep = !vep.empty() && (ours_p == vep || (ours_p == "(Met1?)" && vep == "(Met1?)"));
        if (by_mutalyzer || by_vep) p_ok++;
        else note(f[1] + ": p." + ours_p + " vs Mutalyzer p." + mut + " and VEP " + vep);
        if (!by_mutalyzer) note("Mutalyzer differs on " + f[1] + ": p." + mut + " (ours p." + ours_p + ", VEP " + vep + ")");
      }
      if (!f[10].empty()) {
        vep_n++;
        const std::string most = f[10].substr(0, f[10].find(','));
        if (!pp.consequences.empty() && pp.consequences.front() == most) vep_ok++;
        else note(f[1] + ": " + (pp.consequences.empty() ? std::string("none") : pp.consequences.front()) + " vs VEP " + f[10]);
      }
    }
    if (!n) {
      s.skip("the variant fixture is absent (tests/fixtures/nomenclature-variants.tsv)");
    } else {
      s.require("real variants: the g. and m. descriptions are Mutalyzer's (the 3′ rule)", g_ok == n, std::to_string(g_ok) + " of " + std::to_string(n) + " (" + header + ")");
      s.require("their contextual SPDI is NCBI's (VOCA)", spdi_ok == n, std::to_string(spdi_ok) + " of " + std::to_string(n));
      s.require("their VCF fields are NCBI's", vcf_ok == n, std::to_string(vcf_ok) + " of " + std::to_string(n));
      s.require("projected through their transcripts, the c. descriptions are Mutalyzer's", c_ok == c_n, std::to_string(c_ok) + " of " + std::to_string(c_n));
      s.require("the predicted p. descriptions are Mutalyzer's or VEP's (mitochondrial: VEP's), up to their stated conventions", p_ok == p_n, std::to_string(p_ok) + " of " + std::to_string(p_n));
      s.require("the most severe consequence is VEP's", vep_ok == vep_n, std::to_string(vep_ok) + " of " + std::to_string(vep_n) + (bad.empty() ? "" : "; differences: " + bad));
    }
  }
  return s.done();
}
