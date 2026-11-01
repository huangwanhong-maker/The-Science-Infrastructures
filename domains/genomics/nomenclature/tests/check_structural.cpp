// Checks of structural.hpp on the VCF 4.3 specification's own examples (tests/fixtures/vcf43-*.*,
// made by make_vcf_sv_fixtures.py from VCFv4.3.tex at a pinned commit) and on the HGVS
// Nomenclature's descriptions of rearrangements (DNA/complex.md, in hgvs-examples.tsv):
//
// - §5.3's example: kinds, END, SVLEN, CIPOS, CIEND, CN read; described in HGVS with their uncertain
//   ends, g.(a_b)_(c_d)del (uncertain.md);
// - §5.4's breakend tables: every bracket form read and written back to the same ALT; a breakend and
//   its mate assert one adjacency; the inversion INV0 as <INV> and as its four breakends assert the
//   same two adjacencies, described as one inversion; telomeric breakends at POS 0; single breakends;
// - 2000 random adjacencies of every orientation, with inserted bases, written as breakend records
//   and read back;
// - translocations described as the derivative chromosome's delins: the standard's t(2;11)(p25.1;p15.2)
//   in both derivatives, its unbalanced der(2), and der(14) of t(3;14) with its inverted piece;
// - copy-number segments.

#include <genomics/nomenclature/structural.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace nm = genomics::nomenclature;

namespace {

std::string read_file(const std::string& p) {
  std::ifstream in(p, std::ios::binary);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

bool in_standard(const std::string& description) {
  const std::string t = read_file(std::string(GENOMICS_FIXTURES) + "/hgvs-examples.tsv");
  return t.find("\n" + description + "\t") != std::string::npos;
}

std::string upper(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- §5.3 ----------------------------------------------------------------------------------------------
  {
    std::vector<genomics::sequence::Problem> problems;
    const auto svs = nm::read_structural(read_file(std::string(GENOMICS_FIXTURES) + "/vcf43-structural.vcf"), &problems);
    s.require("VCF 4.3 §5.3: five symbolic records read (the precise deletion is a small variant)", svs.size() == 5 && problems.empty(), std::to_string(svs.size()));
    if (svs.size() == 5) {
      s.require("<DEL>: POS 321682, END 321887, SVLEN 205, CIPOS −56,20, CIEND −10,62",
                svs[0].kind == nm::SvKind::deletion && svs[0].pos == 321682 && svs[0].end == 321887 && svs[0].length == 205 && svs[0].cipos[0] == -56 && svs[0].ciend[1] == 62);
      s.require("<DEL:ME:ALU> is a deletion of subtype ME:ALU; <INS:ME:L1> an insertion of 6027 bases",
                svs[1].kind == nm::SvKind::deletion && svs[1].subtype == "ME:ALU" && svs[2].kind == nm::SvKind::insertion && svs[2].length == 6027);
      s.require("<DUP> with FORMAT CN 3; <DUP:TANDEM> a tandem duplication with CN 5",
                svs[3].kind == nm::SvKind::duplication && svs[3].copy_number == 3 && svs[4].kind == nm::SvKind::tandem_duplication && svs[4].copy_number == 5);
      const std::string d = nm::hgvs_of(svs[0]), t = nm::hgvs_of(svs[4]), u = nm::hgvs_of(svs[3]), i = nm::hgvs_of(svs[2]);
      s.require("the imprecise deletion in HGVS: 2:g.(321627_321703)_(321877_321949)del", d == "2:g.(321627_321703)_(321877_321949)del", d);
      s.require("the tandem duplication: 4:g.(18665119_18665139)_(18665194_18665214)dup", t == "4:g.(18665119_18665139)_(18665194_18665214)dup", t);
      s.require("the duplication with its copy number: 3:g.(12664601_12665601)_(12685700_12686700)[3]", u == "3:g.(12664601_12665601)_(12685700_12686700)[3]", u);
      s.require("the insertion of 6027 unknown bases: 3:g.(9425900_9425939)insN[6027]", i == "3:g.(9425900_9425939)insN[6027]", i);
      bool parse = true;
      for (const auto& x : {d, t, u, i}) parse = parse && nm::parse_hgvs("chr" + x).ok();  // the spec's contigs are 1 to 4; a reference has two characters at least
      s.require("each description (its contig written chrN) is in the HGVS grammar and meets its conditions", parse);
    }
  }

  // ---- §5.4 ---------------------------------------------------------------------------------------------------
  std::map<std::string, std::vector<std::vector<std::string>>> tables;
  {
    std::ifstream in(std::string(GENOMICS_FIXTURES) + "/vcf43-breakends.tsv");
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::vector<std::string> f;
      std::stringstream ss(line);
      std::string x;
      while (std::getline(ss, x, '\t')) f.push_back(x);
      tables[f[0]].push_back(f);
    }
  }
  const auto vcf_of = [](const std::vector<std::vector<std::string>>& rows) {
    std::string t = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    for (const auto& f : rows) t += f[1] + "\t" + f[2] + "\t" + f[3] + "\t" + f[4] + "\t" + f[5] + "\t" + f[6] + "\t" + f[7] + "\t" + f[8] + "\n";
    return t;
  };
  {
    std::size_t records = 0, exact = 0;
    std::string bad;
    for (const auto& [section, rows] : tables)
      for (const auto& f : rows) {
        if (f[5].find('<') != std::string::npos) continue;  // contigs of an assembly and the INS shorthand
        for (const auto& alt : std::vector<std::string>{f[5].substr(0, f[5].find(',')), f[5].find(',') == std::string::npos ? "" : f[5].substr(f[5].find(',') + 1)}) {
          if (alt.empty()) continue;
          const auto b = nm::parse_breakend(alt, f[1], std::stoll(f[2]), f[4]);
          records++;
          const std::string back = b.ok ? nm::breakend_alt(b.adjacency, f[4][0]) : "";
          if (b.ok && upper(back) == upper(alt)) exact++;
          else bad += section + " " + f[3] + " " + alt + " → " + back + (b.ok ? "" : " (" + b.error + ")") + "; ";
        }
      }
    s.require("§5.4: every bracket form of the specification read and written back to the same ALT", records > 30 && exact == records, std::to_string(exact) + " of " + std::to_string(records) + (bad.empty() ? "" : ": " + bad));
  }
  {
    const auto adj = nm::adjacencies(nm::read_structural(vcf_of(tables["Breakends"])));
    s.require("Figure 1's six breakends assert three adjacencies, each with its mate", adj.size() == 3, std::to_string(adj.size()));
    const auto ins = nm::adjacencies(nm::read_structural(vcf_of(tables["Inserted Sequence"])));
    s.require("the inserted sequence AGTNNNNNCA between 13:123456 and 2:321682, read the same from both mates",
              ins.size() == 1 && ins[0].inserted == "AGTNNNNNCA" && ins[0].a.contig == "13" && ins[0].a.position == 123455 && ins[0].b.position == 321681, ins.empty() ? "" : ins[0].inserted);
    std::vector<std::vector<std::string>> sym, bnd;
    for (const auto& f : tables["Inversions"]) (f[5] == "<INV>" ? sym : bnd).push_back(f);
    const auto a = nm::adjacencies(nm::read_structural(vcf_of(sym))), b = nm::adjacencies(nm::read_structural(vcf_of(bnd)));
    std::set<std::string> sa, sb;
    for (const auto& x : a) sa.insert(nm::breakend_alt(x, 'N'));
    for (const auto& x : b) sb.insert(nm::breakend_alt(x, 'N'));
    const auto da = nm::hgvs_of(a), db = nm::hgvs_of(b);
    s.require("INV0 as <INV> and as its four breakends asserts the same two adjacencies", a.size() == 2 && sa == sb, std::to_string(a.size()) + " and " + std::to_string(b.size()));
    s.require("both described as one inversion, 2:g.321682_421681inv", da.size() == 1 && db.size() == 1 && da[0] == "2:g.321682_421681inv" && db[0] == da[0], da.empty() ? "" : da[0]);
    std::vector<genomics::sequence::Problem> problems;
    const auto telo = nm::read_structural(vcf_of(tables["Telomeres"]), &problems);
    const auto ta = nm::adjacencies(telo);
    s.require("telomeric breakends at POS 0 are read (variant/vcf.hpp reads from POS 1), four records two adjacencies",
              telo.size() == 4 && telo[0].pos == 0 && telo[0].adjacency.a.position == -1 && ta.size() == 2 && problems.empty(), std::to_string(telo.size()) + " " + std::to_string(ta.size()));
    const auto single = nm::read_structural(vcf_of(tables["Single breakends"]));
    std::size_t singles = 0;
    for (const auto& x : single) singles += x.kind == nm::SvKind::single_breakend;
    s.require("single breakends (G., .A, .TGCA, TCC.) are read as single", singles == 6, std::to_string(singles));
  }

  // ---- random adjacencies ---------------------------------------------------------------------------------------
  {
    calx::random::Stream r(20261101, 0x42524b4e, 0, 0);
    std::size_t ok = 0, n = 0;
    for (int k = 0; k < 2000; k++) {
      nm::Adjacency x;
      x.a = {std::string("chr") + std::to_string(1 + r.below(3)), static_cast<std::int64_t>(r.below(100000)), r.below(2) ? '+' : '-'};
      x.b = {std::string("chr") + std::to_string(1 + r.below(3)), static_cast<std::int64_t>(r.below(100000)), r.below(2) ? '+' : '-'};
      for (std::size_t i = 0, L = r.below(5); i < L; i++) x.inserted += "ACGT"[r.below(4)];
      const auto recs = nm::breakend_records(x, "bnd" + std::to_string(k), [](const std::string&, std::int64_t p) { return "ACGT"[p % 4]; });
      std::string text = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
      for (const auto& rec : recs) text += genomics::variant::vcf_line(rec);
      const auto back = nm::adjacencies(nm::read_structural(text));
      n++;
      ok += back.size() == 1 && back[0] == nm::canonical(x);
    }
    s.require("2000 random adjacencies of every orientation, with inserted bases, written as two mates and read back as one", ok == n, std::to_string(ok) + " of " + std::to_string(n));
  }

  // ---- HGVS of rearrangements ------------------------------------------------------------------------------------------
  {
    nm::SvHgvsOptions o;
    const std::map<std::string, std::string> acc = {{"chr2", "NC_000002.12"}, {"chr11", "NC_000011.10"}, {"chr3", "NC_000003.12"}, {"chr14", "NC_000014.9"}};
    o.accession = [acc](const std::string& c) { return acc.at(c); };
    // A base within each centromere (GRCh38; chromosome 14 is acrocentric, its centromere near the p-arm's end).
    const std::map<std::string, std::int64_t> cen = {{"chr2", 93000000}, {"chr11", 53000000}, {"chr3", 91000000}, {"chr14", 17500000}};
    o.centromere = [cen](const std::string& c) -> std::optional<std::int64_t> { return cen.at(c); };
    const auto bk = [](const char* contig, std::int64_t one_based, char strand) { return nm::Breakend{contig, one_based - 1, strand}; };
    const std::string der2 = nm::hgvs_of(nm::Adjacency{bk("chr11", 15825272, '+'), bk("chr2", 8247757, '-'), "", false}, o);
    const std::string der11 = nm::hgvs_of(nm::Adjacency{bk("chr2", 8247756, '+'), bk("chr11", 15825273, '-'), "", false}, o);
    s.require("t(2;11)(p25.1;p15.2), der(2): NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825272] (complex.md)",
              der2 == "NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825272]" && in_standard(der2), der2);
    s.require("its der(11): NC_000011.10:g.pter_15825272delins[NC_000002.12:g.pter_8247756]", der11 == "NC_000011.10:g.pter_15825272delins[NC_000002.12:g.pter_8247756]" && in_standard(der11), der11);
    const std::string unb = nm::hgvs_of(nm::Adjacency{bk("chr2", 8247757, '-'), bk("chr11", 15825266, '+'), "", false}, o);
    s.require("the unbalanced der(2)t(2;11): NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825266]", unb == "NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825266]" && in_standard(unb), unb);
    const std::string der14 = nm::hgvs_of(nm::Adjacency{bk("chr14", 29745313, '+'), bk("chr3", 36969141, '+'), "", false}, o);
    s.require("der(14) of t(3;14): the inverted piece of chromosome 3, its ends written 5′ to 3′ (the standard writes 36969141_pterinv)",
              der14 == "NC_000014.9:g.29745314_qterdelins[NC_000003.12:g.pter_36969141inv]" && in_standard("NC_000014.9:g.29745314_qterdelins[NC_000003.12:g.36969141_pterinv]"), der14);
    const std::string del = nm::hgvs_of(nm::Adjacency{bk("chr2", 1000, '+'), bk("chr2", 2001, '-'), "", false}, o);
    const std::string dup = nm::hgvs_of(nm::Adjacency{bk("chr2", 2000, '+'), bk("chr2", 1001, '-'), "", false}, o);
    const std::string ins = nm::hgvs_of(nm::Adjacency{bk("chr2", 1000, '+'), bk("chr2", 1001, '-'), "ACGT", false}, o);
    s.require("facing breakends on one contig: a deletion, a tandem duplication, an insertion",
              del == "NC_000002.12:g.1001_2000del" && dup == "NC_000002.12:g.1001_2000dup" && ins == "NC_000002.12:g.1000_1001insACGT", del + " " + dup + " " + ins);
    bool parse = true;
    for (const auto& x : {der2, der11, unb, der14, del, dup, ins}) parse = parse && nm::parse_hgvs(x).accepted;
    s.require("every description is in the HGVS grammar", parse);
  }

  // ---- copy number ------------------------------------------------------------------------------------------------------
  {
    const auto h = [](int cn) { return nm::hgvs_of(nm::CopyNumberSegment{"chr2", 999, 2000, cn, 2}); };
    s.require("copy-number segments: 0 copies a homozygous deletion, 1 a deletion, 2 no change, 3 three copies",
              h(0) == "chr2:g.[1000_2000del];[1000_2000del]" && h(1) == "chr2:g.1000_2000del" && h(2) == "chr2:g.1000_2000=" && h(3) == "chr2:g.1000_2000[3]", h(0) + " " + h(1) + " " + h(3));
    const auto rec = nm::symbolic_record(nm::SvKind::deletion, "chr2", 999, 2000, 'A');
    std::string text = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + genomics::variant::vcf_line(rec);
    const auto back = nm::read_structural(text);
    s.require("a symbolic <DEL> written and read back: POS 999, END 2000, SVLEN 1001", back.size() == 1 && back[0].pos == 999 && back[0].end == 2000 && back[0].length == 1001,
              genomics::variant::vcf_line(rec));
  }
  return s.done();
}
