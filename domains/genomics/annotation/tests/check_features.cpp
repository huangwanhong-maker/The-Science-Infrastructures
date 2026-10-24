// Checks of features as tracks: a GenBank record of a circular genome written here, with a gene
// and its CDS across the origin, a feature across the origin on the reverse strand, a CDS on the
// complement strand, a spliced CDS whose /translation is wrong, and an RNA, a repeat and a mobile
// element of each kind, read into tracks with their kinds, names, segments, extents, strands
// and translations, and lanes in which nothing overlaps; and GFF3: the specification's canonical
// gene (a fixture transcribed from version 1.26) read with its part-of hierarchy, its
// discontinuous CDSs and phases, percent-encoding decoded, the specification's circular example
// split at the origin, a ##FASTA section, and each fault found at its line. When the data cache
// has NCBI's genomes, every record's CDSs are drawn and their translations agree.

#include <genomics/annotation/features.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::annotation;
using sequence::Feature;
using sequence::GenBank;

namespace {

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

Codes random_cds(calx::random::Stream& r, std::size_t codons, const sequence::GeneticCode& g) {
  Codes x = sequence::encode_dna("ATG");
  while (x.size() < 3 * (codons - 1)) {
    const Code a = static_cast<Code>(r.below(4)), b = static_cast<Code>(r.below(4)), c = static_cast<Code>(r.below(4));
    if (!sequence::is_stop(g, a, b, c)) { x.push_back(a); x.push_back(b); x.push_back(c); }
  }
  const Codes stop = sequence::encode_dna("TAA");
  x.insert(x.end(), stop.begin(), stop.end());
  return x;
}

std::string qualifier(const std::string& name, const std::string& value, bool quoted = true) {
  const std::string text = "/" + name + "=" + (quoted ? "\"" + value + "\"" : value);
  std::string out;
  for (std::size_t i = 0; i < text.size(); i += 58) out += std::string(21, ' ') + text.substr(i, 58) + "\n";
  return out;
}

std::string feature(const std::string& key, const std::string& location) {
  return "     " + key + std::string(16 - key.size(), ' ') + location + "\n";
}

// The intervals a track feature covers.
std::vector<std::pair<std::size_t, std::size_t>> cover(const TrackFeature& f, std::size_t length) {
  if (f.wraps) return {{f.start, length}, {0, f.stop}};
  return {{f.start, f.stop}};
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x66656174, 0, 0);

  // ---- a GenBank record -------------------------------------------------------------------------------------
  {
    const auto& g11 = sequence::genetic_code(11);
    const std::size_t n = 3000;
    Codes genome(n);
    for (auto& c : genome) c = static_cast<Code>(r.below(4));
    const Codes cds1 = random_cds(r, 101, g11);  // join(2901..3000,1..203)
    for (std::size_t i = 0; i < 100; i++) genome[2900 + i] = cds1[i];
    for (std::size_t i = 100; i < cds1.size(); i++) genome[i - 100] = cds1[i];
    const Codes cds2 = random_cds(r, 100, g11), rc2 = sequence::reverse_complement(cds2);  // complement(1001..1300)
    std::copy(rc2.begin(), rc2.end(), genome.begin() + 1000);
    const Codes cds3 = random_cds(r, 100, g11);  // join(1501..1560,1601..1840)
    for (std::size_t i = 0; i < 60; i++) genome[1500 + i] = cds3[i];
    for (std::size_t i = 60; i < cds3.size(); i++) genome[1600 + i - 60] = cds3[i];
    const std::string p1 = sequence::translate(cds1, g11, 0, true), p2 = sequence::translate(cds2, g11, 0, true);
    std::string p3 = sequence::translate(cds3, g11, 0, true);
    p3[10] = p3[10] == 'W' ? 'A' : 'W';  // a /translation that the sequence does not give

    std::string gb = "LOCUS       TEST0002                3000 bp    DNA     circular PHG 24-OCT-2026\n";
    gb += "DEFINITION  A test genome with features across the origin.\nACCESSION   TEST0002\nVERSION     TEST0002.1\n";
    gb += "FEATURES             Location/Qualifiers\n";
    gb += feature("source", "1..3000") + qualifier("organism", "synthetic construct");
    gb += feature("gene", "join(2901..3000,1..203)") + qualifier("gene", "one");
    gb += feature("CDS", "join(2901..3000,1..203)") + qualifier("gene", "one") + qualifier("transl_table", "11", false) + qualifier("translation", p1);
    gb += feature("gene", "complement(1001..1300)") + qualifier("locus_tag", "T0002");
    gb += feature("CDS", "complement(1001..1300)") + qualifier("locus_tag", "T0002") + qualifier("transl_table", "11", false) + qualifier("translation", p2);
    const std::size_t spliced_line = static_cast<std::size_t>(std::count(gb.begin(), gb.end(), '\n')) + 1;
    gb += feature("CDS", "join(1501..1560,1601..1840)") + qualifier("product", "spliced protein") + qualifier("transl_table", "11", false) + qualifier("translation", p3);
    gb += feature("tRNA", "2001..2076") + qualifier("product", "tRNA-Gly");
    gb += feature("rRNA", "2101..2400") + qualifier("product", "16S ribosomal RNA");
    gb += feature("repeat_region", "2450..2500") + qualifier("rpt_family", "REP");
    gb += feature("mobile_element", "2501..2800") + qualifier("mobile_element_type", "insertion sequence:IS1");
    gb += feature("misc_feature", "complement(join(2971..3000,1..30))") + qualifier("note", "a site across the origin on the reverse strand");
    gb += "ORIGIN\n";
    const std::string letters = sequence::decode_dna(genome);
    for (std::size_t i = 0; i < n; i += 60) {
      std::string line = std::to_string(i + 1);
      line = std::string(9 - line.size(), ' ') + line;
      for (std::size_t j = i; j < std::min(n, i + 60); j += 10) line += " " + letters.substr(j, 10);
      gb += line + "\n";
    }
    gb += "//\n";

    std::vector<sequence::Problem> read_problems, problems;
    const auto recs = sequence::read_genbank(gb, &read_problems);
    s.require("the record is read without problems", recs.size() == 1 && read_problems.empty(), read_problems.empty() ? "" : read_problems[0].message);
    if (recs.size() == 1) {
      const Tracks t = feature_tracks(recs[0], {}, &problems);
      const auto find = [&](Kind k, const std::string& name) -> const TrackFeature* {
        for (const auto& f : t.features)
          if (f.kind == k && f.name == name) return &f;
        return nullptr;
      };
      s.require("ten features (the source left out), of their kinds", t.features.size() == 10 && t.counts[static_cast<std::size_t>(Kind::gene)] == 2 && t.counts[static_cast<std::size_t>(Kind::cds)] == 3 &&
                                                                     t.counts[static_cast<std::size_t>(Kind::trna)] == 1 && t.counts[static_cast<std::size_t>(Kind::rrna)] == 1 && t.counts[static_cast<std::size_t>(Kind::repeat)] == 1 &&
                                                                     t.counts[static_cast<std::size_t>(Kind::mobile)] == 1 && t.counts[static_cast<std::size_t>(Kind::misc)] == 1);
      const TrackFeature* one = find(Kind::cds, "one");
      s.require("a CDS across the origin: it wraps, from 2900 round to 203, in two segments of 303 bases, its translation agreeing",
                one && one->wraps && one->start == 2900 && one->stop == 203 && one->bases == 303 && one->segments.size() == 2 && one->segments[0].begin == 0 &&
                    one->segments[0].end == 203 && one->segments[1].begin == 2900 && one->segments[1].end == 3000 && !one->reverse && one->translation == Translation::agrees);
      const TrackFeature* two = find(Kind::cds, "T0002");
      s.require("a CDS on the complement strand, named by its locus tag, its translation agreeing",
                two && two->reverse && !two->wraps && two->start == 1000 && two->stop == 1300 && two->translation == Translation::agrees);
      const TrackFeature* three = find(Kind::cds, "spliced protein");
      s.require("a spliced CDS named by its product: two segments, its wrong /translation found and reported at its line",
                three && three->segments.size() == 2 && three->segments[0].end == 1560 && three->segments[1].begin == 1600 && three->translation == Translation::differs &&
                    problems.size() == 1 && problems[0].line == spliced_line,
                problems.empty() ? "no problem reported" : problems[0].message + " (line " + std::to_string(problems[0].line) + ")");
      const TrackFeature* site = nullptr;
      for (const auto& f : t.features)
        if (f.kind == Kind::misc) site = &f;
      s.require("a feature across the origin on the reverse strand wraps, from 2970 round to 30", site && site->wraps && site->reverse && site->start == 2970 && site->stop == 30 && site->bases == 60);
      s.require("names from /product, /rpt_family and /mobile_element_type", find(Kind::trna, "tRNA-Gly") && find(Kind::rrna, "16S ribosomal RNA") && find(Kind::repeat, "REP") &&
                                                                                find(Kind::mobile, "insertion sequence:IS1") && find(Kind::gene, "one") && find(Kind::gene, "T0002"));
      bool apart = true;
      for (std::size_t a = 0; a < t.features.size(); a++)
        for (std::size_t b = a + 1; b < t.features.size(); b++) {
          const auto& x = t.features[a];
          const auto& y = t.features[b];
          if (x.kind != y.kind || x.lane != y.lane) continue;
          for (const auto& [p, q] : cover(x, t.length))
            for (const auto& [u, v] : cover(y, t.length))
              if (std::min(q, v) > std::max(p, u)) apart = false;
        }
      s.require("features of a kind in one lane never overlap, across the origin included", apart);
      s.require("the CDSs between the ends of the one across the origin share its lane: one lane for the three",
                one && two && three && one->lane == 0 && two->lane == 0 && three->lane == 0 && t.lanes[static_cast<std::size_t>(Kind::cds)] == 1);
    }
  }

  // ---- GFF3: the canonical gene ------------------------------------------------------------------------------
  {
    const std::string text = read_file(std::filesystem::path(GENOMICS_FIXTURES) / "gff3-canonical-gene.gff3");
    std::vector<sequence::Problem> problems;
    const Gff g = read_gff3(text, &problems);
    s.require("the canonical gene of GFF3 1.26 is read without problems: 23 lines, version 3.1.26", !text.empty() && problems.empty() && g.features.size() == 23 && g.version == "3.1.26",
              problems.empty() ? "" : problems[0].message + " at line " + std::to_string(problems[0].line));
    s.require("##sequence-region ctg123 1 1497228 is [0, 1497228)", g.regions.size() == 1 && g.regions[0].seqid == "ctg123" && g.regions[0].begin == 0 && g.regions[0].end == 1497228);
    const auto entity = [&](const std::string& id) -> const GffEntity* {
      for (const auto& e : g.entities)
        if (e.id == id) return &e;
      return nullptr;
    };
    const GffEntity* gene = entity("gene00001");
    s.require("fourteen features: the gene, a binding site, three mRNAs, five exons, four CDSs; the gene the only root",
              g.entities.size() == 14 && g.roots.size() == 1 && gene && g.roots[0] == static_cast<std::size_t>(gene - g.entities.data()));
    s.require("the gene's parts are the binding site and the three mRNAs; EDEN is its name", gene && gene->children.size() == 4 && g.features[gene->lines[0]].name == "EDEN");
    const GffEntity* m3 = entity("mRNA00003");
    const GffEntity* e4 = entity("exon00004");
    s.require("mRNA00003 has four exons and two CDSs (its two starts); exon00004 is part of all three mRNAs", m3 && m3->children.size() == 6 && e4 && e4->parents.size() == 3);
    const GffEntity* c1 = entity("cds00001");
    const GffEntity* c3 = entity("cds00003");
    bool phases = c3 && c3->lines.size() == 3 && g.features[c3->lines[0]].phase == 0 && g.features[c3->lines[1]].phase == 1 && g.features[c3->lines[2]].phase == 1;
    s.require("a discontinuous CDS is one feature of four lines over [1200, 7600); cds00003's phases are 0, 1, 1", c1 && c1->lines.size() == 4 && c1->begin == 1200 && c1->end == 7600 && phases);
    s.require("coordinates are read as [start − 1, end): the binding site 1000..1012 is [999, 1012)", g.features[1].begin == 999 && g.features[1].end == 1012 && g.features[1].strand == '+');
  }

  // ---- GFF3: encoding, circular landmarks, FASTA and faults -------------------------------------------------------
  {
    const std::string text =
        "##gff-version 3\n"
        "chr%201\tsrc\tgene\t10\t20\t0.5\t-\t.\tID=g1;Note=a%3Bb%2Cc%25,second;Name=tab%09here\n"
        "##FASTA\n>chr 1\nACGT\nACGT\n>other\nGG\n";
    std::vector<sequence::Problem> problems;
    const Gff g = read_gff3(text, &problems);
    const auto* note = g.features.empty() ? nullptr : g.features[0].attribute("Note");
    s.require("percent-encoding decoded in every column, commas separating values and encoded ones not",
              problems.empty() && g.features.size() == 1 && g.features[0].seqid == "chr 1" && note && note->size() == 2 && (*note)[0] == "a;b,c%" && g.features[0].name == "tab\there" &&
                  g.features[0].score == 0.5 && g.features[0].strand == '-');
    s.require("a ##FASTA section read as FASTA", g.fasta.size() == 2 && g.fasta[0].sequence == "ACGTACGT" && g.fasta[1].id == "other");

    // The specification's circular example (bacteriophage f1, gene II across the origin).
    const std::string f1 =
        "##gff-version 3.1.26\n# organism Enterobacteria phage f1\n##sequence-region J02448 1 6407\n"
        "J02448\tGenBank\tregion\t1\t6407\t.\t+\t.\tID=J02448;Name=J02448;Is_circular=true;\n"
        "J02448\tGenBank\tCDS\t6006\t7238\t.\t+\t0\tID=geneII;Name=II;Note=protein II;\n";
    problems.clear();
    const Gff c = read_gff3(f1, &problems);
    const auto seg = c.entities.size() == 2 ? segments(c, c.entities[1], 6407) : std::vector<Segment>{};
    s.require("a circular landmark's feature past its end is no fault, and is split at the origin into [6005, 6407) and [0, 831)",
              problems.empty() && c.circular.size() == 1 && seg.size() == 2 && seg[0].begin == 0 && seg[0].end == 831 && seg[1].begin == 6005 && seg[1].end == 6407,
              problems.empty() ? "" : problems[0].message);

    const std::string bad =
        "##gff-version 3\n"                                  // 1
        "##sequence-region s 1 100\n"                        // 2
        "s\t.\tgene\t1\t50\t.\t+\t.\tID=a;Parent=b\n"        // 3
        "s\t.\tgene\t1\t50\t.\t+\t.\tID=b;Parent=a\n"        // 4: a cycle a → b → a, found back at a (line 3)
        "s\t.\tgene\t60\t50\t.\t+\t.\tID=c\n"                // 5: start after end
        "s\t.\tCDS\t1\t30\t.\t+\t.\tID=d\n"                  // 6: a CDS without phase
        "s\t.\tgene\t1\t30\t.\t*\t.\tID=e\n"                 // 7: a bad strand
        "s\t.\tgene\t1\t30\t.\t+\t.\tID=f;Parent=nobody\n"   // 8: a Parent that is no ID
        "s\t.\tgene\t1\t30\t.\t+\t.\tID=g%2\n"               // 9: a broken percent-encoding
        "s\t.\tgene\t90\t130\t.\t+\t.\tID=h\n"               // 10: outside the sequence region
        "s\t.\tgene\t1\t30\n";                               // 11: four columns
    problems.clear();
    read_gff3(bad, &problems);
    std::vector<std::size_t> lines;
    for (const auto& p : problems) lines.push_back(p.line);
    std::sort(lines.begin(), lines.end());
    const std::vector<std::size_t> expected{3, 5, 6, 7, 8, 9, 10, 11};
    std::string detail;
    for (const auto& p : problems) detail += std::to_string(p.line) + ": " + p.message + "; ";
    s.require("faults found at their lines: a Parent cycle, start after end, a CDS without phase, a bad strand, an unknown Parent, broken encoding, outside the region, too few columns",
              lines == expected, detail);
  }

  // ---- the genomes of NCBI, when the data library has fetched them ------------------------------------------------
  {
    const char* env = std::getenv("SCIENCE_CACHE");
    const std::filesystem::path cache = env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
    const std::filesystem::path dir = cache / "ncbi-genomes";
    std::vector<std::filesystem::path> files;
    if (std::filesystem::is_directory(dir))
      for (const auto& e : std::filesystem::directory_iterator(dir))
        if (e.path().extension() == ".gb") files.push_back(e.path());
    std::sort(files.begin(), files.end());
    if (files.empty()) s.skip("the NCBI genomes are not in the data cache (data/cache/ncbi-genomes): their tracks are not drawn");
    for (const auto& file : files) {
      const auto recs = sequence::read_genbank(read_file(file));
      std::size_t cds = 0, agree = 0, wraps = 0, features = 0;
      std::vector<sequence::Problem> problems;
      for (const auto& rec : recs) {
        const Tracks t = feature_tracks(rec, {}, &problems);
        features += t.features.size();
        for (const auto& f : t.features) {
          wraps += f.wraps;
          if (f.kind != Kind::cds || f.translation == Translation::none || f.translation == Translation::unchecked) continue;
          cds++;
          agree += f.translation == Translation::agrees;
        }
      }
      s.require("NCBI " + file.stem().string() + ": " + std::to_string(features) + " features drawn (" + std::to_string(wraps) + " across the origin), every checked CDS's translation agreeing",
                cds > 0 && agree == cds, std::to_string(agree) + " of " + std::to_string(cds));
    }
  }
  return s.done();
}
