// Checks of the readers and writers: FASTA with comments, wrapped lines, descriptions and
// CRLF endings, read back after writing; FASTQ in four lines and in many, its faults found
// where they are; the Phred scale; the location grammar of the feature table (spans, partial
// ends, between-base sites, complement, join, order, remote parts, and its errors); and a
// GenBank record of a circular genome written here with genes across the origin and on the
// complement strand, read whole: LOCUS, a two-line DEFINITION, the lineage, a reference whose
// continuation lines stay out of the definition, every qualifier (quoted over lines, with
// doubled quotation marks, bare), the ORIGIN, and every CDS translated as its /translation.
// When the data library has fetched the genomes of NCBI (data/cache/ncbi-genomes), every CDS
// of each is translated against GenBank's own /translation.

#include <genomics/sequence/fasta.hpp>
#include <genomics/sequence/genbank.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics::sequence;

namespace {

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

// A random coding sequence of `codons` codons: ATG, codons without a stop, TAA.
Codes random_cds(calx::random::Stream& r, std::size_t codons, const GeneticCode& g) {
  Codes x = encode_dna("ATG");
  while (x.size() < 3 * (codons - 1)) {
    const Code a = static_cast<Code>(r.below(4)), b = static_cast<Code>(r.below(4)), c = static_cast<Code>(r.below(4));
    if (!is_stop(g, a, b, c)) { x.push_back(a); x.push_back(b); x.push_back(c); }
  }
  const Codes stop = encode_dna("TAA");
  x.insert(x.end(), stop.begin(), stop.end());
  return x;
}

// A qualifier's value wrapped as GenBank wraps it, from column 22, 58 characters a line.
std::string qualifier(const std::string& name, const std::string& value, bool quoted = true) {
  std::string text = "/" + name + "=" + (quoted ? "\"" + value + "\"" : value);
  std::string out;
  std::size_t i = 0;
  while (i < text.size()) {
    // A translation is cut anywhere, text at a space (dropped, as the reader restores it).
    std::size_t len = std::min<std::size_t>(58, text.size() - i), next = i + len;
    if (name != "translation" && i + len < text.size()) {
      const std::size_t sp = text.rfind(' ', i + len);
      if (sp != std::string::npos && sp > i) { len = sp - i; next = sp + 1; }
    }
    out += std::string(21, ' ') + text.substr(i, len) + "\n";
    i = next;
  }
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x666d7421, 0, 0);

  // ---- FASTA ----------------------------------------------------------------------------------
  {
    const std::string text = ";a comment\r\n>seq1 the first record\r\nACGT\r\nacgtNN\r\n\r\n>seq2\nMKV LL\n>seq3 empty\n";
    std::vector<Problem> problems;
    const auto rec = read_fasta(text, &problems);
    s.require("FASTA: three records, their identifiers and descriptions", rec.size() == 3 && rec[0].id == "seq1" && rec[0].description == "the first record" && rec[1].id == "seq2" && rec[1].description.empty());
    s.require("FASTA: wrapped lines joined, white space and CRLF removed, comments skipped", rec[0].sequence == "ACGTacgtNN" && rec[1].sequence == "MKVLL" && rec[2].sequence.empty() && problems.empty());
    const auto back = read_fasta(write_fasta(rec, 4));
    s.require("FASTA: written and read again, the same records", back.size() == 3 && back[0].sequence == rec[0].sequence && back[0].description == rec[0].description);
    read_fasta("ACGT\n>x\nAC\n", &problems);
    s.require("FASTA: a sequence before the first header is a problem on its line", problems.size() == 1 && problems[0].line == 1);
  }

  // ---- FASTQ ----------------------------------------------------------------------------------
  {
    const std::string text = "@r1 one\nACGT\n+\nIIII\n@r2\nAC\nGT\n+r2\nII\nI#\n@r3\n\n+\n\n";
    std::vector<Problem> problems;
    const auto rec = read_fastq(text, &problems);
    s.require("FASTQ: records in four lines and in many, and an empty one", rec.size() == 3 && rec[0].quality == "IIII" && rec[1].sequence == "ACGT" && rec[1].quality == "III#" && rec[2].sequence.empty() && problems.empty());
    const auto back = read_fastq(write_fastq(rec));
    s.require("FASTQ: written and read again", back.size() == 3 && back[1].quality == "III#");
    problems.clear();
    read_fastq("@a\nACGT\n+b\nIIIIII\n@c\nAC\n+\nI", &problems);
    s.require("FASTQ: a '+' naming another record, a quality too long, a record cut short: three problems", problems.size() == 3, std::to_string(problems.size()));
    s.require("Phred: Q 30 is an error of 0.001, written '?'", std::fabs(error_probability(30) - 0.001) < 1e-15 && phred_char(30) == '?' && phred_of('?') == 30 && phred_from_probability(0.001) == 30);
  }

  // ---- a line at a time ------------------------------------------------------------------------
  {
    // Each text fed to the readers line by line as std::getline gives them (a \r\n line keeps its
    // '\r'), against the whole-text readers: the same records and the same problems on the same lines.
    const auto same = [](const std::vector<Record>& a, const std::vector<Record>& b) {
      if (a.size() != b.size()) return false;
      for (std::size_t i = 0; i < a.size(); i++)
        if (a[i].id != b[i].id || a[i].description != b[i].description || a[i].sequence != b[i].sequence || a[i].quality != b[i].quality) return false;
      return true;
    };
    const auto same_problems = [](const std::vector<Problem>& a, const std::vector<Problem>& b) {
      if (a.size() != b.size()) return false;
      for (std::size_t i = 0; i < a.size(); i++)
        if (a[i].line != b[i].line || a[i].message != b[i].message) return false;
      return true;
    };
    const auto by_lines = [](const std::string& text, bool fastq, std::vector<Problem>& problems) {
      std::vector<Record> out;
      std::istringstream in(text);
      std::string line;
      Record r;
      FastaReader fa(&problems);
      FastqReader fq(&problems);
      while (std::getline(in, line))
        if (fastq ? fq.line(line, r) : fa.line(line, r)) out.push_back(std::move(r));
      if (fastq) fq.finish();
      else if (fa.finish(r)) out.push_back(std::move(r));
      return out;
    };
    // Random records, wrapped, with comments, CRLF endings, empty records and faults.
    std::string fasta = "ACGT\n", fastq;
    for (int k = 0; k < 200; k++) {
      std::string seq;
      const std::size_t n = r.below(300);
      for (std::size_t i = 0; i < n; i++) seq += "ACGTN"[r.below(5)];
      const std::string end = r.below(4) == 0 ? "\r\n" : "\n";
      fasta += (r.below(10) == 0 ? ";comment" + end : "") + ">r" + std::to_string(k) + (r.below(2) ? " a description" : "") + end;
      for (std::size_t i = 0; i < seq.size(); i += 60) fasta += seq.substr(i, 60) + end;
      std::string q;
      for (std::size_t i = 0; i < seq.size(); i++) q += phred_char(static_cast<int>(r.below(41)));
      const bool many = r.below(3) == 0;
      fastq += "@r" + std::to_string(k) + end + (many ? seq.substr(0, seq.size() / 2) + end + seq.substr(seq.size() / 2) : seq) + end + (r.below(5) == 0 ? "+r" + std::to_string(k) : "+") + end;
      fastq += (many ? q.substr(0, q.size() / 3) + end + q.substr(q.size() / 3) : q) + end;
      if (k == 70) fastq += "not a header" + end;
      if (k == 120) fastq += "@x\nACGT\n+y\nIIIIII" + end;
    }
    fastq += "@cut\nACGT\n+\nII";
    std::vector<Problem> p1, p2, q1, q2;
    const auto whole_a = read_fasta(fasta, &p1), lines_a = by_lines(fasta, false, p2);
    const auto whole_q = read_fastq(fastq, &q1), lines_q = by_lines(fastq, true, q2);
    s.require("FASTA a line at a time: the same 200 records as the whole text, and the same problem on the same line", whole_a.size() == 200 && same(whole_a, lines_a) && same_problems(p1, p2) && p1.size() == 1,
              std::to_string(lines_a.size()));
    s.require("FASTQ a line at a time: the same records as the whole text, in four lines and in many, and the same problems on the same lines", whole_q.size() == 201 && same(whole_q, lines_q) && same_problems(q1, q2) && q1.size() == 4,
              std::to_string(lines_q.size()) + " records, " + std::to_string(q1.size()) + " problems");
  }

  // ---- locations ----------------------------------------------------------------------------------
  {
    Location L;
    std::string e;
    s.require("a span 10..20 is [9, 20) forward", parse_location("10..20", L, e) && L.parts.size() == 1 && L.parts[0].begin == 9 && L.parts[0].end == 20 && !L.parts[0].reverse && L.kind == "span");
    s.require("a single base 7 is [6, 7)", parse_location("7", L, e) && L.parts[0].begin == 6 && L.parts[0].end == 7);
    s.require("partial ends <1..>200", parse_location("<1..>200", L, e) && L.parts[0].partial_lower && L.parts[0].partial_upper && L.parts[0].end == 200);
    s.require("complement(join(10..20,30..40)) reads 30..40 then 10..20, both reversed",
              parse_location("complement(join(10..20,30..40))", L, e) && L.parts.size() == 2 && L.parts[0].begin == 29 && L.parts[0].reverse && L.parts[1].begin == 9 && L.parts[1].reverse && L.kind == "join" && L.reverse());
    s.require("join(complement(1..5),10..20) mixes strands in order", parse_location("join(complement(1..5),10..20)", L, e) && L.parts[0].reverse && !L.parts[1].reverse && L.length() == 16);
    s.require("order(1..3, 7..9) keeps its kind", parse_location("order(1..3, 7..9)", L, e) && L.kind == "order" && L.parts.size() == 2);
    s.require("123^124 lies between two bases: no bases, uncertain", parse_location("123^124", L, e) && L.parts[0].begin == L.parts[0].end && L.uncertain);
    s.require("a remote part, J00194.1:100..202, is noted and left out", parse_location("join(1..10,J00194.1:100..202)", L, e) && L.remote && L.parts.size() == 1);
    s.require("errors: an unclosed join, a reversed span, a span beyond the sequence",
              !parse_location("join(1..5", L, e) && !parse_location("9..3", L, e) && !parse_location("1..600", L, e, 500));
  }

  // ---- a GenBank record of a circular genome -------------------------------------------------------------
  {
    const auto& g = genetic_code(11);
    const std::size_t n = 3000;
    Codes genome(n);
    for (auto& c : genome) c = static_cast<Code>(r.below(4));
    // CDS 1 across the origin: 2901..3000 then 1..203 (303 bases, 101 codons).
    const Codes cds1 = random_cds(r, 101, g);
    for (std::size_t i = 0; i < 100; i++) genome[2900 + i] = cds1[i];
    for (std::size_t i = 100; i < cds1.size(); i++) genome[i - 100] = cds1[i];
    // CDS 2 on the complement strand, 1001..1300.
    const Codes cds2 = random_cds(r, 100, g);
    const Codes rc2 = reverse_complement(cds2);
    for (std::size_t i = 0; i < rc2.size(); i++) genome[1000 + i] = rc2[i];
    // CDS 3 starting with GTG, spliced: join(1501..1560,1601..1840).
    Codes cds3 = random_cds(r, 100, g);
    cds3[0] = G;
    for (std::size_t i = 0; i < 60; i++) genome[1500 + i] = cds3[i];
    for (std::size_t i = 60; i < cds3.size(); i++) genome[1600 + i - 60] = cds3[i];
    const std::string p1 = translate(cds1, g, 0, true), p2 = translate(cds2, g, 0, true), p3 = translate(cds3, g, 0, true);

    std::string gb;
    gb += "LOCUS       TEST0001                3000 bp    DNA     circular PHG 24-OCT-2026\n";
    gb += "DEFINITION  A test genome, written by the checks of the genomics domain, with\n";
    gb += "            three genes.\n";
    gb += "ACCESSION   TEST0001\nVERSION     TEST0001.1\nKEYWORDS    synthetic; test\n            more keywords.\n";
    gb += "SOURCE      synthetic construct\n  ORGANISM  synthetic construct\n            other sequences; artificial sequences.\n";
    gb += "REFERENCE   1  (bases 1 to 3000)\n  AUTHORS   Nobody,A. and\n            Somebody,B.\n  TITLE     Direct Submission\n";
    gb += "FEATURES             Location/Qualifiers\n";
    gb += "     source          1..3000\n" + qualifier("organism", "synthetic construct") + qualifier("mol_type", "genomic DNA");
    gb += "     CDS             join(2901..3000,1..203)\n" + qualifier("gene", "one") + qualifier("codon_start", "1", false) + qualifier("transl_table", "11", false) + qualifier("translation", p1);
    gb += "     CDS             complement(1001..1300)\n" + qualifier("gene", "two") + qualifier("note", "a note that says \"\"quoted\"\" words and runs over more than one line of the table") + qualifier("transl_table", "11", false) + qualifier("translation", p2);
    gb += "     CDS             join(1501..1560,\n                     1601..1840)\n" + qualifier("gene", "three") + qualifier("transl_table", "11", false) + qualifier("translation", p3);
    gb += "                     /pseudo\n";
    gb += "ORIGIN      \n";
    const std::string letters = decode_dna(genome);
    for (std::size_t i = 0; i < n; i += 60) {
      std::string line = std::to_string(i + 1);
      line = std::string(9 - line.size(), ' ') + line;
      for (std::size_t j = i; j < std::min(n, i + 60); j += 10) {
        std::string block = letters.substr(j, 10);
        for (auto& c : block) c = static_cast<char>(c + 32);
        line += " " + block;
      }
      gb += line + "\n";
    }
    gb += "//\n";

    std::vector<Problem> problems;
    const auto recs = read_genbank(gb, &problems);
    s.require("GenBank: one record, read without problems", recs.size() == 1 && problems.empty(), problems.empty() ? "" : problems[0].message);
    if (recs.size() == 1) {
      const GenBank& b = recs[0];
      s.require("GenBank: the LOCUS line (name, length, molecule, topology, division, date)",
                b.locus == "TEST0001" && b.length == 3000 && b.molecule == "DNA" && b.circular() && b.division == "PHG" && b.date == "24-OCT-2026");
      s.require("GenBank: a DEFINITION of two lines, and continuations of other keywords kept out of it",
                b.definition == "A test genome, written by the checks of the genomics domain, with three genes.", b.definition);
      s.require("GenBank: ACCESSION, VERSION, ORGANISM and its lineage",
                b.accession == "TEST0001" && b.version == "TEST0001.1" && b.organism == "synthetic construct" && b.lineage == std::vector<std::string>{"other sequences", "artificial sequences"});
      s.require("GenBank: four features, keys and lines", b.features.size() == 4 && b.features[1].key == "CDS" && b.features[0].key == "source");
      s.require("GenBank: the ORIGIN is the genome", encode_dna(b.sequence) == genome);
      if (b.features.size() == 4) {
        const Feature& two = b.features[2];
        s.require("GenBank: a quoted value over lines, its doubled quotation marks single",
                  two.qualifier("note") == "a note that says \"quoted\" words and runs over more than one line of the table", two.qualifier("note"));
        s.require("GenBank: a location continued over lines, and a bare qualifier", b.features[3].location.parts.size() == 2 && b.features[3].has("pseudo") && b.features[3].qualifier("pseudo").empty());
        bool same = true;
        for (std::size_t f = 1; f < 4; f++) same = same && translate_cds(genome, b.features[f]) == b.features[f].qualifier("translation");
        s.require("GenBank: every CDS translated as its /translation (across the origin, on the complement strand, from GTG)", same);
      }
    }
  }

  // ---- the genomes of NCBI, when the data library has fetched them ------------------------------------
  {
    const char* env = std::getenv("SCIENCE_CACHE");
    const std::filesystem::path cache = env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
    const std::filesystem::path dir = cache / "ncbi-genomes";
    std::vector<std::filesystem::path> files;
    if (std::filesystem::is_directory(dir))
      for (const auto& e : std::filesystem::recursive_directory_iterator(dir))
        if (e.path().extension() == ".gb") files.push_back(e.path());
    if (files.empty()) s.skip("the NCBI genomes are not in the data cache (data/cache/ncbi-genomes): their CDSs are not translated");
    for (const auto& file : files) {
      std::vector<Problem> problems;
      const auto recs = read_genbank(read_file(file), &problems);
      std::size_t cds = 0, agree = 0;
      std::string first_wrong;
      for (const auto& rec : recs) {
        const Codes genome = encode_dna(rec.sequence);
        for (const auto& f : rec.features) {
          if (f.key != "CDS" || !f.has("translation") || f.location.remote || f.has("pseudo") || f.has("ribosomal_slippage") || f.has("exception")) continue;
          cds++;
          const std::string t = translate_cds(genome, f);
          if (t == f.qualifier("translation")) agree++;
          else if (first_wrong.empty()) first_wrong = f.qualifier("locus_tag") + " " + f.location.text;
        }
      }
      s.require("NCBI " + file.stem().string() + ": read without problems, and every CDS translated as GenBank's /translation",
                problems.empty() && cds > 0 && agree == cds,
                std::to_string(agree) + " of " + std::to_string(cds) + (first_wrong.empty() ? "" : "; first different: " + first_wrong));
    }
  }

  return s.done();
}
