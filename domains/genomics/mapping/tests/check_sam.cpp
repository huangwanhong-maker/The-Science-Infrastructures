// Checks of SAM: the example of the SAM specification (§1.1) read field by field, its CIGARs
// consistent with SEQ, its NM recomputed from the reference; mappings written and read back
// to the same mappings (FLAG, POS, MAPQ, CIGAR with soft clips, SEQ and QUAL on the strand
// mapped, NM, AS, tp); unmapped reads; and lines that are not records reported as problems
// with their line numbers while the others are read.

#include <genomics/align/types.hpp>
#include <genomics/mapping/mapper.hpp>
#include <genomics/mapping/sam.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <string>
#include <vector>

using namespace genomics;
using mapping::Mapping;
using sequence::Codes;

int main() {
  calx::verification::Suite s;

  // ---- the specification's example -----------------------------------------------------------------------
  {
    // The SAM/BAM Format Specification (version 1.6), §1.1, "An example": the reference and six
    // records, columns separated by tabs as the format requires (the specification displays them
    // aligned with spaces). The unpadded reference is the padded one without its two '*'.
    const std::string ref = "AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT";
    const std::string text =
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:ref\tLN:45\n"
        "r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t*\n"
        "r002\t0\tref\t9\t30\t3S6M1P1I4M\t*\t0\t0\tAAAAGATAAGGATA\t*\n"
        "r003\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCCTAAGCTAA\t*\tSA:Z:ref,29,-,6H5M,17,0;\n"
        "r004\t0\tref\t16\t30\t6M14N5M\t*\t0\t0\tATAGCTTCAGC\t*\n"
        "r003\t2064\tref\t29\t17\t6H5M\t*\t0\t0\tTAGGC\t*\tSA:Z:ref,9,+,5S6M,30,1;\n"
        "r001\t147\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGGCAT\t*\tNM:i:1\n";
    std::vector<sequence::Problem> problems;
    const auto sam = mapping::read_sam(text, &problems);
    s.require("the specification's example reads without problems: 2 header lines, 1 reference of 45 bases, 6 records",
              problems.empty() && sam.header.size() == 2 && sam.references.size() == 1 && sam.references[0].second == 45 && sam.records.size() == 6);
    bool positions = sam.records.size() == 6;
    const std::int64_t pos[6] = {6, 8, 8, 15, 28, 36};
    for (std::size_t i = 0; positions && i < 6; i++) positions = sam.records[i].position == pos[i];
    s.require("POS is read 1-based into 0-based positions (7, 9, 9, 16, 29, 37)", positions);
    const char* written[6] = {"8M2I4M1D3M", "3S6M1P1I4M", "5S6M", "6M14N5M", "6H5M", "9M"};
    bool cigars = sam.records.size() == 6;
    for (std::size_t i = 0; cigars && i < 6; i++)
      cigars = align::query_length(sam.records[i].cigar) == sam.records[i].seq.size() && align::cigar_string(sam.records[i].cigar) == written[i];
    s.require("every CIGAR reads back as written, its query length SEQ's (P and N consume no query, H none)", cigars);
    s.require("FLAG 99 is paired and first, 147 reverse, 2064 supplementary and reverse",
              sam.records.size() == 6 && !sam.records[0].reverse() && sam.records[5].reverse() && sam.records[4].supplementary() && sam.records[4].reverse() &&
                  sam.records[5].pnext == 6 && sam.records[5].tlen == -39);
    // NM of r001 from the reference: 2 inserted, 1 deleted; of its mate, one mismatch (NM:i:1 in the example).
    const Codes reference = sequence::encode_dna(ref);
    const auto nm = [&](const mapping::SamRecord& r) {
      align::Alignment a;
      a.cigar = r.cigar;
      a.query_end = align::query_length(r.cigar);
      a.target_begin = static_cast<std::size_t>(r.position);
      a.target_end = a.target_begin + align::target_length(r.cigar);
      return align::edit_distance(a, sequence::encode_dna(r.seq), reference);
    };
    s.require("NM recomputed from the reference: 3 for r001 (8M2I4M1D3M), 1 for its mate, as its NM:i:1 says",
              sam.records.size() == 6 && nm(sam.records[0]) == 3 && nm(sam.records[5]) == 1 && sam.records[5].tag_int("NM") == 1);
    s.require("optional fields keep their type and value (SA:Z of r003)", sam.records.size() == 6 && sam.records[2].tag("SA") && sam.records[2].tag("SA")->type == 'Z' &&
                                                                                 sam.records[2].tag("SA")->value == "ref,29,-,6H5M,17,0;");
  }

  // ---- mappings written and read back ----------------------------------------------------------------------
  {
    calx::scheduler::ThreadPool pool(4);
    simulate::GenomeSpec spec;
    spec.length = 120000;
    spec.repeats = {3, 3, 800, 0.0};
    const std::vector<Codes> contigs{simulate::make_genome(spec, 5).sequence, simulate::random_sequence(30000, 0.4, 6, 1)};
    const std::vector<std::string> names{"chrA", "chrB"};
    const mapping::Mapper mapper(pool, contigs, mapping::MapperOptions::short_reads());
    auto reads = simulate::short_reads(simulate::Source::haploid(contigs), 600, simulate::ShortReadModel{}, 9);
    // Two reads that map nowhere.
    reads[0].bases = simulate::random_sequence(150, 0.5, 99, 7);
    reads[1].bases = Codes(20, sequence::A);
    reads[0].quality = reads[1].quality = "";
    const auto all = mapper.map_all(pool, reads);
    std::vector<mapping::SamRecord> records;
    std::vector<std::pair<std::size_t, std::size_t>> which;  // (read, mapping)
    for (std::size_t i = 0; i < reads.size(); i++) {
      if (all[i].empty()) records.push_back(mapping::sam_unmapped(reads[i].name, reads[i].bases, reads[i].quality));
      for (std::size_t j = 0; j < all[i].size(); j++) {
        records.push_back(mapping::sam_record(all[i][j], reads[i].name, reads[i].bases, reads[i].quality, names));
        which.push_back({i, j});
      }
    }
    const std::string text = mapping::write_sam(mapping::sam_header(names, {contigs[0].size(), contigs[1].size()}, "genomics-mapping", "check_sam"), records);
    std::vector<sequence::Problem> problems;
    const auto sam = mapping::read_sam(text, &problems);
    s.require("written SAM reads back without problems, record for record", problems.empty() && sam.records.size() == records.size() && sam.references.size() == 2,
              std::to_string(problems.size()) + " problems");
    bool same = true, strand = true;
    std::size_t k = 0, unmapped = 0, secondary = 0;
    for (const auto& r : sam.records) {
      if (r.unmapped()) { unmapped++; same = same && r.position == -1 && r.cigar.empty() && r.rname == "*"; continue; }
      const auto [i, j] = which[k++];
      const Mapping& m = all[i][j];
      Mapping back;
      if (!mapping::mapping_of(r, names, back)) { same = false; continue; }
      secondary += !m.primary;
      same = same && back.contig == m.contig && back.target_begin == m.target_begin && back.target_end == m.target_end && back.reverse == m.reverse &&
             back.mapq == m.mapq && back.primary == m.primary && back.supplementary == m.supplementary && back.alignment.cigar == m.alignment.cigar &&
             back.alignment.query_begin == m.alignment.query_begin && back.alignment.query_end == m.alignment.query_end && back.alignment.score == m.alignment.score &&
             back.nm == m.nm && back.query_length == m.query_length && r.tag("tp") && r.tag("tp")->value == (m.primary ? "P" : "S");
      const Codes onstrand = m.reverse ? sequence::reverse_complement(reads[i].bases) : reads[i].bases;
      std::string qual = reads[i].quality;
      if (m.reverse) std::reverse(qual.begin(), qual.end());
      strand = strand && r.seq == sequence::decode_dna(onstrand) && (qual.empty() ? r.qual == "*" : r.qual == qual);
    }
    s.require("every mapping reads back the same: contig, POS, strand, MAPQ, flags, CIGAR with clips, AS, NM, tp", same && k == which.size(),
              std::to_string(k) + " mappings, " + std::to_string(secondary) + " secondary");
    s.require("SEQ and QUAL are written on the strand mapped", strand);
    s.require("unmapped reads are written with FLAG 4, RNAME '*', POS 0 and CIGAR '*'", unmapped == 2, std::to_string(unmapped));
  }

  // ---- problems ----------------------------------------------------------------------------------------
  {
    const std::string text =
        "@HD\tVN:1.6\n"
        "@SQ\tSN:x\n"
        "q1\t0\tx\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n"
        "q2\tzero\tx\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n"
        "q3\t0\tx\t5\t60\t10Q\t*\t0\t0\tACGTACGTAC\t*\n"
        "q4\t0\tx\t5\n"
        "q5\t0\tx\t5\t60\t8M\t*\t0\t0\tACGTACGTAC\t*\tNM\n"
        "q6\t16\tx\t9\t3\t4S6M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tAS:i:12\n";
    std::vector<sequence::Problem> problems;
    const auto sam = mapping::read_sam(text, &problems);
    std::vector<std::size_t> lines;
    for (const auto& p : problems) lines.push_back(p.line);
    s.require("lines that are not records are problems with their line numbers (2: @SQ without LN; 4: FLAG; 5: CIGAR; 6: fields; 7: CIGAR length and a malformed tag)",
              lines == std::vector<std::size_t>{2, 4, 5, 6, 7, 7}, std::to_string(problems.size()) + " problems");
    s.require("the other records are read, never throwing", sam.records.size() == 3 && sam.records[2].tag_int("AS") == 12 && sam.records[2].reverse() && sam.records[2].position == 8);
  }

  return s.done();
}
