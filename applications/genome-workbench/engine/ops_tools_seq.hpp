// The toolbox's operations on sequences and their annotation (src/tools/sequence.js and
// annotation.js): composition, the other strand, translation in one, three or six frames, open
// reading frames, motifs in IUPAC letters on both strands, the k-mer spectrum, codon usage; and
// genes, CpG islands and repeats found from a sequence alone.
//
// The sequences are the user's, given as text or as a file (inputs.hpp). Every computation is the
// domain's (genomics::sequence, genomics::index, genomics::kmers, genomics::annotation,
// genomics::hmm): the engine reads the inputs, holds them within bounds and reports.

#pragma once

#include "genomes.hpp"
#include "inputs.hpp"
#include "shared.hpp"
#include "shared_annotation.hpp"

#include <genomics/assembly/assemble.hpp>
#include <genomics/index/patterns.hpp>
#include <genomics/kmers/counter.hpp>
#include <genomics/kmers/spectrum.hpp>
#include <genomics/sequence/composition.hpp>
#include <genomics/sequence/genetic_code.hpp>
#include <genomics/sequence/kmer.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <optional>
#include <queue>
#include <string>
#include <utility>
#include <vector>

namespace workbench {

namespace detail {

// The protein of an ORF or a gene, [begin, end) of the forward strand read on its own strand:
// its first codon M when it is a start, its final stop dropped (translate's coding form).
inline std::string protein_of(const seq::Codes& x, std::size_t begin, std::size_t end, bool reverse, const seq::GeneticCode& code) {
  seq::Codes part(x.begin() + static_cast<std::ptrdiff_t>(begin), x.begin() + static_cast<std::ptrdiff_t>(end));
  if (reverse) part = seq::reverse_complement(part);
  return seq::translate(part, code, 0, true);
}

// The bases of a first record of DNA, at most `limit` of them.
inline seq::Codes dna_of(const seq::Record& r, std::size_t limit, const char* what) {
  require_dna(r, what);
  if (r.sequence.size() > limit) throw Fault(r.id + " has " + std::to_string(r.sequence.size()) + " bases; " + what + " here takes up to " + std::to_string(limit));
  return seq::encode_dna(r.sequence);
}

}  // namespace detail

// Each sequence's length, GC, GC skew, CpG observed over expected, entropy of its 2-mers and
// ambiguous letters; and the whole input's (N50 by the assembly module's statistics). A file is
// read a record at a time, so a genome of gigabytes is measured without being held.
BENCH_OP("tool.seq.stats") {
  detail::Bounds b;
  b.bytes = std::size_t{4} << 30;
  b.records = 5'000'000;
  b.whole = false;
  Json rows = Json::array();
  std::vector<std::size_t> lengths;
  std::array<std::size_t, 5> all{};
  std::size_t protein = 0, ambiguous = 0;
  const auto reading = each_record(request, "fasta", b, [&](seq::Record&& r) {
    lengths.push_back(r.sequence.size());
    Json row = Json::object({{"name", r.id}, {"length", r.sequence.size()}});
    if (detail::looks_protein(r.sequence)) {
      protein++;
      const auto x = static_cast<std::size_t>(std::count_if(r.sequence.begin(), r.sequence.end(), [](char c) { return c == 'X' || c == 'x'; }));
      ambiguous += x;
      row.set("alphabet", "protein");
      for (const char* k : {"gc", "gc_skew", "cpg_oe", "entropy"}) row.set(k, Json());
      row.set("ambiguous", x);
    } else {
      const seq::Codes codes = seq::encode_dna(r.sequence);
      const auto k = seq::base_counts(codes);
      for (int i = 0; i < 5; i++) all[static_cast<std::size_t>(i)] += k[static_cast<std::size_t>(i)];
      ambiguous += k[4];
      const double acgt = static_cast<double>(k[0] + k[1] + k[2] + k[3]), gc = static_cast<double>(k[1] + k[2]);
      row.set("alphabet", "dna");
      row.set("gc", acgt > 0 ? Json(detail::round_to(gc / acgt, 1e5)) : Json());
      row.set("gc_skew", gc > 0 ? Json(detail::round_to((static_cast<double>(k[2]) - static_cast<double>(k[1])) / gc, 1e4)) : Json());
      row.set("cpg_oe", k[1] && k[2] ? Json(detail::round_to(seq::cpg_observed_expected(codes), 1e3)) : Json());
      row.set("entropy", codes.size() >= 2 ? Json(detail::round_to(seq::kmer_entropy(codes, 2), 1e4)) : Json());
      row.set("ambiguous", k[4]);
    }
    if (rows.size() < 20000) rows.push(std::move(row));
  });
  const auto s = genomics::assembly::statistics(lengths);
  const double acgt = static_cast<double>(all[0] + all[1] + all[2] + all[3]);
  Json total = Json::object({{"records", lengths.size()}, {"length", s.total}, {"gc", acgt > 0 ? Json(detail::round_to(static_cast<double>(all[1] + all[2]) / acgt, 1e5)) : Json()},
                             {"n50", s.n50}, {"l50", s.l50}, {"longest", s.largest}, {"ambiguous", ambiguous}, {"protein", protein},
                             {"bases", Json::object({{"A", all[0]}, {"C", all[1]}, {"G", all[2]}, {"T", all[3]}, {"N", all[4]}})}});
  return Json::object({{"records", rows}, {"listed", rows.size()}, {"total", total}, {"reading", reading.json()}});
}

// The other strand read 5′ to 3′ (IUPAC letters complemented, case kept), the complement alone,
// or the letters reversed.
BENCH_OP("tool.seq.revcomp") {
  const std::string what = str(request, "what", "revcomp");
  if (what != "revcomp" && what != "complement" && what != "reverse") throw Fault("give revcomp, complement or reverse");
  detail::Bounds b;
  b.bytes = std::size_t{32} << 20;
  Json records = Json::array();
  const auto reading = each_record(request, "fasta", b, [&](seq::Record&& r) {
    std::string out;
    if (what == "reverse") {
      out.assign(r.sequence.rbegin(), r.sequence.rend());
    } else {
      require_dna(r, "the complement");
      out = seq::reverse_complement(std::string_view(r.sequence));
      if (what == "complement") std::reverse(out.begin(), out.end());
    }
    records.push(Json::object({{"name", r.id}, {"seq", std::move(out)}}));
  });
  return Json::object({{"records", records}, {"what", what}, {"reading", reading.json()}});
}

// Translations under a genetic code of NCBI's: frame 1, the three forward frames, or all six
// (frames −1 to −3 read the reverse complement from its start). Every codon is translated, stops
// as *; a codon with an ambiguous base gives the amino acid every reading of it gives, else X.
BENCH_OP("tool.seq.translate") {
  const auto& code = genetic_code_of(request, 1);
  const auto frames = bounded(request, "frames", 1, 1, 6);
  if (frames != 1 && frames != 3 && frames != 6) throw Fault("translate one, three or six frames");
  detail::Bounds b;
  b.bytes = std::size_t{16} << 20;
  Json records = Json::array();
  const auto reading = each_record(request, "fasta", b, [&](seq::Record&& r) {
    require_dna(r, "translation");
    const seq::Codes x = seq::encode_dna(r.sequence);
    for (int f = 0; f < std::min<int>(3, static_cast<int>(frames)); f++)
      records.push(Json::object({{"name", r.id}, {"frame", f + 1}, {"protein", seq::translate(x, code, static_cast<std::size_t>(f))}}));
    if (frames == 6) {
      const seq::Codes y = seq::reverse_complement(x);
      for (int f = 0; f < 3; f++) records.push(Json::object({{"name", r.id}, {"frame", -(f + 1)}, {"protein", seq::translate(y, code, static_cast<std::size_t>(f))}}));
    }
  });
  return Json::object({{"records", records}, {"table", code.id}, {"code", std::string(code.name)}, {"reading", reading.json()}});
}

// The ORFs of both strands of the first sequence: from a start (ATG, or every start of the code)
// to the first stop in frame, at least so many codons (the stop included), each with its protein.
BENCH_OP("tool.seq.orfs") {
  const auto& code = genetic_code_of(request, 11);
  const auto min_codons = static_cast<std::size_t>(bounded(request, "min_codons", 60, 10, 100000));
  const bool any_start = flag(request, "any_start", false);
  detail::Bounds b;
  b.bytes = std::size_t{64} << 20;
  b.whole = false;
  detail::Reading reading;
  const auto first = first_record(request, "fasta", b, reading);
  const seq::Codes x = detail::dna_of(first.record, 50'000'000, "the search of ORFs");
  const auto orfs = seq::open_reading_frames(x, code, min_codons, any_start);
  Json list = Json::array();
  for (const auto& o : orfs) {
    if (list.size() >= 20000) break;
    list.push(Json::object({{"begin", o.begin}, {"end", o.end}, {"strand", o.reverse ? -1 : 1}, {"frame", o.frame + 1}, {"codons", o.codons()},
                            {"protein", detail::protein_of(x, o.begin, o.end, o.reverse, code)}}));
  }
  return Json::object({{"orfs", list}, {"count", orfs.size()}, {"length", x.size()}, {"record", first.record.id}, {"table", code.id},
                       {"starts", any_start ? "every start codon of the code" : "ATG"}, {"reading", reading.json()}});
}

// Every occurrence of motifs in IUPAC letters on both strands, by one Aho–Corasick automaton of
// their words (genomics::index::PatternSet); a palindromic site once.
BENCH_OP("tool.seq.find") {
  std::vector<std::string> patterns;
  const Json& p = request["patterns"];
  if (p.is_array())
    for (const auto& x : p.as_array())
      if (!x.text().empty()) patterns.push_back(x.text());
  if (patterns.empty()) throw Fault("give motifs to find");
  if (patterns.size() > 200) throw Fault("at most 200 motifs at once");
  for (const auto& x : patterns)
    if (x.size() > 100) throw Fault("a motif has at most 100 letters");
  const bool both = flag(request, "both", true);
  std::optional<genomics::index::PatternSet> set;
  try {
    set.emplace(patterns, both, std::size_t{1} << 20);
  } catch (const std::invalid_argument& e) {
    throw Fault(detail::message_of(e));
  }
  detail::Bounds b;
  b.bytes = std::size_t{1} << 30;
  b.whole = false;
  const std::size_t cap = static_cast<std::size_t>(bounded(request, "max_hits", 5000, 1, 100000));
  std::vector<std::array<std::size_t, 2>> counts(patterns.size(), {0, 0});
  std::size_t total = 0, searched = 0;
  Json hits = Json::array();
  const auto reading = each_record(request, "fasta", b, [&](seq::Record&& r) {
    require_dna(r, "the search of motifs");
    const seq::Codes x = seq::encode_dna(r.sequence);
    searched += x.size();
    for (const auto& h : set->find(x)) {
      total++;
      counts[h.pattern][h.reverse ? 1 : 0]++;
      if (hits.size() >= cap) continue;
      seq::Codes w(x.begin() + static_cast<std::ptrdiff_t>(h.begin), x.begin() + static_cast<std::ptrdiff_t>(h.end));
      if (h.reverse) w = seq::reverse_complement(w);
      hits.push(Json::object({{"record", r.id}, {"motif", patterns[h.pattern]}, {"begin", h.begin}, {"end", h.end}, {"strand", h.reverse ? -1 : 1}, {"match", seq::decode_dna(w)}}));
    }
  });
  Json by = Json::array();
  std::size_t words = 0;
  for (std::size_t k = 0; k < patterns.size(); k++) {
    words += set->forward_words(k);
    by.push(Json::object({{"motif", patterns[k]}, {"words", set->forward_words(k)}, {"count", counts[k][0] + counts[k][1]}, {"forward", counts[k][0]}, {"reverse", counts[k][1]}}));
  }
  return Json::object({{"count", total}, {"patterns", words}, {"automaton_words", set->words()}, {"searched", searched}, {"both", both}, {"by_motif", by}, {"hits", hits},
                       {"hits_listed", hits.size()}, {"reading", reading.json()}});
}

// Canonical k-mers counted (genomics::kmers, partitioned and deterministic) over sequences or reads,
// a batch of records at a time: the spectrum of their counts and the commonest. Reading stops at a
// bound of the file or of distinct k-mers, so that the tables fit in a worker's memory.
BENCH_OP("tool.seq.kmers") {
  const int k = static_cast<int>(bounded(request, "k", 21, 3, 31));
  const auto top = static_cast<std::size_t>(bounded(request, "top", 25, 1, 1000));
  constexpr std::size_t DISTINCT = 80'000'000, BATCH = std::size_t{16} << 20;
  genomics::kmers::CounterOptions co;
  co.k = k;
  genomics::kmers::KmerCounter counter(co);
  std::vector<seq::Codes> batch;
  std::size_t in_batch = 0;
  bool full = false;
  const auto flush = [&] {
    if (batch.empty()) return;
    counter.add(batch, worker.pool);
    batch.clear();
    in_batch = 0;
    full = counter.distinct() > DISTINCT;
  };
  detail::Bounds b;
  b.bytes = std::size_t{1} << 30;
  b.records = 50'000'000;
  b.whole = false;
  const std::string name = !str(request, "reads_file").empty() ? "reads" : "fasta";
  auto reading = each_record(request, name, b, [&](seq::Record&& r) {
    require_dna(r, "counting k-mers");
    in_batch += r.sequence.size();
    batch.push_back(seq::encode_dna(r.sequence));
    if (in_batch >= BATCH) flush();
    return !full;
  });
  flush();
  if (full) reading.stopped = "reading stopped at " + std::to_string(counter.distinct()) + " distinct k-mers, the most a worker holds here";
  counter.release_buffers();
  const auto spec = genomics::kmers::spectrum(counter, worker.pool, 10000);
  Json spectrum = Json::array();
  for (std::size_t c = 1; c < spec.histogram.size() && spectrum.size() < 2000; c++)
    if (spec.histogram[c]) spectrum.push(Json::array({c, spec.histogram[c]}));
  // The commonest: a heap of the `top` largest counts, then sorted, larger first and ties by k-mer.
  using Entry = std::pair<std::uint32_t, seq::Kmer>;
  std::priority_queue<Entry, std::vector<Entry>, std::greater<>> heap;
  counter.for_each([&](seq::Kmer x, std::uint32_t c) {
    if (heap.size() < top) heap.push({c, x});
    else if (Entry{c, x} > heap.top()) { heap.pop(); heap.push({c, x}); }
  });
  std::vector<Entry> best;
  while (!heap.empty()) { best.push_back(heap.top()); heap.pop(); }
  std::sort(best.begin(), best.end(), [](const Entry& a, const Entry& b) { return a.first != b.first ? a.first > b.first : a.second < b.second; });
  Json common = Json::array();
  for (const auto& [c, x] : best) common.push(Json::array({seq::kmer_string(x, k), c}));
  return Json::object({{"k", k}, {"distinct", counter.distinct()}, {"total", counter.total()}, {"spectrum", spectrum}, {"top", common},
                       {"memory", counter.peak_bytes()}, {"reading", reading.json()}});
}

// Codon usage over coding sequences read from their first base in whole codons (codons with an
// ambiguous base skipped): each codon's count and relative synonymous usage (Sharp, Tuohy and
// Mosurski 1986), and Wright's effective number of codons (1990), under a genetic code.
BENCH_OP("tool.seq.codons") {
  const auto& code = genetic_code_of(request, 1);
  detail::Bounds b;
  b.bytes = std::size_t{512} << 20;
  b.whole = false;
  std::array<std::size_t, 64> counts{};
  std::size_t partial = 0;
  const auto reading = each_record(request, "fasta", b, [&](seq::Record&& r) {
    require_dna(r, "codon usage");
    partial += r.sequence.size() % 3 != 0;
    const auto n = seq::codon_counts({seq::encode_dna(r.sequence)});
    for (std::size_t i = 0; i < 64; i++) counts[i] += n[i];
  });
  const auto rscu = seq::rscu(counts, code);
  std::size_t counted = 0;
  Json codons = Json::array();
  for (int i = 0; i < 64; i++) {
    const auto c = seq::codon_of(i);
    const auto u = static_cast<std::size_t>(i);
    counted += counts[u];
    codons.push(Json::object({{"codon", seq::decode_dna({c[0], c[1], c[2]})}, {"aa", std::string(1, code.amino_acids[u])}, {"count", counts[u]}, {"rscu", detail::round_to(rscu[u], 1e4)}}));
  }
  return Json::object({{"codons", codons}, {"codons_counted", counted}, {"enc", counted ? Json(detail::round_to(seq::effective_number_of_codons(counts, code), 100)) : Json()},
                       {"table", code.id}, {"code", std::string(code.name)}, {"partial", partial}, {"reading", reading.json()}});
}

// Genes of the first sequence found from it alone: periodic Markov chains of coding trained on the
// sequence itself, round after round (genomics::annotation::find_genes).
BENCH_OP("tool.annot.genes") {
  const auto& code = genetic_code_of(request, 11);
  detail::Bounds b;
  b.bytes = std::size_t{32} << 20;
  b.whole = false;
  detail::Reading reading;
  const auto first = first_record(request, "fasta", b, reading);
  const seq::Codes x = detail::dna_of(first.record, 20'000'000, "the gene finder");
  genomics::annotation::GeneReport report;
  const auto genes = detail::predicted_genes(x, code.id, &report);
  static constexpr const char* STARTS[3] = {"ATG", "GTG", "TTG"};
  Json list = Json::array();
  std::vector<std::pair<std::size_t, std::size_t>> spans;
  for (const auto& g : genes) {
    spans.push_back({g.begin, g.end});
    if (list.size() >= 20000) continue;
    list.push(Json::object({{"begin", g.begin}, {"end", g.end}, {"strand", g.reverse ? -1 : 1}, {"codons", (g.end - g.begin) / 3}, {"score", detail::round_to(g.score, 10)},
                            {"start", g.start >= 0 && g.start < 3 ? STARTS[g.start] : "?"}, {"protein", detail::protein_of(x, g.begin, g.end, g.reverse, code)}}));
  }
  // The share of the sequence some gene covers: the union of their intervals.
  std::sort(spans.begin(), spans.end());
  std::size_t covered = 0, reach = 0;
  for (const auto& [lo, hi] : spans) {
    if (hi <= reach) continue;
    covered += hi - std::max(lo, reach);
    reach = hi;
  }
  return Json::object({{"genes", list}, {"count", genes.size()}, {"length", x.size()}, {"coding", x.empty() ? 0.0 : detail::round_to(static_cast<double>(covered) / static_cast<double>(x.size()), 1e4)},
                       {"rounds", report.rounds}, {"converged", report.converged}, {"candidates", report.candidates}, {"record", first.record.id}, {"table", code.id},
                       {"reading", reading.json()}});
}

// CpG islands of the first sequence by the two-state HMM of Durbin et al. (1998, ch. 3), Viterbi
// decoded, at least 200 bases.
BENCH_OP("tool.annot.cpg") {
  detail::Bounds b;
  b.bytes = std::size_t{64} << 20;
  b.whole = false;
  detail::Reading reading;
  const auto first = first_record(request, "fasta", b, reading);
  const seq::Codes x = detail::dna_of(first.record, 60'000'000, "the search of CpG islands");
  Json list = Json::array();
  const auto islands = detail::cpg_islands_of(x);
  for (const auto& i : islands) {
    if (list.size() >= 20000) break;
    list.push(Json::object({{"begin", i.begin}, {"end", i.end}, {"gc", detail::round_to(i.gc, 1e4)}, {"oe", detail::round_to(i.cpg_observed_expected, 1e3)},
                            {"log_odds", detail::round_to(i.log_odds, 10)}}));
  }
  return Json::object({{"islands", list}, {"count", islands.size()}, {"length", x.size()}, {"record", first.record.id}, {"reading", reading.json()}});
}

// Maximal repeats of the first sequence on both strands, at least so long, by the suffix array of
// the sequence and its reverse complement and their LCP (Gusfield 1997, §7.12); the longest first.
BENCH_OP("tool.annot.repeats") {
  const auto min_length = static_cast<std::size_t>(bounded(request, "min_length", 30, 10, 100000));
  detail::Bounds b;
  b.bytes = std::size_t{32} << 20;
  b.whole = false;
  detail::Reading reading;
  const auto first = first_record(request, "fasta", b, reading);
  const seq::Codes x = detail::dna_of(first.record, 20'000'000, "the search of repeats");
  auto repeats = detail::repeats_of(x, min_length);
  std::stable_sort(repeats.begin(), repeats.end(), [](const auto& a, const auto& c) { return a.length > c.length; });
  Json list = Json::array();
  std::size_t longest = 0;
  for (const auto& r : repeats) {
    longest = std::max(longest, r.length);
    if (list.size() >= 2000) continue;
    std::string where;
    for (std::size_t i = 0; i < r.occurrences.size() && i < 20; i++) {
      const auto& o = r.occurrences[i];
      where += (i ? ", " : "") + std::string(o.reverse ? "−" : "") + std::to_string(o.position + 1) + ".." + std::to_string(o.position + r.length);
    }
    if (r.occurrences.size() > 20) where += ", and " + std::to_string(r.occurrences.size() - 20) + " more";
    const auto& o = r.occurrences.front();
    seq::Codes w(x.begin() + static_cast<std::ptrdiff_t>(o.position), x.begin() + static_cast<std::ptrdiff_t>(o.position + std::min<std::size_t>(r.length, 2000)));
    if (o.reverse) w = seq::reverse_complement(w);
    list.push(Json::object({{"length", r.length}, {"copies", r.occurrences.size()}, {"supermaximal", r.supermaximal}, {"where", where}, {"word", seq::decode_dna(w)}}));
  }
  return Json::object({{"repeats", list}, {"count", repeats.size()}, {"longest", longest}, {"min_length", min_length}, {"length", x.size()}, {"record", first.record.id},
                       {"reading", reading.json()}});
}

}  // namespace workbench
