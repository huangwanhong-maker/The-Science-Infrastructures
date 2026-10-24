// The toolbox's operations of alignment (src/tools/alignment.js and the statistics of data.js):
// many sequences aligned together, progressively, by the centre star or by partial orders; one
// query searched against many by Smith–Waterman, one target to each vector lane, ranked by
// Karlin–Altschul E-values; and what a local score means under a scoring.
//
// Every computation is the domain's (genomics::msa, genomics::align): the engine reads the inputs
// (inputs.hpp), holds them within bounds and reports.

#pragma once

#include "genomes.hpp"
#include "inputs.hpp"

#include <genomics/align/inter.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/align/statistics.hpp>
#include <genomics/msa/centre_star.hpp>
#include <genomics/msa/msa.hpp>
#include <genomics/msa/partial_order.hpp>
#include <genomics/msa/progressive.hpp>
#include <genomics/msa/scores.hpp>

#include <calx/cpu/target.hpp>
#include <calx/scheduler/parallel_for.hpp>

#include <algorithm>
#include <chrono>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace workbench {

namespace detail {

// The alphabet a request names, or the one its letters read as.
inline bool protein_of_request(const Json& request, const std::vector<seq::Record>& records) {
  const std::string a = str(request, "alphabet", "auto");
  if (a == "protein") return true;
  if (a == "dna") return false;
  if (a != "auto") throw Fault("the letters are dna or protein");
  std::string sample;
  for (const auto& r : records) {
    sample += r.sequence.substr(0, 5000 - std::min<std::size_t>(5000, sample.size()));
    if (sample.size() >= 5000) break;
  }
  return looks_protein(sample);
}

// The statistics of local scores under a scoring: NCBI BLAST's gapped λ, K and H where its table
// has the matrix with these gaps (Altschul and Gish 1996), else the ungapped ones computed from the
// matrix under the background (Robinson and Robinson's amino acids, or uniform bases).
struct ScoreStatistics {
  std::optional<genomics::align::KarlinAltschul> ka;
  bool gapped = false;
  std::optional<genomics::align::GappedParameters> row;
  std::string note;
};
inline ScoreStatistics score_statistics(const genomics::align::Scoring& s, bool protein, const std::string& matrix) {
  namespace aln = genomics::align;
  ScoreStatistics out;
  if (protein) {
    if (const auto g = aln::gapped_parameters(matrix, s.gap_open, s.gap_extend)) {
      out.ka = g->statistics();
      out.gapped = true;
      out.row = g;
      return out;
    }
    out.note = "NCBI BLAST tabulates no gapped values for " + matrix + " with gaps " + std::to_string(s.gap_open) + " + L·" + std::to_string(s.gap_extend) +
               ": these are the ungapped values, which make a gapped score look more significant than it is";
  } else {
    out.note = "ungapped values for the scores of DNA under uniform bases; NCBI's gapped values for DNA are not in the domain";
  }
  try {
    const auto bg_p = aln::robinson_frequencies();
    const auto bg_d = aln::uniform_dna_frequencies();
    out.ka = protein ? aln::karlin_altschul(s.matrix, bg_p, bg_p) : aln::karlin_altschul(s.matrix, bg_d, bg_d);
  } catch (const std::invalid_argument& e) {
    out.note = std::string("no Karlin–Altschul statistics: ") + message_of(e);
  }
  return out;
}

}  // namespace detail

// Many sequences aligned together: progressively along a guide tree (MUSCLE's three stages: a
// k-mer tree, a tree re-estimated from the first alignment, refinement), by the centre star
// (Gusfield 1993), or by partial orders (Lee, Grasso and Sharlow 2002); with the sum-of-pairs
// score, the consensus and Clustal's marks of conservation.
BENCH_OP("tool.align.msa") {
  namespace msa = genomics::msa;
  const std::string method = str(request, "method", "progressive");
  if (method != "progressive" && method != "centre_star" && method != "poa") throw Fault("align progressively, by the centre star, or by partial orders (poa)");
  detail::Bounds b;
  b.bytes = std::size_t{8} << 20;
  b.records = 500;
  detail::Reading reading;
  const auto records = all_records(request, "fasta", b, reading);
  if (records.size() < 2) throw Fault("give at least two sequences");
  std::size_t total = 0, longest = 0;
  for (const auto& r : records) {
    if (r.sequence.empty()) throw Fault(r.id + " is empty");
    total += r.sequence.size();
    longest = std::max(longest, r.sequence.size());
  }
  if (total > 400'000 || longest > 30'000) throw Fault("a multiple alignment here holds up to 400,000 residues, none longer than 30,000");
  const bool protein = detail::protein_of_request(request, records);
  const auto A = protein ? msa::Alphabet::protein : msa::Alphabet::dna;
  const auto names = distinct_names(records);
  std::vector<seq::Codes> seqs;
  for (const auto& r : records) seqs.push_back(protein ? seq::encode_protein(r.sequence) : seq::encode_dna(r.sequence));
  const auto scoring = msa::default_scoring(A);
  const auto t0 = std::chrono::steady_clock::now();
  msa::Msa m;
  Json tree;
  if (method == "progressive") {
    auto res = msa::progressive(seqs, names, A, msa::default_options(A), worker.pool);
    tree = res.tree.newick(names);
    m = std::move(res.msa);
  } else if (method == "centre_star") {
    m = msa::centre_star(seqs, names, scoring, worker.pool).msa;
  } else {
    m = msa::partial_order_alignment(seqs, names, scoring);
  }
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  const std::string conservation = msa::clustal_conservation(m, A);
  Json rows = Json::array();
  for (std::size_t i = 0; i < m.size(); i++) rows.push(Json::object({{"name", m.names[i]}, {"row", msa::row_string(m.rows[i], A)}}));
  return Json::object({{"rows", rows}, {"columns", m.columns()}, {"sp", msa::sp_score(m, scoring, false, worker.pool)},
                       {"identical_columns", static_cast<std::size_t>(std::count(conservation.begin(), conservation.end(), '*'))},
                       {"alphabet", protein ? "protein" : "dna"}, {"method", method}, {"consensus", msa::row_string(msa::consensus(m, msa::letters(A)), A)},
                       {"conservation", conservation}, {"tree", tree},
                       {"scoring", scoring.matrix.name + ", gaps " + std::to_string(scoring.gap_open) + " + L·" + std::to_string(scoring.gap_extend)},
                       {"seconds", seconds}, {"reading", reading.json()}});
}

// One query against every sequence of a set by Smith–Waterman with affine gaps, one target to each
// vector lane (Rognes 2011), the batches on the worker's threads; the best ranked by score, with
// bit scores and E-values in a search of this size (Karlin and Altschul 1990).
BENCH_OP("tool.align.search") {
  namespace aln = genomics::align;
  const bool protein = str(request, "alphabet", "dna") == "protein";
  const std::string q = str(request, "query");
  if (q.empty()) throw Fault("give a query");
  if (q.size() > 10000) throw Fault("a query has at most 10,000 letters here");
  const seq::Codes query = protein ? seq::encode_protein(q) : seq::encode_dna(q);
  const std::string matrix = protein ? str(request, "matrix", "BLOSUM62") : "dna";
  aln::Scoring s;
  try {
    s = aln::scoring(protein ? aln::matrix(matrix) : aln::Matrix::dna(static_cast<int>(bounded(request, "match", 2, 1, 100)), static_cast<int>(bounded(request, "mismatch", -4, -100, -1)), -1),
                     static_cast<int>(bounded(request, "gap_open", protein ? 11 : 5, 0, 100)), static_cast<int>(bounded(request, "gap_extend", protein ? 1 : 2, 1, 100)));
  } catch (const std::invalid_argument& e) {
    throw Fault(detail::message_of(e));
  }
  const auto top = static_cast<std::size_t>(bounded(request, "top", 50, 1, 1000));
  detail::Bounds b;
  b.bytes = std::size_t{256} << 20;
  b.records = 5'000'000;
  std::vector<std::string> names;
  std::vector<seq::Codes> targets;
  std::size_t residues = 0;
  const auto reading = each_record(request, "fasta", b, [&](seq::Record&& r) {
    names.push_back(std::move(r.id));
    targets.push_back(protein ? seq::encode_protein(r.sequence) : seq::encode_dna(r.sequence));
    residues += targets.back().size();
  });
  if (targets.empty()) throw Fault("give the targets to search");
  const std::size_t cells = query.size() * residues;
  if (cells > 500'000'000'000ull) throw Fault("a search here computes at most 5 × 10¹¹ cells");
  std::vector<std::int64_t> scores(targets.size(), 0);
  const auto t0 = std::chrono::steady_clock::now();
  calx::scheduler::parallel_for(worker.pool, targets.size(), 512, [&](std::size_t lo, std::size_t hi, std::size_t) {
    const std::vector<seq::Codes> part(targets.begin() + static_cast<std::ptrdiff_t>(lo), targets.begin() + static_cast<std::ptrdiff_t>(hi));
    const auto got = aln::inter_scores(query, part, s);
    std::copy(got.begin(), got.end(), scores.begin() + static_cast<std::ptrdiff_t>(lo));
  });
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  const auto st = detail::score_statistics(s, protein, matrix);
  std::vector<std::size_t> order(targets.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  const std::size_t shown = std::min(top, order.size());
  std::partial_sort(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(shown), order.end(), [&](std::size_t a, std::size_t c) { return scores[a] != scores[c] ? scores[a] > scores[c] : a < c; });
  Json hits = Json::array();
  for (std::size_t k = 0; k < shown; k++) {
    const std::size_t i = order[k];
    Json h = Json::object({{"name", names[i]}, {"length", targets[i].size()}, {"score", scores[i]}});
    h.set("bits", st.ka ? Json(aln::bit_score(static_cast<double>(scores[i]), *st.ka)) : Json());
    h.set("evalue", st.ka ? Json(aln::evalue(static_cast<double>(scores[i]), static_cast<double>(query.size()), static_cast<double>(residues), *st.ka)) : Json());
    hits.push(std::move(h));
  }
  Json statistics;
  if (st.ka) statistics = Json::object({{"lambda", st.ka->lambda}, {"K", st.ka->K}, {"H", st.ka->H}, {"gapped", st.gapped}, {"note", st.note}});
  const std::string kernel = "one target to each of " + std::to_string(calx::cpu::LANE_BYTES / 2) + " lanes of 16 bits (" + calx::cpu::ISA + "), " + std::to_string(worker.pool.size()) + " threads";
  return Json::object({{"searched", targets.size()}, {"query_length", query.size()}, {"database_length", residues}, {"cells", cells}, {"seconds", seconds},
                       {"gcups", seconds > 0 ? static_cast<double>(cells) / seconds / 1e9 : 0}, {"kernel", kernel}, {"statistics", statistics},
                       {"note", st.note}, {"hits", hits}, {"reading", reading.json()}});
}

// A local score's bits, E-value and P-value in a search of m by n residues, under a protein matrix
// with affine gaps: λ, K and H from NCBI BLAST's table where it has these gaps, else the ungapped
// values computed from the matrix (Karlin and Altschul 1990).
BENCH_OP("tool.stats.evalue") {
  namespace aln = genomics::align;
  const double score = num(request, "score", 60), m = num(request, "m", 300), n = num(request, "n", 1e8);
  if (!(score >= 0)) throw Fault("a score is at least 0");
  if (!(m >= 1 && n >= 1)) throw Fault("the lengths of the query and of the database are at least 1");
  const std::string matrix = str(request, "matrix", "BLOSUM62");
  aln::Scoring s;
  try {
    s = aln::scoring(aln::matrix(matrix), static_cast<int>(bounded(request, "gap_open", 11, 0, 100)), static_cast<int>(bounded(request, "gap_extend", 1, 1, 100)));
  } catch (const std::invalid_argument& e) {
    throw Fault(detail::message_of(e));
  }
  const auto st = detail::score_statistics(s, true, matrix);
  if (!st.ka) throw Fault(st.note);
  const double e = aln::evalue(score, m, n, *st.ka);
  Json out = Json::object({{"bits", aln::bit_score(score, *st.ka)}, {"evalue", e}, {"pvalue", aln::pvalue(e)}, {"gapped", st.gapped}, {"lambda", st.ka->lambda},
                           {"K", st.ka->K}, {"H", st.ka->H}, {"matrix", matrix}, {"gap_open", s.gap_open}, {"gap_extend", s.gap_extend}, {"m", m}, {"n", n}, {"note", st.note}});
  if (st.row) out.set("blast", Json::object({{"alpha", st.row->alpha}, {"beta", st.row->beta}}));
  return out;
}

}  // namespace workbench
