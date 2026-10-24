// The operations of the Align view: two sequences aligned by one schedule, with what the page
// draws of it (the field of best scores and the path, the band, the wavefronts); every schedule
// on one pair, timed and held to the reference; the same edit graph in other semirings; the
// parametric polygon and its fan; and a pair of random relatives to try them on.
//
// Every schedule is the domain's (genomics::align): the engine only chooses, times and reports.

#pragma once

#include "genomes.hpp"

#include <genomics/align/banded.hpp>
#include <genomics/align/bitparallel.hpp>
#include <genomics/align/gotoh.hpp>
#include <genomics/align/hirschberg.hpp>
#include <genomics/align/inter.hpp>
#include <genomics/align/pairhmm.hpp>
#include <genomics/align/poa.hpp>
#include <genomics/align/statistics.hpp>
#include <genomics/align/parametric.hpp>
#include <genomics/align/reference.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/align/striped.hpp>
#include <genomics/align/wavefront.hpp>
#include <genomics/simulate/family.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/scheduler/parallel_for.hpp>

#include <chrono>
#include <functional>
#include <string>
#include <vector>

namespace workbench {

namespace aln = genomics::align;

namespace detail {

struct Pair {
  seq::Codes a, b;
  bool dna = true;
  aln::Scoring scoring;
  aln::Mode mode = aln::Mode::global;
  seq::Code limit() const { return dna ? 4 : 20; }  // codes at or above it match nothing
};

inline aln::Mode mode_of(const std::string& m) {
  if (m == "local") return aln::Mode::local;
  if (m == "semiglobal") return aln::Mode::semiglobal;
  if (m == "overlap") return aln::Mode::overlap;
  if (m == "global" || m.empty()) return aln::Mode::global;
  throw Fault("a mode is global, local, semiglobal or overlap");
}

inline Pair pair_of(const Json& r, std::size_t limit = 200000) {
  Pair p;
  p.dna = str(r, "alphabet", "dna") != "protein";
  const std::string a = str(r, "a"), b = str(r, "b");
  if (a.empty() || b.empty()) throw Fault("two sequences are needed");
  if (a.size() > limit || b.size() > limit) throw Fault("a sequence has at most " + std::to_string(limit) + " letters here");
  p.a = p.dna ? seq::encode_dna(a) : seq::encode_protein(a);
  p.b = p.dna ? seq::encode_dna(b) : seq::encode_protein(b);
  const std::string m = str(r, "matrix", p.dna ? "dna" : "BLOSUM62");
  aln::Matrix matrix = m == "dna" ? aln::Matrix::dna(static_cast<int>(integer(r, "match", 2)), static_cast<int>(integer(r, "mismatch", -4)), -1) : aln::matrix(m);
  p.scoring = aln::scoring(matrix, static_cast<int>(bounded(r, "gap_open", 4, 0, 100)), static_cast<int>(bounded(r, "gap_extend", 2, 1, 100)));
  p.mode = mode_of(str(r, "mode", "global"));
  return p;
}

// Seconds of one call, from as many calls as fill about 40 ms (at most 200).
inline double timed(const std::function<void()>& f) {
  using clock = std::chrono::steady_clock;
  const auto t0 = clock::now();
  f();
  double one = std::chrono::duration<double>(clock::now() - t0).count();
  if (one > 0.04) return one;
  const int reps = static_cast<int>(std::min(200.0, 0.04 / std::max(one, 1e-7)));
  const auto t1 = clock::now();
  for (int k = 0; k < reps; k++) f();
  return std::chrono::duration<double>(clock::now() - t1).count() / std::max(1, reps);
}

inline Json rows_json(const aln::Alignment& al, const Pair& p) {
  const auto r = aln::render(al, p.a, p.b, p.dna ? seq::dna_letter : seq::protein_letter);
  return Json::object({{"query", r.query}, {"middle", r.middle}, {"target", r.target}});
}

inline Json alignment_json(const aln::Alignment& al, const Pair& p) {
  return Json::object({{"score", al.score}, {"cigar", aln::cigar_string(al.cigar)}, {"query_begin", al.query_begin}, {"query_end", al.query_end},
                       {"target_begin", al.target_begin}, {"target_end", al.target_end}, {"identity", std::round(aln::identity(al, p.a, p.b, p.limit()) * 1e5) / 1e5},
                       {"distance", aln::edit_distance(al, p.a, p.b, p.limit())}, {"rows", rows_json(al, p)}});
}

// The cells of a path, from its start, as (i, j) of the edit graph (i down the query).
inline Json path_json(const aln::Alignment& al) {
  Json path = Json::array();
  std::size_t i = al.query_begin, j = al.target_begin;
  path.push(Json::array({i, j}));
  for (const auto& o : al.cigar)
    for (std::uint32_t k = 0; k < o.length; k++) {
      if (aln::consumes_query(o.op) && o.op != 'S') i++;
      if (aln::consumes_target(o.op)) j++;
      path.push(Json::array({i, j}));
    }
  return path;
}

// The penalties of the pair's scoring (Eizenga and Paten 2022), for the WFA: exact only for one
// match and one mismatch score on A, C, G, T.
inline aln::PenaltyForm penalties_of(const Pair& p) {
  if (!p.dna) throw Fault("the WFA minimises penalties: one match and one mismatch score on DNA");
  if (seq::count_ambiguous(p.a) + seq::count_ambiguous(p.b) > 0) throw Fault("the WFA's penalties equal the scores only on A, C, G and T, and these sequences hold N");
  try { return aln::penalty_form(p.scoring); } catch (const std::invalid_argument& e) { throw Fault(e.what()); }
}

}  // namespace detail

// Two random relatives: a root and a copy diverged by substitutions and indels.
BENCH_OP("align.random") {
  const auto n = static_cast<std::size_t>(bounded(request, "length", 300, 10, 100000));
  const double d = num(request, "divergence", 0.1);
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  const auto root = sim::random_sequence(n, 0.5, seed);
  const auto fam = sim::evolve_family(root, 2, d / 2, d / 10, seed + 1);
  return Json::object({{"a", seq::decode_dna(fam.sequences[0])}, {"b", seq::decode_dna(fam.sequences[1])}});
}

BENCH_OP("align.pair") {
  const auto p = detail::pair_of(request);
  const std::string algorithm = str(request, "algorithm", "gotoh");
  const auto n = p.a.size(), m = p.b.size();
  aln::Alignment al;
  double seconds = 0;
  std::size_t cells = n * m;
  Json extra = Json::object();
  Json fronts_of_wfa;
  if (algorithm == "gotoh") {
    seconds = detail::timed([&] { al = aln::gotoh(p.a, p.b, p.scoring, p.mode); });
  } else if (algorithm == "hirschberg") {
    seconds = detail::timed([&] { al = aln::hirschberg(p.a, p.b, p.scoring, p.mode); });
  } else if (algorithm == "banded") {
    const auto half = bounded(request, "band", 8, 0, 100000);
    const auto d = static_cast<std::int64_t>(m) - static_cast<std::int64_t>(n);
    const std::int64_t lo = std::min<std::int64_t>(0, d) - half, hi = std::max<std::int64_t>(0, d) + half;
    seconds = detail::timed([&] { al = aln::banded(p.a, p.b, p.scoring, p.mode, lo, hi); });
    cells = 0;
    for (std::size_t i = 0; i <= n; i++) {
      const std::int64_t a0 = std::max<std::int64_t>(0, static_cast<std::int64_t>(i) + lo), a1 = std::min<std::int64_t>(static_cast<std::int64_t>(m), static_cast<std::int64_t>(i) + hi);
      if (a1 >= a0) cells += static_cast<std::size_t>(a1 - a0 + 1);
    }
    extra.set("band", Json::array({lo, hi}));
  } else if (algorithm == "myers") {
    if (!p.dna) throw Fault("Myers' bit vectors here count unit edits on DNA");
    std::optional<aln::Alignment> got;
    aln::MyersStats st;
    const auto mode = p.mode == aln::Mode::global ? aln::Mode::global : aln::Mode::semiglobal;
    seconds = detail::timed([&] { st = {}; got = aln::myers_align(p.a, p.b, mode, -1, 4, &st); });
    if (!got) throw Fault("Myers' algorithm found no alignment");
    al = *got;
    cells = st.cells();
    extra.set("unit_costs", true);
    extra.set("distance_unit", -got->score);
    extra.set("passes", st.passes);
  } else if (algorithm == "wfa" || algorithm == "biwfa") {
    // Penalties minimised, in the order of their values; the score follows from the penalty.
    if (p.mode != aln::Mode::global && !(algorithm == "wfa" && p.mode == aln::Mode::semiglobal))
      throw Fault(algorithm == "wfa" ? "the WFA here is global or semi-global" : "BiWFA here is global");
    const auto f = detail::penalties_of(p);
    aln::WavefrontStats st;
    seconds = detail::timed([&] {
      st = {};
      al = algorithm == "wfa" ? *aln::wavefront(p.a, p.b, f.penalties, p.mode, {}, &st) : aln::biwfa(p.a, p.b, f.penalties, &st);
    });
    const std::int64_t pen = -al.score;
    al.score = p.mode == aln::Mode::global ? aln::score_from_penalty(f, pen, n, m) : aln::rescore(al, p.a, p.b, p.scoring);
    cells = st.cells;
    extra.set("penalty", pen);
    extra.set("penalties", Json::object({{"mismatch", f.penalties.mismatch}, {"gap_open", f.penalties.gap_open}, {"gap_extend", f.penalties.gap_extend}}));
    extra.set("memory", st.peak_bytes);
    extra.set("scores_visited", st.scores);
    if (algorithm == "wfa" && flag(request, "with_matrix", false) && (n + 1) * (m + 1) <= 400000) {
      // Every wavefront, as the cells (i, j) = (h − k, h) of its furthest points, by penalty.
      aln::WavefrontRecord rec;
      aln::wavefront(p.a, p.b, f.penalties, p.mode, {}, nullptr, &rec);
      Json fronts = Json::array();
      for (const auto& w : rec.fronts) {
        Json points = Json::array();
        for (std::int64_t k = w.lo; k <= w.hi; k++) {
          const auto h = w.offsets[static_cast<std::size_t>(k - w.lo)];
          points.push(h < 0 ? Json() : Json::array({h - k, h}));
        }
        fronts.push(Json::object({{"score", w.score}, {"points", points}}));
      }
      fronts_of_wfa = fronts;
      extra.set("record", Json::object({{"offsets", rec.offsets_computed}, {"extended", rec.cells_extended}}));
    }
  } else if (algorithm == "striped") {
    if (p.mode != aln::Mode::local) throw Fault("Farrar's striped kernel is Smith–Waterman's: local");
    const aln::StripedProfile profile(p.a, p.scoring);
    aln::StripedStats st;
    seconds = detail::timed([&] { st = {}; al = aln::striped_local(profile, p.b, &st); });
    extra.set("lanes", Json::object({{"bits", st.bits}, {"passes", st.passes}, {"lazy_f", st.lazy_f}, {"width", static_cast<int>(calx::cpu::LANE_BYTES * 8) / std::max(1, st.bits)},
                                     {"isa", calx::cpu::ISA}}));
  } else {
    throw Fault("the schedule " + algorithm + " is not in this engine yet");
  }
  Json out = detail::alignment_json(al, p);
  if (algorithm == "myers") out.set("score", aln::rescore(al, p.a, p.b, p.scoring));
  out.set("algorithm", algorithm);
  out.set("mode", aln::mode_name(p.mode));
  out.set("seconds", seconds);
  out.set("cells", cells);
  out.set("gcups", seconds > 0 ? static_cast<double>(cells) / seconds / 1e9 : 0);
  // The reference: Gotoh's score (equal to morphe's generic sum, which the domain's checks hold it to).
  const auto ref = aln::gotoh_score(p.a, p.b, p.scoring, p.mode);
  const auto mine = aln::rescore(al, p.a, p.b, p.scoring);
  // Penalties and scores order the alignments alike only end to end: an ends-free WFA minimises
  // the penalty, which need not maximise the score (Eizenga and Paten 2022).
  const bool comparable = algorithm != "myers" && !(algorithm == "wfa" && p.mode != aln::Mode::global);
  out.set("reference", Json::object({{"score", ref}, {"same", comparable ? Json(mine == ref) : Json()}}));
  if (flag(request, "with_matrix", false) && (n + 1) * (m + 1) <= 400000) {
    const auto field = aln::gotoh_field(p.a, p.b, p.scoring, p.mode);
    Json values = Json::array();
    for (auto v : field.best) values.push(v < -(std::int64_t(1) << 40) ? Json() : Json(v));
    Json matrix = Json::object({{"rows", n}, {"cols", m}, {"values", values}, {"path", detail::path_json(al)}});
    if (extra.has("band")) matrix.set("band", extra["band"]);
    if (!fronts_of_wfa.is_null()) matrix.set("fronts", fronts_of_wfa);
    out.set("matrix", matrix);
  }
  for (const auto& [k, v] : extra.as_object()) out.set(k, v);
  return out;
}

// Every schedule on one pair: its value, cells, time, and whether it agrees with the reference.
BENCH_OP("align.compare") {
  const auto p = detail::pair_of(request);
  const auto n = p.a.size(), m = p.b.size();
  Json runs = Json::array();
  const auto ref = aln::gotoh_score(p.a, p.b, p.scoring, p.mode);
  const auto add = [&](const char* name, const char* computes, Json value, std::size_t cells, double seconds, Json agrees) {
    runs.push(Json::object({{"name", name}, {"computes", computes}, {"value", value}, {"cells", cells}, {"seconds", seconds},
                            {"gcups", seconds > 0 && cells ? static_cast<double>(cells) / seconds / 1e9 : 0}, {"agrees", agrees}}));
  };
  if (n * m <= 4'000'000) {
    std::int64_t generic = 0;
    const double t = detail::timed([&] { generic = aln::reference_score(p.a, p.b, p.scoring, p.mode); });
    add("morphe’s generic sum", "the definition: (max, +) over every path", generic, n * m, t, generic == ref);
  }
  {
    aln::Alignment al;
    const double t = detail::timed([&] { al = aln::gotoh(p.a, p.b, p.scoring, p.mode); });
    add("Gotoh", "best score, with its path", al.score, n * m, t, al.score == ref);
    std::int64_t s = 0;
    const double t2 = detail::timed([&] { s = aln::gotoh_score(p.a, p.b, p.scoring, p.mode); });
    add("Gotoh, score only", "best score in linear memory", s, n * m, t2, s == ref);
  }
  {
    aln::Alignment al;
    const double t = detail::timed([&] { al = aln::hirschberg(p.a, p.b, p.scoring, p.mode); });
    add("Hirschberg", "best score and path in linear space", al.score, 2 * n * m, t, al.score == ref);
  }
  if (p.mode == aln::Mode::local) {
    const aln::StripedProfile profile(p.a, p.scoring);
    aln::StripedStats st;
    aln::Endpoint e;
    const double t = detail::timed([&] { st = {}; e = profile.end(p.b, &st); });
    const std::string lanes = "best local score on " + std::to_string(calx::cpu::LANE_BYTES * 8 / static_cast<std::size_t>(std::max(1, st.bits))) + " lanes of " +
                              std::to_string(st.bits) + " bits";
    add("Farrar, striped", lanes.c_str(), e.score, n * m, t, e.score == ref);
    aln::Alignment al;
    const double t2 = detail::timed([&] { al = aln::striped_local(profile, p.b); });
    add("Farrar, with its path", "the end, the start from the reversed prefixes, the path in a band", al.score, n * m, t2, al.score == ref);
  }
  const bool acgt = p.dna && seq::count_ambiguous(p.a) + seq::count_ambiguous(p.b) == 0;
  if (acgt && (p.mode == aln::Mode::global || p.mode == aln::Mode::semiglobal)) {
    std::optional<aln::PenaltyForm> f;
    try { f = aln::penalty_form(p.scoring); } catch (const std::invalid_argument&) {}
    if (f) {
      const bool global = p.mode == aln::Mode::global;
      aln::WavefrontStats st;
      std::optional<aln::Alignment> al;
      const double t = detail::timed([&] { st = {}; al = aln::wavefront(p.a, p.b, f->penalties, p.mode, {}, &st); });
      const auto score = global ? aln::score_from_penalty(*f, -al->score, n, m) : aln::rescore(*al, p.a, p.b, p.scoring);
      add("WFA", global ? "least penalty, with its path; the score from it" : "least penalty, ends free (a different order)", score, st.cells, t, global ? Json(score == ref) : Json());
      std::int64_t pen = 0;
      const double t2 = detail::timed([&] { pen = aln::wavefront_penalty(p.a, p.b, f->penalties, p.mode); });
      const auto score2 = global ? aln::score_from_penalty(*f, pen, n, m) : -pen;
      add("WFA, penalty only", global ? "the last wavefronts kept: O(s) memory" : "the least penalty (negated)", score2, st.cells, t2, global ? Json(score2 == ref) : Json());
      if (global) {
        aln::Alignment b;
        aln::WavefrontStats sb;
        const double t3 = detail::timed([&] { sb = {}; b = aln::biwfa(p.a, p.b, f->penalties, &sb); });
        const auto sc = aln::score_from_penalty(*f, -b.score, n, m);
        add("BiWFA", "both ends towards the middle, O(s) memory", sc, sb.cells, t3, sc == ref);
      }
    }
  }
  if (p.dna) {
    std::int64_t d = 0;
    aln::MyersStats st;
    const auto mode = p.mode == aln::Mode::global ? aln::Mode::global : aln::Mode::semiglobal;
    const double t = detail::timed([&] { st = {}; d = aln::myers_distance(p.a, p.b, mode, -1, 4, &st); });
    add("Myers", "unit edit distance (a different scoring)", d, st.cells(), t, Json());
  }
  return Json::object({{"runs", runs}, {"reference", ref}, {"cells", n * m}});
}

// The same edit graph in other semirings.
BENCH_OP("align.semiring") {
  const auto p = detail::pair_of(request, 3000);
  Json out = Json::object();
  const auto count = aln::reference_count(p.a, p.b, p.scoring, p.mode);
  out.set("score", count.score);
  out.set("optimal", count.count.str());
  out.set("alignments", aln::reference_alignments(p.a, p.b, p.mode).str());
  if (p.mode == aln::Mode::global) out.set("delannoy", aln::reference_alignments(p.a, p.b, aln::Mode::global).str());
  // The probability semiring: the pair HMM (Durbin et al. 1998, ch. 4), its forward sum, the
  // posterior of every pair, and the alignment of maximum expected accuracy.
  if (p.a.size() * p.b.size() <= 250'000 && seq::count_ambiguous(p.dna ? p.a : seq::Codes{}) == 0 && (p.dna ? seq::count_ambiguous(p.b) == 0 : true)) {
    aln::PairHmm h;
    if (p.dna) {
      h = aln::PairHmm::dna(0.9, 0.05, 0.4, 0.01);
    } else {
      const auto bg = aln::robinson_frequencies();
      const auto ka = aln::karlin_altschul(p.scoring.matrix, bg, bg);
      h = aln::PairHmm::from_matrix(p.scoring.matrix, bg, ka.lambda, 0.05, 0.4, 0.01);
    }
    const auto post = aln::posterior(h, p.a, p.b);
    out.set("forward", aln::forward(h, p.a, p.b).log_probability);
    out.set("posterior", Json::object({{"rows", p.a.size()}, {"cols", p.b.size()}, {"values", numbers(post, 4)}}));
    const auto mea = aln::mea_alignment(post, p.a.size(), p.b.size());
    out.set("mea", Json::object({{"rows", detail::rows_json(mea.alignment, p)}, {"accuracy", mea.accuracy}, {"expected_pairs", mea.expected_pairs}}));
  }
  return out;
}

// Reads of a genome to score against one of them: the input of the GPU's comparison with the CPU.
BENCH_OP("align.batch_input") {
  const auto g = worker.genome(request);
  std::vector<seq::Codes> contigs;
  for (const auto& c : g->contigs) contigs.push_back(c.codes);
  sim::ShortReadModel model;
  model.length = static_cast<std::size_t>(bounded(request, "length", 150, 20, 1000));
  const auto n = static_cast<std::size_t>(bounded(request, "targets", 4000, 1, 200000));
  const auto reads = short_reads(worker, sim::Source::haploid(contigs), n, model, static_cast<std::uint64_t>(integer(request, "seed", 1)));
  Json targets = Json::array();
  for (const auto& r : reads) targets.push(seq::decode_dna(r.bases));
  const auto q = static_cast<std::size_t>(bounded(request, "query", 120, 10, 256));
  std::string query = seq::decode_dna(reads.front().bases).substr(0, std::min<std::size_t>(q, reads.front().bases.size()));
  return Json::object({{"query", query}, {"targets", targets}});
}

// Local scores of one query against many targets, on the worker's threads: the CPU's side of the
// comparison with the GPU's kernel (src/gpu/smithwaterman.js).
BENCH_OP("align.batch") {
  const seq::Codes query = seq::encode_dna(str(request, "query"));
  std::vector<seq::Codes> targets;
  for (const auto& t : request["targets"].as_array()) targets.push_back(seq::encode_dna(t.text()));
  const auto s = aln::scoring(aln::Matrix::dna(2, -4, -1), static_cast<int>(bounded(request, "gap_open", 4, 0, 100)), static_cast<int>(bounded(request, "gap_extend", 2, 1, 100)));
  std::vector<std::int64_t> scores(targets.size());
  std::size_t cells = 0;
  for (const auto& t : targets) cells += t.size() * query.size();
  // Three schedules of the same scores: Gotoh's scalar rows; Farrar's stripes along the query;
  // one target to a lane (Rognes 2011). Each on every thread of the worker.
  const std::string kernel = str(request, "kernel", "inter");
  const auto threads = std::to_string(worker.pool.size()) + " threads";
  std::string described;
  const auto t0 = std::chrono::steady_clock::now();
  if (kernel == "gotoh") {
    calx::scheduler::parallel_for(worker.pool, targets.size(), 64, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t k = b; k < e; k++) scores[k] = aln::gotoh_score(query, targets[k], s, aln::Mode::local);
    });
    described = "Gotoh's scores in linear memory, " + threads;
  } else if (kernel == "striped") {
    const aln::StripedProfile profile(query, s);
    calx::scheduler::parallel_for(worker.pool, targets.size(), 64, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t k = b; k < e; k++) scores[k] = profile.end(targets[k]).score;
    });
    described = std::string("Farrar's stripes along the query (") + calx::cpu::ISA + "), " + threads;
  } else if (kernel == "inter") {
    calx::scheduler::parallel_for(worker.pool, targets.size(), 512, [&](std::size_t b, std::size_t e, std::size_t) {
      const std::vector<seq::Codes> part(targets.begin() + static_cast<std::ptrdiff_t>(b), targets.begin() + static_cast<std::ptrdiff_t>(e));
      const auto got = aln::inter_scores(query, part, s);
      std::copy(got.begin(), got.end(), scores.begin() + static_cast<std::ptrdiff_t>(b));
    });
    described = "one target to each of " + std::to_string(calx::cpu::LANE_BYTES / 2) + " lanes of 16 bits (" + calx::cpu::ISA + "), " + threads;
  } else {
    throw Fault("a batch kernel is gotoh, striped or inter");
  }
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return Json::object({{"scores", Json(std::vector<double>(scores.begin(), scores.end()))}, {"cells", cells}, {"seconds", seconds},
                       {"gcups", static_cast<double>(cells) / seconds / 1e9}, {"kernel", described}});
}

// Partial-order alignment: noisy long reads of one sequence fused into a graph, read by read, and
// the heaviest bundle through it as their consensus (Lee, Grasso and Sharlow 2002; Lee 2003).
BENCH_OP("align.poa") {
  const auto n = static_cast<std::size_t>(bounded(request, "length", 150, 20, 5000));
  const auto count = static_cast<std::size_t>(bounded(request, "reads", 10, 2, 60));
  const double error = num(request, "error", 0.08);
  if (!(error >= 0 && error <= 0.3)) throw Fault("an error rate is between 0 and 0.3");
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  const auto truth = sim::random_sequence(n, 0.5, seed);
  sim::LongReadModel model;
  model.error = error;
  model.mean_length = 1e9;  // whole reads: each spans the sequence
  model.min_length = n;
  std::vector<seq::Codes> reads;
  for (const auto& r : sim::long_reads(sim::Source::haploid({truth}), count, model, seed + 1)) {
    seq::Codes x = r.bases;
    if (r.reverse) x = seq::reverse_complement(x);
    reads.push_back(std::move(x));
  }
  const auto s = aln::scoring(aln::Matrix::dna(2, -4, -1), 4, 2);
  const auto t0 = std::chrono::steady_clock::now();
  aln::PoaGraph g(reads.front());
  Json growth = Json::array({Json::array({1, g.size(), g.edges()})});
  for (std::size_t k = 1; k < reads.size(); k++) {
    g.add(reads[k], s, aln::Mode::global);
    growth.push(Json::array({k + 1, g.size(), g.edges()}));
  }
  std::vector<std::uint32_t> path;
  const auto consensus = g.consensus(&path);
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  // The graph, for drawing: each node's base, column and weight (the reads through it); each edge.
  Json nodes = Json::array(), edges = Json::array();
  for (std::uint32_t u = 0; u < g.size(); u++)
    nodes.push(Json::array({std::string(1, seq::dna_letter(g.node(u).base)), g.column_of()[u], g.node(u).sequences.size()}));
  for (std::uint32_t e = 0; e < g.edges(); e++) edges.push(Json::array({g.edge(e).from, g.edge(e).to, g.edge(e).weight}));
  Json rows = Json::array();
  for (const auto& r : g.msa(true)) rows.push(seq::decode_dna(r));
  Json errors = Json::array();
  for (const auto& r : reads) errors.push(aln::myers_distance(r, truth));
  return Json::object({{"truth", seq::decode_dna(truth)}, {"consensus", seq::decode_dna(consensus)}, {"distance", aln::myers_distance(consensus, truth)},
                       {"read_distances", errors}, {"nodes", nodes}, {"edges", edges}, {"columns", g.columns()},
                       {"path", Json(std::vector<double>(path.begin(), path.end()))}, {"msa", rows}, {"growth", growth}, {"seconds", seconds}});
}

// The parametric polygon of (mismatches, gap columns), its vertices with their alignments, and
// which vertex is optimal on a grid of (α, β) in the quadrant biology uses.
BENCH_OP("align.parametric") {
  auto p = detail::pair_of(request, 400);
  const auto par = aln::parametric(p.a, p.b, p.dna);
  Json vertices = Json::array();
  for (const auto& v : par.vertices) {
    Json j = Json::object({{"mismatches", v.mismatches}, {"gaps", v.gaps}, {"matches", v.matches}, {"cigar", aln::cigar_string(v.alignment.cigar)}});
    j.set("rows", detail::rows_json(v.alignment, p));
    vertices.push(j);
  }
  // The fan's cones issue from (−1, −1/2), so its structure spreads over orders of magnitude:
  // (α, β) are sampled on logarithmic axes from `low` to `high`.
  const auto grid = static_cast<std::size_t>(bounded(request, "grid", 140, 10, 400));
  const double lo = std::log(num(request, "low", 0.01)), hi = std::log(num(request, "high", 100));
  const auto at = [&](double k) { return std::exp(lo + (hi - lo) * (k + 0.5) / static_cast<double>(grid)); };
  Json cells = Json::array();
  for (std::size_t r = 0; r < grid; r++)       // β from high to low, so rows read top-down
    for (std::size_t c = 0; c < grid; c++) cells.push(par.optimal(at(static_cast<double>(c)), at(static_cast<double>(grid - 1 - r))));
  return Json::object({{"vertices", vertices}, {"grid", grid}, {"low", std::exp(lo)}, {"high", std::exp(hi)}, {"optimal", cells},
                       {"query_length", par.query_length}, {"target_length", par.target_length}});
}

}  // namespace workbench
