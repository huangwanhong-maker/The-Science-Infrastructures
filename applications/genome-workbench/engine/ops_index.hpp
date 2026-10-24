// The operations of the Index view: a word's suffix array, LCP, sorted rotations and BWT, the
// steps of backward search, its Lyndon factorisation and bijective transform (morphe's formal
// module); and a genome's FMD-index (the genomics index module): patterns counted and located on
// both strands, a read's super-maximal exact matches, and reads to try.

#pragma once

#include "protocol.hpp"
#include "genomes.hpp"

#include <genomics/index/fmd.hpp>
#include <genomics/simulate/reads.hpp>

#include <morphe/formal/bwt.hpp>
#include <morphe/formal/suffix_array.hpp>
#include <morphe/formal/words.hpp>

#include <algorithm>
#include <chrono>
#include <set>
#include <string>
#include <vector>

namespace workbench {

namespace mf = morphe::formal;

BENCH_OP("index.word") {
  const std::string text = str(request, "text");
  if (text.empty() || text.size() > 400) throw Fault("a word of 1 to 400 letters");
  if (text.find('$') != std::string::npos) throw Fault("$ is the sentinel, and cannot be a letter of the word");
  // The alphabet: the word's letters in byte order, as 0 … σ − 1.
  const std::set<char> letters(text.begin(), text.end());
  const std::string alphabet(letters.begin(), letters.end());
  const auto code = [&](char c) { return static_cast<mf::Letter>(alphabet.find(c)); };
  mf::Word w;
  for (char c : text) w.push_back(code(c));
  const std::size_t n = w.size(), sigma = alphabet.size();
  const auto sa = mf::suffix_array(w, sigma);
  const auto lcp = mf::lcp_array(w, sa);
  const mf::Word L = mf::bwt(w, sa);  // of w$, $ = 0
  const auto letter = [&](mf::Letter c) { return c == 0 ? '$' : alphabet[c - 1]; };
  std::string bwt;
  for (auto c : L) bwt += letter(c);
  // The rows of w$'s sorted suffixes: the suffix $ first, then w's in order.
  std::vector<std::size_t> rows{n};
  for (auto s : sa) rows.push_back(s);
  Json sa_json = Json::array(), lcp_json = Json::array(), rotations = Json::array();
  const std::string ws = text + "$";
  for (std::size_t r = 0; r <= n; r++) {
    sa_json.push(rows[r]);
    lcp_json.push(r <= 1 ? 0 : lcp[r - 1]);
    if (n <= 80) rotations.push(ws.substr(rows[r]) + ws.substr(0, rows[r]));
  }
  std::size_t lcp_sum = 0;
  for (auto v : lcp) lcp_sum += v;
  Json out = Json::object({{"text", text}, {"word", text}, {"alphabet", alphabet}, {"sa", sa_json}, {"lcp", lcp_json}, {"bwt", bwt}, {"runs", mf::runs(L)},
                           {"distinct", n * (n + 1) / 2 - lcp_sum}});
  if (n <= 80) out.set("rotations", rotations);
  // The text again, from the BWT alone.
  {
    const mf::Word back = mf::inverse_bwt(L);
    std::string s;
    for (auto c : back) s += alphabet[c];
    out.set("inverse", s);
  }
  // Lyndon factorisation and the bijective BWT.
  const auto starts = mf::lyndon_factorisation(w);
  Json lyndon = Json::array();
  for (std::size_t f = 0; f + 1 < starts.size(); f++) lyndon.push(Json::array({starts[f], starts[f + 1]}));
  out.set("lyndon", lyndon);
  {
    std::string s;
    for (auto c : mf::bijective_bwt(w)) s += alphabet[c];
    out.set("bijective", s);
  }
  // Backward search: C(c) + occ(c, lo) to C(c) + occ(c, hi), one letter at a time from the right.
  const std::string pattern = str(request, "pattern");
  if (!pattern.empty() && pattern.size() <= 200) {
    const auto C = mf::letter_starts(L, sigma + 1);
    const auto occ = [&](mf::Letter c, std::size_t i) { return static_cast<std::size_t>(std::count(L.begin(), L.begin() + static_cast<std::ptrdiff_t>(i), c)); };
    std::size_t lo = 0, hi = L.size();
    Json steps = Json::array();
    for (std::size_t k = pattern.size(); k-- > 0;) {
      const auto at = alphabet.find(pattern[k]);
      if (at == std::string::npos) { lo = hi = 0; }
      else {
        const auto c = static_cast<mf::Letter>(at + 1);
        lo = C[c] + occ(c, lo);
        hi = C[c] + occ(c, hi);
      }
      steps.push(Json::object({{"letter", std::string(1, pattern[k])}, {"lo", lo}, {"hi", hi}}));
      if (lo >= hi) break;
    }
    out.set("pattern", pattern);
    out.set("search", steps);
  }
  return out;
}

namespace detail {

inline std::shared_ptr<genomics::index::FmdIndex> fmd_of(Worker& worker, const Genome& g, double* seconds = nullptr) {
  return worker.derived<genomics::index::FmdIndex>(g, "fmd", [&] {
    std::vector<seq::Codes> contigs;
    for (const auto& c : g.contigs) contigs.push_back(c.codes);
    return genomics::index::FmdIndex(contigs);
  }, seconds);
}

inline Json index_json(const genomics::index::FmdIndex& f, const Genome& g, double build) {
  return Json::object({{"length", g.length()}, {"bwt_length", f.size()}, {"bytes", f.bytes()}, {"sa_sample", f.options().sample}, {"build_seconds", build}});
}

}  // namespace detail

BENCH_OP("index.count") {
  const auto g = worker.genome(request);
  const std::string text = str(request, "pattern");
  if (text.empty() || text.size() > 10000) throw Fault("a pattern of 1 to 10,000 bases");
  const seq::Codes p = seq::encode_dna(text);
  if (seq::count_ambiguous(p)) throw Fault("a pattern is written in A, C, G and T");
  double build = 0;
  const auto f = detail::fmd_of(worker, *g, &build);
  const auto t0 = std::chrono::steady_clock::now();
  const auto x = f->search(p);
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  Json located = Json::array();
  for (const auto& h : f->locate(x, p.size(), 50)) located.push(Json::array({h.contig, h.position, h.reverse}));
  return Json::object({{"count", x.s}, {"steps", p.size()}, {"seconds", seconds}, {"located", located}, {"palindrome", seq::reverse_complement(p) == p},
                       {"index", detail::index_json(*f, *g, build)}});
}

BENCH_OP("index.sample_read") {
  const auto g = worker.genome(request);
  std::vector<seq::Codes> contigs;
  for (const auto& c : g->contigs) contigs.push_back(c.codes);
  sim::ShortReadModel m;
  m.length = static_cast<std::size_t>(bounded(request, "length", 250, 30, 5000));
  m.error_first = m.error_last = num(request, "error", 0.03);
  m.indel = num(request, "indel", 0.002);
  const auto reads = sim::short_reads(sim::Source::haploid(contigs), 1, m, static_cast<std::uint64_t>(integer(request, "seed", 1)));
  if (reads.empty()) throw Fault("the genome is shorter than a read");
  const auto& r = reads.front();
  return Json::object({{"read", seq::decode_dna(r.bases)}, {"contig", r.contig}, {"begin", r.begin}, {"end", r.end}, {"reverse", r.reverse}, {"errors", r.errors}});
}

BENCH_OP("index.smems") {
  const auto g = worker.genome(request);
  const seq::Codes q = seq::encode_dna(str(request, "read"));
  if (q.size() < 10 || q.size() > 100000) throw Fault("a read of 10 to 100,000 bases");
  const auto min_length = static_cast<std::size_t>(bounded(request, "min_length", 19, 1, 1000));
  double build = 0;
  const auto f = detail::fmd_of(worker, *g, &build);
  const auto t0 = std::chrono::steady_clock::now();
  const auto smems = f->smems(q, min_length);
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  Json list = Json::array();
  for (const auto& m : smems) {
    Json hits = Json::array();
    for (const auto& h : f->locate(m.interval, m.length(), 5)) hits.push(Json::array({h.contig, h.position, h.reverse}));
    list.push(Json::object({{"begin", m.begin}, {"end", m.end}, {"hits", m.interval.s}, {"where", hits}}));
  }
  return Json::object({{"smems", list}, {"min_length", min_length}, {"seconds", seconds}, {"index", detail::index_json(*f, *g, build)}});
}

}  // namespace workbench
