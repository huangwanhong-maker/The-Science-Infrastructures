// What the Expression view (expression.quantify, on a simulated sample) and the toolbox
// (tool.expr.quant, on the user's transcripts and reads) do alike to quantify a sample: the
// transcriptome's index of k-mers and their classes, the sample pseudoaligned to it, effective
// lengths from the law of fragment lengths, and the abundances by EM (Bray et al. 2016), each step
// timed. They differ in how the reads reach the index, so that step is the caller's.

#pragma once

#include "genomes.hpp"

#include <genomics/expression/pseudoalign.hpp>
#include <genomics/expression/quantify.hpp>
#include <genomics/expression/transcriptome.hpp>

#include <calx/scheduler/thread_pool.hpp>

#include <chrono>
#include <vector>

namespace workbench {

namespace gex = genomics::expression;

namespace detail {

struct Quantification {
  gex::Pseudoalignment pa;
  std::vector<double> effective;
  gex::Abundances abundances;
  std::size_t kmers = 0, index_bytes = 0;
  double index_seconds = 0, pseudoalign_seconds = 0, em_seconds = 0;
};

// pseudoalign(index) gives the sample's pseudoalignment against the index built here.
template <class Pseudoalign>
Quantification quantify_sample(const gex::Transcriptome& tx, const gex::FragmentLengths& d, calx::scheduler::ThreadPool& pool, Pseudoalign&& pseudoalign, int k = 31) {
  using clock = std::chrono::steady_clock;
  const auto since = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };
  Quantification q;
  auto t = clock::now();
  const gex::Index index(tx, k);
  q.index_seconds = since(t);
  q.kmers = index.kmers();
  q.index_bytes = index.memory_bytes();
  t = clock::now();
  q.pa = pseudoalign(index);
  q.pseudoalign_seconds = since(t);
  q.effective = gex::effective_lengths(tx.lengths(), d);
  t = clock::now();
  q.abundances = gex::quantify(q.pa, q.effective, gex::EmOptions{}, pool);
  q.em_seconds = since(t);
  return q;
}

}  // namespace detail

}  // namespace workbench
