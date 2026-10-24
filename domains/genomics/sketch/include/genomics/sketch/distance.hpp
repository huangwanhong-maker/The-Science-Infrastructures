// The all-against-all distance matrix of many genomes' sketches, as a task graph of blocks.
//
// The n × n matrix is cut into blocks of b × b genomes; the blocks on and above the diagonal are
// the tasks of a calx::scheduler::TaskGraph, each with the cost of its pairs, run by priority of
// upward rank (HEFT; Topcuoglu, Hariri and Wu 2002) on the pool. A task computes the distance of
// each of its pairs once and writes both (i, j) and (j, i); no two tasks write one cell, and the
// distance is a pure function of two sketches, so the matrix is the same, to the bit, on any number
// of threads, whatever the schedule the trace records. O(n² s) for sketches of size s (the Mash
// distance), the tasks independent, so the time is about that over the number of threads.
//
// Any symmetric distance of two sketches may be given; mash_distances() uses the Mash distance of
// minhash.hpp. The diagonal is the distance of each sketch to itself, 0 for a distance.

#pragma once

#include <genomics/sketch/minhash.hpp>

#include <calx/scheduler/task_graph.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::sketch {

struct DistanceMatrix {
  std::size_t n = 0;
  std::vector<double> values;          // row-major, n × n
  std::size_t block = 0, tasks = 0;
  calx::scheduler::Trace trace;        // where and when each block ran

  double operator()(std::size_t i, std::size_t j) const { return values[i * n + j]; }
};

template <class Sketch, class Distance>
DistanceMatrix distance_matrix(calx::scheduler::ThreadPool& pool, const std::vector<Sketch>& sketches, Distance distance, std::size_t block = 32) {
  if (block == 0) throw std::invalid_argument("genomics::sketch::distance_matrix: blocks of at least one genome");
  DistanceMatrix d;
  d.n = sketches.size();
  d.block = block;
  d.values.assign(d.n * d.n, 0);
  const std::size_t blocks = (d.n + block - 1) / block;
  calx::scheduler::TaskGraph graph;
  for (std::size_t bi = 0; bi < blocks; bi++)
    for (std::size_t bj = bi; bj < blocks; bj++) {
      const std::size_t i0 = bi * block, i1 = std::min(d.n, i0 + block), j0 = bj * block, j1 = std::min(d.n, j0 + block);
      const double pairs = bi == bj ? static_cast<double>((i1 - i0) * (i1 - i0 + 1) / 2) : static_cast<double>((i1 - i0) * (j1 - j0));
      graph.add("block " + std::to_string(bi) + "," + std::to_string(bj),
                [&, i0, i1, j0, j1](std::size_t) {
                  for (std::size_t i = i0; i < i1; i++)
                    for (std::size_t j = std::max(j0, i); j < j1; j++) {
                      const double x = distance(sketches[i], sketches[j]);
                      d.values[i * d.n + j] = x;
                      d.values[j * d.n + i] = x;
                    }
                },
                {}, pairs);
    }
  d.tasks = graph.size();
  d.trace = graph.run(pool);
  return d;
}

// The Mash distances of bottom-s sketches.
inline DistanceMatrix mash_distances(calx::scheduler::ThreadPool& pool, const std::vector<MinHash>& sketches, std::size_t block = 32) {
  return distance_matrix(pool, sketches, [](const MinHash& a, const MinHash& b) { return compare(a, b).distance; }, block);
}

}  // namespace genomics::sketch
