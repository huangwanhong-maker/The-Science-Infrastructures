// Multiple alignment by partial orders (Lee, Grasso and Sharlow 2002, Bioinformatics 18,
// 452–464): the sequences are added one by one to a partial-order graph (align/poa.hpp), each
// aligned to the whole graph rather than to a profile, and the alignment is the one the graph
// holds, its columns the aligned rings of nodes in an order of the graph.
//
// Unlike a profile, the graph keeps every sequence's own path, so an indel shared by some
// sequences is a branch that later sequences can follow without paying again; but the result
// depends on the order of addition, and each sequence is aligned once, never revised. The checks
// pass a guide tree's leaves from left to right, so that close relatives enter together.
//
// Complexity: the i-th sequence of length L against a graph of V_i nodes and E_i edges,
// O((V_i + E_i) L) time and 9 bytes per cell of traceback (poa.hpp); V grows with the
// divergence, from about L for identical sequences towards n L.

#pragma once

#include <genomics/align/poa.hpp>
#include <genomics/align/types.hpp>
#include <genomics/msa/msa.hpp>

#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::msa {

// The alignment of the sequences added to one graph in `order` (empty: as given); rows in the
// input's order.
inline Msa partial_order_alignment(const std::vector<Codes>& seqs, const std::vector<std::string>& names_in, const align::Scoring& s,
                                   std::vector<std::size_t> order = {}, align::PoaGraph* graph = nullptr) {
  const std::size_t n = seqs.size();
  if (!names_in.empty() && names_in.size() != n) throw std::invalid_argument("genomics::msa::partial_order_alignment: one name per sequence");
  if (order.empty()) {
    order.resize(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
  }
  std::vector<bool> seen(n, false);
  if (order.size() != n) throw std::invalid_argument("genomics::msa::partial_order_alignment: the order is not a permutation of the sequences");
  for (std::size_t i : order) {
    if (i >= n || seen[i]) throw std::invalid_argument("genomics::msa::partial_order_alignment: the order is not a permutation of the sequences");
    seen[i] = true;
  }
  align::PoaGraph g;
  for (std::size_t i : order) g.add(seqs[i], s, align::Mode::global);
  const std::vector<Codes> rows = g.msa();
  Msa out;
  out.names.resize(n);
  out.rows.resize(n);
  for (std::size_t k = 0; k < n; k++) {
    out.rows[order[k]] = rows[k];
    out.names[order[k]] = names_in.empty() ? "s" + std::to_string(order[k] + 1) : names_in[order[k]];
  }
  if (graph) *graph = std::move(g);
  return out;
}

}  // namespace genomics::msa
