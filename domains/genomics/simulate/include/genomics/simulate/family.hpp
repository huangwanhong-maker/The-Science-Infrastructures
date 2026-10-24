// A family of sequences evolved from a root, with their true alignment: the truth that
// multiple alignment is checked against (the TC and SP scores of Thompson, Plewniak and Poch
// 1999, Nucleic Acids Research 27, 2682–2690, measure how much of it an aligner recovers).
//
// The members descend from the root by a balanced binary tree; on each branch every site is
// substituted with probability `substitution`, and an insertion or a deletion of a length
// geometric with mean 2 starts after each site with probability `indel` (half each). Every base
// keeps the column it was born in: a column is created by the root or by an insertion, placed
// in the global order right after the column of the base before it (insertions on different
// branches at one place get columns of their own, since they are not homologous), so that the
// true alignment is the matrix of members by columns, with a gap where a member lacks a column.

#pragma once

#include <genomics/simulate/genome.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace genomics::simulate {

struct Family {
  std::vector<std::string> names;      // m1, m2, ...
  std::vector<Codes> sequences;        // ungapped
  std::vector<Codes> alignment;        // the true alignment, sequence::GAP in its gaps
  std::string newick;                  // the tree, branch lengths the substitution probability
};

namespace detail {

struct Lineage {
  std::vector<std::size_t> column;  // the column of each base
  Codes base;
};

}  // namespace detail

inline Family evolve_family(const Codes& root, std::size_t members, double substitution, double indel, std::uint64_t seed) {
  calx::random::Stream r(seed, FAMILY, 0, 0);
  std::vector<std::size_t> order;   // the columns, in their global order
  std::size_t next_column = 0;
  detail::Lineage start;
  for (Code b : root) {
    start.column.push_back(next_column);
    start.base.push_back(b);
    order.push_back(next_column++);
  }
  // One branch: substitutions, then indels, on a copy of the parent.
  const auto branch = [&](const detail::Lineage& parent) {
    detail::Lineage child;
    for (std::size_t i = 0; i < parent.base.size(); i++) {
      Code b = parent.base[i];
      if (r.uniform() < substitution) b = substitute(r, b);
      child.base.push_back(b);
      child.column.push_back(parent.column[i]);
    }
    detail::Lineage out;
    for (std::size_t i = 0; i <= child.base.size(); i++) {
      if (i < child.base.size()) {
        out.base.push_back(child.base[i]);
        out.column.push_back(child.column[i]);
      }
      if (r.uniform() >= indel) continue;
      const std::size_t len = 1 + r.geometric(0.5);
      if (r.bernoulli(0.5)) {
        // An insertion after the base just written (or at the very start).
        std::size_t after_pos = 0;  // where in `order` the new columns go
        if (!out.column.empty()) after_pos = static_cast<std::size_t>(std::find(order.begin(), order.end(), out.column.back()) - order.begin()) + 1;
        for (std::size_t k = 0; k < len; k++) {
          order.insert(order.begin() + static_cast<std::ptrdiff_t>(after_pos + k), next_column);
          out.column.push_back(next_column++);
          out.base.push_back(static_cast<Code>(r.below(4)));
        }
      } else {
        i += len;  // the next len bases are deleted
      }
    }
    return out;
  };
  // The balanced tree, members split in halves.
  std::vector<detail::Lineage> leaves;
  std::string newick;
  const auto descend = [&](auto&& self, const detail::Lineage& at, std::size_t m) -> std::string {
    if (m == 1) {
      leaves.push_back(at);
      return "m" + std::to_string(leaves.size());
    }
    const std::size_t left = m / 2;
    const auto a = branch(at);
    const std::string sa = self(self, a, left);
    const auto b = branch(at);
    const std::string sb = self(self, b, m - left);
    char len[32];
    std::snprintf(len, sizeof len, "%.6g", substitution);
    return "(" + sa + ":" + len + "," + sb + ":" + len + ")";
  };
  newick = descend(descend, start, std::max<std::size_t>(1, members)) + ";";
  Family f;
  f.newick = newick;
  // Columns no leaf has are dropped from the alignment.
  std::vector<bool> used(next_column, false);
  for (const auto& l : leaves)
    for (std::size_t c : l.column) used[c] = true;
  std::vector<std::size_t> kept;
  for (std::size_t c : order)
    if (used[c]) kept.push_back(c);
  std::vector<std::size_t> index(next_column, 0);
  for (std::size_t k = 0; k < kept.size(); k++) index[kept[k]] = k;
  for (std::size_t m = 0; m < leaves.size(); m++) {
    f.names.push_back("m" + std::to_string(m + 1));
    f.sequences.push_back(leaves[m].base);
    Codes row(kept.size(), sequence::GAP);
    for (std::size_t i = 0; i < leaves[m].base.size(); i++) row[index[leaves[m].column[i]]] = leaves[m].base[i];
    f.alignment.push_back(std::move(row));
  }
  return f;
}

}  // namespace genomics::simulate
