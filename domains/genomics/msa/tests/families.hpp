// A test helper: families of protein sequences evolved from a root, with their true alignment,
// by the process of genomics::simulate::evolve_family (simulate/family.hpp) on twenty letters:
// a balanced binary tree, on each branch every site substituted with probability
// `substitution` and an indel of geometric length (mean 2) starting after each site with
// probability `indel`. The root and inserted residues are drawn from Robinson and Robinson's
// frequencies; a substituted residue a becomes b ≠ a with probability proportional to BLOSUM62's
// target frequency q_ab = p_a p_b e^{λ s_ab} (Altschul 1991, J. Mol. Biol. 219, 555–565), λ
// from align/statistics.hpp, so that the substitutions are those BLOSUM62 scores as likely.

#pragma once

#include <genomics/align/scoring.hpp>
#include <genomics/align/statistics.hpp>
#include <genomics/msa/msa.hpp>
#include <genomics/simulate/family.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace test_families {

using genomics::sequence::Code;
using genomics::sequence::Codes;

struct ProteinModel {
  std::array<double, 20> background{};          // cumulative
  std::array<std::array<double, 20>, 20> step{};  // cumulative P(b | a), b ≠ a

  ProteinModel() {
    const auto p = genomics::align::robinson_frequencies();
    const auto m = genomics::align::matrix("BLOSUM62");
    const auto d = genomics::align::score_distribution(m, std::span<const double>(p.data(), p.size()), std::span<const double>(p.data(), p.size()));
    const double lambda = genomics::align::karlin_lambda(d);
    double acc = 0;
    for (int a = 0; a < 20; a++) background[static_cast<std::size_t>(a)] = acc += p[static_cast<std::size_t>(a)];
    for (int a = 0; a < 20; a++) {
      double row = 0;
      std::array<double, 20> q{};
      for (int b = 0; b < 20; b++)
        if (b != a) row += q[static_cast<std::size_t>(b)] = p[static_cast<std::size_t>(b)] * std::exp(lambda * m(static_cast<Code>(a), static_cast<Code>(b)));
      double c = 0;
      for (int b = 0; b < 20; b++) step[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = c += q[static_cast<std::size_t>(b)] / row;
    }
  }
  static Code pick(const std::array<double, 20>& cumulative, double u) {
    for (int a = 0; a < 19; a++)
      if (u < cumulative[static_cast<std::size_t>(a)]) return static_cast<Code>(a);
    return 19;
  }
  Code residue(calx::random::Stream& r) const { return pick(background, r.uniform() * background[19]); }
  Code substitute(calx::random::Stream& r, Code a) const { return pick(step[a], r.uniform() * step[a][19]); }
};

inline Codes random_protein(const ProteinModel& model, std::size_t length, std::uint64_t seed) {
  calx::random::Stream r(seed, 0x50524f54, 0, 0);
  Codes x(length);
  for (auto& c : x) c = model.residue(r);
  return x;
}

// evolve_family's process on protein letters (the same structure, line for line).
inline genomics::simulate::Family evolve_protein_family(const ProteinModel& model, const Codes& root, std::size_t members, double substitution, double indel,
                                                       std::uint64_t seed) {
  calx::random::Stream r(seed, genomics::simulate::FAMILY, 1, 0);
  struct Lineage {
    std::vector<std::size_t> column;
    Codes base;
  };
  std::vector<std::size_t> order;
  std::size_t next_column = 0;
  Lineage start;
  for (Code b : root) {
    start.column.push_back(next_column);
    start.base.push_back(b);
    order.push_back(next_column++);
  }
  const auto branch = [&](const Lineage& parent) {
    Lineage child;
    for (std::size_t i = 0; i < parent.base.size(); i++) {
      Code b = parent.base[i];
      if (r.uniform() < substitution) b = model.substitute(r, b);
      child.base.push_back(b);
      child.column.push_back(parent.column[i]);
    }
    Lineage out;
    for (std::size_t i = 0; i <= child.base.size(); i++) {
      if (i < child.base.size()) {
        out.base.push_back(child.base[i]);
        out.column.push_back(child.column[i]);
      }
      if (r.uniform() >= indel) continue;
      const std::size_t len = 1 + r.geometric(0.5);
      if (r.bernoulli(0.5)) {
        std::size_t after_pos = 0;
        if (!out.column.empty()) after_pos = static_cast<std::size_t>(std::find(order.begin(), order.end(), out.column.back()) - order.begin()) + 1;
        for (std::size_t k = 0; k < len; k++) {
          order.insert(order.begin() + static_cast<std::ptrdiff_t>(after_pos + k), next_column);
          out.column.push_back(next_column++);
          out.base.push_back(model.residue(r));
        }
      } else {
        i += len;
      }
    }
    return out;
  };
  std::vector<Lineage> leaves;
  const auto descend = [&](auto&& self, const Lineage& at, std::size_t m) -> void {
    if (m == 1) { leaves.push_back(at); return; }
    const std::size_t left = m / 2;
    const auto a = branch(at);
    self(self, a, left);
    const auto b = branch(at);
    self(self, b, m - left);
  };
  descend(descend, start, std::max<std::size_t>(1, members));
  genomics::simulate::Family f;
  std::vector<bool> used(next_column, false);
  for (const auto& l : leaves)
    for (std::size_t c : l.column) used[c] = true;
  std::vector<std::size_t> kept;
  for (std::size_t c : order)
    if (used[c]) kept.push_back(c);
  std::vector<std::size_t> index(next_column, 0);
  for (std::size_t k = 0; k < kept.size(); k++) index[kept[k]] = k;
  for (std::size_t m = 0; m < leaves.size(); m++) {
    f.names.push_back("p" + std::to_string(m + 1));
    f.sequences.push_back(leaves[m].base);
    Codes row(kept.size(), genomics::sequence::GAP);
    for (std::size_t i = 0; i < leaves[m].base.size(); i++) row[index[leaves[m].column[i]]] = leaves[m].base[i];
    f.alignment.push_back(std::move(row));
  }
  return f;
}

inline genomics::msa::Msa truth(const genomics::simulate::Family& f) { return {f.names, f.alignment}; }

}  // namespace test_families
