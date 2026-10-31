// Resources carried and converted among places (theory/landscape-dynamics.md §9): the balance
// ṁ = (B ⊗ I_s) j + N r + u of s typed stocks at n places, the combinations it conserves exactly, and the
// difference between moving a thing, converting it, copying what is known, and entering a claim in a ledger.
//
// The balance. m the stocks, place-major (place i, type α at i·s + α); j the flows of
// each type along each directed link, link-major; B the incidence of the links (−1 at the tail, +1 at the head);
// N the stoichiometric matrix of the conversion channels, a column each, r their rates; u the inputs from outside.
// A weighting w ∈ ℤ^{ns} is conserved by every transport and conversion exactly when wᵀ(B ⊗ I_s) = 0 and wᵀN = 0,
// and then d(wᵀm)/dt = wᵀu. The conserved weightings form a lattice, the left kernel over ℤ of
// [B ⊗ I_s | N], found here from the Smith normal form of its transpose (morphe), as physis finds a reaction
// network's conservation laws.
//
// Kinds of change. Moving a good along a link conserves the good. A conversion channel conserves what its column
// is orthogonal to: one unit of A made into two of B conserves 2A + B. Copying
// what is known — the source keeps it, the recipient gains it, ν⁻ = (1, 0), ν⁺ = (1, 1) — conserves no positive
// count of holders: knowledge is not a stock. A sale entered in a ledger creates offsetting claims, an asset and a
// liability, and conserves their difference, not their number. The value of any of these to someone is a
// valuation of the state, not a balance, and is left to the caller. Checks: tests/check_landscape.cpp.

#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/algebra/smith.hpp>

namespace urban::landscape {

struct Resources {
  std::size_t places = 0, types = 0;
  std::vector<std::string> type_names;
  std::vector<std::array<std::size_t, 2>> links;  // (tail, head)
  std::vector<std::vector<long long>> channels;   // each a column of N, of length places · types
  std::vector<bool> transported;                  // which types the links carry (all when empty)

  std::size_t index(std::size_t place, std::size_t type) const { return place * types + type; }
  bool carries(std::size_t type) const { return transported.empty() || transported.at(type); }
  // [B ⊗ I_s | N] over ℤ, the columns of the transported types' link moves first.
  morphe::algebra::Matrix moves() const {
    const std::size_t rows = places * types;
    std::vector<std::vector<long long>> cols;
    for (const auto& [t, h] : links)
      for (std::size_t a = 0; a < types; a++) {
        if (!carries(a)) continue;
        std::vector<long long> c(rows, 0);
        c[index(t, a)] = -1, c[index(h, a)] = 1;
        cols.push_back(c);
      }
    for (const auto& c : channels) {
      if (c.size() != rows) throw std::invalid_argument("urban::landscape::Resources: a channel needs a change for each place and type");
      cols.push_back(c);
    }
    morphe::algebra::Matrix M(rows, cols.size());
    for (std::size_t j = 0; j < cols.size(); j++)
      for (std::size_t i = 0; i < rows; i++) M(i, j) = morphe::algebra::Integer(cols[j][i]);
    return M;
  }
  // A basis over ℤ of the conserved weightings {w : wᵀ [B ⊗ I | N] = 0}.
  std::vector<std::vector<long long>> conserved() const {
    const auto M = moves();
    const std::size_t rows = places * types;
    std::vector<std::vector<long long>> out;
    if (M.cols == 0) {
      for (std::size_t i = 0; i < rows; i++) {
        std::vector<long long> w(rows, 0);
        w[i] = 1;
        out.push_back(w);
      }
      return out;
    }
    const auto S = morphe::algebra::smith(M.transpose(), {.V = true});
    for (std::size_t j = S.rank(); j < rows; j++) {
      std::vector<long long> w(rows);
      for (std::size_t i = 0; i < rows; i++) w[i] = S.V(i, j).to_int64();
      out.push_back(std::move(w));
    }
    return out;
  }
  // ṁ for flows j (one for each carried type on each link, in moves() order), channel rates r and inputs u.
  std::vector<double> rate(const std::vector<double>& j, const std::vector<double>& r, const std::vector<double>& u) const {
    const auto M = moves();
    std::vector<double> x(j);
    x.insert(x.end(), r.begin(), r.end());
    if (x.size() != M.cols || u.size() != M.rows) throw std::invalid_argument("urban::landscape::Resources::rate: one flow per carried type and link, one rate per channel, one input per stock");
    std::vector<double> d(u);
    for (std::size_t i = 0; i < M.rows; i++)
      for (std::size_t k = 0; k < M.cols; k++) d[i] += static_cast<double>(M(i, k).to_int64()) * x[k];
    return d;
  }
};

}  // namespace urban::landscape
