// Where value settles: each economy passes its stock to its partners in the shares
// of its exports, and the stocks run to the stationary ones of that rule. A C++
// plug-in that uses the economics domain itself: its stocks under a stochastic
// rule (domains/economics/theory/value-in-space-and-time.md §3.3), the user's
// warning that accumulation is not exploitation made a map.
//
//   M_ji = F_ij / X_i     the share of i's exports that goes to j (column i sums to 1);
//                         an economy that exports nothing keeps its stock, M_ii = 1;
//   M' = (1 − q) M + q/n  damped, so that M' is irreducible and aperiodic and its
//                         stationary stocks unique (q, as PageRank's teleportation);
//   M' π = π, Σ π = n     economics::stationary: π_i above 1 is a place where value settles.

#include "network_plugin.hpp"

#include <economics/circulation/stocks.hpp>

#include <algorithm>
#include <numeric>

namespace net = science::network;
using science::Json;

struct Stationary : net::Analysis {
  Json describe() const override {
    return Json::object({
        {"id", "stationary"},
        {"title", "Where value settles"},
        {"description", "Each economy passes its stock to its partners in the shares of its exports; the stationary "
                        "stocks of that rule (the economics domain's), damped by q as PageRank is. Above 1, value settles."},
        {"params", Json::array({net::param("q", "Damping q", 0.15, "number", 0.0, 1.0, 0.01, {},
                                           "the share of each stock spread evenly over all economies every period")})},
    });
  }

  std::vector<Json> run(const net::Network& network, const net::Params& params, net::Progress& progress) override {
    const std::size_t n = network.size();
    if (n < 2) throw std::invalid_argument("a network of at least two places");
    const double q = params.number("q");
    const auto F = network.matrix();
    const auto X = network.exports();
    calx::linalg::Matrix M(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) {
        const double passed = X[i] > 0 ? F[i][j] / X[i] : (i == j ? 1.0 : 0.0);
        M(j, i) = (1 - q) * passed + q / static_cast<double>(n);
      }
    progress(0.5, "solving M' pi = pi");
    const auto pi = economics::stationary(M, static_cast<double>(n));
    std::map<std::string, double> by_place;
    for (std::size_t i = 0; i < n; i++) by_place[network.places[i].id] = pi[i];
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return pi[a] > pi[b]; });
    std::vector<std::vector<Json>> rows;
    for (std::size_t k = 0; k < std::min<std::size_t>(15, n); k++) {
      const std::size_t i = order[k];
      rows.push_back({network.places[i].name, pi[i], X[i]});
    }
    return {
        net::values("Stationary stocks (mean 1)", by_place, "", "sequential",
                    "Above 1: value passed on in the shares of exports settles here."),
        net::table("Where value settles most", {"economy", "stationary stock", "exports"}, rows,
                   "economics::stationary on M' = (1 − q) M + q/n, q = " + Json::number_text(q) + "."),
        net::note("Stocks drift with no loop and no deficit: where value settles is not where it is extracted "
                  "(the theory note, §3.3)."),
    };
  }
};

SCIENCE_NETWORK_PLUGIN(Stationary)
