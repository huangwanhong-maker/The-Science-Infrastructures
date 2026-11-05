// Checks of capillary networks (capillary.hpp).
//
// - A perfect maze of 12 × 12 cells is a spanning tree: 143 tubes.
// - Driven from cell 0 to cell 143, its flowing tubes are exactly its solution, the unique path found by a
//   breadth-first search; every dead end carries no flow (1e−9 of the peak).
// - The flow through the path is Δp × g / (the path's length in tubes), equal tubes in series.
// - Washburn's filling time 2 μ L²/(γ r) for each tube.

#include <physical_computing/fluids/capillary.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <queue>
#include <vector>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  const std::size_t W = 12;
  const auto maze = fluids::CapillaryNetwork::maze(W, W, 7);
  s.require("maze: a spanning tree of 144 cells, 143 tubes", maze.tubes().size() == W * W - 1);
  const double dp = 100;
  const auto sol = maze.flow(0, W * W - 1, dp);
  auto flowing = maze.flowing(sol);
  // The path by breadth-first search.
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> adj(W * W);
  for (std::size_t k = 0; k < maze.tubes().size(); k++) adj[maze.tubes()[k].a].push_back({maze.tubes()[k].b, k}), adj[maze.tubes()[k].b].push_back({maze.tubes()[k].a, k});
  std::vector<std::size_t> via(W * W, static_cast<std::size_t>(-1));
  std::vector<char> seen(W * W, 0);
  std::queue<std::size_t> q;
  q.push(0), seen[0] = 1;
  while (!q.empty()) {
    const auto v = q.front();
    q.pop();
    for (auto [u, k] : adj[v])
      if (!seen[u]) seen[u] = 1, via[u] = k, q.push(u);
  }
  std::vector<std::size_t> path;
  for (std::size_t v = W * W - 1; v != 0;) {
    const auto k = via[v];
    path.push_back(k);
    v = maze.tubes()[k].a == v ? maze.tubes()[k].b : maze.tubes()[k].a;
  }
  std::sort(path.begin(), path.end()), std::sort(flowing.begin(), flowing.end());
  s.require("maze: the flowing tubes are exactly the solution", flowing == path, std::to_string(flowing.size()) + " flowing, path of " + std::to_string(path.size()));
  const auto net = maze.conductances();
  const double g = net.edges[0].conductance;
  s.check("maze: the flow through the path, Δp g / (its length) (relative)", std::fabs(sol.flow[path[0]]) / (dp * g / path.size()) - 1, 1e-10);
  const auto t = maze.filling_times();
  const auto w = fluids::Fluid::water();
  s.check("Washburn: each tube fills in 2 μ L²/(γ r) (relative)", t[5] / (2 * w.density * w.kinematic_viscosity * 4e-6 / (w.surface_tension * 1e-4)) - 1, 1e-12);
  return s.done();
}
