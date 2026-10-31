// The crowd as a kinematic wave (theory/street-oceanography.md, Representation 6.5): a representation, its
// standing marked by Withholding 2.1 of Street Oceanography (a crowd is not asserted to be a continuum).
//
// A street of blocks, each a link of physis's kinematic network (Godunov's scheme on Lighthill–Whitham–
// Richards), joined by crossings, each a short link whose flux law a signal switches: green, the walkway's
// law; red, a law of negligible speed, which passes nothing. Pedestrians' fundamental diagram is
// Greenshields's form with Weidmann's (1993, Transporttechnik der Fussgänger, ETH Zürich, IVT 90, §4) free
// speed 1.34 m/s and jam density 5.4 m⁻², times the walkable width: q(ρ) = v ρ (1 − ρ/ρ_j), ρ per metre of
// street. A rule that opens a carriageway (Ginza's Chūō-dōri, Saturdays, Sundays and holidays from 12:00 to
// 18:00, or 17:00 from October to March; Huang, Towards Street Oceanography, §11) widens every block while it
// holds and stops the crossings' signals; the walkable widths and the demand are declared.
//
// What it shows: the same arrivals converted into unlike densities, queues and times on the street under the
// two configurations (Street Oceanography §7.2, the reversible alteration of terrain by rule); the waves of
// a queue spilling back from a red crossing at the Rankine–Hugoniot speed; and the conservation of the crowd
// exactly to rounding (the budget).
//
// A street of its own: the blocks' lengths one by one (block_lengths), the cross street ending each, and
// whether a signal controls each crossing (signalled); an unsignalled crossing always passes walkers.
//
// A window of the day can be kept at a fine step, for drawing the waves second by second: its densities, the
// flow through every interface, the walkers' speed there, the walkable width and each crossing's signal. In
// the window the steps end at its snapshots. The walkers' speed at an interface is its flow over the
// upwind cell's density, F/ρ, which is v(ρ) where the flow is the upwind cell's own and zero where a red
// crossing passes nothing. Where the upwind cell is empty it is the speed v(0) of a walker alone. The speed at
// the street's start uses the first cell's density.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

#include <physis/continuum/kinematic_waves.hpp>

namespace urban::flow {

namespace pc = physis::continuum;

struct StreetDay {
  double block_length = 125;      // m
  std::size_t blocks = 8;         // the declared street of [G] §11: eight blocks
  std::vector<double> block_lengths;   // m, block by block (empty: `blocks` blocks of block_length)
  std::vector<char> signalled;         // per crossing, whether a signal controls it (empty: every one)
  double cell = 12.5;             // m
  double footway = 12;            // walkable width on an ordinary day, m (declared)
  double carriageway = 27;        // added while the rule holds, m ([G] §11.2)
  double crossing = 20;           // m, the length of a crossing link
  double cycle = 120, green = 0.5;   // the crossings' signal: period (s) and green share
  double free_speed = 1.34;       // m/s (Weidmann 1993)
  double jam = 5.4;               // m⁻² (Weidmann 1993)
  bool rule = false;              // the carriageway opened
  double rule_from = 12, rule_to = 18;   // hours of the day
  double from = 10, to = 20;      // the simulated hours of the day
  std::vector<double> demand;     // arrivals per second at the street's start, one per hour from 0 to 23
  double step = 0;                // s; 0: 0.9 of the CFL bound
  double snapshot = 300;          // s between stored density profiles
  double window_from = 0, window_to = 0;   // s since midnight: a window of fine snapshots (none if not after window_from)
  double window_step = 3;                  // s between them
};

struct StreetDayResult {
  std::vector<double> times;                      // s since midnight, of the snapshots
  std::vector<double> position;                   // m, the centres of the cells along the street (blocks and crossings)
  std::vector<std::vector<double>> density;       // per snapshot, per cell: people per m² (density / width)
  std::vector<double> entered, left;              // cumulative counts at the snapshots
  double budget_residual = 0;                     // entered − left − stored change, people
  double max_density = 0;                         // m⁻²
  double max_density_in_rule_hours = 0;           // m⁻², over the snapshots within [rule_from, rule_to) (on either day)
  double mean_time = 0;                           // s, from the cumulative curves (FIFO)
  double queue_spill = 0;                         // m: the longest run of the street at half the jam density or more
  double throughput = 0;                          // people who left
  // The window: snapshots every window_step from window_from to window_to inclusive (within the simulated hours).
  std::vector<double> interfaces;                 // m, the interfaces' positions along the street (cells + 1)
  std::vector<double> crossings;                  // m, where each crossing begins
  double crossing_length = 0;                     // m
  std::vector<double> window_times;               // s since midnight
  std::vector<std::vector<double>> window_density;   // per snapshot, per cell: people per m²
  std::vector<std::vector<double>> window_flow;      // per snapshot, per interface: people per second
  std::vector<std::vector<double>> window_speed;     // per snapshot, per interface: the walkers' speed, m/s
  std::vector<double> window_width;               // m, walkable at each snapshot
  std::vector<std::vector<char>> window_green;    // per snapshot, per crossing: whether it passes walkers
};

inline StreetDayResult street_day(const StreetDay& d) {
  if (d.demand.size() != 24) throw std::invalid_argument("urban::flow::street_day: a demand for each hour of the day");
  const double width_day = d.footway, width_rule = d.footway + d.carriageway;
  auto open = [&](double t) { const double h = t / 3600; return d.rule && h >= d.rule_from && h < d.rule_to; };
  auto walkway = [&](double width) { return pc::Flux::greenshields(d.free_speed, d.jam * width); };
  const auto red = pc::Flux::greenshields(1e-9, d.jam * width_day);
  // Nodes 0 (start), then alternately block and crossing links.
  std::vector<pc::Link> links;
  std::size_t node = 0;
  std::vector<bool> is_crossing;
  const std::size_t nb = d.block_lengths.empty() ? d.blocks : d.block_lengths.size();
  if (nb == 0) throw std::invalid_argument("urban::flow::street_day: a block or more");
  if (!d.signalled.empty() && d.signalled.size() != nb) throw std::invalid_argument("urban::flow::street_day: one signal flag for each crossing");
  std::vector<char> signal_of;   // per link: whether a signal controls it
  for (std::size_t b = 0; b < nb; b++) {
    const double length = d.block_lengths.empty() ? d.block_length : d.block_lengths[b];
    if (!(length > 0)) throw std::invalid_argument("urban::flow::street_day: blocks of positive length");
    links.push_back({node, node + 1, length, static_cast<std::size_t>(std::max(1.0, std::round(length / d.cell))), walkway(width_day)});
    signal_of.push_back(0);
    is_crossing.push_back(false);
    node++;
    links.push_back({node, node + 1, d.crossing, static_cast<std::size_t>(std::max(1.0, std::round(d.crossing / d.cell))), walkway(width_day)});
    is_crossing.push_back(true);
    signal_of.push_back(d.signalled.empty() || d.signalled[b] ? 1 : 0);
    node++;
  }
  pc::KinematicNetwork net(node + 1, links);
  net.entry(0, [&](double t) { const auto h = static_cast<std::size_t>(std::clamp(t / 3600, 0.0, 23.999)); return d.demand[h]; });
  net.exit(node);
  for (std::size_t n = 1; n < node; n++) net.junction(n, {{1.0}});
  for (std::size_t l = 0; l < links.size(); l++) {
    if (is_crossing[l] && signal_of[l]) {
      net.schedule(l, [&, w = width_day](double t) {
        if (open(t)) return walkway(width_rule);
        const double phase = std::fmod(t, d.cycle) / d.cycle;
        return phase < d.green ? walkway(w) : red;
      });
    } else {
      net.schedule(l, [&](double t) { return walkway(open(t) ? width_rule : width_day); });
    }
  }
  StreetDayResult out;
  double x = 0;
  for (std::size_t l = 0; l < links.size(); l++) {
    const auto& L = links[l];
    const double dx = L.length / static_cast<double>(L.cells);
    if (is_crossing[l]) out.crossings.push_back(x);
    for (std::size_t c = 0; c < L.cells; c++) {
      out.position.push_back(x + (static_cast<double>(c) + 0.5) * dx);
      out.interfaces.push_back(x + static_cast<double>(c) * dx);
    }
    x += L.length;
  }
  out.interfaces.push_back(x);
  out.crossing_length = d.crossing;
  const bool windowed = d.window_to > d.window_from && d.window_step > 0;
  double next_window = std::max(d.window_from, d.from * 3600);
  auto keep_window = [&](const pc::KinematicNetwork::State& s) {
    const double width = open(s.t) ? width_rule : width_day;
    const auto F = net.link_flows(s);
    std::vector<double> rho, flow, speed;
    for (std::size_t l = 0; l < links.size(); l++)
      for (std::size_t c = 0; c < links[l].cells; c++) rho.push_back(s.density[l][c]);
    for (std::size_t l = 0; l < links.size(); l++) {
      for (std::size_t c = 0; c < links[l].cells; c++) flow.push_back(F[l][c]);
      if (l + 1 == links.size()) flow.push_back(F[l][links[l].cells]);
    }
    for (std::size_t i = 0; i < flow.size(); i++) {
      const double up = rho[i == 0 ? 0 : i - 1];
      speed.push_back(i > 0 && up > 1e-9 ? std::clamp(flow[i] / up, 0.0, d.free_speed) : d.free_speed * std::max(0.0, 1 - up / (d.jam * width)));
    }
    std::vector<char> green;
    for (std::size_t l = 0; l < links.size(); l++)
      if (is_crossing[l]) green.push_back(!signal_of[l] || open(s.t) || std::fmod(s.t, d.cycle) / d.cycle < d.green ? 1 : 0);
    for (double& v : rho) v /= width;
    out.window_times.push_back(s.t);
    out.window_density.push_back(std::move(rho));
    out.window_flow.push_back(std::move(flow));
    out.window_speed.push_back(std::move(speed));
    out.window_width.push_back(width);
    out.window_green.push_back(std::move(green));
  };
  const double h = d.step > 0 ? d.step : 0.9 * net.cfl();
  auto state = net.start(d.from * 3600);
  pc::KinematicNetwork::Budget budget;
  double next_snapshot = state.t;
  while (state.t < d.to * 3600) {
    if (state.t >= next_snapshot - 1e-9) {
      std::vector<double> rho;
      const double hour = state.t / 3600;
      const bool in_rule_hours = hour >= d.rule_from && hour < d.rule_to;
      for (std::size_t l = 0; l < links.size(); l++) {
        const double width = open(state.t) ? width_rule : width_day;
        for (double v : state.density[l]) {
          rho.push_back(v / width);
          out.max_density = std::max(out.max_density, v / width);
          if (in_rule_hours) out.max_density_in_rule_hours = std::max(out.max_density_in_rule_hours, v / width);
        }
      }
      out.times.push_back(state.t);
      out.density.push_back(std::move(rho));
      out.entered.push_back(budget.entered);
      out.left.push_back(budget.left);
      next_snapshot += d.snapshot;
    }
    const bool in_window = windowed && next_window <= d.window_to + 1e-9;
    if (in_window && state.t >= next_window - 1e-9) {
      keep_window(state);
      next_window += d.window_step;
    }
    double hh = std::min(h, d.to * 3600 - state.t);
    if (windowed && next_window <= d.window_to + 1e-9 && state.t < next_window) hh = std::min(hh, next_window - state.t);
    state = net.step(state, hh, &budget);
  }
  if (windowed && next_window <= d.window_to + 1e-9 && state.t >= next_window - 1e-9) keep_window(state);
  out.budget_residual = budget.residual();
  out.throughput = budget.left;
  // Mean time on the street from the cumulative curves (first in, first out): the area between them over the
  // number who left, at the snapshots' resolution.
  double area = 0;
  for (std::size_t k = 1; k < out.times.size(); k++)
    area += 0.5 * ((out.entered[k] - out.left[k]) + (out.entered[k - 1] - out.left[k - 1])) * (out.times[k] - out.times[k - 1]);
  out.mean_time = out.throughput > 0 ? area / out.entered.back() : 0;
  // The longest queue: the longest run of cells at half the jam density or more, in any snapshot.
  for (const auto& rho : out.density) {
    double run = 0;
    for (std::size_t c = 0; c < rho.size(); c++) {
      const double dx = c + 1 < out.position.size() ? out.position[c + 1] - out.position[c] : d.cell;
      run = rho[c] >= 0.5 * d.jam ? run + dx : 0;
      out.queue_spill = std::max(out.queue_spill, run);
    }
  }
  return out;
}

}  // namespace urban::flow
