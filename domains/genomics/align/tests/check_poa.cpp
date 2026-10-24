// Checks of partial-order alignment: aligned to the chain of one sequence, the graph's optimum
// is Gotoh's (global and local, affine and linear gaps); aligned to graphs of several sequences,
// it is the best over every source-to-sink path of the graph (brute force), and the path it
// returns is a walk of the graph that costs what it claims; the graph's columns are a partial
// order and its multiple alignment holds every sequence, in increasing columns; the consensus of
// 20 noisy long reads of a 2 kb sequence recovers it; the induced alignment of a simulated
// family recovers most of its true pairs (the sum-of-pairs score of Thompson, Plewniak and Poch
// 1999, Nucleic Acids Research 27, 2682–2690).

#include <genomics/align/bitparallel.hpp>
#include <genomics/align/gotoh.hpp>
#include <genomics/align/poa.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/simulate/family.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::align;

namespace {

constexpr std::uint32_t NOISE = 0x504f4153;

Codes random_codes(std::size_t n, int alphabet, calx::random::Stream& r) {
  Codes x(n);
  for (auto& c : x) c = static_cast<Code>(r.below(static_cast<std::uint64_t>(alphabet)));
  return x;
}

// The score of a path of steps, its runs of I (a base against no node) and of D (a node against
// no base) each paying one opening: what the alignment claims, recomputed. Local paths are
// scored between their first and last pair.
std::int64_t path_score(const PoaGraph& g, const PoaGraph::Path& path, const Codes& x, const Scoring& s, bool local) {
  std::size_t a = 0, b = path.size();
  if (local) {
    while (a < b && !(path[a].node >= 0 && path[a].position >= 0)) a++;
    while (b > a && !(path[b - 1].node >= 0 && path[b - 1].position >= 0)) b--;
  }
  std::int64_t score = 0;
  char run = 0;
  for (std::size_t k = a; k < b; k++) {
    const auto& st = path[k];
    if (st.node >= 0 && st.position >= 0) {
      score += s.matrix(g.node(static_cast<std::uint32_t>(st.node)).base, x[static_cast<std::size_t>(st.position)]);
      run = 0;
    } else {
      const char kind = st.node < 0 ? 'I' : 'D';
      score -= kind == run ? s.gap_extend : s.gap_open + s.gap_extend;
      run = kind;
    }
  }
  return score;
}

// Whether the path's nodes are a walk of the graph (consecutive nodes joined by an edge), from a
// source to a sink when `whole`, and its positions cover the sequence in order.
bool walk(const PoaGraph& g, const PoaGraph::Path& path, std::size_t n, bool whole) {
  std::vector<std::int32_t> nodes;
  std::int32_t next = 0;
  for (const auto& st : path) {
    if (st.node >= 0) nodes.push_back(st.node);
    if (st.position >= 0 && st.position != next++) return false;
  }
  if (static_cast<std::size_t>(next) != n) return false;
  for (std::size_t k = 1; k < nodes.size(); k++) {
    bool joined = false;
    for (std::uint32_t e : g.node(static_cast<std::uint32_t>(nodes[k - 1])).out) joined = joined || g.edge(e).to == static_cast<std::uint32_t>(nodes[k]);
    if (!joined) return false;
  }
  if (whole && !nodes.empty())
    return g.node(static_cast<std::uint32_t>(nodes.front())).in.empty() && g.node(static_cast<std::uint32_t>(nodes.back())).out.empty();
  return true;
}

// Every source-to-sink path of a small graph, as its sequence of bases.
std::vector<Codes> all_paths(const PoaGraph& g) {
  std::vector<Codes> out;
  Codes cur;
  std::function<void(std::uint32_t)> go = [&](std::uint32_t u) {
    cur.push_back(g.node(u).base);
    if (g.node(u).out.empty()) out.push_back(cur);
    for (std::uint32_t e : g.node(u).out) go(g.edge(e).to);
    cur.pop_back();
  };
  for (std::uint32_t u = 0; u < g.size(); u++)
    if (g.node(u).in.empty()) go(u);
  return out;
}

Codes ungapped(const Codes& row) {
  Codes x;
  for (Code c : row)
    if (c != sequence::GAP) x.push_back(c);
  return x;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  Scoring linear;
  linear.matrix = Matrix::dna(1, -1, -1);
  linear.gap_open = 0;
  linear.gap_extend = 1;
  const std::vector<Scoring> scorings = {minimap2_map_ont(), linear, blastn(), bwa_mem()};

  // ---- a chain is Gotoh's edit graph -------------------------------------------------------------------
  {
    calx::random::Stream r(1, NOISE, 1, 0);
    int ok = 0, total = 0;
    std::string bad;
    for (int trial = 0; trial < 400; trial++) {
      const Scoring& sc = scorings[static_cast<std::size_t>(trial) % scorings.size()];
      const Codes y = random_codes(1 + r.below(trial < 300 ? 30 : 400), 4, r), x = random_codes(r.below(trial < 300 ? 30 : 400), 4, r);
      const PoaGraph g(y);
      for (Mode mode : {Mode::global, Mode::local}) {
        total++;
        const auto res = g.align(x, sc, mode);
        const std::int64_t want = gotoh_score(x, y, sc, mode);
        if (res.score == want && path_score(g, res.path, x, sc, mode == Mode::local) == want && walk(g, res.path, x.size(), false)) ok++;
        else if (bad.empty()) bad = std::string(mode_name(mode)) + " |x|=" + std::to_string(x.size()) + " |y|=" + std::to_string(y.size()) + ": " + std::to_string(res.score) + " vs " + std::to_string(want);
      }
    }
    s.require("aligned to the chain of one sequence, the optimum is Gotoh's, global and local, under 4 scorings (" + std::to_string(total) + "), each path costing what it claims",
              ok == total, bad);
  }

  // ---- graphs of several sequences against every path, by brute force -------------------------------------
  {
    calx::random::Stream r(2, NOISE, 2, 0);
    int ok = 0, total = 0, orders = 0;
    std::string bad;
    for (int trial = 0; trial < 300; trial++) {
      const Scoring& sc = scorings[static_cast<std::size_t>(trial) % scorings.size()];
      const int alphabet = trial % 2 ? 2 : 4;
      PoaGraph g;
      const std::size_t members = 2 + r.below(3);
      for (std::size_t k = 0; k < members; k++) g.add(random_codes(1 + r.below(7), alphabet, r), sc, k % 2 ? Mode::local : Mode::global);
      // The order of the graph is a linear extension: every edge goes forward, columns are runs.
      std::vector<std::size_t> rank(g.size());
      for (std::size_t i = 0; i < g.order().size(); i++) rank[g.order()[i]] = i;
      bool forward = true;
      for (std::uint32_t e = 0; e < g.edges(); e++) forward = forward && rank[g.edge(e).from] < rank[g.edge(e).to] && g.column_of()[g.edge(e).from] < g.column_of()[g.edge(e).to];
      orders += forward;
      const auto paths = all_paths(g);
      const Codes x = random_codes(r.below(9), alphabet, r);
      for (Mode mode : {Mode::global, Mode::local}) {
        std::int64_t want = std::numeric_limits<std::int64_t>::min();
        for (const auto& p : paths) want = std::max(want, gotoh_score(x, p, sc, mode));
        const auto res = g.align(x, sc, mode);
        total++;
        if (res.score == want && path_score(g, res.path, x, sc, mode == Mode::local) == want && walk(g, res.path, x.size(), mode == Mode::global)) ok++;
        else if (bad.empty()) bad = std::string(mode_name(mode)) + " trial " + std::to_string(trial) + ": " + std::to_string(res.score) + " vs " + std::to_string(want);
      }
    }
    s.require("aligned to graphs of 2 to 4 sequences, the optimum is the best over every source-to-sink path (" + std::to_string(total) + "), each path a walk of the graph",
              ok == total, bad);
    s.require("every edge of every graph goes forward in its order and between columns in order", orders == 300);
  }

  // ---- the consensus of noisy long reads -----------------------------------------------------------------
  {
    const Codes truth = simulate::random_sequence(2000, 0.5, 31);
    simulate::LongReadModel model;
    model.mean_length = 2000;
    model.min_length = 2000;  // every read the whole sequence
    model.error = 0.08;
    const auto reads = simulate::long_reads(simulate::Source::haploid({truth}), 20, model, 32);
    PoaGraph g;
    std::vector<Codes> oriented;
    double read_errors = 0;
    for (const auto& rd : reads) {
      const Codes x = rd.reverse ? sequence::reverse_complement(rd.bases) : rd.bases;
      oriented.push_back(x);
      read_errors += static_cast<double>(myers_distance(x, truth));
      g.add(x, minimap2_map_ont(), Mode::global);
    }
    const Codes cons = g.consensus();
    const std::int64_t d = myers_distance(cons, truth);
    s.require("the heaviest bundle of 20 reads of a 2 kb sequence at 8% error (" + std::to_string(read_errors / 20 / 20) + "% measured) is within 4 edits of it",
              reads.size() == 20 && d <= 4, std::to_string(d) + " edits, " + std::to_string(g.size()) + " nodes");
    // The multiple alignment the graph induces.
    const auto rows = g.msa(true);
    bool same_length = true, holds = true, increasing = true;
    for (std::size_t k = 0; k < rows.size(); k++) {
      same_length = same_length && rows[k].size() == g.columns();
      if (k < oriented.size()) holds = holds && ungapped(rows[k]) == oriented[k];
    }
    holds = holds && ungapped(rows.back()) == cons;
    for (std::size_t k = 0; k < g.sequences(); k++) {
      const auto& p = g.path_of(k);
      for (std::size_t i = 1; i < p.size(); i++) increasing = increasing && g.column_of()[p[i - 1]] < g.column_of()[p[i]];
    }
    // A column's nodes have distinct bases, and the rows that share a base in a column share its node.
    std::set<std::pair<std::uint32_t, Code>> cells;
    for (std::uint32_t u = 0; u < g.size(); u++) cells.insert({g.column_of()[u], g.node(u).base});
    s.require("its multiple alignment has a row for each read and the consensus, each the sequence with gaps, in " + std::to_string(g.columns()) + " columns",
              rows.size() == 21 && same_length && holds);
    s.require("each read passes its columns in increasing order, and a column's nodes hold distinct bases", increasing && cells.size() == g.size());
  }

  // ---- the induced alignment of a family against its truth --------------------------------------------------
  {
    const Codes root = simulate::random_sequence(600, 0.5, 41);
    const simulate::Family fam = simulate::evolve_family(root, 8, 0.02, 0.005, 42);
    PoaGraph g;
    for (const auto& x : fam.sequences) g.add(x, minimap2_map_ont(), Mode::global);
    const auto rows = g.msa();
    // The pairs (member a's base i, member b's base j) in one column, in the truth and in the graph.
    const auto pairs = [](const std::vector<Codes>& aln) {
      std::set<std::array<std::size_t, 4>> out;
      std::vector<std::size_t> pos(aln.size(), 0);
      const std::size_t cols = aln.empty() ? 0 : aln[0].size();
      for (std::size_t c = 0; c < cols; c++) {
        for (std::size_t a = 0; a < aln.size(); a++)
          for (std::size_t b = a + 1; b < aln.size(); b++)
            if (aln[a][c] != sequence::GAP && aln[b][c] != sequence::GAP) out.insert({a, pos[a], b, pos[b]});
        for (std::size_t a = 0; a < aln.size(); a++) pos[a] += aln[a][c] != sequence::GAP;
      }
      return out;
    };
    const auto truth = pairs(fam.alignment), found = pairs(rows);
    std::size_t shared = 0;
    for (const auto& p : truth) shared += found.count(p);
    const double sp = truth.empty() ? 0 : static_cast<double>(shared) / static_cast<double>(truth.size());
    s.check("the graph's alignment of a family of 8 (600 bp; per branch 2% substitutions, 0.5% indels) misses under 5% of its true pairs", 1 - sp, 0.05);
  }

  return s.done();
}
