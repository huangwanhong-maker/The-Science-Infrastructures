// Partial-order alignment: a graph of the bases of many sequences, a sequence aligned to the
// graph, the graph grown by it, the consensus of the graph, and the multiple alignment it holds.
//
// The graph (Lee, Grasso and Sharlow 2002, Bioinformatics 18, 452–464). A partial-order
// graph is a DAG whose nodes are bases and whose edges join the consecutive bases of the
// sequences that pass through them, weighted by how many do. Nodes that an alignment put in one
// column, with different bases, are *aligned* to each other: a column is a set of nodes of
// distinct bases, closed under the relation (Lee's "aligned ring"). The columns, ordered by the
// edges between their members, are a partial order whose linear extensions are the columns of
// a multiple alignment; this module keeps one, by Kahn's algorithm on the columns, the smallest
// first, so that a column's nodes are adjacent in the order of the graph.
//
// The alignment of a sequence to the graph is the (max, +) sum over the paths of the edit
// graph of the DAG and the sequence: its vertices are (node, position), and a step is a pair
// (from a predecessor of the node, consuming a base of the sequence), a base of the sequence
// against no node (I), or a node against no base (D, along an edge of the DAG). With Gotoh's
// states (gotoh.hpp), for a node u with predecessors P(u) (the virtual start S when it has
// none) and position j:
//   I(u, j) = max(H(u, j − 1) − (o + e), I(u, j − 1) − e)
//   D(u, j) = max over p ∈ P(u) of max(H(p, j) − (o + e), D(p, j) − e)
//   H(u, j) = max(H(p, j − 1) + s(base u, x_j) over p ∈ P(u), I(u, j), D(u, j) [, 0 local])
// computed in the order of the graph, so that every predecessor's row is complete (Lee et al.
// 2002; this is the edit graph of the product of the DAG and the path of the sequence,
// itself a DAG). Global: the whole sequence against a path from a source (S's successors) to a
// sink; local: any part of each. A gap of length L costs o + L e (types.hpp), so open 0 is a
// linear gap. The traceback prefers, at H, a pair (predecessors in the order of their edges),
// then I, then D; in a gap, its opening before its extension.
//
// Adding the sequence (Lee et al. 2002): a base paired with a node of the same base is
// fused into it; paired with a node of another base, into the node of its base in that node's
// column, or else a new node joins the column; a base against no node is a new node. Edges join
// the sequence's consecutive nodes, adding 1 (its weight) to those that exist. In local mode the
// unaligned ends of the sequence enter as new nodes.
//
// The consensus by the heaviest bundle (Lee 2003, Bioinformatics 19, 999–1008): each node
// takes its heaviest incoming edge (ties to the predecessor of greater score), its score that
// edge's weight plus its predecessor's; the path back from the best node, continued forward
// along the heaviest outgoing edges to a sink, is the consensus. The multiple alignment: one row
// per sequence, its base in the column of each node it passes, sequence::GAP elsewhere.
//
// Complexity: aligning a sequence of length n to a graph of V nodes and E edges costs
// O((V + E) n) time and 9 bytes per (node, position) for the traceback (2 kb against 4000
// nodes: 72 MB); adding it O(n + V + E) with the new order. Scores fit in 32 bits.

#pragma once

#include <genomics/align/types.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace genomics::align {

class PoaGraph {
 public:
  struct Node {
    Code base = 0;
    std::vector<std::uint32_t> in, out;  // edges
    std::vector<std::uint32_t> aligned;  // the other nodes of its column
    std::vector<std::uint32_t> sequences;
  };
  struct Edge {
    std::uint32_t from = 0, to = 0;
    std::int64_t weight = 0;
    std::vector<std::uint32_t> sequences;
  };
  // One step of a sequence's alignment to the graph: a node and a position of the sequence, or
  // −1 for none (a base against no node, a node against no base).
  struct Step {
    std::int32_t node = -1, position = -1;
    friend bool operator==(const Step&, const Step&) = default;
  };
  using Path = std::vector<Step>;
  struct Result {
    std::int64_t score = 0;
    Path path;
  };

  PoaGraph() = default;
  // The graph of one sequence: a chain.
  explicit PoaGraph(const Codes& first, std::int64_t weight = 1) {
    Path p;
    for (std::size_t j = 0; j < first.size(); j++) p.push_back({-1, static_cast<std::int32_t>(j)});
    add(first, p, weight);
  }

  std::size_t size() const { return nodes_.size(); }
  std::size_t edges() const { return edges_.size(); }
  std::size_t sequences() const { return paths_.size(); }
  const Node& node(std::uint32_t u) const { return nodes_[u]; }
  const Edge& edge(std::uint32_t e) const { return edges_[e]; }
  // The nodes in an order of the graph, each column's nodes adjacent.
  const std::vector<std::uint32_t>& order() const { return order_; }
  // The nodes a sequence passes through, in order.
  const std::vector<std::uint32_t>& path_of(std::size_t sequence) const { return paths_[sequence]; }
  // The column of each node, numbered in the order of the graph.
  const std::vector<std::uint32_t>& column_of() const { return column_; }
  std::size_t columns() const { return columns_; }

  // The best alignment of a sequence to the graph, global or local, and its path.
  Result align(const Codes& x, const Scoring& s, Mode mode = Mode::global) const;
  // Grows the graph by a sequence along its path (from align); returns the sequence's index.
  std::size_t add(const Codes& x, const Path& path, std::int64_t weight = 1);
  // Aligns and adds; returns the alignment's score.
  std::int64_t add(const Codes& x, const Scoring& s, Mode mode = Mode::global, std::int64_t weight = 1) {
    if (nodes_.empty()) {
      PoaGraph g(x, weight);
      *this = std::move(g);
      return 0;
    }
    Result r = align(x, s, mode);
    add(x, r.path, weight);
    return r.score;
  }

  // The consensus by the heaviest bundle, and the nodes it passes.
  Codes consensus(std::vector<std::uint32_t>* nodes = nullptr) const;
  // The rows of the multiple alignment the graph holds: one per sequence (and the consensus
  // last, when asked), sequence::GAP where a row has no base.
  std::vector<Codes> msa(bool with_consensus = false) const;

 private:
  std::uint32_t new_node(Code base) {
    nodes_.push_back({});
    nodes_.back().base = base;
    return static_cast<std::uint32_t>(nodes_.size() - 1);
  }
  void link(std::uint32_t from, std::uint32_t to, std::int64_t weight, std::uint32_t sequence) {
    for (std::uint32_t e : nodes_[from].out)
      if (edges_[e].to == to) {
        edges_[e].weight += weight;
        edges_[e].sequences.push_back(sequence);
        return;
      }
    edges_.push_back({from, to, weight, {sequence}});
    const std::uint32_t e = static_cast<std::uint32_t>(edges_.size() - 1);
    nodes_[from].out.push_back(e);
    nodes_[to].in.push_back(e);
  }
  void sort();

  std::vector<Node> nodes_;
  std::vector<Edge> edges_;
  std::vector<std::vector<std::uint32_t>> paths_;
  std::vector<std::uint32_t> order_, rank_, column_;
  std::size_t columns_ = 0;
};

// Kahn's algorithm on the columns: a column is ready when no edge enters it from a column not
// yet placed; the ready column of the smallest node goes first, its nodes in increasing order.
inline void PoaGraph::sort() {
  const std::size_t V = nodes_.size();
  std::vector<std::uint32_t> head(V);  // the smallest node of each node's column
  for (std::uint32_t u = 0; u < V; u++) {
    head[u] = u;
    for (std::uint32_t a : nodes_[u].aligned) head[u] = std::min(head[u], a);
  }
  std::vector<std::uint32_t> waiting(V, 0);  // per column head: edges from other columns not yet placed
  for (const auto& e : edges_)
    if (head[e.from] != head[e.to]) waiting[head[e.to]]++;
    else throw std::logic_error("genomics::align::PoaGraph: an edge inside a column");
  std::priority_queue<std::uint32_t, std::vector<std::uint32_t>, std::greater<>> ready;
  for (std::uint32_t u = 0; u < V; u++)
    if (head[u] == u && waiting[u] == 0) ready.push(u);
  order_.clear();
  column_.assign(V, 0);
  columns_ = 0;
  std::vector<std::uint32_t> members;
  while (!ready.empty()) {
    const std::uint32_t c = ready.top();
    ready.pop();
    members.assign(1, c);
    members.insert(members.end(), nodes_[c].aligned.begin(), nodes_[c].aligned.end());
    std::sort(members.begin(), members.end());
    for (std::uint32_t u : members) {
      order_.push_back(u);
      column_[u] = static_cast<std::uint32_t>(columns_);
    }
    columns_++;
    for (std::uint32_t u : members)
      for (std::uint32_t e : nodes_[u].out)
        if (--waiting[head[edges_[e].to]] == 0) ready.push(head[edges_[e].to]);
  }
  if (order_.size() != V) throw std::logic_error("genomics::align::PoaGraph: the columns are not a partial order");
  rank_.assign(V, 0);
  for (std::uint32_t r = 0; r < V; r++) rank_[order_[r]] = r;
}

inline std::size_t PoaGraph::add(const Codes& x, const Path& path, std::int64_t weight) {
  const std::uint32_t seq = static_cast<std::uint32_t>(paths_.size());
  std::vector<std::uint32_t> visited;
  visited.reserve(x.size());
  std::int64_t prev = -1;
  std::size_t covered = 0;
  for (const Step& st : path) {
    if (st.position < 0) continue;
    if (static_cast<std::size_t>(st.position) != covered || covered >= x.size()) throw std::invalid_argument("genomics::align::PoaGraph::add: the path does not cover the sequence in order");
    covered++;
    const Code c = x[static_cast<std::size_t>(st.position)];
    std::uint32_t id;
    if (st.node < 0) {
      id = new_node(c);
    } else {
      const std::uint32_t u = static_cast<std::uint32_t>(st.node);
      if (u >= nodes_.size()) throw std::invalid_argument("genomics::align::PoaGraph::add: a node not in the graph");
      if (nodes_[u].base == c) id = u;
      else {
        id = std::numeric_limits<std::uint32_t>::max();
        for (std::uint32_t a : nodes_[u].aligned)
          if (nodes_[a].base == c) { id = a; break; }
        if (id == std::numeric_limits<std::uint32_t>::max()) {
          // A new member of u's column, aligned to every node of it.
          id = new_node(c);
          std::vector<std::uint32_t> ring = nodes_[u].aligned;
          ring.push_back(u);
          for (std::uint32_t a : ring) {
            nodes_[a].aligned.push_back(id);
            nodes_[id].aligned.push_back(a);
          }
        }
      }
    }
    nodes_[id].sequences.push_back(seq);
    if (prev >= 0) link(static_cast<std::uint32_t>(prev), id, weight, seq);
    visited.push_back(id);
    prev = id;
  }
  if (covered != x.size()) throw std::invalid_argument("genomics::align::PoaGraph::add: the path does not cover the sequence");
  paths_.push_back(std::move(visited));
  sort();
  return seq;
}

inline PoaGraph::Result PoaGraph::align(const Codes& x, const Scoring& s, Mode mode) const {
  if (mode != Mode::global && mode != Mode::local) throw std::invalid_argument("genomics::align::PoaGraph::align: global or local");
  if (s.gap_open < 0 || s.gap_extend < 0) throw std::invalid_argument("genomics::align::PoaGraph::align: negative gap costs");
  for (Code c : x)
    if (c >= s.matrix.size) throw std::invalid_argument("genomics::align::PoaGraph::align: a code outside the matrix");
  for (const auto& nd : nodes_)
    if (nd.base >= s.matrix.size) throw std::invalid_argument("genomics::align::PoaGraph::align: a node's code outside the matrix");
  const bool local = mode == Mode::local;
  const std::size_t V = nodes_.size(), n = x.size(), W = n + 1;
  const std::int32_t NEG = std::numeric_limits<std::int32_t>::min() / 4;
  const std::int32_t oe = s.gap_open + s.gap_extend, e = s.gap_extend;
  // Row 0 is the virtual start S; row r + 1 the node of rank r.
  std::vector<std::int32_t> H((V + 1) * W), D((V + 1) * W, NEG);
  std::vector<std::uint8_t> dir((V + 1) * W, 0);
  enum : std::uint8_t { FROM_START = 0, FROM_PAIR = 1, FROM_I = 2, FROM_D = 3, I_EXTENDS = 4 };
  // The start row: the sequence's prefix against no node (a gap), or free in local mode.
  H[0] = 0;
  for (std::size_t j = 1; j <= n; j++) H[j] = local ? 0 : -(s.gap_open + static_cast<std::int32_t>(j) * e);
  // The predecessors' rows of each node (0 for S).
  std::vector<std::uint32_t> pred_rows;
  std::vector<std::size_t> pred_begin(V + 1, 0);
  for (std::size_t r = 0; r < V; r++) {
    pred_begin[r] = pred_rows.size();
    const Node& nd = nodes_[order_[r]];
    if (nd.in.empty()) pred_rows.push_back(0);
    for (std::uint32_t ed : nd.in) pred_rows.push_back(rank_[edges_[ed].from] + 1);
  }
  pred_begin[V] = pred_rows.size();
  struct Best {
    std::int32_t score = std::numeric_limits<std::int32_t>::min();
    std::size_t row = 0, j = 0;
  } best;
  if (local) best = {0, 0, 0};
  for (std::size_t r = 0; r < V; r++) {
    const std::size_t row = r + 1;
    const Node& nd = nodes_[order_[r]];
    const int* prof = s.matrix.score.data() + static_cast<std::size_t>(nd.base) * static_cast<std::size_t>(s.matrix.size);
    std::int32_t* h = H.data() + row * W;
    std::int32_t* d = D.data() + row * W;
    std::uint8_t* t = dir.data() + row * W;
    // D first, for every position: it needs only the predecessors' rows.
    for (std::size_t k = pred_begin[r]; k < pred_begin[r + 1]; k++) {
      const std::int32_t* ph = H.data() + static_cast<std::size_t>(pred_rows[k]) * W;
      const std::int32_t* pd = D.data() + static_cast<std::size_t>(pred_rows[k]) * W;
      for (std::size_t j = 0; j <= n; j++) d[j] = std::max(d[j], std::max(ph[j] - oe, pd[j] - e));
    }
    // Position 0: the node against no base, or a free start.
    h[0] = local ? std::max(0, d[0]) : d[0];
    t[0] = (local && h[0] == 0 && d[0] <= 0) ? FROM_START : FROM_D;
    std::int32_t ins = NEG;  // I(u, j)
    for (std::size_t j = 1; j <= n; j++) {
      const std::int32_t io = h[j - 1] - oe, ix = ins - e;
      const bool iext = ix > io;
      ins = iext ? ix : io;
      std::int32_t pair = NEG;
      for (std::size_t k = pred_begin[r]; k < pred_begin[r + 1]; k++) pair = std::max(pair, H[static_cast<std::size_t>(pred_rows[k]) * W + j - 1]);
      pair += prof[x[j - 1]];
      std::int32_t v = local ? 0 : NEG;
      std::uint8_t src = FROM_START;
      if (pair > v) { v = pair; src = FROM_PAIR; }
      if (ins > v) { v = ins; src = FROM_I; }
      if (d[j] > v) { v = d[j]; src = FROM_D; }
      h[j] = v;
      t[j] = static_cast<std::uint8_t>(src | (iext ? I_EXTENDS : 0));
      if (local && v > best.score) best = {v, row, j};
    }
    if (!local && nd.out.empty() && h[n] > best.score) best = {h[n], row, n};
  }
  Result res;
  if (V == 0) {
    // No graph: the sequence against nothing.
    res.score = local ? 0 : (n ? -(s.gap_open + static_cast<std::int32_t>(n) * e) : 0);
    for (std::size_t j = 0; j < n && !local; j++) res.path.push_back({-1, static_cast<std::int32_t>(j)});
    return res;
  }
  res.score = best.score;
  // The traceback.
  Path rev;
  std::size_t row = best.row, j = best.j;
  enum { IN_H, IN_I, IN_D } state = IN_H;
  const auto node_of = [&](std::size_t rw) { return static_cast<std::int32_t>(order_[rw - 1]); };
  for (;;) {
    if (row == 0) {
      // The start row: in global mode the rest of the sequence is a gap before the graph.
      if (!local)
        for (std::size_t k = j; k > 0; k--) rev.push_back({-1, static_cast<std::int32_t>(k - 1)});
      break;
    }
    const std::size_t at = row * W + j;
    const std::size_t r = row - 1;
    if (state == IN_H) {
      const std::uint8_t src = j == 0 ? (dir[at] & 3) : (dir[at] & 3);
      if (src == FROM_START) break;
      if (src == FROM_PAIR) {
        const std::int32_t want = H[at] - s.matrix(nodes_[order_[r]].base, x[j - 1]);
        std::size_t prow = SIZE_MAX;
        for (std::size_t k = pred_begin[r]; k < pred_begin[r + 1]; k++)
          if (H[static_cast<std::size_t>(pred_rows[k]) * W + j - 1] == want) { prow = pred_rows[k]; break; }
        if (prow == SIZE_MAX) throw std::logic_error("genomics::align::PoaGraph::align: the traceback lost its pair");
        rev.push_back({node_of(row), static_cast<std::int32_t>(j - 1)});
        row = prow;
        j--;
        continue;
      }
      state = src == FROM_I ? IN_I : IN_D;
      continue;
    }
    if (state == IN_I) {
      rev.push_back({-1, static_cast<std::int32_t>(j - 1)});
      state = (dir[at] & I_EXTENDS) ? IN_I : IN_H;
      j--;
      continue;
    }
    // In D: the node against no base, from the predecessor that gives D(u, j).
    rev.push_back({node_of(row), -1});
    std::size_t prow = SIZE_MAX;
    bool opens = false;
    for (std::size_t k = pred_begin[r]; k < pred_begin[r + 1] && prow == SIZE_MAX; k++)
      if (H[static_cast<std::size_t>(pred_rows[k]) * W + j] - oe == D[at]) { prow = pred_rows[k]; opens = true; }
    for (std::size_t k = pred_begin[r]; k < pred_begin[r + 1] && prow == SIZE_MAX; k++)
      if (D[static_cast<std::size_t>(pred_rows[k]) * W + j] - e == D[at]) prow = pred_rows[k];
    if (prow == SIZE_MAX) throw std::logic_error("genomics::align::PoaGraph::align: the traceback lost its gap");
    row = prow;
    state = opens ? IN_H : IN_D;
    if (row == 0 && !opens) throw std::logic_error("genomics::align::PoaGraph::align: a gap extended from the start");
  }
  std::reverse(rev.begin(), rev.end());
  if (local) {
    // The unaligned ends of the sequence enter as bases against no node.
    std::int32_t first = static_cast<std::int32_t>(n), last = -1;
    for (const Step& st : rev)
      if (st.position >= 0) { first = std::min(first, st.position); last = std::max(last, st.position); }
    if (last < 0) { first = static_cast<std::int32_t>(n); last = static_cast<std::int32_t>(n) - 1; }
    for (std::int32_t k = 0; k < first; k++) res.path.push_back({-1, k});
    res.path.insert(res.path.end(), rev.begin(), rev.end());
    for (std::int32_t k = last + 1; k < static_cast<std::int32_t>(n); k++) res.path.push_back({-1, k});
  } else {
    res.path = std::move(rev);
  }
  return res;
}

inline Codes PoaGraph::consensus(std::vector<std::uint32_t>* nodes) const {
  const std::size_t V = nodes_.size();
  Codes out;
  if (nodes) nodes->clear();
  if (V == 0) return out;
  std::vector<std::int64_t> score(V, 0);
  std::vector<std::int64_t> pred(V, -1);
  std::int64_t best = -1;
  std::uint32_t at = order_[0];
  for (std::uint32_t u : order_) {
    std::int64_t w = -1;
    for (std::uint32_t e : nodes_[u].in) {
      const Edge& ed = edges_[e];
      if (ed.weight > w || (ed.weight == w && score[ed.from] > score[static_cast<std::size_t>(pred[u])])) {
        w = ed.weight;
        pred[u] = ed.from;
      }
    }
    score[u] = pred[u] >= 0 ? w + score[static_cast<std::size_t>(pred[u])] : 0;
    if (score[u] > best) { best = score[u]; at = u; }
  }
  std::vector<std::uint32_t> path;
  for (std::int64_t u = at; u >= 0; u = pred[static_cast<std::size_t>(u)]) path.push_back(static_cast<std::uint32_t>(u));
  std::reverse(path.begin(), path.end());
  // On to a sink along the heaviest outgoing edges.
  for (std::uint32_t u = at; !nodes_[u].out.empty();) {
    std::uint32_t next = edges_[nodes_[u].out[0]].to;
    std::int64_t w = edges_[nodes_[u].out[0]].weight;
    for (std::uint32_t e : nodes_[u].out)
      if (edges_[e].weight > w) { w = edges_[e].weight; next = edges_[e].to; }
    path.push_back(next);
    u = next;
  }
  for (std::uint32_t u : path) out.push_back(nodes_[u].base);
  if (nodes) *nodes = std::move(path);
  return out;
}

inline std::vector<Codes> PoaGraph::msa(bool with_consensus) const {
  std::vector<Codes> rows(paths_.size() + (with_consensus ? 1 : 0), Codes(columns_, sequence::GAP));
  for (std::size_t k = 0; k < paths_.size(); k++)
    for (std::uint32_t u : paths_[k]) rows[k][column_[u]] = nodes_[u].base;
  if (with_consensus) {
    std::vector<std::uint32_t> cons;
    consensus(&cons);
    for (std::uint32_t u : cons) rows.back()[column_[u]] = nodes_[u].base;
  }
  return rows;
}

}  // namespace genomics::align
