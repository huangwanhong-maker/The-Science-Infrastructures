// The toolbox's games and choice: nomos's two-player games (equilibria, dominance, zero-sum games),
// matching (deferred acceptance) and social choice (scoring rules, runoff, Condorcet methods), on
// games, preferences and ballots typed in the page.

#pragma once

#include <nomos/choice/condorcet.hpp>
#include <nomos/choice/profile.hpp>
#include <nomos/choice/scoring.hpp>
#include <nomos/games/bimatrix.hpp>
#include <nomos/games/dominance.hpp>
#include <nomos/games/equilibria.hpp>
#include <nomos/games/zero_sum.hpp>
#include <nomos/matching/deferred_acceptance.hpp>
#include <nomos/matching/stability.hpp>

#include <numeric>

#include "answer.hpp"

namespace science::toolbox {

namespace ng = nomos::games;

inline std::string strategy(bool row, std::size_t i) { return (row ? "R" : "C") + std::to_string(i + 1); }

inline std::string mixed(const std::vector<double>& p, bool row) {
  std::string s;
  for (std::size_t i = 0; i < p.size(); i++)
    if (p[i] > 1e-9) s += (s.empty() ? "" : ", ") + strategy(row, i) + " " + fixed(p[i], 3);
  return s;
}

// A bimatrix game: one row a line, each cell a pair of payoffs (the row player's, the column player's).
inline Json game(const Request& r) {
  std::vector<std::vector<double>> a, b;
  for (const auto& line : lines_of(r.text("payoffs", ""))) {
    const auto v = numbers_in(line);
    if (v.empty()) continue;
    if (v.size() % 2) throw std::invalid_argument("each cell is a pair of payoffs: \"" + line + "\"");
    std::vector<double> ra, rb;
    for (std::size_t k = 0; k < v.size(); k += 2) {
      ra.push_back(v[k]);
      rb.push_back(v[k + 1]);
    }
    a.push_back(ra);
    b.push_back(rb);
  }
  if (a.empty() || a[0].empty()) throw std::invalid_argument("type the payoffs, one row a line");
  for (const auto& row : a)
    if (row.size() != a[0].size()) throw std::invalid_argument("every row needs as many cells as the first");
  const ng::Bimatrix g(calx::linalg::Rectangular::from_rows(a), calx::linalg::Rectangular::from_rows(b));
  Arr blocks;
  // The game itself, with its pure equilibria marked.
  const auto pure = ng::pure_equilibria(g);
  Arr rows;
  std::vector<std::string> head{""};
  for (std::size_t j = 0; j < g.cols(); j++) head.push_back(strategy(false, j));
  for (std::size_t i = 0; i < g.rows(); i++) {
    Arr row{Json(strategy(true, i))};
    for (std::size_t j = 0; j < g.cols(); j++) {
      const bool eq = std::any_of(pure.begin(), pure.end(), [&](const ng::Cell& c) { return c.row == i && c.col == j; });
      row.push_back(Json(fixed(g.a(i, j), 2) + " , " + fixed(g.b(i, j), 2) + (eq ? "  ★" : "")));
    }
    rows.push_back(Json(std::move(row)));
  }
  blocks.push_back(table("The game: the row player's payoff, the column player's (★ a pure equilibrium)", head, std::move(rows)));
  const auto all = ng::support_enumeration(g);
  Arr eq;
  for (const auto& e : all)
    eq.push_back(Json(Arr{Json(mixed(e.x, true)), Json(mixed(e.y, false)), num(e.row_value), num(e.column_value)}));
  blocks.push_back(table("Nash equilibria, pure and mixed", {"the row player", "the column player", "row's payoff", "column's payoff"}, std::move(eq),
                         all.empty() ? "None found: the game is degenerate, and support enumeration finds the equilibria of nondegenerate games." : ""));
  for (const auto kind : {ng::Dominance::pure, ng::Dominance::mixed}) {
    const auto red = ng::eliminate_dominated(g, kind);
    std::string removed, left;
    for (const auto& s : red.removed) removed += (removed.empty() ? "" : ", ") + strategy(s.row_player, s.strategy);
    for (std::size_t i : red.rows) left += (left.empty() ? "" : " ") + strategy(true, i);
    for (std::size_t j : red.cols) left += " " + strategy(false, j);
    blocks.push_back(words(kind == ng::Dominance::pure ? "Strict dominance by pure strategies" : "Strict dominance by mixed strategies",
                           {{"removed, in order", removed.empty() ? "none" : removed}, {"what survives", left}}));
  }
  blocks.push_back(note("Equilibria by support enumeration (von Stengel 2007); dominance by mixed strategies by linear programming on calx. nomos::games."));
  return answer("A two-player game", std::move(blocks));
}

inline Json zero_sum(const Request& r) {
  std::vector<std::vector<double>> a;
  for (const auto& line : lines_of(r.text("matrix", ""))) {
    const auto v = numbers_in(line);
    if (!v.empty()) a.push_back(v);
  }
  if (a.empty()) throw std::invalid_argument("type the row player's payoffs, one row a line");
  for (const auto& row : a)
    if (row.size() != a[0].size()) throw std::invalid_argument("every row needs as many entries as the first");
  const auto s = ng::solve_zero_sum(calx::linalg::Rectangular::from_rows(a));
  std::vector<std::string> rl, cl;
  for (std::size_t i = 0; i < s.x.size(); i++) rl.push_back(strategy(true, i));
  for (std::size_t j = 0; j < s.y.size(); j++) cl.push_back(strategy(false, j));
  std::vector<double> x = s.x, y = s.y;
  Arr blocks;
  blocks.push_back(figures("The value", {{"the value of the game", s.value}, {"the row player guarantees", s.lower_value}, {"the column player concedes at most", s.upper_value}}));
  blocks.push_back(bars("The row player's optimal strategy", rl, {{"probability", x}}));
  blocks.push_back(bars("The column player's optimal strategy", cl, {{"probability", y}}));
  blocks.push_back(note("max_x min_y xᵀAy = min_y max_x xᵀAy (von Neumann 1928), both by linear programming on calx's simplex. nomos::games."));
  return answer("A zero-sum game", std::move(blocks));
}

// Lines "name: first second third" as names and rankings over the other side's names.
struct Rankings {
  std::vector<std::string> names;
  std::vector<std::vector<std::string>> ranking;
};
inline Rankings rankings(const std::string& text) {
  Rankings R;
  for (const auto& line : lines_of(text)) {
    const auto colon = line.find(':');
    if (colon == std::string::npos) throw std::invalid_argument("each line: a name, a colon, a ranking: \"" + line + "\"");
    const auto who = words_in(line.substr(0, colon));
    if (who.size() != 1) throw std::invalid_argument("one name before the colon: \"" + line + "\"");
    R.names.push_back(who[0]);
    R.ranking.push_back(words_in(line.substr(colon + 1), ">;"));
  }
  return R;
}

inline Json matching(const Request& r) {
  const auto P = rankings(r.text("proposers", "")), Q = rankings(r.text("receivers", ""));
  if (P.names.empty() || Q.names.empty()) throw std::invalid_argument("both sides need their rankings");
  auto index = [](const std::vector<std::string>& names, const std::string& who, const std::string& side) {
    const auto it = std::find(names.begin(), names.end(), who);
    if (it == names.end()) throw std::invalid_argument(who + " is not among the " + side);
    return static_cast<std::size_t>(it - names.begin());
  };
  nomos::matching::Market m;
  for (const auto& rk : P.ranking) {
    std::vector<std::size_t> v;
    for (const auto& w : rk) v.push_back(index(Q.names, w, "receivers"));
    m.applicants.push_back(v);
  }
  for (const auto& rk : Q.ranking) {
    std::vector<std::size_t> v;
    for (const auto& w : rk) v.push_back(index(P.names, w, "proposers"));
    m.institutions.push_back(v);
  }
  using nomos::matching::Proposing;
  const auto a = nomos::matching::deferred_acceptance(m, Proposing::applicants);
  const auto b = nomos::matching::deferred_acceptance(m, Proposing::institutions);
  auto partner = [&](const nomos::matching::Assignment& x, std::size_t i) {
    return x.institution_of[i] == nomos::matching::unmatched ? std::string("—") : Q.names[x.institution_of[i]];
  };
  auto rank_of = [&](std::size_t i, const nomos::matching::Assignment& x) {
    const auto& v = m.applicants[i];
    const auto it = std::find(v.begin(), v.end(), x.institution_of[i]);
    return it == v.end() ? NaN : static_cast<double>(it - v.begin() + 1);
  };
  Arr rows;
  for (std::size_t i = 0; i < P.names.size(); i++)
    rows.push_back(Json(Arr{Json(P.names[i]), Json(partner(a, i)), num(rank_of(i, a)), Json(partner(b, i)), num(rank_of(i, b))}));
  Arr blocks;
  blocks.push_back(table("The stable matchings at the two ends of the lattice", {"proposer", "proposers propose", "its rank", "receivers propose", "its rank"}, std::move(rows),
                         "Every stable matching lies between these: each proposer likes the first best and the second least among all stable matchings (Gale and Shapley 1962; Knuth 1976)."));
  blocks.push_back(words("Checks", {{"blocking pairs, proposers proposing", std::to_string(nomos::matching::blocking_pairs(m, a).size())},
                                    {"blocking pairs, receivers proposing", std::to_string(nomos::matching::blocking_pairs(m, b).size())},
                                    {"one stable matching only", a.institution_of == b.institution_of ? "yes" : "no"}}));
  blocks.push_back(note("Deferred acceptance; proposing truthfully is a dominant strategy for the proposers (Dubins and Freedman 1981; Roth 1982). nomos::matching."));
  return answer("Stable matching", std::move(blocks));
}

// Ballots "8: A > B > C" as a profile over the candidates named.
inline Json voting(const Request& r) {
  std::vector<std::string> names;
  nomos::choice::Profile prof;
  for (const auto& line : lines_of(r.text("ballots", ""))) {
    const auto colon = line.find(':');
    const auto count = colon == std::string::npos ? std::vector<double>{1} : numbers_in(line.substr(0, colon));
    if (count.size() != 1 || !(count[0] > 0)) throw std::invalid_argument("each line: a count, a colon, a ranking: \"" + line + "\"");
    nomos::choice::Ballot b;
    b.weight = count[0];
    for (const auto& w : words_in(colon == std::string::npos ? line : line.substr(colon + 1), ">;")) {
      auto it = std::find(names.begin(), names.end(), w);
      if (it == names.end()) {
        names.push_back(w);
        it = names.end() - 1;
      }
      const auto c = static_cast<std::size_t>(it - names.begin());
      if (std::find(b.ranking.begin(), b.ranking.end(), c) != b.ranking.end()) throw std::invalid_argument(w + " twice on one ballot");
      b.ranking.push_back(c);
    }
    prof.ballots.push_back(b);
  }
  prof.candidates = names.size();
  if (prof.candidates < 2) throw std::invalid_argument("at least two candidates");
  namespace nc = nomos::choice;
  const auto d = nc::pairwise(prof);
  const auto plural = nc::plurality(prof), borda = nc::borda(prof), cope = nc::copeland(d);
  const auto runoff = nc::instant_runoff(prof);
  const auto sch = nc::schulze(d);
  auto top = [&](const std::vector<double>& s) {
    std::string out;
    for (std::size_t w : nc::winners(s)) out += (out.empty() ? "" : " and ") + names[w];
    return out;
  };
  const auto cw = nc::condorcet_winner(d);
  std::string schulze_w;
  for (std::size_t w : sch.winners) schulze_w += (schulze_w.empty() ? "" : " and ") + names[w];
  Arr blocks;
  blocks.push_back(words("The winners", {{"plurality", top(plural)}, {"Borda", top(borda)}, {"instant runoff", names[runoff.winner]},
                                         {"Copeland", top(cope)}, {"Schulze", schulze_w}, {"the Condorcet winner", cw ? names[*cw] : "none: a majority cycle"}}));
  std::vector<std::string> order;
  for (std::size_t c = 0; c < names.size(); c++) order.push_back(names[c]);
  std::vector<std::size_t> elim_round(names.size(), 0);
  for (std::size_t k = 0; k < runoff.eliminated.size(); k++) elim_round[runoff.eliminated[k]] = k + 1;
  Arr rows;
  for (std::size_t c = 0; c < names.size(); c++)
    rows.push_back(Json(Arr{Json(names[c]), num(plural[c]), num(borda[c]), num(cope[c]),
                            Json(c == runoff.winner ? std::string("wins")
                                 : elim_round[c] ? "out in round " + std::to_string(elim_round[c]) : std::string("in the last round"))}));
  blocks.push_back(table("The scores", {"candidate", "plurality", "Borda", "Copeland", "instant runoff"}, std::move(rows)));
  std::vector<std::vector<double>> dm(names.size(), std::vector<double>(names.size(), NaN));
  for (std::size_t a = 0; a < names.size(); a++)
    for (std::size_t b = 0; b < names.size(); b++)
      if (a != b) dm[a][b] = d(a, b);
  blocks.push_back(matrix("Pairwise majorities: voters ranking the row above the column", names, names, dm));
  std::string ranking;
  for (std::size_t c : sch.ranking) ranking += (ranking.empty() ? "" : " > ") + names[c];
  blocks.push_back(note("Schulze's ranking: " + ranking + ". Arrow (1951): no rule of ranked ballots satisfies every fairness condition at once, so the rules may disagree. nomos::choice."));
  return answer("Voting rules", std::move(blocks));
}

}  // namespace science::toolbox
