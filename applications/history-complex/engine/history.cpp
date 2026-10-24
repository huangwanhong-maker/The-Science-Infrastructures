// The engine of the History Complex: one request in JSON on standard input, one answer in JSON
// on standard output (the application's README, §3). It composes the history domain (the
// record, the construction, the readings, the spin foam and the other complexes) with morphe's
// invariants (homology with generators, the fundamental group, persistence) and physis's (the ℤ₂
// gauge theory of the foam, exact and sampled; the recoupling of SU(2)). The canons are read
// from the paths the request names, and nothing else is read. Errors are answered as
// {"error": "..."} with exit status 1.
//
//   {"op": "build", "canon": path}                  the record, the complexes, their invariants, the readings, the foam
//   {"op": "build", "canons": [paths]}              the same for the canons joined into one record
//   {"op": "moment", "canon": path, "date": edtf}   the network of the moment at the date's last day, with its posterior
//   moment, power and field take "couplings": {name: β}, those a build fitted, and do not fit again
//   {"op": "power", "canon": path}                  the power of each being over the reconstruction
//   {"op": "field", "canon": path}                  the history as a Feynman graph and as a sample path

#include <history/construction/complex.hpp>
#include <history/construction/loops.hpp>
#include <history/readings/colours.hpp>
#include <history/readings/signed.hpp>
#include <history/foam/model.hpp>
#include <history/foam/models.hpp>
#include <history/foam/refined.hpp>
#include <history/foam/transfer.hpp>
#include <history/field/process.hpp>
#include <history/field/stranded.hpp>
#include <morphe/topology/invariants.hpp>
#include <morphe/topology/persistence.hpp>
#include <physis/symmetry/recoupling.hpp>

#include <cmath>
#include <exception>
#include <iostream>
#include <iterator>
#include <map>
#include <numbers>
#include <tuple>
#include <utility>
#include <optional>
#include <string>
#include <vector>

namespace {

using history::Json;
using Obj = Json::Object;
using Arr = Json::Array;

Json index(std::size_t i) { return Json(static_cast<double>(i)); }
Json index_or_null(const std::optional<std::size_t>& i) { return i ? index(*i) : Json(nullptr); }
Json text_or_null(const std::string& s) { return s.empty() ? Json(nullptr) : Json(s); }

template <class T, class F>
Json array_of(const std::vector<T>& v, F f) {
  Arr a;
  for (const auto& x : v) a.push_back(f(x));
  return Json(std::move(a));
}

Json indices(const std::vector<std::size_t>& v) { return array_of(v, [](std::size_t i) { return index(i); }); }
Json texts(const std::vector<std::string>& v) { return array_of(v, [](const std::string& s) { return Json(s); }); }

// The Complex Atlas's JSON of the complex: time along x, the beings round a cylinder.
Json atlas_json(const history::HistoryComplex& K) {
  const auto& r = K.record;
  const auto& X = K.cw;
  Arr edges, faces, positions, labels;
  for (const auto& e : X.edges) edges.push_back(Json(Arr{index(e[0]), index(e[1])}));
  for (const auto& f : X.faces) {
    Arr word;
    for (const auto& s : f.boundary) word.push_back(Json(static_cast<double>(s.sign) * static_cast<double>(s.edge + 1)));
    faces.push_back(Json(Obj{{"base", index(f.base)}, {"boundary", Json(std::move(word))}}));
  }
  double first = 0, last = 1;
  if (!K.order.empty()) {
    first = history::historical_year(r.occasions[K.order.front()].date.first);
    last = history::historical_year(r.occasions[K.order.back()].date.first);
  }
  const double span = last > first ? last - first : 1;
  const double n = static_cast<double>(r.beings.size() ? r.beings.size() : 1);
  for (std::size_t o : K.order) {
    const auto& occ = r.occasions[o];
    double y = 0, z = 0;
    for (const auto& p : occ.present) {
      const double angle = 2 * std::numbers::pi * static_cast<double>(p.being) / n;
      y += std::cos(angle);
      z += std::sin(angle);
    }
    const double k = occ.present.empty() ? 1 : static_cast<double>(occ.present.size());
    positions.push_back(Json(Arr{Json(16 * (history::historical_year(occ.date.first) - first) / span - 8), Json(4 * y / k), Json(4 * z / k)}));
  }
  for (const auto& s : K.segments) labels.push_back(Json(r.beings[s.being].name));
  return Json(Obj{{"name", Json(X.name)},
                  {"title", Json(r.canon.title.empty() ? X.name : r.canon.title)},
                  {"vertices", index(X.vertices)},
                  {"edges", Json(std::move(edges))},
                  {"faces", Json(std::move(faces))},
                  {"positions", Json(std::move(positions))},
                  {"labels", Json(Obj{{"edges", Json(std::move(labels))}})}});
}

Json moment_json(const history::Moment& M) {
  return Json(Obj{
      {"jdn", Json(static_cast<double>(M.jdn))},
      {"nodes", indices(M.nodes)},
      {"links", array_of(M.links, [](const history::Link& l) {
         return Json(Obj{{"tie", index(l.tie)}, {"a", index(l.a)}, {"b", index(l.b)}, {"sign", Json(l.sign)}, {"expected", Json(l.expected)},
                         {"since", index(l.since)}, {"lenses_since", index(l.lenses_since)}});
       })},
      {"triangles", array_of(M.triangles, [](const history::MomentTriangle& t) {
         return Json(Obj{{"beings", Json(Arr{index(t.beings[0]), index(t.beings[1]), index(t.beings[2])})},
                         {"ties", Json(Arr{index(t.ties[0]), index(t.ties[1]), index(t.ties[2])})},
                         {"product", Json(t.product)},
                         {"expected", Json(t.expected)}});
       })},
      {"frustrated", index(M.frustrated)},
      {"expected", Json(M.expected)}});
}


// ---- the spin foam at G = ℤ₂ (theory note §6) ----------------------------------------------------

// The couplings in use are the maximum a posteriori of the exact likelihood (transfer.hpp) where
// the contraction reaches, else the KT estimates; and the ensembles are exact where it reaches,
// else sampled. The answer says which.
struct Foam {
  history::RefinedComplex R;
  history::Couplings C;          // the KT estimates, cell by cell independent
  history::Couplings used;       // the couplings of M
  history::HistoryModel M;
  history::FoamSpace space;
  std::optional<history::Fit> fit;
  std::string inexact;           // why the contraction is not exact, or empty
};

using Given = std::optional<std::map<std::string, double>>;

// The couplings a request carries (those a build fitted), if any.
Given couplings_in(const Json& request) {
  const Json* c = request.find("couplings");
  if (!c || !c->is_object()) return std::nullopt;
  std::map<std::string, double> out;
  for (const auto& [k, v] : c->as_object())
    if (v.is_number()) out[k] = v.as_number();
  return out;
}

// The foam with the couplings given, or fitted here: exactly where the contraction reaches (a
// first pass for log Z alone says whether it does), else the KT estimates.
Foam foam_of(const history::HistoryComplex& K, const Given& given = std::nullopt) {
  Foam F;
  F.R = history::refine(K);
  F.C = history::estimate_couplings(K, F.R);
  F.used = F.C;
  F.M = history::history_model(K, F.R, F.C);
  F.space = history::foam_space(F.R);
  const auto P = history::parameters_of(K, F.M);
  try {
    std::size_t largest = 0;
    (void)history::log_z(K, F.M, false, largest);
    (void)history::log_z(K, F.M, true, largest);
  } catch (const history::Intractable& e) {
    F.inexact = e.what();
  }
  if (given) {
    for (std::size_t c = 0; c < F.M.cells.size(); c++) {
      const auto it = given->find(P.names[P.of_cell[c]]);
      if (it != given->end()) F.M.cells[c].beta = it->second;
    }
    for (const auto& [name, b] : *given) {
      if (name == "balance") F.used.balance = b;
      else F.used.lens[name] = b;
    }
    return F;
  }
  if (!F.inexact.empty()) return F;
  try {
    F.fit = history::fit_couplings(K, F.M, F.C);
    for (std::size_t c = 0; c < F.M.cells.size(); c++) F.M.cells[c].beta = F.fit->beta[P.of_cell[c]];
    for (std::size_t p = 0; p < F.fit->names.size(); p++) {
      if (F.fit->names[p] == "balance") F.used.balance = F.fit->beta[p];
      else F.used.lens[F.fit->names[p]] = F.fit->beta[p];
    }
  } catch (const history::Intractable& e) {
    F.inexact = e.what();
  }
  return F;
}

// The cut at the end of a day: after the last occasion that begins on or before it.
std::size_t cut_at(const history::HistoryComplex& K, long jdn) {
  std::size_t v = 0;
  for (std::size_t k = 0; k < K.order.size(); k++)
    if (K.record.occasions[K.order[k]].date.first <= jdn) v = k;
  return v;
}

// The mean and the 5% and 95% points of a law of counts (the mean alone when there is no law).
std::array<double, 3> law_summary(const std::vector<double>& law, double mean) {
  if (law.empty()) return {mean, mean, mean};
  double acc = 0, lo = -1, hi = -1;
  for (std::size_t k = 0; k < law.size(); k++) {
    acc += law[k];
    if (lo < 0 && acc >= 0.05 - 1e-12) lo = static_cast<double>(k);
    if (hi < 0 && acc >= 0.95 - 1e-12) hi = static_cast<double>(k);
  }
  return {mean, lo, hi};
}

// The triangles of the network of the moment at a day, as triples of rungs (the last rung of each
// live tie at or before the day).
history::Counter moment_counter(const history::HistoryComplex& K, const history::HistoryModel& M, const history::Moment& moment) {
  history::Counter counter;
  const auto& occ = K.record.occasions;
  auto last_rung = [&](std::size_t t) {
    std::size_t r = M.index[t].front();
    for (std::size_t k = 0; k < K.rungs[t].size(); k++)
      if (occ[K.rungs[t][k].occasion].date.first <= moment.jdn) r = M.index[t][k];
    return r;
  };
  for (const auto& T : moment.triangles) counter.triangles.push_back({last_rung(T.ties[0]), last_rung(T.ties[1]), last_rung(T.ties[2])});
  return counter;
}

// The mean and the 5% and 95% points of a histogram of counts.
std::array<double, 3> summary(const std::vector<std::size_t>& hist) {
  double n = 0, mean = 0;
  for (std::size_t k = 0; k < hist.size(); k++) { n += static_cast<double>(hist[k]); mean += static_cast<double>(k * hist[k]); }
  if (n == 0) return {0, 0, 0};
  mean /= n;
  double acc = 0, lo = -1, hi = -1;
  for (std::size_t k = 0; k < hist.size(); k++) {
    acc += static_cast<double>(hist[k]);
    if (lo < 0 && acc >= 0.05 * n) lo = static_cast<double>(k);
    if (hi < 0 && acc >= 0.95 * n) hi = static_cast<double>(k);
  }
  return {mean, lo, hi};
}

Json foam_json(const history::HistoryComplex& K, const Foam& F) {
  const auto& M = F.M;
  const history::Memory m = history::memory(history::persistence(K));
  // The posterior, with the frustration of the moment at the end of every year; the prior, its
  // spins and the Wilson loops of the ties' lives.
  const auto years = history::moments_by_year(K, m);
  std::vector<history::Counter> counters;
  for (const auto& y : years) counters.push_back(moment_counter(K, M, history::moment(K, y.jdn, m)));
  std::vector<std::vector<std::size_t>> lives;
  std::vector<std::size_t> lived;
  for (std::size_t t = 0; t < K.rungs.size(); t++)
    if (K.rungs[t].size() >= 2) { lives.push_back({M.index[t].front(), M.index[t].back()}); lived.push_back(t); }
  // Exactly where the contraction reaches; else by the heat bath.
  std::vector<double> post_rung, post_cell, prior_cell, loop, loop_error;
  std::vector<std::array<double, 3>> year_law;
  std::size_t width = 0, largest = 0;
  double log_z = 0;
  bool exact = F.inexact.empty();
  if (exact) {
    try {
      std::vector<history::Cut> cuts;
      for (std::size_t y = 0; y < years.size(); y++) {
        history::Cut c;
        c.after = cut_at(K, years[y].jdn);
        for (const auto& t : counters[y].triangles) c.keep.insert(c.keep.end(), t.begin(), t.end());
        cuts.push_back(std::move(c));
      }
      std::vector<std::array<std::size_t, 2>> pairs;
      for (const auto& l : lives) pairs.push_back({l[0], l[1]});
      const auto post = history::contract(K, M, true, {}, cuts, counters);
      const auto prior = history::contract(K, M, false, pairs);
      post_rung = post.link;
      post_cell = post.plaquette;
      prior_cell = prior.plaquette;
      loop = prior.two_point;
      loop_error.assign(loop.size(), 0);
      for (std::size_t y = 0; y < years.size(); y++) year_law.push_back(law_summary(post.cuts[y].negative, post.cuts[y].mean_negative));
      width = std::max(post.width, prior.width);
      largest = std::max(post.largest, prior.largest);
      log_z = prior.log_z;
    } catch (const history::Intractable&) {
      exact = false;
    }
  }
  history::Sampling S;
  if (!exact) {
    S.burn = 1000;
    S.sweeps = 20000;
    history::Observables at_years, around_lives;
    at_years.counters = counters;
    for (const std::size_t t : lived) around_lives.loops.push_back(history::tie_loop(K, F.R, t, 0, K.rungs[t].size() - 1));
    const history::Ensemble post = history::sample(F.R, M, true, S, at_years);
    S.seed += 1;
    const history::Ensemble prior = history::sample(F.R, M, false, S, around_lives);
    post_rung = post.rung;
    post_cell = post.cell;
    prior_cell = prior.cell;
    loop = prior.loop;
    loop_error = prior.loop_error;
    for (std::size_t y = 0; y < years.size(); y++) year_law.push_back(summary(post.count_hist[y]));
  }
  const auto j = history::spins(M, prior_cell, &F.space);

  // The couplings: the MAP of each parameter beside its KT estimate, or the KT estimates alone.
  Arr couplings;
  auto kt_of = [&](const std::string& of) -> const history::CouplingEstimate* {
    for (const auto& e : F.C.estimates)
      if (e.of == of) return &e;
    return nullptr;
  };
  if (F.fit) {
    const auto& fit = *F.fit;
    for (std::size_t p = 0; p < fit.names.size(); p++) {
      const auto* e = kt_of(fit.names[p]);
      const double b = fit.beta[p], se = fit.error[p];
      const bool known = std::isfinite(se);
      couplings.push_back(Json(Obj{{"of", Json(fit.names[p])}, {"cells", index(fit.cells[p])},
                                   {"n", e ? index(e->n) : index(0)}, {"k", e ? index(e->k) : index(0)},
                                   {"t", Json(std::tanh(b))}, {"beta", Json(b)}, {"error", known ? Json(se) : Json(nullptr)},
                                   {"low", known ? Json(std::tanh(b - 1.96 * se)) : Json(nullptr)},
                                   {"high", known ? Json(std::tanh(b + 1.96 * se)) : Json(nullptr)},
                                   {"kt_t", Json(std::tanh(fit.start[p]))}, {"kt_beta", Json(fit.start[p])}}));
    }
  } else {
    for (const auto& e : F.C.estimates)
      couplings.push_back(Json(Obj{{"of", Json(e.of)}, {"n", index(e.n)}, {"k", index(e.k)}, {"p", Json(e.p)}, {"t", Json(e.t)}, {"beta", Json(e.beta)},
                                   {"kt_t", Json(e.t)}, {"kt_beta", Json(e.beta)}}));
  }
  Arr cells;
  const std::size_t lenses = K.lenses.size();
  for (std::size_t c = 0; c < M.cells.size(); c++) {
    const double t = std::tanh(M.cells[c].beta);
    cells.push_back(Json(Obj{{"kind", Json(M.cells[c].triangle ? "triangle" : "lens")}, {"index", index(M.cells[c].triangle ? c - lenses : c)},
                             {"beta", Json(M.cells[c].beta)}, {"t", Json(t)}, {"prior", Json(prior_cell[c])}, {"spin", Json(j[c])},
                             {"support", Json(static_cast<bool>(F.space.support[c]))}, {"posterior", Json(post_cell[c])}}));
  }
  Arr triangles;
  for (const auto& T : F.R.triangles)
    triangles.push_back(Json(Obj{{"occasion", index(T.occasion)}, {"beings", Json(Arr{index(T.beings[0]), index(T.beings[1]), index(T.beings[2])})},
                                 {"ties", Json(Arr{index(T.ties[0]), index(T.ties[1]), index(T.ties[2])})}}));
  Arr ties;
  for (std::size_t t = 0; t < K.rungs.size(); t++) {
    Arr rungs;
    for (std::size_t k = 0; k < K.rungs[t].size(); k++) rungs.push_back(Json(post_rung[M.index[t][k]]));
    ties.push_back(Json(Obj{{"tie", index(t)}, {"posterior", Json(std::move(rungs))}}));
  }
  Arr loops;
  for (std::size_t l = 0; l < lived.size(); l++) {
    const std::size_t t = lived[l];
    const double tc = std::tanh(F.used.of(K.record.ties[t].colour));
    loops.push_back(Json(Obj{{"tie", index(t)}, {"lenses", index(K.rungs[t].size() - 1)}, {"wilson", Json(loop[l])}, {"error", Json(loop_error[l])},
                             {"chain", Json(std::pow(tc, static_cast<double>(K.rungs[t].size() - 1)))}}));
  }
  Arr moments;
  for (std::size_t y = 0; y < years.size(); y++) {
    const auto [mean, lo, hi] = year_law[y];
    moments.push_back(Json(Obj{{"year", Json(static_cast<double>(years[y].year))}, {"triangles", index(counters[y].triangles.size())},
                               {"mean", Json(mean)}, {"lo", Json(lo)}, {"hi", Json(hi)}}));
  }
  // SU(2) kinematics: at spin ½ on every live tie, the channels of each being-segment.
  std::vector<std::size_t> live(K.segments.size(), 0);
  for (const auto& f : K.cw.faces)
    for (const auto& step : f.boundary) live[step.edge]++;
  Arr channels;
  std::size_t admissible = 0, carrying = 0;
  for (std::size_t e = 0; e < live.size(); e++) {
    const double d = physis::symmetry::invariant_dimension(std::vector<int>(live[e], 1)).to_double();
    channels.push_back(Json(Obj{{"live", index(live[e])}, {"channels", Json(d)}}));
    if (live[e]) { carrying++; admissible += d > 0; }
  }
  return Json(Obj{{"couplings", Json(std::move(couplings))},
                  {"foams", index(F.space.dimension)},
                  {"cells", Json(std::move(cells))},
                  {"triangles", Json(std::move(triangles))},
                  {"ties", Json(std::move(ties))},
                  {"lives", Json(std::move(loops))},
                  {"moments", Json(std::move(moments))},
                  {"sweeps", exact ? Json(nullptr) : index(S.sweeps)},
                  {"exact", Json(exact)},
                  {"inexact", exact ? Json(nullptr) : Json(F.inexact.empty() ? std::string("the contraction is beyond reach") : F.inexact)},
                  {"width", index(width)}, {"largest", index(largest)},
                  {"log_z", exact ? Json(log_z) : Json(nullptr)},
                  {"log_likelihood", F.fit ? Json(F.fit->log_likelihood) : Json(nullptr)},
                  {"su2", Json(Obj{{"segments", Json(std::move(channels))}, {"carrying", index(carrying)}, {"admissible", index(admissible)}})}});
}

// ---- the other complexes (theory note §7) ----------------------------------------------------------

Json homology_json(const std::vector<morphe::algebra::HomologyGroup>& H) {
  return array_of(H, [](const morphe::algebra::HomologyGroup& h) {
    return Json(Obj{{"dimension", index(h.dimension)}, {"rank", index(h.rank)}, {"text", Json(h.str())}});
  });
}

Json counts_json(const morphe::topology::CWComplex& X) {
  Arr c;
  for (std::size_t n : X.counts()) c.push_back(index(n));
  return Json(std::move(c));
}

Json models_json(const history::HistoryComplex& K, const Foam& F, const std::vector<history::Shell>& shells) {
  namespace mt = morphe::topology;
  // K̂.
  Json refined = Json(Obj{{"counts", counts_json(F.R.cw)}, {"euler", Json(static_cast<double>(F.R.cw.euler_characteristic()))},
                          {"homology", homology_json(mt::homology(F.R.cw, false))}, {"triangles", index(F.R.triangles.size())},
                          {"foams", index(F.space.dimension)}});
  // K³.
  const auto K3 = history::shells_filled(K, shells);
  Json filled = Json(Obj{{"counts", counts_json(K3)}, {"euler", Json(static_cast<double>(K3.euler_characteristic()))}, {"homology", homology_json(mt::homology(K3, false))}});
  // The Dowker complex on the beings, and Atkin's Q-analysis.
  const auto D = history::dowker_beings(K);
  Arr simplices;
  for (std::size_t d = 0; d <= D.dimension(); d++) simplices.push_back(index(D.simplices[d].size()));
  Arr q;
  for (std::size_t n : history::q_analysis(K)) q.push_back(index(n));
  Json dowker = Json(Obj{{"simplices", Json(std::move(simplices))}, {"homology", homology_json(morphe::algebra::homology(D.chain_complex(), false))}, {"q", Json(std::move(q))}});
  // The history in time.
  const auto P = history::persistence_in_time(K);
  auto cell_json = [&](std::size_t dim, std::size_t i) -> Json {
    if (dim == 0) return Json(Obj{{"kind", Json("occasion")}, {"occasion", index(K.order[i])}});
    if (dim == 1) {
      const auto& s = K.segments[i];
      return Json(Obj{{"kind", Json("segment")}, {"being", index(s.being)}, {"from", index(s.from)}, {"to", index(s.to)}});
    }
    const auto& L = K.lenses[i];
    return Json(Obj{{"kind", Json("lens")}, {"tie", index(L.tie)}, {"from", index(L.from)}, {"to", index(L.to)}});
  };
  Arr intervals;
  for (const auto& I : P.intervals) {
    if (I.dimension == 0 && I.essential()) continue;  // the history itself, born with its first occasion
    intervals.push_back(Json(Obj{{"dimension", index(I.dimension)}, {"birth", index(K.order[I.birth])},
                                 {"death", I.death ? index(K.order[*I.death]) : Json(nullptr)},
                                 {"creator", cell_json(I.dimension, I.creator)},
                                 {"killer", I.killer ? cell_json(I.dimension + 1, *I.killer) : Json(nullptr)}}));
  }
  Arr betti;
  for (std::size_t stage = 0; stage < K.order.size(); stage++)
    betti.push_back(Json(Arr{index(P.betti(0, stage)), index(P.betti(1, stage)), index(P.betti(2, stage))}));
  Json time = Json(Obj{{"intervals", Json(std::move(intervals))}, {"betti", Json(std::move(betti))}});
  return Json(Obj{{"refined", refined}, {"filled", filled}, {"dowker", dowker}, {"time", time}});
}

Json build(const history::Record& record) {
  if (record.canon.id.empty() && !record.ok()) throw std::invalid_argument(record.problems.front().where + ": " + record.problems.front().message);
  const history::HistoryComplex K = history::construct(record);
  const history::Record& r = K.record;  // with the bereaved derived
  const auto P = history::persistence(K);
  const history::Memory m = history::memory(P);
  const history::Colours C = history::colours_of(K);
  Json answer = Json::object();
  const auto& c = r.canon;
  answer.set("canon", Json(Obj{{"id", Json(c.id)}, {"title", Json(c.title)}, {"about", Json(c.about)}, {"codebook", Json(c.codebook)},
                              {"period", Json(c.period)}, {"status", Json(c.status)}, {"drafted_by", Json(c.drafted_by)}, {"drafted", Json(c.drafted)}}));
  answer.set("sources", array_of(r.sources, [](const history::Source& s) { return Json(Obj{{"id", Json(s.id)}, {"cite", Json(s.cite)}}); }));

  Arr beings;
  for (std::size_t b = 0; b < r.beings.size(); b++) {
    const auto& x = r.beings[b];
    beings.push_back(Json(Obj{{"id", Json(x.id)}, {"name", Json(x.name)}, {"type", Json(x.type)}, {"lived", text_or_null(x.lived)},
                              {"wikidata", text_or_null(x.wikidata)}, {"colour", text_or_null(C.of_being[b])}, {"occasions", indices(K.worldline[b])},
                              {"drafted_by", text_or_null(x.drafted_by)}, {"checked_by", text_or_null(x.checked_by)}, {"note", text_or_null(x.note)}}));
  }
  answer.set("beings", Json(std::move(beings)));

  Arr occasions;
  for (const auto& o : r.occasions) {
    Arr present;
    for (const auto& p : o.present)
      present.push_back(Json(Obj{{"being", index(p.being)}, {"role", Json(p.role)}, {"derived", Json(p.derived)}}));
    occasions.push_back(Json(Obj{{"id", Json(o.id)}, {"date", Json(o.date.text)}, {"calendar", Json(history::calendar_name(o.date.calendar))},
                                 {"place", Json(o.place)}, {"type", Json(o.type)}, {"label", Json(o.label)},
                                 {"first", Json(static_cast<double>(o.date.first))}, {"last", Json(static_cast<double>(o.date.last))},
                                 {"year", Json(history::historical_year(o.date.first))}, {"present", Json(std::move(present))},
                                 {"sources", texts(o.sources)}, {"confidence", Json(o.confidence)}, {"drafted_by", text_or_null(o.drafted_by)},
                                 {"checked_by", text_or_null(o.checked_by)}, {"note", text_or_null(o.note)}}));
  }
  answer.set("occasions", Json(std::move(occasions)));

  Arr ties;
  for (std::size_t t = 0; t < r.ties.size(); t++) {
    const auto& x = r.ties[t];
    Arr attested, rungs;
    for (const auto& a : x.attested) attested.push_back(Json(Obj{{"occasion", index(a.occasion)}, {"sign", Json(a.sign)}}));
    for (const auto& g : K.rungs[t])
      rungs.push_back(Json(Obj{{"occasion", index(g.occasion)}, {"sign", Json(g.sign)}, {"attested", Json(g.attested)},
                               {"expected", Json(history::expected(g, m.of(x.colour)))}}));
    ties.push_back(Json(Obj{{"id", Json(x.id)}, {"a", index(x.a)}, {"b", index(x.b)}, {"type", Json(x.type)}, {"colour", Json(x.colour)},
                            {"confidence", Json(x.confidence)}, {"sources", texts(x.sources)}, {"note", text_or_null(x.note)},
                            {"drafted_by", text_or_null(x.drafted_by)}, {"checked_by", text_or_null(x.checked_by)},
                            {"attested", Json(std::move(attested))}, {"ends", index_or_null(x.ends)},
                            {"life", x.attested.empty() ? Json(nullptr) : Json(Arr{index(K.life[t][0]), index(K.life[t][1])})},
                            {"rungs", Json(std::move(rungs))}}));
  }
  answer.set("ties", Json(std::move(ties)));

  Arr problems;
  for (const auto* list : {&r.problems, &K.problems})
    for (const auto& p : *list)
      problems.push_back(Json(Obj{{"severity", Json(p.error ? "error" : "warning")}, {"where", Json(p.where)}, {"message", Json(p.message)}}));
  answer.set("problems", Json(std::move(problems)));

  // The complex.
  Arr edges, faces;
  for (const auto& s : K.segments) edges.push_back(Json(Obj{{"being", index(s.being)}, {"from", index(s.from)}, {"to", index(s.to)}}));
  for (std::size_t f = 0; f < K.lenses.size(); f++) {
    const auto& L = K.lenses[f];
    const auto &a = K.rungs[L.tie][L.rung_from], &b = K.rungs[L.tie][L.rung_to];
    Arr word;
    for (const auto& s : K.cw.faces[f].boundary) word.push_back(Json(static_cast<double>(s.sign) * static_cast<double>(s.edge + 1)));
    faces.push_back(Json(Obj{{"tie", index(L.tie)}, {"from", index(L.from)}, {"to", index(L.to)}, {"signs", Json(Arr{Json(a.sign), Json(b.sign)})},
                             {"attested", Json(Arr{Json(a.attested), Json(b.attested)})}, {"flip", Json(K.flip(f))}, {"boundary", Json(std::move(word))}}));
  }
  answer.set("complex", Json(Obj{{"vertices", indices(K.order)}, {"edges", Json(std::move(edges))}, {"faces", Json(std::move(faces))}, {"atlas", atlas_json(K)}}));

  // The invariants.
  const auto H = morphe::topology::homology(K.cw, true);
  Arr homology;
  std::vector<std::size_t> betti(3, 0);
  for (const auto& h : H) {
    if (h.dimension < 3) betti[h.dimension] = h.rank;
    homology.push_back(Json(Obj{{"dimension", index(h.dimension)}, {"rank", index(h.rank)},
                                {"torsion", array_of(h.torsion, [](const morphe::algebra::Integer& t) { return Json(t.str()); })}, {"text", Json(h.str())}}));
  }
  Json pi1 = Json(nullptr);
  if (K.cw.vertices) {
    const auto G = morphe::topology::fundamental_group(K.cw, 0);
    const auto& S = G.simplified.presentation;
    std::string text = S.relators.empty() ? (S.generators ? "free of rank " + std::to_string(S.generators) : "trivial")
                                          : std::to_string(S.generators) + " generators, " + std::to_string(S.relators.size()) + " relators";
    pi1 = Json(Obj{{"generators", index(G.presentation.generators)}, {"relators", index(G.presentation.relators.size())},
                   {"simplified_generators", index(S.generators)}, {"simplified_relators", index(S.relators.size())},
                   {"abelianisation", Json(G.abelianisation.str())}, {"text", Json(text)}});
  }
  const auto L = history::loops(K, H);
  const auto S = history::shells(K);
  answer.set("invariants",
             Json(Obj{{"counts", Json(Arr{index(K.cw.vertices), index(K.cw.edges.size()), index(K.cw.faces.size())})},
                      {"euler", Json(static_cast<double>(K.cw.euler_characteristic()))},
                      {"components", index(betti[0])},
                      {"homology", Json(std::move(homology))},
                      {"pi1", pi1},
                      {"loops", array_of(L, [&](const history::Loop& l) {
                         Arr steps;
                         for (const auto& s : l.steps) {
                           const auto& seg = K.segments[s.edge];
                           steps.push_back(Json(Obj{{"edge", index(s.edge)}, {"being", index(seg.being)},
                                                    {"from", index(s.forward ? seg.from : seg.to)}, {"to", index(s.forward ? seg.to : seg.from)},
                                                    {"forward", Json(s.forward)}}));
                         }
                         return Json(Obj{{"steps", Json(std::move(steps))}, {"beings", indices(l.beings)}, {"occasions", indices(l.occasions)}});
                       })},
                      {"shells", array_of(S, [](const history::Shell& s) {
                         auto three = [](const std::array<std::size_t, 3>& v) { return Json(Arr{index(v[0]), index(v[1]), index(v[2])}); };
                         return Json(Obj{{"beings", three(s.beings)}, {"ties", three(s.ties)}, {"faces", three(s.faces)}, {"from", index(s.from)},
                                         {"to", index(s.to)}, {"product_from", Json(s.product_from)}, {"product_to", Json(s.product_to)},
                                         {"flips", index(s.flips)}});
                       })},
                      {"shell_rank", index(betti[2])}}));

  // The readings.
  Json colours_lenses = Json::object();
  for (const auto& [colour, n] : C.lenses) colours_lenses.set(colour, index(n));
  answer.set("readings",
             Json(Obj{{"persistence", array_of(P, [](const history::Persistence& p) {
                         return Json(Obj{{"colour", Json(p.colour)}, {"lenses", index(p.lenses)}, {"flips", index(p.flips)}, {"p", Json(p.p)},
                                         {"memory", Json(p.memory)}, {"beta", p.beta ? Json(*p.beta) : Json(nullptr)}});
                       })},
                      {"flips", indices(history::flips(K))},
                      {"triangles", array_of(history::triangles(K, m), [](const history::Triangle& t) {
                         return Json(Obj{{"occasion", index(t.occasion)}, {"beings", Json(Arr{index(t.beings[0]), index(t.beings[1]), index(t.beings[2])})},
                                         {"ties", Json(Arr{index(t.ties[0]), index(t.ties[1]), index(t.ties[2])})}, {"product", Json(t.product)},
                                         {"expected", Json(t.expected)}, {"attested", Json(t.attested)}});
                       })},
                      {"moments", array_of(history::moments_by_year(K, m), [](const history::YearMoment& y) {
                         return Json(Obj{{"year", Json(static_cast<double>(y.year))}, {"jdn", Json(static_cast<double>(y.jdn))}, {"ties", index(y.ties)},
                                         {"triangles", index(y.triangles)}, {"frustrated", index(y.frustrated)}, {"expected", Json(y.expected)}});
                       })},
                      {"colours", Json(Obj{{"lenses", colours_lenses},
                                           {"carriers", array_of(C.carriers, [](const auto& c) {
                                              Json by = Json::object();
                                              for (const auto& [colour, n] : c.second) by.set(colour, index(n));
                                              return Json(Obj{{"being", index(c.first)}, {"colours", by}});
                                            })},
                                           {"meetings", array_of(C.meetings, [](const auto& c) {
                                              return Json(Obj{{"occasion", index(c.first)}, {"colours", texts(c.second)}});
                                            })}})}}));
  const Foam F = foam_of(K);
  answer.set("foam", foam_json(K, F));
  answer.set("models", models_json(K, F, S));
  return answer;
}

Json moment(const history::Record& record, const std::string& date, const Given& given) {
  const history::HistoryComplex K = history::construct(record);
  const history::Memory m = history::memory(history::persistence(K));
  const history::Date d = history::parse_historical(date);
  const history::Moment M = history::moment(K, d.last, m);
  Json answer = moment_json(M);
  answer.set("date", Json(date));
  // The posterior of the reconstruction at this moment: each link's inferred sign at its last
  // meeting, each triangle's mean product, and the law of the number frustrated; and, from the
  // prior, the spin networks along which the foam crosses this cut, with their weights.
  const Foam F = foam_of(K, given);
  const history::Counter counter = moment_counter(K, F.M, M);
  std::vector<std::size_t> link_rung;
  for (const auto& l : M.links) {
    std::size_t r = F.M.index[l.tie].front();
    for (std::size_t k = 0; k < K.rungs[l.tie].size(); k++)
      if (K.record.occasions[K.rungs[l.tie][k].occasion].date.first <= M.jdn) r = F.M.index[l.tie][k];
    link_rung.push_back(r);
  }
  Arr links, triangles, networks;
  double mean = 0, lo = 0, hi = 0;
  bool exact = F.inexact.empty();
  if (exact) {
    try {
      history::Cut cut;
      cut.after = cut_at(K, M.jdn);
      cut.keep = link_rung;
      for (const auto& t : counter.triangles) cut.keep.insert(cut.keep.end(), t.begin(), t.end());
      const auto post = history::contract(K, F.M, true, {}, {cut}, {counter});
      for (std::size_t r : link_rung) links.push_back(Json(post.link[r]));
      for (double x : post.cuts[0].products) triangles.push_back(Json(x));
      std::tie(mean, lo, hi) = std::tuple_cat(law_summary(post.cuts[0].negative, post.cuts[0].mean_negative));
      history::Cut open;
      open.after = cut.after;
      const auto prior = history::contract(K, F.M, false, {}, {open}, {}, 12);
      for (const auto& [set, weight] : prior.cuts[0].networks) {
        Arr ties;
        for (std::size_t r : set) ties.push_back(index(F.M.rung[r][0]));
        networks.push_back(Json(Obj{{"ties", Json(std::move(ties))}, {"weight", Json(weight)}}));
      }
    } catch (const history::Intractable&) {
      exact = false;
      links = Arr();
      triangles = Arr();
      networks = Arr();
    }
  }
  if (!exact) {
    history::Sampling S;
    S.burn = 1000;
    S.sweeps = 20000;
    history::Observables O;
    for (const auto& t : counter.triangles) O.rungs.push_back({t[0], t[1], t[2]});
    O.counters = {counter};
    const history::Ensemble post = history::sample(F.R, F.M, true, S, O);
    for (std::size_t r : link_rung) links.push_back(Json(post.rung[r]));
    for (double x : post.rung_product) triangles.push_back(Json(x));
    std::tie(mean, lo, hi) = std::tuple_cat(summary(post.count_hist[0]));
  }
  answer.set("posterior", Json(Obj{{"links", Json(std::move(links))}, {"triangles", Json(std::move(triangles))},
                                   {"frustrated", Json(Obj{{"mean", Json(mean)}, {"lo", Json(lo)}, {"hi", Json(hi)}})},
                                   {"exact", Json(exact)}}));
  answer.set("networks", Json(std::move(networks)));
  return answer;
}

// The power of each being over the reconstruction (the social spin foam's Def. 10.3): the act is
// the record without the being's ties, the couplings kept; its power is the divergence of the
// posterior marginals of every other unattested sign. And the power of balance itself: the act
// is β_Δ = 0.
Json power(const history::Record& record, const Given& given) {
  const history::HistoryComplex K = history::construct(record);
  const Foam F = foam_of(K, given);
  history::Sampling S;
  S.burn = 500;
  S.sweeps = 6000;
  // The posterior's marginals, exactly where the contraction reaches.
  bool exact = F.inexact.empty();
  auto marginals = [&](const history::HistoryComplex& KK, const history::HistoryModel& MM) {
    if (exact) {
      try {
        history::Ensemble E;
        E.rung = history::contract(KK, MM, true).link;
        return E;
      } catch (const history::Intractable&) {
        exact = false;
      }
    }
    return history::sample(history::refine(KK), MM, true, S);
  };
  const history::Ensemble base = marginals(K, F.M);
  // The unattested rungs, by (tie id, occasion id).
  std::map<std::pair<std::string, std::string>, double> mean;
  for (std::size_t r = 0; r < F.M.rung.size(); r++)
    if (!F.M.attested[r]) {
      const auto [t, k] = F.M.rung[r];
      mean[{K.record.ties[t].id, K.record.occasions[K.rungs[t][k].occasion].id}] = base.rung[r];
    }
  auto divergence = [&](const history::HistoryComplex& K2, const history::HistoryModel& M2, const history::Ensemble& E) {
    double d = 0;
    std::size_t n = 0;
    for (std::size_t r = 0; r < M2.rung.size(); r++) {
      if (M2.attested[r]) continue;
      const auto [t, k] = M2.rung[r];
      const auto it = mean.find({K2.record.ties[t].id, K2.record.occasions[K2.rungs[t][k].occasion].id});
      if (it == mean.end()) continue;
      const double p = std::clamp((1 + E.rung[r]) / 2, 1e-6, 1 - 1e-6), q = std::clamp((1 + it->second) / 2, 1e-6, 1 - 1e-6);
      d += p * std::log(p / q) + (1 - p) * std::log((1 - p) / (1 - q));
      n++;
    }
    return std::pair<double, std::size_t>{d, n};
  };
  // One reconstruction for each being: on a large record, for the forty with the most ties.
  std::vector<std::size_t> ties_of(record.beings.size(), 0), chosen;
  for (const auto& t : record.ties) { ties_of[t.a]++; ties_of[t.b]++; }
  for (std::size_t b = 0; b < record.beings.size(); b++)
    if (ties_of[b]) chosen.push_back(b);
  std::stable_sort(chosen.begin(), chosen.end(), [&](std::size_t x, std::size_t y) { return ties_of[x] > ties_of[y]; });
  const std::size_t considered = chosen.size();
  if (chosen.size() > 40) chosen.resize(40);
  Arr beings;
  for (std::size_t b : chosen) {
    history::Record without = record;
    std::vector<history::Tie> kept;
    std::size_t removed = 0;
    for (const auto& t : record.ties) {
      if (t.a == b || t.b == b) { removed++; continue; }
      kept.push_back(t);
    }
    if (!removed) continue;
    without.ties = kept;
    const history::HistoryComplex K2 = history::construct(without);
    const history::RefinedComplex R2 = history::refine(K2);
    const history::HistoryModel M2 = history::history_model(K2, R2, F.used);
    const auto [d, n] = divergence(K2, M2, marginals(K2, M2));
    beings.push_back(Json(Obj{{"being", index(b)}, {"ties", index(removed)}, {"divergence", Json(d)}, {"rungs", index(n)}}));
  }
  history::Couplings flat = F.used;
  flat.balance = 0;
  const history::HistoryModel M0 = history::history_model(K, F.R, flat);
  const auto [d0, n0] = divergence(K, M0, marginals(K, M0));
  return Json(Obj{{"beings", Json(std::move(beings))}, {"balance", Json(Obj{{"divergence", Json(d0)}, {"rungs", index(n0)}})},
                  {"sweeps", exact ? Json(nullptr) : index(S.sweeps)}, {"exact", Json(exact)},
                  {"considered", index(considered)}, {"computed", index(chosen.size())}});
}

// The history as a Feynman graph of a group field theory, and in time as a sample path (theory
// note §8): the stranded graph, its bubbles and their types, the check of Proposition 8.1, the
// amplitude, and the rates of the time-ordered process with its martingale checks.
Json field(const history::Record& record, const Given& given) {
  const history::HistoryComplex K = history::construct(record);
  const history::RefinedComplex R = history::refine(K);
  const history::StrandedGraph G = history::stranded(K);
  const auto forms = history::feynman_forms(K, R, G);
  const auto boundary = history::lens_boundaries(R, K.lenses.size());
  bool identity = forms.consistent && forms.free_strands == K.lenses.size() && history::pairing_complete(K, G);
  for (std::size_t f = 0; f < K.lenses.size(); f++) identity = identity && forms.argument[f] == boundary[f];
  // The field's ranks.
  std::map<std::size_t, std::size_t> ranks;
  std::size_t strands = 0;
  for (std::size_t e = 0; e < K.segments.size(); e++) { ranks[G.rank(e)]++; strands += G.rank(e); }
  Arr rank_list;
  for (const auto& [n, count] : ranks) rank_list.push_back(Json(Obj{{"rank", index(n)}, {"segments", index(count)}}));
  // The bubbles, and their types by frequency.
  std::vector<std::pair<std::string, const std::vector<std::size_t>*>> types;
  for (const auto& [t, vs] : G.types) types.push_back({t, &vs});
  std::stable_sort(types.begin(), types.end(), [](const auto& a, const auto& b) { return a.second->size() > b.second->size(); });
  std::map<std::string, std::size_t> type_index;
  for (std::size_t i = 0; i < types.size(); i++) type_index[types[i].first] = i;
  Arr bubbles;
  for (std::size_t v = 0; v < G.bubbles.size(); v++) {
    const auto& B = G.bubbles[v];
    Arr halves, corners;
    for (std::size_t h = 0; h < B.half_edges.size(); h++)
      halves.push_back(Json(Obj{{"being", index(B.owner[h])}, {"out", Json(B.half_edges[h].out)}, {"rank", index(G.rank(B.half_edges[h].segment))}}));
    for (const auto& c : B.corners)
      corners.push_back(Json(Obj{{"lens", index(c.lens)}, {"tie", index(K.lenses[c.lens].tie)},
                                 {"kind", Json(c.kind == history::Corner::passing ? "passing" : c.kind == history::Corner::closing ? "closing" : "opening")},
                                 {"ends", Json(Arr{index(c.ends[0]), index(c.ends[1])})}}));
    bubbles.push_back(Json(Obj{{"occasion", index(B.occasion)}, {"type", index(type_index.at(B.type))}, {"k_in", index(B.k_in)}, {"k_new", index(B.k_new)},
                               {"k_last", index(B.k_last)}, {"meetings", index(B.meetings)}, {"begun", index(B.begun)}, {"ended", index(B.ended)},
                               {"automorphisms", index(B.automorphisms)}, {"canonical", Json(B.canonical)},
                               {"half_edges", Json(std::move(halves))}, {"corners", Json(std::move(corners))}}));
  }
  Arr type_list;
  for (const auto& [t, vs] : types) {
    const auto& B = G.bubbles[vs->front()];
    type_list.push_back(Json(Obj{{"count", index(vs->size())}, {"example", index(vs->front())}, {"beings", index(B.beings.size())},
                                 {"automorphisms", index(B.automorphisms)}, {"canonical", Json(B.canonical)}}));
  }
  // The amplitude: Z of the spin foam, and the history's Feynman amplitude 2^{F + E} Z.
  const Foam F = foam_of(K, given);
  Json amplitude = Json(nullptr);
  if (F.inexact.empty()) {
    try {
      std::size_t largest = 0;
      const double prior = history::log_z(K, F.M, false, largest), post = history::log_z(K, F.M, true, largest);
      const double ln2 = std::log(2.0);
      amplitude = Json(Obj{{"log_z", Json(prior)},
                           {"log_amplitude", Json(prior + ln2 * static_cast<double>(K.lenses.size() + K.segments.size()))},
                           {"log_likelihood", Json(post - prior)},
                           {"largest", index(largest)}});
    } catch (const history::Intractable&) {
    }
  }
  Arr couplings;
  {
    const auto P = history::parameters_of(K, F.M);
    std::vector<bool> seen(P.names.size(), false);
    std::vector<std::size_t> cells(P.names.size(), 0);
    for (std::size_t p : P.of_cell) cells[p]++;
    for (std::size_t c = 0; c < F.M.cells.size(); c++) {
      const std::size_t p = P.of_cell[c];
      if (seen[p]) continue;
      seen[p] = true;
      const double b = F.M.cells[c].beta;
      couplings.push_back(Json(Obj{{"of", Json(P.names[p])}, {"cells", index(cells[p])}, {"beta", Json(b)}, {"t", Json(std::tanh(b))}}));
    }
  }
  // The process.
  const history::Process P = history::process_of(K);
  auto prop = [](const history::Proportion& q) {
    return Json(Obj{{"trials", index(q.trials)}, {"successes", index(q.successes)}, {"p", Json(q.p)}, {"error", Json(q.error)}});
  };
  Json process = Json(Obj{
      {"from", Json(P.from)}, {"to", Json(P.to)},
      {"rates", array_of(P.rates, [](const history::Rate& r) {
         return Json(Obj{{"k", index(r.k)}, {"count", index(r.count)}, {"exposure", Json(r.exposure)}, {"rate", Json(r.rate)},
                         {"error", Json(r.error)}, {"upper", Json(r.upper)}, {"mean_new", Json(r.mean_new)}});
       })},
      {"open", prop(P.open)}, {"close", prop(P.close)},
      {"binding", Json(Obj{{"tied", Json(P.mu_tied)}, {"untied", Json(P.mu_untied)}, {"ratio", Json(P.binding)}, {"low", Json(P.binding_low)},
                            {"high", Json(P.binding_high)}, {"tied_meetings", index(P.tied_meetings)}, {"untied_meetings", index(P.untied_meetings)},
                            {"tied_exposure", Json(P.tied_exposure)}, {"untied_exposure", Json(P.untied_exposure)}})},
      {"eras", array_of(P.eras, [](const history::EraCheck& c) {
         return Json(Obj{{"from", Json(c.from)}, {"to", Json(c.to)}, {"occasions", index(c.occasions)}, {"expected", Json(c.expected)}, {"z", Json(c.z)},
                         {"opened", index(c.opened)}, {"opened_expected", Json(c.opened_expected)}, {"opened_z", Json(c.opened_z)}});
       })},
      {"chi2", Json(P.chi2)}, {"dof", index(P.dof)},
      {"log_likelihood", Json(Obj{{"occasions", Json(P.log_likelihood_occasions)}, {"strands", Json(P.log_likelihood_strands)}})}});
  return Json(Obj{{"segments", index(K.segments.size())}, {"lenses", index(K.lenses.size())}, {"occasions", index(K.order.size())},
                  {"strands", index(strands)}, {"ranks", Json(std::move(rank_list))},
                  {"identity", Json(Obj{{"holds", Json(identity)}, {"free_strands", index(forms.free_strands)}})},
                  {"bubbles", Json(std::move(bubbles))}, {"types", Json(std::move(type_list))},
                  {"amplitude", amplitude}, {"couplings", Json(std::move(couplings))},
                  {"exact", Json(F.inexact.empty())}, {"process", process}});
}

// The record a request names: one canon, or several joined.
history::Record load(const Json& request) {
  if (const Json* list = request.find("canons"); list && list->is_array()) {
    std::vector<std::string> paths;
    for (const auto& x : list->as_array())
      if (x.is_string()) paths.push_back(x.as_string());
    if (paths.empty()) throw std::invalid_argument("\"canons\" names no canon");
    return history::join_records(paths);
  }
  const Json* c = request.find("canon");
  if (!c || !c->is_string() || c->as_string().empty()) throw std::invalid_argument("the request has no \"canon\"");
  return history::read_record_file(c->as_string());
}

std::string text(const Json& request, const char* key) {
  const Json* v = request.find(key);
  if (!v || !v->is_string() || v->as_string().empty()) throw std::invalid_argument(std::string("the request has no \"") + key + "\"");
  return v->as_string();
}

}  // namespace

int main() {
  std::ios::sync_with_stdio(false);
  const std::string input((std::istreambuf_iterator<char>(std::cin)), std::istreambuf_iterator<char>());
  try {
    const Json request = Json::parse(input);
    const std::string op = text(request, "op");
    Json answer;
    if (op == "build") answer = build(load(request));
    else if (op == "moment") answer = moment(load(request), text(request, "date"), couplings_in(request));
    else if (op == "power") answer = power(load(request), couplings_in(request));
    else if (op == "field") answer = field(load(request), couplings_in(request));
    else throw std::invalid_argument("no operation \"" + op + "\"");
    std::cout << answer.dump() << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cout << Json(Obj{{"error", Json(e.what())}}).dump() << "\n";
    return 1;
  }
}
