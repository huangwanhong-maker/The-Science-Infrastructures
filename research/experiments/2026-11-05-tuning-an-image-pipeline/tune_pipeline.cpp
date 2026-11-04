// Tuning an image pipeline (README.md).
//
// eikon's chain for tubes in an image (Frangi's vesselness over a range of scales, hysteresis on its response,
// growth into the smoothed image above an intensity level) is tuned on synthetic phantoms of exact geometry
// (eikon::phantoms: a tree of round cones keeping Murray's law, rasterised with partial volumes, Gaussian noise),
// with eight parameters:
//
//   NSGA-II for two objectives, 1 − Dice (overlap of volumes) and 1 − clDice (overlap of centrelines,
//   Shit et al. 2021), against random search with the same budget;
//   CMA-ES for one score, 1 − (Dice + clDice)/2, against random search;
//
// on two training phantoms, and the chosen chains held to three phantoms not seen in tuning, with Betti numbers
// and surface distances. Every number of the record is written to <results>/ by this program.
//
//   tune_pipeline <results directory> [time|all]

#include <heuresis/multiobjective/multiobjective.hpp>
#include <heuresis/search/search.hpp>
#include <heuresis/strategies/strategies.hpp>
#include <heuresis/variation/variation.hpp>

#include <eikon/features/vesselness.hpp>
#include <eikon/filters/filters.hpp>
#include <eikon/fusion/fusion.hpp>
#include <eikon/phantoms/phantoms.hpp>
#include <eikon/segmentation/segmentation.hpp>
#include <eikon/skeleton/skeleton.hpp>

#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace hs = heuresis::search;
namespace hv = heuresis::variation;
namespace hm = heuresis::multiobjective;
namespace eg = eikon::grid;
namespace ep = eikon::phantoms;
namespace ef = eikon::filters;
namespace ev = eikon::features;
namespace es = eikon::segmentation;
namespace efu = eikon::fusion;
using Real = hv::Real;

namespace {

std::filesystem::path out_dir;
double seconds_since(std::chrono::steady_clock::time_point t0) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }

// The ladder of scales: σ_k = 0.35 mm · 2^{k/4}, k = 0 … 12 (0.35 to 2.8 mm; radii √2σ from 0.5 to 4 mm).
constexpr int K = 13;
double ladder(int k) { return 0.35 * std::pow(2.0, k / 4.0); }

double noise_sd = 0.5;              // on a contrast of 1 (varied in the time mode)
constexpr double smoothing = 0.5;  // mm, the image grown into
constexpr int blobs = 4;           // balls of the vessels' intensity, not vessels

// A phantom: eikon's tree of tubes (as its skeleton's check builds it: depth 3, 7 junctions and 9 ends) on a grid
// of 72³ voxels of 0.5 mm, with four balls of radius 2.5 to 3.5 mm placed at random at least 2 mm from the tree
// (blobs, which a vesselness measure is meant to reject), their occupancy as intensity (background 0, tubes and
// balls 1) with Gaussian noise of deviation 0.5; the truth is the set of voxels whose centres lie inside the
// tubes. The Hessian's eigenvalues at every scale of the ladder are computed once.
struct Phantom {
  std::uint64_t seed = 0;
  eg::Image<float> image, smoothed;
  es::Mask truth, blob;
  std::vector<ev::Eigenvalues> E;
  double truth_voxels = 0;
};

eg::Geometry box(std::size_t n, double h) {
  eg::Geometry g;
  g.size = {n, n, n};
  g.spacing = {h, h, h};
  const double half = h * double(n - 1) / 2;
  g.origin = {-half, -half, -half};
  return g;
}

Phantom make_phantom(std::uint64_t seed) {
  Phantom P;
  P.seed = seed;
  const auto g = box(72, 0.5);
  ep::TreeOptions o;
  o.depth = 3;
  o.angle = 0.7;
  o.length_ratio = 0.9;
  calx::random::Stream shape = hs::Draws(seed)(hs::tag("study.phantom.tree"), 0, 0);
  const auto tree = ep::tree({0, 0, -16}, {0.1, 0.05, 1}, 10, 2.0, o, shape);
  auto shapes = ep::shapes(tree);
  const auto vessels = ep::rasterise(g, shapes, 4);
  // The balls: centres uniform in the inner box (±14 mm), kept when at least r + 2 mm from every tube and ball.
  calx::random::Stream place = hs::Draws(seed)(hs::tag("study.phantom.blobs"), 0, 0);
  std::vector<ep::Tube> balls;
  while (int(balls.size()) < blobs) {
    const ep::Vec3 c{place.uniform(-14, 14), place.uniform(-14, 14), place.uniform(-14, 14)};
    const double r = place.uniform(2.5, 3.5);
    bool clear = true;
    for (const auto& t : tree.tubes) clear = clear && ep::distance(t, c) > r + 2;
    for (const auto& b : balls) clear = clear && ep::distance(b, c) > r + 2;
    if (clear) balls.push_back({c, c, r, r, 100 + int(balls.size())});
  }
  for (const auto& b : balls) shapes.push_back(ep::shape(b));
  const auto raster = ep::rasterise(g, shapes, 4);
  P.image = ep::intensities(raster.occupancy, 0, 1);
  calx::random::Stream noise = hs::Draws(seed)(hs::tag("study.phantom.noise"), 0, 0);
  ep::gaussian_noise(P.image, noise_sd, noise);
  P.truth = es::Mask(g);
  P.blob = es::Mask(g);
  for (std::size_t v = 0; v < g.count(); v++) {
    P.truth.data[v] = vessels.distance.data[v] <= 0, P.truth_voxels += P.truth.data[v];
    P.blob.data[v] = raster.distance.data[v] <= 0 && !P.truth.data[v];
  }
  P.smoothed = ef::smooth(P.image, smoothing);
  for (int k = 0; k < K; k++) P.E.push_back(ev::hessian_eigenvalues(ef::hessian(P.image, ladder(k), true), ev::Polarity::bright));
  return P;
}

// The chain's eight parameters from the unit cube:
//   scales: lo = round(12x₀), hi = lo + round((12 − lo)x₁), on the ladder;
//   Frangi's α = 0.1 + 1.4x₂, β = 0.1 + 1.4x₃, c = 10^{−3 + 3x₄} (0.001 to 1, on a contrast of 1);
//   hysteresis on the response: high = 0.02 + 0.88x₅, low = (0.05 + 0.95x₆)·high;
//   growth from the detected voxels into the smoothed image above level = 0.1 + 0.8x₇.
struct Chain {
  int lo = 0, hi = 0;
  double alpha = 0.5, beta = 0.5, c = 0.1, high = 0.3, low = 0.1, level = 0.5;
};
Chain decode(const Real& x) {
  Chain ch;
  ch.lo = int(std::lround(12 * std::clamp(x[0], 0.0, 1.0)));
  ch.hi = ch.lo + int(std::lround((12 - ch.lo) * std::clamp(x[1], 0.0, 1.0)));
  ch.alpha = 0.1 + 1.4 * x[2];
  ch.beta = 0.1 + 1.4 * x[3];
  ch.c = std::pow(10.0, -3 + 3 * x[4]);
  ch.high = 0.02 + 0.88 * x[5];
  ch.low = (0.05 + 0.95 * x[6]) * ch.high;
  ch.level = 0.1 + 0.8 * x[7];
  return ch;
}

es::Mask reconstruct(const Phantom& P, const Chain& ch, eg::Image<float>* response = nullptr) {
  eg::Image<float> V(P.image.geometry, 0.0f);
  for (int k = ch.lo; k <= ch.hi; k++) {
    const auto v = ev::frangi(P.E[std::size_t(k)], ch.alpha, ch.beta, ch.c);
    for (std::size_t n = 0; n < V.data.size(); n++) V.data[n] = std::max(V.data[n], v.data[n]);
  }
  const auto detected = es::hysteresis(V, ch.low, ch.high, 26);
  std::vector<std::size_t> seeds;
  for (std::size_t n = 0; n < detected.data.size(); n++)
    if (detected.data[n]) seeds.push_back(n);
  const double level = ch.level;
  if (response) *response = std::move(V);
  return es::grow(P.smoothed, seeds, [level](float x) { return double(x) >= level; }, 26);
}

struct Scores {
  double dice = 0, cldice = 0;
};
Scores score(const Phantom& P, const es::Mask& m) { return {efu::dice(m, P.truth), efu::cldice(m, P.truth).cldice}; }

// The objectives on the training phantoms: the means of 1 − Dice and 1 − clDice.
std::vector<const Phantom*> training;
Scores train_scores(const Chain& ch) {
  Scores s;
  for (const auto* P : training) {
    const Scores q = score(*P, reconstruct(*P, ch));
    s.dice += q.dice / double(training.size());
    s.cldice += q.cldice / double(training.size());
  }
  return s;
}

hs::Evaluation two_objectives(const Real& x) {
  const Scores s = train_scores(decode(x));
  return hs::evaluation({1 - s.dice, 1 - s.cldice});
}
hs::Evaluation one_objective(const Real& x) {
  const Scores s = train_scores(decode(x));
  hs::Evaluation e = hs::scalar(1 - 0.5 * (s.dice + s.cldice));
  e.descriptors = {1 - s.dice, 1 - s.cldice};  // kept beside the score, for the record
  return e;
}

const hv::Box unit8 = hv::Box::cube(8, 0, 1);

void write_record(FILE* f, const char* who, std::uint64_t seed, const hs::Record<Real>& r) {
  std::size_t i = 0;
  for (const auto& x : r.told) {
    const Chain ch = decode(x.genotype);
    const auto& e = x.evaluation;
    const double d = e.objectives.size() == 2 ? e.objectives[0] : e.descriptors[0], c = e.objectives.size() == 2 ? e.objectives[1] : e.descriptors[1];
    std::fprintf(f, "%s,%llu,%zu,%.6f,%.6f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", who, static_cast<unsigned long long>(seed), ++i, 1 - d, 1 - c, ladder(ch.lo), ladder(ch.hi), ch.alpha, ch.beta,
                 ch.c, ch.high, ch.low, ch.level);
  }
}

// The hypervolume (reference (1, 1)) of the non-dominated points among the first m evaluations, every `step`.
std::vector<std::pair<std::size_t, double>> anytime_hypervolume(const hs::Record<Real>& r, std::size_t step, bool two) {
  std::vector<std::pair<std::size_t, double>> out;
  std::vector<hm::Point> pts;
  const std::vector<double> ref{1, 1};
  for (std::size_t i = 0; i < r.told.size(); i++) {
    const auto& e = r.told[i].evaluation;
    pts.push_back(two ? hm::Point{e.objectives[0], e.objectives[1]} : hm::Point{e.descriptors[0], e.descriptors[1]});
    if ((i + 1) % step == 0 || i + 1 == r.told.size()) out.push_back({i + 1, hm::hypervolume(pts, ref)});
  }
  return out;
}

struct Chosen {
  std::string name;
  Chain chain;
};

void run_all(calx::scheduler::ThreadPool& pool, const std::vector<Phantom>& tests) {
  const std::size_t mu = 32, generations = 24, nsga_budget = mu * (generations + 1);
  const std::size_t cma_budget = 320, seeds = 3;
  hs::Pool executor(pool);
  FILE* fe = std::fopen((out_dir / "evaluations.csv").string().c_str(), "w");
  std::fprintf(fe, "algorithm,seed,index,dice,cldice,sigma_lo,sigma_hi,alpha,beta,c,high,low,level\n");
  FILE* fh = std::fopen((out_dir / "hypervolume.csv").string().c_str(), "w");
  std::fprintf(fh, "algorithm,seed,evaluations,hypervolume\n");
  FILE* fc = std::fopen((out_dir / "convergence.csv").string().c_str(), "w");
  std::fprintf(fc, "algorithm,seed,evaluations,best\n");
  std::map<std::string, double> walls;

  // NSGA-II and random search with the same budget, for the two objectives.
  hs::Record<Real> rec_nsga, rec_random;
  {
    const auto t0 = std::chrono::steady_clock::now();
    hm::Nsga2<Real> nsga(hm::real_kernels(unit8), mu, hs::Draws(41));
    hs::run(nsga, two_objectives, hs::Budget{double(nsga_budget)}, executor, hs::unlimited, &rec_nsga);
    walls["NSGA-II"] = seconds_since(t0);
    const auto t1 = std::chrono::steady_clock::now();
    hs::RandomSearch<Real> rs([](calx::random::Stream& s) { return hv::uniform_point(unit8, s); }, hs::Draws(42));
    hs::run(rs, two_objectives, hs::Budget{double(nsga_budget)}, executor, 32, &rec_random);
    walls["random (two objectives)"] = seconds_since(t1);
    write_record(fe, "NSGA-II", 1, rec_nsga);
    write_record(fe, "random", 1, rec_random);
    for (const auto& [m, h] : anytime_hypervolume(rec_nsga, 32, true)) std::fprintf(fh, "NSGA-II,1,%zu,%.6f\n", m, h);
    for (const auto& [m, h] : anytime_hypervolume(rec_random, 32, true)) std::fprintf(fh, "random,1,%zu,%.6f\n", m, h);
    std::printf("NSGA-II and random: %zu evaluations each, %.1f s and %.1f s\n", nsga_budget, walls["NSGA-II"], walls["random (two objectives)"]);
  }

  // CMA-ES and random search for the one score, three seeds each.
  std::vector<hs::Record<Real>> rec_cma(seeds), rec_rs1(seeds);
  for (std::uint64_t sd = 1; sd <= seeds; sd++) {
    const auto t0 = std::chrono::steady_clock::now();
    heuresis::strategies::CmaOptions o;
    o.bounds = unit8;
    heuresis::strategies::Cma cma(Real(8, 0.5), 0.3, o, hs::Draws(50 + sd));
    const auto out = hs::run(cma, one_objective, hs::Budget{double(cma_budget)}, executor, hs::unlimited, &rec_cma[sd - 1]);
    walls["CMA-ES " + std::to_string(sd)] = seconds_since(t0);
    for (const auto& p : out.progress) std::fprintf(fc, "CMA-ES,%llu,%.0f,%.6f\n", static_cast<unsigned long long>(sd), p.evaluations, p.value);
    hs::RandomSearch<Real> rs([](calx::random::Stream& s) { return hv::uniform_point(unit8, s); }, hs::Draws(60 + sd));
    const auto out2 = hs::run(rs, one_objective, hs::Budget{double(cma_budget)}, executor, 32, &rec_rs1[sd - 1]);
    for (const auto& p : out2.progress) std::fprintf(fc, "random,%llu,%.0f,%.6f\n", static_cast<unsigned long long>(sd), p.evaluations, p.value);
    write_record(fe, "CMA-ES", sd, rec_cma[sd - 1]);
    write_record(fe, "random (score)", sd, rec_rs1[sd - 1]);
    std::printf("CMA-ES seed %llu: best %.4f; random %.4f\n", static_cast<unsigned long long>(sd), out.best->evaluation.value(), out2.best->evaluation.value());
  }
  // NSGA-II's anytime best of the same score, for comparison.
  {
    double best = 1e9;
    for (std::size_t i = 0; i < rec_nsga.told.size(); i++) {
      const auto& e = rec_nsga.told[i].evaluation;
      const double v = 0.5 * (e.objectives[0] + e.objectives[1]);
      if (v < best) best = v, std::fprintf(fc, "NSGA-II,1,%zu,%.6f\n", i + 1, v);
    }
  }
  std::fclose(fe);
  std::fclose(fh);
  std::fclose(fc);

  // The non-dominated chains of NSGA-II's run and of random search's (on the training phantoms), each held to the
  // test phantoms.
  {
    FILE* ff = std::fopen((out_dir / "front_test.csv").string().c_str(), "w");
    std::fprintf(ff, "algorithm,index,train_dice,train_cldice,phantom,dice,cldice,b0,blob_fraction,sigma_lo,sigma_hi,alpha,beta,c,high,low,level\n");
    for (const auto* rec : {&rec_nsga, &rec_random}) {
      std::vector<hm::Point> pts;
      for (const auto& x : rec->told) pts.push_back(x.evaluation.objectives);
      for (auto i : hm::non_dominated(pts)) {
        const Chain ch = decode(rec->told[i].genotype);
        for (const auto& P : tests) {
          const auto m = reconstruct(P, ch);
          const Scores q = score(P, m);
          double in_blob = 0, blob = 0;
          for (std::size_t v = 0; v < m.data.size(); v++) blob += P.blob.data[v], in_blob += P.blob.data[v] && m.data[v];
          std::fprintf(ff, "%s,%zu,%.6f,%.6f,%llu,%.6f,%.6f,%ld,%.4f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", rec == &rec_nsga ? "NSGA-II" : "random", i + 1, 1 - pts[i][0], 1 - pts[i][1],
                       static_cast<unsigned long long>(P.seed), q.dice, q.cldice, efu::betti(m).b0, blob > 0 ? in_blob / blob : 0.0, ladder(ch.lo), ladder(ch.hi), ch.alpha, ch.beta, ch.c,
                       ch.high, ch.low, ch.level);
        }
      }
    }
    std::fclose(ff);
  }

  // The chains held to the test phantoms: the hand-set chain (Frangi's defaults, α = β = 0.5 and c half the largest
  // Hessian norm over the scales, σ from 0.5 to 2 mm, hysteresis 0.1–0.3, growth at the half-maximum 0.5); NSGA-II's
  // front at its best Dice, its best clDice and its best mean; CMA-ES's best of each seed; random search's best.
  std::vector<Chosen> chosen;
  {
    Chain d;
    d.lo = 2, d.hi = 10, d.alpha = 0.5, d.beta = 0.5, d.high = 0.3, d.low = 0.1, d.level = 0.5;
    std::vector<double> sc;
    for (int k = d.lo; k <= d.hi; k++) sc.push_back(ladder(k));
    d.c = ev::frangi_c(training[0]->image, sc);
    chosen.push_back({"hand-set", d});
  }
  auto pick = [&](const hs::Record<Real>& r, auto key) {
    std::size_t best = 0;
    for (std::size_t i = 1; i < r.told.size(); i++)
      if (key(r.told[i].evaluation) < key(r.told[best].evaluation)) best = i;
    return decode(r.told[best].genotype);
  };
  chosen.push_back({"NSGA-II best Dice", pick(rec_nsga, [](const hs::Evaluation& e) { return e.objectives[0] + 1e-9 * e.objectives[1]; })});
  chosen.push_back({"NSGA-II best clDice", pick(rec_nsga, [](const hs::Evaluation& e) { return e.objectives[1] + 1e-9 * e.objectives[0]; })});
  chosen.push_back({"NSGA-II best mean", pick(rec_nsga, [](const hs::Evaluation& e) { return e.objectives[0] + e.objectives[1]; })});
  for (std::uint64_t sd = 1; sd <= seeds; sd++)
    chosen.push_back({"CMA-ES seed " + std::to_string(sd), pick(rec_cma[sd - 1], [](const hs::Evaluation& e) { return e.value(); })});
  chosen.push_back({"random best mean", pick(rec_random, [](const hs::Evaluation& e) { return e.objectives[0] + e.objectives[1]; })});

  FILE* ft = std::fopen((out_dir / "held_out.csv").string().c_str(), "w");
  std::fprintf(ft, "chain,set,phantom,dice,cldice,b0,b1,b2,surface_mean_mm,surface_p95_mm,voxels,truth_voxels,blob_fraction,sigma_lo,sigma_hi,alpha,beta,c,high,low,level\n");
  std::vector<const Phantom*> all = training;
  for (const auto& P : tests) all.push_back(&P);
  for (const auto& c : chosen)
    for (std::size_t i = 0; i < all.size(); i++) {
      const Phantom& P = *all[i];
      const auto m = reconstruct(P, c.chain);
      const Scores s = score(P, m);
      const auto b = efu::betti(m);
      const auto sd = efu::surface_distances(m, P.truth);
      const double vox = std::accumulate(m.data.begin(), m.data.end(), 0.0);
      double in_blob = 0, blob = 0;
      for (std::size_t v = 0; v < m.data.size(); v++) blob += P.blob.data[v], in_blob += P.blob.data[v] && m.data[v];
      const auto& ch = c.chain;
      std::fprintf(ft, "%s,%s,%llu,%.5f,%.5f,%ld,%ld,%ld,%.4f,%.4f,%.0f,%.0f,%.4f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", c.name.c_str(), i < training.size() ? "training" : "test",
                   static_cast<unsigned long long>(P.seed), s.dice, s.cldice, b.b0, b.b1, b.b2, sd.mean, sd.p95, vox, P.truth_voxels, blob > 0 ? in_blob / blob : 0.0, ladder(ch.lo), ladder(ch.hi), ch.alpha, ch.beta, ch.c,
                   ch.high, ch.low, ch.level);
    }
  std::fclose(ft);
  {
    const Phantom& P = *training[0];
    const auto b = efu::betti(P.truth);
    std::printf("truth of phantom %llu: Betti (%ld, %ld, %ld), %.0f voxels\n", static_cast<unsigned long long>(P.seed), b.b0, b.b1, b.b2, P.truth_voxels);
  }

  // Projections along y of the second test phantom (seed 4): the image's maximum, the truth, and four chains'
  // reconstructions.
  {
    const Phantom& P = tests[1];
    const auto& g = P.image.geometry;
    FILE* fp = std::fopen((out_dir / "projections.csv").string().c_str(), "w");
    std::fprintf(fp, "layer,k,values\n");
    auto project = [&](const char* layer, auto value) {
      for (std::size_t k = 0; k < g.size[2]; k++) {
        std::fprintf(fp, "%s,%zu,", layer, k);
        for (std::size_t i = 0; i < g.size[0]; i++) {
          double v = -1e9;
          for (std::size_t j = 0; j < g.size[1]; j++) v = std::max(v, double(value(g.linear(i, j, k))));
          std::fprintf(fp, i ? " %.3g" : "%.3g", v);
        }
        std::fprintf(fp, "\n");
      }
    };
    project("image", [&](std::size_t n) { return P.image.data[n]; });
    project("truth", [&](std::size_t n) { return P.truth.data[n]; });
    // The middle slice (j = 36) of the image, as one sees it.
    for (std::size_t k = 0; k < g.size[2]; k++) {
      std::fprintf(fp, "slice,%zu,", k);
      for (std::size_t i = 0; i < g.size[0]; i++) std::fprintf(fp, i ? " %.3g" : "%.3g", double(P.image.data[g.linear(i, 36, k)]));
      std::fprintf(fp, "\n");
    }
    for (const auto& c : chosen) {
      if (c.name != "hand-set" && c.name != "NSGA-II best Dice" && c.name != "NSGA-II best mean" && c.name != "CMA-ES seed 2") continue;
      const auto m = reconstruct(P, c.chain);
      project(c.name.c_str(), [&](std::size_t n) { return m.data[n]; });
    }
    std::fclose(fp);
  }
  FILE* fj = std::fopen((out_dir / "summary.json").string().c_str(), "w");
  std::fprintf(fj, "{\n  \"nsga_budget\": %zu, \"nsga_mu\": %zu, \"cma_budget\": %zu, \"seeds\": %zu,\n", nsga_budget, mu, cma_budget, seeds);
  std::fprintf(fj, "  \"hypervolume_nsga\": %.6f, \"hypervolume_random\": %.6f,\n", anytime_hypervolume(rec_nsga, nsga_budget, true).back().second,
               anytime_hypervolume(rec_random, nsga_budget, true).back().second);
  std::fprintf(fj, "  \"walls_s\": {");
  bool first = true;
  for (const auto& [k, v] : walls) std::fprintf(fj, "%s\"%s\": %.1f", first ? "" : ", ", k.c_str(), v), first = false;
  std::fprintf(fj, "}\n}\n");
  std::fclose(fj);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: tune_pipeline <results directory> [time|all]\n");
    return 2;
  }
  out_dir = argv[1];
  std::filesystem::create_directories(out_dir);
  const std::string mode = argc > 2 ? argv[2] : "all";
  calx::scheduler::ThreadPool pool(std::max(1u, std::thread::hardware_concurrency()));
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<Phantom> train_set, tests;
  for (std::uint64_t s : {1, 2}) train_set.push_back(make_phantom(s));
  for (const auto& P : train_set) training.push_back(&P);
  std::printf("two phantoms in %.1f s\n", seconds_since(t0));
  if (mode == "time") {
    {
      const Phantom& Q = train_set[0];
      Chain d;
      d.lo = 2, d.hi = 10;
      std::vector<double> sc;
      for (int k = d.lo; k <= d.hi; k++) sc.push_back(ladder(k));
      d.c = ev::frangi_c(Q.image, sc);
      for (double high : {0.1, 0.3, 0.6}) {
        d.high = high, d.low = high / 3;
        const auto m = reconstruct(Q, d);
        const Scores s = score(Q, m);
        double in_blob = 0, blob = 0;
        for (std::size_t v = 0; v < m.data.size(); v++) blob += Q.blob.data[v], in_blob += Q.blob.data[v] && m.data[v];
        std::printf("noise %.2f, high %.2f: hand-set Dice %.4f clDice %.4f (c %.4f), blobs taken %.3f\n", noise_sd, high, s.dice, s.cldice, d.c, in_blob / blob);
      }
    }
    Chain d;
    d.lo = 2, d.hi = 10;
    const auto t1 = std::chrono::steady_clock::now();
    const Scores s = train_scores(d);
    std::printf("one evaluation (two phantoms): %.2f s; Dice %.4f clDice %.4f\n", seconds_since(t1), s.dice, s.cldice);
    const auto t2 = std::chrono::steady_clock::now();
    const auto m = reconstruct(train_set[0], d);
    std::printf("  reconstruct %.3f s\n", seconds_since(t2));
    const auto t3 = std::chrono::steady_clock::now();
    (void)efu::cldice(m, train_set[0].truth);
    std::printf("  clDice %.3f s\n", seconds_since(t3));
    return 0;
  }
  for (std::uint64_t s : {3, 4, 5}) tests.push_back(make_phantom(s));
  run_all(pool, tests);
  std::printf("total %.1f s\n", seconds_since(t0));
  return 0;
}
