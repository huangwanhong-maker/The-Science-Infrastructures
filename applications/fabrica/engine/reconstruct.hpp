// Fabrica's reconstructions (ADR 0070): vessels reconstructed from images by the medical-imaging domain, given
// to the page as surfaces, centrelines and slices, with their agreement with a reference.
//
// Cases:
//   - the Medical Segmentation Decathlon's hepatic vessels (the data library's product msd-hepatic-vessels):
//     portal-venous CT, reconstructed within the annotated territory and compared with the experts' labels;
//   - "simulated": one body imaged twice by the modalities' forward models, a contrast CT and a TOF-MRA of the
//     patient moved between scans, integrated: registered, fused, arteries told from veins, against the truth.
// Positions are in the patient frame LPS (NIfTI's RAS carried across), in millimetres.
//
//   {"op": "cases"}
//   {"op": "reconstruct", "case": "hepaticvessel_166"}   surfaces, the tree and its measures, the agreement
//   {"op": "slice", "case": "…", "k": 40}                 a slice, windowed, with the masks' codes per pixel

#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <numbers>
#include <string>
#include <vector>

#include <anatomy/vessels/vessels.hpp>
#include <calx/random/stream.hpp>
#include <eikon/fusion/fusion.hpp>
#include <eikon/grid/nifti.hpp>
#include <eikon/phantoms/phantoms.hpp>
#include <medical_imaging/dicom/dicom.hpp>
#include <medical_imaging/integration/integration.hpp>
#include <medical_imaging/modalities/modalities.hpp>
#include <medical_imaging/reconstruction/reconstruction.hpp>
#include <science/json.hpp>

namespace fabrica {

using science::Json;
namespace eg = eikon::grid;
namespace mr = medical_imaging::reconstruction;
namespace mi = medical_imaging::integration;
namespace mm = medical_imaging::modalities;
namespace ep = eikon::phantoms;
namespace av = anatomy::vessels;

struct Case {
  std::string id, kind;           // "ct" or "simulated"
  eg::Image<float> image;         // the reference image (LPS geometry)
  mr::Mask reference, found, arteries, veins, region;
  Json result;
};

struct Reconstructions {
  std::filesystem::path products;
  std::map<std::string, Case> cases;
};

inline double round1(double x) { return std::round(x * 10) / 10; }
inline double round3(double x) { return std::round(x * 1000) / 1000; }

inline std::string base64(const std::vector<std::uint8_t>& b) {
  static const char* A = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  std::string out;
  for (std::size_t i = 0; i < b.size(); i += 3) {
    const unsigned v = unsigned(b[i]) << 16 | (i + 1 < b.size() ? unsigned(b[i + 1]) << 8 : 0) | (i + 2 < b.size() ? b[i + 2] : 0);
    out += A[v >> 18 & 63], out += A[v >> 12 & 63];
    out += i + 1 < b.size() ? A[v >> 6 & 63] : '=';
    out += i + 2 < b.size() ? A[v & 63] : '=';
  }
  return out;
}

// A mask's surface for drawing: closed by a voxel, smoothed, islands dropped; flattened and rounded.
inline Json surface_json(const mr::Mask& m, double min_volume = 20) {
  const auto& g = m.geometry;
  const double h = std::min({g.spacing[0], g.spacing[1], g.spacing[2]});
  const auto closed = eikon::distance::close(m, std::max({g.spacing[0], g.spacing[1], g.spacing[2]}));
  eg::Image<float> f(g);
  for (std::size_t v = 0; v < f.data.size(); v++) f.data[v] = closed.data[v];
  const auto S = eikon::surfaces::without_small_pieces(eikon::surfaces::smoothed(eikon::surfaces::isosurface(eikon::filters::smooth(f, 0.5 * h), 0.5), 10), min_volume);
  Json pts = Json::array(), tri = Json::array();
  for (const auto& p : S.points)
    for (int a = 0; a < 3; a++) pts.push(round1(p[a]));
  for (const auto& t : S.triangles)
    for (int a = 0; a < 3; a++) tri.push(double(t[std::size_t(a)]));
  return Json::object({{"points", pts}, {"triangles", tri}});
}

// The tree of a mask: centrelines with radii, rooted at the widest end, and its measures.
inline Json tree_json(const mr::Mask& m) {
  const auto D = eikon::distance::distance_to(eikon::distance::complement(m));
  const auto G = eikon::skeleton::graph(eikon::skeleton::prune(eikon::skeleton::thin(m), D), D);
  auto t = mr::tree(G);
  std::size_t root = 0;
  double widest = -1;
  for (std::size_t n = 0; n < G.nodes.size(); n++)
    if (G.nodes[n].kind == eikon::skeleton::NodeKind::end && G.nodes[n].radius > widest) widest = G.nodes[n].radius, root = n;
  Json segs = Json::array();
  double total = 0, murray = std::nan("");
  if (!t.nodes.empty()) {
    const auto R = av::root_at(t, root);
    const auto bs = av::bifurcations(t, R);
    if (bs.size() >= 3) murray = av::murray_exponent(bs);
    if (murray < 1.01 || murray > 4.99) murray = std::nan("");  // at the search's bounds: not determined by these bifurcations
    for (std::size_t e = 0; e < t.segments.size(); e++) {
      const auto& s = t.segments[e];
      Json pts = Json::array();
      for (const auto& p : s.points)
        for (int a = 0; a < 3; a++) pts.push(round1(p[a]));
      total += av::length(s);
      segs.push(Json::object({{"points", pts}, {"radius", round3(av::radius(s))}, {"length", round1(av::length(s))},
                              {"tortuosity", round3(av::tortuosity(s))}, {"generation", R.reached[e] ? R.generation[e] : -1},
                              {"strahler", R.reached[e] ? R.strahler[e] : -1}}));
    }
  }
  return Json::object({{"segments", segs}, {"junctions", double(G.count(eikon::skeleton::NodeKind::junction))},
                       {"ends", double(G.count(eikon::skeleton::NodeKind::end))}, {"length", round1(total)},
                       {"components", double(G.components())}, {"cycles", double(G.cycles())},
                       {"murray", std::isfinite(murray) ? Json(round3(murray)) : Json()}});
}

inline Json agreement(const mr::Mask& found, const mr::Mask& reference) {
  const auto cl = eikon::fusion::cldice(found, reference);
  const auto sd = eikon::fusion::surface_distances(found, reference);
  return Json::object({{"dice", round3(eikon::fusion::dice(found, reference))}, {"cldice", round3(cl.cldice)},
                       {"precision", round3(cl.precision)}, {"sensitivity", round3(cl.sensitivity)},
                       {"mean_distance", round3(sd.mean)}, {"p95_distance", round3(sd.p95)}});
}

inline Json geometry_json(const eg::Geometry& g) {
  Json d = Json::array();
  for (int c = 0; c < 3; c++) d.push(Json::array({round3(g.direction[0][c]), round3(g.direction[1][c]), round3(g.direction[2][c])}));
  return Json::object({{"size", Json::array({double(g.size[0]), double(g.size[1]), double(g.size[2])})},
                       {"spacing", Json::array({round3(g.spacing[0]), round3(g.spacing[1]), round3(g.spacing[2])})},
                       {"origin", Json::array({round1(g.origin[0]), round1(g.origin[1]), round1(g.origin[2])})}, {"direction", d}});
}

inline Json list_cases(Reconstructions& R) {
  Json out = Json::array();
  out.push(Json::object({{"id", "simulated"}, {"kind", "simulated"},
                         {"title", "A simulated body: contrast CT and TOF-MRA, integrated"},
                         {"note", "One body imaged twice by the medical-imaging domain's forward models, the patient moved between scans; arteries and veins told apart by what each modality shows."}}));
  const auto dir = R.products / "msd-hepatic-vessels";
  if (std::filesystem::exists(dir)) {
    std::vector<std::string> ids;
    for (const auto& e : std::filesystem::directory_iterator(dir)) {
      const auto n = e.path().filename().string();
      if (n.size() > 10 && n.ends_with("-image.nii")) ids.push_back(n.substr(0, n.size() - 10));
    }
    std::sort(ids.begin(), ids.end());
    for (const auto& id : ids)
      out.push(Json::object({{"id", id}, {"kind", "ct"}, {"title", "Hepatic vessels, CT case " + id.substr(id.find('_') + 1)},
                             {"note", "Portal-venous CT from the Medical Segmentation Decathlon (Task 08, CC BY-SA 4.0), compared with the experts' vessel labels."}}));
  }
  return out;
}

namespace detail {
// The CT of a Decathlon case: read, carried to LPS, reconstructed within the annotated territory.
inline void run_ct(Reconstructions& R, Case& c) {
  const auto dir = R.products / "msd-hepatic-vessels";
  auto image = eg::read_nifti(dir / (c.id + "-image.nii"));
  auto label = eg::read_nifti(dir / (c.id + "-label.nii"));
  image.geometry = medical_imaging::dicom::to_lps(image.geometry);
  label.geometry = image.geometry;
  mr::Mask vessel(image.geometry), annotated(image.geometry), tumour(image.geometry);
  for (std::size_t v = 0; v < vessel.data.size(); v++) {
    vessel.data[v] = label.data[v] > 0.5f && label.data[v] < 1.5f;
    tumour.data[v] = label.data[v] > 1.5f;
    annotated.data[v] = label.data[v] > 0.5f;
  }
  auto region = eikon::distance::dilate(annotated, 10.0);
  for (std::size_t v = 0; v < region.data.size(); v++) region.data[v] = region.data[v] && !tumour.data[v];
  const auto rec = mr::reconstruct(image, mr::portal_venous_ct(1.0, 6.0), &region);
  auto found = mr::embed(rec.mask, image.geometry);
  for (std::size_t v = 0; v < found.data.size(); v++) found.data[v] = found.data[v] && region.data[v];
  c.image = std::move(image), c.reference = vessel, c.found = found, c.region = region;
  c.result = Json::object({{"id", c.id}, {"kind", "ct"}, {"geometry", geometry_json(c.image.geometry)},
                           {"surfaces", Json::object({{"found", surface_json(found)}, {"reference", surface_json(vessel)}})},
                           {"tree", tree_json(found)}, {"agreement", agreement(found, vessel)},
                           {"values", Json::object({{"background", round1(rec.background)}, {"lumen", round1(rec.lumen)}})},
                           {"window", Json::array({100.0, 400.0})},
                           {"credit", "CT and labels: the Medical Segmentation Decathlon, Task 08 (Memorial Sloan Kettering Cancer Center), CC BY-SA 4.0"}});
}

// The simulated body: an organ with an arterial and a venous tree, imaged by CT and TOF-MRA, integrated.
inline void run_simulated(Case& c) {
  calx::random::Stream rng(20261104, 51, 0, 0);
  eg::Geometry fine;
  fine.size = {200, 172, 144};
  fine.spacing = {0.35, 0.35, 0.35};
  fine.origin = {-34.825, -29.925, -25.025};
  eg::Image<float> organ(fine);
  for (std::size_t v = 0; v < fine.count(); v++) {
    const auto ijk = fine.unlinear(v);
    const auto x = fine.world(double(ijk[0]), double(ijk[1]), double(ijk[2]));
    organ.data[v] = x[0] * x[0] / 900 + x[1] * x[1] / 625 + x[2] * x[2] / 484 <= 1 ? 1.0f : 0.0f;
  }
  ep::TreeOptions to;
  to.depth = 2, to.angle = 0.45, to.length_ratio = 0.9;
  const auto arteries = ep::tree({-11, 0, -19}, {-0.1, 0.05, 1}, 13, 2.2, to, rng);
  const auto veins = ep::tree({11, 2, 19}, {0.1, -0.05, -1}, 13, 2.6, to, rng);
  const auto art = ep::rasterise(fine, ep::shapes(arteries), 2), ven = ep::rasterise(fine, ep::shapes(veins), 2);
  mm::CtProtocol cp;
  cp.grid.size = {100, 86, 50};
  cp.grid.spacing = {0.7, 0.7, 1.0};
  cp.grid.origin = {-34.65, -29.75, -24.5};
  cp.background = -100;
  const auto ct = mm::ct({{organ, 60}, {ven.occupancy, 180}, {art.occupancy, 250}}, cp, rng);
  mr::Mask ct_region(cp.grid);
  const auto oc = eg::resample(organ, cp.grid);
  for (std::size_t v = 0; v < oc.data.size(); v++) ct_region.data[v] = oc.data[v] >= 0.5f;
  const auto T = eikon::registration::about({5 * std::numbers::pi / 180, -3 * std::numbers::pi / 180, 6 * std::numbers::pi / 180}, {2, -1.5, 3}, {0, 0, 0});
  const auto Ti = T.inverse();
  mm::TofProtocol tp;
  tp.grid.size = {110, 96, 84};
  tp.grid.spacing = {0.6, 0.6, 0.6};
  tp.grid.origin = {-32.7, -28.5, -24.9};
  tp.slab_entry = {0, 0, -25};
  tp.noise = 0.008;
  const auto in_mra = [&](const eg::Image<float>& im) { return eg::resample(im, tp.grid, [&](const eg::Vec3& y) { return Ti(y); }); };
  const auto organ_m = in_mra(organ), art_m = in_mra(art.occupancy), ven_m = in_mra(ven.occupancy);
  eg::Image<float> tissue(tp.grid);
  for (std::size_t v = 0; v < tissue.data.size(); v++) tissue.data[v] = organ_m.data[v] + 0.4f * (1 - organ_m.data[v]);
  const auto mra = mm::tof(tissue, {{art_m, 1.6, 350, false}, {ven_m, 1.6, 120, true}}, tp, rng);
  mr::Mask mra_region(tp.grid);
  for (std::size_t v = 0; v < organ_m.data.size(); v++) mra_region.data[v] = organ_m.data[v] >= 0.5f;
  std::vector<mi::Acquisition> acq{{"contrast CT", ct, mr::contrast_ct(0.8, 3.5), mi::Shows::all, &ct_region},
                                   {"TOF-MRA", mra, mr::tof_mra(0.8, 3.0), mi::Shows::arteries, &mra_region}};
  const auto I = mi::integrate(acq);
  const auto truth = [&](const ep::Raster& r) {
    const auto o = eg::resample(r.occupancy, cp.grid);
    mr::Mask m(cp.grid);
    for (std::size_t v = 0; v < o.data.size(); v++) m.data[v] = o.data[v] >= 0.5f && ct_region.data[v];
    return m;
  };
  const auto ta = truth(art), tv = truth(ven);
  mr::Mask fa(cp.grid), fv(cp.grid), both(cp.grid);
  for (std::size_t v = 0; v < fa.data.size(); v++) {
    fa.data[v] = I.arteries.data[v] && ct_region.data[v], fv.data[v] = I.veins.data[v] && ct_region.data[v];
    both.data[v] = ta.data[v] || tv.data[v];
  }
  double reg = 0;
  for (int k = 0; k < 8; k++) {
    const auto x = cp.grid.world(k & 1 ? 99 : 0, k & 2 ? 85 : 0, k & 4 ? 49 : 0);
    const auto u = I.to_acquisition[1](x), w = T(x);
    reg = std::max(reg, std::sqrt((u[0] - w[0]) * (u[0] - w[0]) + (u[1] - w[1]) * (u[1] - w[1]) + (u[2] - w[2]) * (u[2] - w[2])));
  }
  c.image = ct, c.reference = both, c.found = I.vessels, c.arteries = fa, c.veins = fv, c.region = ct_region;
  c.result = Json::object({{"id", c.id}, {"kind", "simulated"}, {"geometry", geometry_json(cp.grid)},
                           {"surfaces", Json::object({{"arteries", surface_json(fa)}, {"veins", surface_json(fv)}, {"reference", surface_json(both)}})},
                           {"tree", tree_json(I.vessels)}, {"agreement", agreement(I.vessels, both)},
                           {"arteries", agreement(fa, ta)}, {"veins", agreement(fv, tv)},
                           {"registration", Json::object({{"error", round3(reg)}, {"similarity", round3(I.similarity[1])}})},
                           {"window", Json::array({60.0, 500.0})},
                           {"credit", "Simulated by the medical-imaging domain's forward models: no person's data."}});
}
}  // namespace detail

inline Case& ensure(Reconstructions& R, const std::string& id) {
  auto it = R.cases.find(id);
  if (it != R.cases.end()) return it->second;
  Case c;
  c.id = id;
  if (id == "simulated") c.kind = "simulated", detail::run_simulated(c);
  else {
    if (id.find('/') != std::string::npos || id.find('\\') != std::string::npos || id.find("..") != std::string::npos) throw std::invalid_argument("no case " + id);
    if (!std::filesystem::exists(R.products / "msd-hepatic-vessels" / (id + "-image.nii"))) throw std::invalid_argument("no case " + id);
    c.kind = "ct", detail::run_ct(R, c);
  }
  return R.cases.emplace(id, std::move(c)).first->second;
}

// A slice k along the image's third axis: the values windowed to bytes, and a code per pixel: 1 found, 2
// reference, 4 artery, 8 vein, 16 outside the region searched.
inline Json slice(Reconstructions& R, const Json& q) {
  const auto& c = ensure(R, q["case"].text());
  const auto& g = c.image.geometry;
  const double kq = q["k"].as_number();
  const std::size_t k = std::size_t(std::clamp(std::isfinite(kq) ? kq : double(g.size[2]) / 2, 0.0, double(g.size[2] - 1)));
  const double level = c.result["window"][0].as_number(), width = c.result["window"][1].as_number();
  std::vector<std::uint8_t> px(g.size[0] * g.size[1]), code(g.size[0] * g.size[1]);
  for (std::size_t j = 0; j < g.size[1]; j++)
    for (std::size_t i = 0; i < g.size[0]; i++) {
      const std::size_t v = g.linear(i, j, k), p = i + g.size[0] * j;
      px[p] = std::uint8_t(std::clamp((c.image.data[v] - (level - width / 2)) / width, 0.0, 1.0) * 255);
      std::uint8_t b = 0;
      if (c.found.data.size() == g.count() && c.found.data[v]) b |= 1;
      if (c.reference.data.size() == g.count() && c.reference.data[v]) b |= 2;
      if (c.arteries.data.size() == g.count() && c.arteries.data[v]) b |= 4;
      if (c.veins.data.size() == g.count() && c.veins.data[v]) b |= 8;
      if (c.region.data.size() == g.count() && !c.region.data[v]) b |= 16;
      code[p] = b;
    }
  return Json::object({{"case", c.id}, {"k", double(k)}, {"slices", double(g.size[2])}, {"width", double(g.size[0])}, {"height", double(g.size[1])},
                       {"spacing", Json::array({round3(g.spacing[0]), round3(g.spacing[1]), round3(g.spacing[2])})},
                       {"z", round1(g.world(0, 0, double(k))[2])}, {"geometry", geometry_json(g)},
                       {"pixels", base64(px)}, {"codes", base64(code)}});
}

}  // namespace fabrica
