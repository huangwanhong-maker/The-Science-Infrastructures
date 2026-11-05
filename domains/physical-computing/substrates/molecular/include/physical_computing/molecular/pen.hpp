// Enzymatic DNA networks: the PEN toolbox's predator and prey, and its travelling waves (Fujii and Rondelez 2013;
// Padirac, Fujii, Estévez-Torres and Rondelez 2013) (substrates/molecular/README.md, simulator 3). Provisional.
//
// The toolbox (Montagne, Plasson, Sakai, Fujii and Rondelez 2011): short DNA strands are made on templates by a
// polymerase and a nicking enzyme, and destroyed by an exonuclease. In the molecular ecosystem of Fujii and
// Rondelez (2013, ACS Nano 7, 27), a prey strand N replicates on its template, its growth limited by the template and
// the enzymes (logistic in N with capacity K, set by the template's concentration); a predator strand P replicates by
// using the prey as its template, the polymerase saturating in N (an appetite of Holling's type II); the exonuclease
// destroys the predator at a rate m. The system has Rosenzweig and MacArthur's form (physis::dynamics), whose Hopf
// point K_H = 2x* + 1/(a h) is in closed form: more prey template, and the ecosystem oscillates, as they saw.
//
// Waves: the prey alone, diffusing at D, invades as a pulled front at 2√(D r) less Bramson's correction (Padirac
// et al. 2013, Proc. Natl. Acad. Sci. 110, 7537 saw such fronts, and predator–prey waves).
//
// Limits (declared): scaled parameters, not fitted to the papers' recipes; the enzymes' shared saturation and the
// prey's degradation folded into its logistic term; no binding kinetics. Checks: tests/check_pen.cpp.
#pragma once

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include <physis/chemistry/field_kinetics.hpp>
#include <physis/dynamics/flows.hpp>

#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::molecular {

struct PenPredatorPrey {
  // growth r, capacity K (the prey's template), predation a, handling h, the predator's yield e and destruction m.
  physis::dynamics::RosenzweigMacArthur law{1.0, 2.0, 1.0, 1.0, 0.5, 0.2};

  void operator()(double t, const physis::dynamics::State& x, physis::dynamics::State& dx) const { law(t, x, dx); }
  double hopf_template() const { return law.hopf_capacity(); }

  // The prey alone on a chain of n sites of spacing dx, diffusing at D: its front.
  physis::chemistry::FieldKinetics prey_front(std::size_t n, double spacing, double diffusion) const {
    std::vector<std::tuple<std::size_t, std::size_t, double>> edges;
    for (std::size_t i = 0; i + 1 < n; i++) edges.push_back({i, i + 1, 1 / (spacing * spacing)});
    const double r = law.r, K = law.K;
    physis::chemistry::FieldKinetics F(1, [r, K](const double* c, double* dc) { dc[0] = r * c[0] * (1 - c[0] / K); }, physis::chemistry::SparseLaplacian::graph(n, edges), {diffusion});
    return F;
  }
};

inline runtime::TargetDescription describe_pen() {
  runtime::TargetDescription d;
  d.family = "molecular", d.substrate = "enzymatic DNA networks (PEN)";
  d.accepts = {tasks::Kind::phenomenon, tasks::Kind::flow};
  d.inputs = {{"template concentrations", "scaled"}};
  d.outputs = {{"prey and predator strands", "scaled"}};
  d.materials.materials = {{"DNA", "templates and signals"}, {"enzymes (polymerase, nickase, exonuclease)", "production and destruction"}, {"dNTPs", "the fuel"}};
  d.materials.reuse = computer::Reuse::consumed_per_run;
  d.materials.end_of_life = computer::EndOfLife::biodegradable;
  d.model = "the PEN predator–prey ecosystem (Fujii and Rondelez 2013) in Rosenzweig and MacArthur's form (physis::dynamics); fronts (Padirac et al. 2013)";
  d.limits = {"scaled, illustrative parameters", "the enzymes' saturation folded into the prey's logistic term", "no binding kinetics"};
  return d;
}

}  // namespace physical_computing::molecular
