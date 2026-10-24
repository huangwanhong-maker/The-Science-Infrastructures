// Checks of root-to-tip regression: on a tree whose branch lengths are a strict clock's (rate ×
// time) with tips sampled over thirty years, the best root by each objective is the true root,
// the rate and the root's date are exact, and R² is 1; with lengths estimated by maximum
// likelihood from sequences simulated along it, the rate and the root's date are recovered
// within their error.

#include <evolution/phylogeny/dating.hpp>
#include <evolution/phylogeny/likelihood.hpp>
#include <evolution/phylogeny/simulate.hpp>
#include <evolution/phylogeny/splits.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
#include <vector>

using namespace evolution::phylogeny;

int main() {
  calx::verification::Suite s;

  // A Yule tree in years, each tip's branch lengthened by its own sampling delay of 0 to 30
  // years; the root in 1950; the clock 0.004 substitutions per site per year.
  const double rate = 0.004, root_date = 1950;
  Tree timed = scale_tree(yule_tree(30, 51, 1.0 / 8), 1.0);
  calx::random::Stream r(52, 0x44415445, 0, 0);
  std::vector<std::string> labels;
  std::vector<double> dates;
  for (int v : timed.leaves()) timed[v].length += 30 * r.uniform();
  {
    const auto dep = depths(timed);
    for (int v : timed.leaves()) {
      labels.push_back(timed[v].label);
      dates.push_back(root_date + dep[static_cast<std::size_t>(v)]);
    }
  }
  const Tree clock = scale_tree(timed, rate);
  const Tree unrooted = unroot(clock);
  {
    for (ClockObjective o : {ClockObjective::residual_mean_square, ClockObjective::r_squared, ClockObjective::correlation}) {
      const ClockFit f = root_to_tip(unrooted, labels, dates, o);
      const std::string name = o == ClockObjective::residual_mean_square ? "residual mean square" : o == ClockObjective::r_squared ? "R²" : "correlation";
      s.check("on a strict clock the best root (" + name + ") gives the rate exactly", f.rate / rate - 1, 1e-9);
      s.check("and the root's date, 1950 (" + name + ")", f.root_date - root_date, 1e-6, "years");
      s.check("with R² = 1 (" + name + ")", 1 - f.r_squared, 1e-12);
      s.require("its root splits the leaves as the true root does (" + name + ")", [&] {
        // The two sides of the root, as sets of labels.
        const auto side = [](const Tree& t) {
          std::vector<std::string> a;
          std::vector<int> stack{t[t.root].children[0]};
          while (!stack.empty()) {
            const int v = stack.back();
            stack.pop_back();
            if (t.is_leaf(v)) a.push_back(t[v].label);
            for (int c : t[v].children) stack.push_back(c);
          }
          std::sort(a.begin(), a.end());
          return a;
        };
        auto a = side(f.tree), b = side(clock);
        std::vector<std::string> all = labels;
        std::sort(all.begin(), all.end());
        std::vector<std::string> nb;
        std::set_difference(all.begin(), all.end(), b.begin(), b.end(), std::back_inserter(nb));
        return a == b || a == nb;
      }());
    }
  }
  {
    // Sequences along the clock tree; lengths by maximum likelihood on the true topology.
    const Patterns p = compress(simulate_alignment(clock, hky85(5, {0.3, 0.2, 0.2, 0.3}), Rates{}, 20000, 53));
    Likelihood lik(p, hky85(5, base_frequencies(p)), Rates{});
    Tree fit = unrooted;
    lik.optimise_branch_lengths(fit);
    const ClockFit f = root_to_tip(fit, labels, dates);
    s.check("from 20,000 simulated sites the rate is recovered within 5%", f.rate / rate - 1, 0.05);
    s.check("and the root's date within 5 years", f.root_date - root_date, 5, "years");
    s.require("R² stays high (above 0.9)", f.r_squared > 0.9, std::to_string(f.r_squared));
    const ClockFit at = root_to_tip_at_root(reroot(fit, f.branch, f.position), labels, dates);
    s.check("the fit at the root it chose is the fit it reports", at.rate - f.rate, 1e-12);
  }
  return s.done();
}
