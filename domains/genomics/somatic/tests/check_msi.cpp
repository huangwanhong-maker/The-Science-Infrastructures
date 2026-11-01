// Checks of MSIsensor's pieces: the microsatellites of random sequences against a brute-force scan
// for maximal primitive runs; the copies read between the flanks; Pearson's χ² of 2 × K tables of
// lengths against scipy's chi2_contingency (the fixture somatic-chi-square.txt, made by
// make_somatic_fixtures.py); MSIsensor's rank rule against a hand example, where it rejects fewer
// sites than Benjamini and Hochberg's step-up.

#include <genomics/somatic/msi.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics;

namespace {

std::map<std::size_t, std::size_t> parse(const std::string& s) {
  std::map<std::size_t, std::size_t> m;
  std::stringstream in(s);
  std::string item;
  while (std::getline(in, item, ',')) {
    const auto colon = item.find(':');
    m[std::stoul(item.substr(0, colon))] = std::stoul(item.substr(colon + 1));
  }
  return m;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the scan --------------------------------------------------------------------------------------
  {
    // Random sequence with planted runs.
    sequence::Codes x = simulate::random_sequence(20000, 0.5, 3);
    const auto plant = [&](std::size_t at, const std::string& unit, std::size_t copies) {
      for (std::size_t k = 0; k < copies * unit.size(); k++) x[at + k] = sequence::dna_code(unit[k % unit.size()]);
    };
    plant(1000, "A", 14);
    plant(3000, "CA", 9);
    plant(5000, "GAT", 7);
    plant(7000, "ACGT", 6);
    plant(9000, "AATGC", 5);
    plant(11000, "T", 9);      // too short a homopolymer
    plant(13000, "AT", 4);     // too few copies
    for (std::size_t p : {10999u, 11009u, 12998u, 12999u, 13008u, 13009u}) x[p] = sequence::G;  // walled off from chance extension
    plant(15000, "ATAT", 6);   // not primitive: (AT)^12
    somatic::MicrosatelliteOptions o;
    const auto found = somatic::find_microsatellites({x}, o);
    // Brute force: every maximal run of a primitive unit of 1–5 bases within the bounds.
    std::size_t brute = 0;
    for (std::size_t u = 1; u <= 5; u++)
      for (std::size_t a = 0; a + u <= x.size(); a++) {
        if (a >= u && x[a - 1] == x[a + u - 1]) continue;  // not maximal on the left
        std::size_t e = a + u;
        while (e < x.size() && x[e] == x[e - u]) e++;
        const std::size_t copies = (e - a) / u;
        bool prim = true;
        for (std::size_t d = 1; d < u && prim; d++) {
          if (u % d) continue;
          bool p = true;
          for (std::size_t i = d; i < u; i++) p = p && x[a + i] == x[a + i % d];
          if (p) prim = false;
        }
        const std::size_t lo = u == 1 ? 10 : 5, hi = u == 1 ? 50 : 40;
        if (prim && copies >= lo && copies <= hi && a >= 5 && a + copies * u + 5 <= x.size()) brute++;
      }
    s.require("the scan finds the maximal primitive runs in bounds (a brute-force count)", found.size() == brute, std::to_string(found.size()) + " of " + std::to_string(brute));
    bool planted = true;
    const std::vector<std::pair<std::size_t, std::size_t>> want{{1000, 14}, {3000, 9}, {5000, 7}, {7000, 6}, {9000, 5}, {15000, 12}};
    for (auto [at, copies] : want) {
      bool hit = false;
      // A planted run may extend by chance into the random bases about it.
      for (const auto& m : found) hit = hit || (m.begin <= at && m.begin + 5 >= at && m.copies >= copies);
      planted = planted && hit;
    }
    bool none = true;
    for (const auto& m : found) none = none && m.begin != 11000 && m.begin != 13000;
    s.require("the planted sites are found ((ATAT)^6 as (AT)^12), the short ones not", planted && none);
    // Copies read between the flanks.
    const auto& m = found.front();
    sequence::Codes read(m.left);
    for (int k = 0; k < 3; k++) read.insert(read.end(), m.motif.begin(), m.motif.end());
    sequence::Codes withflank = read;
    withflank.insert(withflank.end(), m.right.begin(), m.right.end());
    s.require("copies between the flanks: 3", somatic::copies_in_read(withflank, m) == 3);
    s.require("no right flank: no reading", somatic::copies_in_read(read, m) == -1);
  }

  // ---- χ² against scipy -------------------------------------------------------------------------------
  {
    std::ifstream f(std::string(GENOMICS_FIXTURES) + "/somatic-chi-square.txt");
    std::string line;
    std::size_t n = 0;
    double worst_stat = 0, worst_p = 0;
    bool dof = true;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::stringstream in(line);
      std::string a, b, stat, d, p;
      std::getline(in, a, '\t');
      std::getline(in, b, '\t');
      std::getline(in, stat, '\t');
      std::getline(in, d, '\t');
      std::getline(in, p, '\t');
      const auto r = somatic::chi_square_two(parse(a), parse(b));
      worst_stat = std::max(worst_stat, calx::verification::relative(r.statistic, std::stod(stat)));
      worst_p = std::max(worst_p, calx::verification::relative(r.p, std::stod(p)));
      dof = dof && r.dof == std::stoul(d);
      n++;
    }
    if (n == 0) s.skip("χ² against scipy: the fixture somatic-chi-square.txt is missing");
    else {
      s.check("χ² statistics against scipy 1.18.1 (" + std::to_string(n) + " tables, relative)", worst_stat, 1e-12);
      s.check("χ² p-values against scipy (relative)", worst_p, 1e-9);
      s.require("degrees of freedom: the lengths seen, less one", dof);
    }
    const auto one = somatic::chi_square_two({{10, 30}}, {{10, 25}});
    s.require("one length seen: p = 1", one.p == 1 && one.dof == 0);
  }

  // ---- the rank rule --------------------------------------------------------------------------------
  {
    std::vector<somatic::MsiSite> sites(3);
    const double p[3] = {0.045, 0.01, 0.04};
    for (int i = 0; i < 3; i++) {
      sites[static_cast<std::size_t>(i)].tested = true;
      sites[static_cast<std::size_t>(i)].test.p = p[i];
    }
    somatic::msisensor_fdr(sites, 0.05);
    // Ranks: 0.01 → 0.03 ≤ 0.05; 0.04 → 0.06 > 0.05; 0.045 → 0.045 ≤ 0.05. BH would reject all three.
    s.require("MSIsensor's rule: ranks 1 and 3 somatic, rank 2 not (BH's step-up would take all three)",
              sites[1].somatic && sites[0].somatic && !sites[2].somatic && std::fabs(sites[2].fdr - 0.06) < 1e-12);
  }
  return s.done();
}
