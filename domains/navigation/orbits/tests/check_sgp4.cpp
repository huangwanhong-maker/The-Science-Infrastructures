// Checks of the navigation bridge's SGP4/SDP4, its TLE reader and the TEME frame (navigation/orbits).
//
// References:
// - Vallado's verification set and his C++'s output on it (tests/fixtures/SGP4-VER.TLE and
//   tcppver.out; Vallado, Crawford, Hujsak and Kelso 2006, AIAA 2006-6753 rev. 2): every listed
//   time of every satellite, the deep-space, resonant, decaying and failing cases among them, run
//   as python-sgp4's own check runs them (gravity WGS 72, the improved operation mode);
// - Vallado's C++ again (python-sgp4's compiled copy) in AFSPC's mode with WGS 72 and the improved
//   mode with WGS 84 (tests/fixtures/sgp4-reference.txt, tests/make_sgp4_fixtures.py);
// - the TEME → ITRF example of AIAA 2006-6753 rev. 2, Appendix C, as published (and as skyfield's
//   tests hold it), and as skyfield computes it at the instant exactly (tests/fixtures/
//   teme-reference.txt);
// - a two-line element set of the ISS (Wikipedia, "Two-line element set", the example of 2008
//   day 264) for the reader's fields and checksums; Alpha-5 catalogue numbers by their definition.

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include <calx/verification/check.hpp>
#include <navigation/orbits/sgp4.hpp>
#include <navigation/orbits/teme.hpp>
#include <navigation/orbits/tle.hpp>

using namespace navigation::orbits;
using calx::verification::Suite;
using metrology::time::JulianDate;

namespace {

std::vector<std::string> read_lines(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  std::vector<std::string> out;
  std::string line;
  while (std::getline(in, line)) {
    while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
    out.push_back(line);
  }
  return out;
}

struct VerificationSet {
  std::string line1, line2;
  double start = 0, stop = 0, step = 0;
};

struct Row {
  double t;
  std::array<double, 6> state;  // km, km/s
};

struct Block {
  long number = 0;
  std::vector<Row> rows;
};

std::vector<double> numbers(const std::string& s) {
  std::istringstream in(s);
  std::vector<double> v;
  double x;
  while (in >> x) v.push_back(x);
  return v;
}

double km_to_m(double x) { return x * 1000.0; }

}  // namespace

int main() {
  Suite s;
  const std::string fixtures = NAVIGATION_FIXTURES;
  constexpr double deg = std::numbers::pi / 180.0, arcsec = std::numbers::pi / 648000.0;

  // ---- the TLE reader ------------------------------------------------------------------------------
  {
    const std::string l1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927";
    const std::string l2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537";
    const Tle t = parse_tle(l1, l2, "ISS (ZARYA)");
    s.require("TLE: the ISS's set of 2008 day 264 is read with its name and number", t.name == "ISS (ZARYA)" && t.number == 25544 && t.designator == "98067A");
    s.require("TLE: both its checksums (7, 7) verify", t.checksums_ok() && t.checksum1 == 7 && t.checksum2 == 7);
    s.check("TLE: inclination 51.6416°", t.inclination / deg - 51.6416, 1e-12, "deg");
    s.check("TLE: eccentricity 0.0006703 (implied decimal point)", t.eccentricity - 0.0006703, 1e-16);
    s.check("TLE: mean motion 15.72125391 rev/day", t.mean_motion / (2 * std::numbers::pi / 86400.0) - 15.72125391, 1e-11, "rev/day");
    s.check("TLE: B* = −0.11606 × 10⁻⁴ per Earth radius (mantissa and exponent)", t.bstar + 0.11606e-4, 1e-18);
    s.check("TLE: ṅ/2 = −0.00002182 rev/day²", t.mean_motion_dot_half / (2 * std::numbers::pi / 86400.0 / 86400.0) + 0.00002182, 1e-15);
    s.require("TLE: revolution 56353 and element set 292", t.revolution == 56353 && t.element_number == 292);
    // 2008 day 264 is September 20 (a leap year): JD of its 0h 2454729.5.
    s.check("TLE: the epoch's day as a Julian date (2008 September 20, 0h UTC)", t.epoch.d1 - 2454729.5, 0, "d");
    s.check("TLE: the epoch's fraction as written", t.epoch.d2 - 0.51782528, 1e-17, "d");
    std::string bad = l1;
    bad[68] = '8';
    s.require("TLE: a wrong checksum is reported, not refused", !parse_tle(bad, l2).checksums_ok());
    bool refused = false;
    try {
      parse_tle(l1.substr(0, 40), l2);
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("TLE: a line too short is refused", refused);
    refused = false;
    try {
      parse_tle(l1, "2 25545" + l2.substr(7));
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("TLE: lines of two different satellites are refused", refused);
    s.require("TLE: Alpha-5 numbers (A0000 = 100000, J0001 = 180001, Z9999 = 339999, I and O skipped)",
              alpha5_number("A0000") == 100000 && alpha5_number("J0001") == 180001 && alpha5_number("P1234") == 231234 &&
                  alpha5_number("Z9999") == 339999 && alpha5_number("00005") == 5);
    const std::string text = "ISS (ZARYA)\r\n" + l1 + "\r\n" + l2 + "\r\n# a comment\n" + l1 + "\n" + l2 + "\n";
    const auto many = read_tles(text);
    s.require("TLE: a file of a three-line and a two-line set", many.size() == 2 && many[0].name == "ISS (ZARYA)" && many[1].name.empty());
  }

  // ---- Vallado's verification set ------------------------------------------------------------------------
  std::vector<VerificationSet> sets;
  {
    const auto lines = read_lines(fixtures + "/SGP4-VER.TLE");
    for (std::size_t i = 0; i + 1 < lines.size(); i++) {
      if (lines[i].rfind("1 ", 0) == 0 && lines[i + 1].rfind("2 ", 0) == 0) {
        VerificationSet v{lines[i], lines[i + 1]};
        const auto run = numbers(lines[i + 1].substr(69));
        if (run.size() == 3) v.start = run[0], v.stop = run[1], v.step = run[2];
        sets.push_back(v);
        i++;
      }
    }
  }
  std::vector<Block> blocks;
  {
    for (const auto& line : read_lines(fixtures + "/tcppver.out")) {
      if (line.find(" xx") != std::string::npos) {
        blocks.push_back({std::stol(line), {}});
      } else if (!blocks.empty()) {
        const auto v = numbers(line.substr(0, std::min<std::size_t>(line.size(), 107)));
        if (v.size() >= 7) blocks.back().rows.push_back({v[0], {v[1], v[2], v[3], v[4], v[5], v[6]}});
      }
    }
  }
  s.require("verification: 33 sets read, each with its run", sets.size() == 33 && std::all_of(sets.begin(), sets.end(), [](const auto& v) { return v.step > 0; }));
  s.require("verification: 33 blocks of Vallado's output read (700 lines)", blocks.size() == 33);
  {
    std::size_t lines = blocks.size();
    for (const auto& b : blocks) lines += b.rows.size();
    s.require("verification: Vallado's output has 700 lines", lines == 700, std::to_string(lines));
  }

  // The 00005 epoch: 2000 day 179.78495062 = June 27, 18:50:19.73357 UTC.
  {
    const Tle t = parse_tle(sets[0].line1, sets[0].line2);
    s.check("verification: the TEME example's epoch, 2000 June 27 (JD 2451722.5) + 0.78495062 d", (t.epoch.d1 - 2451722.5) + (t.epoch.d2 - 0.78495062), 1e-15, "d");
  }

  // Run each set as python-sgp4's check runs it, and hold every line to Vallado's.
  {
    double worst_r = 0, worst_v = 0, worst_r_deep = 0, worst_v_deep = 0, worst_r_long = 0, worst_v_long = 0;
    int compared = 0, mismatched_times = 0, missing = 0, extra = 0, wrong_numbers = 0;
    std::vector<int> errors;
    int stale = 0;
    for (std::size_t k = 0; k < sets.size() && k < blocks.size(); k++) {
      const Tle t = parse_tle(sets[k].line1, sets[k].line2);
      if (t.number != blocks[k].number) wrong_numbers++;
      const Sgp4 sat(t);
      std::vector<Row> mine;
      bool stopped = false;
      auto emit = [&](double ts) {
        if (stopped) return;
        std::array<double, 3> r{}, v{};
        const int e = sat.propagate_km(ts, r, v);
        if (e != 0) {
          errors.push_back(e);
          stopped = true;
          return;
        }
        mine.push_back({ts, {r[0], r[1], r[2], v[0], v[1], v[2]}});
      };
      // The epoch's line: when the epoch itself fails, Vallado's program prints the previous
      // satellite's numbers again (python-sgp4: "Use previous data line"), which are not compared.
      std::array<double, 3> r0{}, v0{};
      const int e0 = sat.propagate_km(0.0, r0, v0);
      std::size_t first = 0;
      if (e0 != 0 && e0 != 6) {
        errors.push_back(e0);
        stale++;
        first = 1;
        stopped = true;
      } else {
        mine.push_back({0.0, {r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]}});
      }
      double ts = sets[k].start;
      while (!stopped && ts <= sets[k].stop) {
        if (!(ts == sets[k].start && ts == 0.0)) emit(ts);
        ts += sets[k].step;
      }
      if (!stopped && ts - sets[k].stop < sets[k].step - 1e-6) emit(sets[k].stop);

      const auto& expect = blocks[k].rows;
      if (mine.size() + first < expect.size()) missing += static_cast<int>(expect.size() - first - mine.size());
      if (mine.size() + first > expect.size()) extra += static_cast<int>(mine.size() + first - expect.size());
      for (std::size_t i = first; i < expect.size() && i - first < mine.size(); i++) {
        const Row& a = mine[i - first];
        const Row& e = expect[i];
        if (std::fabs(a.t - e.t) > 1e-6) mismatched_times++;
        double dr = 0, dv = 0;
        for (int c = 0; c < 3; c++) dr = std::max(dr, std::fabs(a.state[c] - e.state[c]));
        for (int c = 3; c < 6; c++) dv = std::max(dv, std::fabs(a.state[c] - e.state[c]));
        compared++;
        // The last case runs 3.5 years from its epoch (1 844 000 min): kept apart, below.
        if (std::fabs(a.t) > 1e5) {
          worst_r_long = std::max(worst_r_long, dr);
          worst_v_long = std::max(worst_v_long, dv);
          continue;
        }
        worst_r = std::max(worst_r, dr);
        worst_v = std::max(worst_v, dv);
        if (sat.deep_space()) {
          worst_r_deep = std::max(worst_r_deep, dr);
          worst_v_deep = std::max(worst_v_deep, dv);
        }
      }
    }
    s.require("verification: each block is the satellite of its set", wrong_numbers == 0);
    s.require("verification: every line of Vallado's output is reproduced, at its time, and no other",
              missing == 0 && extra == 0 && mismatched_times == 0 && compared == 700 - 33 - stale,
              std::to_string(compared) + " compared, " + std::to_string(missing) + " missing, " + std::to_string(extra) + " extra");
    // tcppver.out prints positions to 1e-8 km and velocities to 1e-9 km/s: half of that is its
    // rounding, and every line within 10⁵ min of its epoch agrees to it.
    s.check("verification: positions against Vallado's, every line within 10⁵ min (printed to 10 µm)", km_to_m(worst_r), 5.1e-6, "m");
    s.check("verification: velocities against Vallado's, every line within 10⁵ min (printed to 1 µm/s)", km_to_m(worst_v), 5.1e-7, "m/s");
    s.check("verification: positions of the deep-space and resonant cases among them", km_to_m(worst_r_deep), 5.1e-6, "m");
    s.check("verification: velocities of the deep-space and resonant cases among them", km_to_m(worst_v_deep), 5.1e-7, "m/s");
    // 3.5 years out, rounding accumulated in 1.8 × 10⁶ minutes of secular terms parts Vallado's file
    // from any recompilation of his code: python-sgp4's pure model and its compiled copy of his C++
    // both differ from it by 0.12 mm there, as this does; python-sgp4 holds it to 0.2 mm.
    s.check("verification: the 3.5-year run (1 844 000 to 1 845 100 min), positions, within python-sgp4's 0.2 mm", km_to_m(worst_r_long), 2e-4, "m");
    s.check("verification: the 3.5-year run, velocities (printed to 1 µm/s)", km_to_m(worst_v_long), 5.1e-7, "m/s");
    std::string got;
    for (int e : errors) got += std::to_string(e) + " ";
    s.require("verification: the cases that fail fail as Vallado's: errors 1, 1, 6, 6, 4, 3, 6 in order", errors == std::vector<int>{1, 1, 6, 6, 4, 3, 6}, got);
    s.require("verification: one case fails at its epoch (error 3, a perturbed eccentricity above 1)", stale == 1);
  }

  // ---- AFSPC's mode and WGS 84, against Vallado's C++ ---------------------------------------------------
  {
    const auto lines = read_lines(fixtures + "/sgp4-reference.txt");
    double worst_r = 0, worst_v = 0, worst_r_decay = 0, worst_v_decay = 0;
    int compared = 0, error_mismatch = 0, rows = 0;
    for (const auto& line : lines) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream in(line);
      std::string catalogue, gravity, mode;
      double ts;
      int code;
      in >> catalogue >> gravity >> mode >> ts >> code;
      std::array<double, 6> want{};
      for (double& x : want) {
        std::string w;
        in >> w;
        x = w == "nan" ? NAN : std::stod(w);
      }
      rows++;
      const auto it = std::find_if(sets.begin(), sets.end(), [&](const auto& v) { return v.line1.substr(2, 5) == catalogue; });
      if (it == sets.end()) continue;
      const Tle t = parse_tle(it->line1, it->line2);
      const Sgp4 sat(t, gravity == "wgs84" ? Gravity::wgs84 : Gravity::wgs72, mode == "a" ? OperationMode::afspc : OperationMode::improved);
      std::array<double, 3> r{}, v{};
      const int e = sat.propagate_km(ts, r, v);
      if (e != code) error_mismatch++;
      if (e == 0 || e == 6) {
        // The two cases of decay (28872, 29141) are propagated past it, where the drag polynomial
        // 1 − C₁t − D₂t² − … nears zero and is squared: there the last bit of a library's pow
        // or sin is magnified (python-sgp4's pure model and its compiled C++ part by 4 µm there).
        const bool decay = catalogue == "28872" || catalogue == "29141";
        double& wr = decay ? worst_r_decay : worst_r;
        double& wv = decay ? worst_v_decay : worst_v;
        for (int c = 0; c < 3; c++) wr = std::max(wr, std::fabs(r[c] - want[c]));
        for (int c = 0; c < 3; c++) wv = std::max(wv, std::fabs(v[c] - want[3 + c]));
        compared++;
      }
    }
    s.require("Vallado's C++, AFSPC mode (WGS 72) and improved mode (WGS 84): every set at five times", rows == 330 && compared > 280, std::to_string(compared) + " states of " + std::to_string(rows));
    s.require("… with the same error codes", error_mismatch == 0, std::to_string(error_mismatch) + " differ");
    // Built with another compiler and its maths library, his code and this agree to the micrometres
    // in which python-sgp4's two implementations agree with each other (1.3 µm on these sets).
    s.check("… positions against Vallado's C++ (live orbits, to 3 days)", km_to_m(worst_r), 5e-6, "m");
    s.check("… velocities against Vallado's C++ (live orbits)", km_to_m(worst_v), 1e-8, "m/s");
    s.check("… positions of the two cases propagated past their decay", km_to_m(worst_r_decay), 5e-5, "m");
    s.check("… velocities of the two cases propagated past their decay", km_to_m(worst_v_decay), 1e-7, "m/s");
  }

  // ---- the interface: SI, TEME, the epoch and instants -------------------------------------------------------
  {
    const Tle t = parse_tle(sets[0].line1, sets[0].line2);
    const Sgp4 sat(t);
    std::array<double, 3> r{}, v{};
    sat.propagate_km(4320.0, r, v);
    const Sgp4State st = sat.at(t.epoch.plus_days(3.0));
    s.check("interface: three days after the epoch by date, as by minutes (m)",
            std::max({std::fabs(st.position[0] - 1000 * r[0]), std::fabs(st.position[1] - 1000 * r[1]), std::fabs(st.position[2] - 1000 * r[2])}), 1e-6, "m");
    s.check("interface: the TEME example at 3 days, Vallado's −9060.47373569 km", st.position[0] + 9060473.73569, 5e-6, "m");
    s.require("interface: a decayed satellite answers its state with error 6", [&] {
      const auto it = std::find_if(sets.begin(), sets.end(), [](const auto& v) { return v.line1.substr(2, 5) == "28872"; });
      const Sgp4 decayed(parse_tle(it->line1, it->line2));
      const Sgp4State d = decayed.at_minutes(60.0);
      return d.error == Sgp4Error::decayed && std::isfinite(d.position[0]);
    }());
    s.check("interface: the period of the TEME example from its mean motion, 133.0 min", 2 * std::numbers::pi / (t.mean_motion * 60.0) - 1440.0 / 10.82419157, 1e-9, "min");
  }

  // ---- TEME → ITRS: Vallado et al. 2006, Appendix C ---------------------------------------------------------
  {
    const physis::geometry::Vec3 r_teme{5094.18016210, 6127.64465950, 6380.34453270}, v_teme{-4.746131487, 0.785818041, 5.531931288};
    const double xp = -0.140682 * arcsec, yp = 0.333309 * arcsec, lod = 0.0015563;
    const physis::geometry::Vec3 r_paper{-1033.4793830, 7901.2952754, 6380.3565958}, v_paper{-3.225636520, -2.872451450, 5.531924446};
    // As published: the instant in one double, as Vallado's code and skyfield's test form it
    // (2004 April 6, 07:51:28.386 UTC, UT1 − UTC = −0.439961 s).
    const double jd_single = 2453101.5 + (7 * 3600 + 51 * 60 + 28.386) / 86400.0 + (-0.439961) / 86400.0;
    {
      const auto x = teme_to_itrs({jd_single, 0.0}, xp, yp, lod);
      const auto [r, v] = x.state(r_teme, v_teme);
      double dr = 0, dv = 0;
      for (int c = 0; c < 3; c++) dr = std::max(dr, std::fabs(r[c] - r_paper[c])), dv = std::max(dv, std::fabs(v[c] - v_paper[c]));
      s.check("TEME → ITRF: the published position (AIAA 2006-6753 rev. 2, App. C; printed to 0.1 mm)", km_to_m(dr), 5e-5, "m");
      s.check("TEME → ITRF: the published velocity, with Vallado's ω⊕ and LOD 1.5563 ms (printed to 1 µm/s)", km_to_m(dv), 1e-6, "m/s");
    }
    // At the instant exactly (two parts), against skyfield, whose rate is dθ_GMST82/dt.
    {
      const auto lines = read_lines(fixtures + "/teme-reference.txt");
      JulianDate ut1;
      physis::geometry::Vec3 r_ref{}, v_ref{};
      for (const auto& line : lines) {
        std::istringstream in(line);
        std::string key;
        in >> key;
        if (key == "ut1") in >> ut1.d1 >> ut1.d2;
        if (key == "r") in >> r_ref[0] >> r_ref[1] >> r_ref[2];
        if (key == "v") in >> v_ref[0] >> v_ref[1] >> v_ref[2];
      }
      const auto x = teme_to_itrs(ut1, xp, yp, 0.0, TemeRate::gmst1982);
      const auto [r, v] = x.state(r_teme, v_teme);
      double dr = 0, dv = 0;
      for (int c = 0; c < 3; c++) dr = std::max(dr, std::fabs(r[c] - r_ref[c])), dv = std::max(dv, std::fabs(v[c] - v_ref[c]));
      s.check("TEME → ITRF at the instant exactly, against skyfield (position)", km_to_m(dr), 1e-6, "m");
      s.check("TEME → ITRF at the instant exactly, against skyfield (velocity, its rate)", km_to_m(dv), 1e-9, "m/s");
      // The exact instant moves the published position by the 40 µs a one-double Julian date resolves.
      const auto [rp, vp] = teme_to_itrs(ut1, xp, yp, lod).state(r_teme, v_teme);
      double shift = 0;
      for (int c = 0; c < 3; c++) shift = std::max(shift, std::fabs(rp[c] - r_paper[c]));
      s.check("TEME → ITRF: the exact instant differs from the published one by the one-double date's rounding (< 2 cm)", km_to_m(shift), 0.02, "m");
    }
    // The two rates differ by the general precession in right ascension of IAU 1976, on which GMST
    // 1982 rests: m = 4612.4362″ per Julian century (Lieske et al. 1977, Astron. Astrophys. 58, 1),
    // 7.09 × 10⁻¹² rad/s.
    const JulianDate t{2453101.5, 0.3274};
    s.check("TEME: dθ_GMST82/dt − ω⊕ is IAU 1976's precession in right ascension, 4612.4362″/century",
            (gmst1982_rate(t) - EARTH_ANGULAR_VELOCITY) / (4612.4362 * arcsec / (36525.0 * 86400.0)) - 1, 1e-3);
    s.check("TEME: θ_GMST82 from one double and from two parts agree to the double's rounding", std::remainder(gmst1982(jd_single) - gmst1982(JulianDate{2453101.5, jd_single - 2453101.5}), 2 * std::numbers::pi), 3e-9, "rad");
    // The rotation is orthogonal and its rate is its derivative (central difference, ±1 s).
    const auto a = teme_to_pef(t, 0.0, TemeRate::gmst1982);
    const auto p = teme_to_pef(t.plus_seconds(1.0), 0.0, TemeRate::gmst1982), m = teme_to_pef(t.plus_seconds(-1.0), 0.0, TemeRate::gmst1982);
    double worst = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) worst = std::max(worst, std::fabs((p.r[i][j] - m.r[i][j]) / 2.0 - a.r_dot[i][j]));
    // The central difference's own error is ω³h²/6 = 6.5 × 10⁻¹⁴ s⁻¹ at h = 1 s.
    s.check("TEME: Ṙ of TEME → PEF is the derivative of R (with dθ_GMST82/dt; central difference ±1 s)", worst, 1e-13, "1/s");
    s.check("TEME: θ_GMST82 at J2000.0 UT1 is 280.46061837°", gmst1982(JulianDate{2451545.0, 0.0}) / deg - 280.46061837, 1e-8, "deg");
  }

  return s.done();
}
