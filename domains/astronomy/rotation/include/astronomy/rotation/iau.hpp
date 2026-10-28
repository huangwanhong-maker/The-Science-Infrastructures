// The rotation of the Sun, the Moon and the planets by the IAU's models: the direction of each
// body's north pole and the angle of its prime meridian.
//
// The IAU Working Group on Cartographic Coordinates and Rotational Elements gives, in the ICRF,
// the right ascension α₀ and declination δ₀ of a body's north pole and the angle W of its prime
// meridian, measured eastward along the body's equator from the ascending node of that equator
// on the ICRF equator; with d days and T Julian centuries of TDB since J2000.0 (Archinal et al.
// 2018, Celest. Mech. Dyn. Astron. 130:22, Tables 1–3). A vector in the ICRF becomes body-fixed
// by the frame rotations
//   v_body = R_z(W) R_x(90° − δ₀) R_z(90° + α₀) v_ICRF,                                       (1)
// with x toward the prime meridian and z toward the north pole. The models are those of the
// astronomy domain's JavaScript (rotation.js), which this replaces for the C++ and against which
// it is checked. For the Moon, the IAU model approximates the mean-Earth/polar-axis frame of the
// ephemeris: the checks find its sub-Earth points within 0.0023° (70 m on the surface) of DE440's
// at the dates of the fixtures; moon.hpp gives the frames of DE440 themselves.

#pragma once

#include <array>
#include <cmath>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include <physis/geometry/kinematics.hpp>
#include <physis/geometry/rotation.hpp>

namespace astronomy::rotation {

using physis::geometry::Axis;
using physis::geometry::frame_rotation;
using physis::geometry::Mat3;
using namespace physis::geometry::operators;

struct PoleAndMeridian {
  double alpha0, delta0, w;  // degrees; w in [0, 360)
  double w_rate;             // degrees per day (the rotation's mean rate, for velocities)
};

namespace detail {
inline double sind(double x) { return std::sin(x * std::numbers::pi / 180); }
inline double cosd(double x) { return std::cos(x * std::numbers::pi / 180); }
}  // namespace detail

// The pole and prime meridian of a body (sun, mercury, venus, earth, moon, mars, jupiter,
// saturn, uranus, neptune, pluto) at a TDB Julian date; nullopt for a body without a model here.
inline std::optional<PoleAndMeridian> pole_and_meridian(std::string_view body, double jd_tdb) {
  using detail::cosd;
  using detail::sind;
  const double d = jd_tdb - 2451545.0, T = d / 36525.0;
  double a = 0, dl = 0, w = 0, rate = 0;
  if (body == "sun") {
    a = 286.13, dl = 63.87, w = 84.176 + 14.1844000 * d, rate = 14.1844000;
  } else if (body == "mercury") {
    const double M1 = 174.7910857 + 4.092335 * d, M2 = 349.5821714 + 8.184670 * d, M3 = 164.3732571 + 12.277005 * d,
                 M4 = 339.1643429 + 16.369340 * d, M5 = 153.9554286 + 20.461675 * d;
    a = 281.0103 - 0.0328 * T, dl = 61.4155 - 0.0049 * T;
    w = 329.5988 + 6.1385108 * d + 0.01067257 * sind(M1) - 0.00112309 * sind(M2) - 0.00011040 * sind(M3) - 0.00002539 * sind(M4) -
        0.00000571 * sind(M5);
    rate = 6.1385108;
  } else if (body == "venus") {
    a = 272.76, dl = 67.16, w = 160.20 - 1.4813688 * d, rate = -1.4813688;
  } else if (body == "earth") {
    a = 0 - 0.641 * T, dl = 90 - 0.557 * T, w = 190.147 + 360.9856235 * d, rate = 360.9856235;
  } else if (body == "moon") {
    const double E[13] = {125.045 - 0.0529921 * d, 250.089 - 0.1059842 * d, 260.008 + 13.0120009 * d, 176.625 + 13.3407154 * d,
                          357.529 + 0.9856003 * d, 311.589 + 26.4057084 * d, 134.963 + 13.0649930 * d, 276.617 + 0.3287146 * d,
                          34.226 + 1.7484877 * d,  15.134 - 0.1589763 * d,  119.743 + 0.0036096 * d, 239.961 + 0.1643573 * d,
                          25.053 + 12.9590088 * d};
    auto s = [&](int i) { return sind(E[i - 1]); };
    auto c = [&](int i) { return cosd(E[i - 1]); };
    a = 269.9949 + 0.0031 * T - 3.8787 * s(1) - 0.1204 * s(2) + 0.0700 * s(3) - 0.0172 * s(4) + 0.0072 * s(6) - 0.0052 * s(10) + 0.0043 * s(13);
    dl = 66.5392 + 0.0130 * T + 1.5419 * c(1) + 0.0239 * c(2) - 0.0278 * c(3) + 0.0068 * c(4) - 0.0029 * c(6) + 0.0009 * c(7) + 0.0008 * c(10) -
         0.0009 * c(13);
    w = 38.3213 + 13.17635815 * d - 1.4e-12 * d * d + 3.5610 * s(1) + 0.1208 * s(2) - 0.0642 * s(3) + 0.0158 * s(4) + 0.0252 * s(5) - 0.0066 * s(6) -
        0.0047 * s(7) - 0.0046 * s(8) + 0.0028 * s(9) + 0.0052 * s(10) + 0.0040 * s(11) + 0.0019 * s(12) - 0.0044 * s(13);
    rate = 13.17635815;
  } else if (body == "mars") {
    a = 317.269202 - 0.10927547 * T + 0.000068 * sind(198.991226 + 19139.4819985 * T) + 0.000238 * sind(226.292679 + 38280.8511281 * T) +
        0.000052 * sind(249.663391 + 57420.7251593 * T) + 0.000009 * sind(266.183510 + 76560.6367950 * T) + 0.419057 * sind(79.398797 + 0.5042615 * T);
    dl = 54.432516 - 0.05827105 * T + 0.000051 * cosd(122.433576 + 19139.9407476 * T) + 0.000141 * cosd(43.058401 + 38280.8753272 * T) +
         0.000031 * cosd(57.663379 + 57420.7517205 * T) + 0.000005 * cosd(79.476401 + 76560.6495004 * T) + 1.591274 * cosd(166.325722 + 0.5042615 * T);
    w = 176.049863 + 350.891982443297 * d + 0.000145 * sind(129.071773 + 19140.0328244 * T) + 0.000157 * sind(36.352167 + 38281.0473591 * T) +
        0.000040 * sind(56.668646 + 57420.9295360 * T) + 0.000001 * sind(67.364003 + 76560.2552215 * T) + 0.000001 * sind(104.792680 + 95700.4387578 * T) +
        0.584542 * sind(95.391654 + 0.5042615 * T);
    rate = 350.891982443297;
  } else if (body == "jupiter") {
    const double J[5] = {99.360714 + 4850.4046 * T, 175.895369 + 1191.9605 * T, 300.323162 + 262.5475 * T, 114.012305 + 6070.2476 * T,
                         49.511251 + 64.3000 * T};
    a = 268.056595 - 0.006499 * T + 0.000117 * sind(J[0]) + 0.000938 * sind(J[1]) + 0.001432 * sind(J[2]) + 0.000030 * sind(J[3]) + 0.002150 * sind(J[4]);
    dl = 64.495303 + 0.002413 * T + 0.000050 * cosd(J[0]) + 0.000404 * cosd(J[1]) + 0.000617 * cosd(J[2]) - 0.000013 * cosd(J[3]) + 0.000926 * cosd(J[4]);
    w = 284.95 + 870.5360000 * d, rate = 870.5360000;
  } else if (body == "saturn") {
    a = 40.589 - 0.036 * T, dl = 83.537 - 0.004 * T, w = 38.90 + 810.7939024 * d, rate = 810.7939024;
  } else if (body == "uranus") {
    a = 257.311, dl = -15.175, w = 203.81 - 501.1600928 * d, rate = -501.1600928;
  } else if (body == "neptune") {
    const double N = 357.85 + 52.316 * T;
    a = 299.36 + 0.70 * sind(N), dl = 43.46 - 0.51 * cosd(N), w = 249.978 + 541.1397757 * d - 0.48 * sind(N), rate = 541.1397757;
  } else if (body == "pluto") {
    a = 132.993, dl = -6.163, w = 302.695 + 56.3625225 * d, rate = 56.3625225;
  } else {
    return std::nullopt;
  }
  w = std::fmod(w, 360.0);
  if (w < 0) w += 360;
  return PoleAndMeridian{a, dl, w, rate};
}

// The change of coordinates (1) from the ICRF to a body's IAU frame, with the derivative of the
// rotation about the pole (the slow motion of the pole is neglected in the derivative).
inline physis::geometry::StateTransform icrf_to_iau_frame(std::string_view body, double jd_tdb) {
  const auto p = pole_and_meridian(body, jd_tdb);
  if (!p) throw std::invalid_argument("astronomy: no IAU rotation model for " + std::string(body));
  const double deg = std::numbers::pi / 180;
  const Mat3 pole = frame_rotation(Axis::x, (90 - p->delta0) * deg) * frame_rotation(Axis::z, (90 + p->alpha0) * deg);
  const Mat3 spin = frame_rotation(Axis::z, p->w * deg);
  const double rate = p->w_rate * deg / 86400.0;  // rad/s
  const Mat3 dspin = (-rate) * (physis::geometry::skew({0, 0, 1}) * spin);
  return {spin * pole, dspin * pole, {0, 0, 0}, {0, 0, 0}};
}

// Bodies whose longitudes the IAU counts positive to the west (prograde rotators other than the
// Sun, the Earth and the Moon, by its convention).
inline bool west_positive(std::string_view body) {
  return body == "mercury" || body == "mars" || body == "jupiter" || body == "saturn" || body == "neptune";
}

}  // namespace astronomy::rotation
