// SGP4/SDP4: the analytical propagator of the mean elements of a two-line element set.
//
// SGP4 (simplified general perturbations 4) is the theory with which the United States' space
// surveillance network fits its element sets, and therefore the only one with which a TLE's mean
// elements give positions: Brouwer's first-order zonal theory (J2, J3, J4) in Lyddane's form, with
// a power-density atmosphere (the B* term) for orbits of periods under 225 minutes; SDP4 adds the
// Sun's and the Moon's long-period and secular effects and the resonances of 12-hour and 24-hour
// orbits with the Earth's tesserals, integrated by a fixed-step Euler–Maclaurin scheme (Hoots and
// Roehrich 1980, Spacetrack Report #3; Hoots 1986, Spacetrack Report #6; Hoots, Schumacher and
// Glover 2004, J. Guid. Control Dyn. 27, 174). The code here is Vallado's of 2006 and after, which
// gathered the two models into one, mended them ("sgp4fix") and fixed their behaviour against the
// Space Force's own (Vallado, Crawford, Hujsak and Kelso 2006, "Revisiting Spacetrack Report #3",
// AIAA 2006-6753, rev. 2; companion code of Vallado 2013, Fundamentals of Astrodynamics and
// Applications, 4th ed.: sgp4unit.cpp), as python-sgp4 2.27 transcribes it line by line
// (sgp4/propagation.py), from which this is translated, keeping Vallado's names and his order of
// operations so that the results agree with his to rounding. Section titles below are his.
//
//   python-sgp4: Copyright (c) 2012–2023 Brandon Rhodes. MIT licence: Permission is hereby
//   granted, free of charge, to any person obtaining a copy of this software and associated
//   documentation files (the "Software"), to deal in the Software without restriction, including
//   without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
//   and/or sell copies of the Software, and to permit persons to whom the Software is furnished to
//   do so, subject to the following conditions: The above copyright notice and this permission
//   notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE
//   IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.
//   Vallado's code is distributed by him and by CelesTrak without restriction.
//
// Units. Inside, the model's: distances in Earth radii, times in minutes, angles in radians, as
// Vallado's; at the interface, metres and metres per second in TEME (teme.hpp), the epoch and the
// instants as two-part Julian dates in UTC. Gravity: WGS 72 by default, the constants with which
// the element sets are made and with which Vallado's verification is run; WGS 84 and the older
// WGS 72 (with its xke) are offered as he offers them. Operation mode: 'i' (improved, the default:
// sidereal time by IAU 1982 at the epoch) or 'a' (AFSPC's own, with its sidereal time from 1970
// and its intrinsic-function quirks).
//
// Errors are those of Vallado's code: 1 the mean eccentricity outside [0, 1) (or a < 0.95 Earth
// radii); 2 the mean motion not positive; 3 the perturbed eccentricity outside [0, 1]; 4 the
// semi-latus rectum negative; 5 (no longer used) sub-orbital elements; 6 the satellite decayed (its
// radius below the Earth's). On errors 1–4 there is no state (NaN is answered); on error 6 the
// state is answered with the error, as Vallado answers it.
//
// The deep-space integrator restarts from the epoch at every call (as python-sgp4 does) rather than
// keep its last step, so that a propagation is a pure function of the elements and the time; with
// its fixed steps the two give the same numbers.

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <string>

#include <metrology/time/julian.hpp>
#include <navigation/orbits/teme.hpp>
#include <navigation/orbits/tle.hpp>
#include <physis/geometry/rotation.hpp>

namespace navigation::orbits {

using physis::geometry::Vec3;

enum class Gravity { wgs72old, wgs72, wgs84 };
enum class OperationMode : char { improved = 'i', afspc = 'a' };

// getgravconst: the Earth's constants of each model (Vallado 2006, sgp4unit.cpp).
struct GravityConstants {
  double tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;
  static GravityConstants of(Gravity g) {
    GravityConstants c{};
    if (g == Gravity::wgs72old) {
      c.mu = 398600.79964;  // km³/s²
      c.radiusearthkm = 6378.135;
      c.xke = 0.0743669161;
      c.j2 = 0.001082616;
      c.j3 = -0.00000253881;
      c.j4 = -0.00000165597;
    } else if (g == Gravity::wgs72) {
      c.mu = 398600.8;
      c.radiusearthkm = 6378.135;
      c.xke = 60.0 / std::sqrt(c.radiusearthkm * c.radiusearthkm * c.radiusearthkm / c.mu);
      c.j2 = 0.001082616;
      c.j3 = -0.00000253881;
      c.j4 = -0.00000165597;
    } else {
      c.mu = 398600.5;
      c.radiusearthkm = 6378.137;
      c.xke = 60.0 / std::sqrt(c.radiusearthkm * c.radiusearthkm * c.radiusearthkm / c.mu);
      c.j2 = 0.00108262998905;
      c.j3 = -0.00000253215306;
      c.j4 = -0.00000161098761;
    }
    c.tumin = 1.0 / c.xke;
    c.j3oj2 = c.j3 / c.j2;
    return c;
  }
};

enum class Sgp4Error : int {
  none = 0,
  mean_eccentricity = 1,
  mean_motion = 2,
  perturbed_eccentricity = 3,
  semi_latus_rectum = 4,
  suborbital = 5,
  decayed = 6,
};

inline const char* describe(Sgp4Error e) {
  switch (e) {
    case Sgp4Error::none: return "no error";
    case Sgp4Error::mean_eccentricity: return "the mean eccentricity is outside the range 0.0 to 1.0";
    case Sgp4Error::mean_motion: return "the mean motion is not positive";
    case Sgp4Error::perturbed_eccentricity: return "the perturbed eccentricity is outside the range 0.0 to 1.0";
    case Sgp4Error::semi_latus_rectum: return "the semi-latus rectum is negative";
    case Sgp4Error::suborbital: return "the epoch elements are sub-orbital";
    default: return "the satellite has decayed (its radius is below the Earth's)";
  }
}

// The elements SGP4 starts from, in its own units: what sgp4init takes.
struct Sgp4Elements {
  double epoch = 0;     // days from 1949 December 31, 0h UTC
  double bstar = 0;     // per Earth radius
  double ndot = 0;      // ṅ/2, rad/min² (kept, not used by the theory)
  double nddot = 0;     // n̈/6, rad/min³ (kept, not used)
  double ecco = 0;      // eccentricity
  double argpo = 0;     // argument of perigee, rad
  double inclo = 0;     // inclination, rad
  double mo = 0;        // mean anomaly, rad
  double no_kozai = 0;  // Kozai mean motion, rad/min
  double nodeo = 0;     // right ascension of the ascending node, rad

  static Sgp4Elements of(const Tle& t) {
    Sgp4Elements e;
    // The epoch as Vallado's twoline2rv forms it: the Julian date in one double, less 2433281.5
    // ((jdsatepoch + jdsatepochF) − 2433281.5). It resolves about 20 µs; the lunar-solar terms
    // take their arguments from it, and in the most eccentric orbits (WIND, e = 0.97, whose
    // periodics carry e to 0.99) that rounding moves the position by millimetres, so it is kept
    // to reproduce his numbers. The minutes since the epoch are counted from the two parts.
    e.epoch = (t.epoch.d1 + t.epoch.d2) - 2433281.5;
    e.bstar = t.bstar;
    e.ndot = t.mean_motion_dot_half * 3600.0;
    e.nddot = t.mean_motion_ddot_sixth * 216000.0;
    e.ecco = t.eccentricity;
    e.argpo = t.argument_of_perigee;
    e.inclo = t.inclination;
    e.mo = t.mean_anomaly;
    e.no_kozai = t.mean_motion * 60.0;
    e.nodeo = t.raan;
    return e;
  }
};

// The singly averaged mean elements at an instant (Vallado's "sgp4fix recover singly averaged
// mean elements"): semi-major axis (Earth radii), eccentricity, inclination, node, argument of
// perigee, mean anomaly (rad) and mean motion (rad/min).
struct MeanElements {
  double am = 0, em = 0, im = 0, Om = 0, om = 0, mm = 0, nm = 0;
};

struct Sgp4State {
  Vec3 position{0, 0, 0};  // TEME, m
  Vec3 velocity{0, 0, 0};  // TEME, m/s
  Sgp4Error error = Sgp4Error::none;
  MeanElements mean;
  bool ok() const { return error == Sgp4Error::none; }
};

class Sgp4 {
 public:
  // From a two-line element set; its epoch keeps the TLE's two parts.
  explicit Sgp4(const Tle& t, Gravity g = Gravity::wgs72, OperationMode mode = OperationMode::improved)
      : Sgp4(Sgp4Elements::of(t), t.epoch, g, mode) {
    catalogue_ = t.catalogue;
  }
  // From elements, with the epoch as a two-part Julian date in UTC (which elements.epoch, days from
  // 1950, must name).
  Sgp4(const Sgp4Elements& e, const metrology::time::JulianDate& epoch, Gravity g = Gravity::wgs72,
       OperationMode mode = OperationMode::improved)
      : epoch_(epoch) {
    init(GravityConstants::of(g), static_cast<char>(mode), e);
  }

  // The state at tsince minutes from the epoch.
  Sgp4State at_minutes(double tsince) const {
    Sgp4State s;
    std::array<double, 3> r{}, v{};
    s.error = static_cast<Sgp4Error>(propagate(tsince, r, v, s.mean));
    if (s.error != Sgp4Error::none && s.error != Sgp4Error::decayed) {
      constexpr double nan = std::numeric_limits<double>::quiet_NaN();
      s.position = s.velocity = {nan, nan, nan};
      return s;
    }
    for (int k = 0; k < 3; k++) {
      s.position[k] = r[k] * 1000.0;
      s.velocity[k] = v[k] * 1000.0;
    }
    return s;
  }
  // The state at an instant of UTC, as a two-part Julian date.
  Sgp4State at(const metrology::time::JulianDate& utc) const { return at_minutes(minutes_since_epoch(utc)); }
  double minutes_since_epoch(const metrology::time::JulianDate& utc) const {
    return (utc.d1 - epoch_.d1) * 1440.0 + (utc.d2 - epoch_.d2) * 1440.0;
  }

  // The model's own answer, in its units (km, km/s), as Vallado's sgp4() gives it; the error code.
  int propagate_km(double tsince, std::array<double, 3>& r, std::array<double, 3>& v) const {
    MeanElements m;
    return propagate(tsince, r, v, m);
  }

  const metrology::time::JulianDate& epoch() const { return epoch_; }
  const std::string& catalogue() const { return catalogue_; }
  Sgp4Error init_error() const { return static_cast<Sgp4Error>(init_error_); }
  bool deep_space() const { return r_.method == 'd'; }
  int resonance() const { return r_.irez; }  // 0 none, 1 one day, 2 half a day
  double semi_major_axis() const { return r_.a * r_.radiusearthkm * 1000.0; }  // m, un-Kozai'd mean motion's
  double perigee_height() const { return r_.altp * r_.radiusearthkm * 1000.0; }  // m above the model's radius
  double apogee_height() const { return r_.alta * r_.radiusearthkm * 1000.0; }
  double gst_at_epoch() const { return r_.gsto; }  // rad
  const GravityConstants& gravity() const { return gravity_; }

 private:
  static constexpr double pi = std::numbers::pi;
  static constexpr double twopi = 2.0 * std::numbers::pi;

  // The satellite record (Vallado's elsetrec), fixed after initialisation.
  struct Record {
    double tumin = 0, mu = 0, radiusearthkm = 0, xke = 0, j2 = 0, j3 = 0, j4 = 0, j3oj2 = 0;
    char operationmode = 'i', method = 'n', init = 'y';
    int isimp = 0, irez = 0;
    double bstar = 0, ndot = 0, nddot = 0, ecco = 0, argpo = 0, inclo = 0, mo = 0, no_kozai = 0, nodeo = 0, no_unkozai = 0;
    double a = 0, alta = 0, altp = 0;
    // near earth
    double aycof = 0, con41 = 0, cc1 = 0, cc4 = 0, cc5 = 0, d2 = 0, d3 = 0, d4 = 0, delmo = 0, eta = 0, argpdot = 0, omgcof = 0,
           sinmao = 0, t2cof = 0, t3cof = 0, t4cof = 0, t5cof = 0, x1mth2 = 0, x7thm1 = 0, mdot = 0, nodedot = 0, xlcof = 0,
           xmcof = 0, nodecf = 0;
    // deep space
    double d2201 = 0, d2211 = 0, d3210 = 0, d3222 = 0, d4410 = 0, d4422 = 0, d5220 = 0, d5232 = 0, d5421 = 0, d5433 = 0;
    double dedt = 0, del1 = 0, del2 = 0, del3 = 0, didt = 0, dmdt = 0, dnodt = 0, domdt = 0;
    double e3 = 0, ee2 = 0, peo = 0, pgho = 0, pho = 0, pinco = 0, plo = 0, se2 = 0, se3 = 0, sgh2 = 0, sgh3 = 0, sgh4 = 0,
           sh2 = 0, sh3 = 0, si2 = 0, si3 = 0, sl2 = 0, sl3 = 0, sl4 = 0, gsto = 0, xfact = 0, xgh2 = 0, xgh3 = 0, xgh4 = 0,
           xh2 = 0, xh3 = 0, xi2 = 0, xi3 = 0, xl2 = 0, xl3 = 0, xl4 = 0, xlamo = 0, zmol = 0, zmos = 0, atime = 0, xli = 0,
           xni = 0;
  };

  Record r_;
  GravityConstants gravity_{};
  metrology::time::JulianDate epoch_;
  std::string catalogue_;
  int init_error_ = 0;

  // ---- procedure dpper: deep space long-period periodics -------------------------------------------
  static void dpper(const Record& s, double t, double inclo, char init, double& ep, double& inclp, double& nodep, double& argpp,
                    double& mp, char opsmode) {
    (void)inclo;  // the GSFC choice: the perturbed inclination decides Lyddane's form (below)
    //  ---------------------- constants -----------------------------
    const double zns = 1.19459e-5, zes = 0.01675, znl = 1.5835218e-4, zel = 0.05490;
    //  --------------- calculate time varying periodics -----------
    double zm = s.zmos + zns * t;
    // be sure that the initial call has time set to zero
    if (init == 'y') zm = s.zmos;
    double zf = zm + 2.0 * zes * std::sin(zm);
    double sinzf = std::sin(zf);
    double f2 = 0.5 * sinzf * sinzf - 0.25;
    double f3 = -0.5 * sinzf * std::cos(zf);
    const double ses = s.se2 * f2 + s.se3 * f3;
    const double sis = s.si2 * f2 + s.si3 * f3;
    const double sls = s.sl2 * f2 + s.sl3 * f3 + s.sl4 * sinzf;
    const double sghs = s.sgh2 * f2 + s.sgh3 * f3 + s.sgh4 * sinzf;
    const double shs = s.sh2 * f2 + s.sh3 * f3;
    zm = s.zmol + znl * t;
    if (init == 'y') zm = s.zmol;
    zf = zm + 2.0 * zel * std::sin(zm);
    sinzf = std::sin(zf);
    f2 = 0.5 * sinzf * sinzf - 0.25;
    f3 = -0.5 * sinzf * std::cos(zf);
    const double sel = s.ee2 * f2 + s.e3 * f3;
    const double sil = s.xi2 * f2 + s.xi3 * f3;
    const double sll = s.xl2 * f2 + s.xl3 * f3 + s.xl4 * sinzf;
    const double sghl = s.xgh2 * f2 + s.xgh3 * f3 + s.xgh4 * sinzf;
    const double shll = s.xh2 * f2 + s.xh3 * f3;
    double pe = ses + sel, pinc = sis + sil, pl = sls + sll, pgh = sghs + sghl, ph = shs + shll;

    if (init == 'n') {
      pe = pe - s.peo;
      pinc = pinc - s.pinco;
      pl = pl - s.plo;
      pgh = pgh - s.pgho;
      ph = ph - s.pho;
      inclp = inclp + pinc;
      ep = ep + pe;
      const double sinip = std::sin(inclp), cosip = std::cos(inclp);
      //  ----------------- apply periodics directly ------------
      //  sgp4fix for lyddane choice: the GSFC version, with the perturbed inclination
      if (inclp >= 0.2) {
        ph = ph / sinip;
        pgh = pgh - cosip * ph;
        argpp = argpp + pgh;
        nodep = nodep + ph;
        mp = mp + pl;
      } else {
        //  ---- apply periodics with lyddane modification ----
        const double sinop = std::sin(nodep), cosop = std::cos(nodep);
        double alfdp = sinip * sinop, betdp = sinip * cosop;
        const double dalf = ph * cosop + pinc * cosip * sinop;
        const double dbet = -ph * sinop + pinc * cosip * cosop;
        alfdp = alfdp + dalf;
        betdp = betdp + dbet;
        nodep = std::fmod(nodep, twopi);
        //  sgp4fix for afspc written intrinsic functions: nodep used without a trigonometric function ahead
        if (nodep < 0.0 && opsmode == 'a') nodep = nodep + twopi;
        const double xls = mp + argpp + pl + pgh + (cosip - pinc * sinip) * nodep;
        const double xnoh = nodep;
        nodep = std::atan2(alfdp, betdp);
        if (nodep < 0.0 && opsmode == 'a') nodep = nodep + twopi;
        if (std::fabs(xnoh - nodep) > pi) {
          if (nodep < xnoh) nodep = nodep + twopi;
          else nodep = nodep - twopi;
        }
        mp = mp + pl;
        argpp = xls - mp - cosip * nodep;
      }
    }
  }

  // ---- procedure dscom: deep space common items ----------------------------------------------------
  struct Dscom {
    double snodm, cnodm, sinim, cosim, sinomm, cosomm, day, em, emsq, gam, rtemsq;
    double s1, s2, s3, s4, s5, s6, s7, ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23,
        sz31, sz32, sz33;
    double nm, z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33;
  };
  static Dscom dscom(double epoch, double ep, double argpp, double tc, double inclp, double nodep, double np, Record& s) {
    //  -------------------------- constants -------------------------
    const double zes = 0.01675, zel = 0.05490, c1ss = 2.9864797e-6, c1l = 4.7968065e-7, zsinis = 0.39785416,
                 zcosis = 0.91744867, zcosgs = 0.1945905, zsings = -0.98088458;
    Dscom d{};
    //  --------------------- local variables ------------------------
    d.nm = np;
    d.em = ep;
    d.snodm = std::sin(nodep);
    d.cnodm = std::cos(nodep);
    d.sinomm = std::sin(argpp);
    d.cosomm = std::cos(argpp);
    d.sinim = std::sin(inclp);
    d.cosim = std::cos(inclp);
    d.emsq = d.em * d.em;
    const double betasq = 1.0 - d.emsq;
    d.rtemsq = std::sqrt(betasq);

    //  ----------------- initialize lunar solar terms ---------------
    s.peo = 0.0;
    s.pinco = 0.0;
    s.plo = 0.0;
    s.pgho = 0.0;
    s.pho = 0.0;
    d.day = epoch + 18261.5 + tc / 1440.0;
    const double xnodce = std::fmod(4.5236020 - 9.2422029e-4 * d.day, twopi);
    const double stem = std::sin(xnodce), ctem = std::cos(xnodce);
    const double zcosil = 0.91375164 - 0.03568096 * ctem;
    const double zsinil = std::sqrt(1.0 - zcosil * zcosil);
    const double zsinhl = 0.089683511 * stem / zsinil;
    const double zcoshl = std::sqrt(1.0 - zsinhl * zsinhl);
    d.gam = 5.8351514 + 0.0019443680 * d.day;
    double zx = 0.39785416 * stem / zsinil;
    const double zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    zx = std::atan2(zx, zy);
    zx = d.gam + zx - xnodce;
    const double zcosgl = std::cos(zx), zsingl = std::sin(zx);

    //  ------------------------- do solar terms ---------------------
    double zcosg = zcosgs, zsing = zsings, zcosi = zcosis, zsini = zsinis, zcosh = d.cnodm, zsinh = d.snodm, cc = c1ss;
    const double xnoi = 1.0 / d.nm;

    for (int lsflg = 1; lsflg <= 2; lsflg++) {
      const double a1 = zcosg * zcosh + zsing * zcosi * zsinh;
      const double a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
      const double a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
      const double a8 = zsing * zsini;
      const double a9 = zsing * zsinh + zcosg * zcosi * zcosh;
      const double a10 = zcosg * zsini;
      const double a2 = d.cosim * a7 + d.sinim * a8;
      const double a4 = d.cosim * a9 + d.sinim * a10;
      const double a5 = -d.sinim * a7 + d.cosim * a8;
      const double a6 = -d.sinim * a9 + d.cosim * a10;

      const double x1 = a1 * d.cosomm + a2 * d.sinomm;
      const double x2 = a3 * d.cosomm + a4 * d.sinomm;
      const double x3 = -a1 * d.sinomm + a2 * d.cosomm;
      const double x4 = -a3 * d.sinomm + a4 * d.cosomm;
      const double x5 = a5 * d.sinomm;
      const double x6 = a6 * d.sinomm;
      const double x7 = a5 * d.cosomm;
      const double x8 = a6 * d.cosomm;

      d.z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
      d.z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
      d.z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
      d.z1 = 3.0 * (a1 * a1 + a2 * a2) + d.z31 * d.emsq;
      d.z2 = 6.0 * (a1 * a3 + a2 * a4) + d.z32 * d.emsq;
      d.z3 = 3.0 * (a3 * a3 + a4 * a4) + d.z33 * d.emsq;
      d.z11 = -6.0 * a1 * a5 + d.emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
      d.z12 = -6.0 * (a1 * a6 + a3 * a5) + d.emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
      d.z13 = -6.0 * a3 * a6 + d.emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
      d.z21 = 6.0 * a2 * a5 + d.emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
      d.z22 = 6.0 * (a4 * a5 + a2 * a6) + d.emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
      d.z23 = 6.0 * a4 * a6 + d.emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
      d.z1 = d.z1 + d.z1 + betasq * d.z31;
      d.z2 = d.z2 + d.z2 + betasq * d.z32;
      d.z3 = d.z3 + d.z3 + betasq * d.z33;
      d.s3 = cc * xnoi;
      d.s2 = -0.5 * d.s3 / d.rtemsq;
      d.s4 = d.s3 * d.rtemsq;
      d.s1 = -15.0 * d.em * d.s4;
      d.s5 = x1 * x3 + x2 * x4;
      d.s6 = x2 * x3 + x1 * x4;
      d.s7 = x2 * x4 - x1 * x3;

      //  ----------------------- do lunar terms -------------------
      if (lsflg == 1) {
        d.ss1 = d.s1;
        d.ss2 = d.s2;
        d.ss3 = d.s3;
        d.ss4 = d.s4;
        d.ss5 = d.s5;
        d.ss6 = d.s6;
        d.ss7 = d.s7;
        d.sz1 = d.z1;
        d.sz2 = d.z2;
        d.sz3 = d.z3;
        d.sz11 = d.z11;
        d.sz12 = d.z12;
        d.sz13 = d.z13;
        d.sz21 = d.z21;
        d.sz22 = d.z22;
        d.sz23 = d.z23;
        d.sz31 = d.z31;
        d.sz32 = d.z32;
        d.sz33 = d.z33;
        zcosg = zcosgl;
        zsing = zsingl;
        zcosi = zcosil;
        zsini = zsinil;
        zcosh = zcoshl * d.cnodm + zsinhl * d.snodm;
        zsinh = d.snodm * zcoshl - d.cnodm * zsinhl;
        cc = c1l;
      }
    }

    s.zmol = std::fmod(4.7199672 + 0.22997150 * d.day - d.gam, twopi);
    s.zmos = std::fmod(6.2565837 + 0.017201977 * d.day, twopi);

    //  ------------------------ do solar terms ----------------------
    s.se2 = 2.0 * d.ss1 * d.ss6;
    s.se3 = 2.0 * d.ss1 * d.ss7;
    s.si2 = 2.0 * d.ss2 * d.sz12;
    s.si3 = 2.0 * d.ss2 * (d.sz13 - d.sz11);
    s.sl2 = -2.0 * d.ss3 * d.sz2;
    s.sl3 = -2.0 * d.ss3 * (d.sz3 - d.sz1);
    s.sl4 = -2.0 * d.ss3 * (-21.0 - 9.0 * d.emsq) * zes;
    s.sgh2 = 2.0 * d.ss4 * d.sz32;
    s.sgh3 = 2.0 * d.ss4 * (d.sz33 - d.sz31);
    s.sgh4 = -18.0 * d.ss4 * zes;
    s.sh2 = -2.0 * d.ss2 * d.sz22;
    s.sh3 = -2.0 * d.ss2 * (d.sz23 - d.sz21);

    //  ------------------------ do lunar terms ----------------------
    s.ee2 = 2.0 * d.s1 * d.s6;
    s.e3 = 2.0 * d.s1 * d.s7;
    s.xi2 = 2.0 * d.s2 * d.z12;
    s.xi3 = 2.0 * d.s2 * (d.z13 - d.z11);
    s.xl2 = -2.0 * d.s3 * d.z2;
    s.xl3 = -2.0 * d.s3 * (d.z3 - d.z1);
    s.xl4 = -2.0 * d.s3 * (-21.0 - 9.0 * d.emsq) * zel;
    s.xgh2 = 2.0 * d.s4 * d.z32;
    s.xgh3 = 2.0 * d.s4 * (d.z33 - d.z31);
    s.xgh4 = -18.0 * d.s4 * zel;
    s.xh2 = -2.0 * d.s2 * d.z22;
    s.xh3 = -2.0 * d.s2 * (d.z23 - d.z21);
    return d;
  }

  // ---- procedure dsinit: deep space contributions to mean motion dot from resonance -------------------
  static void dsinit(double xke, double cosim, double emsq, double argpo, double s1, double s2, double s3, double s4, double s5,
                     double sinim, double ss1, double ss2, double ss3, double ss4, double ss5, double sz1, double sz3, double sz11,
                     double sz13, double sz21, double sz23, double sz31, double sz33, double t, double tc, double gsto, double mo,
                     double mdot, double no, double nodeo, double nodedot, double xpidot, double z1, double z3, double z11,
                     double z13, double z21, double z23, double z31, double z33, double ecco, double eccsq, double& em,
                     double& argpm, double& inclm, double& mm, double& nm, double& nodem, Record& s, double& dndt) {
    const double q22 = 1.7891679e-6, q31 = 2.1460748e-6, q33 = 2.2123015e-7, root22 = 1.7891679e-6, root44 = 7.3636953e-9,
                 root54 = 2.1765803e-9, rptim = 4.37526908801129966e-3,  // equates to 7.29211514668855e-5 rad/sec
        root32 = 3.7393792e-7, root52 = 1.1428639e-7, x2o3 = 2.0 / 3.0, znl = 1.5835218e-4, zns = 1.19459e-5;

    //  -------------------- deep space initialization ------------
    s.irez = 0;
    if (nm < 0.0052359877 && nm > 0.0034906585) s.irez = 1;
    if (nm >= 8.26e-3 && nm <= 9.24e-3 && em >= 0.5) s.irez = 2;

    //  ------------------------ do solar terms -------------------
    const double ses = ss1 * zns * ss5;
    const double sis = ss2 * zns * (sz11 + sz13);
    const double sls = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
    const double sghs = ss4 * zns * (sz31 + sz33 - 6.0);
    double shs = -zns * ss2 * (sz21 + sz23);
    //  sgp4fix for 180 deg incl
    if (inclm < 5.2359877e-2 || inclm > pi - 5.2359877e-2) shs = 0.0;
    if (sinim != 0.0) shs = shs / sinim;
    const double sgs = sghs - cosim * shs;

    //  ------------------------- do lunar terms ------------------
    s.dedt = ses + s1 * znl * s5;
    s.didt = sis + s2 * znl * (z11 + z13);
    s.dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
    const double sghl = s4 * znl * (z31 + z33 - 6.0);
    double shll = -znl * s2 * (z21 + z23);
    //  sgp4fix for 180 deg incl
    if (inclm < 5.2359877e-2 || inclm > pi - 5.2359877e-2) shll = 0.0;
    s.domdt = sgs + sghl;
    s.dnodt = shs;
    if (sinim != 0.0) {
      s.domdt = s.domdt - cosim / sinim * shll;
      s.dnodt = s.dnodt + shll / sinim;
    }

    //  ----------- calculate deep space resonance effects --------
    dndt = 0.0;
    const double theta = std::fmod(gsto + tc * rptim, twopi);
    em = em + s.dedt * t;
    inclm = inclm + s.didt * t;
    argpm = argpm + s.domdt * t;
    nodem = nodem + s.dnodt * t;
    mm = mm + s.dmdt * t;
    //   sgp4fix for negative inclinations: the correction stays commented out, as Vallado's

    //  -------------- initialize the resonance terms -------------
    if (s.irez != 0) {
      const double aonv = std::pow(nm / xke, x2o3);

      //  ---------- geopotential resonance for 12 hour orbits ------
      if (s.irez == 2) {
        const double cosisq = cosim * cosim;
        const double emo = em;
        em = ecco;
        const double emsqo = emsq;
        const double emsq2 = eccsq;  // emsq = eccsq within this block, as Vallado's
        const double eoc = em * emsq2;
        const double g201 = -0.306 - (em - 0.64) * 0.440;
        double g211, g310, g322, g410, g422, g520, g533, g521, g532;
        if (em <= 0.65) {
          g211 = 3.616 - 13.2470 * em + 16.2900 * emsq2;
          g310 = -19.302 + 117.3900 * em - 228.4190 * emsq2 + 156.5910 * eoc;
          g322 = -18.9068 + 109.7927 * em - 214.6334 * emsq2 + 146.5816 * eoc;
          g410 = -41.122 + 242.6940 * em - 471.0940 * emsq2 + 313.9530 * eoc;
          g422 = -146.407 + 841.8800 * em - 1629.014 * emsq2 + 1083.4350 * eoc;
          g520 = -532.114 + 3017.977 * em - 5740.032 * emsq2 + 3708.2760 * eoc;
        } else {
          g211 = -72.099 + 331.819 * em - 508.738 * emsq2 + 266.724 * eoc;
          g310 = -346.844 + 1582.851 * em - 2415.925 * emsq2 + 1246.113 * eoc;
          g322 = -342.585 + 1554.908 * em - 2366.899 * emsq2 + 1215.972 * eoc;
          g410 = -1052.797 + 4758.686 * em - 7193.992 * emsq2 + 3651.957 * eoc;
          g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq2 + 12422.520 * eoc;
          if (em > 0.715) g520 = -5149.66 + 29936.92 * em - 54087.36 * emsq2 + 31324.56 * eoc;
          else g520 = 1464.74 - 4664.75 * em + 3763.64 * emsq2;
        }
        if (em < 0.7) {
          g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq2 + 5542.21 * eoc;
          g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq2 + 5337.524 * eoc;
          g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq2 + 5341.4 * eoc;
        } else {
          g533 = -37995.780 + 161616.52 * em - 229838.20 * emsq2 + 109377.94 * eoc;
          g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq2 + 146349.42 * eoc;
          g532 = -40023.880 + 170470.89 * em - 242699.48 * emsq2 + 115605.82 * eoc;
        }

        const double sini2 = sinim * sinim;
        const double f220 = 0.75 * (1.0 + 2.0 * cosim + cosisq);
        const double f221 = 1.5 * sini2;
        const double f321 = 1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq);
        const double f322 = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq);
        const double f441 = 35.0 * sini2 * f220;
        const double f442 = 39.3750 * sini2 * sini2;
        const double f522 = 9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq) + 0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
        const double f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim + 10.0 * cosisq) + 6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
        const double f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim + cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
        const double f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim + cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq));
        const double xno2 = nm * nm;
        const double ainv2 = aonv * aonv;
        double temp1 = 3.0 * xno2 * ainv2;
        double temp = temp1 * root22;
        s.d2201 = temp * f220 * g201;
        s.d2211 = temp * f221 * g211;
        temp1 = temp1 * aonv;
        temp = temp1 * root32;
        s.d3210 = temp * f321 * g310;
        s.d3222 = temp * f322 * g322;
        temp1 = temp1 * aonv;
        temp = 2.0 * temp1 * root44;
        s.d4410 = temp * f441 * g410;
        s.d4422 = temp * f442 * g422;
        temp1 = temp1 * aonv;
        temp = temp1 * root52;
        s.d5220 = temp * f522 * g520;
        s.d5232 = temp * f523 * g532;
        temp = 2.0 * temp1 * root54;
        s.d5421 = temp * f542 * g521;
        s.d5433 = temp * f543 * g533;
        s.xlamo = std::fmod(mo + nodeo + nodeo - theta - theta, twopi);
        s.xfact = mdot + s.dmdt + 2.0 * (nodedot + s.dnodt - rptim) - no;
        em = emo;
        (void)emsqo;  // emsq is restored, as Vallado's, though it is not used after
      }

      //  ---------------- synchronous resonance terms --------------
      if (s.irez == 1) {
        const double g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
        const double g310 = 1.0 + 2.0 * emsq;
        const double g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
        const double f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
        const double f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
        double f330 = 1.0 + cosim;
        f330 = 1.875 * f330 * f330 * f330;
        s.del1 = 3.0 * nm * nm * aonv * aonv;
        s.del2 = 2.0 * s.del1 * f220 * g200 * q22;
        s.del3 = 3.0 * s.del1 * f330 * g300 * q33 * aonv;
        s.del1 = s.del1 * f311 * g310 * q31 * aonv;
        s.xlamo = std::fmod(mo + nodeo + argpo - theta, twopi);
        s.xfact = mdot + xpidot - rptim + s.dmdt + s.domdt + s.dnodt - no;
      }

      //  ------------ for sgp4, initialize the integrator ----------
      s.xli = s.xlamo;
      s.xni = no;
      s.atime = 0.0;
      nm = no + dndt;
    }
  }

  // ---- procedure dspace: deep space contributions to the mean elements --------------------------------
  // The resonance is integrated from the epoch at each call (atime, xli and xni are this call's).
  static void dspace(const Record& s, double t, double tc, double& em, double& argpm, double& inclm, double& mm, double& nodem,
                     double& dndt, double& nm) {
    const double fasx2 = 0.13130908, fasx4 = 2.8843198, fasx6 = 0.37448087, g22 = 5.7686396, g32 = 0.95240898,
                 g44 = 1.8014998, g52 = 1.0508330, g54 = 4.4108898, rptim = 4.37526908801129966e-3, stepp = 720.0,
                 stepn = -720.0, step2 = 259200.0;

    //  ----------- calculate deep space resonance effects -----------
    dndt = 0.0;
    const double theta = std::fmod(s.gsto + tc * rptim, twopi);
    em = em + s.dedt * t;
    inclm = inclm + s.didt * t;
    argpm = argpm + s.domdt * t;
    nodem = nodem + s.dnodt * t;
    mm = mm + s.dmdt * t;

    //  - update resonances : numerical (euler-maclaurin) integration -
    double ft = 0.0;
    if (s.irez != 0) {
      //  sgp4fix streamline check: the integration starts at the epoch
      double atime = 0.0, xni = s.no_unkozai, xli = s.xlamo;
      const double delt = t > 0.0 ? stepp : stepn;
      int iretn = 381;
      double xndt = 0, xldot = 0, xnddt = 0;
      while (iretn == 381) {
        //  ------------------- dot terms calculated -------------
        //  ----------- near - synchronous resonance terms -------
        if (s.irez != 2) {
          xndt = s.del1 * std::sin(xli - fasx2) + s.del2 * std::sin(2.0 * (xli - fasx4)) + s.del3 * std::sin(3.0 * (xli - fasx6));
          xldot = xni + s.xfact;
          xnddt = s.del1 * std::cos(xli - fasx2) + 2.0 * s.del2 * std::cos(2.0 * (xli - fasx4)) +
                  3.0 * s.del3 * std::cos(3.0 * (xli - fasx6));
          xnddt = xnddt * xldot;
        } else {
          // --------- near - half-day resonance terms --------
          const double xomi = s.argpo + s.argpdot * atime;
          const double x2omi = xomi + xomi;
          const double x2li = xli + xli;
          xndt = s.d2201 * std::sin(x2omi + xli - g22) + s.d2211 * std::sin(xli - g22) + s.d3210 * std::sin(xomi + xli - g32) +
                 s.d3222 * std::sin(-xomi + xli - g32) + s.d4410 * std::sin(x2omi + x2li - g44) + s.d4422 * std::sin(x2li - g44) +
                 s.d5220 * std::sin(xomi + xli - g52) + s.d5232 * std::sin(-xomi + xli - g52) +
                 s.d5421 * std::sin(xomi + x2li - g54) + s.d5433 * std::sin(-xomi + x2li - g54);
          xldot = xni + s.xfact;
          xnddt = s.d2201 * std::cos(x2omi + xli - g22) + s.d2211 * std::cos(xli - g22) + s.d3210 * std::cos(xomi + xli - g32) +
                  s.d3222 * std::cos(-xomi + xli - g32) + s.d5220 * std::cos(xomi + xli - g52) +
                  s.d5232 * std::cos(-xomi + xli - g52) +
                  2.0 * (s.d4410 * std::cos(x2omi + x2li - g44) + s.d4422 * std::cos(x2li - g44) +
                         s.d5421 * std::cos(xomi + x2li - g54) + s.d5433 * std::cos(-xomi + x2li - g54));
          xnddt = xnddt * xldot;
        }

        //  ----------------------- integrator -------------------
        //  sgp4fix move end checks to end of routine
        if (std::fabs(t - atime) >= stepp) {
          iretn = 381;
        } else {
          ft = t - atime;
          iretn = 0;
        }
        if (iretn == 381) {
          xli = xli + xldot * delt + xndt * step2;
          xni = xni + xndt * delt + xnddt * step2;
          atime = atime + delt;
        }
      }

      nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
      const double xl = xli + xldot * ft + xndt * ft * ft * 0.5;
      if (s.irez != 1) {
        mm = xl - 2.0 * nodem + 2.0 * theta;
        dndt = nm - s.no_unkozai;
      } else {
        mm = xl - nodem - argpm + theta;
        dndt = nm - s.no_unkozai;
      }
      nm = s.no_unkozai + dndt;
    }
  }

  // ---- procedure initl and sgp4init -------------------------------------------------------------------
  void init(const GravityConstants& g, char opsmode, const Sgp4Elements& e) {
    Record& s = r_;
    gravity_ = g;
    const double temp4 = 1.5e-12;
    //  ------------------------ earth constants -----------------------
    s.tumin = g.tumin;
    s.mu = g.mu;
    s.radiusearthkm = g.radiusearthkm;
    s.xke = g.xke;
    s.j2 = g.j2;
    s.j3 = g.j3;
    s.j4 = g.j4;
    s.j3oj2 = g.j3oj2;
    s.operationmode = opsmode;
    s.bstar = e.bstar;
    s.ndot = e.ndot;
    s.nddot = e.nddot;
    s.ecco = e.ecco;
    s.argpo = e.argpo;
    s.inclo = e.inclo;
    s.mo = e.mo;
    s.no_kozai = e.no_kozai;
    s.nodeo = e.nodeo;
    const double epoch = e.epoch;

    const double ss = 78.0 / s.radiusearthkm + 1.0;
    //  sgp4fix use multiply for speed instead of pow
    const double qzms2ttemp = (120.0 - 78.0) / s.radiusearthkm;
    const double qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
    const double x2o3 = 2.0 / 3.0;
    s.init = 'y';
    const double t = 0.0;

    //  ------------- initl: calculate auxillary epoch quantities ----------
    const double eccsq = s.ecco * s.ecco;
    const double omeosq = 1.0 - eccsq;
    const double rteosq = std::sqrt(omeosq);
    const double cosio = std::cos(s.inclo);
    const double cosio2 = cosio * cosio;
    //  ------------------ un-kozai the mean motion -----------------
    const double ak = std::pow(s.xke / s.no_kozai, x2o3);
    const double d1 = 0.75 * s.j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
    double del = d1 / (ak * ak);
    const double adel = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
    del = d1 / (adel * adel);
    s.no_unkozai = s.no_kozai / (1.0 + del);
    const double ao = std::pow(s.xke / s.no_unkozai, x2o3);
    const double sinio = std::sin(s.inclo);
    const double po = ao * omeosq;
    const double con42 = 1.0 - 5.0 * cosio2;
    s.con41 = -con42 - cosio2 - cosio2;
    const double posq = po * po;
    const double rp = ao * (1.0 - s.ecco);
    s.method = 'n';
    //  sgp4fix modern approach to finding sidereal time
    if (opsmode == 'a') {
      //  sgp4fix use old way of finding gst: count integer number of days from 0 jan 1970
      const double ts70 = epoch - 7305.0;
      const double ds70 = std::floor(ts70 + 1.0e-8);
      const double tfrac = ts70 - ds70;
      //  find greenwich location at epoch
      const double c1 = 1.72027916940703639e-2, thgr70 = 1.7321343856509374, fk5r = 5.07551419432269442e-15;
      const double c1p2p = c1 + twopi;
      s.gsto = std::fmod(thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r, twopi);
      if (s.gsto < 0.0) s.gsto = s.gsto + twopi;
    } else {
      s.gsto = gmst1982(epoch + 2433281.5);
    }

    s.a = std::pow(s.no_unkozai * s.tumin, -2.0 / 3.0);
    s.alta = s.a * (1.0 + s.ecco) - 1.0;
    s.altp = s.a * (1.0 - s.ecco) - 1.0;

    if (omeosq >= 0.0 || s.no_unkozai >= 0.0) {
      s.isimp = 0;
      if (rp < 220.0 / s.radiusearthkm + 1.0) s.isimp = 1;
      double sfour = ss, qzms24 = qzms2t;
      const double perige = (rp - 1.0) * s.radiusearthkm;
      //  - for perigees below 156 km, s and qoms2t are altered -
      if (perige < 156.0) {
        sfour = perige - 78.0;
        if (perige < 98.0) sfour = 20.0;
        //  sgp4fix use multiply for speed instead of pow
        const double qzms24temp = (120.0 - sfour) / s.radiusearthkm;
        qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
        sfour = sfour / s.radiusearthkm + 1.0;
      }
      const double pinvsq = 1.0 / posq;
      const double tsi = 1.0 / (ao - sfour);
      s.eta = ao * s.ecco * tsi;
      const double etasq = s.eta * s.eta;
      const double eeta = s.ecco * s.eta;
      const double psisq = std::fabs(1.0 - etasq);
      const double coef = qzms24 * std::pow(tsi, 4.0);
      const double coef1 = coef / std::pow(psisq, 3.5);
      const double cc2 = coef1 * s.no_unkozai *
                         (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
                          0.375 * s.j2 * tsi / psisq * s.con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
      s.cc1 = s.bstar * cc2;
      double cc3 = 0.0;
      if (s.ecco > 1.0e-4) cc3 = -2.0 * coef * tsi * s.j3oj2 * s.no_unkozai * sinio / s.ecco;
      s.x1mth2 = 1.0 - cosio2;
      s.cc4 = 2.0 * s.no_unkozai * coef1 * ao * omeosq *
              (s.eta * (2.0 + 0.5 * etasq) + s.ecco * (0.5 + 2.0 * etasq) -
               s.j2 * tsi / (ao * psisq) *
                   (-3.0 * s.con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
                    0.75 * s.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * std::cos(2.0 * s.argpo)));
      s.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
      const double cosio4 = cosio2 * cosio2;
      const double temp1 = 1.5 * s.j2 * pinvsq * s.no_unkozai;
      const double temp2 = 0.5 * temp1 * s.j2 * pinvsq;
      const double temp3 = -0.46875 * s.j4 * pinvsq * pinvsq * s.no_unkozai;
      s.mdot = s.no_unkozai + 0.5 * temp1 * rteosq * s.con41 + 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
      s.argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                  temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
      const double xhdot1 = -temp1 * cosio;
      s.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
      const double xpidot = s.argpdot + s.nodedot;
      s.omgcof = s.bstar * cc3 * std::cos(s.argpo);
      s.xmcof = 0.0;
      if (s.ecco > 1.0e-4) s.xmcof = -x2o3 * coef * s.bstar / eeta;
      s.nodecf = 3.5 * omeosq * xhdot1 * s.cc1;
      s.t2cof = 1.5 * s.cc1;
      //  sgp4fix for divide by zero with xinco = 180 deg
      if (std::fabs(cosio + 1.0) > 1.5e-12) s.xlcof = -0.25 * s.j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
      else s.xlcof = -0.25 * s.j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
      s.aycof = -0.5 * s.j3oj2 * sinio;
      //  sgp4fix use multiply for speed instead of pow
      const double delmotemp = 1.0 + s.eta * std::cos(s.mo);
      s.delmo = delmotemp * delmotemp * delmotemp;
      s.sinmao = std::sin(s.mo);
      s.x7thm1 = 7.0 * cosio2 - 1.0;

      //  --------------- deep space initialization -------------
      if (2.0 * pi / s.no_unkozai >= 225.0) {
        s.method = 'd';
        s.isimp = 1;
        const double tc = 0.0;
        double inclm = s.inclo;
        const Dscom d = dscom(epoch, s.ecco, s.argpo, tc, s.inclo, s.nodeo, s.no_unkozai, s);
        double em = d.em, nm = d.nm;
        // dpper at the epoch (init 'y') leaves the elements as they are.
        dpper(s, t, inclm, s.init, s.ecco, s.inclo, s.nodeo, s.argpo, s.mo, s.operationmode);
        double argpm = 0.0, nodem = 0.0, mm = 0.0, dndt = 0.0;
        dsinit(s.xke, d.cosim, d.emsq, s.argpo, d.s1, d.s2, d.s3, d.s4, d.s5, d.sinim, d.ss1, d.ss2, d.ss3, d.ss4, d.ss5, d.sz1,
               d.sz3, d.sz11, d.sz13, d.sz21, d.sz23, d.sz31, d.sz33, t, tc, s.gsto, s.mo, s.mdot, s.no_unkozai, s.nodeo,
               s.nodedot, xpidot, d.z1, d.z3, d.z11, d.z13, d.z21, d.z23, d.z31, d.z33, s.ecco, eccsq, em, argpm, inclm, mm, nm,
               nodem, s, dndt);
      }

      //----------- set variables if not deep space -----------
      if (s.isimp != 1) {
        const double cc1sq = s.cc1 * s.cc1;
        s.d2 = 4.0 * ao * tsi * cc1sq;
        const double temp = s.d2 * tsi * s.cc1 / 3.0;
        s.d3 = (17.0 * ao + sfour) * temp;
        s.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * s.cc1;
        s.t3cof = s.d2 + 2.0 * cc1sq;
        s.t4cof = 0.25 * (3.0 * s.d3 + s.cc1 * (12.0 * s.d2 + 10.0 * cc1sq));
        s.t5cof = 0.2 * (3.0 * s.d4 + 12.0 * s.cc1 * s.d3 + 6.0 * s.d2 * s.d2 + 15.0 * cc1sq * (2.0 * s.d2 + cc1sq));
      }
    }
    //  finally propogate to zero epoch to initialize all others.
    std::array<double, 3> r{}, v{};
    MeanElements m;
    init_error_ = propagate(0.0, r, v, m);
    s.init = 'n';
  }

  // ---- procedure sgp4 -----------------------------------------------------------------------------------
  int propagate(double tsince, std::array<double, 3>& r, std::array<double, 3>& v, MeanElements& mean) const {
    const Record& s = r_;
    const double temp4 = 1.5e-12;
    const double x2o3 = 2.0 / 3.0;
    const double vkmpersec = s.radiusearthkm * s.xke / 60.0;
    const double t = tsince;
    // The coefficients the deep-space branch recomputes at every call (Vallado's writes them into
    // the record; here they are this call's).
    double aycof = s.aycof, xlcof = s.xlcof, con41 = s.con41, x1mth2 = s.x1mth2, x7thm1 = s.x7thm1;
    double mrt = 0.0;

    //  ------- update for secular gravity and atmospheric drag -----
    const double xmdf = s.mo + s.mdot * t;
    const double argpdf = s.argpo + s.argpdot * t;
    const double nodedf = s.nodeo + s.nodedot * t;
    double argpm = argpdf;
    double mm = xmdf;
    const double t2 = t * t;
    double nodem = nodedf + s.nodecf * t2;
    double tempa = 1.0 - s.cc1 * t;
    double tempe = s.bstar * s.cc4 * t;
    double templ = s.t2cof * t2;

    if (s.isimp != 1) {
      const double delomg = s.omgcof * t;
      //  sgp4fix use mutliply for speed instead of pow
      const double delmtemp = 1.0 + s.eta * std::cos(xmdf);
      const double delm = s.xmcof * (delmtemp * delmtemp * delmtemp - s.delmo);
      const double temp = delomg + delm;
      mm = xmdf + temp;
      argpm = argpdf - temp;
      const double t3 = t2 * t;
      const double t4 = t3 * t;
      tempa = tempa - s.d2 * t2 - s.d3 * t3 - s.d4 * t4;
      tempe = tempe + s.bstar * s.cc5 * (std::sin(mm) - s.sinmao);
      templ = templ + s.t3cof * t3 + t4 * (s.t4cof + t * s.t5cof);
    }

    double nm = s.no_unkozai;
    double em = s.ecco;
    double inclm = s.inclo;
    if (s.method == 'd') {
      const double tc = t;
      double dndt = 0;
      dspace(s, t, tc, em, argpm, inclm, mm, nodem, dndt, nm);
    }

    if (nm <= 0.0) return 2;

    const double am = std::pow(s.xke / nm, x2o3) * tempa * tempa;
    nm = s.xke / std::pow(am, 1.5);
    em = em - tempe;

    //  fix tolerance for error recognition; sgp4fix am is fixed from the previous nm check
    if (em >= 1.0 || em < -0.001) return 1;
    //  sgp4fix fix tolerance to avoid a divide by zero
    if (em < 1.0e-6) em = 1.0e-6;
    mm = mm + s.no_unkozai * templ;
    double xlm = mm + argpm + nodem;

    nodem = std::fmod(nodem, twopi);
    argpm = std::fmod(argpm, twopi);
    xlm = std::fmod(xlm, twopi);
    mm = std::fmod(xlm - argpm - nodem, twopi);

    // sgp4fix recover singly averaged mean elements
    mean = {am, em, inclm, nodem, argpm, mm, nm};

    //  ----------------- compute extra mean quantities -------------
    const double sinim = std::sin(inclm), cosim = std::cos(inclm);

    //  -------------------- add lunar-solar periodics --------------
    double ep = em, xincp = inclm, argpp = argpm, nodep = nodem, mp = mm;
    double sinip = sinim, cosip = cosim;
    if (s.method == 'd') {
      dpper(s, t, s.inclo, 'n', ep, xincp, nodep, argpp, mp, s.operationmode);
      if (xincp < 0.0) {
        xincp = -xincp;
        nodep = nodep + pi;
        argpp = argpp - pi;
      }
      if (ep < 0.0 || ep > 1.0) return 3;
    }

    //  -------------------- long period periodics ------------------
    if (s.method == 'd') {
      sinip = std::sin(xincp);
      cosip = std::cos(xincp);
      aycof = -0.5 * s.j3oj2 * sinip;
      //  sgp4fix for divide by zero for xincp = 180 deg
      if (std::fabs(cosip + 1.0) > 1.5e-12) xlcof = -0.25 * s.j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
      else xlcof = -0.25 * s.j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
    }
    const double axnl = ep * std::cos(argpp);
    double temp = 1.0 / (am * (1.0 - ep * ep));
    const double aynl = ep * std::sin(argpp) + temp * aycof;
    const double xl = mp + argpp + nodep + temp * xlcof * axnl;

    //  --------------------- solve kepler's equation ---------------
    const double u = std::fmod(xl - nodep, twopi);
    double eo1 = u, tem5 = 9999.9, sineo1 = 0, coseo1 = 0;
    int ktr = 1;
    //    sgp4fix for kepler iteration: the following iteration needs better limits on corrections
    while (std::fabs(tem5) >= 1.0e-12 && ktr <= 10) {
      sineo1 = std::sin(eo1);
      coseo1 = std::cos(eo1);
      tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
      tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
      if (std::fabs(tem5) >= 0.95) tem5 = tem5 > 0.0 ? 0.95 : -0.95;
      eo1 = eo1 + tem5;
      ktr = ktr + 1;
    }

    //  ------------- short period preliminary quantities -----------
    const double ecose = axnl * coseo1 + aynl * sineo1;
    const double esine = axnl * sineo1 - aynl * coseo1;
    const double el2 = axnl * axnl + aynl * aynl;
    const double pl = am * (1.0 - el2);
    if (pl < 0.0) return 4;

    const double rl = am * (1.0 - ecose);
    const double rdotl = std::sqrt(am) * esine / rl;
    const double rvdotl = std::sqrt(pl) / rl;
    const double betal = std::sqrt(1.0 - el2);
    temp = esine / (1.0 + betal);
    const double sinu = am / rl * (sineo1 - aynl - axnl * temp);
    const double cosu = am / rl * (coseo1 - axnl + aynl * temp);
    double su = std::atan2(sinu, cosu);
    const double sin2u = (cosu + cosu) * sinu;
    const double cos2u = 1.0 - 2.0 * sinu * sinu;
    temp = 1.0 / pl;
    const double temp1 = 0.5 * s.j2 * temp;
    const double temp2 = temp1 * temp;

    //  -------------- update for short period periodics ------------
    if (s.method == 'd') {
      const double cosisq = cosip * cosip;
      con41 = 3.0 * cosisq - 1.0;
      x1mth2 = 1.0 - cosisq;
      x7thm1 = 7.0 * cosisq - 1.0;
    }
    mrt = rl * (1.0 - 1.5 * temp2 * betal * con41) + 0.5 * temp1 * x1mth2 * cos2u;
    su = su - 0.25 * temp2 * x7thm1 * sin2u;
    const double xnode = nodep + 1.5 * temp2 * cosip * sin2u;
    const double xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
    const double mvt = rdotl - nm * temp1 * x1mth2 * sin2u / s.xke;
    const double rvdot = rvdotl + nm * temp1 * (x1mth2 * cos2u + 1.5 * con41) / s.xke;

    //  --------------------- orientation vectors -------------------
    const double sinsu = std::sin(su), cossu = std::cos(su);
    const double snod = std::sin(xnode), cnod = std::cos(xnode);
    const double sini = std::sin(xinc), cosi = std::cos(xinc);
    const double xmx = -snod * cosi;
    const double xmy = cnod * cosi;
    const double ux = xmx * sinsu + cnod * cossu;
    const double uy = xmy * sinsu + snod * cossu;
    const double uz = sini * sinsu;
    const double vx = xmx * cossu - cnod * sinsu;
    const double vy = xmy * cossu - snod * sinsu;
    const double vz = sini * cossu;

    //  --------- position and velocity (in km and km/sec) ----------
    const double mr = mrt * s.radiusearthkm;
    r = {mr * ux, mr * uy, mr * uz};
    v = {(mvt * ux + rvdot * vx) * vkmpersec, (mvt * uy + rvdot * vy) * vkmpersec, (mvt * uz + rvdot * vz) * vkmpersec};

    //  sgp4fix for decaying satellites
    if (mrt < 1.0) return 6;
    return 0;
  }
};

}  // namespace navigation::orbits
