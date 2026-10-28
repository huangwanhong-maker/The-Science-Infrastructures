// Time scales and the conversions among them: a groupoid whose arrows are affine maps, some of
// them read from data.
//
// The scales (IAU Resolutions 1991 A4, 2000 B1.3–B1.9, 2006 B3; IERS Conventions 2010, ch. 10;
// Explanatory Supplement to the Astronomical Almanac, 3rd ed., 2013, ch. 3):
//   TAI  International Atomic Time, the realised atomic scale.
//   TT   Terrestrial Time, TAI + 32.184 s exactly: the time of geocentric ephemerides.
//   TDB  Barycentric Dynamical Time, TT plus a periodic term of 1.7 ms amplitude (Fairhead and
//        Bretagnon's series here, at the geocentre): the time of the JPL ephemerides.
//   TCG, TCB  Geocentric and Barycentric Coordinate Time, the coordinate times of the IAU's
//        relativistic frames: TT and TDB rescaled by L_G = 6.969290134 × 10⁻¹⁰ and
//        L_B = 1.550519768 × 10⁻⁸, with TDB₀ = −6.55 × 10⁻⁵ s, all at 1977 January 1, 0h TAI.
//   UTC  Coordinated Universal Time: TAI less the leap seconds (leap_seconds.hpp).
//   UT1  The Earth's rotation as a time: UTC + (UT1 − UTC), measured (earth_orientation).
//   GPS  GPS time, TAI − 19 s exactly.
//
// Each conversion goes through TAI, by SOFA's algorithms: UTC days that contain a leap second
// are 86 401 s long, and a UTC Julian date's fraction of such a day is scaled accordingly
// (iauUtctai, iauTaiutc). UT1 needs the IERS's measurements; where they do not reach, a policy
// says whether to refuse or to approximate (UT1 ≈ UTC after 1972, within 0.9 s; TT − ΔT by
// Espenak and Meeus's polynomials before), and every instant says which was done.

#pragma once

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include <metrology/earth_orientation/eop.hpp>
#include <metrology/time/julian.hpp>
#include <metrology/time/leap_seconds.hpp>
#include <metrology/time/series_tdb.hpp>

namespace metrology::time {

enum class Scale { TAI, TT, TDB, TCG, TCB, UTC, UT1, GPS };

inline std::string_view name(Scale s) {
  switch (s) {
    case Scale::TAI: return "TAI";
    case Scale::TT: return "TT";
    case Scale::TDB: return "TDB";
    case Scale::TCG: return "TCG";
    case Scale::TCB: return "TCB";
    case Scale::UTC: return "UTC";
    case Scale::UT1: return "UT1";
    default: return "GPS";
  }
}
inline Scale parse_scale(std::string_view s) {
  for (Scale x : {Scale::TAI, Scale::TT, Scale::TDB, Scale::TCG, Scale::TCB, Scale::UTC, Scale::UT1, Scale::GPS})
    if (name(x) == s) return x;
  if (s == "Z") return Scale::UTC;
  throw std::invalid_argument("metrology::time: unknown time scale " + std::string(s));
}

inline constexpr double TT_MINUS_TAI = 32.184;        // s, exact
inline constexpr double TAI_MINUS_GPS = 19.0;         // s, exact
inline constexpr double L_G = 6.969290134e-10;        // IAU 2000 B1.9
inline constexpr double L_B = 1.550519768e-8;         // IAU 2006 B3
inline constexpr double TDB0 = -6.55e-5;              // s, IAU 2006 B3
inline constexpr double MJD_1977 = 43144.0;           // 1977 January 1, 0h TAI

// TDB − TT at the geocentre, in seconds: Fairhead and Bretagnon's series in Julian millennia,
// with the adjustments to JPL's planetary masses; SOFA's iauDtdb with u = v = 0, whose
// topocentric terms (≤ 2 µs) then vanish. Accurate to a few nanoseconds within a few
// millennia of J2000.
inline double tdb_minus_tt(const JulianDate& t) {
  const double T = t.since_j2000() / JULIAN_MILLENNIUM;
  const auto& f = series::fairhead;
  auto sum = [&](int hi, int lo) {
    double w = 0;
    for (int j = hi; j >= lo; j--) w += f[j].amplitude * std::sin(f[j].frequency * T + f[j].phase);
    return w;
  };
  const double w0 = sum(473, 0), w1 = sum(678, 474), w2 = sum(763, 679), w3 = sum(783, 764), w4 = sum(786, 784);
  const double wf = T * (T * (T * (T * w4 + w3) + w2) + w1) + w0;
  const double wj = 0.00065e-6 * std::sin(6069.776754 * T + 4.021194) + 0.00033e-6 * std::sin(213.299095 * T + 5.543132) -
                    0.00196e-6 * std::sin(6208.294251 * T + 5.696701) - 0.00173e-6 * std::sin(74.781599 * T + 2.435900) +
                    0.03638e-6 * T * T;
  return wf + wj;
}

// ΔT = TT − UT1 in seconds from 1800 to 1972 (Espenak and Meeus 2006, NASA/TP-2006-214141,
// "Polynomial expressions for Delta T"; the metrology domain's JavaScript has the same).
inline double delta_t_espenak_meeus(double year) {
  if (year < 1800 || year >= 1973) throw std::domain_error("metrology::time: the ΔT polynomials here cover 1800 to 1972");
  if (year < 1860) {
    const double t = year - 1800;
    return 13.72 - 0.332447 * t + 0.0068612 * t * t + 0.0041116 * t * t * t - 0.00037436 * std::pow(t, 4) +
           0.0000121272 * std::pow(t, 5) - 0.0000001699 * std::pow(t, 6) + 0.000000000875 * std::pow(t, 7);
  }
  if (year < 1900) {
    const double t = year - 1860;
    return 7.62 + 0.5737 * t - 0.251754 * t * t + 0.01680668 * t * t * t - 0.0004473624 * std::pow(t, 4) + std::pow(t, 5) / 233174;
  }
  if (year < 1920) {
    const double t = year - 1900;
    return -2.79 + 1.494119 * t - 0.0598939 * t * t + 0.0061966 * t * t * t - 0.000197 * std::pow(t, 4);
  }
  if (year < 1941) {
    const double t = year - 1920;
    return 21.20 + 0.84493 * t - 0.076100 * t * t + 0.0020936 * t * t * t;
  }
  if (year < 1961) {
    const double t = year - 1950;
    return 29.07 + 0.407 * t - t * t / 233 + t * t * t / 2547;
  }
  const double t = year - 1975;
  return 45.45 + 1.067 * t - t * t / 260 - t * t * t / 718;
}

// How UT1 was obtained for an instant.
enum class Ut1Source {
  measured,     // the IERS's Bulletin A, a measurement
  predicted,    // the IERS's Bulletin A, a prediction
  approximate,  // beyond the IERS's series after 1972: UT1 taken as UTC, within 0.9 s
  delta_t,      // before the series: TT − ΔT by Espenak and Meeus, within a second or so after 1800
};
inline std::string_view name(Ut1Source s) {
  switch (s) {
    case Ut1Source::measured: return "IERS Bulletin A, measured";
    case Ut1Source::predicted: return "IERS Bulletin A, predicted";
    case Ut1Source::approximate: return "approximate: UT1 = UTC, within 0.9 s";
    default: return "approximate: TT - Delta T (Espenak and Meeus 2006)";
  }
}

// What to do for UT1 beyond the Earth orientation series.
enum class Ut1Policy { refuse, approximate };

// An instant in every scale, with the Earth orientation that holds at it.
struct Instant {
  JulianDate tai, tt, tdb, ut1, utc;
  bool has_utc = false;  // false before 1960
  earth_orientation::Values eop;  // interpolated; polar motion and pole offsets zero when approximated
  Ut1Source ut1_source = Ut1Source::measured;

  // Seconds past J2000 TDB: the time argument of the JPL ephemerides.
  double ephemeris_seconds() const { return tdb.since_j2000() * DAY; }
};

class TimeScales {
 public:
  explicit TimeScales(LeapSeconds leap = LeapSeconds::built_in(), const earth_orientation::Series* eop = nullptr,
                      Ut1Policy policy = Ut1Policy::approximate)
      : leap_(std::move(leap)), eop_(eop), policy_(policy) {}

  const LeapSeconds& leap_seconds() const { return leap_; }
  const earth_orientation::Series* earth_orientation() const { return eop_; }

  // ---- UTC and TAI (SOFA's iauUtctai and iauTaiutc) ------------------------------------------
  JulianDate utc_to_tai(const JulianDate& utc) const {
    const bool big1 = std::fabs(utc.d1) >= std::fabs(utc.d2);
    const double u1 = big1 ? utc.d1 : utc.d2, u2 = big1 ? utc.d2 : utc.d1;
    const CalendarDate today = jd_to_calendar({u1, u2});
    const double dat0 = leap_.tai_minus_utc(today.year, today.month, today.day, 0.0);
    const double dat12 = leap_.tai_minus_utc(today.year, today.month, today.day, 0.5);
    const CalendarDate tomorrow = jd_to_calendar({u1 + 1.5, u2 - today.fraction});
    const double dat24 = leap_.tai_minus_utc(tomorrow.year, tomorrow.month, tomorrow.day, 0.0);
    const double dlod = 2 * (dat12 - dat0), dleap = dat24 - (dat0 + dlod);
    double fd = today.fraction * (DAY + dleap) / DAY;
    fd *= (DAY + dlod) / DAY;
    const JulianDate z = calendar_to_jd(today.year, today.month, today.day);
    double a2 = z.d1 - u1;
    a2 += z.d2;
    a2 += fd + dat0 / DAY;
    return big1 ? JulianDate{u1, a2} : JulianDate{a2, u1};
  }
  JulianDate tai_to_utc(const JulianDate& tai) const {
    const bool big1 = std::fabs(tai.d1) >= std::fabs(tai.d2);
    const double a1 = big1 ? tai.d1 : tai.d2, a2 = big1 ? tai.d2 : tai.d1;
    double u1 = a1, u2 = a2;
    for (int i = 0; i < 3; i++) {
      const JulianDate g = utc_to_tai({u1, u2});
      u2 += a1 - g.d1;
      u2 += a2 - g.d2;
    }
    return big1 ? JulianDate{u1, u2} : JulianDate{u2, u1};
  }
  // The length of the UTC day containing a UTC date, in SI seconds (86 401 on a leap second day).
  double utc_day_length(const JulianDate& utc) const {
    const CalendarDate today = jd_to_calendar(utc);
    const double dat0 = leap_.tai_minus_utc(today.year, today.month, today.day, 0.0);
    const double dat12 = leap_.tai_minus_utc(today.year, today.month, today.day, 0.5);
    const CalendarDate tomorrow = jd_to_calendar({utc.d1 + 1.5, utc.d2 - today.fraction});
    const double dat24 = leap_.tai_minus_utc(tomorrow.year, tomorrow.month, tomorrow.day, 0.0);
    return DAY + (dat24 - (2 * dat12 - dat0));
  }

  // ---- TT, TDB, TCG, TCB ------------------------------------------------------------------------
  static JulianDate tt_to_tdb(const JulianDate& tt) { return tt.plus_seconds(tdb_minus_tt(tt)); }
  static JulianDate tdb_to_tt(const JulianDate& tdb) {
    JulianDate tt = tdb.plus_seconds(-tdb_minus_tt(tdb));
    return tdb.plus_seconds(-tdb_minus_tt(tt));
  }
  static JulianDate tt_to_tcg(const JulianDate& tt) {
    const double t77t = MJD_1977 + TT_MINUS_TAI / DAY, elgg = L_G / (1 - L_G);
    if (std::fabs(tt.d1) > std::fabs(tt.d2)) return {tt.d1, tt.d2 + ((tt.d1 - MJD_ZERO) + (tt.d2 - t77t)) * elgg};
    return {tt.d1 + ((tt.d2 - MJD_ZERO) + (tt.d1 - t77t)) * elgg, tt.d2};
  }
  static JulianDate tcg_to_tt(const JulianDate& tcg) {
    const double t77t = MJD_1977 + TT_MINUS_TAI / DAY;
    if (std::fabs(tcg.d1) > std::fabs(tcg.d2)) return {tcg.d1, tcg.d2 - ((tcg.d1 - MJD_ZERO) + (tcg.d2 - t77t)) * L_G};
    return {tcg.d1 - ((tcg.d2 - MJD_ZERO) + (tcg.d1 - t77t)) * L_G, tcg.d2};
  }
  static JulianDate tdb_to_tcb(const JulianDate& tdb) {
    const double t77td = MJD_ZERO + MJD_1977, t77tf = TT_MINUS_TAI / DAY, tdb0 = TDB0 / DAY, elbb = L_B / (1 - L_B);
    if (std::fabs(tdb.d1) > std::fabs(tdb.d2)) {
      const double d = t77td - tdb.d1, f = tdb.d2 - tdb0;
      return {tdb.d1, f - (d - (f - t77tf)) * elbb};
    }
    const double d = t77td - tdb.d2, f = tdb.d1 - tdb0;
    return {f - (d - (f - t77tf)) * elbb, tdb.d2};
  }
  static JulianDate tcb_to_tdb(const JulianDate& tcb) {
    const double t77td = MJD_ZERO + MJD_1977, t77tf = TT_MINUS_TAI / DAY, tdb0 = TDB0 / DAY;
    if (std::fabs(tcb.d1) > std::fabs(tcb.d2)) {
      const double d = tcb.d1 - t77td;
      return {tcb.d1, tcb.d2 + tdb0 - (d + (tcb.d2 - t77tf)) * L_B};
    }
    const double d = tcb.d2 - t77td;
    return {tcb.d1 + tdb0 - (d + (tcb.d1 - t77tf)) * L_B, tcb.d2};
  }

  // ---- UT1 --------------------------------------------------------------------------------------
  struct Dut1 {
    double seconds = 0;
    Ut1Source source = Ut1Source::measured;
    earth_orientation::Values eop;
  };
  // UT1 − UTC at a UTC MJD, with where it came from.
  Dut1 dut1(double mjd_utc) const {
    if (eop_ && eop_->covers(mjd_utc)) {
      const auto v = eop_->at(mjd_utc);
      return {v.dut1, v.predicted ? Ut1Source::predicted : Ut1Source::measured, v};
    }
    if (policy_ == Ut1Policy::refuse)
      throw std::domain_error("metrology::time: no Earth orientation for UT1 at MJD " + std::to_string(mjd_utc) +
                              (eop_ ? " (outside the IERS series)" : " (no series given)"));
    earth_orientation::Values v;
    v.mjd = mjd_utc;
    if (mjd_utc >= 41317.0) return {0.0, Ut1Source::approximate, v};  // from 1972 January 1
    // 1960–1971: UT1 − UTC = (TT − UTC) − ΔT.
    const CalendarDate c = jd_to_calendar(JulianDate::from_mjd(mjd_utc));
    const double year = c.year + (c.month - 0.5) / 12;
    const double tt_minus_utc = leap_.tai_minus_utc(c.year, c.month, c.day, c.fraction) + TT_MINUS_TAI;
    v.dut1 = tt_minus_utc - delta_t_espenak_meeus(year);
    return {v.dut1, Ut1Source::delta_t, v};
  }

  // ---- any scale to any other ---------------------------------------------------------------------
  JulianDate to_tai(const JulianDate& t, Scale s) const {
    switch (s) {
      case Scale::TAI: return t;
      case Scale::TT: return t.plus_seconds(-TT_MINUS_TAI);
      case Scale::TDB: return tdb_to_tt(t).plus_seconds(-TT_MINUS_TAI);
      case Scale::TCG: return tcg_to_tt(t).plus_seconds(-TT_MINUS_TAI);
      case Scale::TCB: return tdb_to_tt(tcb_to_tdb(t)).plus_seconds(-TT_MINUS_TAI);
      case Scale::UTC: return utc_to_tai(t);
      case Scale::GPS: return t.plus_seconds(TAI_MINUS_GPS);
      case Scale::UT1: return ut1_to_tai(t);
    }
    throw std::logic_error("unreachable");
  }
  JulianDate from_tai(const JulianDate& tai, Scale s) const {
    switch (s) {
      case Scale::TAI: return tai;
      case Scale::TT: return tai.plus_seconds(TT_MINUS_TAI);
      case Scale::TDB: return tt_to_tdb(tai.plus_seconds(TT_MINUS_TAI));
      case Scale::TCG: return tt_to_tcg(tai.plus_seconds(TT_MINUS_TAI));
      case Scale::TCB: return tdb_to_tcb(tt_to_tdb(tai.plus_seconds(TT_MINUS_TAI)));
      case Scale::UTC: return tai_to_utc(tai);
      case Scale::GPS: return tai.plus_seconds(-TAI_MINUS_GPS);
      case Scale::UT1: return tai_to_ut1(tai);
    }
    throw std::logic_error("unreachable");
  }
  JulianDate convert(const JulianDate& t, Scale from, Scale to) const {
    if (from == to) return t;
    // The two dynamical scales are related directly, without UTC's range.
    if (from == Scale::TT && to == Scale::TDB) return tt_to_tdb(t);
    if (from == Scale::TDB && to == Scale::TT) return tdb_to_tt(t);
    return from_tai(to_tai(t, from), to);
  }

  // Every scale at an instant given in one, with the Earth orientation there.
  Instant instant(const JulianDate& t, Scale s) const {
    Instant i;
    i.tai = to_tai(t, s);
    i.tt = i.tai.plus_seconds(TT_MINUS_TAI);
    i.tdb = tt_to_tdb(i.tt);
    i.has_utc = utc_defined(i.tai);
    if (i.has_utc) {
      i.utc = s == Scale::UTC ? t : tai_to_utc(i.tai);
      const Dut1 d = dut1(i.utc.mjd());
      i.ut1 = i.utc.plus_seconds(d.seconds);
      i.eop = d.eop;
      i.ut1_source = d.source;
    } else {
      i.ut1 = s == Scale::UT1 ? t : tt_to_ut1_by_delta_t(i.tt);
      i.ut1_source = Ut1Source::delta_t;
      i.eop.mjd = i.tt.mjd();
    }
    return i;
  }

 private:
  static bool utc_defined(const JulianDate& tai) { return tai.mjd() >= 36934.0; }  // 1960 January 1
  static double decimal_year(const JulianDate& t) {
    const CalendarDate c = jd_to_calendar(t);
    return c.year + (c.month - 0.5) / 12;
  }
  JulianDate tt_to_ut1_by_delta_t(const JulianDate& tt) const {
    if (policy_ == Ut1Policy::refuse) throw std::domain_error("metrology::time: no UT1 before 1960 without a ΔT model");
    return tt.plus_seconds(-delta_t_espenak_meeus(decimal_year(tt)));
  }
  JulianDate ut1_to_tai(const JulianDate& ut1) const {
    // Before 1960: TT = UT1 + ΔT.
    if (ut1.mjd() < 36934.0) {
      if (policy_ == Ut1Policy::refuse) throw std::domain_error("metrology::time: no UT1 before 1960 without a ΔT model");
      return ut1.plus_seconds(delta_t_espenak_meeus(decimal_year(ut1)) - TT_MINUS_TAI);
    }
    // UTC = UT1 − (UT1 − UTC)(UTC), by iteration.
    JulianDate utc = ut1;
    for (int k = 0; k < 3; k++) utc = ut1.plus_seconds(-dut1(utc.mjd()).seconds);
    return utc_to_tai(utc);
  }
  JulianDate tai_to_ut1(const JulianDate& tai) const {
    if (!utc_defined(tai)) return tt_to_ut1_by_delta_t(tai.plus_seconds(TT_MINUS_TAI));
    const JulianDate utc = tai_to_utc(tai);
    return utc.plus_seconds(dut1(utc.mjd()).seconds);
  }

  LeapSeconds leap_;
  const earth_orientation::Series* eop_ = nullptr;
  Ut1Policy policy_ = Ut1Policy::approximate;
};

}  // namespace metrology::time
