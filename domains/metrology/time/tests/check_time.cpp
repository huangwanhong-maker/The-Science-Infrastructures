// Checks of the time scales, the calendar, leap seconds and the Earth orientation series
// (metrology/time, metrology/earth_orientation).
//
// References: ERFA's own test values (t_erfa_c.c, the SOFA test suite) written in below, and
// tests/fixtures/erfa-reference.txt, ERFA evaluated through pyerfa at dates from 1850 to 2200;
// the IERS's finals2000A for the rows read here; IANA's leap-seconds.list.

#include <cmath>
#include <string>

#include <calx/verification/check.hpp>
#include <metrology/earth_orientation/eop.hpp>
#include <metrology/time/datetime.hpp>
#include <metrology/time/julian.hpp>
#include <metrology/time/leap_seconds.hpp>
#include <metrology/time/scales.hpp>

#include "../../tests/reference.hpp"

using namespace metrology::time;
namespace eo = metrology::earth_orientation;
using calx::verification::Suite;

namespace {

// Rows of finals2000A around the leap second of 2017 January 1, as the IERS distributed them on
// 2026-10-28 (https://datacenter.iers.org/data/9/finals2000A.all), and its last two rows' forms.
const char* FINALS_2017 =
    "161230 57752.00 I  0.082883 0.000052  0.263539 0.000038  I-0.4069180 0.0000079  0.8286 0.0056  I     0.022    0.119    -0.155    0.024  0.082924  0.263518 -0.4069106    -0.023    -0.038\n"
    "161231 57753.00 I  0.081400 0.000052  0.263094 0.000039  I-0.4077601 0.0000078  0.8842 0.0055  I     0.025    0.119    -0.169    0.024  0.081318  0.262990 -0.4077600    -0.021    -0.048\n"
    "17 1 1 57754.00 I  0.080504 0.000028  0.263145 0.000028  I 0.5912821 0.0000077  1.0342 0.0050  I     0.012    0.119    -0.168    0.018  0.080450  0.263074  0.5912975    -0.019    -0.057\n"
    "17 1 2 57755.00 I  0.080285 0.000029  0.263605 0.000032  I 0.5901752 0.0000063  1.1730 0.0048  I    -0.001    0.040    -0.156    0.033  0.080275  0.263595  0.5902149    -0.027    -0.065\n"
    "2610 1 61314.00 I  0.174596 0.000026  0.325269 0.000026  I-0.0225232 0.0000134  0.6645 0.0074  I     0.365    0.050    -0.104    0.044\n"
    "271217 61756.00\n";

// Lines of IANA's leap-seconds.list (retrieved 2026-10-28), with its expiry.
const char* IANA_LIST =
    "#\tUpdated through IERS Bulletin C (https://hpiers.obspm.fr/iers/bul/bulc/bulletinc.dat)\n"
    "#@\t4023129600\n"
    "2272060800\t10\t# 1 Jan 1972\n"
    "2287785600\t11\t# 1 Jul 1972\n"
    "2303683200\t12\t# 1 Jan 1973\n"
    "3644697600\t36\t# 1 Jul 2015\n"
    "3692217600\t37\t# 1 Jan 2017\n";

double days_apart(const JulianDate& a, const JulianDate& b) { return std::fabs(a.minus(b)); }

}  // namespace

int main() {
  Suite s;
  auto ref = metrology::checks::read_reference(METROLOGY_FIXTURES "/erfa-reference.txt");
  const TimeScales ts;

  // ---- the calendar ----------------------------------------------------------------------------
  {
    const JulianDate j = calendar_to_jd(2003, 6, 1);
    s.require("cal2jd 2003-06-01 is MJD 52791 (SOFA test)", j.mjd() == 52791.0);
    const CalendarDate c = jd_to_calendar({2400000.5, 50123.9999});
    s.require("jd2cal MJD 50123.9999 is 1996-02-10 (SOFA test)", c.year == 1996 && c.month == 2 && c.day == 10);
    s.check("jd2cal fraction 0.9999 (SOFA test)", c.fraction - 0.9999, 1e-7);
    double worst = 0;
    bool dates = true;
    for (const auto& r : ref["cal2jd"]) {
      const JulianDate x = calendar_to_jd(r.integer(0), r.integer(1), r.integer(2));
      worst = std::max(worst, std::fabs(x.mjd() - r.number(4)));
    }
    for (const auto& r : ref["jd2cal"]) {
      const CalendarDate x = jd_to_calendar({r.number(0), r.number(1)});
      dates = dates && x.year == r.integer(2) && x.month == r.integer(3) && x.day == r.integer(4);
      worst = std::max(worst, std::fabs(x.fraction - r.number(5)));
    }
    s.check("the calendar against ERFA at 23 dates, 1850–2200", worst, 1e-15, "day");
    s.require("the calendar's dates agree with ERFA", dates && !ref["jd2cal"].empty());
    bool refused = false;
    try {
      calendar_to_jd(2023, 2, 29);
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("2023-02-29 is refused", refused);
  }

  // ---- leap seconds ----------------------------------------------------------------------------
  {
    const LeapSeconds leap = LeapSeconds::built_in();
    s.check("TAI − UTC on 2003-06-01 is 32 s (SOFA test)", leap.tai_minus_utc(2003, 6, 1) - 32, 0);
    s.check("TAI − UTC on 2008-01-17 is 33 s (SOFA test)", leap.tai_minus_utc(2008, 1, 17) - 33, 0);
    s.check("TAI − UTC on 2017-09-01 is 37 s (SOFA test)", leap.tai_minus_utc(2017, 9, 1) - 37, 0);
    double worst = 0;
    for (const auto& r : ref["dat"]) worst = std::max(worst, std::fabs(leap.tai_minus_utc(r.integer(0), r.integer(1), r.integer(2), r.number(3)) - r.number(4)));
    s.check("TAI − UTC against ERFA, with the drifting offsets of 1960–1971", worst, 1e-12, "s");
    const LeapSeconds iana = LeapSeconds::parse_iana(IANA_LIST);
    s.check("IANA's list read: TAI − UTC on 2016-12-31 is 36 s", iana.tai_minus_utc(2016, 12, 31) - 36, 0);
    s.check("IANA's list read: 37 s from 2017-01-01", iana.tai_minus_utc(2017, 1, 1) - 37, 0);
    s.check("IANA's list read: the drifting offsets kept", iana.tai_minus_utc(1965, 3, 1, 0.25) - leap.tai_minus_utc(1965, 3, 1, 0.25), 0);
    s.check("IANA's list valid until 2027-06-28 (MJD 61584)", iana.valid_until_mjd() - 61584, 1e-9, "day");
    bool refused = false;
    try {
      (void)leap.tai_minus_utc(1959, 12, 31);
    } catch (const std::domain_error&) {
      refused = true;
    }
    s.require("UTC before 1960 is refused", refused);
  }

  // ---- UTC and TAI -----------------------------------------------------------------------------
  {
    const JulianDate tai = ts.utc_to_tai({2453750.5, 0.892100694});
    s.check("UTC → TAI (SOFA test)", tai.d2 - 0.8924826384444444444, 1e-12, "day");
    const JulianDate utc = ts.tai_to_utc({2453750.5, 0.892482639});
    s.check("TAI → UTC (SOFA test)", utc.d2 - 0.8921006945555555556, 1e-12, "day");
    double worst = 0;
    for (const auto& r : ref["utctai"]) worst = std::max(worst, days_apart(ts.utc_to_tai({r.number(0), r.number(1)}), {r.number(2), r.number(3)}));
    for (const auto& r : ref["taiutc"]) worst = std::max(worst, days_apart(ts.tai_to_utc({r.number(0), r.number(1)}), {r.number(2), r.number(3)}));
    s.check("UTC ↔ TAI against ERFA across leap seconds and the 1960s", worst * DAY, 1e-9, "s");
    s.check("the UTC day of 2016-12-31 lasts 86 401 s", ts.utc_day_length(calendar_to_jd(2016, 12, 31)) - 86401, 1e-9, "s");
    s.check("the UTC day of 2016-12-30 lasts 86 400 s", ts.utc_day_length(calendar_to_jd(2016, 12, 30)) - 86400, 1e-9, "s");
  }

  // ---- the dynamical and coordinate scales -----------------------------------------------------
  {
    const JulianDate g = TimeScales::tt_to_tcg({2453750.5, 0.892482639});
    s.check("TT → TCG (SOFA test)", g.d2 - 0.8924900312508587113, 1e-12, "day");
    const JulianDate b = TimeScales::tdb_to_tcb({2453750.5, 0.892855137});
    s.check("TDB → TCB (SOFA test)", b.d2 - 0.8930195997253656716, 1e-12, "day");
    double tdb = 0, tcg = 0, tcb = 0, back = 0;
    for (const auto& r : ref["dtdb"]) tdb = std::max(tdb, std::fabs(tdb_minus_tt({r.number(0), r.number(1)}) - r.number(2)));
    for (const auto& r : ref["tttcg"]) {
      const JulianDate x = TimeScales::tt_to_tcg({r.number(0), r.number(1)});
      tcg = std::max(tcg, days_apart(x, {r.number(2), r.number(3)}));
      back = std::max(back, days_apart(TimeScales::tcg_to_tt(x), {r.number(0), r.number(1)}));
    }
    for (const auto& r : ref["tdbtcb"]) {
      const JulianDate x = TimeScales::tdb_to_tcb({r.number(0), r.number(1)});
      tcb = std::max(tcb, days_apart(x, {r.number(2), r.number(3)}));
      back = std::max(back, days_apart(TimeScales::tcb_to_tdb(x), {r.number(0), r.number(1)}));
    }
    s.check("TDB − TT (Fairhead–Bretagnon, 787 terms) against ERFA, 1850–2200", tdb, 1e-15, "s");
    s.check("TT → TCG against ERFA", tcg * DAY, 1e-9, "s");
    s.check("TDB → TCB against ERFA", tcb * DAY, 1e-9, "s");
    s.check("TCG → TT and TCB → TDB invert them", back * DAY, 1e-9, "s");
    const JulianDate tt{2451545.0, 0.123};
    s.check("TT → TDB → TT", days_apart(TimeScales::tdb_to_tt(TimeScales::tt_to_tdb(tt)), tt) * DAY, 1e-12, "s");
  }

  // ---- conversions among every scale -----------------------------------------------------------
  {
    const Scale all[] = {Scale::TAI, Scale::TT, Scale::TDB, Scale::TCG, Scale::TCB, Scale::UTC, Scale::UT1, Scale::GPS};
    double worst = 0;
    for (const double mjd : {41317.25, 51544.5, 57753.999, 57754.0001, 61341.5}) {
      const JulianDate t = JulianDate::from_mjd(mjd);
      for (Scale a : all)
        for (Scale b : all) worst = std::max(worst, days_apart(ts.convert(ts.convert(t, a, b), b, a), t));
    }
    s.check("every scale to every other and back, at five dates (approximate UT1)", worst * DAY, 2e-9, "s");
    const JulianDate tai = calendar_to_jd(2026, 10, 28).plus_seconds(43200);
    s.check("TT − TAI = 32.184 s", ts.convert(tai, Scale::TAI, Scale::TT).minus(tai) * DAY - 32.184, 1e-9, "s");
    s.check("TAI − GPS = 19 s", tai.minus(ts.convert(tai, Scale::TAI, Scale::GPS)) * DAY - 19, 1e-9, "s");
    s.check("TT − UTC = 69.184 s in 2026", ts.convert(tai, Scale::TAI, Scale::TT).minus(ts.convert(tai, Scale::TAI, Scale::UTC)) * DAY - 69.184, 1e-9, "s");
  }

  // ---- date and time, with leap seconds --------------------------------------------------------
  {
    double worst = 0;
    bool fields = true;
    for (const auto& r : ref["dtf2d"]) {
      const Scale sc = parse_scale(r.fields[0]);
      const DateTime d{r.integer(1), r.integer(2), r.integer(3), r.integer(4), r.integer(5), r.number(6)};
      worst = std::max(worst, days_apart(from_datetime(d, sc, ts), {r.number(7), r.number(8)}));
    }
    for (const auto& r : ref["d2dtf"]) {
      const Scale sc = parse_scale(r.fields[0]);
      const DateTime d = to_datetime({r.number(2), r.number(3)}, sc, ts, 5);
      const double sec = r.integer(9) + r.integer(10) * 1e-5;
      fields = fields && d.year == r.integer(4) && d.month == r.integer(5) && d.day == r.integer(6) && d.hour == r.integer(7) &&
               d.minute == r.integer(8) && std::fabs(d.second - sec) < 5e-6;
    }
    s.check("calendar date and time to JD against ERFA, with 23:59:60", worst * DAY, 1e-9, "s");
    s.require("JD to calendar date and time against ERFA, with 23:59:60", fields && !ref["d2dtf"].empty());
    const DateTime leap = to_datetime({2400000.5, 49533.99999}, Scale::UTC, ts, 5);
    s.require("1994-06-30 23:59:60.13599 UTC (SOFA test)", leap.hour == 23 && leap.minute == 59 && std::fabs(leap.second - 60.13599) < 5e-6);
    const auto [t1, sc1] = parse_time("2016-12-31T23:59:60.5Z", ts);
    s.require("ISO 8601 with Z is UTC", sc1 == Scale::UTC);
    s.require("the leap second written back as 23:59:60.500Z", format_iso(t1, Scale::UTC, ts) == "2016-12-31T23:59:60.500Z");
    const JulianDate tai = ts.convert(t1, Scale::UTC, Scale::TAI);
    s.require("23:59:60.5 UTC is 00:00:36.5 TAI the next day", format_iso(tai, Scale::TAI, ts, 1) == "2017-01-01T00:00:36.5 TAI");
    const auto [t2, sc2] = parse_time("2000-01-01 12:00 TT", ts);
    s.require("\"2000-01-01 12:00 TT\" is J2000.0", sc2 == Scale::TT && std::fabs(t2.value() - J2000) < 1e-12);
    const auto [t3, sc3] = parse_time("MJD 60000.5 UTC", ts);
    s.require("\"MJD 60000.5 UTC\" read", sc3 == Scale::UTC && std::fabs(t3.mjd() - 60000.5) < 1e-12);
    bool refused = false;
    try {
      parse_time("2026-10-28T12:00:00", ts);
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("a time without a scale is refused", refused);
  }

  // ---- Earth orientation ------------------------------------------------------------------------
  {
    const eo::Series eop = eo::Series::parse_finals2000a(FINALS_2017);
    s.require("finals2000A: five rows read, the empty prediction skipped", eop.rows().size() == 5);
    const eo::Values v = eop.at(57753.0);
    s.check("finals2000A columns: x_p on 2016-12-31", v.xp - 0.081400, 1e-12, "arcsec");
    s.check("finals2000A columns: y_p", v.yp - 0.263094, 1e-12, "arcsec");
    s.check("finals2000A columns: UT1 − UTC", v.dut1 + 0.4077601, 1e-12, "s");
    s.check("finals2000A columns: LOD", v.lod - 0.8842e-3, 1e-12, "s");
    s.check("finals2000A columns: dX", v.dx - 0.025e-3, 1e-12, "arcsec");
    s.check("finals2000A columns: dY", v.dy + 0.169e-3, 1e-12, "arcsec");
    // Halfway through the leap-second day, UT1 − UTC is interpolated as UT1 − TAI.
    const eo::Values half = eop.at(57753.5);
    s.check("UT1 − UTC halfway through the leap second's day, without the jump",
            half.dut1 - 0.5 * (-0.4077601 + (0.5912821 - 1.0)), 1e-12, "s");
    s.check("UT1 − UTC after the leap second", eop.at(57754.0).dut1 - 0.5912821, 1e-12, "s");
    const TimeScales with(LeapSeconds::built_in(), &eop, Ut1Policy::refuse);
    const JulianDate utc = calendar_to_jd(2016, 12, 31).plus_days(0.25);
    const Instant i = with.instant(utc, Scale::UTC);
    s.check("UT1 from the series at 2016-12-31 06:00 UTC",
            i.ut1.minus(utc) * DAY - (-0.4077601 + 0.25 * (0.5912821 - 1.0 + 0.4077601)), 1e-9, "s");
    s.require("an instant says its UT1 was measured", i.ut1_source == Ut1Source::measured);
    s.check("UT1 → UTC → UT1 with the series", days_apart(with.convert(with.convert(i.ut1, Scale::UT1, Scale::UTC), Scale::UTC, Scale::UT1), i.ut1) * DAY, 1e-9, "s");
    bool refused = false;
    try {
      (void)with.instant(calendar_to_jd(2030, 1, 1), Scale::UTC);
    } catch (const std::domain_error&) {
      refused = true;
    }
    s.require("UT1 beyond the series is refused under the policy of refusal", refused);
    const Instant far = ts.instant(calendar_to_jd(1900, 6, 1), Scale::TT);
    s.require("an instant in 1900 has no UTC and takes UT1 from ΔT", !far.has_utc && far.ut1_source == Ut1Source::delta_t);
    s.check("ΔT in 1900 by Espenak and Meeus is −2.79 s", far.tt.minus(far.ut1) * DAY - delta_t_espenak_meeus(1900 + 5.5 / 12), 1e-9, "s");
  }
  return s.done();
}
