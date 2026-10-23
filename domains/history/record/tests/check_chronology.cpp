// Checks of the chronology: Julian day numbers of both calendars against Meeus's worked
// examples and table (Astronomical Algorithms, 2nd ed., 1998, ch. 7) and against Python's
// proleptic Gregorian ordinals (retrieved 2026-10-24); the reform of 1582, the German change of
// 1700 and the Russian of 1918; Old Style's legal year from 25 March (Newton's birth and
// death); the inverses over the whole range; and the EDTF subset of codebook §4.

#include <history/record/chronology.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

using namespace history;

namespace {

bool throws(const char* text, Calendar c = Calendar::gregorian) {
  try {
    parse_date(text, c);
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  auto same = [&](const char* name, long got, long want) { s.require(name, got == want, std::to_string(got) + (got == want ? "" : ", expected " + std::to_string(want))); };

  // Julian day numbers (Meeus ch. 7; Python's date.toordinal() + 1721425 for the Gregorian).
  same("JDN: 1 January 2000, Gregorian (J2000.0 is its noon)", jdn_gregorian(2000, 1, 1), 2451545);
  same("JDN: 4 October 1957, Gregorian (Meeus ex. 7.a, JD 2436116.31)", jdn_gregorian(1957, 10, 4), 2436116);
  same("JDN: 27 January 333, Julian (Meeus ex. 7.b, JD 1842713.0)", jdn_julian(333, 1, 27), 1842713);
  same("JDN: 10 April 837, Julian (Meeus, JD 2026871.8)", jdn_julian(837, 4, 10), 2026872);
  same("JDN: 12 July −1000, Julian (Meeus, JD 1356001.0)", jdn_julian(-1000, 7, 12), 1356001);
  same("JDN: 1 January −4712, Julian, is day 0", jdn_julian(-4712, 1, 1), 0);
  same("JDN: 1 January 1600, Gregorian (Meeus, JD 2305447.5)", jdn_gregorian(1600, 1, 1), 2305448);
  same("JDN: 19 June 1987, Gregorian (Meeus, JD 2446966.0)", jdn_gregorian(1987, 6, 19), 2446966);
  same("JDN: 15 March 44 BC, Julian (the Ides of March)", jdn_julian(-43, 3, 15), 1705426);

  // The reforms.
  same("1582: Thursday 4 October (Julian) is followed by Friday 15 October (Gregorian)", jdn_gregorian(1582, 10, 15) - jdn_julian(1582, 10, 4), 1);
  same("1700: the Protestant German states pass from 18 February (Julian) to 1 March (Gregorian)", jdn_gregorian(1700, 3, 1) - jdn_julian(1700, 2, 18), 1);
  same("1700 is a leap year of the Julian calendar", days_in_month(Calendar::julian, 1700, 2), 29);
  same("1700 is not a leap year of the Gregorian calendar", days_in_month(Calendar::gregorian, 1700, 2), 28);
  same("1917: the October Revolution, 25 October Julian, was 7 November Gregorian", jdn_julian(1917, 10, 25), jdn_gregorian(1917, 11, 7));

  // Old Style: Julian, the legal year from 25 March.
  same("Old Style: Newton born 25 December 1642 is 4 January 1643 Gregorian", jdn_of(Calendar::julian_march, 1642, 12, 25), jdn_gregorian(1643, 1, 4));
  same("Old Style: Newton died 20 March 1726 is 31 March 1727 Gregorian", jdn_of(Calendar::julian_march, 1726, 3, 20), jdn_gregorian(1727, 3, 31));
  same("Old Style: 25 March 1726 begins the legal year 1726", jdn_of(Calendar::julian_march, 1726, 3, 25), jdn_julian(1726, 3, 25));

  // The historians' convention for a reader's dates (Meeus ch. 7): Julian before the reform, Gregorian from it.
  same("convention: 4 October 1582 is followed by 15 October 1582", jdn_historical(1582, 10, 15) - jdn_historical(1582, 10, 4), 1);
  same("convention: the Ides of March, 44 BC, are Julian", jdn_historical(-43, 3, 15), 1705426);
  same("convention: 1 January 2000 is Gregorian", jdn_historical(2000, 1, 1), 2451545);
  same("convention: the end of 338 BC is the day before the Julian 337 BC begins", jdn_historical(-337, 12, 31) + 1, jdn_julian(-336, 1, 1));
  same("convention: read in the proleptic Gregorian, 338 BC would end 5 days into the Julian 337 BC", jdn_gregorian(-337, 12, 31) - jdn_julian(-336, 1, 1), 4);
  {
    bool ok = true;
    for (long j = 1500000; j <= 2500000 && ok; j += 41) {
      const Ymd h = historical_of(j);
      ok = jdn_historical(h.year, h.month, h.day) == j;
    }
    s.require("convention: historical_of inverts jdn_historical, every 41st day from 1500000 to 2500000", ok);
  }
  s.check("convention: the decimal year of 1 July 338 BC (Julian) is −337 + 181/365", historical_year(jdn_julian(-337, 7, 1)) - (-337 + 181.0 / 365), 1e-12);
  same("convention: a reader's −0337-12-31 is Julian", parse_historical("-0337-12-31").last, jdn_julian(-337, 12, 31));
  same("convention: a reader's 1814-11-15 is Gregorian", parse_historical("1814-11-15").last, jdn_gregorian(1814, 11, 15));

  // The inverses, over the whole range.
  {
    bool ok = true;
    long bad = -1;
    for (long j = 0; j <= 2600000 && ok; j += 37) {
      const Ymd g = gregorian_of(j), u = julian_of(j);
      if (jdn_gregorian(g.year, g.month, g.day) != j || jdn_julian(u.year, u.month, u.day) != j) { ok = false; bad = j; }
    }
    s.require("the inverses of both calendars, every 37th day from 0 to 2,600,000", ok, ok ? "" : "fails at JDN " + std::to_string(bad));
  }
  s.require("Gregorian of JDN 2451545 is 2000-01-01", iso_gregorian(2451545) == "2000-01-01", iso_gregorian(2451545));
  s.check("decimal year of 2 October 1808 is 1808 + 275/366", decimal_year(jdn_gregorian(1808, 10, 2)) - (1808 + 275.0 / 366), 1e-12);

  // EDTF, the subset of codebook §4.
  {
    const Date d = parse_date("1785-02-12");
    s.require("EDTF: a day", d.precision == 3 && d.first == d.last && d.first == jdn_gregorian(1785, 2, 12) && !d.interval);
    const Date m = parse_date("1790-12~");
    s.require("EDTF: an approximate month", m.precision == 2 && m.approximate && !m.uncertain && m.first == jdn_gregorian(1790, 12, 1) && m.last == jdn_gregorian(1790, 12, 31));
    const Date y = parse_date("1795?");
    s.require("EDTF: an uncertain year", y.precision == 1 && y.uncertain && !y.approximate && y.first == jdn_gregorian(1795, 1, 1) && y.last == jdn_gregorian(1795, 12, 31));
    const Date b = parse_date("1797-03/1797-05");
    s.require("EDTF: an interval of months", b.interval && b.first == jdn_gregorian(1797, 3, 1) && b.last == jdn_gregorian(1797, 5, 31));
    const Date x = parse_date("190X");
    s.require("EDTF: unspecified digits of a year", x.first == jdn_gregorian(1900, 1, 1) && x.last == jdn_gregorian(1909, 12, 31));
    const Date both = parse_date("1804%");
    s.require("EDTF: % is uncertain and approximate", both.uncertain && both.approximate);
    const Date ides = parse_date("-0043-03-15", Calendar::julian);
    s.require("EDTF: a negative year, Julian", ides.first == 1705426, std::to_string(ides.first));
    const Date legal = parse_date("1715", Calendar::julian_march);
    s.require("EDTF: an Old Style year runs from 25 March to 24 March", legal.first == jdn_julian(1715, 3, 25) && legal.last == jdn_julian(1716, 3, 24));
    const Date feb = parse_date("1715-02", Calendar::julian_march);
    s.require("EDTF: February of the legal year 1715 is February 1716, a leap year", feb.first == jdn_julian(1716, 2, 1) && feb.last == jdn_julian(1716, 2, 29));
    s.require("EDTF: refused: month 13", throws("1795-13"));
    s.require("EDTF: refused: 30 February", throws("1795-02-30"));
    s.require("EDTF: accepted: 29 February 1700 Julian, refused Gregorian", !throws("1700-02-29", Calendar::julian) && throws("1700-02-29"));
    s.require("EDTF: refused: not a date", throws("about 1800") && throws("") && throws("180"));
    s.require("EDTF: refused: an interval that ends before it begins", throws("1800/1799"));
  }
  return s.done();
}
