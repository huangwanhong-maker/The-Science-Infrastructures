"""Reference values for the geodesy domain's normal gravity (gravity/normal.hpp).

Writes tests/fixtures/gravity-reference.txt with three kinds of line:

  published <ellipsoid> <name> <value> <half-unit>
      The derived constants as the defining documents print them: for WGS 84, NGA.STND.0036_1.0.0_WGS84
      (2014), Tables 3.5 and 3.6; for GRS 80, Moritz's "Geodetic Reference System 1980" (Bull. Geod. 54,
      395, 1980; J. Geodesy 74, 128, 2000), section 4, in the copy at MIT. Both are fetched, and every
      value below is required to appear, digit for digit, in the document's own text layer before it is
      written; half-unit is half a unit of its last printed digit.
  boule <ellipsoid> <name> <value>
      The same constants from boule (Fatiando a Terra, BSD-3-Clause), an independent implementation of the
      level ellipsoid; boule's GRS 80 is built from the printed 1/f = 298.257 222 101, not from J2.
  field <ellipsoid> <lat (deg)> <lon (deg)> <h (m)> <gamma> <U> <u> <beta (deg)> <gamma_beta> <down> <north>
      Normal gravity off the ellipsoid: gamma = -gamma_u, boule's normal_gravity (the closed formula of
      Lakshmanan 1991 and Li and Goetze 2001, Geophysics 66, 1660, which is the u-component, NGA eq. 4-5);
      U, boule's normal_gravity_potential (Hofmann-Wellenhof and Moritz 2006, eq. 2-126); u and beta,
      boule's ellipsoidal-harmonic coordinates of the point; gamma_beta = (1/(w sqrt(u^2 + E^2))) dU/dbeta,
      the gradient of boule's U along beta (w of NGA eq. 4-10, the metric of the coordinates), by a
      fourth-order central difference with an arc of 1 km, at fixed u, where boule's U is a smooth
      function of beta (error below 1e-10 m s-2); and the components down the geodetic normal and north
      in the meridian, down = -dU/dh and north = (1/(M + h)) dU/dphi, by the same differences of boule's
      U in geodetic coordinates with steps of 5 km. boule sums q(u) in closed form, whose two terms
      cancel to a part in 10^5, so its U carries rounding of some 1e-13 of itself as u varies; these two
      are good to some 2e-9 m s-2, and serve to check the projection on the local frame.

    python domains/geodesy/tests/make_gravity_fixtures.py
"""

from __future__ import annotations

import datetime
import hashlib
import tempfile
import urllib.request
import warnings
from decimal import Decimal
from importlib.metadata import version
from pathlib import Path

import boule
import numpy as np
import pypdf

OUT = Path(__file__).resolve().parent / "fixtures" / "gravity-reference.txt"
NGA_URL = "https://earth-info.nga.mil/php/download.php?file=coord-wgs84"
GRS80_URL = "http://geoweb.mit.edu/~tah/12.221_2005/grs80_corr.pdf"

# As NGA.STND.0036_1.0.0_WGS84 prints them (Tables 3.5, 3.6).
WGS84_PUBLISHED = [
    ("f", "3.3528106647475 × 10-03"), ("b", "6356752.314 2"), ("e2", "6.694379990141 × 10-03"),
    ("second_e", "8.2094437949696 × 10-02"), ("second_e2", "6.739496742276 × 10-03"), ("E", "5.2185400842339 × 10+05"),
    ("c", "6399593.6258"), ("C20", "-4.84166774985 × 10-04"), ("J2", "1.082629821313 × 10-03"),
    ("U0", "6.26368517146 × 10+07"), ("gamma_e", "9.7803253359"), ("gamma_p", "9.8321849379"),
    ("gamma_mean", "9.7976432223"), ("k", "1.931852652458 × 10-03"), ("m", "3.449786506841 × 10-03"),
]
# As Moritz prints them (section 4).
GRS80_PUBLISHED = [
    ("b", "6 356 752.3141"), ("E", "521 854.0097"), ("c", "6 399 593.6259"), ("e2", "0.006 694 380 022 90"),
    ("second_e2", "0.006 739 496 775 48"), ("f", "0.003 352 810 681 18"), ("inverse_f", "298.257 222 101"),
    ("U0", "6 263 686.0850 x 10"), ("J4", "-0.000 002 370 912 22"), ("J6", "0.000 000 006 083 47"),
    ("J8", "-0.000 000 000 014 27"), ("m", "0.003 449 786 003 08"), ("gamma_e", "9.780 326 7715"),
    ("gamma_p", "9.832 186 3685"), ("f_star", "0.005 302 440 112"), ("k", "0.001 931 851 353"),
    ("gamma_mean", "9.797 644 656"), ("gamma_45", "9.806 199 203"),
]


def squeeze(text: str) -> str:
    return "".join(text.replace("/nobreakspace", " ").split())


def printed(text: str) -> tuple[float, float]:
    """The value of a printed number ('6.26368517146 × 10+07', '6 263 686.0850 x 10') and half its last unit."""
    s = squeeze(text)
    exponent = 0
    for times in ("×10", "x10"):
        if times in s:
            s, rest = s.split(times)
            exponent = int(rest) if rest else 1
    d = Decimal(s)
    return float(d.scaleb(exponent)), 0.5 * float(Decimal(1).scaleb(d.as_tuple().exponent + exponent))


def document(url: str, folder: Path) -> tuple[str, str]:
    """A PDF's text layer, whitespace removed, and the SHA-256 of the file."""
    path = folder / "document.pdf"
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (science-infrastructures fixture generator)"})
    with urllib.request.urlopen(request, timeout=120) as response:
        path.write_bytes(response.read())
    text = "".join(page.extract_text() or "" for page in pypdf.PdfReader(str(path)).pages)
    return squeeze(text), hashlib.sha256(path.read_bytes()).hexdigest()


def published(name: str, url: str, entries: list[tuple[str, str]], rows: list[str]) -> str:
    with tempfile.TemporaryDirectory() as folder:
        text, digest = document(url, Path(folder))
    for key, as_printed in entries:
        if squeeze(as_printed) not in text:
            raise SystemExit(f"{name} {key}: '{as_printed}' is not in the document at {url}")
        value, half = printed(as_printed)
        rows.append(f"published {name} {key} {value!r} {half!r}")
    return digest


def main() -> None:
    rows: list[str] = []
    nga = published("WGS84", NGA_URL, WGS84_PUBLISHED, rows)
    grs = published("GRS80", GRS80_URL, GRS80_PUBLISHED, rows)
    ellipsoids = {"WGS84": boule.WGS84, "GRS80": boule.GRS80}
    for name, e in ellipsoids.items():
        for key, value in [("b", e.semiminor_axis), ("E", e.linear_eccentricity), ("e2", e.first_eccentricity**2),
                           ("second_e2", e.second_eccentricity**2), ("m", e._emm), ("gamma_e", e.gravity_equator),
                           ("gamma_p", e.gravity_pole), ("U0", e.reference_normal_gravity_potential)]:
            rows.append(f"boule {name} {key} {float(value)!r}")

    rng = np.random.default_rng(20261028)
    latitudes = [-90.0, -80.0, -60.0, -45.0, -30.0, -10.0, 0.0, 5.0, 15.0, 30.0, 40.0, 45.0, 50.0, 60.0, 75.0, 85.0, 90.0]
    latitudes += [float(x) for x in np.degrees(np.arcsin(rng.uniform(-1, 1, 20)))]
    heights = [0.0, 10.0, 100.0, 1000.0, 3000.0, 8848.0, 20000.0, 50000.0, 100000.0]

    def derivative(f, x, step):
        return (-f(x + 2 * step) + 8 * f(x + step) - 8 * f(x - step) + f(x - 2 * step)) / (12 * step)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # boule warns at negative heights, which the stencils at h = 0 reach
        for name, e in ellipsoids.items():
            big_e = e.linear_eccentricity
            U_geodetic = lambda lat, h: float(e.normal_gravity_potential((0.0, lat, h)))  # noqa: E731
            U_harmonic = lambda beta, u: float(e.normal_gravity_potential((0.0, beta, u), coordinate_system="ellipsoidal harmonic"))  # noqa: E731
            for lat in latitudes:
                for h in heights:
                    lon = float(rng.uniform(-180, 180))
                    gamma = float(e.normal_gravity((lon, lat, h), si_units=True))
                    potential = float(e.normal_gravity_potential((lon, lat, h)))
                    _, beta, u = (float(v) for v in e.geodetic_to_ellipsoidal_harmonic((lon, lat, h)))
                    sb2 = np.sin(np.radians(beta)) ** 2
                    w = np.sqrt((u * u + big_e * big_e * sb2) / (u * u + big_e * big_e))
                    arc = w * np.sqrt(u * u + big_e * big_e)  # metres per radian of beta
                    gamma_beta = float(derivative(lambda x: U_harmonic(x, u), beta, np.degrees(1000.0 / arc)) / np.radians(1.0) / arc)
                    down = float(-derivative(lambda x: U_geodetic(lat, x), h, 5000.0))
                    s2 = np.sin(np.radians(lat)) ** 2
                    radius = e.semimajor_axis * (1 - e.first_eccentricity**2) / (1 - e.first_eccentricity**2 * s2) ** 1.5 + h
                    north = float(derivative(lambda x: U_geodetic(x, h), lat, np.degrees(5000.0 / radius)) / np.radians(1.0) / radius)
                    rows.append(f"field {name} {lat!r} {lon!r} {h!r} {gamma!r} {potential!r} {u!r} {beta!r} {gamma_beta!r} {down!r} {north!r}")

    header = [
        f"# Normal gravity references; computed {datetime.date.today().isoformat()} by domains/geodesy/tests/make_gravity_fixtures.py.",
        f"# published WGS84: NGA.STND.0036_1.0.0_WGS84 (2014), Tables 3.5, 3.6, {NGA_URL}, SHA-256 {nga}.",
        f"# published GRS80: H. Moritz, Geodetic Reference System 1980, sec. 4, {GRS80_URL}, SHA-256 {grs}.",
        f"# boule {version('boule')} (numpy {np.__version__}): derived constants, and field rows:",
        "# field ellipsoid lat lon (deg) h (m) gamma=-gamma_u (m s-2) U (m2 s-2) u (m) beta (deg) gamma_beta down north (m s-2; the last three by differences of boule's U).",
    ]
    OUT.write_text("\n".join(header + rows) + "\n", encoding="utf-8")
    print(f"wrote {OUT}: {sum(r.startswith('published') for r in rows)} published, {sum(r.startswith('field') for r in rows)} field rows")


if __name__ == "__main__":
    main()
