"""Reference values for the meteorology domain's U.S. Standard Atmosphere 1976 (atmosphere/standard.hpp).

Three references, written to tests/fixtures/atmosphere-reference.txt:

1. The document's own tables. U.S. Standard Atmosphere, 1976 (NOAA-S/T 76-1562; NASA-TM-X-74335,
   https://ntrs.nasa.gov/citations/19770009539), Part 4, Tables I-III, geometric-altitude pages, at
   thirty heights from -1 km to 86 km. The PDF is a scan without a usable text layer, so the values
   below were transcribed by eye from its pages (the copy whose SHA-256 is given in SOURCE_SHA256,
   rendered at 150 dpi); each entry keeps the printed digits, and the fixture carries half a unit in
   the last printed digit as its rounding. Below 86 km the tables were computed without the ratios
   M/M0 of Table 8 (Sec. 1.2.4, p. 9), and the document says how to correct them: T and the mean free
   path are multiplied by M/M0, the number density and the collision frequency divided by it, and the
   viscosities recomputed. This script applies that correction at 85 km (M/M0 = 0.999694, Table 8)
   and omits the viscosities there; pressure, density, gravity, the scale height, the particle speed
   and the speed of sound need none.
2. fluids (Caleb Bell, MIT licence), its ATMOSPHERE_1976 class, an implementation of the 1976 below
   86 km with the document's constants; it keeps T = T_M (no Table 8), so its T is the
   molecular-scale temperature, and its viscosity and conductivity are compared below 80 km only.
3. ambiance (Aaron Dettmann, Apache-2.0), an implementation of the ICAO Standard Atmosphere 1993 to
   80 km' (81.02 km), which shares the 1976's temperature profile below 80 km but rounds its base
   pressures to six figures and uses ICAO's constants (R = 287.052 87 J kg-1 K-1,
   N_A = 6.022 57 x 10^23 mol-1); compared to 80 km.

    python domains/meteorology/tests/make_atmosphere_fixtures.py
"""

from __future__ import annotations

import datetime
from decimal import Decimal
from importlib.metadata import version
from pathlib import Path

import ambiance
import numpy as np
from fluids.atmosphere import ATMOSPHERE_1976

OUT = Path(__file__).resolve().parent / "fixtures" / "atmosphere-reference.txt"
SOURCE_URL = "https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf"
SOURCE_SHA256 = "23eb6b68346d4abc6071fc598872d51c40cf8c613d1a5d1006c6261bcd09dcea"

# Table I (pp. 53-68): Z (m) -> T (K), P (mb), rho (kg m-3), with the document page.
TABLE_I = [
    (53, -1000, "294.651", "1.1393e3", "1.3470e0"), (53, 0, "288.150", "1.01325e3", "1.2250e0"),
    (53, 1000, "281.651", "8.9876e2", "1.1117e0"), (53, 2000, "275.154", "7.9501e2", "1.0066e0"),
    (55, 3000, "268.659", "7.0121e2", "9.0925e-1"), (55, 5000, "255.676", "5.4048e2", "7.3643e-1"),
    (57, 8000, "236.215", "3.5651e2", "5.2579e-1"), (57, 10000, "223.252", "2.6499e2", "4.1351e-1"),
    (59, 11000, "216.774", "2.2699e2", "3.6480e-1"), (59, 15000, "216.650", "1.2111e2", "1.9476e-1"),
    (59, 18000, "216.650", "7.5652e1", "1.2165e-1"), (61, 20000, "216.650", "5.5293e1", "8.8910e-2"),
    (61, 25000, "221.552", "2.5492e1", "4.0084e-2"), (63, 30000, "226.509", "1.1970e1", "1.8410e-2"),
    (63, 32000, "228.490", "8.8906e0", "1.3555e-2"), (63, 36000, "239.282", "4.9852e0", "7.2579e-3"),
    (65, 40000, "250.350", "2.8714e0", "3.9957e-3"), (65, 44000, "261.403", "1.6949e0", "2.2589e-3"),
    (65, 48000, "270.650", "1.0229e0", "1.3167e-3"), (65, 50000, "270.650", "7.9779e-1", "1.0269e-3"),
    (65, 55000, "260.771", "4.2525e-1", "5.6810e-4"), (67, 60000, "247.021", "2.1958e-1", "3.0968e-4"),
    (67, 65000, "233.292", "1.0929e-1", "1.6321e-4"), (67, 70000, "219.585", "5.2209e-2", "8.2829e-5"),
    (67, 75000, "208.399", "2.3881e-2", "3.9921e-5"), (67, 80000, "198.639", "1.0524e-2", "1.8458e-5"),
    (67, 85000, "188.893", "4.4568e-3", "8.2196e-6"), (68, 86000, "186.87", "3.7338e-3", "6.958e-6"),
]
# T_M at 86 km (Table I, p. 68); below 86 km the printed T_M equals the printed T.
TABLE_I_TM_86 = (68, 86000, "186.95")
# Table I's geopotential-altitude pages (pp. 52, 60): H (m') -> T (K), P (mb), rho (kg m-3). Remark: the
# printed pressure in mb is truncated, not rounded (226.9996 mb is printed 226.99 at Z = 11 km), as the
# P/P0 column beside it shows; the checks hold it accordingly.
TABLE_I_GEOPOTENTIAL = [
    (52, -1000, "294.650", "1.1392e3", "1.3470e0"), (52, 500, "284.900", "9.5460e2", "1.1673e0"),
    (52, 1000, "281.650", "8.9874e2", "1.1116e0"), (52, 1500, "278.400", "8.4555e2", "1.0581e0"),
    (52, 2000, "275.150", "7.9495e2", "1.0065e0"), (52, 2500, "271.900", "7.4682e2", "9.5686e-1"),
    (60, 19000, "216.650", "6.4100e1", "1.0307e-1"), (60, 23000, "219.650", "3.4224e1", "5.4280e-2"),
]
# Table II (pp. 77-92): Z -> g (m s-2), H_P (m), n (m-3), V (m s-1), nu (s-1), L (m).
TABLE_II = [
    (77, -1000, "9.8097", "8622.1", "2.8007e25", "464.09", "7.6934e9", "6.0324e-8"),
    (77, 0, "9.8066", "8434.5", "2.5470e25", "458.94", "6.9189e9", "6.6332e-8"),
    (77, 1000, "9.8036", "8246.9", "2.3113e25", "453.74", "6.2075e9", "7.3095e-8"),
    (77, 2000, "9.8005", "8059.2", "2.0928e25", "448.48", "5.5554e9", "8.0728e-8"),
    (83, 11000, "9.7728", "6367.2", "7.5848e24", "398.07", "1.7871e9", "2.2274e-7"),
    (83, 15000, "9.7605", "6371.6", "4.0493e24", "397.95", "9.5380e8", "4.1723e-7"),
    (83, 18000, "9.7513", "6377.6", "2.5292e24", "397.95", "5.9576e8", "6.6797e-7"),
    (85, 20000, "9.7452", "6381.6", "1.8486e24", "397.95", "4.3543e8", "9.1393e-7"),
    (85, 25000, "9.7300", "6536.2", "8.3341e23", "402.43", "1.9852e8", "2.0272e-6"),
    (87, 30000, "9.7147", "6692.9", "3.8278e23", "406.91", "9.2192e7", "4.4137e-6"),
    (87, 32000, "9.7087", "6755.7", "2.8183e23", "408.68", "6.8175e7", "5.9945e-6"),
    (87, 36000, "9.6965", "7083.7", "1.5090e23", "418.22", "3.7356e7", "1.1196e-5"),
    (89, 40000, "9.6844", "7420.6", "8.3077e22", "427.78", "2.1036e7", "2.0336e-5"),
    (89, 50000, "9.6542", "8047.4", "2.1351e22", "444.79", "5.6210e6", "7.9130e-5"),
    (89, 55000, "9.6391", "7765.8", "1.1812e22", "436.60", "3.0524e6", "1.4303e-4"),
    (91, 60000, "9.6241", "7367.8", "6.4387e21", "424.93", "1.6195e6", "2.6239e-4"),
    (91, 70000, "9.5942", "6569.9", "1.7222e21", "400.64", "4.0839e5", "9.8102e-4"),
    (91, 80000, "9.5644", "5961.7", "3.8378e20", "381.05", "8.6559e4", "4.4022e-3"),
    (91, 85000, "9.5496", "5678.0", "1.7090e20", "371.59", "3.7588e4", "9.8858e-3"),
    (92, 86000, "9.5466", "5621.", "1.447e20", "369.7", "3.17e4", "1.17e-2"),
]
# Table III (pp. 101-115): Z -> C_s (m s-1), mu (N s m-2), eta (m2 s-1).
TABLE_III = [
    (101, -1000, "344.11", "1.8206e-5", "1.3516e-5"), (101, 0, "340.29", "1.7894e-5", "1.4607e-5"),
    (101, 1000, "336.43", "1.7579e-5", "1.5813e-5"), (101, 2000, "332.53", "1.7260e-5", "1.7147e-5"),
    (107, 11000, "295.15", "1.4223e-5", "3.8988e-5"), (107, 15000, "295.07", "1.4216e-5", "7.2995e-5"),
    (107, 18000, "295.07", "1.4216e-5", "1.1686e-4"), (109, 20000, "295.07", "1.4216e-5", "1.5989e-4"),
    (109, 25000, "298.39", "1.4484e-5", "3.6135e-4"), (111, 30000, "301.71", "1.4753e-5", "8.0134e-4"),
    (111, 32000, "303.02", "1.4859e-5", "1.0962e-3"), (111, 36000, "310.10", "1.5433e-5", "2.1264e-3"),
    (113, 40000, "317.19", "1.6009e-5", "4.0066e-3"), (113, 44000, "324.12", "1.6573e-5", "7.3371e-3"),
    (113, 48000, "329.80", "1.7037e-5", "1.2939e-2"), (113, 50000, "329.80", "1.7037e-5", "1.6591e-2"),
    (113, 55000, "323.72", "1.6541e-5", "2.9117e-2"), (115, 60000, "315.07", "1.5837e-5", "5.1141e-2"),
    (115, 65000, "306.19", "1.5116e-5", "9.2617e-2"), (115, 70000, "297.06", "1.4377e-5", "1.7357e-1"),
    (115, 75000, "289.40", "1.3759e-5", "3.4465e-1"), (115, 80000, "282.54", "1.3208e-5", "7.1557e-1"),
    (115, 85000, "275.52", "1.2647e-5", "1.5386e0"),
]
RATIO_85KM = 0.999694  # M/M0 at Z = 85.0 km, Table 8 (p. 9)


def printed(text: str) -> tuple[float, float]:
    """A printed number and half a unit in its last printed digit."""
    d = Decimal(text)
    return float(d), 0.5 * float(Decimal(1).scaleb(d.as_tuple().exponent))


def table_rows() -> list[str]:
    rows = []

    def put(quantity, page, z, text, factor=1.0):
        v, h = printed(text)
        rows.append(f"table {quantity} {z} {v * factor!r} {h * factor!r} {page}")

    for page, z, t, p, rho in TABLE_I:
        put("T", page, z, t, RATIO_85KM if z == 85000 else 1.0)
        put("P", page, z, p, 100.0)  # mb -> Pa
        put("rho", page, z, rho)
    put("TM", TABLE_I_TM_86[0], TABLE_I_TM_86[1], TABLE_I_TM_86[2])
    for page, z, g, hp, n, v, nu, ell in TABLE_II:
        c = RATIO_85KM if z == 85000 else 1.0
        put("g", page, z, g)
        put("Hp", page, z, hp)
        put("n", page, z, n, 1 / c)
        put("V", page, z, v)
        put("nu", page, z, nu, 1 / c)
        put("L", page, z, ell, c)
    for page, z, cs, mu, eta in TABLE_III:
        put("Cs", page, z, cs)
        if z < 85000:
            put("mu", page, z, mu)
            put("eta", page, z, eta)
    for page, h, t, p, rho in TABLE_I_GEOPOTENTIAL:
        for quantity, text, factor in (("T", t, 1.0), ("P", p, 100.0), ("rho", rho, 1.0)):
            v, half = printed(text)
            rows.append(f"tableH {quantity} {h} {v * factor!r} {half * factor!r} {page}")
    return rows


def main() -> None:
    rng = np.random.default_rng(20261028)
    r0 = 6356766.0
    boundaries = [r0 * h / (r0 - h) for h in (11000, 20000, 32000, 47000, 51000, 71000)]
    heights = sorted(set([float(z) for z in np.arange(-4900.0, 86000.0 + 1, 500.0)] + [86000.0]
                         + [b + d for b in boundaries for d in (-0.5, 0.5)]
                         + [float(z) for z in rng.uniform(-4900, 86000, 120)]))
    rows = [
        f"# U.S. Standard Atmosphere 1976 references; computed {datetime.date.today().isoformat()} by domains/meteorology/tests/make_atmosphere_fixtures.py.",
        f"# table: the document's Tables I-III (NOAA-S/T 76-1562, {SOURCE_URL}, SHA-256 {SOURCE_SHA256}), transcribed by eye;",
        "#   quantity Z (m) value half-unit-of-the-last-printed-digit document-page; at 85 km corrected by Table 8 as Sec. 1.2.4 directs.",
        "# tableH: Table I's geopotential-altitude pages, the same with H (m') in place of Z.",
        f"# fluids: fluids {version("fluids")} ATMOSPHERE_1976 (T is T_M): Z T_M P rho C_s mu k g (SI).",
        f"# ambiance: ambiance {version("ambiance")} (ICAO 1993): Z H T P rho C_s mu eta n L V nu g H_P (SI).",
    ]
    rows += table_rows()
    for z in heights:
        a = ATMOSPHERE_1976(z)
        rows.append(f"fluids {z!r} {a.T!r} {a.P!r} {a.rho!r} {a.v_sonic!r} {a.mu!r} {a.k!r} {a.g!r}")
    for z in heights:
        if z > 80000:
            continue
        s = ambiance.Atmosphere(z)
        vals = [s.H, s.temperature, s.pressure, s.density, s.speed_of_sound, s.dynamic_viscosity, s.kinematic_viscosity,
                s.number_density, s.mean_free_path, s.mean_particle_speed, s.collision_frequency, s.grav_accel, s.pressure_scale_height]
        rows.append("ambiance " + repr(z) + " " + " ".join(repr(float(v[0])) for v in vals))
    OUT.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"wrote {OUT}: {len(table_rows())} table entries, {len(heights)} heights")


if __name__ == "__main__":
    main()
