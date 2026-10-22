"""Checks of observations: the units and the humidity formula against reference
values, the decoding of METAR groups, the products of METAR, NDBC and IBTrACS
when they have been built, the rolling window, and the sensor store, run in a
temporary directory so that the user's readings and index are untouched.

Run with the data library's other checks: node data/tests/run-tests.mjs
"""

from __future__ import annotations

import datetime as dt
import json
import math
import shutil
import sys
import tempfile
from pathlib import Path

from checks import PRODUCTS, ROOT, assert_, check, skip

import observation_set as obs  # noqa: E402
import product_index  # noqa: E402
import ibtracs_tracks  # noqa: E402

sys.path.insert(0, str(ROOT / "observation" / "sensors"))
import store  # noqa: E402

UTC = dt.timezone.utc


# ---- 1. The humidity formula -------------------------------------------------
# Saturation vapour pressure of water by IAPWS-95 (Wagner and Pruss 2002), as
# tabulated in the CRC Handbook of Chemistry and Physics ('Vapor pressure of
# water from 0 to 370 C'): 0.6117 kPa at 0.01 C, 1.2282 kPa at 10 C, 2.3393 kPa
# at 20 C, 4.2470 kPa at 30 C. Alduchov and Eskridge (1996) give the Magnus
# form's error as under about 0.4 per cent from -40 to +50 C.

for celsius, kpa in [(0.01, 0.6117), (10.0, 1.2282), (20.0, 2.3393), (30.0, 4.2470)]:
    e = obs.saturation_vapour_pressure(celsius)
    check(f"humidity: Magnus (AERK) saturation vapour pressure at {celsius:g} C against IAPWS-95, relative error",
          abs(e / (kpa * 1000.0) - 1.0) * 100.0, 0.4, "%")
# At 20 C with a dew point of 10 C, the reference pressures give 100 x 1.2282 / 2.3393 = 52.50 per cent.
rh = obs.relative_humidity(293.15, 283.15)
check("humidity: relative humidity at 20 C, dew point 10 C, against 52.50 % from IAPWS-95", abs(rh - 52.50), 0.15, "%")
check("humidity: equal temperature and dew point give 100 %", abs(obs.relative_humidity(280.0, 280.0) - 100.0), 1e-9, "%")


# ---- 2. Units ----------------------------------------------------------------
# By definition (BIPM SI Brochure, 2019; NIST SP 811, 2008, B.8); the standard
# atmosphere's 1013.25 hPa is 29.9213 inHg (NIST SP 811: 1 inHg = 3386.389 Pa).

check("units: 15 C is 288.15 K", abs(obs.kelvin(15.0) - 288.15), 1e-12, "K")
check("units: 10 kt is 5.144 m s-1 (1852 m per hour)", abs(10 * obs.KNOT - 5.144444), 1e-6, "m s-1")
check("units: 10 statute miles are 16 093.44 m", abs(10 * obs.STATUTE_MILE - 16093.44), 1e-9, "m")
check("units: 29.9213 inHg is the standard atmosphere, 101 325 Pa", abs(29.9213 * obs.INCH_OF_MERCURY - 101325.0), 1.0, "Pa")
check("units: 1 in of precipitation is 25.4 kg m-2", abs(obs.INCH * obs.WATER_DENSITY - 25.4), 1e-12, "kg m-2")
check("units: a store's 59 F is 288.15 K", abs(store.convert(59.0, "degree_F", "K") - 288.15), 1e-9, "K")
check("units: a store's 1013.25 hPa is 101 325 Pa", abs(store.convert(1013.25, "hPa", "Pa") - 101325.0), 1e-9, "Pa")
check("units: a store's 36 km/h is 10 m s-1", abs(store.convert(36.0, "km/h", "m s-1") - 10.0), 1e-12, "m s-1")
check("units: a store's 25 ug m-3 is 2.5e-8 kg m-3", abs(store.convert(25.0, "ug m-3", "kg m-3") - 2.5e-8), 1e-20, "kg m-3")


# ---- 3. METAR groups ---------------------------------------------------------

d, s, g = obs.metar_wind("METAR UUWW 221500Z 12005MPS CAVOK 08/01 Q1032 R06/000062 NOSIG")
assert_("metar: 12005MPS is 5 m s-1 from 120 degrees, without gusts", (d, s, g) == (120.0, 5.0, None), (d, s, g))
d, s, g = obs.metar_wind("METAR KJFK 221451Z 36016G23KT 10SM FEW250 15/M01 A3001 RMK AO2 PK WND 22028/1500 SLP163")
check("metar: 36016G23KT is 16 kt from the north with gusts of 23 kt",
      abs(s - 16 * obs.KNOT) + abs(g - 23 * obs.KNOT) + abs(d - 360.0), 1e-9, "")
assert_("metar: variable (VRB03KT) and calm (00000KT) winds have no direction",
        obs.metar_wind("METAR KRME 221504Z VRB03KT 1SM BR OVC002 06/05 A3018")[0] is None
        and obs.metar_wind("METAR KMTJ 221504Z AUTO 00000KT 2SM RA BR 09/09 A3034")[:2] == (None, 0.0))
assert_("metar: 9999 and CAVOK are 10 km; 0800 is 800 m; statute miles are left to the column",
        obs.metar_visibility("METAR RJTT 221500Z 04010KT 9999 FEW045 17/09 Q1026") == 10000.0
        and obs.metar_visibility("METAR UUWW 221500Z 12005MPS CAVOK 08/01 Q1032") == 10000.0
        and obs.metar_visibility("METAR EGLL 221450Z 24012KT 210V270 0800 FG 15/06 Q1009") == 800.0
        and obs.metar_visibility("METAR KJFK 221451Z 36016G23KT 10SM FEW250 15/M01 A3001") is None)
check("metar: Q1026 is 102 600 Pa", abs(obs.metar_qnh("METAR RJTT 221500Z 04010KT 9999 17/09 Q1026 NOSIG") - 102600.0), 1e-9, "Pa")
check("metar: A3001 is 30.01 inHg, 101 625 Pa",
      abs(obs.metar_qnh("METAR KJFK 221451Z 36016G23KT 10SM 15/M01 A3001 RMK SLP163") - 30.01 * 3386.389), 1e-6, "Pa")


# ---- 4. The rolling window ---------------------------------------------------

temporary = Path(tempfile.mkdtemp(prefix="science-observations-"))
real_out = product_index.OUT
try:
    product_index.OUT = temporary / "products"
    variables = [obs.variable("air_temperature", "Air temperature", "K", "2 m", 2)]

    def reading(ident, hour, value):
        t = dt.datetime(2026, 10, 22, 0, tzinfo=UTC) + dt.timedelta(hours=hour)
        return {"id": ident, "name": ident, "lon": 1.0, "lat": 2.0, "elevation": 3.0, "time": obs.iso(t),
                "values": {"air_temperature": value}}

    end = dt.datetime(2026, 10, 22, 0, tzinfo=UTC)
    obs.write("test", "Test", "meteorology", {"dataset": "test"}, variables,
              [reading("A", -50, 270.0), reading("A", -1, 280.0), reading("B", -1, 281.0)], hours=48, end=end)
    path = obs.write("test", "Test", "meteorology", {"dataset": "test"}, variables,
                     [reading("A", 0, 282.0), reading("B", -1, 285.0), reading("C", 0, 400.0)], hours=48,
                     end=end + dt.timedelta(hours=1))
    product = json.loads(path.read_text(encoding="utf-8"))
    rows = [(r[0], r[5], r[6]) for r in product["rows"]]
    assert_("window: rows merged once per station and time, the newer retrieval winning, by station then time",
            rows == [("A", "2026-10-21T23:00:00Z", 280.0), ("A", "2026-10-22T00:00:00Z", 282.0),
                     ("B", "2026-10-21T23:00:00Z", 285.0), ("C", "2026-10-22T00:00:00Z", None)], rows)
    assert_("window: a row older than the window's 48 hours is dropped, and 400 K is nulled as a gross error",
            product["window"] == {"start": "2026-10-20T01:00:00Z", "end": "2026-10-22T01:00:00Z"}
            and product["screened"]["nulled"] == {"air_temperature": 1}, product["window"])
    entry = next(e for e in json.loads((product_index.OUT / "index.json").read_text(encoding="utf-8"))["observations"]
                 if e["name"] == "test")
    assert_("window: registered with its count of stations, rows and variables",
            entry["count"] == 3 and entry["rows"] == 4 and entry["variables"] == ["air_temperature"], entry)
finally:
    product_index.OUT = real_out


# ---- 5. The sensor store -----------------------------------------------------

real_store = store.STORE
try:
    store.STORE = temporary / "observations"
    product_index.OUT = temporary / "products"
    refused = []
    for bad in ["../garden", "Garden", "garden_station", "a--b", ""]:
        try:
            store.define(bad, "x", ["air_temperature"])
            refused.append(bad)
        except ValueError:
            pass
    assert_("store: names that could leave the store, or are not lowercase and hyphenated, are refused",
            not refused, refused)
    desc = store.define("garden", "Garden station", ["air_temperature:degree_C", "relative_humidity:%",
                                                     {"name": "battery_voltage", "units": "V"}],
                        [{"id": "g1", "name": "Lawn", "lon": 139.70, "lat": 35.68, "elevation": 40}])
    assert_("store: define keeps a CF variable in SI with its input units, and a custom one as declared",
            desc["variables"][0] == {"name": "air_temperature", "units": "K", "long_name": "Air temperature",
                                     "standard_name": "air_temperature", "input_units": "degree_C"}
            and desc["variables"][2]["standard_name"] is None, desc["variables"])
    try:
        store.define("garden", "Garden station", ["air_temperature:furlong"])
        assert_("store: units without a known conversion are refused", False)
    except ValueError:
        assert_("store: units without a known conversion are refused", True)

    heard = []
    unsubscribe = store.subscribe(lambda name, readings: heard.append((name, len(readings))))
    now = dt.datetime.now(UTC).replace(microsecond=0)
    t = lambda minutes: obs.iso(now - dt.timedelta(minutes=minutes))  # noqa: E731
    good = [
        {"station": "g1", "time": t(30), "values": {"air_temperature": 15.0, "relative_humidity": 60.0}},
        {"station": "g2", "lon": 139.71, "lat": 35.69, "time": int((now - dt.timedelta(minutes=20)).timestamp() * 1000),
         "values": {"air_temperature": 16.5, "battery_voltage": 3.7}},
        {"station": "boat", "lon": -170.0, "lat": 10.0, "time": t(10), "values": {"air_temperature": 27.0}},
    ]
    bad = [
        {"station": "g1", "time": t(5), "values": {"wind_speed": 3.0}},                     # not a variable of the source
        {"station": "g1", "time": t(5), "values": {"air_temperature": float("nan")}},       # not finite
        {"station": "g1", "time": t(5), "values": {"air_temperature": 95.0}},               # 368 K: implausible
        {"station": "g1", "time": 1761000000, "values": {"air_temperature": 15.0}},         # seconds, read as ms: 1970
        {"station": "g1", "time": obs.iso(now + dt.timedelta(hours=2)), "values": {"air_temperature": 15.0}},
        {"station": "g3", "time": t(5), "values": {"air_temperature": 15.0}},               # a new station without a place
        {"station": "g4", "lon": 200.0, "lat": 10.0, "time": t(5), "values": {"air_temperature": 15.0}},
        {"station": "../x", "lon": 1.0, "lat": 1.0, "time": t(5), "values": {"air_temperature": 15.0}},
    ]
    rejected = []
    accepted = store.append("garden", good + bad, rejected)
    assert_("store: append accepts the 3 good readings and refuses the 8 bad ones",
            accepted == 3 and [i for i, _ in rejected] == list(range(3, 11)), rejected)
    kept = store.read("garden")
    check("store: 15 C arrives as 288.15 K", abs(kept[0]["values"]["air_temperature"] - 288.15), 1e-9, "K")
    assert_("store: a subscriber hears each append", heard == [("garden", 3)], heard)
    unsubscribe()
    store.append("garden", [{"station": "g1", "time": t(1), "values": {"air_temperature": 14.0}}])
    assert_("store: an unsubscribed callback hears no more", heard == [("garden", 3)], heard)
    stations = {s["id"]: s for s in store.sources()[0]["stations"]}
    assert_("store: new stations are recorded with their places", set(stations) == {"g1", "g2", "boat"}
            and stations["g2"]["lon"] == 139.71, stations)
    assert_("store: read with a time window", [r["station"] for r in store.read("garden", since=t(25), until=t(5))]
            == ["g2", "boat"])
    assert_("store: read inside a box, across the antimeridian", [r["station"] for r in store.read("garden", bbox=(170, 0, -160, 20))]
            == ["boat"])
    assert_("store: read inside a box given as text, and the last n", [r["station"] for r in store.read(
            "garden", bbox="139,35,140,36", limit=2)] == ["g2", "g1"])
    path = store.build_product("garden", hours=48)
    product = json.loads(path.read_text(encoding="utf-8"))
    ok = (product["schema"] == "observations/1" and product["name"] == "sensor-garden"
          and product["fields"][:6] == ["id", "name", "lon", "lat", "elevation", "time"]
          and product["fields"][6:] == [v["name"] for v in product["variables"]]
          and all(len(r) == len(product["fields"]) for r in product["rows"]) and len(product["rows"]) == 4
          and all(k in product["source"] for k in ("dataset", "provider", "credit", "licence", "retrieved_at"))
          and all(r[5].endswith("Z") for r in product["rows"]))
    assert_("store: build_product writes a valid observations/1 product of the 4 readings", ok, product["fields"])
    first = next(r for r in product["rows"] if r[0] == "g1")
    assert_("store: the product names the station and keeps SI values",
            first[1] == "Lawn" and abs(first[6] - 288.15) < 1e-9, first)
    entry = next(e for e in json.loads((product_index.OUT / "index.json").read_text(encoding="utf-8"))["observations"]
                 if e["name"] == "sensor-garden")
    assert_("store: the product is registered as observations sensor-garden", entry["count"] == 3, entry)
finally:
    store.STORE = real_store
    product_index.OUT = real_out
    shutil.rmtree(temporary, ignore_errors=True)


# ---- 6. IBTrACS parsing ------------------------------------------------------
# An excerpt of two rows made in the form of ibtracs.*.list.v04r01.csv: the WMO
# agency's wind is preferred; where it gives none, outside its span of times,
# the US agency's; longitudes beyond 180 are reduced.

names = ["SID", "SEASON", "NUMBER", "BASIN", "SUBBASIN", "NAME", "ISO_TIME", "NATURE", "LAT", "LON", "WMO_WIND",
         "WMO_PRES", "WMO_AGENCY", "TRACK_TYPE", "DIST2LAND", "LANDFALL", "IFLAG", "USA_AGENCY", "USA_ATCF_ID",
         "USA_WIND", "USA_PRES", "TOKYO_WIND", "TOKYO_PRES"]
excerpt = [
    ["2024001N10170", "2024", "1", "WP", "MM", "TEST", "2024-01-01 00:00:00", "TS", "10.0", "185.0", " ", " ", " ",
     "main", "500", "500", "OP_____________", "jtwc_wp", "WP012024", "40", "995", " ", " "],
    ["2024001N10170", "2024", "1", "WP", "MM", "TEST", "2024-01-01 03:00:00", "TS", "10.2", "184.0", " ", " ", " ",
     "main", "500", "500", "PP_____________", "jtwc_wp", "WP012024", "42", "994", " ", " "],
    ["2024001N10170", "2024", "1", "WP", "MM", "TEST", "2024-01-01 06:00:00", "TS", "10.5", "183.0", "45", "990",
     "tokyo", "main", "500", "500", "OO_____________", "jtwc_wp", "WP012024", "55", "985", "45", "990"],
    ["2024001N10170", "2024", "1", "WP", "MM", "TEST", "2024-01-01 06:00:00", "TS", "11.5", "183.0", " ", " ", " ",
     "spur", "500", "500", "O______________", "jtwc_wp", "WP012024", "55", "985", " ", " "],
]
with tempfile.TemporaryDirectory() as folder:
    csv_path = Path(folder) / "ibtracs.csv"
    csv_path.write_text("\n".join(",".join(r) for r in [names, [" "] * len(names)] + excerpt) + "\n", encoding="utf-8")
    tracks = ibtracs_tracks.read(csv_path)
points = tracks[0]["points"] if tracks else []
assert_("ibtracs: interpolated times and spur tracks are left out", len(tracks) == 1 and len(points) == 2, points)
check("ibtracs: before the WMO agency reports, the US agency's 40 kt (1-minute) and 995 hPa",
      abs(points[0][5] - 40 * obs.KNOT) + abs(points[0][4] - 995.0), 0.005, "")   # winds kept to 0.01 m s-1
assert_("ibtracs: the WMO agency's 45 kt and 990 hPa are preferred, and named",
        abs(points[1][5] - 45 * obs.KNOT) < 0.01 and points[1][4] == 990.0 and points[1][8] == "tokyo"
        and points[0][8] == "jtwc_wp", points)
check("ibtracs: longitude 185 E is reduced to -175", abs(points[0][1] + 175.0), 1e-9, "degree")


# ---- 7. The products, when built ---------------------------------------------

def load(name: str) -> dict | None:
    path = PRODUCTS / name
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def column(product: dict, name: str) -> int:
    return product["fields"].index(name)


def times_ordered(product: dict) -> bool:
    by_station = {}
    for r in product["rows"]:
        by_station.setdefault(r[0], []).append(r[5])
    return all(ts == sorted(ts) and len(ts) == len(set(ts)) for ts in by_station.values())


metar = load("observations-metar.json")
if metar is None:
    skip("metar: observations-metar.json not built (fetch.py metar; observation_set.py metar)")
else:
    rows = metar["rows"]
    stations = {r[0] for r in rows}
    assert_("metar: more than 1000 stations report", len(stations) > 1000, f"{len(stations)} stations")
    ta = [r[column(metar, "air_temperature")] for r in rows if r[column(metar, "air_temperature")] is not None]
    assert_("metar: every temperature lies within 180-340 K", ta and all(180 <= x <= 340 for x in ta),
            f"{min(ta)} to {max(ta)} K")
    rh = [r[column(metar, "relative_humidity")] for r in rows if r[column(metar, "relative_humidity")] is not None]
    assert_("metar: every relative humidity lies within 0-102 %", rh and all(0 <= x <= 102 for x in rh))
    assert_("metar: each station's reports in time order, once per time", times_ordered(metar))
    known = [r for r in rows if r[0] in ("RJTT", "KJFK")]
    if not known:
        assert_("metar: Tokyo Haneda (RJTT) or New York JFK (KJFK) reports", False)
    else:
        r = known[-1]
        lon, lat = (139.78, 35.55) if r[0] == "RJTT" else (-73.76, 40.64)
        check(f"metar: {r[0]} is where it is (distance in degrees)", math.hypot(r[2] - lon, r[3] - lat), 0.05, "degree")
        t, p = r[column(metar, "air_temperature")], r[column(metar, "air_pressure_at_sea_level")]
        assert_(f"metar: {r[0]} reports a plausible temperature ({t} K) and pressure ({p} Pa)",
                t is not None and 250 <= t <= 320 and p is not None and 95000 <= p <= 105000)

ndbc = load("observations-ndbc.json")
if ndbc is None:
    skip("ndbc: observations-ndbc.json not built (fetch.py ndbc; observation_set.py ndbc)")
else:
    rows = ndbc["rows"]
    assert_("ndbc: more than 300 stations report", len({r[0] for r in rows}) > 300, f"{len({r[0] for r in rows})} stations")
    sst = [r[column(ndbc, "sea_surface_temperature")] for r in rows if r[column(ndbc, "sea_surface_temperature")] is not None]
    assert_("ndbc: sea-surface temperatures, 100 or more, lie within 269-312 K",
            len(sst) >= 100 and all(269 <= x <= 312 for x in sst), f"{len(sst)}: {min(sst, default=None)} to {max(sst, default=None)}")
    hs = [r[column(ndbc, "sea_surface_wave_significant_height")] for r in rows
          if r[column(ndbc, "sea_surface_wave_significant_height")] is not None]
    assert_("ndbc: significant wave heights lie within 0-25 m", all(0 <= x <= 25 for x in hs), f"{len(hs)} heights")
    assert_("ndbc: each station's reports in time order, once per time", times_ordered(ndbc))

for profile in ("active", "last3years"):
    product = load(f"tracks-ibtracs-{profile}.json")
    if product is None:
        skip(f"ibtracs: tracks-ibtracs-{profile}.json not built (fetch.py ibtracs; ibtracs_tracks.py)")
        continue
    tracks = product["tracks"]
    ordered = all(all(b[0] > a[0] for a, b in zip(t["points"], t["points"][1:])) for t in tracks)
    assert_(f"ibtracs {profile}: every track has positions, in time order", tracks and ordered
            and all(t["points"] for t in tracks), f"{len(tracks)} tracks")
    inside = all(-180 <= p[1] < 180 and -90 <= p[2] <= 90 for t in tracks for p in t["points"])
    assert_(f"ibtracs {profile}: positions within [-180, 180) and [-90, 90]", inside)
    plausible = all((p[4] is None or 850 <= p[4] <= 1050) and (p[5] is None or 0 < p[5] < 100)
                    for t in tracks for p in t["points"])
    assert_(f"ibtracs {profile}: pressures within 850-1050 hPa and winds within 0-100 m s-1", plausible)
    if profile == "last3years":
        # Hurricane Milton, 2024: 155 kt and 895 hPa at its peak (NHC Tropical Cyclone Report, AL142024).
        milton = next((t for t in tracks if t["name"] == "MILTON" and t["season"] == 2024), None)
        if milton is None:
            skip("ibtracs: Milton (2024) has left the last three years")
        else:
            check("ibtracs: Milton's peak wind against NHC's 155 kt", abs(milton["max_wind_ms"] / obs.KNOT - 155.0), 0.5, "kt")
            check("ibtracs: Milton's lowest pressure against NHC's 895 hPa", abs(milton["min_pressure_hpa"] - 895.0), 0.5, "hPa")
