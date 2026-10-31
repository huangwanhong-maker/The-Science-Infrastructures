"""Checks of the street records (data/formats/opening_hours.py, data/products/street_records.py), offline:
the reading of OpenStreetMap's opening_hours over the cases of its grammar the subset covers and those it
leaves out; Melbourne's civil time (the rule against the IANA database where Python has one); the uses of
premises from ANZSIC classes and OpenStreetMap tags; the street axis and projection onto it; the exact
east-north-up projection (closed forms on a parallel and a meridian), ring assembly, winding and clipping
of the city model; lengths, heights and widths read from tags. Where the products have been built, their
shape, and a few numbers seen by eye. Nothing is fetched.
"""

from __future__ import annotations

import datetime as dt
import json
import math

import numpy as np

from checks import PRODUCTS, assert_, check, skip

import opening_hours as oh  # noqa: E402
import street_records as sr  # noqa: E402

# ---- opening_hours: what the subset reads ------------------------------------------------------------------
P = oh.parse
WEEK = lambda a, b, days=range(7): [[d, a, b] for d in days]  # noqa: E731
cases = [
    ("Mo-Fr 10:00-19:00; Sa 10:00-17:00; Su,PH off", WEEK(600, 1140, range(5)) + [[5, 600, 1020]], [], "ok"),
    ("24/7", WEEK(0, 1440), None, "ok"),
    ("Mo-Su 00:00-24:00", WEEK(0, 1440), None, "ok"),
    ("08:00-18:00", WEEK(480, 1080), None, "ok"),
    ("Mo,We,Fr", [[0, 0, 1440], [2, 0, 1440], [4, 0, 1440]], None, "ok"),
    ("Fr-Mo 10:00-12:00", [[0, 600, 720], [4, 600, 720], [5, 600, 720], [6, 600, 720]], None, "ok"),
    ("Mo-Fr 08:00-12:00,13:00-17:30", [x for d in range(5) for x in ([d, 480, 720], [d, 780, 1050])], None, "ok"),
    ("Mo-Su 22:00-02:00; Tu off", [[d, 1320, 1560] for d in (0, 2, 3, 4, 5, 6)], None, "ok"),
    ("Mo-Fr 18:00-26:00", WEEK(1080, 1560, range(5)), None, "ok"),
    ("Mo-Fr 08:00-12:00, We 14:00-18:00", [[0, 480, 720], [1, 480, 720], [2, 480, 720], [2, 840, 1080], [3, 480, 720], [4, 480, 720]], None, "ok"),
    ("Mo-Fr 08:00-18:00, We off", [[d, 480, 1080] for d in (0, 1, 3, 4)], None, "ok"),
    ("Mo-Fr 10:00-20:00; PH 11:00-15:00", WEEK(600, 1200, range(5)), [[660, 900]], "ok"),
    ("off", [], None, "ok"),
    ("Mo-Fr: 9:00 – 17:00", WEEK(540, 1020, range(5)), None, "ok"),
    ("11:00～20:00", WEEK(660, 1200), None, "ok"),
    ("mo-fr 10:00-18:00", WEEK(600, 1080, range(5)), None, "ok"),
    ('Mo-Fr 10:00-18:00 "ring first"', WEEK(600, 1080, range(5)), None, "ok"),
    ("Mo-Fr 10:00-18:00; Dec 25 off", WEEK(600, 1080, range(5)), None, "partial"),
    ("Apr-Oct Mo-Fr 10:00-18:00; Nov-Mar Mo-Fr 10:00-16:00", [], None, "partial"),
    ("week 01-10 Mo 10:00-12:00; Tu 09:00-10:00", [[1, 540, 600]], None, "partial"),
    ("2024 Mo 10:00-12:00", [], None, "partial"),
    ("Th[1] 10:00-12:00; Mo 09:00-10:00", [[0, 540, 600]], None, "partial"),
    ("sunrise-sunset", [], None, "partial"),
    ("Mo-Sa 10:00+", [], None, "partial"),
    ("Mo-Su 12:00-18:00+", WEEK(720, 1080), None, "partial"),
    ("We,Th 12:00-15:00, Fr-Su 11:00-15:00, We-Sa 17:00-20:00+", [[2, 720, 900], [2, 1020, 1200], [3, 720, 900], [3, 1020, 1200],
     [4, 660, 900], [4, 1020, 1200], [5, 660, 900], [5, 1020, 1200], [6, 660, 900]], None, "partial"),
    ("Mo-Fr 10:00-16:00/01:00", WEEK(600, 960, range(5)), None, "partial"),
    ("Mo-Fr 09:00-17:00 || \"by appointment\"", WEEK(540, 1020, range(5)), None, "partial"),
    ("Mo-Fr,SH 08:00-09:00", WEEK(480, 540, range(5)), None, "partial"),
    ('"by appointment"', [], None, "partial"),
    ("Mo-Fr 10:00-18:00; Sat 10:00-14:00", WEEK(600, 1080, range(5)), None, "partial"),
    ("Mo-Fr 09:00-17:00 unknown", [], None, "partial"),
    ("Mon-Fri 9am-5pm", [], None, "unparsed"),
    ("10am-6pm", [], None, "unparsed"),
    (None, [], None, "none"),
    ("", [], None, "none"),
]
for text, hours, ph, status in cases:
    r = P(text)
    assert_(f"opening_hours: {text!r} reads as {status}", r["hours"] == hours and r["ph"] == ph and r["status"] == status,
            json.dumps(r, ensure_ascii=False))
assert_("opening_hours: what is left out is said, with its reason", "a year, month, date or week selector" in P("Mo 10:00-12:00; Dec 25 off")["dropped"][0])
h = P("Fr 22:00-02:00")["hours"]
assert_("opening_hours: open past midnight counts on the next day", oh.is_open(h, 5, 60) and not oh.is_open(h, 5, 120) and oh.is_open(h, 4, 1380))
assert_("opening_hours: a week's minutes", oh.weekly_minutes(P("Mo-Fr 09:00-17:00")["hours"]) == 5 * 480 and oh.weekly_minutes(P("24/7")["hours"]) == 7 * 1440)

# ---- Melbourne's civil time ------------------------------------------------------------------------------------
UTC = dt.timezone.utc
end, start = sr.transitions(2025)
assert_("civil time: 2025's changes at 16:00 UTC on the Saturdays before 6 April and 5 October",
        end == dt.datetime(2025, 4, 5, 16, tzinfo=UTC) and start == dt.datetime(2025, 10, 4, 16, tzinfo=UTC))
assert_("civil time: AEDT in January, AEST in July", sr.offset_hours(dt.datetime(2024, 1, 15, tzinfo=UTC)) == 11
        and sr.offset_hours(dt.datetime(2024, 7, 15, tzinfo=UTC)) == 10)
u, kind = sr.local_to_utc(dt.date(2025, 4, 6), 2)
assert_("civil time: the repeated 02:00 of 6 April 2025 is its first, 15:00 UTC", kind == "repeated" and sr.iso(u) == "2025-04-05T15:00:00Z")
assert_("civil time: 02:00 on 5 October 2025 does not exist", sr.local_to_utc(dt.date(2025, 10, 5), 2) == (None, "skipped"))
assert_("civil time: 03:00 on 5 October 2025 is 16:00 UTC", sr.iso(sr.local_to_utc(dt.date(2025, 10, 5), 3)[0]) == "2025-10-04T16:00:00Z")
try:
    from zoneinfo import ZoneInfo
    mel = ZoneInfo("Australia/Melbourne")
    bad = 0
    t = dt.datetime(2009, 1, 1, tzinfo=UTC)
    while t.year < 2027:
        bad += sr.offset_hours(t) != t.astimezone(mel).utcoffset().total_seconds() / 3600
        t += dt.timedelta(hours=1)
    check("civil time: the rule against the IANA database, every hour 2009-2026", bad, 0, "hours")
except Exception:  # noqa: BLE001 (no time-zone database in this Python)
    skip("civil time: no IANA time-zone database in this Python (pip install tzdata)")
assert_("civil time: the offsets over 2025 change twice", [o[1] for o in sr.offset_changes(sr.hour_of(dt.datetime(2025, 1, 1, tzinfo=UTC)),
                                                                                          sr.hour_of(dt.datetime(2025, 12, 31, tzinfo=UTC)))] == [11, 10, 11])

# ---- uses -----------------------------------------------------------------------------------------------------
assert_("uses: ANZSIC classes by their longest prefix",
        [sr.anzsic_use(c) for c in ("4511", "4253", "4255", "4110", "6931", "6221", "0000", "8102", "9533", "6970", "5513", "2591")]
        == ["food-and-drink", "jewellery", "clothing-and-fashion", "other-retail", "offices-and-professional", "finance", "vacant",
            "education", "transport", "health", "entertainment-and-culture", "wholesale-and-industry"])
assert_("uses: the bridal heuristic on names, as words, never for jewellers",
        sr.anzsic_use("4251", "Brides of Melbourne") == "bridal-and-wedding" and sr.anzsic_use("7211", "Cleveland McBride") == "offices-and-professional"
        and sr.anzsic_use("4253", "Bridal Bizzarr") == "jewellery" and sr.anzsic_use("4400", "Stamford Plaza Wedding & Conference Centre") == "bridal-and-wedding")
assert_("uses: OpenStreetMap tags", [sr.osm_use(t) for t in (
    {"shop": "jewelry"}, {"shop": "watches"}, {"shop": "wedding"}, {"shop": "clothes", "clothes": "wedding;women"},
    {"shop": "clothes", "name": "銀座ブライダル"}, {"amenity": "restaurant"}, {"amenity": "bank"}, {"office": "lawyer"},
    {"tourism": "hotel"}, {"healthcare": "dentist"}, {"shop": "vacant"}, {"craft": "brewery"}, {"shop": "books"})] ==
    ["jewellery", "jewellery", "bridal-and-wedding", "bridal-and-wedding", "bridal-and-wedding", "food-and-drink", "finance",
     "offices-and-professional", "accommodation", "health", "vacant", "wholesale-and-industry", "other-retail"])
assert_("uses: furniture is not a holding", not sr.is_premises({"amenity": "bench"}) and not sr.is_premises({"leisure": "park"})
        and sr.is_premises({"amenity": "cafe"}) and sr.is_premises({"amenity": "bench", "shop": "kiosk"}))

# ---- a street's axis and the projection onto it -------------------------------------------------------------------
# Two one-way carriageways 10 m apart along x, a holding 25 m north of the centreline.
segs = [((0.0, -5.0), (200.0, -5.0)), ((200.0, 5.0), (0.0, 5.0))]
line = sr.centreline(segs, (0.0, 0.0), (200.0, 0.0))
check("streets: the centreline lies midway between the carriageways", float(np.max(np.abs(line[:, 1]))), 1e-6, "m")
check("streets: the axis runs from one end to the other", abs(float(np.hypot(*np.diff(line, axis=0).T).sum()) - 200.0), 1e-6, "m")
s, off, beside = sr.project(line, (120.0, 25.0))
assert_("streets: s along the axis, offset positive on the left", abs(s - 120) < 1e-9 and abs(off - 25) < 1e-9 and beside)
assert_("streets: a point beyond an end is not beside the axis", not sr.project(line, (-30.0, 2.0))[2])

# ---- the exact east-north-up projection --------------------------------------------------------------------------------
for name, (lat0, lon0) in (("Melbourne", (-37.8095, 144.9525)), ("35.672° N", (35.672, 139.7665))):
    e, n, u = sr.enu(lat0, lon0, lat0, lon0)
    assert_(f"ENU: the origin of {name} is (0, 0, 0)", abs(float(e)) < 1e-9 and abs(float(n)) < 1e-9 and abs(float(u)) < 1e-9)
    e, n, u = sr.enu(lat0 + 0.001, lon0, lat0, lon0)
    M = sr.meridian_radius(lat0 + 0.0005)
    check(f"ENU: 0.001° north of {name} is M(φ) Δφ = {M * math.radians(0.001):.3f} m north", abs(float(n) - M * math.radians(0.001)), 1e-4, "m")
    check(f"ENU: ... and no east", abs(float(e)), 1e-9, "m")
    # On the origin's parallel: e = N cos φ0 sin Δλ, u = N cos² φ0 (cos Δλ − 1) exactly.
    N = sr.A_WGS84 / math.sqrt(1 - sr.E2 * math.sin(math.radians(lat0)) ** 2)
    d = math.radians(0.02)
    e, n, u = sr.enu(lat0, lon0 + 0.02, lat0, lon0)
    check(f"ENU: 0.02° east of {name} on its parallel (closed form)", abs(float(e) - N * math.cos(math.radians(lat0)) * math.sin(d))
          + abs(float(u) - N * math.cos(math.radians(lat0)) ** 2 * (math.cos(d) - 1)), 1e-6, "m")
    lat = lat0 + np.array([-0.015, 0.015]); lon = lon0 + np.array([-0.02, 0.02])
    e, n, _ = sr.enu(lat, lon, lat0, lon0)
    x, y = sr.Plane(lat0, lon0).xy(lat, lon)
    check(f"ENU: the equirectangular plane within 2 km of {name}'s centre", float(np.max(np.hypot(x - e, y - n))), 0.5, "m")

# ---- rings ------------------------------------------------------------------------------------------------------------
sq = [(0, 0), (0, 1), (1, 1), (1, 0)]
r = sr.assemble_rings([[(0, 0), (0, 1)], [(1, 0), (1, 1)], [(0, 1), (1, 1)], [(1, 0), (0, 0)]])
assert_("rings: ways joined end to end, reversed where they must be, into one closed ring", r is not None and len(r) == 1
        and r[0][0] == r[0][-1] and len(r[0]) == 5 and set(r[0]) == set(sq))
assert_("rings: ways that do not close are refused", sr.assemble_rings([[(0, 0), (0, 1)], [(0, 1), (1, 1)]]) is None)
two = sr.assemble_rings([[(0, 0), (0, 1), (1, 1)], [(1, 1), (1, 0), (0, 0)], [(5, 5), (5, 6), (6, 6), (5, 5)]])
assert_("rings: two rings from three ways", two is not None and len(two) == 2)
ccw = sr.oriented([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)], True)
assert_("rings: footprints counter-clockwise and not closed", sr.signed_area(ccw) > 0 and len(ccw) == 4 and ccw[0] != ccw[-1])
assert_("rings: holes clockwise", sr.signed_area(sr.oriented(sq, False)) < 0)
mp = {"type": "relation", "members": [
    {"type": "way", "role": "outer", "geometry": [{"lat": 0, "lon": 0}, {"lat": 0, "lon": 4}, {"lat": 4, "lon": 4}]},
    {"type": "way", "role": "outer", "geometry": [{"lat": 4, "lon": 4}, {"lat": 4, "lon": 0}, {"lat": 0, "lon": 0}]},
    {"type": "way", "role": "inner", "geometry": [{"lat": 1, "lon": 1}, {"lat": 1, "lon": 2}, {"lat": 2, "lon": 2}, {"lat": 2, "lon": 1}, {"lat": 1, "lon": 1}]}]}
pg = sr.polygons(mp)
assert_("rings: a multipolygon's hole goes to the outer that holds it", pg is not None and len(pg) == 1 and len(pg[0][1]) == 1)
clipped = sr.clip_to_box([(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)], [1, -1, 3, 3])
check("rings: a square half outside its box keeps half its area", abs(abs(sr.signed_area(clipped)) - 2.0), 1e-12)
assert_("lines: a road leaving and re-entering its box is cut into pieces",
        sr.clip_line_to_box([(0.5, -1), (0.5, 0.5), (1.5, 0.5), (0.5, 0.8)], [0, 0, 1, 1]) == [[(0.5, 0.0), (0.5, 0.5), (1.0, 0.5)], [(1.0, 0.65), (0.5, 0.8)]]
        and sr.clip_line_to_box([(2, 2), (3, 3)], [0, 0, 1, 1]) == [])

# ---- tags read as lengths ----------------------------------------------------------------------------------------------
assert_("lengths: units stripped, feet converted, the first of a list",
        [sr.parse_length(v) for v in ("12", "12 m", "12.5m", "40'", "10'6\"", "12;15", "~20", "abc", None)] ==
        [12.0, 12.0, 12.5, 40 * 0.3048, 10 * 0.3048 + 6 * 0.0254, 12.0, 20.0, None, None])
assert_("heights: the height tag, then levels (and roof levels), then the kind's default",
        sr.building_heights({"height": "45 m"}, "office")[::3] == (45.0, "height")
        and sr.building_heights({"building:levels": "10", "roof:levels": "1"}, "yes")[0] == 10 * 3.2 + 1.6
        and sr.building_heights({}, "house")[::3] == (6.4, "default") and sr.building_heights({}, "yes")[0] == sr.DEFAULT_HEIGHT_OTHER)
assert_("heights: the base from min_height or min_level, below the top",
        sr.building_heights({"height": "30", "min_height": "12"}, "yes")[1] == 12.0
        and sr.building_heights({"building:levels": "4", "building:min_level": "2"}, "yes")[1] == 6.4
        and sr.building_heights({"height": "3", "min_height": "5"}, "yes")[0] > 5)
assert_("widths: the width tag, then lanes and sidewalks, then the class's default",
        sr.road_width({"highway": "primary", "width": "18"}) == (18.0, "width")
        and sr.road_width({"highway": "secondary", "lanes": "2", "sidewalk": "both"}) == (2 * 3.2 + 3.0, "lanes")
        and sr.road_width({"highway": "primary", "oneway": "yes"}) == (6.5, "default") and sr.road_width({"highway": "footway"}) == (2.0, "default"))
assert_("spaces: kinds from tags", [sr.space_kind(t) for t in ({"leisure": "park"}, {"natural": "water"}, {"place": "square"},
                                                                {"landuse": "grass"}, {"landuse": "railway"}, {"building": "yes"})]
        == ["park", "water", "plaza", "green", "rail", None])

# ---- the products, where built ---------------------------------------------------------------------------------------------
path = PRODUCTS / "street-sensors.json"
if path.exists():
    sensors = json.loads(path.read_text(encoding="utf-8"))
    keys = {"id", "code", "name", "lat", "lon", "installed", "status", "location_type", "directions", "series", "first_utc", "last_utc", "hours", "missing"}
    assert_("street-sensors: every sensor has the contract's fields", all(keys <= set(s) for s in sensors["sensors"]) and sensors["timezone"] == "Australia/Melbourne")
    bad = []
    for s in sensors["sensors"][:200]:
        if not s["series"]:
            continue
        c = json.loads((PRODUCTS / s["series"]).read_text(encoding="utf-8"))
        nulls = sum(v is None for v in c["counts"])
        if len(c["counts"]) != s["hours"] or nulls != s["missing"] or c["t0_utc"] != s["first_utc"] or c["step_s"] != 3600 \
                or c["counts"][0] is None or c["counts"][-1] is None or any(v is not None and v < 0 for v in c["counts"]):
            bad.append(s["id"])
    assert_("street-counts: each series spans its sensor's first to last record, its nulls counted", not bad, str(bad))
    one = json.loads((PRODUCTS / "street-counts" / "1.json").read_text(encoding="utf-8"))
    t0 = dt.datetime.fromisoformat(one["t0_utc"].replace("Z", "+00:00"))
    gap = [one["counts"][int((dt.datetime(2023, 6, 1, tzinfo=UTC) - t0).total_seconds() // 3600) + k] for k in range(24)]
    assert_("street-counts: the City's gap (November 2022 to October 2024) is null", all(v is None for v in gap))
    by_hour = [[] for _ in range(24)]
    for i, v in enumerate(one["counts"]):
        t = t0 + dt.timedelta(hours=i)
        if v is not None and t.year == 2019:
            local = t + dt.timedelta(hours=sr.offset_hours(t))
            if local.weekday() < 5:
                by_hour[local.hour].append(v)
    mean = [sum(x) / len(x) if x else 0 for x in by_hour]
    assert_("street-counts: Bourke Street Mall (North) is busiest at lunch on weekdays (2019)", 12 <= mean.index(max(mean)) <= 14, str([round(m) for m in mean]))
    def day(date):
        k0 = int((dt.datetime(date.year, date.month, date.day, tzinfo=UTC) - dt.timedelta(hours=11) - t0).total_seconds() // 3600)
        return sum(v for v in one["counts"][k0:k0 + 24] if v)
    assert_("street-counts: Christmas Day 2019 is quieter in the Mall than the Wednesday before", day(dt.date(2019, 12, 25)) < 0.9 * day(dt.date(2019, 12, 18)))
    # The directions: d1 + d2 is the count at every hour where both are recorded, on every series that has them;
    # their names are the locations dataset's, absent for a sensor it no longer lists.
    bad, n, unnamed = [], 0, []
    for s in sensors["sensors"]:
        if not s["series"]:
            continue
        c = json.loads((PRODUCTS / s["series"]).read_text(encoding="utf-8"))
        d = c.get("directions")
        if not d:
            continue
        n += 1
        k0 = int((dt.datetime.fromisoformat(d["t0_utc"].replace("Z", "+00:00")) - dt.datetime.fromisoformat(c["t0_utc"].replace("Z", "+00:00"))).total_seconds() // 3600)
        if len(d["d1"]) != len(d["d2"]) or any(a is not None and b is not None and c["counts"][k0 + i] != a + b for i, (a, b) in enumerate(zip(d["d1"], d["d2"]))):
            bad.append(s["id"])
        if len(s["directions"]) != 2:
            unnamed.append(s["id"])
    assert_(f"street-counts: each direction's counts sum to the count, at {n} sensors ({len(unnamed)} unlisted, their directions unnamed)",
            n >= 90 and not bad and all(next(x for x in sensors["sensors"] if x["id"] == u)["status"] == "unlisted" for u in unnamed), str(bad))
else:
    skip("street-sensors: not built (python data/products/street_records.py --only sensors)")
path = PRODUCTS / "street-weather-melbourne.json"
if path.exists():
    w = json.loads(path.read_text(encoding="utf-8"))
    n = len(w["temperature_c"])
    assert_("street-weather: three hourly series of one length", len(w["precipitation_mm"]) == n == len(w["wind_m_s"]) and w["step_s"] == 3600)
    assert_("street-weather: plausible values", all(v is None or -10 < v < 50 for v in w["temperature_c"]) and all(v is None or 0 <= v < 200 for v in w["precipitation_mm"])
            and all(v is None or 0 <= v < 60 for v in w["wind_m_s"]))
else:
    skip("street-weather-melbourne: not built")
path = PRODUCTS / "street-holidays-victoria.json"
if path.exists():
    hol = {h["date"]: h for h in json.loads(path.read_text(encoding="utf-8"))["holidays"]}
    assert_("street-holidays: Melbourne Cup Day 2019 (5 November), Victorian, not national", "2019-11-05" in hol and not hol["2019-11-05"]["national"])
    assert_("street-holidays: Christmas Day 2025 and Anzac Day 2024", "2025-12-25" in hol and "2024-04-25" in hol)
else:
    skip("street-holidays-victoria: not built")
path = PRODUCTS / "street-terrain-melbourne.json"
if path.exists():
    ter = json.loads(path.read_text(encoding="utf-8"))
    uses = {u["id"] for u in ter["uses"]}
    rows = [e for s in ter["sensors"].values() for e in s["establishments"]]
    assert_("street-terrain: every establishment within 200 m, its use one of the classes", all(math.hypot(e["x_m"], e["y_m"]) <= 200.05 and e["use"] in uses for e in rows))
    assert_("street-terrain: jewellers about the Mall", any(e["use"] == "jewellery" for e in ter["sensors"]["1"]["establishments"]))
else:
    skip("street-terrain-melbourne: not built")
path = PRODUCTS / "street-holdings.json"
if path.exists():
    hs = json.loads(path.read_text(encoding="utf-8"))
    ok = all(0 <= d <= 6 and 0 <= a < b <= a + 1440 and a < 1440 for s in hs["streets"] for h in s["holdings"] for d, a, b in h["hours"])
    assert_("street-holdings: hours are [day, open, close] within a day and its night", ok)
    assert_("street-holdings: holdings within 40 m and on the axis", all(abs(h["offset_m"]) <= 40 and 0 <= h["s_m"] <= s["length_m"] + 0.1
                                                                         for s in hs["streets"] for h in s["holdings"]))
    g = next((s for s in hs["streets"] if s["id"] == "melbourne-swanston"), None)
    if g:
        check("street-holdings: Swanston Street from Flinders Street to La Trobe Street is about 930 m", abs(g["length_m"] - 930), 40, "m")
        assert_("street-holdings: shops with declared hours on Swanston Street", sum(bool(h["hours"]) for h in g["holdings"]) >= 10)
else:
    skip("street-holdings: not built")
path = PRODUCTS / "city-model.json"
if path.exists():
    cm = json.loads(path.read_text(encoding="utf-8"))
    for a in cm["areas"]:
        fp = [b for b in a["buildings"]]
        assert_(f"city-model {a['id']}: footprints counter-clockwise, holes clockwise, not closed",
                all(sr.signed_area(b["footprint"]) > 0 and b["footprint"][0] != b["footprint"][-1] and all(sr.signed_area(h) < 0 for h in b["holes"]) for b in fp))
        assert_(f"city-model {a['id']}: heights above their bases", all(b["height_m"] > b["min_height_m"] >= 0 for b in fp))
        assert_(f"city-model {a['id']}: roads have widths and lines", all(r["width_m"] > 0 and len(r["line"]) >= 2 for r in a["roads"]))
        lat0, lon0, lat1, lon1 = a["bbox"]
        e, n, _ = sr.enu(np.array([lat0, lat0, lat1, lat1]), np.array([lon0, lon1, lon0, lon1]), a["origin"]["lat"], a["origin"]["lon"])
        inside = sum(1 for s in a["spaces"] for x, y in s["ring"] if min(e) - 2 <= x <= max(e) + 2 and min(n) - 2 <= y <= max(n) + 2)
        assert_(f"city-model {a['id']}: open spaces clipped to the area", inside == sum(len(s["ring"]) for s in a["spaces"]))
else:
    skip("city-model: not built (python data/products/street_records.py --only city)")
