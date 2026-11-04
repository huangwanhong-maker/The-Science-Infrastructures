"""Checks of the record of a district's events (data/products/city_events.py), offline: ANZSIC 2006's divisions against
classes the ABS publishes; points joined to footprints with holes, parts and neighbours; time labels at their
resolutions and the window; streets found in free text; the frame's inverse. Where the product has been built, its
schema, its references between tables, and its agreement with the city model and the workbench's district. Nothing is
fetched.
"""

from __future__ import annotations

import datetime as dt
import json
import math
import re

from checks import PRODUCTS, ROOT, assert_, check, skip

import city_events as ce  # noqa: E402

# ---- ANZSIC 2006: classes and their divisions, as the ABS lists them -----------------------------------------------
# Australian Bureau of Statistics, ANZSIC 2006 (Revision 2.0), cat. no. 1292.0, "Division, subdivision, group and class
# codes and titles" (retrieved 2026-11-04): one class of every division.
PUBLISHED = [("0111", "Nursery Production (Under Cover)", "A"), ("0600", "Coal Mining", "B"), ("1111", "Meat Processing", "C"),
             ("2611", "Fossil Fuel Electricity Generation", "D"), ("3011", "House Construction", "E"), ("3311", "Wool Wholesaling", "F"),
             ("4110", "Supermarket and Grocery Stores", "G"), ("4511", "Cafes and Restaurants", "H"), ("4610", "Road Freight Transport", "I"),
             ("5411", "Newspaper Publishing", "J"), ("6221", "Banking", "K"), ("6611", "Passenger Car Rental and Hiring", "L"),
             ("6931", "Legal Services", "M"), ("7000", "Computer System Design and Related Services", "M"),
             ("7211", "Employment Placement and Recruitment Services", "N"), ("7520", "State Government Administration", "O"),
             ("8021", "Primary Education", "P"), ("8401", "Hospitals (Except Psychiatric Hospitals)", "Q"), ("8910", "Museum Operation", "R"),
             ("9411", "Automotive Electrical Services", "S"), ("9601", "Private Households Employing Staff", "S")]
wrong = [(c, d, ce.anzsic_division(int(c[:2]))) for c, _, d in PUBLISHED if ce.anzsic_division(int(c[:2])) != d]
assert_("anzsic: the division of a class of each of the 19 divisions, as the ABS lists it", not wrong and len({d for _, _, d in PUBLISHED}) == 19, str(wrong))
assert_("anzsic: the numbers between divisions are no subdivision (61, 65, 68, 71, 74, 78, 79, 83, 88, 93; and 00, 97)",
        all(ce.anzsic_division(n) is None for n in (0, 61, 65, 68, 71, 74, 78, 79, 83, 88, 93, 97)))
assert_("anzsic: the divisions cover 01-96 in order without overlap",
        all(a[2] < b[1] for a, b in zip(ce.DIVISIONS, ce.DIVISIONS[1:])) and ce.DIVISIONS[0][1] == 1 and ce.DIVISIONS[-1][2] == 96)
lv = ce.anzsic_levels("4511")
assert_("anzsic: a class's levels (4511: H, 45, 451, 4511)", lv == {"division": "H", "subdivision": "45", "group": "451", "class": "4511"}, str(lv))
assert_("anzsic: CLUE's vacant space (0000) has no division", ce.anzsic_levels("0000")["division"] is None)

# ---- points joined to footprints ---------------------------------------------------------------------------------------
square = [[0, 0], [10, 0], [10, 10], [0, 10]]
court = [[4, 4], [4, 6], [6, 6], [6, 4]]                      # a hole, clockwise
assert_("join: a point in a footprint is inside it; in its hole, outside; beyond it, outside",
        ce.point_in_footprint((2, 2), square, [court]) and not ce.point_in_footprint((5, 5), square, [court])
        and not ce.point_in_footprint((12, 5), square, [court]))
model = [{"id": "w/1", "kind": "yes", "footprint": square, "holes": [court]},
         {"id": "w/2", "kind": "yes", "footprint": [[20, 0], [40, 0], [40, 20], [20, 20]], "holes": []},
         {"id": "w/3", "kind": "yes", "footprint": [[25, 5], [30, 5], [30, 10], [25, 10]], "holes": []},
         {"id": "p/4", "kind": "part:yes", "footprint": [[21, 1], [24, 1], [24, 4], [21, 4]], "holes": []},
         {"id": "p/5", "kind": "part:yes", "footprint": [[100, 0], [110, 0], [110, 10], [100, 10]], "holes": []}]
F = ce.Footprints(model, [-50, -50, 150, 50])
got = {p: F.join(p) for p in [(2, 2), (5, 5), (13, 5), (60, 5), (27, 7), (22, 2), (105, 5)]}
assert_("join: inside the footprint that holds the point", got[(2, 2)][:2] == ("w/1", "inside"))
assert_("join: a point in a court is nearest its building, 1 m from the court's edge", got[(5, 5)][:2] == ("w/1", "nearest") and abs(got[(5, 5)][2] - 1.0) < 1e-9, str(got[(5, 5)]))
assert_("join: 3 m outside, the nearest footprint", got[(13, 5)][:2] == ("w/1", "nearest") and abs(got[(13, 5)][2] - 3.0) < 1e-9, str(got[(13, 5)]))
assert_("join: 20 m from every footprint, none", got[(60, 5)] == (None, "none", None), str(got[(60, 5)]))
assert_("join: of two outlines that hold a point, the smaller", got[(27, 7)][:2] == ("w/3", "inside"), str(got[(27, 7)]))
assert_("join: an outline before a part that holds the point too", got[(22, 2)][:2] == ("w/2", "inside"), str(got[(22, 2)]))
assert_("join: a part when no outline holds the point", got[(105, 5)][:2] == ("p/5", "inside"), str(got[(105, 5)]))
assert_("join: the labels (inside, nearest:<m to 0.1>, none)", ce.join_label("nearest", 3.04) == "nearest:3.0" and ce.join_label("inside", None) == "inside")

# ---- time labels -------------------------------------------------------------------------------------------------------
assert_("labels: a year is 1 January to the next 1 January", ce.year_label(2002) == ("2002-01-01", "2003-01-01") and ce.year_label(1999) == ("1999-01-01", "2000-01-01"))
assert_("labels: a day is that day to the next, across a leap day and a year's end",
        ce.day_label(dt.date(2016, 2, 28), dt.date(2016, 2, 29)) == ("2016-02-28", "2016-03-01")
        and ce.day_label(dt.date(2017, 12, 31)) == ("2017-12-31", "2018-01-01"))
assert_("labels: the window keeps 1952 to 2024, from fifty years before 2002 to the last census",
        [ce.in_window(ce.year_label(y)) for y in (1951, 1952, 2001, 2002, 2024, 2025)] == ["before", "in", "in", "in", "in", "after"])
assert_("labels: the Heritage Database's dates read", ce.heritage_date("October 9, 1974") == dt.date(1974, 10, 9) and ce.heritage_date("") is None)

# ---- streets in free text -----------------------------------------------------------------------------------------------
S = ce.StreetNames(["Bourke Street", "Bourke Street Mall", "Collins Street", "Little Collins Street", "Flinders Street", "Flinders Lane"])
found, rest = S.find("Bourke St Mall - next to Tourism booth, Little Collins St, Flinders St Station, flinders lane")
assert_("streets: abbreviations, the longest name first, a station not its street, any case",
        found == ["Bourke Street Mall", "Little Collins Street", "Flinders Lane"], str(found))
assert_("streets: Collins Street is not found inside Little Collins Street", "Collins Street" not in S.find("Little Collins Street")[0])
assert_("streets: other street-like names counted with their word written out",
        ce.other_street_names("St Kilda Rd, Spring St, Queensbridge Square") == ["St Kilda Road", "Spring Street"], str(ce.other_street_names("St Kilda Rd, Spring St, Queensbridge Square")))
assert_("permits: categories' variants made one", ce.permit_kind("Pubilc Event - Non Ticketed") == "Public Event - Non-ticketed"
        and ce.permit_kind("Public Event - Run Walk") == "Public Event - Run/Walk" and ce.permit_kind("Public Event - Low Impact Activity,") == "Public Event - Low Impact Activity")

# ---- the frame -----------------------------------------------------------------------------------------------------------
fr = ce.Frame({"lat": -37.8085, "lon": 144.9525})
x, y = fr.xy(-37.8173, 144.9697)
lat, lon = fr.latlon(float(x), float(y))
check("frame: the inverse returns the place (a district's corner, 1.8 km away)", max(abs(lat + 37.8173), abs(lon - 144.9697)) * 111e3, 1e-3, "m")

# ---- the product ---------------------------------------------------------------------------------------------------------
path = PRODUCTS / "city-events-melbourne.json"
if path.exists():
    p = json.loads(path.read_text(encoding="utf-8"))
    keys = {"product", "id", "built", "frame", "origin", "box", "window", "sources", "properties", "censuses", "establishments", "anzsic", "events", "report"}
    assert_("city-events: its keys", keys <= set(p) and p["product"] == "city-events" and p["id"] == "melbourne", str(sorted(keys - set(p))))
    cols = {"properties": ["id", "base_property_id", "property_ids", "block", "at", "building", "join"],
            "censuses": ["year", "property", "property_id", "construction_year", "refurbished_year", "floors_above", "predominant_use"],
            "establishments": ["year", "property", "property_id", "name", "anzsic4"],
            "events": ["id", "type", "kind", "from", "to", "match", "source", "provenance", "detail"]}
    assert_("city-events: the tables' columns", all(p[t]["columns"][:len(c)] == c for t, c in cols.items())
            and p["events"]["columns"] == cols["events"] and p["properties"]["columns"] == cols["properties"])
    assert_("city-events: every row as long as its columns", all(len(r) == len(p[t]["columns"]) for t in cols for r in p[t]["rows"]))
    props = {r[0]: r for r in p["properties"]["rows"]}
    assert_("city-events: properties named base:<n> after their base property, once each",
            len(props) == len(p["properties"]["rows"]) and all(k == f"base:{r[1]}" and r[2] for k, r in props.items()))
    B = p["box"]["local"]
    assert_("city-events: every property's point in the box", all(B[0] <= r[4][0] <= B[2] and B[1] <= r[4][1] <= B[3] for r in props.values()))
    assert_("city-events: joins are inside, nearest within 15 m, or none (with no building)",
            all((r[6] == "inside" and r[5]) or (r[6] == "none" and r[5] is None)
                or (re.fullmatch(r"nearest:\d+\.\d", r[6]) and r[5] and float(r[6][8:]) <= ce.JOIN_M) for r in props.values()))
    assert_("city-events: censuses and establishments name properties of the record",
            all(r[1] in props for r in p["censuses"]["rows"]) and all(r[1] in props for r in p["establishments"]["rows"]))
    assert_("city-events: census years within the window", all(2002 <= r[0] <= 2024 for t in ("censuses", "establishments") for r in p[t]["rows"]))
    codes = {r[4] for r in p["establishments"]["rows"] if r[4]}
    assert_("city-events: every class that occurs is described, with its division by the ABS's ranges",
            codes <= set(p["anzsic"]) and all(v["division"] == ce.anzsic_levels(c)["division"] for c, v in p["anzsic"].items()))
    ev = [dict(zip(p["events"]["columns"], r)) for r in p["events"]["rows"]]
    types = {"birth", "death", "alteration", "institution", "regime", "occurrence"}
    iso = re.compile(r"\d{4}-\d{2}-\d{2}")
    assert_("city-events: events typed, recorded, labelled by intervals [from, to) within the record",
            all(e["type"] in types and e["provenance"] == "recorded" and iso.fullmatch(e["from"]) and iso.fullmatch(e["to"])
                and e["from"] < e["to"] and "1952-01-01" <= e["from"] < p["window"]["to"] for e in ev))
    assert_("city-events: before_window marks exactly the events before 2002",
            all(bool((e["detail"] or {}).get("before_window")) == (e["from"] < p["window"]["from"]) for e in ev))
    assert_("city-events: event ids unique", len({e["id"] for e in ev}) == len(ev))
    allowed = {"property", "building", "rule", "streets", "point"}
    assert_("city-events: matches name properties of the record, rules and streets",
            all(e["match"] and set(e["match"]) <= allowed and (e["match"].get("property") is None or e["match"]["property"] in props)
                and isinstance(e["match"].get("streets", []), list) for e in ev))
    names = set(p["report"]["district_street_names"])
    assert_("city-events: permits match streets of the district", all(set(e["match"].get("streets", [])) <= names for e in ev))
    same = {e["detail"]["same_as"] for e in ev if e["detail"] and e["detail"].get("same_as")}
    assert_("city-events: a completion said to be the same fact names an event of the record", same <= {e["id"] for e in ev}, str(same - {e["id"] for e in ev}))
    rc = p["report"]["counts"]
    assert_("city-events: the report counts its tables", rc == {t: len(p[t]["rows"]) for t in ("properties", "censuses", "establishments", "events")}, str(rc))
    j = p["report"]["joins"]
    assert_("city-events: the joins' outcomes sum to the properties", j["inside"] + j["nearest"] + j["none"] == rc["properties"])
    model_path = PRODUCTS / p["city_model"]["file"]
    if model_path.exists():
        cm = json.loads(model_path.read_text(encoding="utf-8"))
        ids = {b["id"] for b in cm["buildings"]}
        used = {r[5] for r in props.values() if r[5]} | {e["match"]["building"] for e in ev if e["match"].get("building")}
        assert_("city-events: buildings are the city model's, and the model is the one the record was built on",
                used <= ids and cm["origin"] == p["origin"] and cm.get("built") == p["city_model"]["built"], f"{len(used - ids)} unknown")
    areas = ROOT / "applications" / "urban-science-workbench" / "areas.json"
    if areas.exists():
        d = next(a for a in json.loads(areas.read_text(encoding="utf-8"))["areas"] if a["id"] == "melbourne")["district"]
        assert_("city-events: the district is Landscape Dynamics' (areas.json)", all(abs(d[k] - ce.DISTRICT[k]) < 1e-9 for k in ("south", "west", "north", "east")))
    # Seen by eye on 2026-11-04: the Town Hall's registration on the Victorian Heritage Register.
    th = next((e for e in ev if e["id"] == "vhr-registered:H0001"), None)
    assert_("city-events: Melbourne Town Hall (H0001) registered on 9 October 1974, on its property",
            th is not None and th["from"] == "1974-10-09" and th["type"] == "institution" and "property" in th["match"], str(th))
else:
    skip("city-events: not built (python data/products/city_events.py)")
