"""Checks of the data library's naming policy (data/products/naming.json, naming.py; ADR 0069), offline:
every entry says its source, every name of the policy and every alias resolves to its code, a
territory whose status is in dispute is not labelled, a name in dispute is replaced in every language
the policy gives (Natural Earth's own Chinese and Japanese for Taiwan name the Republic of China), and
an unlisted place keeps its names; that Natural Earth's areas in dispute are found where they are, a
place in one losing its country and none being listed as a national capital there or in a territory
not labelled; and, where the products have been built, that they name places by the policy: the
trade flows, the geography (its countries, seas and places), and the gazetteer.
"""

from __future__ import annotations

import json

from checks import PRODUCTS, assert_, skip

from naming import AREAS, POLICY, Naming  # noqa: E402

names = POLICY.names
assert_("naming: every entry says its source", all(e.get("source") for e in names.values()), str([c for c, e in names.items() if not e.get("source")]))
# (Kosovo and Western Sahara have two codes each, of their data and of Natural Earth's polygon, with one name.)
wrong = [(c, a) for c, e in names.items() for a in [e.get("name"), *e.get("aliases", [])] if a and names[POLICY.resolve(a)].get("name") != e.get("name")]
assert_("naming: every name and alias resolves to an entry of the same name", not wrong, str(wrong))
assert_("naming: Taiwan and Turkey resolve, as aliases", POLICY.resolve("Taiwan") == "TWN" and POLICY.resolve("turkey") == "TUR")
for code in ("TWN", "PSE", "XKX", "KOS", "ESH", "FLK", "CYN", "SOL", "KAS"):
    assert_(f"naming: {code} is not labelled on a map", not POLICY.labelled(code))
assert_("naming: France and China are labelled", POLICY.labelled("FRA") and POLICY.labelled("CHN"))
ne_taiwan = {"default": "Taiwan", "en": "Taiwan", "zh": "中华民国", "zht": "中華民國", "ja": "中華民国", "fr": "Taïwan", "long": "Taiwan"}
local = POLICY.local("TWN", ne_taiwan)
assert_("naming: Chinese Taipei in English, 中华台北 and 中華台北 in Chinese, and no language keeps Natural Earth's name",
        local["default"] == "Chinese Taipei" and local["zh"] == "中华台北" and local["zht"] == "中華台北"
        and "中华民国" not in local.values() and "中華民國" not in local.values() and "Taiwan" not in local.values(), json.dumps(local, ensure_ascii=False))
assert_("naming: a territory the policy names nothing has no names", POLICY.local("CYN", {"default": "N. Cyprus", "en": "N. Cyprus"}) == {})
assert_("naming: Türkiye keeps Natural Earth's other languages", POLICY.local("TUR", {"default": "Turkey", "en": "Turkey", "fr": "Turquie"}) == {"default": "Türkiye", "en": "Türkiye", "fr": "Turquie"})
assert_("naming: an unlisted place keeps its names", POLICY.local("FRA", {"default": "France"}) == {"default": "France"} and POLICY.name("FRA", "France") == "France")
assert_("naming: the code of a feature without ISO codes is Natural Earth's own", Naming.code({"ISO_A3": "-99", "ISO_A3_EH": "-99", "ADM0_A3": "KOS"}) == "KOS")
for sea in ("Sea of Japan", "Persian Gulf", "Gulf of Mexico", "South China Sea"):
    assert_(f"naming: the {sea} is not labelled", not POLICY.sea_labelled(sea))
assert_("naming: the Atlantic is", POLICY.sea_labelled("North Atlantic Ocean"))
assert_("naming: the places of the West Bank and Gaza (PSX) are named as PSE", POLICY.name("PSX", "Palestine") == "West Bank and Gaza" and not POLICY.labelled("PSX"))

# ---- areas in dispute and places in them ------------------------------------------------------------
if AREAS is None:
    skip("naming: Natural Earth's areas in dispute are not fetched")
else:
    inside = {"Crimea": (34.1, 45.0), "the Golan Heights": (35.75, 33.0), "Aksai Chin": (79.3, 35.2), "Arunachal Pradesh": (93.5, 28.0), "the southern Kurils": (147.9, 45.0)}
    outside = {"Kyiv": (30.52, 50.45), "Paris": (2.35, 48.86), "Beijing": (116.4, 39.9), "Delhi": (77.2, 28.6), "Tokyo": (139.7, 35.7)}
    assert_("naming: the areas in dispute hold Crimea, the Golan, Aksai Chin, Arunachal Pradesh and the southern Kurils",
            all(AREAS.contains(*q) for q in inside.values()), str([k for k, q in inside.items() if not AREAS.contains(*q)]))
    assert_("naming: and not Kyiv, Paris, Beijing, Delhi or Tokyo", not any(AREAS.contains(*q) for q in outside.values()), str([k for k, q in outside.items() if AREAS.contains(*q)]))
    assert_("naming: a place in an area in dispute has no country, and is no national capital",
            POLICY.place("Somewhere", "Admin-0 capital", ("RUS",), 34.1, 45.0, AREAS) == ("Populated place", False))
    assert_("naming: a capital of a territory not labelled is listed as a populated place, keeping its country",
            POLICY.place("Taipei", "Admin-0 capital", ("TWN",), 121.57, 25.04, AREAS) == ("Populated place", True))
    assert_("naming: Jerusalem has no country and is no capital; Nicosia is the capital of Cyprus",
            POLICY.place("Jerusalem", "Admin-0 capital", ("ISR",), 35.21, 31.78, AREAS) == ("Populated place", False)
            and POLICY.place("Nicosia", "Admin-0 capital", ("CYP",), 33.37, 35.17, AREAS) == ("Admin-0 capital", True))
    assert_("naming: Paris is unchanged", POLICY.place("Paris", "Admin-0 capital", ("FRA",), 2.35, 48.86, AREAS) == ("Admin-0 capital", True))

# ---- the products, where built ---------------------------------------------------------------------
trade = PRODUCTS / "trade-flows.json"
if trade.exists():
    by = {e["iso3"]: e["name"] for e in json.loads(trade.read_text(encoding="utf-8"))["economies"]}
    assert_("naming: the trade flows name TWN Chinese Taipei and PSE West Bank and Gaza", by.get("TWN") == "Chinese Taipei" and by.get("PSE") == "West Bank and Gaza", str({k: by.get(k) for k in ("TWN", "PSE")}))
    assert_("naming: no economy of the trade flows is named Taiwan", "Taiwan" not in by.values())
else:
    skip("naming: the trade flows are not built")
geo = PRODUCTS / "geography-ne50m.json"
if geo.exists():
    g = json.loads(geo.read_text(encoding="utf-8"))
    layers = g.get("meta", g).get("layers", {})
    countries = {c["id"]: c for c in layers.get("countries", [])}
    names_all = [v for c in countries.values() for v in c["names"].values()]
    assert_("naming: the geography names no country 中华民国, 中華民國 or Taiwan", not {"中华民国", "中華民國", "中華民国", "Taiwan"} & set(names_all))
    assert_("naming: its disputed territories are not labelled", all(countries[c]["labelled"] is False for c in ("TWN", "KOS", "CYN", "SOL", "SAH") if c in countries))
    assert_("naming: its seas in dispute are not labelled", all(not m["labelled"] for m in layers.get("marine", []) if m["names"].get("default") in POLICY.seas))
    places = {q["names"].get("default"): q for q in layers.get("places", [])}
    assert_("naming: its places: no capital in a territory not labelled, Jerusalem without a country and no capital, Nicosia Cyprus's",
            not [q for q in places.values() if q["capital"] and not POLICY.labelled(q["country"])]
            and places["Jerusalem"]["country"] is None and not places["Jerusalem"]["capital"] and places["Nicosia"]["capital"],
            str([n for n, q in places.items() if q["capital"] and not POLICY.labelled(q["country"])]))
else:
    skip("naming: the geography is not built")
gaz = PRODUCTS / "places" / "earth.json"
if gaz.exists():
    rows = json.loads(gaz.read_text(encoding="utf-8"))["places"]
    assert_("naming: the gazetteer lists no country Taiwan, Somaliland or N. Cyprus", not [r for r in rows if r[3] == "Country" and r[0] in ("Taiwan", "Somaliland", "Northern Cyprus", "N. Cyprus")])
    assert_("naming: and places Taipei in Chinese Taipei", any(r[0] == "Taipei" and r[4] == "Chinese Taipei" for r in rows))
    by_name = {r[0]: r for r in rows if r[3] != "Country"}
    assert_("naming: the gazetteer has no country row for a territory not labelled",
            not [r for r in rows if r[3] == "Country" and r[0] in ("Chinese Taipei", "West Bank and Gaza", "Kosovo (UNSCR 1244)", "Western Sahara", "Falkland Islands (Malvinas)")])
    assert_("naming: its national capitals are none of Taipei, Pristina, Hargeisa, Jerusalem or Laayoune",
            not [n for n in ("Taipei", "Pristina", "Hargeisa", "Jerusalem", "Laayoune") if by_name[n][3].startswith("Admin-0")])
    assert_("naming: Sevastopol, Srinagar and Jerusalem have no country", all(by_name[n][4] == "" for n in ("Sevastopol", "Srinagar", "Jerusalem")),
            str({n: by_name[n][4] for n in ("Sevastopol", "Srinagar", "Jerusalem")}))
else:
    skip("naming: the gazetteer is not built")
