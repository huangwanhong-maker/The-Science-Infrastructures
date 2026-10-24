"""Checks of the economics products (data/products/trade_flows.py, currency_record.py), offline: the
codes and places of Natural Earth's countries (France and Norway through ISO_A3_EH, Australia's
territories yielding to Australia), the reading of BIMTS rows (exports of all goods only, in
millions), the reading of the ECB's and the BIS's rows, the months between two; and, where the
products have been built, their shape. Nothing is fetched.
"""

from __future__ import annotations

import json

from checks import PRODUCTS, assert_, skip

import currency_record as cr  # noqa: E402
import trade_flows as tf  # noqa: E402

# ---- Natural Earth's countries ---------------------------------------------------------------------
assert_("trade flows: France's code through ISO_A3_EH", tf.economy_code({"ISO_A3": "-99", "ISO_A3_EH": "FRA"}) == "FRA")
assert_("trade flows: Somaliland has no code", tf.economy_code({"ISO_A3": "-99", "ISO_A3_EH": "-99"}) is None)
features = {"features": [
    {"properties": {"ISO_A3": "-99", "ISO_A3_EH": "AUS", "NAME": "Indian Ocean Ter.", "CONTINENT": "Asia", "LABEL_X": 105.0, "LABEL_Y": -10.5}},
    {"properties": {"ISO_A3": "AUS", "ISO_A3_EH": "AUS", "NAME": "Australia", "NAME_LONG": "Australia", "CONTINENT": "Oceania", "LABEL_X": 134.5, "LABEL_Y": -25.7}},
    {"properties": {"ISO_A3": "-99", "ISO_A3_EH": "AUS", "NAME": "Ashmore and Cartier Is.", "CONTINENT": "Oceania", "LABEL_X": 123.6, "LABEL_Y": -12.4}},
]}
p = tf.places(features)
assert_("trade flows: Australia's place is Australia's, not its territories'", p["AUS"]["name"] == "Australia" and p["AUS"]["longitude"] == 134.5, json.dumps(p))

# ---- BIMTS rows ---------------------------------------------------------------------------------------
rows = [
    {"TRADE_FLOW": "X", "PRODUCT_HS": "_T", "REF_AREA": "USA", "COUNTERPART_AREA": "CHN", "TIME_PERIOD": "2022", "OBS_VALUE": "157125.78", "UNIT_MULT": "6"},
    {"TRADE_FLOW": "X", "PRODUCT_HS": "HS17_88", "REF_AREA": "USA", "COUNTERPART_AREA": "CHN", "TIME_PERIOD": "2022", "OBS_VALUE": "4352.9", "UNIT_MULT": "6"},
    {"TRADE_FLOW": "M", "PRODUCT_HS": "_T", "REF_AREA": "USA", "COUNTERPART_AREA": "CHN", "TIME_PERIOD": "2022", "OBS_VALUE": "536307.1", "UNIT_MULT": "6"},
    {"TRADE_FLOW": "X", "PRODUCT_HS": "_T", "REF_AREA": "CHN", "COUNTERPART_AREA": "USA", "TIME_PERIOD": "2022", "OBS_VALUE": "500", "UNIT_MULT": "9"},
    {"TRADE_FLOW": "X", "PRODUCT_HS": "_T", "REF_AREA": "CHN", "COUNTERPART_AREA": "DEU", "TIME_PERIOD": "2022", "OBS_VALUE": "", "UNIT_MULT": "6"},
]
f = tf.read_flows(rows)
assert_("trade flows: only exports of all goods are read", set(f[2022]) == {("USA", "CHN"), ("CHN", "USA")}, str(f))
assert_("trade flows: values in millions of dollars, whatever UNIT_MULT says", f[2022][("CHN", "USA")] == 500000 and f[2022][("USA", "CHN")] == 157125.78)

# ---- the ECB's and the BIS's rows -------------------------------------------------------------------
r = cr.read_rates([{"CURRENCY": "USD", "CURRENCY_DENOM": "EUR", "EXR_SUFFIX": "A", "TIME_PERIOD": "2026-08", "OBS_VALUE": "1.1593"},
                   {"CURRENCY": "USD", "CURRENCY_DENOM": "EUR", "EXR_SUFFIX": "E", "TIME_PERIOD": "2026-08", "OBS_VALUE": "1.2"}])
assert_("currency record: the ECB's monthly averages, not its end-of-period rates", r == {"USD": {"2026-08": 1.1593}}, str(r))
b = cr.read_policy([{"FREQ": "M", "REF_AREA": "US", "TIME_PERIOD": "2024-01", "OBS_VALUE": "5.375"}, {"FREQ": "M", "REF_AREA": "JP", "TIME_PERIOD": "2024-01", "OBS_VALUE": ""}])
assert_("currency record: the BIS's policy rates, a missing value left out", b == {"US": {"2024-01": 5.375}}, str(b))
dfr = cr.read_deposit_rate([{"TIME_PERIOD": "2019-09-18", "OBS_VALUE": "-0.5"}, {"TIME_PERIOD": "2016-03-16", "OBS_VALUE": "-0.4"},
                            {"TIME_PERIOD": "2022-07-27", "OBS_VALUE": "0"}])
assert_("currency record: the deposit facility rate in force at the end of each month, from its dates of change",
        [cr.in_force(dfr, m) for m in ("2016-02", "2016-03", "2019-08", "2019-09", "2022-07", "2026-01")] == [None, -0.4, -0.4, -0.5, 0.0, 0.0],
        str([cr.in_force(dfr, m) for m in ("2016-02", "2016-03", "2019-08", "2019-09", "2022-07", "2026-01")]))
assert_("currency record: the months between two, across a year", cr.months_between("2025-11", "2026-02") == ["2025-11", "2025-12", "2026-01", "2026-02"])
assert_("currency record: every currency has its central bank's area", all(c in cr.AREA for c in cr.NAME))

# ---- the products, where built --------------------------------------------------------------------
trade = PRODUCTS / "trade-flows.json"
if trade.exists():
    t = json.loads(trade.read_text(encoding="utf-8"))
    n = len(t["economies"])
    ok = all(0 <= i < n and 0 <= j < n and i != j and v > 0 for flows in t["flows"].values() for i, j, v in flows)
    assert_("trade-flows: every flow joins two distinct economies with a positive value", ok, f"{n} economies")
    assert_("trade-flows: every economy has a place on the sphere", all(-90 <= e["latitude"] <= 90 and -180 <= e["longitude"] <= 180 for e in t["economies"]))
else:
    skip("trade-flows: not built (python data/products/trade_flows.py)")
cur = PRODUCTS / "currency-spacetime.json"
if cur.exists():
    c = json.loads(cur.read_text(encoding="utf-8"))
    m = len(c["months"])
    assert_("currency-spacetime: one rate and one policy rate a month for each currency",
            all(len(v) == m for v in c["rates"].values()) and all(len(v) == m for v in c["policy"].values()) and "EUR" not in c["rates"])
    assert_("currency-spacetime: rates positive where given", all(x is None or x > 0 for v in c["rates"].values() for x in v))
else:
    skip("currency-spacetime: not built (python data/products/currency_record.py)")
