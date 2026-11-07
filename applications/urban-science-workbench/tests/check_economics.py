"""Checks of Quantitative economics (ADR 0080), run by tests/run-tests.mjs: every question of the plug-in answered by its engine on the
tables of the data library, as the page sends them; a tool for each place held to the same tool run on each part alone; residuals by
row against the rows the tool kept; location quotients by hand; and the tables against their sources. One JSON object a line.

    python applications/urban-science-workbench/tests/check_economics.py
"""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = HERE.parent
ROOT = APP.parents[1]
sys.path.insert(0, str(APP))
sys.path.insert(0, str(ROOT / "data" / "acquisition"))
import engine  # noqa: E402

TABLES = ROOT / "data" / "cache" / "products" / "tables"
CACHE = ROOT / "data" / "cache"
SHOW = "--show" in sys.argv


def emit(record: dict) -> None:
    print(json.dumps(record, ensure_ascii=False), flush=True)


def check(name: str, value: float, limit: float) -> None:
    emit({"kind": "check", "name": name, "value": float(value), "limit": float(limit), "unit": ""})


def assert_(name: str, condition: bool, detail: str = "") -> None:
    emit({"kind": "assert", "name": name, "ok": bool(condition), "detail": str(detail)})


def skip(what: str) -> None:
    emit({"kind": "skip", "name": what})


status = engine.status("economics")
if not status["current"]:
    skip("the engine of Quantitative economics is not built from the present sources (python applications/urban-science-workbench/engine.py build economics)")
    sys.exit(0)
worker = engine.Worker("economics")


def ask(request: dict) -> dict:
    a = worker.ask({**request, "op": "tool"})
    if not a.get("ok"):
        raise RuntimeError(a.get("error"))
    return a["result"]


def table(name: str) -> dict | None:
    p = TABLES / f"{name}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def positions(t: dict) -> tuple[list, list]:
    P = t.get("places") or {}
    if P.get("kind") == "cells":
        centre = {k: (sum(q[0] for q in r) / len(r), sum(q[1] for q in r) / len(r)) for k, r in P["rings"].items()}
    elif P.get("kind") == "points":
        centre = {k: tuple(v) for k, v in P["at"].items()}
    else:
        return [], []
    ids = t["data"][P["column"]]
    return [centre.get(str(i), (None, None))[0] for i in ids], [centre.get(str(i), (None, None))[1] for i in ids]


def request_of(q: dict, t: dict | None) -> tuple[dict, list | None]:
    """What the page sends for a question: the columns its roles name, in the rows its filter keeps."""
    body = {"tool": q["tool"], "roles": q.get("roles", {}), "params": q.get("params", {})}
    if t is None:
        return body, None
    data = dict(t["data"])
    lon, lat = positions(t)
    data["place_lon"], data["place_lat"] = lon, lat
    n = t["rows"]
    filt = q.get("filter", {})
    kept = [i for i in range(n) if all(data[c][i] == v for c, v in filt.items())]
    wanted = set(q.get("each", []))
    for v in body["roles"].values():
        wanted |= set(v if isinstance(v, list) else [v])
    body["columns"] = {c: [data[c][i] for i in kept] for c in wanted}
    if q.get("each"):
        body["each"] = q["each"]
    if q.get("byRow"):
        body["by_row"] = True
        if lon:
            body["places"] = {"lon": [lon[i] for i in kept], "lat": [lat[i] for i in kept]}
    return body, kept


def figure(answer: dict, label: str) -> float | None:
    for b in answer["blocks"]:
        if b["kind"] == "figures":
            for item in b["items"]:
                if item[0] == label:
                    return item[1]
    return None


def block(answer: dict, kind: str, title: str | None = None) -> dict | None:
    for b in answer["blocks"]:
        if b["kind"] == kind and (title is None or title in b.get("title", "")):
            return b
    return None


def questions() -> list[dict]:
    code = "import('" + (APP / "plugins" / "economics" / "questions.js").as_uri() + "').then((m) => console.log(JSON.stringify(m.QUESTIONS)))"
    r = subprocess.run(["node", "-e", code], capture_output=True, text=True, encoding="utf-8")
    if r.returncode:
        raise RuntimeError(r.stderr)
    return json.loads(r.stdout)


# ---- every question answers, on the tables held ------------------------------------------------------------------------------------
kept_tiles = sorted((CACHE / "space2stats").glob("*.json")) if (CACHE / "space2stats").exists() else []
for q in questions():
    name = f"economics: \"{q['ask']}\""
    t = None
    if q.get("table"):
        if q["table"].startswith("cells"):
            skip(f"{name}: its table is made in the page, from the tile looked at")
            continue
        t = table(q["table"])
        if t is None:
            skip(f"{name}: the table {q['table']} is not made")
            continue
    body, kept = request_of(q, t)
    try:
        a = ask(body)
    except RuntimeError as error:
        assert_(f"{name} is answered", False, str(error))
        continue
    drawable = [b for b in a["blocks"] if b["kind"] in ("rows", "map", "each")]
    names = [b.get("title", "") for b in drawable] + [m["name"] for b in drawable if b["kind"] == "each" for m in b["measures"]]
    wanted = q.get("draw")
    assert_(f"{name} is answered, and what it draws first is there", not wanted or any(wanted in n for n in names), f"wanted {wanted!r} among {names[:12]}")
    if SHOW:
        print("\n==", q["ask"], file=sys.stderr)
        for b in a["blocks"]:
            if b["kind"] == "figures":
                print("  ", b["title"], {i[0]: i[1] for i in b["items"]}, file=sys.stderr)
            elif b["kind"] == "table":
                print("  ", b["title"], b["rows"][:6], file=sys.stderr)
            elif b["kind"] == "note":
                print("   note:", b["text"][:170], file=sys.stderr)
            elif b["kind"] == "each":
                m = next((m for m in b["measures"] if wanted and wanted in m["name"]), b["measures"][0])
                vals = sorted((v, k) for v, k in zip(m["values"], b["keys"]) if v is not None)
                print("   each", m["name"], "n", len(vals), "low", vals[:3], "high", vals[-3:], "left", b["left"][:2], file=sys.stderr)
            elif b["kind"] == "rows":
                vals = [v for v in b["values"] if v is not None]
                print("   rows", b["title"][:60], "n", len(vals), "range", min(vals), max(vals), file=sys.stderr)
            elif b["kind"] == "map":
                print("   map", b["title"][:70], "n", len(b["values"]), file=sys.stderr)

# ---- a tool for each place is that tool on each part ---------------------------------------------------------------------------------
blocks_t = table("melbourne-blocks")
if blocks_t:
    d = blocks_t["data"]
    years = sorted(set(d["census_year"]))
    whole = ask({"tool": "inequality", "roles": {"value": "establishments"}, "params": {}, "each": ["census_year"], "columns": {"establishments": d["establishments"], "census_year": d["census_year"]}})
    each_b = whole["blocks"][0]
    assert_("for each: one part a year, in the order they first appear", each_b["kind"] == "each" and sorted(int(k[0]) for k in each_b["keys"]) == years, each_b["keys"][:3])
    gini = next(m for m in each_b["measures"] if m["name"].startswith("Gini"))
    worst = 0.0
    for key, value in zip(each_b["keys"], gini["values"]):
        y = int(key[0])
        one = ask({"tool": "inequality", "roles": {"value": "establishments"}, "params": {}, "columns": {"establishments": [v for v, c in zip(d["establishments"], d["census_year"]) if c == y]}})
        worst = max(worst, abs(figure(one, gini["name"]) - value))
    check("for each: each year's Gini is the tool's own on that year's rows", worst, 1e-12)

    # What cannot be answered is left out and counted, not given a zero.
    few = ask({"tool": "ols", "roles": {"y": "y", "x": ["x"]}, "params": {"covariance": "classical"}, "each": ["g"], "columns": {"g": ["a", "a", "a", "a", "b", "b"], "x": [1, 2, 3, 4, 1, 2], "y": [1.5, 2.1, 2.9, 4.2, 1, 2]}})
    fb = few["blocks"][0]
    estimate = next(m for m in fb["measures"] if m["name"] == "x · estimate")
    assert_("for each: a part with too few rows has no value, and is counted with its reason", fb["keys"][1] == ["b"] and estimate["values"][1] is None and fb["left"] and fb["left"][0]["count"] == 1 and "too few" in fb["left"][0]["reason"], fb["left"])
    # x̄ = 2.5, ȳ = 2.675: Σ(x − x̄)(y − ȳ) = 4.45 over Σ(x − x̄)² = 5.
    check("for each: the slope of the part that can be answered, by hand: 4.45/5 = 0.89", estimate["values"][0] - 0.89, 1e-9)

    # Residuals by row name the rows kept.
    x = [1.0, 2.0, None, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    y = [2.1, 3.9, 5.0, 8.2, 9.8, 12.1, 14.2, 15.8, 18.1, None]
    lon = [144.0 + 0.01 * i for i in range(10)]
    lat = [-37.0 - 0.01 * (i % 3) for i in range(10)]
    r = ask({"tool": "ols", "roles": {"y": "y", "x": ["x"]}, "params": {"covariance": "classical"}, "by_row": True, "columns": {"x": x, "y": y}, "places": {"lon": lon, "lat": lat}})
    res = block(r, "rows", "Residuals")
    assert_("residuals by row: the rows with a missing value have none", res["rows"] == [0, 1, 3, 4, 5, 6, 7, 8], res["rows"])
    check("residuals by row: they sum to zero with a constant", sum(res["values"]), 1e-10)
    fit = block(r, "rows", "Fitted")
    worst = max(abs(y[i] - f - e) for i, f, e in zip(res["rows"], fit["values"], res["values"]))
    check("residuals by row: outcome = fitted + residual, row by row", worst, 1e-12)
    assert_("residuals are tested across space where rows are places", figure(r, "Moran's I of the residuals") is not None)
    plain = ask({"tool": "ols", "roles": {"y": "y", "x": ["x"]}, "params": {"covariance": "classical"}, "columns": {"x": x, "y": y}})
    assert_("a request that does not ask for rows is answered as before, with none", block(plain, "rows") is None and figure(plain, "Moran's I of the residuals") is None)
    repeated = ask({"tool": "ols", "roles": {"y": "y", "x": ["x"]}, "params": {"covariance": "classical"}, "by_row": True, "columns": {"x": x, "y": y}, "places": {"lon": [144.0] * 10, "lat": [-37.0] * 10}})
    assert_("a place repeated: its residuals are not tested across space, and it says so", figure(repeated, "Moran's I of the residuals") is None and any(b["kind"] == "note" and "more than once" in b["text"] for b in repeated["blocks"]))

    # The panel's effects by row.
    panel = ask({"tool": "panel", "roles": {"unit": "u", "time": "t", "y": "y", "x": ["x"]}, "params": {"effects": "one-way", "random": False}, "by_row": True,
                 "columns": {"u": ["a", "a", "a", "b", "b", "b"], "t": [1, 2, 3, 1, 2, 3], "x": [1.0, 2.0, 3.0, 1.0, 2.0, 3.0], "y": [11.0, 12.0, 13.0, 21.0, 22.0, 23.0]}})
    effect = block(panel, "rows", "effect")
    check("a panel's unit effects by row, about their mean: −5 and +5 for y = 10 + x and 20 + x", max(abs(effect["values"][k] - v) for k, v in enumerate([-5, -5, -5, 5, 5, 5])), 1e-9)

# ---- location quotients by hand -------------------------------------------------------------------------------------------------------
lq = ask({"tool": "specialisation", "roles": {"columns": ["a", "b"]}, "params": {"focus": "a"}, "columns": {"a": [10, 30, 0, None], "b": [30, 30, 0, 40]}})
rows = block(lq, "rows")
# Totals: a 40, b 100, all 140: place 0 holds a's share 10/40 over 40/140; place 1 30/60 over 40/140; the third holds nothing, and has no share.
check("location quotients by hand: (10/40)/(40/140) = 0.875", rows["values"][0] - (10 / 40) / (40 / 140), 1e-12)
check("location quotients by hand: (30/60)/(40/140) = 1.75", rows["values"][1] - (30 / 60) / (40 / 140), 1e-12)
assert_("location quotients: a place with nothing of these kinds has no quotient, an empty cell holds none", rows["rows"] == [0, 1, 3] and rows["values"][2] == 0, rows["rows"])

# ---- Moran's I through the tool ------------------------------------------------------------------------------------------------------------
lons = [144.0 + 0.01 * c for r_ in range(6) for c in range(6)]
lats = [-37.0 + 0.01 * r_ for r_ in range(6) for c in range(6)]
smooth = [lons[i] * 100 + lats[i] * 50 for i in range(36)]
m = ask({"tool": "moran", "roles": {"value": "v", "lon": "x", "lat": "y"}, "params": {"k": 4, "permutations": 199}, "columns": {"v": smooth, "x": lons, "y": lats}})
assert_("Moran's tool: a smooth surface has I well above zero, and p at its least", figure(m, "I") > 0.5 and figure(m, "p, by permutation") <= 1 / 200 + 1e-12, figure(m, "I"))
assert_("Moran's tool: clusters for each place", len(block(m, "rows", "Clusters")["values"]) == 36)
try:
    ask({"tool": "moran", "roles": {"value": "v", "lon": "x", "lat": "y"}, "params": {}, "columns": {"v": list(range(10)), "x": [1.0] * 10, "y": [2.0] * 10}})
    assert_("Moran's tool: a place appearing more than once is refused", False)
except RuntimeError as error:
    assert_("Moran's tool: a place appearing more than once is refused", "more than once" in str(error), error)

# ---- the tables against their sources -------------------------------------------------------------------------------------------------
if blocks_t:
    d = blocks_t["data"]
    floor = CACHE / "melbourne-clue-floor-space" / "all" / "floor-space.csv"
    with floor.open(encoding="utf-8-sig", newline="") as f:
        rows_src = list(csv.DictReader(f, delimiter=";"))
    pick = next(r for r in rows_src if r["census_year"] == "2024" and r["block_id"] == "4")
    at = next(i for i in range(blocks_t["rows"]) if d["census_year"][i] == 2024 and d["block_id"][i] == "4")
    check("melbourne-blocks: a block's retail is the census's four kinds summed", d["retail"][at] - sum(float(pick[c] or 0) for c in ("retail_shop", "retail_showroom", "retail_stall", "retail_cars")), 1e-9)
    check("melbourne-blocks: its total is the census's", d["floor_total"][at] - float(pick["total_floor_space_in_block"]), 1e-9)
    withheld = sum(1 for r in rows_src if r["office"] == "")
    assert_("melbourne-blocks: a withheld cell stays null, not zero", sum(1 for v in d["office"] if v is None) >= withheld * 0.99 and withheld > 0, withheld)
    with (CACHE / "melbourne-clue-establishments" / "all" / "establishments.csv").open(encoding="utf-8-sig", newline="") as f:
        count = sum(1 for r in csv.DictReader(f, delimiter=";") if r["census_year"] == "2024" and r["block_id"] == "4")
    check("melbourne-blocks: its establishments are the census's rows for it", d["establishments"][at] - count, 0)

day = table("melbourne-counters-day")
if day:
    counts = json.loads((ROOT / "data/cache/products/street-counts/1.json").read_text(encoding="utf-8"))
    # A civil day from the hours by hand: AEDT, UTC+11, on 30 October 2024 (the offsets the counts carry).
    from datetime import datetime, timezone
    t0 = datetime.fromisoformat(counts["t0_utc"].replace("Z", "+00:00"))
    start = datetime(2024, 10, 29, 13, tzinfo=timezone.utc)
    i0 = int((start - t0).total_seconds() // 3600)
    total = sum(counts["counts"][i0 + k] for k in range(24))
    at = next(i for i in range(day["rows"]) if day["data"]["counter"][i] == "1" and day["data"]["date"][i] == "2024-10-30")
    check("melbourne-counters-day: a day is its 24 civil hours summed by hand", day["data"]["walkers"][at] - total, 1e-9)
    assert_("melbourne-counters-day: no day of nobody is kept as a count", min(day["data"]["walkers"]) > 0)
    weekend = [w for w, ds in zip(day["data"]["weekend"], day["data"]["date"]) if ds == "2024-11-02"]
    assert_("melbourne-counters-day: 2 November 2024 was a Saturday", weekend and all(w == 1 for w in weekend))

month = table("melbourne-counters-month")
if month:
    gap = [m_ for m_ in set(month["data"]["month"]) if "2022-12" <= m_ <= "2024-09"]
    assert_("melbourne-counters-month: nothing stands in the months the City did not publish", not gap, sorted(gap)[:3])

wards = table("freetown-wards")
if wards:
    source = json.loads((ROOT / "data/cache/products/freetown-flood.json").read_text(encoding="utf-8"))
    check("freetown-wards: the wards' losses sum to the assessment's", sum(v or 0 for v in wards["data"]["annual_loss"]) - sum(w.get("aal") or 0 for w in source["wards"] if w.get("polygons")), 1e-6)

worker.stop()
