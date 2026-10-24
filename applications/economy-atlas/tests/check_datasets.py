"""Checks of the Economy Atlas's datasets of one's own (datasets.py, ADR 0046) and of what the engine
makes of them: CSV and JSON read into the library's schemas, every spelling of an economy (CHN,
chn, China) one place, places of one's own with semicolons and a byte-order mark, years kept apart,
and the refusals (an unknown id, a negative flow, missing columns, a word for a number, two places at
one point); a record of currencies on a base of its own. Where the engine is built from the present
sources: the map of a dataset of one's own (its totals adding up, its attraction summing to zero),
one of three places, which makes no map's faces, the spacetime of a record on the US dollar (its
lens's curvature in closed form), and the strip of years of the library's trade (every year's shares
adding to 1). One JSON object a line, as api/tests/check_client.py reports.

    python applications/economy-atlas/tests/check_datasets.py
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
ROOT = APP.parents[1]
sys.path.insert(0, str(APP))
import datasets as ds  # noqa: E402
import engine  # noqa: E402


def report(kind: str, name: str, **fields) -> None:
    print(json.dumps({"kind": kind, "name": name, **fields}, ensure_ascii=False), flush=True)


def check(name: str, value: float, limit: float) -> None:
    report("check", name, value=value, limit=limit)


def require(name: str, ok: bool, detail: str = "") -> None:
    report("assert", name, ok=bool(ok), detail=detail)


def refused(text_files: list[dict], kind: str = "flows") -> str:
    try:
        (ds.make_flows if kind == "flows" else ds.make_currencies)("x", text_files)
    except ds.DatasetError as error:
        return str(error)
    return ""


def ask(request: dict) -> dict:
    done = subprocess.run([sys.executable, str(APP / "engine.py"), "run", "-"], input=json.dumps(request).encode("utf-8"),
                          capture_output=True, cwd=ROOT)
    try:
        return json.loads(done.stdout.decode("utf-8"))
    except ValueError:
        return {"error": (done.stdout + done.stderr).decode("utf-8", "replace").replace("\0", "")[-300:]}


def main() -> int:
    library = ds.TRADE.exists()
    # Codes, names and cases of one economy are one place.
    if library:
        p = ds.make_flows("G3", [{"name": "f.csv", "text": "from,to,value,year\nUSA,CHN,100,2020\nCHN,USA,400,2020\nDEU,USA,90,2020\nusa,deu,60,2021\nChina,Germany,30,2021\n"}])
        require("datasets: CHN, chn and China are one place, and the years are kept apart",
                [e["iso3"] for e in p["economies"]] == ["CHN", "DEU", "USA"] and p["years"] == [2020, 2021] and len(p["flows"]["2021"]) == 2,
                json.dumps(p["flows"]))
    else:
        report("skip", "datasets: the library's economies are not built: ids by code and name are not checked")
    p = ds.make_flows("Cities", [
        {"name": "places.csv", "text": "﻿id;name;latitude;longitude;group\nA;Alpha;10;20;x\nB;Beta;-5;40;x\nC;Gamma;30;-60;y\nD;Delta;50;120;y\n"},
        {"name": "flows.csv", "text": "from;to;value\nA;B;5\nB;C;3\nC;A;4\nD;A;7\nA;A;9\nB;A;2\n"}])
    require("datasets: places of one's own (semicolons, a byte-order mark), flows to oneself left out, two flows of a pair kept",
            len(p["economies"]) == 4 and p["years"] == [0] and sorted(map(tuple, p["flows"]["0"])) == [(0, 1, 5.0), (1, 0, 2.0), (1, 2, 3.0), (2, 0, 4.0), (3, 0, 7.0)],
            json.dumps(p["flows"]))
    j = ds.make_flows("", [{"name": "x.json", "text": json.dumps({"title": "J", "unit": "t", "places": [{"id": "p", "name": "P", "latitude": 0, "longitude": 0},
                                                                                                         {"id": "q", "name": "Q", "latitude": 0, "longitude": 90}],
                                                                  "flows": [{"from": "p", "to": "q", "value": 2, "year": 2001}]})}])
    require("datasets: JSON read, with its title, unit and year", j["title"] == "J" and j["unit"] == "t" and j["flows"] == {"2001": [[0, 1, 2.0]]})
    cases = [("an unknown id", [{"name": "f.csv", "text": "from,to,value\nNOWHERE,ELSEWHERE,1\n"}], "no place"),
             ("a negative flow", [{"name": "f.csv", "text": "from,to,value\nA,B,-1\n"}, {"name": "p.csv", "text": "id,name,latitude,longitude\nA,A,0,0\nB,B,0,10\n"}], "negative"),
             ("missing columns", [{"name": "f.csv", "text": "a,b\n1,2\n"}], "columns"),
             ("a word for a number", [{"name": "f.csv", "text": "from,to,value\nA,B,ten\n"}], "not a number"),
             ("two places at one point", [{"name": "p.csv", "text": "id,name,latitude,longitude\nA,A,1,1\nB,B,1,1\n"}, {"name": "f.csv", "text": "from,to,value\nA,B,1\n"}], "same point")]
    for what, files, expected in cases:
        message = refused(files)
        require(f"datasets: {what} is refused, and the refusal says so", expected in message, message)
    c = ds.make_currencies("Toy", [{"name": "r.csv", "text": "month,currency,rate,policy\n2020-01,USD,,1.5\n2020-02,USD,,1.5\n2020-03,USD,,1.0\n"
                                                             "2020-01,JPY,108,-0.1\n2020-02,JPY,110,-0.1\n2020-03,JPY,107,-0.1\n"}])
    require("datasets: a record of currencies finds its base (the one with no rate), its months and rates",
            c["base"] == "USD" and c["months"] == ["2020-01", "2020-02", "2020-03"] and c["rates"] == {"JPY": [108.0, 110.0, 107.0]})

    if not engine.status()["current"]:
        report("skip", "the engine is not built from the present sources: datasets through it are not checked")
        return 0
    kind = (engine.read_build() or {}).get("toolchain") or "native"
    with tempfile.TemporaryDirectory() as tmp:
        store = ds.Datasets(Path(tmp))
        six = store.add("flows", "Six", [
            {"name": "p.csv", "text": "id,name,latitude,longitude\nA,A,0,0\nB,B,0,60\nC,C,0,120\nD,D,0,180\nE,E,0,-120\nF,F,0,-60\n"},
            {"name": "f.csv", "text": "from,to,value\nA,B,10\nB,C,9\nC,D,8\nD,E,7\nE,F,6\nF,A,5\nA,D,20\nD,A,4\nB,E,3\n"}])
        m = ask({"op": "map", "flows": engine.path_for(kind, store.product(six["id"], "flows"))})
        T = m.get("totals", {})
        require("datasets: the map of six places of one's own", "error" not in m and len(m.get("economies", [])) == 6, m.get("error", ""))
        if "error" not in m:
            check("... exchange and net transfer add up to the whole", (T["reciprocal"] + T["net"] - T["total"]) / T["total"], 1e-12)
            check("... the shares of the net flows add to 1", T["driven"] + T["local"] + T["global"] - 1, 1e-9)
            check("... the attraction sums to zero", sum(e["attraction"] for e in m["economies"]), 1e-5)
        three = store.add("flows", "Three", [{"name": "p.csv", "text": "id,name,latitude,longitude\nA,A,0,0\nB,B,0,60\nC,C,30,30\n"},
                                             {"name": "f.csv", "text": "from,to,value\nA,B,10\nB,C,9\nC,A,8\n"}])
        m3 = ask({"op": "map", "flows": engine.path_for(kind, store.product(three["id"], "flows"))})
        require("datasets: three places make no map's faces, and their circulation is global",
                "error" not in m3 and m3.get("faces", {}).get("triangles") == 0 and m3["totals"]["global"] > 0.99, json.dumps(m3)[:200])
        usd = store.add("currencies", "Dollar", [{"name": "r.csv", "text": "month,currency,rate,policy\n2020-01,USD,,6\n2020-02,USD,,6\n2020-01,JPY,100,0\n2020-02,JPY,101,0\n"}])
        s = ask({"op": "spacetime", "record": engine.path_for(kind, store.product(usd["id"], "currencies"))})
        if "error" in s:
            require("datasets: a record on the US dollar in spacetime", False, s["error"])
        else:
            dt = (2020 + 1.5 / 12) - (2020 + 0.5 / 12)
            kappa = (math.log1p(0.06) - 0.0) * dt + math.log(101) - math.log(100)
            require("datasets: a record on the US dollar in spacetime, its base first", s["base"] == "USD" and s["currencies"] == ["USD", "JPY"])
            check("... its lens's curvature (r_USD − r_JPY) Δt + ln 101 − ln 100", s["lenses"][0][3] - kappa, 1e-12)
    if library:
        t = ask({"op": "timeline", "flows": engine.path_for(kind, ds.TRADE)})
        years = t.get("years", [])
        product_years = json.loads(ds.TRADE.read_text(encoding="utf-8"))["years"]
        require("the strip of years: one entry for each year of the library's trade", [y["year"] for y in years] == product_years, t.get("error", ""))
        if years:
            check("... every year's shares of the net flows add to 1", max(abs(y["driven"] + y["local"] + y["global"] - 1) for y in years), 1e-9)
    else:
        report("skip", "the trade flows are not built: the strip of years is not checked")
    return 0


if __name__ == "__main__":
    sys.exit(main())
