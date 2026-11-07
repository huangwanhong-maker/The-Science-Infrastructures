"""Checks of the plug-ins over data in space and time (science/spacetime.py and cpp/spacetime_plugin.hpp, ADR 0079), run by api/tests/run-tests.mjs: a record of each form
read into columns, a missing value NaN; results checked, and refused when their columns differ in length; bins of time against closed forms; a plug-in written to a directory of its
own, described and run through Host, in Python and, where a compiler is found, in C++, with the same results on the same records; a plug-in that raises reported with its reason; a
file outside the directory refused. One JSON object a line, as check_client.py reports.

    python api/tests/check_spacetime.py
"""

from __future__ import annotations

import json
import math
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "api" / "python"))
import numpy as np  # noqa: E402

from science import hosting, spacetime as st  # noqa: E402


def emit(record: dict) -> None:
    print(json.dumps(record, ensure_ascii=False), flush=True)


def check(name: str, value: float, limit: float) -> None:
    emit({"kind": "check", "name": name, "value": float(value), "limit": float(limit), "unit": ""})


def assert_(name: str, condition: bool, detail: str = "") -> None:
    emit({"kind": "assert", "name": name, "ok": bool(condition), "detail": str(detail)})


def skip(what: str) -> None:
    emit({"kind": "skip", "name": what})


# ---- records in columns ------------------------------------------------------------------------------------------------------------------------------------------------
e = st.Record({"kind": "events", "title": "e", "time": [1000, 2000], "lat": [1, None], "lon": [3, 4], "values": {"v": [5, None]}, "units": {"v": "m"}})
assert_("spacetime api: events read into columns of floats", e.time.dtype == float and e.time[1] == 2000 and e.units["v"] == "m")
assert_("spacetime api: a missing value is NaN", math.isnan(e.lat[1]) and math.isnan(e.values["v"][1]) and e.values["v"][0] == 5)
t = st.Record({"kind": "tracks", "tracks": [{"id": "a", "time": [1, 2], "lat": [0, 1], "lon": [0, 1], "values": {"w": [7, 8]}}]})
assert_("spacetime api: tracks, each with its columns", t.tracks[0].id == "a" and t.tracks[0].values["w"][1] == 8)
c = st.Record({"kind": "cells", "time": [1, 2, 3], "cells": [{"id": "x", "ring": [[0, 0], [1, 0], [0, 1]]}, {"id": "y", "ring": [[1, 1], [2, 1], [1, 2]]}], "values": {"n": [[1, 2], [3, None], [5, 6]]}})
assert_("spacetime api: cells, values by time and cell", c.values["n"].shape == (3, 2) and c.values["n"][2, 1] == 6 and math.isnan(c.values["n"][1, 1]) and c.ids == ["x", "y"])
g = st.Record({"kind": "grid", "time": [1, 2], "box": [0, 0, 1, 1], "width": 3, "height": 2, "values": {"z": [[0, 1, 2, 3, 4, 5], [5, 4, 3, 2, 1, 0]]}})
assert_("spacetime api: a grid, values by time, row from the south, column", g.values["z"].shape == (2, 2, 3) and g.values["z"][0, 1, 2] == 5)
f = st.Record({"kind": "flows", "places": [{"id": "a", "lat": 0, "lon": 0}], "flows": []})
assert_("spacetime api: flows keep their places", f.places[0]["id"] == "a")
try:
    st.Record({"kind": "clouds"})
    assert_("spacetime api: a form that does not exist is refused", False)
except ValueError:
    assert_("spacetime api: a form that does not exist is refused", True)
index, starts = st.bins_of([0.0, 5.0, 10.0, 25.0, float("nan")], 10.0)
assert_("spacetime api: bins of time: each time's bin and the bins' starts; no time, no bin", index.tolist() == [0, 0, 1, 2, -1] and starts.tolist() == [0.0, 10.0, 20.0], (index.tolist(), starts.tolist()))
back = st.grid("g", [1, 2], [0, 0, 1, 1], {"z": np.arange(12.0).reshape(2, 2, 3)})
assert_("spacetime api: a grid answered is flattened, its size said", back["width"] == 3 and back["height"] == 2 and back["values"]["z"].shape == (2, 6))
done = st.finish([st.events("a", [1.0], [2.0], [float("nan")]), st.note("n")])
assert_("spacetime api: results made plain for JSON, NaN as null", done[0]["lon"] == [None] and done[1]["text"] == "n")
for bad, why in ([st.events("a", [1, 2], [1], [1])], "columns of different lengths"), ([{"kind": "smoke"}], "a kind that does not exist"):
    try:
        st.finish(bad)
        assert_(f"spacetime api: {why} is refused", False)
    except ValueError:
        assert_(f"spacetime api: {why} is refused", True)

# ---- through the host ---------------------------------------------------------------------------------------------------------------------------------------------------
PY = '''
import numpy as np
def create(api):
    class Late(api.SpaceTimeAnalysis):
        id = "late"; title = "The last ones"; description = "Keeps the last n events, and counts them by a step."
        params = [api.param("n", "How many", 2, kind="integer", minimum=1)]
        def run(self, records, params, context):
            r = records[0]
            order = np.argsort(r.time)[-int(params["n"]):]
            context.progress(0.5, "halfway")
            return [api.events("Late", r.time[order], r.lat[order], r.lon[order], values={"v": r.values["v"][order]}), api.note(f"at {context.now_ms}")]
    return Late()
'''
CPP = r'''
#include "spacetime_plugin.hpp"
#include <algorithm>
#include <numeric>
namespace st = science::spacetime;
struct Late : st::Analysis {
  science::Json describe() const override {
    return science::Json::object({{"id", "late"}, {"title", "The last ones"}, {"description", "Keeps the last n events."}, {"params", science::Json::array({st::param("n", "How many", 2, "integer", 1)})}});
  }
  std::vector<science::Json> run(const std::vector<st::Record>& records, const st::Params& p, const st::Context& context, st::Progress& progress) override {
    const st::Events& e = records.at(0).events;
    std::vector<std::size_t> order(e.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return e.time[a] < e.time[b]; });
    progress(0.5, "halfway");
    st::Events out;
    const auto n = static_cast<std::size_t>(p.integer("n", 2));
    for (std::size_t k = order.size() - std::min(n, order.size()); k < order.size(); k++) {
      out.push(e.time[order[k]], e.lat[order[k]], e.lon[order[k]]);
      out.values["v"].push_back(e.values.at("v")[order[k]]);
    }
    return {out.to_json("Late"), st::note("at " + std::to_string(static_cast<long long>(context.now_ms)))};
  }
};
SCIENCE_SPACETIME_PLUGIN(Late)
'''
work = Path(tempfile.mkdtemp(prefix="spacetime-check-"))
try:
    scripts = work / "scripts"
    scripts.mkdir()
    (scripts / "late.py").write_text(PY, encoding="utf-8")
    (scripts / "raises.py").write_text("def create(api):\n    class R(api.SpaceTimeAnalysis):\n        id = 'r'\n        def run(self, records, params, context):\n            raise ValueError('it failed on purpose')\n    return R()\n", encoding="utf-8")
    host = st.Host(scripts, work / "work")
    d = host.describe("late.py")
    assert_("spacetime host: a script describes itself", d["id"] == "late" and d["kind"] == "spacetime-analysis" and d["params"][0]["name"] == "n", d)
    record = {"kind": "events", "title": "r", "time": [30, 10, 20, 40], "lat": [3, 1, 2, 4], "lon": [3, 1, 2, 4], "values": {"v": [300, 100, None, 400]}}
    seen = []
    out = host.run("late.py", {"records": [record], "context": {"now_ms": 77}}, {"n": 3}, progress=lambda f, x: seen.append((f, x)))
    assert_("spacetime host: a script runs on the records, with its parameter", out[0]["time"] == [20, 30, 40] and out[0]["values"]["v"] == [None, 300, 400], out[0])
    assert_("spacetime host: and the view, and its progress is followed", out[1]["text"] == "at 77" and seen and seen[0][1] == "halfway", (out[1], seen))
    defaulted = host.run("late.py", {"records": [record], "context": {}}, {})
    assert_("spacetime host: a parameter not given takes its default", len(defaulted[0]["time"]) == 2)
    try:
        host.run("raises.py", {"records": [], "context": {}}, {})
        assert_("spacetime host: a script's failure is reported with its reason", False)
    except RuntimeError as error:
        assert_("spacetime host: a script's failure is reported with its reason", "on purpose" in str(error), str(error)[-200:])
    try:
        host.plugin_path("../outside.py")
        assert_("spacetime host: a file outside the directory is refused", False)
    except ValueError:
        assert_("spacetime host: a file outside the directory is refused", True)
    if hosting.compiler():
        (scripts / "late.cpp").write_text(CPP, encoding="utf-8")
        d = host.describe("late.cpp")
        assert_("spacetime host: a C++ program describes itself", d["id"] == "late" and d["language"] == "cpp" and d["kind"] == "spacetime-analysis", d)
        same = host.run("late.cpp", {"records": [record], "context": {"now_ms": 77}}, {"n": 3})
        assert_("spacetime host: C++ answers what Python answers, on the same records", same[0]["time"] == out[0]["time"] and same[0]["lat"] == out[0]["lat"] and same[0]["values"]["v"] == out[0]["values"]["v"] and same[1]["text"] == out[1]["text"], same)
    else:
        skip("no C++ compiler: the C++ plug-in is not run")
finally:
    shutil.rmtree(work, ignore_errors=True)
