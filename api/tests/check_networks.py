"""Checks of the plug-ins over networks of places (science/networks.py and
cpp/network_plugin.hpp, ADR 0046), run by api/tests/run-tests.mjs: a network's
matrix, exports, imports, balance and distances against closed forms; values
completed by the places' ids and NaN left out; parameters defaulted and typed,
and refused out of bounds; a plug-in written to a plugins/ directory of its own
described and run through Host, in Python and, where a compiler is found, in C++
with the same results; a plug-in that raises reported with its reason; a file
outside plugins/ refused; a run stopped by its time limit. One JSON object a
line, as check_client.py reports.

    python api/tests/check_networks.py
"""

from __future__ import annotations

import json
import math
import shutil
import sys
import tempfile
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "api" / "python"))
import numpy as np  # noqa: E402

from science import hosting, networks  # noqa: E402


def emit(record: dict) -> None:
    print(json.dumps(record, ensure_ascii=False), flush=True)


def check(name: str, value: float, limit: float, unit: str = "") -> None:
    emit({"kind": "check", "name": name, "value": float(value), "limit": float(limit), "unit": unit})


def assert_(name: str, condition: bool, detail: str = "") -> None:
    emit({"kind": "assert", "name": name, "ok": bool(condition), "detail": str(detail)})


def skip(what: str) -> None:
    emit({"kind": "skip", "name": what})


# Four places: the north pole, two on the equator a quarter-turn apart, and one in the south.
NETWORK = {
    "title": "four places", "year": 2024, "unit": "m $",
    "places": [{"id": "N", "name": "North", "latitude": 90, "longitude": 0},
               {"id": "E", "name": "Equator", "latitude": 0, "longitude": 0},
               {"id": "F", "name": "Far", "latitude": 0, "longitude": 90},
               {"id": "S", "name": "South", "latitude": -30, "longitude": 45, "group": "south"}],
    "flows": [{"from": a, "to": b, "value": v} for a, b, v in
              [("E", "F", 10), ("F", "E", 4), ("N", "E", 3), ("S", "N", 2), ("E", "S", 1), ("E", "F", 0.5)]],
}

TINY_PY = '''
def create(api):
    class Tiny(api.NetworkAnalysis):
        id = "tiny"
        title = "Tiny"
        description = "balances, scaled; the largest flow"
        params = [api.param("scale", "Scale", 2.0, minimum=0), api.param("count", "Count", 3, kind="integer"),
                  api.param("mode", "Mode", "a", kind="choice", options=["a", "b"]), api.param("flag", "Flag", False, kind="boolean")]

        def run(self, network, params, context):
            context.progress(0.5, "half way")
            b = network.balance() * params["scale"]
            k = int(network.value.argmax())
            return [api.values("balance", list(b), scale="diverging"),
                    api.arcs("largest", [{"from": network.ids[network.source[k]], "to": network.ids[network.target[k]],
                                           "value": float(network.value[k])}]),
                    api.table("params", ["name", "value", "type"], [[n, v, type(v).__name__] for n, v in params.items()]),
                    api.note("done")]
    return Tiny()
'''

TINY_CPP = r'''
#include "network_plugin.hpp"
namespace net = science::network;
using science::Json;

struct Tiny : net::Analysis {
  Json describe() const override {
    return Json::object({{"id", "tiny"}, {"title", "Tiny"}, {"description", "balances, scaled; the largest flow"},
                         {"params", Json::array({net::param("scale", "Scale", 2.0, "number", 0.0),
                                                 net::param("count", "Count", 3, "integer")})}});
  }
  std::vector<Json> run(const net::Network& n, const net::Params& p, net::Progress& progress) override {
    progress(0.5, "half way");
    std::map<std::string, double> b;
    const auto balance = n.balance();
    for (std::size_t i = 0; i < n.size(); i++) b[n.places[i].id] = balance[i] * p.number("scale");
    std::size_t k = 0;
    for (std::size_t f = 0; f < n.flows.size(); f++) if (n.flows[f].value > n.flows[k].value) k = f;
    const auto& f = n.flows[k];
    return {net::values("balance", b), net::arcs("largest", {{n.places[f.from].id, n.places[f.to].id, f.value}}),
            net::note("count " + std::to_string(p.integer("count")) + ", d(N, E) " + Json::number_text(n.distance_km(0, 1)))};
  }
};
SCIENCE_NETWORK_PLUGIN(Tiny)
'''

RAISES_PY = '''
def create(api):
    class Raises(api.NetworkAnalysis):
        id = "raises"
        def run(self, network, params, context):
            raise ValueError("the network is too small for this")
    return Raises()
'''

SLOW_PY = '''
import time
def create(api):
    class Slow(api.NetworkAnalysis):
        id = "slow"
        def run(self, network, params, context):
            time.sleep(30)
            return []
    return Slow()
'''


def main() -> int:
    R = networks.EARTH_RADIUS_KM
    n = networks.Network(NETWORK)
    F = n.matrix()
    assert_("networks: the matrix adds flows between the same two places, row from, column to",
            F[n.index["E"], n.index["F"]] == 10.5 and F[n.index["F"], n.index["E"]] == 4 and F.sum() == 20.5, F.tolist())
    want_x, want_m = [3, 11.5, 4, 2], [2, 7, 10.5, 1]
    check("networks: exports, imports and balance (imports − exports)",
          max(np.abs(n.exports() - want_x).max(), np.abs(n.imports() - want_m).max(),
              np.abs(n.balance() - (np.array(want_m) - want_x)).max()), 0)
    D = n.distances_km()
    check("networks: the pole to the equator and a quarter of the equator are πR/2", max(abs(D[0, 1] - math.pi * R / 2),
                                                                                         abs(D[1, 2] - math.pi * R / 2)), 1e-9, "km")
    check("networks: distances are symmetric with a zero diagonal", max(np.abs(D - D.T).max(), np.abs(np.diag(D)).max()), 1e-9, "km")
    try:
        networks.Network({**NETWORK, "flows": [{"from": "E", "to": "X", "value": 1}]})
        assert_("networks: a flow to an unknown place is refused", False)
    except ValueError as error:
        assert_("networks: a flow to an unknown place is refused, naming it", "'X'" in str(error), str(error))
    done = networks.finish([networks.values("v", [1.0, float("nan"), None, 4.0])], n)[0]
    assert_("networks: values given as a list are keyed by the places' ids, NaN and None left out",
            done["values"] == {"N": 1.0, "S": 4.0}, done["values"])
    described = [networks.param("k", "k", 3, kind="integer", minimum=1), networks.param("x", "x", 0.5),
                 networks.param("m", "m", "a", kind="choice", options=["a", {"value": "b", "label": "B"}]),
                 networks.param("f", "f", False, kind="boolean")]
    typed = networks.typed(described, {"k": "40", "m": "b", "f": "true"})
    assert_("networks: parameters defaulted and typed by kind", typed == {"k": 40, "x": 0.5, "m": "b", "f": True}
            and isinstance(typed["k"], int), typed)
    refused = []
    for given in ({"k": 0}, {"k": 2.5}, {"m": "c"}, {"x": "nan"}):
        try:
            networks.typed(described, given)
        except ValueError as error:
            refused.append(str(error))
    assert_("networks: a parameter below its least, not whole, not an option or not finite is refused", len(refused) == 4, refused)

    # Through a host, on a plugins/ directory of its own.
    tmp = Path(tempfile.mkdtemp(prefix="science-networks-"))
    try:
        plugins = tmp / "plugins"
        plugins.mkdir()
        (plugins / "tiny.py").write_text(TINY_PY, encoding="utf-8")
        (plugins / "raises.py").write_text(RAISES_PY, encoding="utf-8")
        (plugins / "slow.py").write_text(SLOW_PY, encoding="utf-8")
        (plugins / "_hidden.py").write_text(TINY_PY, encoding="utf-8")
        (tmp / "outside.py").write_text(TINY_PY, encoding="utf-8")
        host = networks.Host(plugins, tmp / "work")
        listed = [f["file"] for f in host.files()]
        assert_("networks host: the plug-ins listed, those starting with _ left out", listed == ["plugins/raises.py", "plugins/slow.py", "plugins/tiny.py"], listed)
        d = host.describe("plugins/tiny.py")
        assert_("networks host: a Python plug-in describes itself", d["id"] == "tiny" and d["language"] == "python"
                and [p["name"] for p in d["params"]] == ["scale", "count", "mode", "flag"] and d["file"] == "plugins/tiny.py", d)
        seen = []
        results = host.run("plugins/tiny.py", NETWORK, {"scale": 3}, progress=lambda f, t: seen.append((f, t)))
        by = {r["kind"]: r for r in results}
        want = {pid: 3 * b for pid, b in zip(n.ids, n.balance())}
        check("networks host: a Python run's values, by id", max(abs(by["values"]["values"][k] - want[k]) for k in want), 1e-12)
        types = {row[0]: row[2] for row in by["table"]["rows"]}
        assert_("networks host: its parameters reached it defaulted and typed",
                types == {"scale": "float", "count": "int", "mode": "str", "flag": "bool"}
                and by["arcs"]["arcs"] == [{"from": "E", "to": "F", "value": 10.0}], by["table"]["rows"])
        assert_("networks host: its progress followed", (0.5, "half way") in seen, seen)
        try:
            host.run("plugins/raises.py", NETWORK, {})
            assert_("networks host: a plug-in that raises is reported", False)
        except RuntimeError as error:
            assert_("networks host: a plug-in that raises is reported with its reason", "too small for this" in str(error), str(error)[-200:])
        refused = []
        for file in ("../outside.py", "plugins/../outside.py", str(tmp / "outside.py"), "plugins/missing.py"):
            try:
                host.plugin_path(file)
            except ValueError:
                refused.append(file)
        assert_("networks host: a file outside plugins/, or none, is refused", len(refused) == 4, refused)
        started = time.time()
        try:
            host.run("plugins/slow.py", NETWORK, {}, seconds=2)
            assert_("networks host: a run past its time limit is stopped", False)
        except RuntimeError as error:
            assert_("networks host: a run past its time limit is stopped", "stopped after" in str(error) and time.time() - started < 15,
                    f"{error} after {time.time() - started:.1f} s")
        cancel = threading.Event()
        threading.Timer(1.0, cancel.set).start()
        try:
            host.run("plugins/slow.py", NETWORK, {}, cancel=cancel)
            assert_("networks host: a cancelled run stops", False)
        except RuntimeError as error:
            assert_("networks host: a cancelled run stops", "cancelled" in str(error), str(error))

        # The same plug-in in C++, where a compiler is found.
        if hosting.compiler() is None:
            skip("networks host: no C++ compiler (c++, g++ or clang++, or g++ in WSL), so no C++ plug-in")
        else:
            (plugins / "tiny_cpp.cpp").write_text(TINY_CPP, encoding="utf-8")
            started = time.time()
            dc = host.describe("plugins/tiny_cpp.cpp")
            compiled = time.time() - started
            assert_("networks host: a C++ plug-in compiles and describes itself", dc["id"] == "tiny" and dc["language"] == "cpp"
                    and dc["kind"] == "network-analysis" and [p["name"] for p in dc["params"]] == ["scale", "count"],
                    f"{dc} ({compiled:.0f} s to compile)")
            seen = []
            rc = host.run("plugins/tiny_cpp.cpp", NETWORK, {"scale": 3}, progress=lambda f, t: seen.append((f, t)))
            byc = {r["kind"]: r for r in rc}
            check("networks host: the C++ run's values equal the Python run's",
                  max(abs(byc["values"]["values"][k] - by["values"]["values"][k]) for k in want), 1e-12)
            assert_("networks host: ... and its arc, its progress and its parameters' defaults",
                    byc["arcs"]["arcs"] == by["arcs"]["arcs"] and (0.5, "half way") in seen
                    and byc["note"]["text"].startswith("count 3, d(N, E) 10007.5"), (byc["arcs"], byc["note"], seen))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
