"""Checks of the Economy Atlas's example plug-ins (plugins/, ADR 0046), run through
the API's host of plug-ins over networks of places (api/python/science/networks.py)
with the include directories of the atlas's engine: the concentration of exports
against its closed form, the gravity equation's coefficients against a least-squares
fit of the same equation made here, the stationary stocks of a two-place network
against their closed form (q, 2 − q) and summing to n on any; and, where the trade
flows are built, all three on the world's trade of 2024 (a distance elasticity
between 0.8 and 1.5, indices between 0 and 1, stocks summing to the economies).
One JSON object a line, as api/tests/check_client.py reports; a C++ plug-in is
skipped, and said, where no compiler is found.

    python applications/economy-atlas/tests/check_plugins.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
ROOT = APP.parents[1]
sys.path.insert(0, str(ROOT / "api" / "python"))
sys.path.insert(0, str(APP))
import numpy as np  # noqa: E402

import engine  # noqa: E402
from science import hosting, networks  # noqa: E402


def emit(record: dict) -> None:
    print(json.dumps(record, ensure_ascii=False), flush=True)


def check(name: str, value: float, limit: float, unit: str = "") -> None:
    emit({"kind": "check", "name": name, "value": float(value), "limit": float(limit), "unit": unit})


def assert_(name: str, condition: bool, detail: str = "") -> None:
    emit({"kind": "assert", "name": name, "ok": bool(condition), "detail": str(detail)})


def skip(what: str) -> None:
    emit({"kind": "skip", "name": what})


def place(pid, lat, lon):
    return {"id": pid, "name": pid, "latitude": lat, "longitude": lon}


FIVE = {"title": "five places", "unit": "m $",
        "places": [place("N", 90, 0), place("E", 0, 0), place("F", 0, 90), place("S", -30, 45), place("W", 20, -60)],
        "flows": [{"from": a, "to": b, "value": v} for a, b, v in
                  [("E", "F", 10), ("F", "E", 4), ("N", "E", 3), ("S", "N", 2), ("E", "S", 1), ("W", "E", 5), ("F", "W", 2),
                   ("N", "F", 1.5), ("S", "W", 3), ("W", "N", 1.2)]]}
TWO = {"title": "two places", "places": [place("a", 0, 0), place("b", 0, 10)], "flows": [{"from": "a", "to": "b", "value": 1}]}


def by_kind(results, kind, title=None):
    return [r for r in results if r["kind"] == kind and (title is None or r["title"] == title)][0]


def main() -> int:
    host = networks.Host(APP / "plugins", ROOT / "data" / "cache" / "economy-atlas" / "plugins", engine.include_dirs())
    files = [f["file"] for f in host.files()]
    assert_("plug-ins: the three examples are listed", {"plugins/gravity.py", "plugins/concentration.py", "plugins/stationary.cpp"} <= set(files), files)
    net = networks.Network(FIVE)
    F, X = net.matrix(), net.exports()

    # Concentration: H_i = Σ_j (F_ij / X_i)².
    r = host.run("plugins/concentration.py", FIVE, {"threshold": 0})
    H = by_kind(r, "values")["values"]
    want = {pid: float(((F[i] / X[i]) ** 2).sum()) for i, pid in enumerate(net.ids) if X[i] > 0}
    check("concentration: every economy's index is Σ (F_ij / X_i)²", max(abs(H[k] - want[k]) for k in want), 1e-12)
    check("concentration: the Equator's, (10² + 1²) / 11², by hand", H["E"] - 101 / 121, 1e-12)

    # Gravity: the coefficients of ln F = c + a ln X + b ln M − d ln D, fitted here independently.
    r = host.run("plugins/gravity.py", FIVE, {"minimum": 0, "arcs": 3})
    M, D = net.imports(), np.maximum(net.distances_km(), 100.0)
    i, j = np.nonzero(F > 0)
    A = np.column_stack([np.ones(len(i)), np.log(X[i]), np.log(M[j]), -np.log(D[i, j])])
    coef = np.linalg.lstsq(A, np.log(F[i, j]), rcond=None)[0]
    got = [row[1] for row in by_kind(r, "table", "The gravity equation")["rows"]]
    check("gravity: its coefficients equal a least-squares fit of the same equation made here", np.abs(np.array(got) - coef).max(), 1e-9)
    residual = np.log(F[i, j]) - A @ coef
    mean = {net.ids[k]: residual[i == k].mean() for k in set(i.tolist())}
    vals = by_kind(r, "values")["values"]
    check("gravity: each economy's value is the mean residual of its exports", max(abs(vals[k] - mean[k]) for k in mean), 1e-9)
    arcs = by_kind(r, "arcs")["arcs"]
    assert_("gravity: at most the arcs asked for, each above the line", len(arcs) <= 3 and all(
        math.log(a["value"]) - A[[n for n in range(len(i)) if net.ids[i[n]] == a["from"] and net.ids[j[n]] == a["to"]][0]] @ coef > 0
        for a in arcs), arcs)

    # Stationary stocks: a → b only, b keeping its stock; every column of M' is (q/2, 1 − q/2), so π = 2 (q/2, 1 − q/2).
    if hosting.compiler() is None:
        skip("stationary: no C++ compiler (c++, g++ or clang++, or g++ in WSL), so no C++ plug-in")
    else:
        for q in (0.15, 0.4):
            r = host.run("plugins/stationary.cpp", TWO, {"q": q})
            pi = by_kind(r, "values")["values"]
            check(f"stationary: two places, a passing everything to b, q = {q}: π = (q, 2 − q)", abs(pi["a"] - q) + abs(pi["b"] - (2 - q)), 1e-12)
        r = host.run("plugins/stationary.cpp", FIVE, {})
        pi = by_kind(r, "values")["values"]
        check("stationary: on five places the stocks sum to n", sum(pi.values()) - 5, 1e-9)

    # The world's trade of 2024, where the product is built.
    trade = ROOT / "data" / "cache" / "products" / "trade-flows.json"
    if not trade.exists():
        skip("plug-ins: the trade flows are not built, so the world's trade is not run")
        return 0
    product = json.loads(trade.read_text(encoding="utf-8"))
    E = product["economies"]
    world = {"title": "World trade in goods, 2024", "year": 2024, "unit": product.get("unit", ""),
             "places": [{"id": e["iso3"], "name": e["name"], "latitude": e["latitude"], "longitude": e["longitude"],
                         "group": e.get("continent", "")} for e in E],
             "flows": [{"from": E[a]["iso3"], "to": E[b]["iso3"], "value": v} for a, b, v in product["flows"]["2024"]]}
    r = host.run("plugins/gravity.py", world, {})
    d = by_kind(r, "table", "The gravity equation")["rows"][3][1]
    assert_("gravity, 2024: the distance elasticity d lies between 0.8 and 1.5, as the literature finds", 0.8 < d < 1.5, d)
    H = by_kind(host.run("plugins/concentration.py", world, {}), "values")["values"]
    assert_("concentration, 2024: every index between 0 and 1", all(0 < h <= 1 for h in H.values()), (min(H.values()), max(H.values())))
    if hosting.compiler() is not None:
        pi = by_kind(host.run("plugins/stationary.cpp", world, {}), "values")["values"]
        check("stationary, 2024: the stocks sum to the economies", sum(pi.values()) - len(E), 1e-6)
    return 0


if __name__ == "__main__":
    sys.exit(main())
