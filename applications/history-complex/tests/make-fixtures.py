"""Writes the fixtures of the History Complex from the built engine: for every canon of the history
domain, the answer to ``build`` and the moments at two dates, so that the page can be seen with
``?fixtures`` and its checks run without the engine. Also the list of canons, as the server gives it.

    python applications/history-complex/tests/make-fixtures.py
"""

from __future__ import annotations


import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = HERE.parent
sys.path.insert(0, str(APP))
import engine  # noqa: E402
import serve  # noqa: E402

DATES = ["1808-10-01", "1814-11-15"]
OUT = HERE / "fixtures"


def compact(x):
    """Eight significant digits for the fixtures' reals (a day is 0.0027 of a year); integers exact."""
    if isinstance(x, float):
        return int(x) if x.is_integer() else float(f"{x:.8g}")
    if isinstance(x, list):
        return [compact(v) for v in x]
    if isinstance(x, dict):
        return {k: compact(v) for k, v in x.items()}
    return x


def main() -> int:
    sys.stdout.reconfigure(encoding="utf-8")
    if not engine.status()["current"] and not engine.build():
        print("The engine is not built from the present sources.", file=sys.stderr)
        return 1
    kind = (engine.read_build() or {}).get("toolchain") or "native"
    OUT.mkdir(parents=True, exist_ok=True)
    for old in OUT.glob("*.json"):
        old.unlink()
    canons = serve.canons()
    (OUT / "canons.json").write_text(json.dumps({"canons": canons}, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    for c in canons:
        paths = [engine.path_for(kind, p) for p in serve.canon_paths(c["id"])]
        where = {"canons": paths} if c.get("integrated") else {"canon": paths[0]}
        answer = engine.ask(kind, {"op": "build", **where})
        if "error" in answer:
            print(f"{c['id']}: {answer['error']}", file=sys.stderr)
            return 1
        (OUT / f"build-{c['id']}.json").write_text(json.dumps(compact(answer), ensure_ascii=False) + "\n", encoding="utf-8")
        for date in DATES:
            m = engine.ask(kind, {"op": "moment", **where, "date": date})
            (OUT / f"moment-{c['id']}-{date}.json").write_text(json.dumps(compact(m), ensure_ascii=False) + "\n", encoding="utf-8")
        if c is canons[0]:
            power = engine.ask(kind, {"op": "power", **where})
            (OUT / f"power-{c['id']}.json").write_text(json.dumps(compact(power), ensure_ascii=False) + "\n", encoding="utf-8")
        inv = answer["invariants"]
        print(f"{c['id']}: {len(answer['beings'])} beings, {len(answer['occasions'])} occasions, {len(answer['ties'])} ties; "
              f"V E F = {inv['counts']}, χ = {inv['euler']}, {len(inv['loops'])} loops, {len(inv['shells'])} shells")
    return 0


if __name__ == "__main__":
    sys.exit(main())
