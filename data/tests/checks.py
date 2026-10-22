"""Checks of the data library written in Python report through this module, one
JSON object per line, which data/tests/run-tests.mjs reads into the shared
harness (development/tools/checks.mjs):

    from checks import check, assert_, skip
    check("metar: Haneda's temperature is plausible", abs(t - 288), 25, "K")   # passes when value <= limit
    assert_("ibtracs: every track has a name", ok, "3 unnamed")
    skip("ecmwf: no eccodes in this environment")

A check script is data/tests/check_<subject>.py; it may read the cache
(data/cache/products/) and skips what is missing, as the JavaScript suites do.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# UTF-8 whatever the console: names carry °, ρ and − (on Windows a pipe defaults to cp1252).
sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[2]
PRODUCTS = ROOT / "data" / "cache" / "products"
sys.path.insert(0, str(ROOT / "data" / "acquisition"))
sys.path.insert(0, str(ROOT / "data" / "products"))
sys.path.insert(0, str(ROOT / "data" / "formats"))


def _emit(record: dict) -> None:
    print(json.dumps(record, ensure_ascii=False), flush=True)


def check(name: str, value: float, limit: float, unit: str = "") -> None:
    _emit({"kind": "check", "name": name, "value": float(value), "limit": float(limit), "unit": unit})


def assert_(name: str, condition: bool, detail: str = "") -> None:
    _emit({"kind": "assert", "name": name, "ok": bool(condition), "detail": str(detail)})


def skip(what: str) -> None:
    _emit({"kind": "skip", "name": what})
