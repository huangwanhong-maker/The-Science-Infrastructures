"""Tables for the Economy Atlas's toolbox (ADR 0056): the data library's tables (data/products/economic_tables.py, moved there by ADR 0080)
and tables of one's own.

A table is named columns of equal length, each of numbers (null where missing) or of text:

    {"id", "title", "source", "rows", "columns": [{"name", "kind": "number"|"text", "label"?}],
     "data": {name: [values]}}

Tables of one's own are CSV files (comma, semicolon or tab separated, a header row), kept as the
datasets of kind "table" (datasets.py); a column is numbers when every filled cell is a number.

    python applications/economy-atlas/tables.py list
    python applications/economy-atlas/tables.py make     (the recipe's step: fetch, then make)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

APP = Path(__file__).resolve().parent
sys.path.insert(0, str(APP))
import datasets  # noqa: E402  (which puts data/products on the path)
from economic_tables import KEPT, library, make  # noqa: E402,F401


def own() -> list[dict]:
    """The tables of one's own, kept as datasets of kind "table"."""
    out = []
    store = datasets.Datasets()
    for d in store.list():
        if d.get("kind") == "table":
            try:
                t = json.loads(store.product(d["id"], "table").read_text(encoding="utf-8"))
                out.append({**t, "id": d["id"], "own": True})
            except (OSError, ValueError, datasets.DatasetError):
                continue
    return out


def listing() -> list[dict]:
    """Every table, without its data."""
    return [{k: v for k, v in t.items() if k not in ("data", "stamp")} for t in library() + own()]


def table(ident: str) -> dict:
    for t in library() + own():
        if t["id"] == ident:
            return {k: v for k, v in t.items() if k != "stamp"}
    raise datasets.DatasetError(f"no table {ident!r}")


if __name__ == "__main__":
    if sys.argv[1:2] == ["list"]:
        for t in listing():
            print(t["id"], t["rows"], [c["name"] for c in t["columns"]])
    elif sys.argv[1:2] == ["make"]:
        sys.exit(make())
    else:
        print(__doc__)
