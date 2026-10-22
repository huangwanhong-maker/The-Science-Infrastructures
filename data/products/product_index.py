"""The index of data products, data/cache/products/index.json, that applications read
to learn which products exist. Each kind of product (stars, galaxies, ...) has a
list of entries; an entry is replaced when a product of the same name is rebuilt.
"""

from __future__ import annotations

import json
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "cache" / "products"


def register(kind: str, entry: dict, first: str | None = None) -> None:
    """Adds or replaces entry (by its 'name') in the list for kind."""
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "index.json"
    index = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    items = [e for e in index.get(kind, []) if e["name"] != entry["name"]] + [entry]
    items.sort(key=lambda e: (e["name"] != first, e["name"]))
    index[kind] = items
    path.write_text(json.dumps(index, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
