"""The areas in dispute, for what colours the world by economy: Natural Earth's breakaway and disputed areas as polygons, without their names, since
those are not drawn (ADR 0069). The regions product assigns every point of the land to a country by Natural Earth's polygons, which give each such
area to whoever holds it; a work that draws countries' values asks this as well, so that no claimant's colour falls on them.

    python data/products/disputed_areas.py

writes data/cache/products/areas-in-dispute.json: {polygons: [[ring, ...], ...]}, each ring of [lon, lat], read by the even-odd rule.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from naming import Areas, DISPUTED  # noqa: E402

OUT = HERE.parent / "cache" / "products" / "areas-in-dispute.json"


def main() -> int:
    areas = Areas.load()
    if areas is None:
        print(f"{DISPUTED} is not fetched: python data/acquisition/fetch.py natural-earth-50m", file=sys.stderr)
        return 1
    polygons = [[[[round(x, 5), round(y, 5)] for x, y in ring] for ring in poly] for _, poly in areas.polygons]
    OUT.write_text(json.dumps({"schema": "areas-in-dispute/1", "source": "Natural Earth 1:50m, ne_50m_admin_0_breakaway_disputed_areas (public domain); names not kept", "polygons": polygons}, separators=(",", ":")), encoding="utf-8")
    print(f"{len(polygons)} polygons, to {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
