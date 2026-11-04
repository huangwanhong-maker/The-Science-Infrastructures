"""Write the maintenance module's fixture from the data library's product road-condition-melbourne.json: each rated
surface of the Hoddle Grid with its type, age (empty if unknown), class, rating, area and evaluation date.

    python data/acquisition/fetch.py melbourne-road-surfaces
    python data/products/road_condition.py
    python domains/urban-science/maintenance/tests/make_fixtures.py
"""

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
PRODUCT = ROOT / "data" / "cache" / "products" / "road-condition-melbourne.json"
OUT = Path(__file__).resolve().parent / "fixtures" / "road-condition-hoddle-grid.csv"


def main() -> int:
    d = json.loads(PRODUCT.read_text(encoding="utf-8"))
    c = d["columns"]
    idx = {k: c.index(k) for k in ("type", "age", "class", "rating", "area", "evaluated")}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as f:
        src = d["source"]
        f.write(f"# {src['title']}\n# {src['credit']}; licence: {src['licence']}; retrieved {src['retrieved']}\n")
        f.write(f"# {d['rated']} rated surfaces of {d['rows']}; ages withheld: {json.dumps(d['without_age'])}\n")
        w = csv.writer(f)
        w.writerow(["type", "age", "class", "rating", "area", "evaluated"])
        for s in d["surfaces"]:
            w.writerow([s[idx["type"]], "" if s[idx["age"]] is None else s[idx["age"]], s[idx["class"]], s[idx["rating"]], s[idx["area"]], s[idx["evaluated"]]])
    print(f"wrote {OUT} ({d['rated']} surfaces)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
