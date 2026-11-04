"""Build the condition of the Hoddle Grid's road surfaces for the urban-science domain's maintenance model (ADR 0071):
each surface's latest condition rating with the age at which it was rated.

The City of Melbourne publishes, for every road surface it manages, its condition rating at its last evaluation (a real
number from 1, failed, to 5, good), the evaluation's date, the installation date and the design life. One rating per
surface, at one age, is current-status data: the deterioration of a class of surfaces is fitted from the spread of their
ages (the domain's maintenance module, on tyche's panel likelihood).

1. **Rated surfaces.** Rows with a rating and an evaluation date. Rows without are counted, not kept.
2. **Ages.** The evaluation date less the installation date, in years of 365.2425 days. An installation date that is
   empty, on 1 January or on 25 December of its year (the dataset's placeholders: 1990-01-01, 1980-12-25, 1990-12-25
   hold about 2,100 of the rated surfaces), or after the evaluation, gives no age; such surfaces are kept with their
   age empty and the reason counted, since their ratings still describe the stock now.
3. **Classes.** The rating rounded to the nearest whole class and held to [1, 5].
4. **Groups.** The surface's type (Road Kerb, Road Footway, ...) with its responsible authority; the domain fits one
   chain a type.

Output: data/cache/products/road-condition-melbourne.json.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sys
from collections import Counter
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
CACHE = DATA / "cache"
SOURCE = CACHE / "melbourne-road-surfaces" / "hoddle-grid"
REGISTRY = DATA / "registry"
OUT = CACHE / "products" / "road-condition-melbourne.json"
YEAR = 365.2425


def day(text: str) -> dt.date | None:
    try:
        return dt.date.fromisoformat(text[:10]) if text else None
    except ValueError:
        return None


def build() -> dict:
    rows = list(csv.DictReader((SOURCE / "surfaces.csv").open(encoding="utf-8-sig"), delimiter=";"))
    retrieval = json.loads((SOURCE / "retrieval.json").read_text(encoding="utf-8"))
    manifest = json.loads((REGISTRY / "melbourne-road-surfaces.json").read_text(encoding="utf-8"))
    dropped, no_age = Counter(), Counter()
    columns = ["id", "type", "material", "authority", "area", "rating", "class", "evaluated", "installed", "age", "design_life", "street", "lat", "lon"]
    surfaces = []
    for r in rows:
        rating = r["conditionrating"].strip()
        evaluated = day(r["evaluationdate"])
        if not rating or evaluated is None:
            dropped["no rating or no evaluation date"] += 1
            continue
        installed = day(r["installationdate"])
        age = None
        if installed is None:
            no_age["no installation date"] += 1
        elif (installed.month, installed.day) in ((1, 1), (12, 25)):
            no_age["installation date a placeholder (1 January or 25 December)"] += 1
        elif installed > evaluated:
            no_age["installed after its evaluation"] += 1
        else:
            age = round((evaluated - installed).days / YEAR, 4)
        x = float(rating)
        lat = lon = None
        if r["geo_point_2d"]:
            lat, lon = (float(v) for v in r["geo_point_2d"].split(","))
        surfaces.append([r["assetid"], r["type"], r["surfacematerial"], r["responsibleauthoritycommonname"], float(r["area"]) if r["area"] else None,
                         x, min(5, max(1, round(x))), evaluated.isoformat(), installed.isoformat() if installed else None, age,
                         int(r["designlife"]) if r["designlife"] else None, r["roadsegmentdescription"], lat, lon])
    return {
        "title": "The condition of the Hoddle Grid's road surfaces: each surface's latest rating, with the age at which it was rated",
        "source": {"id": manifest["id"], "title": manifest["title"], "credit": manifest["credit"], "licence": manifest["licence"],
                   "retrieved": retrieval.get("retrieved_at")},
        "rows": len(rows),
        "rated": len(surfaces),
        "dropped": dict(dropped),
        "without_age": dict(no_age),
        "columns": columns,
        "surfaces": surfaces,
    }


def main(argv: list[str]) -> int:
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args(argv)
    if not (SOURCE / "surfaces.csv").exists():
        print("road-condition: missing the surfaces; run: python data/acquisition/fetch.py melbourne-road-surfaces", file=sys.stderr)
        return 1
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    aged = sum(1 for s in out["surfaces"] if s[9] is not None)
    print(f"road-condition: {out['rated']} rated surfaces of {out['rows']}, {aged} with an age; without one: {out['without_age']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
