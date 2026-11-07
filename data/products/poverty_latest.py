"""The latest survey of poverty and inequality of each economy, from the World Bank's Poverty and Inequality Platform.

From every survey-year estimate fetched by the manifest `worldbank-pip`, each economy keeps the latest national one, keyed by its
ISO 3166-1 alpha-3 code so that it joins the regions product. Where an economy has only urban and rural estimates in its latest
year (three do in their early years), it is skipped for that year rather than averaged. Names are the data library's policy's
(naming.py). Values are as the Platform gives them: the headcount a share at the poverty line asked for, the Gini an index on [0, 1].

    python data/products/poverty_latest.py

writes data/cache/products/poverty-latest.json.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import naming  # noqa: E402

CACHE = HERE.parent / "cache"


def build(path: Path = CACHE / "worldbank-pip" / "international" / "pip.json") -> dict:
    rows = json.loads(path.read_text(encoding="utf-8"))
    names = naming.Naming.load()
    latest: dict[str, dict] = {}
    for r in rows:
        if r.get("reporting_level") != "national" or r.get("headcount") is None:
            continue
        code, year = r["country_code"], r["reporting_year"]
        if code not in latest or year > latest[code]["year"]:
            latest[code] = {"name": names.name(code, r["country_name"]), "year": year, "headcount": r["headcount"], "gini": r.get("gini"), "mean": r.get("mean"), "median": r.get("median"),
                            "welfare": r.get("welfare_type"), "survey": r.get("survey_acronym"), "people": r.get("reporting_pop")}
    line = rows[0]["poverty_line"] if rows else None
    return {"schema": "poverty-latest/1", "poverty_line": line, "unit": "2021 PPP dollars a day",
            "source": {"dataset": "worldbank-pip", "credit": "World Bank, Poverty and Inequality Platform, CC BY 4.0", "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},
            "economies": latest}


def main() -> int:
    out = CACHE / "products" / "poverty-latest.json"
    product = build()
    out.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print(f"{len(product['economies'])} economies, to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
