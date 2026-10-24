"""Build the trade-flows product: the goods each economy exports to each other, by year, with a place for each economy.

The flows are the OECD's Balanced International Merchandise Trade Statistics (BIMTS, total of
all goods, in millions of US dollars; balanced, and by default adjusted for re-exports). The
places are Natural Earth's label points of the economies (a point inside each country where its
name is drawn, which serves a map of flows better than a centroid that may fall in the sea).
Economies are keyed by ISO 3166-1 alpha-3. Flows to or from codes with no Natural Earth country
(the OECD's aggregates, such as W, the world, and former economies) are left out and listed.

Output, in data/cache/products/: trade-flows.json

Usage::

    python data/products/trade_flows.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

DATASET = "oecd-bimts"
PROFILE = "total-goods"
COUNTRIES = DATA / "cache" / "natural-earth-50m" / "ne_50m_admin_0_countries.geojson"


def economy_code(properties: dict) -> str | None:
    """The ISO alpha-3 code of a Natural Earth country: ISO_A3, else ISO_A3_EH (France and Norway
    carry -99 in ISO_A3 for their overseas parts), else none."""
    for key in ("ISO_A3", "ISO_A3_EH"):
        code = properties.get(key)
        if code and code != "-99" and len(code) == 3:
            return code
    return None


def places(geojson: dict) -> dict[str, dict]:
    """Each economy's name, continent and label point (degrees). Where two features share a code
    (Australia's external territories carry AUS), the one whose own ISO_A3 is the code wins."""
    out: dict[str, dict] = {}
    primary: dict[str, bool] = {}
    for f in geojson["features"]:
        p = f["properties"]
        code = economy_code(p)
        if not code:
            continue
        own = p.get("ISO_A3") == code
        if code in out and (primary[code] or not own):
            continue
        out[code] = {"iso3": code, "name": p.get("NAME_LONG") or p.get("NAME"), "continent": p.get("CONTINENT"),
                     "latitude": round(float(p["LABEL_Y"]), 4), "longitude": round(float(p["LABEL_X"]), 4)}
        primary[code] = own
    return out


def read_flows(rows) -> dict[int, dict[tuple[str, str], float]]:
    """Exports by year: (exporter, importer) -> millions of US dollars."""
    flows: dict[int, dict[tuple[str, str], float]] = {}
    for r in rows:
        if r.get("TRADE_FLOW") != "X" or r.get("PRODUCT_HS") != "_T":
            continue
        value = r.get("OBS_VALUE")
        if value in (None, ""):
            continue
        scale = 10 ** (int(r.get("UNIT_MULT") or 6) - 6)
        year = int(r["TIME_PERIOD"])
        flows.setdefault(year, {})[(r["REF_AREA"], r["COUNTERPART_AREA"])] = float(value) * scale
    return flows


def build() -> dict:
    folder = DATA / "cache" / DATASET / PROFILE
    manifest = json.loads((DATA / "registry" / f"{DATASET}.json").read_text(encoding="utf-8"))
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    with (folder / "bimts.csv").open(encoding="utf-8", newline="") as fh:
        flows = read_flows(csv.DictReader(fh))
    known = places(json.loads(COUNTRIES.read_text(encoding="utf-8")))
    seen = sorted({c for year in flows.values() for pair in year for c in pair})
    economies = [known[c] for c in seen if c in known]
    index = {e["iso3"]: i for i, e in enumerate(economies)}
    unmatched = [c for c in seen if c not in known]
    years = sorted(flows)
    product = {
        "schema": "trade-flows/1",
        "title": "Exports of goods between economies, by year",
        "unit": "millions of US dollars (current, at market exchange rates)",
        "domain": "economics",
        "key_system": "ISO 3166-1 alpha-3",
        "places": "Natural Earth 1:50m label points, degrees of geocentric latitude and longitude",
        "years": years,
        "economies": economies,
        # Kept to the thousand dollars; smaller flows, which round to nothing, are left out.
        "flows": {str(y): [[index[a], index[b], round(v, 3)] for (a, b), v in sorted(flows[y].items())
                           if a in index and b in index and a != b and round(v, 3) > 0] for y in years},
        "unmatched": unmatched,
        "source": {"dataset": DATASET, "profile": PROFILE, "provider": manifest["provider"], "credit": manifest["credit"],
                   "licence": manifest["licence"], "citation": manifest["citation"],
                   "retrieved_at": retrieval["retrieved_at"], "request": retrieval["requests"][0]["url"],
                   "places": "Natural Earth 1:50m admin-0 countries, public domain"},
    }
    path = OUT / "trade-flows.json"
    OUT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    register("trade", {"name": "trade-flows", "file": path.name, "title": product["title"], "unit": product["unit"],
                       "domain": "economics", "years": [years[0], years[-1]], "economies": len(economies),
                       "flows": {str(y): len(product["flows"][str(y)]) for y in years},
                       "source": {k: product["source"][k] for k in ("dataset", "profile", "provider", "credit", "retrieved_at")}})
    return product


def main(argv: list[str]) -> int:
    if not (DATA / "cache" / DATASET / PROFILE / "retrieval.json").exists():
        print(f"nothing fetched; run: python data/acquisition/fetch.py {DATASET}", file=sys.stderr)
        return 1
    if not COUNTRIES.exists():
        print("no Natural Earth countries; run: python data/acquisition/fetch.py natural-earth-50m", file=sys.stderr)
        return 1
    p = build()
    last = str(p["years"][-1])
    print(f"trade-flows: {len(p['economies'])} economies, {p['years'][0]}-{p['years'][-1]}, "
          f"{len(p['flows'][last])} flows in {last}; left out: {', '.join(p['unmatched']) or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
