"""Build an indicator product: one value per economy and year, from the World Bank.

The product is keyed by the World Bank's three-letter economy codes, which are
ISO 3166-1 alpha-3 codes except for a few, such as XKX for Kosovo; each economy
also carries its ISO alpha-2 code. It holds no geometry: the application joins
it to the countries of a geography product by these keys. Regional and income
aggregates are kept apart from economies.

Output, in data/cache/products/: indicator-<profile>.json

Usage::

    python data/products/country_indicators.py                       # every fetched profile
    python data/products/country_indicators.py --profile life-expectancy
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

DATASET = "worldbank-indicators"


def build(profile: str) -> dict:
    folder = DATA / "cache" / DATASET / profile
    manifest = json.loads((DATA / "registry" / f"{DATASET}.json").read_text(encoding="utf-8"))
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    meta, rows = json.loads((folder / "values.json").read_text(encoding="utf-8"))
    if meta.get("pages", 1) != 1:
        raise SystemExit(f"{profile}: the response has {meta['pages']} pages; raise per_page in the manifest")
    economies_meta, economies_rows = json.loads((folder / "economies.json").read_text(encoding="utf-8"))
    if economies_meta.get("pages", 1) != 1:
        raise SystemExit("the list of economies has more than one page")
    aggregate = {e["id"]: e["region"]["value"] == "Aggregates" for e in economies_rows}
    iso2 = {e["id"]: e["iso2Code"] for e in economies_rows}
    from_iso2 = {e["iso2Code"]: e["id"] for e in economies_rows}

    years = sorted({int(r["date"]) for r in rows})
    first = years[0]
    indicator = rows[0]["indicator"]
    economies, aggregates = {}, {}
    for r in rows:
        # A few aggregates (income groups) come with no three-letter code, only the
        # two-letter one, which the list of economies translates.
        code = r["countryiso3code"] or from_iso2.get(r["country"]["id"], r["country"]["id"])
        if not code:
            continue
        target = aggregates if aggregate.get(code, False) else economies
        entry = target.setdefault(code, {"name": r["country"]["value"], "iso2": iso2.get(code),
                                         "values": [None] * len(years)})
        if r["value"] is not None:
            entry["values"][int(r["date"]) - first] = r["value"]
    unknown = [c for c in economies if c not in aggregate]
    if unknown:
        print(f"  note: codes not in the list of economies, kept as economies: {', '.join(unknown)}")

    name_text = indicator["value"]
    m = re.search(r"\(([^()]*)\)\s*$", name_text)
    coverage = {str(y): sum(e["values"][i] is not None for e in economies.values()) for i, y in enumerate(years)}
    product = {
        "schema": "indicator/1",
        "name": profile,
        "code": indicator["id"],
        "title": name_text,
        "quantity": name_text[:m.start()].strip().rstrip(",") if m else name_text,
        "unit": m.group(1) if m else "",
        "domain": manifest.get("domain", "geography"),
        "key_system": "World Bank economy codes (ISO 3166-1 alpha-3 with exceptions such as XKX); iso2 is ISO 3166-1 alpha-2",
        "time": "calendar year",
        "years": years,
        "economies": economies,
        "aggregates": aggregates,
        "source": {"dataset": DATASET, "profile": profile, "provider": manifest["provider"],
                   "credit": manifest.get("credit"),
                   "licence": manifest["licence"], "citation": manifest["citation"],
                   "retrieved_at": retrieval["retrieved_at"], "last_updated": meta.get("lastupdated"),
                   "request": retrieval["requests"][0]["url"]},
    }
    path = OUT / f"indicator-{profile}.json"
    OUT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    register("indicators", {
        "name": profile, "file": path.name, "title": name_text, "code": indicator["id"],
        "quantity": product["quantity"], "unit": product["unit"], "domain": product["domain"],
        "years": [years[0], years[-1]], "economies": len(economies), "coverage": coverage,
        "source": {k: product["source"][k] for k in ("dataset", "profile", "provider", "credit", "retrieved_at")},
    }, first="population")
    return product


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", action="append", help="profile(s) to build (default: every fetched one)")
    args = parser.parse_args(argv)
    root = DATA / "cache" / DATASET
    profiles = args.profile or (sorted(p.name for p in root.iterdir() if (p / "retrieval.json").exists())
                                if root.exists() else [])
    if not profiles:
        print(f"nothing fetched; run: python data/acquisition/fetch.py {DATASET} --profile population",
              file=sys.stderr)
        return 1
    for profile in profiles:
        p = build(profile)
        print(f"indicator-{profile}: {p['code']}, {len(p['economies'])} economies, "
              f"{p['years'][0]}-{p['years'][-1]}, unit {p['unit']!r}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
