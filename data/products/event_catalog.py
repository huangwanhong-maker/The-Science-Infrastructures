"""Build an event product: things that happen at a place and an instant, such as earthquakes.

Each event is a row: origin time, longitude, latitude, depth, magnitude and its
type, a description of the place, and the catalogue's identifier. Times are in
milliseconds since 1970-01-01 UTC; positions in degrees, WGS 84; depth in km.

Output, in data/cache/products/: events-<source>-<profile>.json

Usage::

    python data/products/event_catalog.py                       # every fetched profile
    python data/products/event_catalog.py --profile month
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

DATASET = "usgs-earthquakes"
FIELDS = ["time_ms", "lon", "lat", "depth_km", "mag", "mag_type", "place", "id"]


def iso(ms: float) -> str:
    return dt.datetime.fromtimestamp(ms / 1000, dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build(profile: str) -> dict:
    folder = DATA / "cache" / DATASET / profile
    manifest = json.loads((DATA / "registry" / f"{DATASET}.json").read_text(encoding="utf-8"))
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    collection = json.loads((folder / "events.geojson").read_text(encoding="utf-8"))
    rows = []
    for f in collection["features"]:
        p, (lon, lat, depth) = f["properties"], f["geometry"]["coordinates"][:3]
        if p.get("type") != "earthquake" or p.get("mag") is None:
            continue
        rows.append([p["time"], round(lon, 4), round(lat, 4), round(depth, 2), p["mag"], p.get("magType"),
                     p.get("place"), f["id"]])
    rows.sort(key=lambda r: r[0])
    name = f"usgs-earthquakes-{profile}"
    resolved = retrieval.get("resolved", {})
    product = {
        "schema": "events/1",
        "name": name,
        "title": "Earthquakes, USGS ComCat",
        "kind": "earthquake",
        "domain": manifest.get("domain", "geophysics"),
        "window": {"start": resolved.get("start"), "end": resolved.get("end")},
        "selection": {"min_magnitude": retrieval["parameters"]["min_magnitude"]},
        "fields": FIELDS,
        "units": {"time_ms": "ms since 1970-01-01T00:00:00Z (UTC)", "lon": "degree_east", "lat": "degree_north",
                  "depth_km": "km", "mag": "magnitude (type in mag_type)"},
        "events": rows,
        "source": {"dataset": DATASET, "profile": profile, "provider": manifest["provider"],
                   "credit": manifest.get("credit"),
                   "licence": manifest["licence"], "citation": manifest["citation"],
                   "retrieved_at": retrieval["retrieved_at"], "request": retrieval["requests"][0]["url"]},
    }
    path = OUT / f"events-{name}.json"
    OUT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    mags = [r[4] for r in rows]
    register("events", {
        "name": name, "file": path.name, "title": product["title"], "kind": product["kind"],
        "domain": product["domain"], "count": len(rows), "window": product["window"],
        "first": iso(rows[0][0]) if rows else None, "last": iso(rows[-1][0]) if rows else None,
        "magnitudes": [min(mags), max(mags)] if mags else None,
        "source": {k: product["source"][k] for k in ("dataset", "profile", "provider", "credit", "retrieved_at")},
    }, first="usgs-earthquakes-week")
    return product


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", action="append", help="profile(s) to build (default: every fetched one)")
    args = parser.parse_args(argv)
    root = DATA / "cache" / DATASET
    profiles = args.profile or (sorted(p.name for p in root.iterdir() if (p / "retrieval.json").exists())
                                if root.exists() else [])
    if not profiles:
        print(f"nothing fetched; run: python data/acquisition/fetch.py {DATASET}", file=sys.stderr)
        return 1
    for profile in profiles:
        p = build(profile)
        print(f"events-{p['name']}: {len(p['events'])} earthquakes from {p['window']['start']} to {p['window']['end']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
