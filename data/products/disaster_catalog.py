"""Build disaster products from GDACS: events with their alert levels, and the
tracks of tropical cyclones, observed and forecast.

Two products:

- an events product (kind `disaster`): one row per event of the period, a
  tropical cyclone (TC), earthquake (EQ), flood (FL), volcano (VO), drought (DR)
  or wildfire (WF), with its alert level (Green, Orange, Red, as GDACS assesses
  the potential need for international assistance), its severity in words, the
  countries concerned, its dates and a point;
- a tracks product (kind `tracks`): for each tropical cyclone, its positions
  every 6 or 12 hours, observed up to the last observation of the latest
  episode, and forecast beyond it by the responsible centre (JTWC, NHC and
  others, as GDACS relays them).

GDACS notifications for cyclones are automatic and do not replace official
warnings; the products say so.

Output, in data/cache/products/: events-gdacs-<profile>.json, tracks-gdacs-<profile>.json

Usage::

    python data/products/disaster_catalog.py
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

DATASET = "gdacs"
FIELDS = ["time_ms", "end_ms", "lon", "lat", "type", "alert", "name", "severity", "country", "id", "report"]
TYPES = {"TC": "tropical cyclone", "EQ": "earthquake", "FL": "flood", "VO": "volcano", "DR": "drought", "WF": "wildfire"}
CAVEAT = ("GDACS notifications for earthquakes and tropical cyclones are automatic and not reviewed by experts; "
          "they do not replace official warnings from national authorities.")


def ms(text: str) -> int:
    t = dt.datetime.fromisoformat(text.replace("Z", "")).replace(tzinfo=dt.timezone.utc)
    return int(t.timestamp() * 1000)


def track_of(path: Path) -> dict | None:
    """Positions of a cyclone from a GDACS geometry file, with times and whether observed."""
    collection = json.loads(path.read_text(encoding="utf-8"))
    points = []
    meta = None
    for f in collection.get("features", []):
        p = f["properties"]
        if not str(p.get("Class", "")).startswith("Point_Polygon_Point_"):
            continue
        meta = meta or p
        ring = f["geometry"]["coordinates"][0][:-1] or f["geometry"]["coordinates"][0]
        # The point is the centre of the small circle GDACS draws round it.
        lon = sum(x for x, _ in ring) / len(ring)
        lat = sum(y for _, y in ring) / len(ring)
        m = re.match(r"(\d{2})/(\d{2}) (\d{2}):(\d{2})", p.get("polygonlabel", ""))
        if not m:
            continue
        points.append((int(p["Class"].rsplit("_", 1)[1]), m, lon, lat))
    if not meta or not points:
        return None
    points.sort()
    start = dt.datetime.fromisoformat(meta["fromdate"]).replace(tzinfo=dt.timezone.utc)
    last_observed = dt.datetime.fromisoformat(meta["todate"]).replace(tzinfo=dt.timezone.utc)
    year, prev, rows = start.year, None, []
    for _, m, lon, lat in points:
        day, month, hour, minute = (int(g) for g in m.groups())
        t = dt.datetime(year, month, day, hour, minute, tzinfo=dt.timezone.utc)
        if prev and t < prev - dt.timedelta(days=180):          # across the new year
            year += 1
            t = t.replace(year=year)
        prev = t
        rows.append([int(t.timestamp() * 1000), round(lon, 3), round(lat, 3), 1 if t <= last_observed else 0])
    return {"id": f"{meta['eventid']}-{meta['episodeid']}", "name": meta["eventname"], "alert": meta.get("alertlevel"),
            "last_observed_ms": int(last_observed.timestamp() * 1000), "points": rows,
            "report": (meta.get("url") or {}).get("report")}


def build(profile: str) -> tuple[dict, dict]:
    folder = DATA / "cache" / DATASET / profile
    manifest = json.loads((DATA / "registry" / f"{DATASET}.json").read_text(encoding="utf-8"))
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    collection = json.loads((folder / "events.json").read_text(encoding="utf-8"))
    rows = []
    for f in collection.get("features", []):
        p = f["properties"]
        lon, lat = f["geometry"]["coordinates"][:2]
        rows.append([ms(p["fromdate"]), ms(p["todate"]), round(lon, 3), round(lat, 3), p["eventtype"], p.get("alertlevel"),
                     p.get("eventname") or p.get("name"), (p.get("severitydata") or {}).get("severitytext"),
                     p.get("country"), p.get("eventid"), (p.get("url") or {}).get("report")])
    rows.sort(key=lambda r: r[0])
    resolved = retrieval.get("resolved", {})
    source = {"dataset": DATASET, "profile": profile, "provider": manifest["provider"], "credit": manifest.get("credit"),
              "licence": manifest["licence"], "citation": manifest["citation"], "retrieved_at": retrieval["retrieved_at"]}
    name = f"gdacs-{profile}"
    events = {
        "schema": "events/1", "name": name, "title": "Disasters, GDACS", "kind": "disaster",
        "domain": manifest.get("domain", "geophysics"),
        "window": {"start": resolved.get("start"), "end": resolved.get("end")},
        "fields": FIELDS, "types": TYPES, "caveat": CAVEAT,
        "units": {"time_ms": "ms since 1970-01-01T00:00:00Z (UTC), start of the event", "end_ms": "ms, its last update",
                  "lon": "degree_east", "lat": "degree_north"},
        "count": len(rows), "events": rows, "source": source,
    }
    (OUT / f"events-{name}.json").write_text(json.dumps(events, ensure_ascii=False) + "\n", encoding="utf-8")
    register("events", {"name": name, "file": f"events-{name}.json", "title": events["title"], "kind": "disaster",
                        "domain": events["domain"], "count": len(rows), "window": events["window"],
                        "source": {k: source[k] for k in ("dataset", "profile", "credit", "retrieved_at")}})

    tracks = [t for t in (track_of(p) for p in sorted(folder.glob("tc-*.json"))) if t]
    product = {
        "schema": "tracks/1", "name": name, "title": "Tropical cyclones, observed and forecast, GDACS",
        "kind": "cyclone", "domain": "meteorology", "caveat": CAVEAT,
        "fields": ["time_ms", "lon", "lat", "observed"],
        "units": {"time_ms": "ms since 1970-01-01T00:00:00Z (UTC)", "lon": "degree_east", "lat": "degree_north",
                  "observed": "1 up to the last observation of the latest episode, 0 for the official forecast"},
        "tracks": tracks, "source": source,
    }
    (OUT / f"tracks-{name}.json").write_text(json.dumps(product, ensure_ascii=False) + "\n", encoding="utf-8")
    register("tracks", {"name": name, "file": f"tracks-{name}.json", "title": product["title"], "kind": "cyclone",
                        "domain": "meteorology", "count": len(tracks),
                        "source": {k: source[k] for k in ("dataset", "profile", "credit", "retrieved_at")}})
    return events, product


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", default="month")
    args = parser.parse_args(argv)
    if not (DATA / "cache" / DATASET / args.profile / "retrieval.json").exists():
        print(f"not fetched; run: python data/acquisition/fetch.py {DATASET}", file=sys.stderr)
        return 1
    events, tracks = build(args.profile)
    by_type = {}
    for r in events["events"]:
        by_type[r[4]] = by_type.get(r[4], 0) + 1
    print(f"events-{events['name']}: {events['count']} events {by_type}")
    for t in tracks["tracks"]:
        obs = sum(p[3] for p in t["points"])
        print(f"  {t['name']}: {obs} observed and {len(t['points']) - obs} forecast positions, alert {t['alert']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
