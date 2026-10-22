"""Build tracks products from IBTrACS: the best tracks of tropical cyclones.

IBTrACS (Knapp et al. 2010; version 4r01, Gahtan et al. 2024) merges the
tracks of every warning centre into one record per storm. For each storm the
product keeps its positions and times and, at each, its central pressure and
maximum sustained wind, each taken first from the agency WMO makes responsible
for the basin (WMO_PRES, WMO_WIND), then, where it gives none, from the first
of the other agencies that does, in this order: the US agencies (NHC, CPHC,
JTWC), Tokyo, CMA, Hong Kong, KMA, New Delhi, La Réunion, the Australian
centres, Nadi, Wellington. An agency is passed over only outside the span of
times at which it reports the storm (intensities()).

The winds are not homogeneous: each agency's maximum sustained wind is a mean
over its own period (1 minute for the US agencies, 2 for CMA, 3 for New Delhi,
10 for the others), and IBTrACS makes no adjustment. Each position therefore
names the agency its wind is from, and the product lists the averaging period
of each. To compare them, Harper et al. (2010, WMO/TD-No. 1555) recommend a
factor of 0.93 from a 1-minute to a 10-minute wind at sea (the older 0.88 is
still met); no factor is applied here.

The positions are IBTrACS's combined ones (LAT, LON). IBTrACS also fills
3-hourly times by interpolation; those whose every agency's position is
interpolated (IFLAG: no O, I or V) are left out, so the agencies' own fixes
remain, usually every 6 hours, with extra ones at landfall. Spur tracks
(alternative positions early in a storm's life, TRACK_TYPE *spur) are left out.
Storms of the last weeks are provisional (TRACK_TYPE PROVISIONAL), from
real-time sources such as TC Vitals, until the agencies' reanalysis.

The product is a tracks product (schema tracks/1, as written by
disaster_catalog.py): each track has its positions as rows of `fields`, here
time_ms, lon, lat, observed (always 1: these are analyses, not forecasts),
pressure_hpa, wind_ms, nature (IBTrACS's storm type: TS tropical, ET
extratropical, SS subtropical, DS disturbance, MX mixed, NR not reported),
pressure_agency and wind_agency. Longitudes are reduced to [-180, 180).

Output, in data/cache/products/: tracks-ibtracs-<profile>.json, for the profiles
`active` (storms of the last 7 days) and `last3years` (for verifying forecasts).

Usage::

    python data/acquisition/fetch.py ibtracs --refresh                 # active storms
    python data/acquisition/fetch.py ibtracs --profile last3years
    python data/products/ibtracs_tracks.py                              # every fetched profile
    python data/products/ibtracs_tracks.py --profile active
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

DATASET = "ibtracs"
KNOT = 1852.0 / 3600.0            # m s-1, exactly (BIPM SI Brochure, 2019, table 8)
FIELDS = ["time_ms", "lon", "lat", "observed", "pressure_hpa", "wind_ms", "nature", "pressure_agency", "wind_agency"]

# The agencies after the WMO one, in the order of preference, by their column prefix.
AGENCIES = ["USA", "TOKYO", "CMA", "HKO", "KMA", "NEWDELHI", "REUNION", "BOM", "NADI", "WELLINGTON"]

# Averaging period of each agency's maximum sustained wind, in minutes (IBTrACS v04
# column documentation; the WMO_AGENCY names and the agencies' prefixes in lower case).
AVERAGING_MINUTES = {
    "hurdat_atl": 1, "hurdat_epa": 1, "atcf": 1, "cphc": 1, "usa": 1, "tcvitals": 1, "nhc_working_bt": 1,
    "jtwc_wp": 1, "jtwc_io": 1, "jtwc_sh": 1, "jtwc_ep": 1, "jtwc_cp": 1,
    "tokyo": 10, "cma": 2, "hko": 10, "kma": 10, "newdelhi": 3, "reunion": 10, "bom": 10, "nadi": 10, "wellington": 10,
}
CAVEAT = ("Best tracks of recent storms are provisional until the warning centres' reanalysis; maximum winds are "
          "means over each agency's own period (see agencies), not adjusted to one another.")


def value(row: dict, name: str, lo: float, hi: float) -> float | None:
    """A number of the row, or None where blank or outside [lo, hi] (as CMA's 0 for an unknown wind)."""
    text = (row.get(name) or "").strip()
    if not text:
        return None
    try:
        x = float(text)
    except ValueError:
        return None
    return x if lo <= x <= hi else None


# What each quantity is called in the columns, and the range outside which a value
# is missing: pressures outside 850-1050 hPa and winds of 0 or less (the US
# agencies' stated ranges, IBTrACS documentation; CMA codes an unknown wind as 0).
QUANTITIES = {"pressure": ("PRES", 850.0, 1050.0), "wind": ("WIND", 1e-6, 300.0)}


def sources(row: dict) -> list[tuple[str, str]]:
    """(column prefix, agency's name) in the order of preference, for one row."""
    usa = (row.get("USA_AGENCY") or "").strip().lower() or "usa"
    return [("WMO", (row.get("WMO_AGENCY") or "").strip().lower() or "wmo")] + \
        [(a, usa if a == "USA" else a.lower()) for a in AGENCIES]


def intensities(rows: list[tuple[int, dict]]) -> list[dict]:
    """The pressure (hPa) and wind (m s-1) of one storm at each of its times, each
    with its agency.

    At each time a quantity comes from the first agency in the order of
    preference that gives it. An agency that gives it at earlier and later times
    of the storm but not at this one ends the search there, and the value is
    missing: otherwise a 10-minute wind at the agency's own times would alternate
    with another agency's 1-minute wind at the times between them.
    """
    spans = {}
    for q, (suffix, lo, hi) in QUANTITIES.items():
        for t, row in rows:
            for prefix, _ in sources(row):
                if value(row, f"{prefix}_{suffix}", lo, hi) is not None:
                    a, b = spans.get((q, prefix), (t, t))
                    spans[(q, prefix)] = (min(a, t), max(b, t))
    out = []
    for t, row in rows:
        point = {}
        for q, (suffix, lo, hi) in QUANTITIES.items():
            point[q] = point[f"{q}_agency"] = None
            for prefix, name in sources(row):
                x = value(row, f"{prefix}_{suffix}", lo, hi)
                if x is not None:
                    point[q], point[f"{q}_agency"] = x, name
                    break
                a, b = spans.get((q, prefix), (None, None))
                if a is not None and a <= t <= b:
                    break
        out.append(point)
    return out


def read(path: Path) -> list[dict]:
    """Tracks from an IBTrACS CSV file (two header lines: names, then units)."""
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.reader(handle)
        names = next(rows)
        next(rows)                                   # the units
        storms = {}
        for r in rows:
            row = dict(zip(names, r))
            if "spur" in row.get("TRACK_TYPE", "").lower():
                continue
            if not any(c in "OIV" for c in row.get("IFLAG", "O")):
                continue
            if value(row, "LAT", -90.0, 90.0) is None or value(row, "LON", -360.0, 360.0) is None:
                continue
            t = dt.datetime.strptime(row["ISO_TIME"].strip(), "%Y-%m-%d %H:%M:%S").replace(tzinfo=dt.timezone.utc)
            storms.setdefault(row["SID"].strip(), {})[int(t.timestamp() * 1000)] = row
    tracks = []
    for sid, by_time in storms.items():
        rows = sorted(by_time.items())
        points = []
        for (t, row), q in zip(rows, intensities(rows)):
            lon, lat = value(row, "LON", -360.0, 360.0), value(row, "LAT", -90.0, 90.0)
            points.append([t, round((lon + 180.0) % 360.0 - 180.0, 2), round(lat, 2), 1,
                           None if q["pressure"] is None else round(q["pressure"], 1),
                           None if q["wind"] is None else round(q["wind"] * KNOT, 2),
                           row.get("NATURE", "").strip() or None, q["pressure_agency"], q["wind_agency"]])
        row = rows[0][1]
        name = row.get("NAME", "").strip()
        if not name or name in ("NOT_NAMED", "UNNAMED"):
            name = row.get("USA_ATCF_ID", "").strip() or sid
        winds = [p[5] for p in points if p[5] is not None]
        pressures = [p[4] for p in points if p[4] is not None]
        wmo = [(r.get("WMO_AGENCY") or "").strip() for _, r in rows]
        wmo = [a for a in wmo if a]
        tracks.append({
            "id": sid, "name": name, "basin": row.get("BASIN", "").strip(), "season": int(row["SEASON"]),
            "atcf_id": row.get("USA_ATCF_ID", "").strip() or None,
            "provisional": any("PROVISIONAL" in r.get("TRACK_TYPE", "").upper() for _, r in rows),
            "wmo_agency": max(set(wmo), key=wmo.count) if wmo else None,
            "max_wind_ms": max(winds) if winds else None, "min_pressure_hpa": min(pressures) if pressures else None,
            "points": points,
        })
    tracks.sort(key=lambda tr: tr["points"][0][0])
    return tracks


def build(profile: str) -> dict:
    folder = DATA / "cache" / DATASET / profile
    manifest = json.loads((DATA / "registry" / f"{DATASET}.json").read_text(encoding="utf-8"))
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    tracks = read(folder / "ibtracs.csv")
    times = [p[0] for tr in tracks for p in tr["points"]]
    iso = lambda ms: dt.datetime.fromtimestamp(ms / 1000, dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")  # noqa: E731
    window = {"start": iso(min(times)), "end": iso(max(times))} if times else {"start": None, "end": None}
    source = {"dataset": DATASET, "profile": profile, "provider": manifest["provider"], "credit": manifest.get("credit"),
              "licence": manifest["licence"], "citation": manifest["citation"], "retrieved_at": retrieval["retrieved_at"],
              "subset": retrieval["parameters"].get("subset")}
    used = sorted({a for tr in tracks for p in tr["points"] for a in p[7:9] if a})
    name = f"ibtracs-{profile}"
    product = {
        "schema": "tracks/1", "name": name, "title": "Tropical cyclones, best tracks, IBTrACS",
        "kind": "cyclone", "domain": "meteorology", "caveat": CAVEAT, "window": window,
        "fields": FIELDS,
        "units": {"time_ms": "ms since 1970-01-01T00:00:00Z (UTC)", "lon": "degree_east", "lat": "degree_north",
                  "observed": "1: a best-track (or provisional) analysis, not a forecast",
                  "pressure_hpa": "hPa, minimum central pressure", "wind_ms": "m s-1, maximum sustained wind, "
                  "a mean over the period of the agency named in `wind_agency`",
                  "nature": "TS tropical, ET extratropical, SS subtropical, DS disturbance, MX mixed, NR not reported",
                  "pressure_agency": "the agency the pressure is from",
                  "wind_agency": "the agency the wind is from, whose averaging period `agencies` gives"},
        "agencies": {a: {"averaging_minutes": AVERAGING_MINUTES.get(a)} for a in used},
        "tracks": tracks, "source": source,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"tracks-{name}.json"
    path.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    register("tracks", {"name": name, "file": path.name, "title": product["title"], "kind": "cyclone",
                        "domain": "meteorology", "count": len(tracks), "window": window,
                        "source": {k: source[k] for k in ("dataset", "profile", "credit", "retrieved_at")}})
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
        product = build(profile)
        tracks = product["tracks"]
        print(f"tracks-{product['name']}: {len(tracks)} storms, {sum(len(t['points']) for t in tracks)} positions, "
              f"{product['window']['start']} to {product['window']['end']}")
        for t in tracks[-8:]:
            wind = f"{t['max_wind_ms']:.0f} m/s" if t["max_wind_ms"] else "no wind"
            pres = f"{t['min_pressure_hpa']:.0f} hPa" if t["min_pressure_hpa"] else "no pressure"
            print(f"  {t['name']} ({t['basin']}, {t['season']}): {len(t['points'])} positions, {wind}, {pres}"
                  + (", provisional" if t["provisional"] else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
