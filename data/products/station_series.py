"""Build a series product: the daily positions of one GPS station, east, north
and up, from a Nevada Geodetic Laboratory tenv3 file.

Positions are given in the file as "graticule distance coordinates" (Blewitt
2024): easting, northing and height in metres, each split into an integer and
a fractional part. Here each is the sum of its parts less its value on the
first day, in millimetres, so that the series shows how far the station has
moved since it began, with steps where it jumped (an earthquake, a new
antenna) and seasons where the ground breathes with water and heat.

Output: data/cache/products/series/series-ngl-<station>.json (schema series/1),
made on demand (data/acquisition/on_demand.py), not listed in the index.

Usage::

    python data/acquisition/fetch.py ngl-station --set station=P395
    python data/products/station_series.py --profile custom-<hash>
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT  # noqa: E402

SERIES = OUT / "series"
MONTHS = {m: i + 1 for i, m in enumerate(["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"])}


def build(folder: Path, manifest: dict) -> Path:
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    station = retrieval["parameters"]["station"]
    days, east, north, up, se, sn, su = [], [], [], [], [], [], []
    place = None
    first = None
    for line in (folder / "station.tenv3").read_text(encoding="latin-1").splitlines():
        f = line.split()
        if len(f) < 23 or f[0] == "site":
            continue
        y, mo, d = int(f[1][:2]), MONTHS[f[1][2:5]], int(f[1][5:7])
        day = dt.date(2000 + y if y < 80 else 1900 + y, mo, d)
        e, n, u = float(f[7]) + float(f[8]), float(f[9]) + float(f[10]), float(f[11]) + float(f[12])
        if first is None:
            first = (e, n, u)
            place = {"latitude": float(f[20]), "longitude": ((float(f[21]) + 180) % 360) - 180, "height_m": float(f[22])}
        days.append(day.isoformat())
        east.append(round((e - first[0]) * 1000, 2)); north.append(round((n - first[1]) * 1000, 2)); up.append(round((u - first[2]) * 1000, 2))
        se.append(round(float(f[14]) * 1000, 2)); sn.append(round(float(f[15]) * 1000, 2)); su.append(round(float(f[16]) * 1000, 2))
    if not days:
        raise SystemExit(f"{station}: no positions in the file")
    name = f"series-ngl-{station.lower()}"
    product = {
        "schema": "series/1", "name": name, "title": f"GPS station {station}: daily positions", "kind": "gnss",
        "domain": "geophysics", "station": station, "place": place, "frame": "IGS20 (ITRF2020)",
        "days": days,
        "components": [
            {"name": "east", "units": "mm", "values": east, "sigma": se},
            {"name": "north", "units": "mm", "values": north, "sigma": sn},
            {"name": "up", "units": "mm", "values": up, "sigma": su},
        ],
        "note": "displacement from the first day; days in UTC",
        "source": {"dataset": manifest["id"], "profile": retrieval["profile"], "provider": manifest["provider"],
                   "credit": manifest["credit"], "licence": manifest["licence"], "citation": manifest["citation"],
                   "retrieved_at": retrieval["retrieved_at"], "parameters": retrieval["parameters"]},
    }
    SERIES.mkdir(parents=True, exist_ok=True)
    path = SERIES / f"{name}.json"
    path.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    return path


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", default="p123", help="the retrieval's cache folder")
    args = parser.parse_args(argv)
    folder = DATA / "cache" / "ngl-station" / args.profile
    if not (folder / "retrieval.json").exists():
        print(f"not fetched: {folder}", file=sys.stderr)
        return 1
    manifest = json.loads((DATA / "registry" / "ngl-station.json").read_text(encoding="utf-8"))
    path = build(folder, manifest)
    p = json.loads(path.read_text(encoding="utf-8"))
    print(f"{path.relative_to(DATA.parent)}: {len(p['days'])} days, {p['days'][0]} to {p['days'][-1]}; "
          f"moved {p['components'][0]['values'][-1]:.0f} mm east, {p['components'][1]['values'][-1]:.0f} mm north")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
