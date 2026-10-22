"""Build the stations product: GPS station velocities from the Nevada Geodetic
Laboratory's MIDAS estimates (Blewitt et al. 2016).

MIDAS takes the median of the slopes between pairs of positions a year apart,
which makes it robust to steps (earthquakes, changes of antenna) and to
seasonal motion without detecting either. Its uncertainties are realistic for
series of a few years or more.

Kept: stations with at least two years of data and horizontal uncertainties
of at most 2 mm/yr. Velocities in mm/yr, east, north and up, in the frame of
the file fetched (IGS20 by default, aligned with ITRF2020); longitudes are
brought from NGL's -360..0 to -180..180.

Output, in data/cache/products/: stations-ngl-midas.json (schema stations/1),
one row per station: [id, lon, lat, height_m, ve, vn, vu, se, sn, su, first, last]
with first and last epochs in decimal years.

Usage::

    python data/acquisition/fetch.py ngl-midas
    python data/products/station_velocities.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

FIELDS = ["id", "lon", "lat", "height_m", "ve", "vn", "vu", "se", "sn", "su", "first", "last"]
MIN_YEARS, MAX_SIGMA = 2.0, 2.0


def build(profile: str) -> dict:
    folder = DATA / "cache" / "ngl-midas" / profile
    manifest = json.loads((DATA / "registry" / "ngl-midas.json").read_text(encoding="utf-8"))
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    rows, seen, dropped = [], 0, 0
    for line in (folder / "midas.txt").read_text(encoding="latin-1").splitlines():
        f = line.split()
        if len(f) < 27:
            continue
        seen += 1
        first, last, years = float(f[2]), float(f[3]), float(f[4])
        ve, vn, vu, se, sn, su = (float(x) * 1000 for x in f[8:14])
        lat, lon, h = float(f[24]), float(f[25]), float(f[26])
        if years < MIN_YEARS or se > MAX_SIGMA or sn > MAX_SIGMA:
            dropped += 1
            continue
        lon = ((lon + 180.0) % 360.0) - 180.0
        rows.append([f[0], round(lon, 5), round(lat, 5), round(h, 1), round(ve, 2), round(vn, 2), round(vu, 2),
                     round(se, 2), round(sn, 2), round(su, 2), round(first, 3), round(last, 3)])
    frame = retrieval["parameters"]["frame"]
    name = "ngl-midas" if frame == "IGS" else f"ngl-midas-{frame.lower()}"
    source = {"dataset": "ngl-midas", "profile": profile, "provider": manifest["provider"], "credit": manifest["credit"],
              "licence": manifest["licence"], "citation": manifest["citation"], "retrieved_at": retrieval["retrieved_at"]}
    product = {
        "schema": "stations/1", "name": name, "title": "GPS station velocities (MIDAS)", "kind": "gnss",
        "domain": "geophysics", "frame": "IGS20 (ITRF2020)" if frame == "IGS" else f"{frame} plate fixed (NGL)",
        "selection": {"min_years": MIN_YEARS, "max_horizontal_sigma_mm_per_yr": MAX_SIGMA, "kept": len(rows), "of": seen},
        "fields": FIELDS,
        "units": {"lon": "degree_east", "lat": "degree_north", "height_m": "m (ellipsoidal)", "ve": "mm/yr", "vn": "mm/yr", "vu": "mm/yr",
                  "se": "mm/yr", "sn": "mm/yr", "su": "mm/yr", "first": "decimal year", "last": "decimal year"},
        "count": len(rows), "stations": rows, "source": source,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"stations-{name}.json").write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    register("stations", {"name": name, "file": f"stations-{name}.json", "title": product["title"], "kind": "gnss",
                          "domain": "geophysics", "count": len(rows), "frame": product["frame"],
                          "source": {k: source[k] for k in ("dataset", "profile", "credit", "retrieved_at")}})
    return product


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", default="igs20")
    args = parser.parse_args(argv)
    if not (DATA / "cache" / "ngl-midas" / args.profile / "retrieval.json").exists():
        print("not fetched; run: python data/acquisition/fetch.py ngl-midas", file=sys.stderr)
        return 1
    p = build(args.profile)
    print(f"stations-{p['name']}: {p['count']} of {p['selection']['of']} stations kept ({p['frame']})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
