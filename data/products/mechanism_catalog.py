"""Build the mechanisms product: the moment tensors of the Global CMT catalogue.

A moment tensor M describes the forces of an earthquake's source as a symmetric
3 x 3 tensor. For a slip on a fault it is a double couple, M0 (n d^T + d n^T),
with n the fault's normal and d the direction of slip; its eigenvectors are the
T (tension), N (null) and P (pressure) axes, and the two planes normal to
T + P and T - P are the fault and the auxiliary plane, which waves alone cannot
tell apart. The atlas draws each mechanism as a "beach ball" and computes its
axes and planes from the six elements kept here; the catalogue's own planes
are kept beside them to check that computation.

Kept for each earthquake: the centroid's time, place and depth; the moment
magnitude Mw = (2/3)(log10 M0 - 16.1), M0 in dyne-cm (Hanks and Kanamori 1979);
the six elements in (r up, t south, p east), divided by M0; and the strike, dip
and rake of the catalogue's two nodal planes.

Output, in data/cache/products/: events-gcmt.json (schema events/1, kind
"mechanism"), registered under "events".

Usage::

    python data/acquisition/fetch.py gcmt
    python data/products/mechanism_catalog.py
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import math
import re
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

FIELDS = ["time_ms", "lon", "lat", "depth_km", "mw", "mrr", "mtt", "mpp", "mrt", "mrp", "mtp",
          "strike1", "dip1", "rake1", "strike2", "dip2", "rake2", "name", "region", "quick"]
NUMBER = re.compile(r"-?\d+(?:\.\d+)?")


def records(text: str):
    lines = [l for l in text.splitlines() if l.strip()]
    for k in range(0, len(lines) - 4, 5):
        yield lines[k:k + 5]


def parse(rec: list[str]) -> list | None:
    l1, l2, l3, l4, l5 = (l.ljust(80) for l in rec)
    try:
        y, mo, d = (int(x) for x in l1[5:15].split("/"))
        hh, mm, ss = l1[16:26].split(":")
        t = dt.datetime(y, mo, d, tzinfo=dt.timezone.utc) + dt.timedelta(hours=int(hh), minutes=int(mm), seconds=float(ss))
        centroid = [float(x) for x in NUMBER.findall(l3[9:58])]
        shift, lat, lon, depth = centroid[0], centroid[2], centroid[4], centroid[6]
        exponent = int(l4[0:2])
        tensor = [float(x) for x in NUMBER.findall(l4[2:80])][0::2]
        scalar = float(l5[49:56])
        planes = [int(x) for x in NUMBER.findall(l5[56:80])]
    except (ValueError, IndexError):
        return None
    if len(tensor) != 6 or len(planes) != 6 or scalar <= 0:
        return None
    m0 = scalar * 10 ** exponent                                      # dyne-cm
    mw = (2 / 3) * (math.log10(m0) - 16.1)
    t += dt.timedelta(seconds=shift)
    unit = [round(x / scalar, 4) for x in tensor]                      # elements over M0
    return [int(t.timestamp() * 1000), round(lon, 2), round(lat, 2), round(depth, 1), round(mw, 2), *unit,
            *planes, l2[0:16].strip(), l1[56:80].strip().title(), 1 if "Q-" in l3[59:80] else 0]


def build(profile: str) -> dict:
    folder = DATA / "cache" / "gcmt" / profile
    manifest = json.loads((DATA / "registry" / "gcmt.json").read_text(encoding="utf-8"))
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    final = gzip.decompress((folder / "final.ndk.gz").read_bytes()).decode("latin-1")
    quick = (folder / "quick.ndk").read_text(encoding="latin-1")
    rows, names, bad = [], set(), 0
    for text in (final, quick):
        for rec in records(text):
            row = parse(rec)
            if row is None:
                bad += 1
                continue
            if row[FIELDS.index("name")] in names:           # a final solution supersedes a quick one
                continue
            names.add(row[FIELDS.index("name")])
            rows.append(row)
    rows.sort(key=lambda r: r[0])
    source = {"dataset": "gcmt", "profile": profile, "provider": manifest["provider"], "credit": manifest["credit"],
              "licence": manifest["licence"], "citation": manifest["citation"], "retrieved_at": retrieval["retrieved_at"]}
    start = dt.datetime.fromtimestamp(rows[0][0] / 1000, dt.timezone.utc).date().isoformat()
    end = dt.datetime.fromtimestamp(rows[-1][0] / 1000, dt.timezone.utc).date().isoformat()
    product = {
        "schema": "events/1", "name": "gcmt", "title": "Earthquake mechanisms, Global CMT", "kind": "mechanism",
        "domain": "geophysics", "window": {"start": start, "end": end},
        "fields": FIELDS,
        "units": {"time_ms": "ms since 1970-01-01T00:00:00Z (UTC), centroid time", "lon": "degree_east", "lat": "degree_north",
                  "depth_km": "km, centroid", "mw": "moment magnitude",
                  "mrr": "moment-tensor element over the scalar moment; r up, t south, p east",
                  "strike1": "degrees, nodal plane 1 of the catalogue (Aki and Richards convention)"},
        "unparsed": bad, "count": len(rows), "events": rows, "source": source,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "events-gcmt.json").write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    register("events", {"name": "gcmt", "file": "events-gcmt.json", "title": product["title"], "kind": "mechanism",
                        "domain": "geophysics", "count": len(rows), "window": product["window"],
                        "source": {k: source[k] for k in ("dataset", "profile", "credit", "retrieved_at")}})
    return product


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", default="catalogue")
    args = parser.parse_args(argv)
    if not (DATA / "cache" / "gcmt" / args.profile / "retrieval.json").exists():
        print("not fetched; run: python data/acquisition/fetch.py gcmt", file=sys.stderr)
        return 1
    p = build(args.profile)
    print(f"events-gcmt: {p['count']} mechanisms, {p['window']['start']} to {p['window']['end']}; {p['unparsed']} records not read")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
