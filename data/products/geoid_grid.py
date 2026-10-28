"""The EGM96 geoid grid as a product: the PGM taken out of GeographicLib's archive.

The dataset egm96-geoid is a bzip2 tar archive; this builder extracts its grid,
geoids/egm96-15.pgm (1440 x 721 big-endian 16-bit values at 15 arc-minutes, N = -108 + 0.003 x
value metres), to data/cache/products/geoid/egm96-15.pgm, where the tile service and the geodesy
domain's C++ checks read it, and records where it came from beside it.

    python data/products/geoid_grid.py
"""

from __future__ import annotations

import hashlib
import json
import sys
import tarfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "cache"
SOURCE = CACHE / "egm96-geoid" / "egm96-15.tar.bz2"
OUT = CACHE / "products" / "geoid"


def main() -> int:
    if not SOURCE.exists():
        print(f"{SOURCE} is missing: python data/acquisition/fetch.py egm96-geoid", file=sys.stderr)
        return 1
    with tarfile.open(SOURCE, "r:bz2") as tar:
        member = tar.getmember("geoids/egm96-15.pgm")
        data = tar.extractfile(member).read()
    if not data.startswith(b"P5"):
        print("the archive's grid is not a binary PGM", file=sys.stderr)
        return 1
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "egm96-15.pgm").write_bytes(data)
    (OUT / "egm96-15.json").write_text(json.dumps({
        "schema": "geoid-grid/1",
        "name": "egm96-15",
        "title": "EGM96 geoid undulation, 15 arc-minute grid",
        "file": "egm96-15.pgm",
        "format": "binary PGM (P5), 1440 x 721 big-endian uint16, point-registered from 90 N and 0 E at 0.25 degrees; N = Offset + Scale x value (Offset -108, Scale 0.003, in its comments), metres above the WGS 84 ellipsoid",
        "source": {"dataset": "egm96-geoid", "file": SOURCE.name, "sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest()},
        "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }, indent=1), encoding="utf-8")
    print(f"wrote {OUT / 'egm96-15.pgm'} ({len(data)} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
