"""A field as a file for other software: the GFS temperature at 2 m over Japan at
the stored time nearest now, as a GeoTIFF (float32, EPSG:4326), which QGIS and
GDAL open with its coordinates, units and no-data value; and the same box at
every time as netCDF.

    python api/examples/python/export_geotiff.py [out-directory] [--url ...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
import science  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("out", nargs="?", default=".", help="the directory to write into")
parser.add_argument("--url")
args = parser.parse_args()

api = science.Client(args.url)
gfs = api.field("gfs-forecast")
japan = (122, 24, 150, 46)                       # west, south, east, north, degrees
tif = gfs.export("tmp2m", time="now", format="tif", bbox=japan, path=Path(args.out))
nc = gfs.export("tmp2m", time="all", format="nc", bbox=japan, path=Path(args.out))
print(f"{tif} ({tif.stat().st_size:,} bytes)\n{nc} ({nc.stat().st_size:,} bytes)")
