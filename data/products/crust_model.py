"""Build the crust products from CRUST1.0 (Laske et al. 2013).

The crust is the Earth's outer rind of light rock, some 7 km thick under the
oceans and 30 to 70 km under the continents, above the denser mantle; the
boundary between them is the Mohorovičić discontinuity (the Moho), where the
speed of P waves jumps from about 7 to about 8 km/s. CRUST1.0 describes it on
1-degree cells as nine layers: water, ice, three of sediment, three of
crystalline crust, and the top of the mantle.

Two products:

- crust-crust1.json and .bin (schema crust/1): every layer's top (km above sea
  level), P and S velocities (km/s) and density (g/cm3), for the column
  beneath a place. Values are stored as 16-bit integers: tops in units of
  10 m, velocities in m/s, densities in kg/m3.
- fields-crust1.json and .bin (a field product): the depth of the Moho below
  sea level, the thickness of the crust (from the solid surface to the Moho),
  and the thickness of the sediments (km), for maps.

Usage::

    python data/acquisition/fetch.py crust1
    python data/products/crust_model.py
"""

from __future__ import annotations

import json
import sys
import tarfile
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(DATA / "formats"))
import netcdf3  # noqa: E402
import field_grid  # noqa: E402
from product_index import OUT, register  # noqa: E402

LAYERS = ["water", "ice", "upper sediments", "middle sediments", "lower sediments",
          "upper crystalline crust", "middle crystalline crust", "lower crystalline crust", "mantle"]
QUANTITIES = [("top", "crust1.bnds", 100.0, "km above sea level"), ("vp", "crust1.vp", 1000.0, "km/s"),
              ("vs", "crust1.vs", 1000.0, "km/s"), ("rho", "crust1.rho", 1000.0, "g/cm3")]


def read(folder: Path) -> dict:
    with tarfile.open(folder / "crust1.0.tar.gz", "r:gz") as tar:
        members = {Path(m.name).name: m for m in tar.getmembers() if m.isfile()}
        out = {}
        for key, name, _, _ in QUANTITIES:
            text = tar.extractfile(members[name]).read().decode("ascii")
            a = np.array(text.split(), dtype=np.float64).reshape(180, 360, 9)
            out[key] = a
    return out


def build() -> tuple[dict, dict]:
    folder = DATA / "cache" / "crust1"
    manifest = json.loads((DATA / "registry" / "crust1.json").read_text(encoding="utf-8"))
    m = read(folder)
    source = {"dataset": "crust1", "provider": manifest["provider"], "credit": manifest["credit"],
              "licence": manifest["licence"], "citation": manifest["citation"]}
    grid = {"lon0": -179.5, "dlon": 1.0, "nlon": 360, "lat0": 89.5, "dlat": 1.0, "nlat": 180, "periodic": True}
    # The column product: per quantity and layer, a block of 180 x 360 int16.
    blocks, offsets = [], {}
    for key, _, scale, _ in QUANTITIES:
        offsets[key] = []
        for k in range(9):
            offsets[key].append(sum(b.nbytes for b in blocks))
            blocks.append(np.round(m[key][:, :, k] * scale).astype("<i2"))
    name = "crust1"
    (OUT / f"crust-{name}.bin").write_bytes(b"".join(b.tobytes() for b in blocks))
    column = {
        "schema": "crust/1", "name": name, "title": "The Earth's crust, CRUST1.0", "domain": "geophysics",
        "file": f"crust-{name}.bin", "grid": grid, "layers": LAYERS,
        "quantities": [{"name": key, "units": units, "scale": 1 / scale, "offsets": offsets[key]} for key, _, scale, units in QUANTITIES],
        "note": "cells of 1 degree, values at their centres; a layer of zero thickness is absent; velocities and densities of absent layers are zero",
        "source": source,
    }
    (OUT / f"crust-{name}.json").write_text(json.dumps(column, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    register("crust", {"name": name, "file": f"crust-{name}.json", "title": column["title"], "domain": "geophysics",
                       "source": {"dataset": "crust1", "credit": manifest["credit"]}})

    # Maps: the depth of the Moho, and the thicknesses of the crust and of its sediments.
    top = m["top"]
    solid = top[:, :, 2]                                   # the top of the upper sediments: the solid surface
    moho = top[:, :, 8]
    lat = 89.5 - np.arange(180, dtype=np.float64)
    lon = -179.5 + np.arange(360, dtype=np.float64)
    variables = {
        "lat": (["lat"], lat, {"units": "degrees_north", "standard_name": "latitude"}),
        "lon": (["lon"], lon, {"units": "degrees_east", "standard_name": "longitude"}),
        "moho_depth": (["lat", "lon"], (-moho).astype(np.float32), {"units": "km", "long_name": "depth of the Moho below sea level"}),
        "crust_thickness": (["lat", "lon"], (solid - moho).astype(np.float32), {"units": "km", "long_name": "thickness of the crust"}),
        "sediment_thickness": (["lat", "lon"], (solid - top[:, :, 5]).astype(np.float32), {"units": "km", "long_name": "thickness of the sediments"}),
    }
    nc = folder / "surfaces.nc"
    netcdf3.write(nc, {"lat": 180, "lon": 360}, variables, {"title": "CRUST1.0 surfaces", "source": manifest["citation"]})
    fields = field_grid.build([(nc, None)], name, "The Earth's crust, CRUST1.0", "geophysics", source, {}, "CRUST1.0")
    return column, fields


def main() -> int:
    if not (DATA / "cache" / "crust1" / "crust1.0.tar.gz").exists():
        print("not fetched; run: python data/acquisition/fetch.py crust1", file=sys.stderr)
        return 1
    column, fields = build()
    names = ", ".join(v["name"] for v in fields["variables"])
    print(f"crust-crust1: 9 layers x 4 quantities on 360 x 180 cells; fields-crust1: {names}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
