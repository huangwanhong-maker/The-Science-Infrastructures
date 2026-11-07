"""Freetown's flood risk, made ready to draw: the depth maps at each return period, the modelled buildings, the wards with their annual
losses, the places the City Council records as flooding, and the loss tables, from the Freetown Flood Risk Assessment (World Bank, JBA
Risk Management, Freetown City Council; CC BY 4.0).

The maps are of two perils, river floods (FLRF) and surface water (FLSW), at return periods of 20 to 1500 years, on 1 m grids in
WGS 84, published with overviews. The first overview but one, about 4 m, is taken: finer than anything drawn at a city's scale. An
overview's file holds no georeference, so that of its map is read from the map's own TIFF, whose directory stands at its end (the whole
member is read, once, for each peril). Depths are reduced to a byte: depth / 0.05 m, to 12.75 m; 0 is dry or no data.

    python data/products/freetown_flood.py

writes data/cache/products/freetown-flood.json, with its files under data/cache/products/freetown-flood/:
    <PERIL>-rp<N>.bin    a depth map, one byte a cell, rows from the south
    buildings.bin        the 153,436 buildings: float64 lat, lon, value (USD) and float64 occupancy code each
"""

from __future__ import annotations

import csv
import io
import json
import sys
import time
import zipfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "formats"))
sys.path.insert(0, str(HERE.parent / "acquisition"))
import remotezip  # noqa: E402
import shapefile  # noqa: E402
from tiff import GeoTIFF  # noqa: E402

CACHE = HERE.parent / "cache"
SOURCE = CACHE / "freetown-flood-risk"
OUT = CACHE / "products" / "freetown-flood"
MAPS_ZIP = "https://datacatalogfiles.worldbank.org/ddh-published/0066934/DR0095923/2023s0512_Freetown_Flood_Hazard_Maps.zip"
MEMBER = "2023s0512_Freetown_Flood_Hazard_Maps/2023s0512_JBA_Freetown_Flood_Hazard_Maps/Freetown_202304_{peril}_U_RP100_RD_1m_4326.tif"
PERILS = {"FLRF": "river", "FLSW": "surface water"}
RETURN_PERIODS = (20, 50, 100, 200, 500, 1500)
LEVEL = 1           # of the overview file: 10382, 5191, 2596, ... cells across
CATALOGUE_YEARS = 10000     # the event set's length: the sum of its event losses over its annual average is 10,000.0 (checked in the plug-in's checks)
RESULTS = "2023s0512_SLE_Freetown_03_gul_{group}_{table}.csv"


def grid_of(peril: str) -> dict:
    """A map's georeference: its corner and cell size, read from the map itself, kept once read."""
    keep = SOURCE / f"grid-{peril}.json"
    if keep.exists():
        return json.loads(keep.read_text())
    data = zipfile.ZipFile(remotezip.RemoteFile(MAPS_ZIP)).read(MEMBER.format(peril=peril))
    t = GeoTIFF(io.BytesIO(data))
    x0, dx, y0, dy = t.transform
    grid = {"x0": x0, "dx": dx, "y0": y0, "dy": dy, "width": t.main.width, "height": t.main.height, "epsg": t.epsg}
    keep.write_text(json.dumps(grid))
    return grid


def depth_map(peril: str, rp: int, grid: dict) -> dict:
    out = OUT / f"{peril}-rp{rp}.bin"
    tif = GeoTIFF.open(SOURCE / f"{peril}-RP{rp}.tif.ovr")
    image = tif.images[LEVEL]
    w, h = image.width, image.height
    if not out.exists():
        a = image.read().astype(np.float32)
        a = np.nan_to_num(a, nan=0, posinf=0, neginf=0)
        a[a < 0] = 0
        codes = np.clip(np.round(a / 0.05), 0, 255).astype(np.uint8)
        np.flipud(codes).tofile(out)
    dx, dy = grid["dx"] * grid["width"] / w, grid["dy"] * grid["height"] / h
    return {"peril": peril, "return_period": rp, "file": f"freetown-flood/{peril}-rp{rp}.bin", "width": w, "height": h,
            "box": [grid["y0"] - h * dy, grid["x0"], grid["y0"], grid["x0"] + w * dx]}


def table(group: str, name: str) -> list[dict]:
    with open(SOURCE / RESULTS.format(group=group, table=name), encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    maps = []
    for peril in PERILS:
        grid = grid_of(peril)
        for rp in RETURN_PERIODS:
            m = depth_map(peril, rp, grid)
            maps.append(m)
            print(f"{peril} RP{rp}: {m['width']} x {m['height']}", flush=True)
    # The buildings.
    rows = []
    with open(SOURCE / "Freetown_ptf_v2.csv", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f):
            try:
                rows.append((float(r["Latitude"]), float(r["Longitude"]), float(r["BuildingTIV"]) + float(r["ContentsTIV"]) + float(r["BITIV"]) + float(r["OtherTIV"]), float(r["OccupancyCode"])))
            except (TypeError, ValueError, KeyError):
                continue
    np.array(rows, dtype=np.float64).tofile(OUT / "buildings.bin")
    # The wards, and what the Council records.
    wards = []
    for ft in shapefile.read(SOURCE / "Freetown_Admin4_Combined_Flood_AAL.shp"):
        a = ft["attributes"]
        wards.append({"name": a.get("admin4Name"), "code": a.get("admin4Pcod"), "aal": a.get("AALoss_los"), "polygons": [{"outer": [[round(x, 6), round(y, 6)] for x, y in p["outer"]], "holes": [[[round(x, 6), round(y, 6)] for x, y in h] for h in p["holes"]]} for p in shapefile.polygons_of(ft)]})
    points = [{"at": [round(ft["parts"][0][0][1], 6), round(ft["parts"][0][0][0], 6)], "name": ft["attributes"].get("Location") or None} for ft in shapefile.read(SOURCE / "Freetown_City_Council_flash_flood_points.shp")]
    areas = [{"polygons": [{"outer": [[round(x, 6), round(y, 6)] for x, y in p["outer"]], "holes": []} for p in shapefile.polygons_of(ft)]} for ft in shapefile.read(SOURCE / "Freetown_City_Council_flash_flood_areas.shp")]
    # The loss tables.
    occupancy = {"1050": "residential", "1100": "commercial", "1150": "industrial"}
    number = lambda rows, key="loss": [float(r[key]) for r in rows]  # noqa: E731
    curves = {kind: {"return_period": number(table("S2", f"leccalc_wheatsheaf_mean_{kind}"), "return_period"), "loss": number(table("S2", f"leccalc_wheatsheaf_mean_{kind}"))} for kind in ("oep", "aep")}
    by = {}
    for kind in ("oep", "aep"):
        for r in table("S1", f"leccalc_wheatsheaf_mean_{kind}"):
            by.setdefault(occupancy.get(r["occupancycode"], r["occupancycode"]), {}).setdefault(kind, {"return_period": [], "loss": []})
            e = by[occupancy.get(r["occupancycode"], r["occupancycode"])][kind]
            e["return_period"].append(float(r["return_period"]))
            e["loss"].append(float(r["loss"]))
    events = sorted(number(table("S2", "eltcalc")), reverse=True)
    product = {
        "schema": "freetown-flood/1",
        "credit": "Freetown Flood Risk Assessment, World Bank, JBA Risk Management, Freetown City Council; CC BY 4.0",
        "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "perils": PERILS, "return_periods": RETURN_PERIODS, "maps": maps,
        "depth_code": {"metres_per_unit": 0.05, "note": "0 is dry or no data"},
        "buildings": {"file": "freetown-flood/buildings.bin", "count": len(rows), "fields": ["lat", "lon", "value_usd", "occupancy"], "occupancy": occupancy},
        "wards": wards, "reported_points": points, "reported_areas": areas,
        "reported_aal": {"all": float(table("S2", "aalcalc")[0]["loss"]), "by_occupancy": {occupancy.get(r["occupancycode"], r["occupancycode"]): float(r["loss"]) for r in table("S1", "aalcalc")}},
        "curves": curves, "curves_by_occupancy": by,
        "events": {"catalogue_years": CATALOGUE_YEARS, "count": len(events), "losses": [round(v, 2) for v in events]},
    }
    (CACHE / "products" / "freetown-flood.json").write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print(f"{len(rows)} buildings, {len(wards)} wards, {len(events)} events")
    return 0


if __name__ == "__main__":
    sys.exit(main())
