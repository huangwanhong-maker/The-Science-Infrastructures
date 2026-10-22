"""Vegetation seen by Sentinel-2: the normalised difference vegetation index,

    NDVI = (ρ_NIR − ρ_red) / (ρ_NIR + ρ_red)

(Rouse, Haas, Schell and Deering 1974, "Monitoring vegetation systems in the
Great Plains with ERTS", NASA SP-351, 309-317; Tucker 1979, Remote Sensing of
Environment 8, 127-150), of the least
cloudy Level-2A acquisition of the last days over a box around the chosen
place, at 20 or 40 m, on a regular longitude–latitude grid (WGS 84). Green
leaves absorb red light for photosynthesis and scatter the near infrared, so
the index runs from about −0.2 over water through 0.1–0.2 over bare soil and
towns to 0.8–0.9 over dense, healthy canopy.

The scenes come from the STAC API of Earth Search (data/acquisition/stac.py),
and their bands from Cloud-Optimised GeoTIFFs read by ranges
(data/formats/tiff.py): red (B04, 665 nm) and near infrared (B08, 842 nm) from
the overview at the resolution asked, and the scene classification (SCL, 20 m)
of Sen2Cor (ESA, Sentinel-2 MSI Level-2A Algorithm Theoretical Basis Document,
S2-PDGS-MPC-ATBD-L2A), which masks no data, saturated pixels, cloud shadows,
clouds of medium and high probability, and thin cirrus. An acquisition is the
tiles of one datatake that meet the box; they are joined, each in its own UTM
zone, by the transverse Mercator projection (Karney 2011). Reflectance is
ρ = DN × scale + offset, with the scale and offset that each item states for
its band (the raster extension of STAC), the offset since processing baseline
04.00 being −0.1 (BOA_ADD_OFFSET = −1000). Earth Search states, in
earthsearch:boa_offset_applied, when it has already removed that offset from
the pixels; the declared offset is then not applied again. The darkest pixels
check the choice: no surface reflects less than −0.05, so an offset that would
make the darkest 1 % of a band darker than that is not in the data, and the
note says so.

Optionally, the acquisition of about a year before, over the same box, gives
the change of NDVI: the same season, so that the change is of the land rather
than of the phenological cycle.

Credit: Contains modified Copernicus Sentinel data (the years of the scenes).
"""

import datetime as dt
import math
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "data" / "acquisition"))
sys.path.insert(0, str(ROOT / "data" / "formats"))
import stac  # noqa: E402
import tiff  # noqa: E402

UTC = dt.timezone.utc
MASKED = {0: "no data", 1: "saturated or defective", 3: "cloud shadow", 8: "cloud, medium probability",
          9: "cloud, high probability", 10: "thin cirrus"}           # classes of the SCL band
CANDIDATES = 3            # acquisitions whose classification is read before one is chosen
DARKEST = -0.05           # no surface reflects less; L2A's negative values are of order -0.01
READERS = 4               # files read at once


def box_grid(lon: float, lat: float, size_m: float, resolution_m: float) -> dict:
    """A regular longitude–latitude grid of about `resolution_m` covering a box of
    `size_m` around a place. Distances are measured on a transverse Mercator
    projection centred on the place with unit scale there (Karney 2011), exact
    along its meridian and within 10⁻⁶ across a box of 40 km."""
    local = tiff.TransverseMercator(central_meridian=lon, k0=1.0)
    _, y = local.forward(lon, lat)
    half = size_m / 2
    east, _ = local.inverse(half, y)
    _, north = local.inverse(0.0, y + half)
    _, south = local.inverse(0.0, y - half)
    dlat = float(north - south) / size_m * resolution_m
    dlon = float(east - lon) / half * resolution_m
    nlat = int(round(float(north - south) / dlat)) + 1
    nlon = int(round(2 * float(east - lon) / dlon)) + 1
    return {"lon0": round(lon - (nlon - 1) / 2 * dlon, 9), "dlon": round(dlon, 12), "nlon": nlon,
            "lat0": round(float(north), 9), "dlat": round(dlat, 12), "nlat": nlat, "periodic": False,
            "crs": "WGS 84 geographic coordinates (EPSG:4326), degrees", "order": "rows north to south, columns eastward"}


def nodes(grid: dict) -> tuple[np.ndarray, np.ndarray]:
    lons = grid["lon0"] + np.arange(grid["nlon"]) * grid["dlon"]
    lats = grid["lat0"] - np.arange(grid["nlat"]) * grid["dlat"]
    return np.meshgrid(lons, lats)


def bounds(grid: dict) -> tuple[float, float, float, float]:
    """(west, south, east, north) of a grid's nodes."""
    return (grid["lon0"], grid["lat0"] - (grid["nlat"] - 1) * grid["dlat"],
            grid["lon0"] + (grid["nlon"] - 1) * grid["dlon"], grid["lat0"])


def at_nodes(values: np.ndarray, transform, x: np.ndarray, y: np.ndarray, nearest: bool = False) -> np.ndarray:
    """Values of a window (rows down, transform (x0, dx, y0, dy) of its outer corner)
    at places (x, y) of its CRS: the pixel containing each place, or bilinear between
    the four pixel centres around it, with the weights of missing pixels (NaN) given
    to the others when they carry at least half (as foundations/physis/js/grid.js does).
    NaN outside the window."""
    x0, dx, y0, dy = transform
    rows, cols = values.shape
    c, r = (x - x0) / dx, (y0 - y) / dy
    outside = (c < 0) | (r < 0) | (c >= cols) | (r >= rows)
    if nearest:
        out = values[np.clip(r.astype(int), 0, rows - 1), np.clip(c.astype(int), 0, cols - 1)].astype(np.float64)
        out[outside] = np.nan
        return out
    c, r = c - 0.5, r - 0.5
    c0, r0 = np.floor(c).astype(int), np.floor(r).astype(int)
    fc, fr = c - c0, r - r0
    total, weight = np.zeros(x.shape), np.zeros(x.shape)
    for dr, dc, w in ((0, 0, (1 - fc) * (1 - fr)), (0, 1, fc * (1 - fr)), (1, 0, (1 - fc) * fr), (1, 1, fc * fr)):
        v = values[np.clip(r0 + dr, 0, rows - 1), np.clip(c0 + dc, 0, cols - 1)].astype(np.float64)
        ok = np.isfinite(v)
        total += np.where(ok, v * w, 0)
        weight += np.where(ok, w, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(weight >= 0.5, total / weight, np.nan)
    out[outside] = np.nan
    return out


def reflectance(dn: np.ndarray, item: dict, key: str) -> tuple[np.ndarray, str, float]:
    """The reflectance of an L2A band's digital numbers (0 is no data); the rule
    used, in words, with its reason; and the reflectance of the darkest 1 %."""
    band = stac.raster_band(item, key)
    if "scale" not in band:
        raise ValueError(f"{item['id']}: the item states no scale for {key}")
    p = item["properties"]
    scale, declared = float(band["scale"]), float(band.get("offset", 0.0))
    applied = p.get("earthsearch:boa_offset_applied") is True
    offset = 0.0 if applied else declared
    reason = ("declared, but already removed from the pixels, as earthsearch:boa_offset_applied says" if applied and declared
              else "declared and applied")
    valid = dn > 0
    darkest = float(np.percentile(dn[valid], 1)) * scale if valid.sum() > 100 else math.nan
    if offset and darkest + offset < DARKEST:
        reason = f"declared, but it would make the darkest 1 % of pixels {darkest + offset:.3f}: not in the data, so not applied"
        offset = 0.0
    rho = np.where(valid, dn * scale + offset, np.nan)
    sign = "+" if offset >= 0 else "−"
    rule = f"ρ = DN × {scale:g} {sign} {abs(offset):g} (offset {declared:+g} {reason}; baseline {p.get('s2:processing_baseline')})"
    return rho, rule, darkest + offset


def rules_note(name: str, found: list[dict]) -> str:
    """One sentence for the reflectance rules of an acquisition's tiles: {band: (rule, darkest)} per tile."""
    rules = {band: sorted({f[band][0] for f in found}) for band in ("red", "nir")}
    if rules["red"] == rules["nir"] and len(rules["red"]) == 1:
        said = f"red and near infrared, {rules['red'][0]}"
    else:
        said = "; ".join(f"{label}, {' or '.join(rules[band])}" for band, label in (("red", "red"), ("nir", "near infrared")))
    dark = {band: [f[band][1] for f in found] for band in ("red", "nir")}
    span = lambda xs: f"{min(xs):.3f}" if max(xs) - min(xs) < 5e-4 else f"{min(xs):.3f} to {max(xs):.3f}"  # noqa: E731
    return (f"Reflectance of {name}: {said}. "
            f"Darkest 1 % of pixels: red {span(dark['red'])}, near infrared {span(dark['nir'])}.")


class Acquisition:
    """The tiles of one datatake that meet the box, the newest processing of each."""

    def __init__(self, items: list[dict]):
        best = {}
        for item in items:
            p = item["properties"]
            key = p.get("grid:code") or item["id"]
            if key not in best or str(p.get("s2:processing_baseline")) > str(best[key]["properties"].get("s2:processing_baseline")):
                best[key] = item
        self.items = sorted(best.values(), key=lambda it: it["properties"].get("eo:cloud_cover", 100))
        self.time = min(stac.when(it) for it in self.items)
        self.cloud = float(np.mean([it["properties"].get("eo:cloud_cover", 100) for it in self.items]))
        self.classes = None

    @property
    def name(self) -> str:
        return " + ".join(it["id"] for it in self.items)

    @property
    def label(self) -> str:
        n = len(self.items)
        return f"the scene of {self.time:%Y-%m-%d %H:%M} UTC ({n} tile{'s' if n > 1 else ''})"

    def classify(self, LON: np.ndarray, LAT: np.ndarray) -> np.ndarray:
        """The scene class at each node (0 where no tile has data): the first tile with a class there."""
        out = np.zeros(LON.shape, dtype=np.int16)
        self.windows = {}
        with ThreadPoolExecutor(READERS) as pool:
            windows = list(pool.map(lambda item: classes_of(item, LON, LAT), self.items))
        for item, window in zip(self.items, windows):
            if window is None:
                continue
            proj, classes, scl, scl_transform = window
            take = (out == 0) & np.isfinite(classes) & (classes > 0)
            out[take] = classes[take].astype(np.int16)
            self.windows[item["id"]] = (proj, take, scl, scl_transform)
        self.classes = out
        return out

    def masked_share(self) -> tuple[float, float]:
        """The shares of the box masked by the classification and outside every tile."""
        outside = float((self.classes == 0).mean())
        masked = float(np.isin(self.classes, [k for k in MASKED if k]).mean())
        return masked, outside

    def ndvi(self, LON: np.ndarray, LAT: np.ndarray, resolution: float, notes: list[str]) -> np.ndarray:
        """NDVI at the nodes, each from the tile that classified it; NaN where masked or outside."""
        out = np.full(LON.shape, np.nan)
        used = [it for it in self.items if it["id"] in self.windows and self.windows[it["id"]][1].any()]
        with ThreadPoolExecutor(READERS) as pool:
            done = list(pool.map(lambda item: ndvi_of(item, self.windows[item["id"]], LON, LAT, resolution), used))
        for item, (values, _) in zip(used, done):
            take = self.windows[item["id"]][1]
            out[take] = values[take]
        if done:
            notes.append(rules_note(self.label, [said for _, said in done]))
        out[np.isin(self.classes, list(MASKED))] = np.nan
        return out


def classes_of(item: dict, LON: np.ndarray, LAT: np.ndarray):
    """A tile's scene classes at the nodes (NaN outside it), with its projection and
    window of the SCL band; None if the box misses the tile's raster."""
    scl = tiff.GeoTIFF.url(stac.asset_href(item, "scl"))
    proj = tiff.utm(scl.epsg)
    x, y = proj.forward(LON, LAT)
    try:
        values, transform = scl.read_box(float(x.min()), float(y.min()), float(x.max()), float(y.max()), level=0)
    except tiff.TIFFError:
        return None
    return proj, at_nodes(values, transform, x, y, nearest=True), values, transform


def ndvi_of(item: dict, window, LON: np.ndarray, LAT: np.ndarray, resolution: float) -> tuple[np.ndarray, dict]:
    """A tile's NDVI at the nodes it classified, and how its reflectances were found: {band: (rule, darkest)}."""
    proj, take, scl, scl_transform = window
    x, y = proj.forward(LON, LAT)
    box = (float(x[take].min()) - 2 * resolution, float(y[take].min()) - 2 * resolution,
           float(x[take].max()) + 2 * resolution, float(y[take].max()) + 2 * resolution)
    red_tif = tiff.GeoTIFF.url(stac.asset_href(item, "red"))
    nir_tif = tiff.GeoTIFF.url(stac.asset_href(item, "nir"))
    dn_red, transform = red_tif.read_box(*box, level=red_tif.level_for(resolution))
    dn_nir, transform_nir = nir_tif.read_box(*box, level=nir_tif.level_for(resolution))
    if dn_red.shape != dn_nir.shape or not np.allclose(transform, transform_nir):
        raise ValueError(f"{item['id']}: the red and near-infrared windows differ")
    red, rule_red, dark_red = reflectance(dn_red.astype(np.float64), item, "red")
    nir, rule_nir, dark_nir = reflectance(dn_nir.astype(np.float64), item, "nir")
    with np.errstate(invalid="ignore", divide="ignore"):
        index = np.where(nir + red > 1e-4, (nir - red) / (nir + red), np.nan)
    # Masked pixels are left out before interpolating, so that a cloud's edge does not reach a clear node.
    x0, dx, y0, dy = transform
    PX, PY = np.meshgrid(x0 + (np.arange(index.shape[1]) + 0.5) * dx, y0 - (np.arange(index.shape[0]) + 0.5) * dy)
    index[np.isin(at_nodes(scl, scl_transform, PX, PY, nearest=True), [k for k in MASKED if k])] = np.nan
    return at_nodes(index, transform, x, y), {"red": (rule_red, dark_red), "nir": (rule_nir, dark_nir)}


def acquisitions(grid: dict, start: dt.datetime, end: dt.datetime) -> list[Acquisition]:
    """The acquisitions of [start, end] meeting the box, least cloudy first."""
    west, south, east, north = bounds(grid)
    body = stac.search_body([stac.SENTINEL_2_L2A["collection"]],
                            intersects={"type": "Polygon", "coordinates": [[[west, south], [east, south], [east, north], [west, north], [west, south]]]},
                            start=start, end=end, query={"eo:cloud_cover": {"lt": 95}}, sortby=[("eo:cloud_cover", "asc")], limit=100)
    groups: dict[str, list[dict]] = {}
    for item in stac.search(stac.EARTH_SEARCH, body, max_items=300):
        p = item["properties"]
        # One datatake, whatever its processing: GS2A_20211024T012651_033055_N05.00 and its N03.01.
        key = re.sub(r"_N\d\d\.\d\d$", "", p.get("s2:datatake_id") or f"{p.get('platform')}-{p['datetime'][:10]}")
        groups.setdefault(key, []).append(item)
    return sorted((Acquisition(items) for items in groups.values()), key=lambda a: a.cloud)


def choose(found: list[Acquisition], LON, LAT, progress) -> Acquisition | None:
    """The acquisition with the least of the box masked or outside, among the least cloudy few."""
    best, best_share = None, 2.0
    for k, acq in enumerate(found[:CANDIDATES]):
        progress(f"classification of {acq.label}, {k + 1} of {min(CANDIDATES, len(found))}")
        acq.classify(LON, LAT)
        masked, outside = acq.masked_share()
        if masked + outside < best_share:
            best, best_share = acq, masked + outside
        if best_share < 0.03:
            break
    return best


def create(api):
    class NDVI(api.Plugin):
        id = "ndvi"
        title = "Vegetation from Sentinel-2 (NDVI)"
        description = ("The normalised difference vegetation index of the least cloudy Sentinel-2 scene of the last days, "
                       "over a box around the chosen place, at 20 or 40 m; clouds and shadows masked by the scene classification.")
        domain = "geography"
        parameters = {
            "days": {"label": "Search the last", "default": 30, "min": 5, "max": 120, "step": 5, "unit": "d"},
            "size": {"label": "Box", "default": 20, "min": 4, "max": 40, "step": 2, "unit": "km"},
            "resolution": {"label": "Resolution", "default": 20, "choices": [[20, "20 m"], [40, "40 m"]]},
            "compare": {"label": "Change since", "default": 0, "choices": [[0, "none"], [1, "a year before"]]},
        }

        def run(self, inputs, params, context):
            if not context.pinned:
                return [api.note("Choose a place first (click the Earth), then run again: the scene is read for a box around it.")]
            lon, lat = float(context.pinned[0]), float(context.pinned[1])
            days, size = float(params.get("days", 30)), float(params.get("size", 20)) * 1000
            resolution = float(params.get("resolution", 20))
            compare = float(params.get("compare", 0)) >= 1
            end = min(context.time, dt.datetime.now(UTC))
            grid = box_grid(lon, lat, size, resolution)
            LON, LAT = nodes(grid)
            say = lambda fraction, text: self.breathe(context, fraction, text)  # noqa: E731

            say(0.05, "searching Earth Search")
            found = acquisitions(grid, end - dt.timedelta(days=days), end)
            if not found:
                return [api.note(f"No Sentinel-2 Level-2A scene meets the box around {lon:.3f}°, {lat:.3f}° "
                                 f"in the {days:g} days before {api.iso(end)}; search a longer period.")]
            now = choose(found, LON, LAT, lambda text: say(0.15, text))
            notes: list[str] = []
            say(0.4, f"red and near infrared of {now.label}")
            ndvi = now.ndvi(LON, LAT, resolution, notes)
            masked, outside = now.masked_share()
            scenes = [now]

            variables = [{"name": "ndvi", "long_name": "Normalised difference vegetation index", "standard_name": None,
                          "units": "1", "cell_methods": "time: point", "style": {"colormap": "speed", "range": [-0.2, 0.9]}, "slices": [ndvi]}]
            rows = [[now.name, api.iso(now.time), f"{now.cloud:.1f} %", f"{masked:.1%}", f"{outside:.1%}"]]
            if compare:
                then = now.time - dt.timedelta(days=365)
                say(0.6, "searching a year before")
                before = choose(acquisitions(grid, then - dt.timedelta(days=days / 2), then + dt.timedelta(days=days / 2)),
                                LON, LAT, lambda text: say(0.65, text))
                if before is None:
                    notes.append(f"No scene meets the box within {days / 2:g} days of {api.iso(then)[:10]}: no change computed.")
                else:
                    say(0.75, f"red and near infrared of {before.label}")
                    earlier = before.ndvi(LON, LAT, resolution, notes)
                    b_masked, b_outside = before.masked_share()
                    rows.append([before.name, api.iso(before.time), f"{before.cloud:.1f} %", f"{b_masked:.1%}", f"{b_outside:.1%}"])
                    variables.append({"name": "ndvi_change", "long_name": f"Change of NDVI since {api.iso(before.time)[:10]}",
                                      "standard_name": None, "units": "1", "cell_methods": "time: point", "difference": True, "slices": [ndvi - earlier]})
                    scenes.append(before)

            years = sorted({s.time.year for s in scenes})
            credit = stac.SENTINEL_2_L2A["credit"].format(year=", ".join(map(str, years)))
            ids = [it["id"] for s in scenes for it in s.items]
            made = api.make_product(
                self.product_name, self.title, grid, [now.time], variables, model=f"Sentinel-2, {api.iso(now.time)[:10]}",
                domain=self.domain,
                source={"dataset": stac.SENTINEL_2_L2A["collection"], "provider": stac.SENTINEL_2_L2A["provider"],
                        "credit": credit, "licence": stac.SENTINEL_2_L2A["licence"], "citation": stac.SENTINEL_2_L2A["citation"],
                        "scenes": ids, "date": api.iso(now.time), "request": f"{stac.EARTH_SEARCH}/search",
                        "retrieved_at": api.iso(dt.datetime.now(UTC))})
            say(1.0, "done")
            finite = ndvi[np.isfinite(ndvi)]
            summary = (f"NDVI of {now.label}, cloud cover of the tiles {now.cloud:.1f} %: "
                       f"{masked:.1%} of the box masked by clouds, shadows or bad pixels, {outside:.1%} outside the scene; "
                       + (f"median {np.median(finite):.2f}, 2–98 % {np.percentile(finite, 2):.2f} to {np.percentile(finite, 98):.2f}."
                          if finite.size else "no clear pixel."))
            return [api.field(made, {"surface": "ndvi"}),
                    api.table("Scenes", ["Scene", "Time (UTC)", "Tile cloud", "Box masked", "Outside"], rows),
                    api.note(summary), *[api.note(n) for n in notes], api.note(credit + ".")]

    return NDVI()
