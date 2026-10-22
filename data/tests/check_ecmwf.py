"""Checks of ECMWF's open forecasts in the data library: the client
(data/acquisition/ecmwf_open.py) on a synthetic index and synthetic responses,
the builder's conversions (data/products/ecmwf_fields.py), and, when they are in
the cache, the products fields-ecmwf-ifs-forecast, fields-ecmwf-aifs-forecast and
fields-ecmwf-ens-probabilities.

The bounds on the products are those of the present climate, not of one day:
the global mean of the 2 m temperature lies near 287 K, that of the pressure at
sea level near 1011 hPa (Trenberth and Smith 2005, J. Climate 18, 864-875), and
the global mean precipitation near 2.7 mm a day (GPCP, Adler et al. 2018,
Atmosphere 9, 138). At the start of a run, AIFS's state is the IFS analysis it
was initialised from, interpolated through its own grid, so the two differ at
+0 h by little more than the interpolation.
"""

from __future__ import annotations

import datetime as dt
import json
import math
import re
import struct

import numpy as np

from checks import PRODUCTS, assert_, check, skip

import ecmwf_open as eo
import ecmwf_fields as ef

UTC = dt.timezone.utc


# ---- the client, on synthetic data ---------------------------------------------

def line(param, step, offset, length, levtype="sfc", levelist=None, kind="fc"):
    e = {"domain": "g", "date": "20261022", "time": "0000", "type": kind, "stream": "oper", "step": step,
         "levtype": levtype, "param": param, "_offset": offset, "_length": length}
    if levelist is not None:
        e["levelist"] = str(levelist)
    return json.dumps(e)


index = "\n".join([
    line("2t", "6", 0, 100), line("msl", "6", 100, 50), line("tp", "6", 150, 25),
    line("z", "6", 400, 80, "pl", 500), line("u", "6", 175, 60, "pl", 500), line("u", "6", 900, 10, "pl", 250),
    "",
    line("tpg10", "96-120", 1000, 30, kind="ep"), line("tpg10", "120", 1030, 30, kind="ep"),
]) + "\n"
entries = eo.parse_index(index)
assert_("ecmwf: an index is read line by line, blank lines skipped", len(entries) == 8, f"{len(entries)} entries")
try:
    eo.parse_index('{"param": "2t", "_offset": 3}\n')
    assert_("ecmwf: an index line without _length is refused", False)
except ValueError:
    assert_("ecmwf: an index line without _length is refused", True)

fields = [{"param": "2t", "levtype": "sfc"}, {"param": "msl", "levtype": "sfc"}, {"param": "tp", "levtype": "sfc"},
          {"param": ["gh", "z"], "levtype": "pl", "levelist": 500}, {"param": "u", "levtype": "pl", "levelist": 500}]
chosen = eo.select(entries, fields, [6], "fc")
assert_("ecmwf: the fields asked for are selected in file order",
        [e["_offset"] for e in chosen] == [0, 100, 150, 175, 400], str([e["_offset"] for e in chosen]))
assert_("ecmwf: where gh is absent its alternative z is taken", any(e["param"] == "z" for e in chosen))
try:
    eo.select(entries, [{"param": "2t", "levtype": "sfc"}], [6, 12], "fc")
    assert_("ecmwf: a field missing at a step is named", False)
except ValueError as error:
    assert_("ecmwf: a field missing at a step is named", "2t at step 12 h" in str(error), str(error))
window = eo.select(entries, [{"param": "tpg10", "levtype": "sfc", "window": 24}], [120], "ep")
assert_("ecmwf: a quantity over 24 h ending at 120 h is the index's step 96-120",
        [e["step"] for e in window] == ["96-120"], str([e["step"] for e in window]))

ranges = eo.plan(chosen)
assert_("ecmwf: adjacent messages merge into one byte range, inclusive",
        ranges == [(0, 234), (400, 479)], str(ranges))
assert_("ecmwf: the Range header lists the runs", eo.range_header(ranges) == "bytes=0-234,400-479", eo.range_header(ranges))
assert_("ecmwf: a file's URL follows ECMWF's layout",
        eo.file_url("https://data.ecmwf.int/forecasts", dt.datetime(2026, 10, 22, tzinfo=UTC),
                    {"model": "ifs", "resolution": "0p25", "stream": "oper", "type": "fc"}, 6)
        == "https://data.ecmwf.int/forecasts/20261022/00z/ifs/0p25/oper/20261022000000-6h-oper-fc.grib2")
assert_("ecmwf: the index sits beside its file",
        eo.index_url("https://x/20261022000000-6h-oper-fc.grib2") == "https://x/20261022000000-6h-oper-fc.index")
assert_("ecmwf: probabilities to 240 h are in the file of step 240, beyond in that of 360",
        (eo.file_step(120, [240, 360]), eo.file_step(246, [240, 360]), eo.file_step(6, None)) == (240, 360, 6))
runs = eo.candidate_cycles(dt.datetime(2026, 10, 22, 15, 7, tzinfo=UTC), (0, 12), 3)
assert_("ecmwf: runs are tried from the latest at 00 or 12 UTC backwards",
        [f"{t:%d %H}" for t in runs] == ["22 12", "22 00", "21 12"], str([f"{t:%d %H}" for t in runs]))


def grib(payload: bytes) -> bytes:
    """A message with GRIB2's framing: 'GRIB', two reserved bytes, discipline, edition 2, length, ..., '7777'."""
    n = 16 + len(payload) + 4
    return b"GRIB\x00\x00\x00\x02" + struct.pack(">Q", n) + payload + b"7777"


a, b = grib(b"first --BOUNDARY inside the data"), grib(b"second")
body = (b"\r\n--BOUNDARY\r\nContent-Type: application/grib\r\nContent-Range: bytes 1000-%d/99\r\n\r\n" % (999 + len(a)) + a
        + b"\r\n--BOUNDARY\r\nContent-Type: application/grib\r\nContent-Range: bytes 5000-%d/99\r\n\r\n" % (4999 + len(b)) + b
        + b"\r\n--BOUNDARY--\r\n")
parts = eo.multipart_parts(body, "BOUNDARY")
assert_("ecmwf: a multipart/byteranges answer is cut by its Content-Range, not by the boundary inside the data",
        parts == {1000: a, 5000: b}, str({k: len(v) for k, v in parts.items()}))
got = eo.grib_message(parts, {"_offset": 5000, "_length": len(b), "param": "x"})
assert_("ecmwf: a message is taken from the part that holds it", got == b)
try:
    eo.grib_message({0: a[:-1] + b"8"}, {"_offset": 0, "_length": len(a), "param": "x"})
    assert_("ecmwf: bytes that are not one GRIB2 message are refused", False)
except eo.RangeError:
    assert_("ecmwf: bytes that are not one GRIB2 message are refused", True)


# ---- the builder's conversions -------------------------------------------------

check("ecmwf: 1 m of precipitation is 1000 kg m-2", abs(ef.factor("m", "kg m-2") - 1000.0), 1e-12, "kg m-2")
check("ecmwf: geopotential 5500 g0 m2 s-2 is 5500 m", abs(5500 * ef.STANDARD_GRAVITY * ef.factor("m**2 s**-2", "m") - 5500), 1e-9, "m")
rates = ef.mean_rates([np.array([0.0]), np.array([6.0]), np.array([5.99]), np.array([12.0])], [0, 6, 12, 18])
assert_("ecmwf: no rate at the start of the run", math.isnan(rates[0][0]))
check("ecmwf: 6 kg m-2 in 6 h is a mean rate of 1/3600 kg m-2 s-1 (1 mm an hour)", abs(rates[1][0] - 1 / 3600), 1e-15, "kg m-2 s-1")
assert_("ecmwf: a total that falls by the packing's rounding gives a rate of zero", rates[2][0] == 0.0, str(rates[2][0]))
lat, lon = ef.sampled_axes({"gridType": "regular_ll", "Ni": 1440, "Nj": 721, "latitudeOfFirstGridPointInDegrees": 90.0,
                            "longitudeOfFirstGridPointInDegrees": 180.0, "iDirectionIncrementInDegrees": 0.25,
                            "jDirectionIncrementInDegrees": 0.25, "iScansNegatively": 0, "jScansPositively": 0}, 4)
assert_("ecmwf: every 4th point of the 0.25° grid is the 1° grid from 90 N and 180 W",
        (len(lat), len(lon), lat[0], lat[-1], lon[0], lon[-1]) == (181, 360, 90.0, -90.0, -180.0, 179.0),
        f"{len(lat)} x {len(lon)}, {lat[0]}..{lat[-1]}, {lon[0]}..{lon[-1]}")


# ---- the products, when they are in the cache -------------------------------------

def product(name: str):
    path = PRODUCTS / f"fields-{name}.json"
    if not path.exists():
        return None
    meta = json.loads(path.read_text(encoding="utf-8"))
    meta["_codes"] = np.fromfile(PRODUCTS / meta["file"], "<i2")
    return meta


def values(meta: dict, name: str, t: int) -> np.ndarray:
    v = next(v for v in meta["variables"] if v["name"] == name)
    g = meta["grid"]
    n = g["nlon"] * g["nlat"]
    q = meta["_codes"][v["offsets"][t] // 2: v["offsets"][t] // 2 + n].reshape(g["nlat"], g["nlon"]).astype(np.float64)
    out = q * v["scale"] + v["offset"]
    out[q == -32768] = np.nan
    return out


def weights(meta: dict) -> np.ndarray:
    g = meta["grid"]
    lat = g["lat0"] - np.arange(g["nlat"]) * g["dlat"]
    return np.cos(np.radians(lat))[:, None] * np.ones((1, g["nlon"]))


def mean(meta, name, t):
    x, w = values(meta, name, t), weights(meta)
    ok = np.isfinite(x)
    return float((x[ok] * w[ok]).sum() / w[ok].sum())


def level_of(v: dict) -> str:
    """The catalogue's level (applications/earth-atlas/src/layers/catalogue.js, levelOf)."""
    if v.get("level"):
        return v["level"]
    m = re.search(r"\bat\s*(\d+(?:\.\d+)?)\s*m\b", v.get("long_name") or "", re.IGNORECASE)
    return f"{m.group(1)} m" if m else ""


QUANTITY = {"air_temperature": "air-temperature", "air_pressure_at_sea_level": "pressure",
            "geopotential_height": "height", "precipitation_flux": "precipitation"}


def keys(meta: dict) -> set[str]:
    """The catalogue's keys for the quantities checked here: speed:wind|<level> for
    vectors, <quantity>|<level>| for the scalars."""
    out = {f"speed:wind|{level_of(next(v for v in meta['variables'] if v['name'] == vec['u']))}"
           for vec in meta["vectors"] if vec["standard_name"] == "wind"}
    out |= {f"{QUANTITY[v['standard_name']]}|{level_of(v)}|" for v in meta["variables"] if v.get("standard_name") in QUANTITY}
    return out


EXPECTED = {"speed:wind|10 m", "speed:wind|500 hPa", "speed:wind|250 hPa", "air-temperature|2 m|", "air-temperature|850 hPa|",
            "pressure||", "height|500 hPa|", "precipitation||"}
gfs = [m for m in (product("gfs-forecast"), product("gfs-levels-forecast")) if m]
gfs_keys = set().union(*(keys(m) for m in gfs)) if gfs else set()
forecasts = {}
for label, name in (("IFS", "ecmwf-ifs-forecast"), ("AIFS", "ecmwf-aifs-forecast")):
    meta = product(name)
    if meta is None:
        skip(f"ecmwf: no product fields-{name} in the cache")
        continue
    forecasts[label] = meta
    n = len(meta["times"])
    t2 = float(np.mean([mean(meta, "t@2m", t) for t in range(n)]))
    check(f"ecmwf {label}: global mean 2 m temperature ({t2:.1f} K) within 280-292 K", abs(t2 - 286), 6, "K")
    p = float(np.mean([mean(meta, "msl", t) for t in range(n)])) / 100
    check(f"ecmwf {label}: global mean pressure at sea level ({p:.1f} hPa) within 1005-1015 hPa", abs(p - 1010), 5, "hPa")
    z = np.stack([values(meta, "gh@500hPa", t) for t in range(n)])
    lo, hi = np.percentile(z, 2), np.percentile(z, 98)
    assert_(f"ecmwf {label}: 500 hPa height between 4800 and 6000 m (2nd and 98th percentiles), extremes within 4500-6100 m",
            4800 <= lo and hi <= 6000 and 4500 <= z.min() and z.max() <= 6100,
            f"p02 {lo:.0f} m, p98 {hi:.0f} m, min {z.min():.0f} m, max {z.max():.0f} m")
    rain = np.stack([values(meta, "mtpr", t) for t in range(1, n)])
    assert_(f"ecmwf {label}: the precipitation rate is never negative", np.nanmin(rain) >= 0, f"min {np.nanmin(rain):.3g}")
    assert_(f"ecmwf {label}: no precipitation rate at the start of the run, one at every later time",
            not np.isfinite(values(meta, "mtpr", 0)).any() and np.isfinite(rain).all())
    daily = float(np.mean([mean(meta, "mtpr", t) for t in range(1, n)])) * 86400
    check(f"ecmwf {label}: global mean precipitation ({daily:.2f} mm a day) within 1.5-4.5 mm a day", abs(daily - 3.0), 1.5, "mm/day")
    g = meta["grid"]
    assert_(f"ecmwf {label}: the 1° grid of the other products, from 90 N and 180 W",
            (g["lon0"], g["dlon"], g["nlon"], g["lat0"], g["dlat"], g["nlat"], g["periodic"]) == (-180.0, 1.0, 360, 90.0, 1.0, 181, True),
            str(g))
    if gfs:
        assert_(f"ecmwf {label}: the grid is the GFS products' grid", all(m["grid"] == g for m in gfs))
    run = dt.datetime.fromisoformat(meta["source"]["run"].replace("Z", "+00:00"))
    expected_times = [f"{run + dt.timedelta(hours=h):%Y-%m-%dT%H:%M:%SZ}" for h in range(0, 121, 6)]
    assert_(f"ecmwf {label}: 21 times every 6 h from the run ({meta['source']['run']}) to +120 h, in UTC",
            meta["times"] == expected_times and meta["time_scale"] == "UTC", f"{n} times, {meta['times'][0]} to {meta['times'][-1]}")
    assert_(f"ecmwf {label}: the model is named for the atlas", meta.get("model") == f"ECMWF {label}", str(meta.get("model")))
    have = keys(meta)
    assert_(f"ecmwf {label}: wind at 10, 500, 250 hPa, temperature at 2 m and 850 hPa, pressure, height and precipitation as the catalogue keys them",
            EXPECTED <= have, f"missing {sorted(EXPECTED - have)}")
    if gfs:
        assert_(f"ecmwf {label}: each of these quantities is also the GFS's, so the atlas puts the models side by side",
                EXPECTED <= gfs_keys, f"not in the GFS products: {sorted(EXPECTED - gfs_keys)}")
    assert_(f"ecmwf {label}: CC BY 4.0 and ECMWF's attribution travel with the product",
            "CC BY 4.0" in meta["source"]["licence"] and "based on data and products of the European Centre" in meta["source"]["credit"])

if len(forecasts) == 2:
    ifs, aifs = forecasts["IFS"], forecasts["AIFS"]
    if ifs["source"]["run"] != aifs["source"]["run"]:
        skip(f"ecmwf: IFS ({ifs['source']['run']}) and AIFS ({aifs['source']['run']}) are of different runs; not compared")
    else:
        w = weights(ifs)
        def rms(name):
            d = values(ifs, name, 0) - values(aifs, name, 0)
            return float(np.sqrt((d * d * w).sum() / w.sum()))
        check("ecmwf: IFS and AIFS at +0 h, root-mean-square difference of 500 hPa height", rms("gh@500hPa"), 30, "m")
        check("ecmwf: IFS and AIFS at +0 h, root-mean-square difference of pressure at sea level", rms("msl"), 300, "Pa")

ens = product("ecmwf-ens-probabilities")
if ens is None:
    skip("ecmwf: no product fields-ecmwf-ens-probabilities in the cache")
else:
    n = len(ens["times"])
    names = [v["name"] for v in ens["variables"]]
    every = np.stack([values(ens, v, t) for v in names for t in range(n)])
    assert_("ecmwf ENS: every probability lies within 0-100 %", np.nanmin(every) >= -1e-3 and np.nanmax(every) <= 100 + 1e-3,
            f"{np.nanmin(every):.3f} to {np.nanmax(every):.3f}")
    assert_("ecmwf ENS: every probability is in percent", all(v["units"] == "%" for v in ens["variables"]))
    # An exceedance probability does not grow with its threshold: P(X >= 50) <= P(X >= 20) <= P(X >= 10).
    # ECMWF packs each probability in steps of a power of two, from 2^-5 % (tpg10, tpg20) to 2^-1 %
    # (tpg50, 10fgg10), rounding each value to within half its step; two values compared may then
    # disagree by up to half the sum of their steps, at most 0.5 %.
    for chain in (["tpg10", "tpg20", "tpg50"], ["fgg10", "fgg15", "fgg25"]):
        worst = max(float(np.nanmax(values(ens, hi, t) - values(ens, lo, t))) for lo, hi in zip(chain, chain[1:]) for t in range(n))
        check(f"ecmwf ENS: {' >= '.join(chain)} at every point and time, within the packing (largest excess)", max(worst, 0.0), 0.5, "%")
    both = max(float(np.nanmax(values(ens, "ptsa_gt_1stdev@850hPa", t) + values(ens, "ptsa_lt_1stdev@850hPa", t))) for t in range(n))
    check("ecmwf ENS: a temperature is not both above and below normal: P(above) + P(below) at most 100 %, within the packing (2^-2 %)",
          max(both - 100, 0.0), 0.25, "%")
    run = dt.datetime.fromisoformat(ens["source"]["run"].replace("Z", "+00:00"))
    expected_times = [f"{run + dt.timedelta(hours=h):%Y-%m-%dT%H:%M:%SZ}" for h in range(24, 121, 12)]
    assert_("ecmwf ENS: 24 h periods ending every 12 h from +24 to +120 h", ens["times"] == expected_times,
            f"{n} times, {ens['times'][0]} to {ens['times'][-1]}")
    assert_("ecmwf ENS: the model is named for the atlas", ens.get("model") == "ECMWF ENS", str(ens.get("model")))
    if forecasts.get("IFS"):
        assert_("ecmwf ENS: on the IFS product's grid", ens["grid"] == forecasts["IFS"]["grid"])
