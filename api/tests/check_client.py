"""Checks of the API's Python side, run by api/tests/run-tests.mjs against a
service it has started: the client against the service, the service against
science.Local reading the same cache (they must agree exactly, being one
implementation), the stream of readings, the write policy's rule, the
interpolation of science/fields.py against cases known in closed form, and the
plug-in context's readers. One JSON object a line, as data/tests/checks.py
reports: {"kind": "check" | "assert" | "skip", ...}.

    python api/tests/check_client.py http://127.0.0.1:<port> <sensor source>
"""

from __future__ import annotations

import datetime as dt
import json
import math
import sys
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "api" / "python"))
sys.path.insert(0, str(ROOT / "api" / "service"))
import numpy as np  # noqa: E402

import science  # noqa: E402
from science import fields, plugins  # noqa: E402

UTC = dt.timezone.utc


def emit(record: dict) -> None:
    print(json.dumps(record, ensure_ascii=False), flush=True)


def check(name: str, value: float, limit: float, unit: str = "") -> None:
    emit({"kind": "check", "name": name, "value": float(value), "limit": float(limit), "unit": unit})


def assert_(name: str, condition: bool, detail: str = "") -> None:
    emit({"kind": "assert", "name": name, "ok": bool(condition), "detail": str(detail)})


def skip(what: str) -> None:
    emit({"kind": "skip", "name": what})


def max_difference(a, b) -> float:
    """The largest difference of two arrays, NaN matching NaN (else infinite)."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != b.shape or not np.array_equal(np.isnan(a), np.isnan(b)):
        return math.inf
    ok = ~np.isnan(a)
    return float(np.max(np.abs(a[ok] - b[ok]))) if ok.any() else 0.0


def interpolation() -> None:
    """science/fields.py against closed forms: bilinear interpolation reproduces a
    function linear in longitude and latitude; missing neighbours are renormalised
    while they carry less than half the weight; places beyond the grid have no value."""
    grid = {"lon0": 10.0, "dlon": 2.0, "nlon": 6, "lat0": 50.0, "dlat": 1.0, "nlat": 5, "periodic": False}
    lons = grid["lon0"] + np.arange(6) * 2.0
    lats = grid["lat0"] - np.arange(5) * 1.0
    LON, LAT = np.meshgrid(lons, lats)
    values = (3 + 0.5 * LON - 2 * LAT).astype(np.float32)
    x, y = np.array([11.3, 17.9, 10.0, 20.0, 9.2]), np.array([48.6, 46.1, 50.0, 46.4, 47.0])
    got = fields.bilinear(values, grid, x, y)
    want = 3 + 0.5 * np.clip(x, 10, 20) - 2 * y
    check("fields: bilinear reproduces a linear function (and clamps within half a cell of the edge)", float(np.max(np.abs(got - want))), 1e-5)
    holed = values.copy()
    holed[1, 0] = np.nan                     # one corner of the cell around (10.5, 48.6) lost
    a = float(fields.bilinear(holed, grid, 10.5, 48.6))
    w = {(1, 0): 0.75 * 0.6, (1, 1): 0.25 * 0.6, (2, 0): 0.75 * 0.4, (2, 1): 0.25 * 0.4}
    kept = {k: v for k, v in w.items() if k != (1, 0)}
    expect = sum(v * float(values[k]) for k, v in kept.items()) / sum(kept.values())
    check("fields: a missing neighbour is renormalised away", abs(a - expect), 1e-5)
    holed[2, 0] = np.nan                     # now the valid corners carry 0.25 of the weight
    assert_("fields: no value where the valid neighbours carry less than half the weight", math.isnan(float(fields.bilinear(holed, grid, 10.5, 48.6))))
    assert_("fields: no value beyond half a cell of a grid that does not wrap",
            bool(np.all(np.isnan(fields.bilinear(values, grid, [8.9, 21.1, 15.0], [48.0, 48.0, 51.6])))))
    ring = {"lon0": -180.0, "dlon": 90.0, "nlon": 4, "lat0": 0.0, "dlat": 1.0, "nlat": 1, "periodic": True}
    v = np.array([[0.0, 1.0, 2.0, 3.0]], dtype=np.float32)
    check("fields: a periodic grid wraps round the antimeridian", float(fields.bilinear(v, ring, 135.0, 0.0)) - 1.5, 1e-12)
    times = [0.0, 3600e3, 7200e3]
    check("fields: the time bracket weighs the later time", abs(fields.time_bracket(times, 900e3)[2] - 0.25), 1e-15)
    assert_("fields: the time bracket clamps to the span", fields.time_bracket(times, -5.0)[:2] == (0, 0) and fields.time_bracket(times, 1e9)[:2] == (2, 2))
    assert_("fields: the nearest time for a time mean", fields.time_bracket(times, 2700e3, nearest=True) == (1, 1, 0.0))


def main(argv: list[str]) -> int:
    url, source = argv[0], argv[1]
    interpolation()
    api = science.Client(url)
    local = science.Local(ROOT)

    info = api.info()
    assert_("python client: the service's index", info.get("api") == "v1" and info.get("version") == science.__version__, info.get("version"))

    # The export of files answers 501 where data/formats/export.py is missing.
    import tempfile
    import handler
    with tempfile.TemporaryDirectory() as empty:
        try:
            science.local.export_module(Path(empty))
            assert_("export: a computer without data/formats/export.py says so", False)
        except science.local.Unimplemented as error:
            assert_("export: a computer without data/formats/export.py answers 501", handler.Api._status_of(error)[0] == 501, str(error))

    # The write policy's rule, as the standalone service applies it.
    import server
    assert_("policy: a service on another address without a token accepts no writes", server.write_policy("0.0.0.0", None, False)[0] is False)
    assert_("policy: with a token it accepts writes (with the token)", server.write_policy("0.0.0.0", "t0ken", False) == (True, None))
    assert_("policy: read-only refuses writes, token or not", server.write_policy("127.0.0.1", "t0ken", True)[0] is False)

    # Fields: the service and Local are one implementation, so they agree exactly.
    names = [f["name"] for f in local.fields()]
    if "gfs-forecast" in names:
        remote, here = api.field("gfs-forecast"), local.field("gfs-forecast")
        lon = [139.69, -0.13, 179.8, -179.9, 0.0, -73.97]
        lat = [35.69, 51.51, -0.2, 10.5, 89.9, 40.78]
        t = here.times[2] + dt.timedelta(hours=2)
        a = remote.sample("tmp2m", lon, lat, time=t)
        b = here.sample("tmp2m", lon, lat, t)
        check("python client: samples equal Local's", max_difference(a, b), 0.0, "K")
        assert_("python client: one place gives a float", isinstance(remote.sample("tmp2m", 139.69, 35.69, time=t), float))
        s1, s2 = remote.series("tmp2m", 139.69, 35.69), here.series("tmp2m", 139.69, 35.69)
        check("python client: a series equals Local's", max_difference(s1.values, s2.values), 0.0, "K")
        assert_("python client: a series has every time", len(s1.times) == len(here.times) and s1.times[0] == here.times[0])
        g1 = remote.grid("tmp2m", time=3, bbox=(170, -10, -170, 10), stride=2)
        g2 = here.grid("tmp2m", 3, (170, -10, -170, 10), 2)
        check("python client: a grid across the antimeridian equals Local's", max_difference(g1.values, g2.values), 0.0, "K")
        assert_("python client: the grid is a numpy array of its shape", isinstance(g1.values, np.ndarray) and g1.values.shape == (g1.grid["nlat"], g1.grid["nlon"]),
                str(getattr(g1.values, "shape", None)))
    else:
        skip("python client: fields-gfs-forecast is not in the cache")

    events = [e["name"] for e in local.events()]
    if events:
        e1, e2 = api.events(events[0], limit=5), local.events(events[0], limit=5)
        assert_(f"python client: the event set {events[0]}, its latest five", e1["events"] == e2["events"] and len(e1["events"]) == min(5, e2["of"]))
    else:
        skip("python client: no event set in the cache")

    # The compute bridge: through the service, and through Local's own workers.
    jd = 2460000.5
    try:
        r1 = api.compute("astronomy", "sunPosition", jd)
        r2 = local.compute("astronomy", "sunPosition", jd)
        assert_("python client: the Sun's position through the service equals Local's", r1 == r2, json.dumps(r1))
        check("python client: wind from the west is 270 degrees", api.compute("meteorology", "windFromDirection", 5, 0) - 270, 1e-12, "degree")
    except science.ApiError as error:
        skip(f"python client: compute: {error}")
    finally:
        local.close()

    try:
        api.field("no-such-field")
        assert_("python client: a missing field is an error", False)
    except science.ApiError as error:
        assert_("python client: a missing field is a 404", error.status == 404, error.message)

    # Readings posted, and received by a stream in another thread within a second.
    got, opened = [], threading.Event()

    def listen():
        try:
            stream = api.stream(source)
            opened.set()
            for event in stream:
                got.append((time.perf_counter(), event))
                return
        except Exception as error:      # reported below
            got.append((time.perf_counter(), {"error": repr(error)}))
    thread = threading.Thread(target=listen, daemon=True)
    thread.start()
    opened.wait(5)
    time.sleep(0.8)                          # the generator connects on its first step
    start = time.perf_counter()
    now = dt.datetime.now(UTC).replace(microsecond=0)
    posted = api.post_observations(source, [{"station": "py1", "lon": 139.76, "lat": 35.68, "time": now, "values": {"air_temperature": 291.0}}],
                                   define={"title": "API check (Python)", "variables": ["air_temperature"]})
    thread.join(5)
    assert_("python client: readings posted are accepted", posted.get("accepted") == 1, json.dumps(posted))
    if got and "readings" in got[0][1]:
        check("python client: a reading reaches the stream", got[0][0] - start, 1.0, "s")
        assert_("python client: the stream carries the reading", got[0][1]["readings"][0]["station"] == "py1")
    else:
        assert_("python client: a reading reaches the stream", False, got[0][1] if got else "nothing within 5 s")
    readings = api.observations(source, limit=1)["readings"]
    assert_("python client: the latest reading reads back", readings and readings[-1]["station"] == "py1", json.dumps(readings)[:200])

    # The plug-ins' context reads the library's other products.
    from science import hosting
    library = hosting.Library()
    context = plugins.Context(time=dt.datetime.now(UTC), listing=library.list, events=library.events, observations=library.observations)
    assert_("plug-in context: list('fields') names the library's fields", len(context.list("fields")) == len(local.list("fields")))
    if events:
        e = context.events(events[0])
        assert_("plug-in context: events(name) reads an event set", bool(e) and "fields" in e and "events" in e)
    observed = local.list("observations")
    if observed:
        o = context.observations(observed[0]["name"])
        assert_("plug-in context: observations(name) reads an observation set", bool(o) and o.get("schema", "").startswith("observations"))
    else:
        skip("plug-in context: no observation set in the cache")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
