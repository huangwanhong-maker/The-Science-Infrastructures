"""The client of the API (api/README.md, api/openapi.json), on the standard
library alone; numpy, when present, makes grids and samples at many places arrays.

    api = Client("http://127.0.0.1:8800")                 # or $SCIENCE_API_URL; token from $SCIENCE_API_TOKEN
    gfs = api.field("gfs-forecast")
    gfs.sample("tmp2m", 139.69, 35.69, time="now+1d")      # one place: a float (NaN where there is no value)
    gfs.series("tmp2m", 139.69, 35.69)                     # Series(times, values, units, variable)
    gfs.grid("tmp2m", time=0, bbox=(120, 20, 150, 50))     # Grid(values, grid, time, variable, units)
    api.compute("astronomy", "sunPosition", 2460000.5)     # a domain's function, run by the service
    for event in api.stream("garden"):                     # readings as they arrive (Server-Sent Events)
        print(event["readings"])

The service may be the API's own (api/service/server.py, port 8800) or any
application's server that mounts it (the Earth Atlas's, at port 8790): give
either's address, with or without /api/v1. Times are UTC: a datetime, ISO 8601
text, `now`, or `now+6h`.
"""

from __future__ import annotations

import datetime as dt
import email.message
import http.client
import json
import math
import os
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from array import array
from pathlib import Path
from typing import Iterator

from .local import Grid, Series

try:
    import numpy as np
except ImportError:          # numpy is optional
    np = None

UTC = dt.timezone.utc
DEFAULT_URL = "http://127.0.0.1:8800"


class ApiError(Exception):
    """An answer that is an error: its HTTP status (None when the service could not be reached) and message."""

    def __init__(self, status: int | None, message: str, url: str = "", job: dict | None = None):
        super().__init__(f"{status or 'unreachable'}: {message}" + (f" ({url})" if url else ""))
        self.status, self.message, self.url, self.job = status, message, url, job


def _text(value) -> str:
    """A query value: datetimes as ISO 8601, sequences comma-separated, numbers as written."""
    if isinstance(value, dt.datetime):
        value = value.replace(tzinfo=UTC) if value.tzinfo is None else value.astimezone(UTC)
        return value.strftime("%Y-%m-%dT%H:%M:%SZ")
    if isinstance(value, (list, tuple)) or (np is not None and isinstance(value, np.ndarray)):
        return ",".join(_text(v) for v in (value.tolist() if hasattr(value, "tolist") else value))
    if isinstance(value, float):
        return repr(value)
    return str(value)


def _values(raw: list, many: bool):
    """Numbers from JSON, null as NaN: a float for one place, an array (or list) for many."""
    out = [math.nan if v is None else float(v) for v in raw]
    if not many:
        return out[0]
    return np.array(out) if np is not None else out


def _time(text):
    return dt.datetime.fromisoformat(text.replace("Z", "+00:00")) if text else None


class Client:
    """A service of the API: fields, observations, events, tracks, regions,
    verification, the compute bridge and plug-ins."""

    def __init__(self, url: str | None = None, token: str | None = None, *, timeout: float = 60.0):
        base = (url or os.environ.get("SCIENCE_API_URL") or DEFAULT_URL).rstrip("/")
        if base.endswith("/api"):
            base += "/v1"
        elif not base.endswith("/api/v1"):
            base += "/api/v1"
        self.url = base
        self.token = token if token is not None else os.environ.get("SCIENCE_API_TOKEN")
        self.timeout = timeout

    # ---- requests ---------------------------------------------------------------------------------

    def _address(self, path: str, query: dict | None = None) -> str:
        q = {k: _text(v) for k, v in (query or {}).items() if v is not None and v != ""}
        return self.url + path + (("?" + urllib.parse.urlencode(q)) if q else "")

    def _open(self, method: str, path: str, query: dict | None = None, body=None, headers: dict | None = None,
              timeout: float | None = None):
        address = self._address(path, query)
        data = json.dumps(body, allow_nan=False).encode("utf-8") if body is not None else None
        request = urllib.request.Request(address, data=data, method=method, headers={"Accept": "application/json", **(headers or {})})
        if data is not None:
            request.add_header("Content-Type", "application/json")
        if self.token:
            request.add_header("Authorization", f"Bearer {self.token}")
        try:
            return urllib.request.urlopen(request, timeout=timeout or self.timeout)
        except urllib.error.HTTPError as error:
            raw = error.read()
            try:
                message = json.loads(raw.decode("utf-8")).get("error") or raw.decode("utf-8", "replace")
            except ValueError:
                message = raw.decode("utf-8", "replace")[:500] or error.reason
            raise ApiError(error.code, message, address) from None
        except (urllib.error.URLError, OSError) as error:
            raise ApiError(None, f"the service cannot be reached: {getattr(error, 'reason', error)}", address) from None

    def _json(self, method: str, path: str, query: dict | None = None, body=None):
        with self._open(method, path, query, body) as response:
            return json.loads(response.read().decode("utf-8"))

    def get(self, path: str, **query):
        """Any GET of the API, by its path under /api/v1, as JSON."""
        return self._json("GET", path, query)

    # ---- the service and its catalogue ----------------------------------------------------------------

    def info(self) -> dict:
        """The version, the endpoints and the write policy."""
        return self._json("GET", "/")

    def catalog(self) -> dict:
        """The registry's datasets (id, kind, title, licence, parameters) and the index of products."""
        return self._json("GET", "/catalog")

    # ---- fields -------------------------------------------------------------------------------------------

    def fields(self) -> list[dict]:
        return self._json("GET", "/fields")["fields"]

    def field(self, name: str) -> "RemoteField":
        return RemoteField(self, self._json("GET", f"/fields/{urllib.parse.quote(name)}"))

    # ---- observations and their stream --------------------------------------------------------------------

    def observations(self, source: str | None = None, *, since=None, until=None, bbox=None, limit=None):
        """The observation sets and sensor sources (no source), or one source's readings or rows."""
        if source is None:
            return self._json("GET", "/observations")
        return self._json("GET", f"/observations/{urllib.parse.quote(source)}",
                          {"since": since, "until": until, "bbox": bbox, "limit": limit})

    def post_observations(self, source: str, readings: list[dict], define: dict | None = None) -> dict:
        """Sends readings [{station, lon, lat, time, values: {variable: value}}] to a sensor
        source, defining it first when `define` is {title, variables, stations?}."""
        body = {"readings": [_reading(r) for r in readings]}
        if define:
            body["define"] = define
        return self._json("POST", f"/observations/{urllib.parse.quote(source)}", body=body)

    def stream(self, source: str | None = None, *, last_event_id: int | None = None, reconnect: bool = True,
               retry: float = 3.0, silence: float = 45.0) -> Iterator[dict]:
        """New readings as they arrive, as {id, source, readings}: a generator over the
        service's Server-Sent Events. Heartbeats come every 15 s; after `silence` seconds
        without one, or when the connection drops, it reconnects (when `reconnect`) and
        asks for what it missed by the last event's id."""
        path = f"/observations/{urllib.parse.quote(source)}/stream" if source else "/stream"
        last = last_event_id
        while True:
            headers = {"Accept": "text/event-stream"}
            if last is not None:
                headers["Last-Event-ID"] = str(last)
            try:
                response = self._open("GET", path, headers=headers, timeout=silence)
            except ApiError as error:
                if error.status is not None or not reconnect:
                    raise
                time.sleep(retry)
                continue
            try:
                for event in _server_sent_events(response):
                    if event.get("retry"):
                        retry = event["retry"] / 1000
                        continue
                    if event["id"] is not None:
                        last = event["id"]
                    if event["event"] == "readings":
                        yield {"id": event["id"], **json.loads(event["data"])}
            except (OSError, socket.timeout, http.client.HTTPException, ValueError):
                if not reconnect:
                    raise
            finally:
                response.close()
            if not reconnect:
                return
            time.sleep(retry)

    # ---- events, tracks, regions, verification --------------------------------------------------------------

    def events(self, name: str | None = None, *, since=None, until=None, bbox=None, limit=None):
        """The event sets (no name), or one: {fields, events: [[row]], window, source, ...}."""
        if name is None:
            return self._json("GET", "/events")["events"]
        return self._json("GET", f"/events/{urllib.parse.quote(name)}", {"since": since, "until": until, "bbox": bbox, "limit": limit})

    def tracks(self, name: str | None = None):
        return self._json("GET", "/tracks")["tracks"] if name is None else self._json("GET", f"/tracks/{urllib.parse.quote(name)}")

    def regions(self, name: str | None = None):
        return self._json("GET", "/regions")["regions"] if name is None else self._json("GET", f"/regions/{urllib.parse.quote(name)}")

    def verification(self, name: str | None = None):
        if name is None:
            return self._json("GET", "/verification")["verification"]
        return self._json("GET", f"/verification/{urllib.parse.quote(name)}")

    # ---- the compute bridge ---------------------------------------------------------------------------------

    def units(self, unit: str | None = None):
        """The domains and foundations that can be called, with their functions; or one unit."""
        if unit is None:
            return self._json("GET", "/compute")["units"]
        return self._json("GET", f"/compute/{urllib.parse.quote(unit)}")

    def compute(self, unit: str, function: str, *args):
        """A pure exported function of a domain or foundation, called with JSON arguments; its result."""
        return self._json("POST", f"/compute/{urllib.parse.quote(unit)}/{urllib.parse.quote(function)}",
                          body={"args": [_plain(a) for a in args]})["result"]

    # ---- plug-ins and jobs --------------------------------------------------------------------------------------

    def plugins(self) -> dict:
        """The plug-ins the service can run, by application."""
        return self._json("GET", "/plugins")["applications"]

    def describe_plugin(self, app: str, file: str) -> dict:
        return self._json("GET", "/plugins/describe", {"app": app, "file": file})

    def job(self, ident: str) -> dict:
        return self._json("GET", f"/jobs/{urllib.parse.quote(ident)}")

    def run_plugin(self, app: str, file: str, chosen: dict | None = None, params: dict | None = None,
                   context: dict | None = None, *, wait: bool = True, poll: float = 1.0, timeout: float | None = None) -> dict:
        """Runs a server plug-in: its sources `chosen` ({need: "product|level"}; the
        described defaults otherwise) and parameters; returns the job, when `wait` once
        it is done. A failed run raises ApiError with the job."""
        body = {"app": app, "file": file, "chosen": chosen or {}, "params": params or {}, "context": context or {}}
        job = self._json("POST", "/plugins/run", body=body)
        if not wait:
            return job
        start = time.monotonic()
        while job["state"] in ("starting", "running"):
            if timeout is not None and time.monotonic() - start > timeout:
                raise ApiError(None, f"job {job['id']} is still running after {timeout:g} s", job=job)
            time.sleep(poll)
            job = self.job(job["id"])
        if job["state"] == "failed":
            raise ApiError(None, job.get("message") or "the run failed", job=job)
        return job


class RemoteField:
    """A field product of a service: its description, and its values at places,
    over time and on its grid, computed by the service as the Earth Atlas computes them."""

    def __init__(self, client: Client, meta: dict):
        self.client, self.meta = client, meta
        self.name = meta["name"]
        self.variables = {v["name"]: v for v in meta["variables"]}
        self.times = [_time(t) for t in meta.get("times", [])]
        self._path = f"/fields/{urllib.parse.quote(self.name)}"

    def __repr__(self) -> str:
        return f"<field {self.name}: {', '.join(self.variables)}; {len(self.times)} times>"

    def sample(self, var: str, lon, lat, time=None):
        """Values at places, bilinear in space and linear in time; a float for one place, NaN where there is none."""
        many = isinstance(lon, (list, tuple)) or isinstance(lat, (list, tuple)) or (np is not None and (np.ndim(lon) or np.ndim(lat)))
        r = self.client._json("GET", self._path + "/sample", {"var": var, "lon": lon, "lat": lat, "time": time})
        return _values(r["values"], bool(many))

    def series(self, var: str, lon: float, lat: float) -> Series:
        """The value at a place at every time of the product."""
        r = self.client._json("GET", self._path + "/series", {"var": var, "lon": lon, "lat": lat})
        return Series([_time(t) for t in r["times"]], _values(r["values"], True), r.get("units"), var)

    def grid(self, var: str, time=None, bbox=None, stride: int = 1) -> Grid:
        """The variable on its grid, or a box of it, every `stride`-th point, at a stored
        time (an index) or an instant: rows north to south, NaN where there is none; a
        numpy array (nlat, nlon) when numpy is present, else a list of rows."""
        query = {"var": var, "time": time, "bbox": bbox, "stride": stride if stride != 1 else None, "format": "f32"}
        with self.client._open("GET", self._path + "/grid", query, headers={"Accept": "application/octet-stream"}) as response:
            body = response.read()
            described = json.loads(response.headers.get("X-Grid") or "{}")
        g = described["grid"]
        if np is not None:
            values = np.frombuffer(body, dtype="<f4").reshape(g["nlat"], g["nlon"]).copy()
        else:
            flat = array("f")
            flat.frombytes(body)
            if sys.byteorder == "big":
                flat.byteswap()
            values = [list(flat[r * g["nlon"]:(r + 1) * g["nlon"]]) for r in range(g["nlat"])]
        return Grid(values, g, described.get("time"), var, described.get("units"))

    def export(self, var: str, time=None, format: str = "tif", bbox=None, stride: int = 1, path=None):
        """The variable at a stored time as a file: GeoTIFF (tif), netCDF (nc), CSV or JSON.
        The bytes, or, when `path` is a file or directory, the path written."""
        query = {"var": var, "time": time, "format": format, "bbox": bbox, "stride": stride if stride != 1 else None}
        with self.client._open("GET", self._path + "/export", query, headers={"Accept": "*/*"}) as response:
            data = response.read()
            disposition = response.headers.get("Content-Disposition", "")
        if path is None:
            return data
        target = Path(path)
        if target.is_dir():
            message = email.message.Message()
            message["Content-Disposition"] = disposition
            target = target / (message.get_filename() or f"{self.name}-{var}.{format}")
        target.write_bytes(data)
        return target


def _server_sent_events(response) -> Iterator[dict]:
    """The events of a text/event-stream: {id, event, data} and {retry}; comments (heartbeats) are skipped."""
    ident, kind, data = None, "message", []
    while True:
        raw = response.readline()
        if not raw:
            return
        line = raw.decode("utf-8").rstrip("\r\n")
        if not line:
            if data:
                yield {"id": ident, "event": kind, "data": "\n".join(data)}
            kind, data = "message", []
            continue
        if line.startswith(":"):
            continue
        name, _, value = line.partition(":")
        value = value[1:] if value.startswith(" ") else value
        if name == "data":
            data.append(value)
        elif name == "event":
            kind = value
        elif name == "id":
            ident = int(value) if value.isdigit() else value
        elif name == "retry" and value.isdigit():
            yield {"retry": int(value)}


def _reading(r: dict) -> dict:
    out = dict(r)
    if isinstance(out.get("time"), dt.datetime):
        out["time"] = _text(out["time"])
    return out


def _plain(value):
    """An argument as JSON: numpy arrays and scalars as lists and numbers, datetimes as ISO 8601."""
    if np is not None:
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
    if isinstance(value, dt.datetime):
        return _text(value)
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value
