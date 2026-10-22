"""The API of The Science Infrastructures: the low layers of the repository (the
data library's products, the sensors' observations, the domains' algorithms and
the servers' plug-ins) answered over HTTP, for any program in any language, as
they happen. `Api` mounts on any standard-library server:

    api = Api(ROOT, allow_writes=True, token=None)
    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            if not api.handle(self): ...                 # paths outside /api/ are the server's own
        do_POST = do_OPTIONS = do_GET

The standalone service is server.py beside this file (port 8800); an
application's server mounts the same object (api/README.md §6). Every endpoint is
under /api/v1 and described in api/openapi.json; answers are JSON, errors are
{"error": "..."} with the status that says why.

The interface is ours and stable; what answers it may change (ADR 0003, and the
README's "own the interfaces, borrow the implementations"). Today the reads go
through science.Local (api/python/science/local.py), which the Python client's
Local also is; calls into the domains go to their JavaScript through a pool of
Node workers (compute.mjs); plug-ins run through science.hosting.

The write policy. Reads are open to every origin (CORS *). The writes are the
posts of observations: refused when the service is read-only; when it has a
token, allowed with `Authorization: Bearer <token>` only; without a token,
allowed from this computer's own programs and pages, and refused to a page of
another origin (a browser's request carries its Origin), so that no web page can
write to a service on this computer. Computing is not writing: calls into the
domains, samples asked by POST and plug-in runs are open to every caller, each
bounded in time, plug-ins one at a time.
"""

from __future__ import annotations

import collections
import datetime as dt
import hmac
import json
import math
import queue
import re
import sys
import threading
import time
import traceback
import urllib.parse
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "api" / "python"))
from science import bridge  # noqa: E402
from science.local import BadRequest, Busy, Local, NotFound, Unavailable, Unimplemented, instant, iso  # noqa: E402

VERSION = "0.1.0"
PREFIX = "/api/v1"
MAX_BODY = 8 << 20                # bytes of a request's body
MAX_PLACES = 100_000              # places sampled at once
MAX_JSON_CELLS = 2_000_000        # grid cells answered as JSON; more as float32
MAX_F32_CELLS = 40_000_000        # grid cells answered as float32
MAX_READINGS = 10_000             # readings posted at once
HEARTBEAT = 15.0                  # seconds between a stream's heartbeats
KEEP_EVENTS = 1000                # events of the stream kept for clients that reconnect
LOOPBACK = {"127.0.0.1", "localhost", "::1"}

# (method, path, summary): the index of the service, and the routes' order.
ENDPOINTS = [
    ("GET", "/", "this index: the version, the endpoints and the write policy"),
    ("GET", "/openapi.json", "the description of every endpoint (OpenAPI 3.1)"),
    ("GET", "/catalog", "the registry's datasets and the index of products"),
    ("GET", "/fields", "the field products: forecasts, analyses, relief, climate normals, plug-ins' results"),
    ("GET", "/fields/{name}", "a field's description: grid, times, variables"),
    ("GET", "/fields/{name}/sample", "values at places and an instant, as the Earth Atlas samples them (?var=&lon=&lat=&time=)"),
    ("POST", "/fields/{name}/sample", "the same, for many places: {var, lon: [...], lat: [...], time}"),
    ("GET", "/fields/{name}/series", "the value at a place at every time (?var=&lon=&lat=)"),
    ("GET", "/fields/{name}/grid", "the grid at a time, or a box of it (?var=&time=&bbox=&stride=&format=json|f32)"),
    ("GET", "/fields/{name}/export", "a file: GeoTIFF, netCDF, CSV or JSON (?var=&time=&format=&bbox=&stride=)"),
    ("GET", "/observations", "the observation sets and the sensors' sources"),
    ("GET", "/observations/{source}", "a source's readings or a set's rows (?since=&until=&bbox=&limit=)"),
    ("POST", "/observations/{source}", "readings from a device: {readings: [...], define?} (a write)"),
    ("GET", "/observations/{source}/stream", "a source's new readings as they arrive (Server-Sent Events)"),
    ("GET", "/stream", "every source's new readings as they arrive (Server-Sent Events)"),
    ("GET", "/events", "the event sets: earthquakes, mechanisms, disasters"),
    ("GET", "/events/{name}", "an event set (?since=&until=&bbox=&limit=)"),
    ("GET", "/tracks", "the track sets: tropical cyclones"),
    ("GET", "/tracks/{name}", "a track set"),
    ("GET", "/regions", "the sets of regions"),
    ("GET", "/regions/{name}", "a set of regions"),
    ("GET", "/regions/{name}/at", "the region at each place (?lon=&lat=)"),
    ("GET", "/verification", "the verifications of forecasts"),
    ("GET", "/verification/{name}", "a verification"),
    ("GET", "/compute", "the domains and foundations whose functions can be called"),
    ("GET", "/compute/{unit}", "one unit's functions, constants and modules"),
    ("POST", "/compute/{unit}/{function}", "a pure function called with {args: [...]}"),
    ("GET", "/plugins", "the plug-ins the service runs, by application"),
    ("GET", "/plugins/describe", "a plug-in's parameters and needs, with the sources for each (?app=&file=)"),
    ("POST", "/plugins/run", "a run of a plug-in: {app, file, chosen, params, context}; returns its job"),
    ("GET", "/jobs/{id}", "a job's state, progress and results"),
]


class HttpError(Exception):
    def __init__(self, status: int, message: str, headers: dict | None = None):
        super().__init__(message)
        self.status, self.message, self.headers = status, message, headers or {}


def finite(value):
    """JSON-safe: NaN and infinities as null."""
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {k: finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite(v) for v in value]
    return value


def numbers(text, name: str) -> list[float]:
    """Comma-separated numbers from a query, or a JSON list from a body."""
    if text is None or text == "":
        raise BadRequest(f"{name} is missing")
    parts = text if isinstance(text, list) else str(text).split(",")
    try:
        out = [float(p) for p in parts]
    except (TypeError, ValueError):
        raise BadRequest(f"{name} is numbers separated by commas, not {str(text)[:80]!r}") from None
    if not all(math.isfinite(x) for x in out):
        raise BadRequest(f"{name} holds a number that is not finite")
    return out


def values_list(array) -> list:
    """A numpy array as nested lists, NaN as None."""
    import numpy as np
    a = np.asarray(array, dtype=np.float64)
    out = a.tolist()
    if not np.isfinite(a).all():
        if a.ndim == 1:
            out = [x if x == x and abs(x) != math.inf else None for x in out]
        else:
            out = [[x if x == x and abs(x) != math.inf else None for x in row] for row in out]
    return out


class Hub:
    """New readings fanned out to every stream that listens, numbered, with the last
    few kept for clients that reconnect (they send the last id they saw)."""

    def __init__(self, store):
        self.lock = threading.Lock()
        self.seq = int(time.time() * 1000)       # ids that grow across restarts of the service
        self.recent: collections.deque = collections.deque(maxlen=KEEP_EVENTS)
        self.queues: set[queue.Queue] = set()
        self.unsubscribe = store.subscribe(self.publish)

    def publish(self, name: str, readings) -> None:
        with self.lock:
            self.seq += 1
            event = (self.seq, name, list(readings))
            self.recent.append(event)
            for q in list(self.queues):
                try:
                    q.put_nowait(event)
                except queue.Full:            # a client too slow to read: it loses events, not the service
                    pass

    def listen(self, after: int | None):
        q: queue.Queue = queue.Queue(maxsize=10_000)
        with self.lock:
            backlog = [e for e in self.recent if after is not None and e[0] > after]
            self.queues.add(q)
        return q, backlog

    def forget(self, q: queue.Queue) -> None:
        with self.lock:
            self.queues.discard(q)


class Api:
    """The API's endpoints, for a server that passes its requests to handle()."""

    def __init__(self, root: Path, *, allow_writes: bool = True, token: str | None = None, store=None,
                 read_only_reason: str | None = None, compute_seconds: float = 10.0, plugin_seconds: float = 600.0,
                 product_seconds: float = 60.0, heartbeat: float = HEARTBEAT):
        self.root = Path(root).resolve()
        self.local = Local(self.root, store=store, compute_seconds=compute_seconds)
        self.allow_writes = allow_writes
        self.token = token or None
        self.read_only_reason = read_only_reason or "this service is read-only"
        self.plugin_seconds = plugin_seconds
        self.product_seconds = product_seconds    # 0: posted readings are not built into products
        self.heartbeat = heartbeat                # seconds between a stream's heartbeats
        self._builds: dict[str, threading.Timer] = {}
        self._built: dict[str, float] = {}
        self.closing = threading.Event()
        self._runs = None
        self._hub = None
        self._lock = threading.Lock()
        self.routes = [(m, self._pattern(p), p) for m, p, _ in ENDPOINTS]

    # ---- the entry point ------------------------------------------------------------------------

    def handle(self, request) -> bool:
        """Answers a request under /api/; False for any other path, which the server answers itself."""
        url = urllib.parse.urlsplit(request.path)
        if not (url.path == "/api" or url.path.startswith("/api/")):
            return False
        try:
            if request.command == "OPTIONS":
                self._preflight(request)
                return True
            if not (url.path == PREFIX or url.path.startswith(PREFIX + "/")):
                if url.path.rstrip("/") == "/api":
                    self.reply(request, {"versions": {"v1": PREFIX + "/"}})
                    return True
                raise HttpError(404, f"no endpoint {url.path}: the API is at {PREFIX}/")
            sub = url.path[len(PREFIX):] or "/"
            allowed = []
            for method, pattern, template in self.routes:
                m = pattern.fullmatch(sub)
                if not m:
                    continue
                allowed.append(method)
                if method != request.command:
                    continue
                params = {k: urllib.parse.unquote(v) for k, v in m.groupdict().items()}
                query = dict(urllib.parse.parse_qsl(url.query, keep_blank_values=True))
                # An endpoint returns its body, or (body, status), or None when it has answered itself.
                answer = getattr(self, "on_" + self._name(method, template))(request, query, **params)
                if isinstance(answer, tuple):
                    self.reply(request, *answer)
                elif answer is not None:
                    self.reply(request, answer)
                return True
            if allowed:
                raise HttpError(405, f"{sub} answers {', '.join(sorted(set(allowed)))}", {"Allow": ", ".join(sorted(set(allowed)) + ["OPTIONS"])})
            raise HttpError(404, f"no endpoint {url.path}: see {PREFIX}/ for the list")
        except (BrokenPipeError, ConnectionResetError, ConnectionAbortedError):
            return True
        except Exception as error:     # every failure becomes an answer that says why
            status, message, headers = self._status_of(error)
            try:
                self.reply(request, {"error": message}, status, headers)
            except (BrokenPipeError, ConnectionResetError, ConnectionAbortedError):
                pass
            return True

    @staticmethod
    def _pattern(template: str):
        """A path template, /fields/{name}/sample, as a pattern with a group per parameter."""
        parts = re.split(r"\{(\w+)\}", template.rstrip("/"))
        text = "".join(f"(?P<{x}>[^/]+)" if k % 2 else re.escape(x) for k, x in enumerate(parts))
        return re.compile(text + "/?")

    @staticmethod
    def _name(method: str, template: str) -> str:
        """The method answering a route: GET /fields/{name}/sample -> on_get_fields_name_sample."""
        words = [w for w in re.sub(r"\{\w+\}", "name", template).replace(".", "_").split("/") if w]
        return (method.lower() + "_" + "_".join(words)) if words else method.lower() + "_index"

    @staticmethod
    def _status_of(error: Exception) -> tuple[int, str, dict]:
        if isinstance(error, HttpError):
            return error.status, error.message, error.headers
        if isinstance(error, (NotFound, LookupError)):
            return 404, str(error.args[0]) if error.args else "not found", {}
        if isinstance(error, Busy):
            return 409, str(error), {"Retry-After": "10"}
        if isinstance(error, bridge.Timeout):
            return 504, str(error), {}
        if isinstance(error, bridge.Failed):
            return 422, str(error), {}
        if isinstance(error, Unimplemented):
            return 501, str(error), {}
        if isinstance(error, Unavailable):
            return 503, str(error), {}
        if isinstance(error, ImportError):
            return 503, f"this needs a module the computer lacks: {error}", {}
        if isinstance(error, (BadRequest, ValueError, TypeError)):
            return 400, str(error), {}
        if isinstance(error, RuntimeError):
            return 503, str(error), {}
        traceback.print_exc()
        return 500, f"{type(error).__name__}: {error}", {}

    # ---- answers ---------------------------------------------------------------------------------

    @staticmethod
    def _common(request) -> None:
        request.send_header("Access-Control-Allow-Origin", "*")
        request.send_header("Access-Control-Expose-Headers", "X-Grid, Content-Disposition")
        request.send_header("X-Content-Type-Options", "nosniff")

    def reply(self, request, body, status: int = 200, headers: dict | None = None) -> None:
        try:
            data = json.dumps(body, ensure_ascii=False, allow_nan=False)
        except ValueError:
            data = json.dumps(finite(body), ensure_ascii=False)
        self.reply_bytes(request, data.encode("utf-8"), "application/json; charset=utf-8", status, headers)

    def reply_bytes(self, request, data: bytes, content_type: str, status: int = 200, headers: dict | None = None) -> None:
        request.send_response(status)
        request.send_header("Content-Type", content_type)
        request.send_header("Content-Length", str(len(data)))
        self._common(request)
        for k, v in (headers or {}).items():
            request.send_header(k, v)
        request.end_headers()
        if request.command != "HEAD":
            request.wfile.write(data)

    def _preflight(self, request) -> None:
        request.send_response(204)
        self._common(request)
        request.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        request.send_header("Access-Control-Allow-Headers", "Authorization, Content-Type, Last-Event-ID")
        request.send_header("Access-Control-Max-Age", "600")
        request.send_header("Content-Length", "0")
        request.end_headers()

    @staticmethod
    def body(request):
        length = request.headers.get("Content-Length")
        if length is None:
            raise HttpError(411, "a body needs its Content-Length")
        try:
            n = int(length)
        except ValueError:
            raise HttpError(400, "Content-Length is not a number") from None
        if n > MAX_BODY:
            raise HttpError(413, f"a body is at most {MAX_BODY >> 20} MB")
        raw = request.rfile.read(n) if n else b""
        if not raw:
            return {}
        try:
            return json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, ValueError) as error:
            raise HttpError(400, f"the body is not JSON: {error}") from None

    # ---- the write policy -----------------------------------------------------------------------------

    def policy(self, request=None) -> dict:
        refusal = self.write_refusal(request) if request is not None else None
        return {"policy": "closed" if not self.allow_writes else ("token" if self.token else "open"),
                "writable": refusal is None, "reason": refusal[1] if refusal else None,
                "writes": ["POST /observations/{source}"],
                "open": ["POST /compute/{unit}/{function}", "POST /fields/{name}/sample", "POST /plugins/run"]}

    def write_refusal(self, request) -> tuple[int, str] | None:
        """None when this request may write; else (status, reason)."""
        if not self.allow_writes:
            return 403, self.read_only_reason
        if self.token:
            given = request.headers.get("Authorization", "")
            if not hmac.compare_digest(given.encode("utf-8"), f"Bearer {self.token}".encode("utf-8")):
                return 401, "writing to this service needs its token: Authorization: Bearer <token>"
            return None
        origin = request.headers.get("Origin")
        if origin and urllib.parse.urlsplit(origin).netloc != request.headers.get("Host", ""):
            return 403, "a page of another origin may not write to a service without a token: start the service with --token"
        return None

    def _require_write(self, request) -> None:
        refusal = self.write_refusal(request)
        if refusal:
            raise HttpError(refusal[0], refusal[1], {"WWW-Authenticate": "Bearer"} if refusal[0] == 401 else {})

    # ---- the index ------------------------------------------------------------------------------------

    def on_get_index(self, request, query):
        bridge_reason = self.local.bridge().available()
        return {
            "name": "The Science Infrastructures API", "version": VERSION, "api": "v1", "base": PREFIX + "/",
            "documentation": PREFIX + "/openapi.json", "time": iso(dt.datetime.now(dt.timezone.utc)),
            "endpoints": [{"method": m, "path": PREFIX + p.rstrip("/"), "summary": s} for m, p, s in ENDPOINTS],
            "writes": self.policy(request),
            "compute": {"available": bridge_reason is None, "reason": bridge_reason},
            "sensors": {"available": self.local.store is not None},
        }

    def on_get_openapi_json(self, request, query):
        path = self.root / "api" / "openapi.json"
        if not path.is_file():
            raise NotFound("api/openapi.json is missing")
        self.reply_bytes(request, path.read_bytes(), "application/json; charset=utf-8")

    def on_get_catalog(self, request, query):
        return self.local.catalog()

    # ---- fields ---------------------------------------------------------------------------------------

    def on_get_fields(self, request, query):
        return {"fields": self.local.fields()}

    def on_get_fields_name(self, request, query, name):
        return self.local.field(name).meta

    def _sample(self, name: str, var, lon, lat, time):
        f = self.local.field(name)
        v = f.variable(var)
        lons, lats = numbers(lon, "lon"), numbers(lat, "lat")
        if max(len(lons), len(lats)) > MAX_PLACES:
            raise HttpError(413, f"at most {MAX_PLACES} places at once")
        when = instant(time)
        values = f.sample(var, lons, lats, when)
        b = f.bracket(var, when)
        return {"field": name, "variable": var, "units": v.get("units"), "interpolation": f.interpolation(var),
                "requested": iso(when), "time": b["time"], "times": b["times"], "weight": b["weight"],
                "lon": lons if len(lons) == len(values) else lons * len(values), "lat": lats if len(lats) == len(values) else lats * len(values),
                "values": values_list(values)}

    def on_get_fields_name_sample(self, request, query, name):
        return self._sample(name, query.get("var"), query.get("lon"), query.get("lat"), query.get("time"))

    def on_post_fields_name_sample(self, request, query, name):
        b = self.body(request)
        if not isinstance(b, dict):
            raise BadRequest("the body is {var, lon, lat, time}")
        return self._sample(name, b.get("var"), b.get("lon"), b.get("lat"), b.get("time"))

    def on_get_fields_name_series(self, request, query, name):
        f = self.local.field(name)
        lon, lat = numbers(query.get("lon"), "lon"), numbers(query.get("lat"), "lat")
        if len(lon) != 1 or len(lat) != 1:
            raise BadRequest("a series is at one place: one lon and one lat")
        s = f.series(query.get("var"), lon[0], lat[0])
        return {"field": name, "variable": s.variable, "units": s.units, "lon": lon[0], "lat": lat[0],
                "times": [iso(t) for t in s.times], "values": values_list(s.values)}

    def on_get_fields_name_grid(self, request, query, name):
        f = self.local.field(name)
        fmt = (query.get("format") or "json").lower()
        if fmt not in ("json", "f32"):
            raise BadRequest("format is json or f32")
        g = f.grid(query.get("var"), query.get("time"), query.get("bbox"), query.get("stride") or 1)
        cells = g.grid["nlat"] * g.grid["nlon"]
        described = {"field": name, "variable": g.variable, "units": g.units, "time": g.time, "grid": g.grid}
        if fmt == "f32":
            if cells > MAX_F32_CELLS:
                raise HttpError(413, f"{cells} cells: at most {MAX_F32_CELLS} at once; ask for a box or a stride")
            header = json.dumps({**described, "encoding": "float32 little-endian, row-major, rows north to south, NaN for no value"},
                                ensure_ascii=True, allow_nan=False)
            self.reply_bytes(request, g.values.astype("<f4").tobytes(), "application/octet-stream", headers={"X-Grid": header})
            return None
        if cells > MAX_JSON_CELLS:
            raise HttpError(413, f"{cells} cells: at most {MAX_JSON_CELLS} as JSON; ask for format=f32, a box or a stride")
        return {**described, "values": values_list(g.values)}

    def on_get_fields_name_export(self, request, query, name):
        f = self.local.field(name)
        data, content_type, filename = f.export(query.get("var"), query.get("time"), query.get("format") or "tif",
                                                query.get("bbox"), query.get("stride") or 1)
        safe = re.sub(r"[^\w.\-]", "_", filename)
        self.reply_bytes(request, data, content_type, headers={"Content-Disposition": f'attachment; filename="{safe}"'})
        return None

    # ---- observations ------------------------------------------------------------------------------------

    def on_get_observations(self, request, query):
        out = self.local.observations()
        return {**out, "streams": {"every source": PREFIX + "/stream", "one source": PREFIX + "/observations/{source}/stream"},
                "sensors_available": self.local.store is not None}

    def _geojson(self, request, query, kind: str, name: str) -> bool:
        """Answers a set as GeoJSON when ?format=geojson asks for it."""
        fmt = (query.get("format") or "json").lower()
        if fmt == "json":
            return False
        if fmt != "geojson":
            raise BadRequest("format is json or geojson")
        data, content_type, filename = self.local.geojson(kind, name, since=query.get("since"), until=query.get("until"),
                                                          bbox=query.get("bbox"), limit=query.get("limit"))
        self.reply_bytes(request, data, content_type, headers={"Content-Disposition": f'inline; filename="{filename}"'})
        return True

    def on_get_observations_name(self, request, query, source):
        is_sensor = any(s.get("name") == source for s in self.local.sensors())
        if not is_sensor and self._geojson(request, query, "observations", source):
            return None
        return self.local.observations(source, since=query.get("since"), until=query.get("until"),
                                       bbox=query.get("bbox"), limit=query.get("limit"))

    def on_post_observations_name(self, request, query, source):
        self._require_write(request)
        body = self.body(request)
        define, readings = None, body
        if isinstance(body, dict):
            define, readings = body.get("define"), body.get("readings", [])
            # A device may also describe its source beside its readings ({title, variables, readings}),
            # as observation/sensors/simulate.py does: it then needs no separate step to define it.
            if define is None and "variables" in body:
                define = {k: body[k] for k in ("title", "variables", "stations") if k in body}
            if "values" in body and "readings" not in body:
                readings = [body]                  # a single reading posted alone
        if not isinstance(readings, list) or not all(isinstance(r, dict) for r in readings):
            raise BadRequest("readings are a list of {station, lon, lat, time, values}")
        if len(readings) > MAX_READINGS:
            raise HttpError(413, f"at most {MAX_READINGS} readings at once")
        self._hub_of_store()                      # listen before writing, so that the stream sees these readings
        if define:
            if not isinstance(define, dict) or "variables" not in define:
                raise BadRequest("define is {title, variables, stations?}")
            self.local.define(source, define.get("title") or source, define["variables"], define.get("stations"))
        rejected: list = []
        accepted = self.local.append(source, readings, rejected)      # an unknown source is refused even without readings
        if accepted:
            self._build_soon(source)
        return {"source": source, "accepted": accepted, "defined": bool(define),
                "rejected": [{"index": i, "reason": why} for i, why in rejected]}

    def _build_soon(self, source: str) -> None:
        """Rebuilds the source's product (sensor-<source>, which the Earth Atlas draws)
        at most every `product_seconds`, after the readings that asked for it."""
        if not self.product_seconds:
            return
        with self._lock:
            if source in self._builds:
                return
            timer = threading.Timer(self.product_seconds if self._built.get(source) else 1.0, self._build, args=(source,))
            timer.daemon = True
            self._builds[source] = timer
        timer.start()

    def _build(self, source: str) -> None:
        try:
            self.local.build(source)
        except Exception:          # the readings are kept; only the product waits
            traceback.print_exc()
        finally:
            with self._lock:
                self._builds.pop(source, None)
                self._built[source] = time.time()

    def _hub_of_store(self) -> Hub:
        if self._hub is None:
            store = self.local.store
            if store is None:
                raise Unavailable("there is no sensors' store on this computer (observation/sensors/store.py)")
            with self._lock:
                if self._hub is None:
                    self._hub = Hub(store)
        return self._hub

    def on_get_observations_name_stream(self, request, query, source=None):
        hub = self._hub_of_store()
        last = request.headers.get("Last-Event-ID") or query.get("after")
        after = int(last) if last and str(last).isdigit() else None
        q, backlog = hub.listen(after)
        try:
            request.send_response(200)
            request.send_header("Content-Type", "text/event-stream; charset=utf-8")
            request.send_header("Cache-Control", "no-cache")
            request.send_header("X-Accel-Buffering", "no")
            self._common(request)
            request.end_headers()
            request.wfile.write(f": readings of {source or 'every source'}, a heartbeat every {self.heartbeat:g} s\nretry: 3000\n\n".encode("utf-8"))
            request.wfile.flush()
            for event in backlog:
                self._send_event(request, event, source)
            while not self.closing.is_set():
                try:
                    event = q.get(timeout=self.heartbeat)
                except queue.Empty:
                    request.wfile.write(f": heartbeat {iso(dt.datetime.now(dt.timezone.utc))}\n\n".encode("utf-8"))
                    request.wfile.flush()
                    continue
                self._send_event(request, event, source)
        except (BrokenPipeError, ConnectionResetError, ConnectionAbortedError, OSError):
            pass
        finally:
            hub.forget(q)
        request.close_connection = True
        return None

    def on_get_stream(self, request, query):
        return self.on_get_observations_name_stream(request, query, None)

    @staticmethod
    def _send_event(request, event, source) -> None:
        ident, name, readings = event
        if source and name != source:
            return
        data = json.dumps(finite({"source": name, "readings": readings}), ensure_ascii=False)
        request.wfile.write(f"id: {ident}\nevent: readings\ndata: {data}\n\n".encode("utf-8"))
        request.wfile.flush()

    # ---- events, tracks, regions, verification -----------------------------------------------------------------

    def on_get_events(self, request, query):
        return {"events": self.local.events()}

    def on_get_events_name(self, request, query, name):
        if self._geojson(request, query, "events", name):
            return None
        return self.local.events(name, since=query.get("since"), until=query.get("until"),
                                 bbox=query.get("bbox"), limit=query.get("limit"))

    def on_get_tracks(self, request, query):
        return {"tracks": self.local.tracks()}

    def on_get_tracks_name(self, request, query, name):
        if self._geojson(request, query, "tracks", name):
            return None
        return self.local.tracks(name)

    def on_get_regions(self, request, query):
        return {"regions": self.local.regions()}

    def on_get_regions_name(self, request, query, name):
        return self.local.regions(name)

    def on_get_regions_name_at(self, request, query, name):
        lons, lats = numbers(query.get("lon"), "lon"), numbers(query.get("lat"), "lat")
        if max(len(lons), len(lats)) > MAX_PLACES:
            raise HttpError(413, f"at most {MAX_PLACES} places at once")
        return {"regions": name, "lon": lons, "lat": lats, "at": self.local.regions_at(name, lons, lats)}

    def on_get_verification(self, request, query):
        return {"verification": self.local.verification()}

    def on_get_verification_name(self, request, query, name):
        return self.local.verification(name)

    # ---- the compute bridge -------------------------------------------------------------------------------

    def on_get_compute(self, request, query):
        return {"units": self.local.units(), "call": f"POST {PREFIX}/compute/{{unit}}/{{function}} with {{\"args\": [...]}}",
                "limit_seconds": self.local.compute_seconds}

    def on_get_compute_name(self, request, query, unit):
        return self.local.bridge().unit(unit)

    def on_post_compute_name_name(self, request, query, unit, function):
        body = self.body(request)
        args = body if isinstance(body, list) else (body.get("args", []) if isinstance(body, dict) else None)
        if not isinstance(args, list):
            raise BadRequest("the body is {\"args\": [...]} or a list of arguments")
        start = time.perf_counter()
        result = self.local.bridge().call(unit, function, args)
        return {"unit": unit, "function": function, "result": result, "seconds": round(time.perf_counter() - start, 6)}

    # ---- plug-ins and jobs ------------------------------------------------------------------------------------

    def runs(self):
        if self._runs is None:
            from science import hosting           # needs numpy
            with self._lock:
                if self._runs is None:
                    self._runs = hosting.Runs(self.root, seconds=self.plugin_seconds)
        return self._runs

    def on_get_plugins(self, request, query):
        return {"applications": {name: host.files() for name, host in self.runs().applications().items()},
                "runs": f"one at a time, at most {self.plugin_seconds:g} s each; open to every caller (they compute, and write only their own products)"}

    def on_get_plugins_describe(self, request, query):
        app, file = query.get("app") or "", query.get("file") or ""
        if not app or not file:
            raise BadRequest("name the application and the file: ?app=earth-atlas&file=plugins/vorticity.py")
        return {"app": app, "file": file, **self.runs().host(app).describe(file)}

    def on_post_plugins_run(self, request, query):
        body = self.body(request)
        if not isinstance(body, dict) or not body.get("app") or not body.get("file"):
            raise BadRequest("the body is {app, file, chosen?, params?, context?}")
        job = self.runs().start(body["app"], body["file"], body.get("chosen"), body.get("params"), body.get("context"))
        return self._job(job), 202

    def on_get_jobs_name(self, request, query, id):
        return self._job(self.runs().status(id))

    @staticmethod
    def _job(job: dict) -> dict:
        out = dict(job)
        out["results"] = [{**r, "api": f"{PREFIX}/fields/{r['field']}"} if r.get("kind") == "field" and r.get("field") else r
                          for r in job.get("results", [])]
        return out

    def close(self) -> None:
        """Ends the streams, the store's subscription and the compute workers."""
        self.closing.set()
        if self._hub is not None:
            self._hub.unsubscribe()
        self.local.close()
