# api: The Interface to the Low Layers

*Status: draft (2026-10-23).* An edge of the architecture. The service and the Python package use the standard library and numpy; the compute bridge needs Node 18 or later; the clients need nothing beyond their language.

The low layers of the infrastructure (the algorithms of the domains, the data library's products, the observations of our own sensors, the plug-ins of the servers) are reachable here from any program, in Python, C and C++ or JavaScript, on this computer or across a network, as they change. The API is the stable interface; what stands behind it may change.

## 1. Purpose

The applications reach the low layers through their own pages and servers. A study in a notebook, a program in C, a microcontroller in a garden or another person's web page cannot. This edge gives them one interface:

- **reads** of the library's products: fields sampled at places and instants, series, grids and files; events, tracks, regions, verification; observations;
- **writes** of readings by sensors, and their **stream** to every listener as they arrive;
- **calls** of the domains' pure functions, and **runs** of the servers' plug-ins.

*Rationale.* The repository's rule is to own the interfaces and borrow the implementations ([architecture §1.5](../docs/architecture.md)). The paths and shapes of `/api/v1` are ours and are kept; behind them the domains are provisional JavaScript today and will be C++ ([ADR 0003](../development/decisions/0003-cpp-computational-core.md)), and the reader of products may become another. A client written against the API does not notice the change.

## 2. Concepts

### 2.1 Resources

Every resource is a product of the data library ([data/README.md §4](../data/README.md)), named as the index of products names it, or a source of the sensors' store.

| Resource | What it is | Read as |
|---|---|---|
| **Fields** | Variables on a regular longitude–latitude grid over time: forecasts, analyses, relief, climate normals, and the products plug-ins have made | a description; values at places and an instant; a series at a place; a grid or a box of it, as JSON or float32; a file (GeoTIFF, netCDF, CSV, JSON) |
| **Observations** | Readings at places: aerodromes (METAR), buoys (NDBC), and our own sensors | a set's rows, or a source's readings, within a span, a box and a limit; GeoJSON; a stream |
| **Events** | Things at a place and an instant: earthquakes, mechanisms, disasters | the set, filtered; GeoJSON |
| **Tracks** | Moving things as positions over time: tropical cyclones | the set; GeoJSON |
| **Regions** | Countries and provinces rasterised on a grid | the set; the region at each place |
| **Verification** | Scores of the forecasts against analyses and observations | the set |

A field is sampled exactly as the Earth Atlas's page samples it ([fields.js](../applications/earth-atlas/src/data/fields.js) over physis's [grid.js](../foundations/physis/js/grid.js)): bilinear between the four grid points around a place, the weights of those with a value renormalised when they carry at least half of the total; linear between the two stored times around an instant, clamped to the product's span; the nearest time for a variable whose cell method is a mean, sum, maximum or minimum over time; the instant's month for a climatology; directions through their unit vectors. The accuracy is that of the product: half a quantum (the `quantum` of each variable, the step of its stored codes), plus the error of linear interpolation.

*Remark.* Positions are longitude and latitude in degrees on WGS 84, as the products write them; times are UTC, given as ISO 8601 (an instant without an offset is UTC), as `now`, or as `now+6h`; values are in the units the product states, SI unless it says otherwise. Degrees are the boundary's convention; the domains work in radians inside.

### 2.2 The Compute Bridge to the Domains

`POST /api/v1/compute/{unit}/{function}` calls an exported function of a domain or foundation (`domains/<d>/js/index.js`, `foundations/<f>/js/index.js`) with JSON arguments, and answers its result as JSON. The calls run in a pool of two persistent Node workers ([service/compute.mjs](service/compute.mjs)), so that Python and C reach the algorithms while they are JavaScript. Only the exported functions of those public interfaces are callable, by name; classes and constants are listed but not called; nothing is evaluated. A call has a time limit (10 s by default), after which its worker is killed and replaced.

`GET /api/v1/compute` lists every unit with its functions, each with its parameters as written and the comment above it in the source, so that the listing documents itself.

*Rationale.* The domains are pure: no network, no page, no storage, which [check-dependencies.mjs](../development/tools/check-dependencies.mjs) verifies. A call therefore reads and changes nothing beyond its arguments, which is what makes it safe to open to every caller. When a domain's C++ exists, the same path may be answered by it.

*Caveat.* JSON has no NaN: NaN and infinities come back as `null`. A function that needs a typed array takes `{"$typed": "Float64Array", "values": [...]}`. Functions whose arguments are large in-memory objects (a parsed ephemeris) are listed, but are of little use through the bridge.

### 2.3 Plug-ins and Jobs

The servers' plug-ins ([ADR 0029](../development/decisions/0029-plugins-in-python-and-cpp.md)), Python and C++ files under each application's `plugins/`, run through the API as through the Earth Atlas's page: describe one (its parameters, its needs and the sources in the library that serve each), run it, and follow its job. A field it makes is written as a product under `data/cache/plugins/<id>/` and is then a field of the API like any other.

The hosting of plug-ins belongs to no application: it is [python/science/hosting.py](python/science/hosting.py), with the classes in [plugins.py](python/science/plugins.py) and the products in [fields.py](python/science/fields.py), moved here from the Earth Atlas's `py/`. The atlas binds it to its own `plugins/` ([applications/earth-atlas/py/host.py](../applications/earth-atlas/py/host.py)). A plug-in's context reads the library: `context.fields()`, `field(name)`, `list(kind)`, `events(name)` and `observations(name)`.

*Rationale.* A run computes; it writes only its own products. The API therefore opens runs to every caller, and bounds them instead: one at a time (a second answers 409), in a process of its own, stopped at a time limit (600 s by default). Only files under an application's `plugins/` are run, as ADR 0029 requires.

### 2.4 Real Time: Server-Sent Events

`GET /api/v1/stream` and `GET /api/v1/observations/{source}/stream` answer `text/event-stream`: each reading posted is sent to every listener as it arrives (within milliseconds on one computer), as an event `readings` with an `id` and the data `{source, readings}`; a comment line every 15 s keeps the connection alive and tells the client it is. A client that reconnects with `Last-Event-ID` receives the events it missed, of the last thousand.

*Rationale.* Server-Sent Events are plain HTTP, one way, and reconnect by design: a browser's `EventSource` does it alone, and the clients here do it in a few lines. Readings flow one way, from the store to the listeners; WebSockets would add a protocol and a handshake for a direction that is not used.

### 2.5 The Write Policy and Tokens

Reads are open to every origin (CORS `*`). The writes are the posts of readings, `POST /observations/{source}`:

| The service is | Writes are |
|---|---|
| on this computer (127.0.0.1), without a token | open to this computer's programs and pages; refused to a page of another origin (a browser's request carries its `Origin`), so that no web page can write to a service on this computer |
| started with a token (`--token`, or `$SCIENCE_API_TOKEN`) | accepted only with `Authorization: Bearer <token>` |
| on another address without a token | refused: secure by default |
| started with `--read-only` | refused |

Computing is not writing: calls into the domains, samples asked by POST (for many places) and plug-in runs are open to every caller, each bounded in time. `GET /api/v1/` states the policy, and whether the request as made could write.

### 2.6 One Implementation, Two Doors

The service answers through `science.Local` ([python/science/local.py](python/science/local.py)), the reader of the cache. The Python client offers the same object, so that a program on the computer holding the cache reads it without a server:

```python
import science
here = science.Local()                                   # the repository's data/cache; SCIENCE_CACHE moves it
here.field("gfs-forecast").sample("tmp2m", 139.69, 35.69, "now+1d")
```

The two cannot drift apart, and the checks hold them to each other exactly.

## 3. Quick Start

### 3.1 The Service

From the repository root, with the data library's products in the cache (the [Quick Start](../README.md#quick-start) fetches them):

```
python api/service/server.py                       # http://127.0.0.1:8800/api/v1/
python api/service/server.py --host 0.0.0.0 --token <a long random string>      # for other computers
```

An application's server may mount the same API (§6); the Earth Atlas's then answers it at http://127.0.0.1:8790/api/v1/.

### 3.2 Python

```
pip install -e api/python          # or nothing: the examples find the package in the repository
```

```python
import science
gfs = science.Client("http://127.0.0.1:8800").field("gfs-forecast")
print(gfs.sample("tmp2m", lon=139.69, lat=35.69, time="now+1d") - 273.15, "°C in Tokyo tomorrow")
```

`grid()` gives a numpy array when numpy is present, lists otherwise; `stream()` is a generator; `run_plugin()` waits for its job. Examples: [examples/python/](examples/python/).

| Example | Shows |
|---|---|
| [forecast_at_a_place.py](examples/python/forecast_at_a_place.py) | GFS's temperature and wind at a place over the coming days, the wind's direction by the meteorology domain, and tomorrow's temperature from every model |
| [sensor_post.py](examples/python/sensor_post.py), [sensor_stream.py](examples/python/sensor_stream.py) | A sensor that posts readings, and a second process that streams them as they arrive |
| [sun_and_wind.py](examples/python/sun_and_wind.py) | The Sun's place over a point now (astronomy), a distance (geodesy), wind directions (meteorology) |
| [export_geotiff.py](examples/python/export_geotiff.py) | A GeoTIFF of a field over Japan, and the same box at every time as netCDF |
| [run_plugin.py](examples/python/run_plugin.py) | The Earth Atlas's vorticity plug-in run through the API, and its result read as a field |

### 3.3 JavaScript

In Node 18 or later, or a browser, with nothing to install:

```js
import { Client } from './api/js/science.mjs';
const gfs = await new Client('http://127.0.0.1:8800').field('gfs-forecast');
console.log(await gfs.sample('tmp2m', 139.69, 35.69, { time: 'now+1d' }));
```

Grids come as `Float32Array`; `stream()` is an async iterator, through `EventSource` in a browser and a reader of fetch's body in Node. Examples: [forecast.mjs](examples/js/forecast.mjs), [stream.mjs](examples/js/stream.mjs), [compute.mjs](examples/js/compute.mjs).

### 3.4 curl

```
curl "http://127.0.0.1:8800/api/v1/fields/gfs-forecast/sample?var=tmp2m&lon=139.69&lat=35.69&time=now"
```

More in [examples/curl.md](examples/curl.md).

### 3.5 C and C++

*Being written.* A client in C, usable from C++, over the same endpoints, with its examples.

## 4. The Endpoints in Brief

Every path is under `/api/v1`; answers are JSON unless stated; an error is `{"error": "..."}` with the status that says why (400 a bad request, 401 a missing token, 403 a refused write, 404 nothing there, 405 a method the path does not answer, 409 a run already running, 413 too much at once, 422 a function that failed, 501 not written yet, 503 what this computer cannot do, 504 a call past its limit). The full description, with every parameter and every shape, is [openapi.json](openapi.json) (OpenAPI 3.1), also served at `/api/v1/openapi.json`.

| Method and path | Answer |
|---|---|
| `GET /` | the version, the endpoints, the write policy |
| `GET /catalog` | the registry's datasets (id, kind, title, licence, parameters, profiles) and the index of products |
| `GET /fields`, `GET /fields/{name}` | the field products; one's description (grid, times, variables with their quanta), without byte offsets |
| `GET /fields/{name}/sample?var=&lon=a,b&lat=c,d&time=` | values at places and an instant; `POST` with `{var, lon, lat, time}` for many places |
| `GET /fields/{name}/series?var=&lon=&lat=` | the value at a place at every time |
| `GET /fields/{name}/grid?var=&time=&bbox=w,s,e,n&stride=&format=json\|f32` | the grid at a stored time (an index) or an instant; `f32` is little-endian float32, row-major, rows north to south, NaN for none, with the grid in the header `X-Grid` |
| `GET /fields/{name}/export?var=&time=&format=tif\|nc\|csv\|json&bbox=&stride=` | a file, by [data/formats/export.py](../data/formats/export.py); `time=all` for every time in netCDF or JSON |
| `GET /observations`, `GET /observations/{source}?since=&until=&bbox=&limit=&format=` | the sets and sensor sources; one's rows or readings |
| `POST /observations/{source}` | readings `{readings: [...], define?}` from a device: a write (§2.5) |
| `GET /stream`, `GET /observations/{source}/stream` | new readings as they arrive (§2.4) |
| `GET /events`, `/events/{name}?since=&until=&bbox=&limit=&format=geojson` | event sets |
| `GET /tracks`, `/tracks/{name}?format=geojson` | track sets |
| `GET /regions`, `/regions/{name}`, `/regions/{name}/at?lon=&lat=` | sets of regions; the region at each place |
| `GET /verification`, `/verification/{name}` | verifications of forecasts |
| `GET /compute`, `GET /compute/{unit}`, `POST /compute/{unit}/{function}` | the units and their functions; a call with `{args: [...]}` (§2.2) |
| `GET /plugins`, `GET /plugins/describe?app=&file=`, `POST /plugins/run`, `GET /jobs/{id}` | the plug-ins, one's description, a run, its job (§2.3) |

## 5. Sensors

A sensor, or the program that reads it, posts its readings to `POST /api/v1/observations/{source}`. The first post may define the source: its title, its variables with the units the device measures in (the store converts them to SI on arrival), and its stations:

```
curl -X POST http://127.0.0.1:8800/api/v1/observations/garden -H "Content-Type: application/json" -d '{
  "define": {"title": "Garden station", "variables": ["air_temperature:degree_C", "relative_humidity:%"]},
  "readings": [{"station": "g1", "lon": 139.69, "lat": 35.69, "time": "2026-10-23T06:00:00Z",
                "values": {"air_temperature": 17.5, "relative_humidity": 64}}]}'
```

The answer counts the readings accepted and lists those refused, each with its reason: the store checks every reading against its source's variables, the gross-error limits, the place and the time ([observation/sensors/store.py](../observation/sensors/store.py); the store, the ingestion of files and the simulator are described in [observation/sensors/README.md](../observation/sensors/README.md)). A client that streams `GET /api/v1/observations/garden/stream` receives each accepted reading as it arrives, and the service rebuilds the source's product `sensor-garden` at most once a minute, so that the Earth Atlas draws it among the observations.

## 6. Mounting the API in an Application's Server

`Api` ([service/handler.py](service/handler.py)) mounts on any standard-library server: it answers the paths under `/api/` and returns `False` for any other, which the server answers itself.

```python
sys.path.insert(0, str(ROOT / "api" / "service"))
from handler import Api
API = Api(ROOT, allow_writes=WRITES["allowed"], token=WRITES["token"])     # created once the server knows its policy

class Handler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if API.handle(self):
            return
        ...                                                                # the application's own paths
    def do_OPTIONS(self):
        if not API.handle(self):
            self.send_error(404)
```

The mounted API keeps the application's write policy and token, so that one server has one policy.

## 7. Security

- **Bound to this computer by default.** The service binds 127.0.0.1 unless `--host` says otherwise; on another address it refuses writes until it has a token.
- **The token** is compared in constant time; give it on the command line or in `$SCIENCE_API_TOKEN`, choose a long random string, and serve over a TLS proxy when it crosses a network, since HTTP carries it in the clear.
- **No page of another origin can write** to a service without a token; with a token, a page can write only if it holds it.
- **Nothing is evaluated.** The compute bridge calls exported functions of the domains' public interfaces by name; the plug-ins that run are the repository's own files under `plugins/`; the store checks every name it is given, so that none leaves its directory.
- **Everything is bounded:** a body at most 8 MB, at most 100 000 places a sample, 2 000 000 cells a JSON grid, 10 000 readings a post; a call 10 s, a run 600 s and one run at a time; a client that stops reading is let go after 120 s.

*Caveat.* A run and a call cost the serving computer time, and they are open to every caller. A public service should sit behind a proxy that limits the rate of requests.

## 8. Versioning

The version is in the path: `/api/v1`. Within v1, endpoints and fields are added and none is removed or changes its meaning; a change that would break a client makes `/api/v2`, served beside v1 for as long as v1 has users. The service's own version (`GET /` → `version`, 0.1.0) follows the package's, and is a draft's until it is declared stable.

## 9. Modules

| Path | Scope | Status |
|---|---|---|
| [`service/handler.py`](service/handler.py) | `Api`: the endpoints, the write policy, the stream; mountable by any standard-library server | draft |
| [`service/server.py`](service/server.py) | The standalone service: `ThreadingHTTPServer`, port 8800, `--host`, `--port`, `--token`, `--read-only` | draft |
| [`service/compute.mjs`](service/compute.mjs) | The Node worker of the compute bridge: JSON lines in and out | draft |
| [`python/science/`](python/science/) | The package `science`: `Client` ([client.py](python/science/client.py)), `Local` ([local.py](python/science/local.py)), the bridge's pool ([bridge.py](python/science/bridge.py)), field products ([fields.py](python/science/fields.py)), plug-ins ([plugins.py](python/science/plugins.py)) and their hosting ([hosting.py](python/science/hosting.py)); [pyproject.toml](python/pyproject.toml) | draft |
| [`js/science.mjs`](js/science.mjs) | The JavaScript client, for Node and browsers | draft |
| `c/` | The C client, usable from C++ | planned |
| [`openapi.json`](openapi.json) | The description of every endpoint (OpenAPI 3.1) | draft |
| [`examples/`](examples/) | Python, JavaScript and curl | draft |
| [`tests/`](tests/) | The checks (§10) | draft |

**Dependencies.** The service reads the data library's cache and registry (`data/`), the sensors' store (`observation/sensors/`), the export of files (`data/formats/export.py`), the domains' and foundations' public interfaces, and the applications' `plugins/`. Nothing in the foundations or the domains depends on it. As an edge it may be used by works; an application that mounts it uses it as it uses the data library.

## 10. Checks

`node api/tests/run-tests.mjs` (run by `node development/tools/run-tests.mjs` with every suite) starts the service on a free port, with the sensors' store and the plug-ins' products in a directory of its own, and checks:

- samples (every variable, five instants, ten places including the poles and both sides of the antimeridian), directions and series equal to the Earth Atlas's own `FieldProduct` on the same product, within 1e-4 relative (they agree exactly); grids equal to the stored values, as JSON and float32, and between stored times their blend; exports as GeoTIFF and netCDF;
- calls of astronomy, meteorology, geodesy and physis equal to the same functions called directly; refusals of what is not a callable function;
- readings posted and received on a stream within a second, implausible readings refused with their reasons, a reconnecting stream receiving what it missed;
- the write policy with and without a token, read-only, and from another origin; the errors (400, 404, 405, 409, 422);
- the plug-ins: listed, described, run one at a time, their field served;
- the Python client (`check_client.py`) against the service and against `Local`, and the interpolation of `fields.py` against cases known in closed form; the JavaScript client.

What the computer cannot do is skipped and said: no Python with numpy, no products in the cache, no export module.

## 11. Open Questions

- **The C++ behind the bridge.** When a domain has C++, whether the bridge calls it through its Python bindings (nanobind, roadmap Phase 0) or a process of its own.
- **Authorisation beyond one token:** tokens per source, so that a sensor can write only its own readings.
- **Streams of fields:** a stream that announces a product rebuilt (a new forecast), as readings are announced now.
