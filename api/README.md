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

### 2.7 Plug-ins over Networks of Places

Some applications draw flows between places rather than fields on a grid: the Economy Atlas draws the world's trade between economies ([ADR 0046](../development/decisions/0046-the-economy-atlas-data-and-plugins.md)). Their plug-ins are given a **network**, places (an id, a name, a latitude and longitude, a group) and flows between them (a value, optionally a kind), and answer **results** the application draws: `values` (a number for each place), `arcs` (flows between places), and the sketches, tables and notes of §2.3. The classes are [python/science/networks.py](python/science/networks.py) (`Network`, `NetworkAnalysis`, `values`, `arcs`, `param`) and [cpp/network_plugin.hpp](cpp/network_plugin.hpp) (`science::network::Analysis`, `SCIENCE_NETWORK_PLUGIN`); an application binds a `networks.Host` to its own `plugins/`, with the include directories its C++ plug-ins may use.

Both languages run as processes behind one file contract: the host writes `input.json` (the network and the parameters, defaulted and typed) into a directory of the run's own, runs the plug-in with `--run <dir>`, follows its progress lines and reads `output.json`. A script runs as `python -m science.networks run <script> <dir>`; a C++ plug-in is compiled as hosting.py compiles its own (natively, or with g++ in WSL), again when it or a header it includes changes.

```python
def create(api):
    class Balance(api.NetworkAnalysis):
        id, title = "balance", "Net inflow"
        def run(self, network, params, context):
            return [api.values("Net inflow", network.balance(), scale="diverging")]
    return Balance()
```

*Rationale.* A plug-in over a network is written by the person using the application, and runs whatever it contains; a process of its own bounds it in time (600 s by default, or as the host is told), lets it be cancelled, and keeps a failure from reaching the server. A network of a few hundred places and tens of thousands of flows passes as JSON in a fraction of a second, so the contract needs no binary files, unlike the fields of §2.3.

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

A client in C++17, header-only ([cpp/science.hpp](cpp/science.hpp), with its JSON in [cpp/science/json.hpp](cpp/science/json.hpp) and its HTTP in [cpp/science/http.hpp](cpp/science/http.hpp)), and a C interface over it ([c/science.h](c/science.h), [c/science.cpp](c/science.cpp)), on the standard library and the system's sockets alone: so that a model, a research code or a sensor's program on a small Linux computer reads fields and observations, posts readings, follows the stream and calls the domains' functions. The C++ has the surface of the Python client in C++'s spelling; the C has the calls a program in C, or any language that calls C, most needs.

**Building.** The C++ needs no build: include the header.

```
g++ -std=c++17 -I api/cpp mine.cpp -o mine                    # add -pthread for threads of its own; -lws2_32 on Windows
make -C api                                                   # api/build/: libscience.a, libscience.so, the examples
cc -std=c99 -I api/c mine.c api/build/libscience.a -lstdc++ -lm -o mine
```

The [Makefile](Makefile) is GNU make's (`CXX=clang++ CC=clang`, `CXXSTD=c++20` and `BUILD=` choose otherwise); a C program links the library with the C++ runtime, as above, or links with `c++`.

```cpp
#include "science.hpp"
science::Client api("http://127.0.0.1:8800");                     // or $SCIENCE_API_URL; the token from $SCIENCE_API_TOKEN
science::Field gfs = api.field("gfs-forecast");
double kelvin = gfs.sample("tmp2m", 139.69, 35.69, "now+1d");       // NaN where there is no value
science::Grid g = gfs.grid("tmp2m", "0", science::Box{120, 20, 150, 50});    // float32, rows north to south
science::Json sun = api.call("astronomy", "sunPosition", 2460000.5);
api.stream("garden", [](const science::Event& e) { std::puts(e.readings.dump().c_str()); return true; });
```

| C++ | What it is |
|---|---|
| `Client(url, token, timeout)` | a service: `info()`, `catalog()`, `fields()`, `field(name)`, `observations(source, since, until, bbox, limit)`, `post_observations(source, readings, define)`, `stream(source, on_event, options)`, `events(name, ...)`, `tracks`, `regions`, `verification`, `units`, `compute(unit, function, args)` with the arguments as a JSON array, `call(unit, function, args...)` with them as they are, `plugins()`, `describe_plugin`, `job`, `run_plugin(app, file, chosen, params, wait)`; `get`, `post` and `request` for any path |
| `Field` | a field: `sample(var, lons, lats, time)` → `std::vector<double>`, or one place → `double`; `series(var, lon, lat)` → `Series {times, values, units}`; `grid(var, time, bbox, stride)` → `Grid {lon0, dlon, nlon, lat0, dlat, nlat, values}`, read as float32 (`format=f32`, described by `X-Grid`); `download(...)` and `export_to(path, var, time, format)` |
| `Json` | the value type: `Json::parse`, `dump()`, `j["key"]`, `j[i]`, `as_number()` (null as NaN), `as_string()`, `as_array()`, `Json::object({...})`, `Json::array({...})` |
| `ApiError` | an error: `status` (0 when no answer came), `message`, `url`, and for a failed run its `job` |
| `Reading`, `Event`, `StreamOptions`, `Box`, `utc_now()`, `unix_seconds(iso)` | a reading to post; an event of the stream; how a stream reconnects; a box (west, south, east, north); instants |

| C | What it does |
|---|---|
| `science_connect(url, token, timeout)`, `science_close` | a client; no request is made until the first call, so it may be made before the network is up |
| `science_sample`, `science_series`, `science_grid_read`, `science_grid_free` | values into the caller's buffer; a series (times as Unix seconds) and a grid allocated for the caller |
| `science_define`, `science_post_reading`, `science_post` | writes: a source defined, one reading, or a body of JSON text |
| `science_stream(client, source, last_id, callback, user)` | the stream: the callback receives each event's id, source and readings as JSON text, and NULLs at each heartbeat; it returns 0 to stop |
| `science_compute(client, unit, function, args, &result)` | a call, with JSON text in and out |
| `science_get`, `science_info` | any GET, as JSON text |
| `science_json_number`, `science_json_numbers`, `science_json_text`, `science_json_length` | reading an answer by a path of keys and indices (`"0.values.air_temperature"`), so that C needs no JSON library |
| `science_last_error()`, `science_last_status()`, `science_free` | every call returns `SCIENCE_OK` or a code (argument, unreachable, HTTP, answer, refused, failed); the message and the HTTP status are kept per thread, as `errno` is |

Examples, each built by `make -C api` as `api/build/cpp-<name>` and `api/build/c-<name>`:

| C++ | C | Shows |
|---|---|---|
| [forecast.cpp](examples/cpp/forecast.cpp) | [forecast.c](examples/c/forecast.c) | GFS's temperature (and in C++ the wind, its direction by the meteorology domain) at a place at every time |
| [sensor.cpp](examples/cpp/sensor.cpp) | [sensor.c](examples/c/sensor.c) | a sensor posting readings in a loop; the C one is the template for a device running Linux |
| [stream.cpp](examples/cpp/stream.cpp) | [stream.c](examples/c/stream.c) | the stream followed, each reading printed as it arrives |
| [sun.cpp](examples/cpp/sun.cpp) | [sun.c](examples/c/sun.c) | where the Sun stands now: `julianDate`, `sunPosition`, `subsolarPoint`, `sunAltAz` and `daylight` of the astronomy domain, the result of one call an argument of the next |
| [grid_mean.cpp](examples/cpp/grid_mean.cpp) | [grid_mean.c](examples/c/grid_mean.c) | a grid as float32 and its mean, over the points and over the area (weights cos φ) |

**A device that posts readings.** A Raspberry Pi, or any small computer running Linux, builds the library and the examples itself (`sudo apt install build-essential`, then `make -C api`), and runs [sensor.c](examples/c/sensor.c) or a program grown from it: `c-sensor --url http://<the service>:8800 --token <its token> --w1 /sys/bus/w1/devices/28-…/temperature` reads a DS18B20 probe on the 1-Wire bus every minute and posts it. The service it posts to is started for other computers with `--host 0.0.0.0 --token <a long random string>` (§2.5). The program defines its source once, waits for its clock to be set (the store refuses times before 1990 or a quarter of an hour ahead), keeps up to 1024 readings while the service is away and sends them, with their own times, when it returns; a reading the store refuses (`SCIENCE_ERROR_REFUSED`) is reported with its reason and not sent again, and a wrong token stops it. Its comment shows the systemd unit that starts it with the device.

*Limits.*

- **Plain HTTP only.** The client brings no TLS library. Across a network that is not trusted, put the service behind a TLS proxy and run the program where the hop to the proxy is trusted (the proxy's computer, or a tunnel), or use the Python client; a token over plain HTTP crosses the network in the clear.
- **One request per connection** (`Connection: close`). A client holds no connection between calls, and so outlives a network that comes and goes; a program that makes thousands of small requests a second would want keep-alive, which is not written. Many places are sampled in one request: more than 500 go by `POST`, up to 100 000.
- **The stream's reconnection.** After a dropped connection, or 45 s without a byte (heartbeats come every 15 s), the stream connects again after 3 s (or the delay the service asks for, `retry:`), sending the last event's id as `Last-Event-ID`; the service replays what was missed of its last thousand events, and what is older is lost. An answer of an error (404, 503) is thrown, not retried.
- **Time limits** bound connecting and each wait for bytes (60 s by default), not a whole answer, which for a large grid may be long. Redirects are not followed.
- **Numbers are doubles**: integers are exact to 2^53, which holds every count and event id of the API; NaN and the infinities go as `null`. JSON nests at most 512 levels.
- **Threads.** A `Client` holds only its address, token and time limit, so one may serve several threads at once; the C error state is per thread.
- **Where it is checked.** Linux, with g++ 11.5 (in WSL), as C++17, C++20 and C99. The Winsock code for Windows is written but not compiled here (there is no compiler on this computer's Windows side), nor is MSVC, Clang or macOS tried; there is no CMake file, since CMake is not installed here.

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
| [`python/science/`](python/science/) | The package `science`: `Client` ([client.py](python/science/client.py)), `Local` ([local.py](python/science/local.py)), the bridge's pool ([bridge.py](python/science/bridge.py)), field products ([fields.py](python/science/fields.py)), plug-ins ([plugins.py](python/science/plugins.py)) and their hosting ([hosting.py](python/science/hosting.py)), plug-ins over networks of places and their host ([networks.py](python/science/networks.py)); [pyproject.toml](python/pyproject.toml) | draft |
| [`js/science.mjs`](js/science.mjs) | The JavaScript client, for Node and browsers | draft |
| [`cpp/science.hpp`](cpp/science.hpp), [`cpp/science/`](cpp/science/) | The C++ client, header-only: `Client`, `Field`, the stream; its JSON ([json.hpp](cpp/science/json.hpp)) and its HTTP over the system's sockets ([http.hpp](cpp/science/http.hpp)) | draft |
| [`cpp/plugin.hpp`](cpp/plugin.hpp), [`cpp/network_plugin.hpp`](cpp/network_plugin.hpp) | The classes of C++ plug-ins: over fields (§2.3), and over networks of places (§2.7) | draft |
| [`c/`](c/) | The C client ([science.h](c/science.h)), implemented over the C++ ([science.cpp](c/science.cpp)) | draft |
| [`Makefile`](Makefile) | Builds `libscience.a`, `libscience.so`, the C and C++ examples and the check programs into `build/` | draft |
| [`openapi.json`](openapi.json) | The description of every endpoint (OpenAPI 3.1) | draft |
| [`examples/`](examples/) | Python, JavaScript, C, C++ and curl | draft |
| [`tests/`](tests/) | The checks (§10) | draft |

**Dependencies.** The service reads the data library's cache and registry (`data/`), the sensors' store (`observation/sensors/`), the export of files (`data/formats/export.py`), the domains' and foundations' public interfaces, and the applications' `plugins/`. Nothing in the foundations or the domains depends on it. As an edge it may be used by works; an application that mounts it uses it as it uses the data library.

## 10. Checks

`node api/tests/run-tests.mjs` (run by `node development/tools/run-tests.mjs` with every suite) starts the service on a free port, with the sensors' store and the plug-ins' products in a directory of its own, and checks:

- samples (every variable, five instants, ten places including the poles and both sides of the antimeridian), directions and series equal to the Earth Atlas's own `FieldProduct` on the same product, within 1e-4 relative (they agree exactly); grids equal to the stored values, as JSON and float32, and between stored times their blend; exports as GeoTIFF and netCDF;
- calls of astronomy, meteorology, geodesy and physis equal to the same functions called directly; refusals of what is not a callable function;
- readings posted and received on a stream within a second, implausible readings refused with their reasons, a reconnecting stream receiving what it missed;
- the write policy with and without a token, read-only, and from another origin; the errors (400, 404, 405, 409, 422);
- the plug-ins: listed, described, run one at a time, their field served;
- the plug-ins over networks of places (`check_networks.py`): a network's matrix, balances and distances against closed forms; parameters defaulted, typed and refused out of bounds; a plug-in described and run through a host in Python and, where a compiler is found, the same in C++ with the same results; a failure reported with its reason, a file outside `plugins/` refused, a run stopped at its time limit and when cancelled;
- the Python client (`check_client.py`) against the service and against `Local`, and the interpolation of `fields.py` against cases known in closed form; the JavaScript client;
- the C and C++ clients, where a compiler is found ([check_cpp.mjs](tests/check_cpp.mjs)): the Makefile builds everything without a warning, and the check programs ([check_cpp.cpp](tests/check_cpp.cpp), [check_c.c](tests/check_c.c)) hold their samples, series, grids and calls to the Python client's for the same requests ([check_cpp.py](tests/check_cpp.py)), exactly; receive a posted reading on the stream within a second, and after reconnecting what they missed; run a plug-in to its field (C++); and check, offline, their JSON (exponents, the limits of doubles, escapes, surrogate pairs, nesting, malformed texts) and their HTTP (chunked bodies, framing, the events of a stream, instants). On Windows the compiler is WSL's: the programs run in WSL, and reach the service in WSL's mirrored network, or else a service started inside WSL with WSL's Python (`SCIENCE_CHECK_CPP_SERVICE=wsl` takes that way always).

What the computer cannot do is skipped and said: no Python with numpy, no products in the cache, no export module, no C and C++ compilers.

## 11. Open Questions

- **The C++ behind the bridge.** When a domain has C++, whether the bridge calls it through its Python bindings (nanobind, roadmap Phase 0) or a process of its own.
- **Authorisation beyond one token:** tokens per source, so that a sensor can write only its own readings.
- **Streams of fields:** a stream that announces a product rebuilt (a new forecast), as readings are announced now.
