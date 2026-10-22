# 0030. The API: an Edge Through Which Other Programs Reach the Infrastructure

- **Status:** Accepted, provisionally: built at the user's request; the placement agreed by the user on 2026-10-23; open to revision on review. Amends [ADR 0029](0029-plugins-in-python-and-cpp.md).
- **Date:** 2026-10-23
- **Decided by:** the user asked that the low layers (algorithms, data reading, processing, spatial data, observations) be usable independently of any application, through an API, in real time, from Python and C/C++ (and JavaScript), with documents. The user chose a new top-level edge, `api/`, over a module of `data/`. Claude chose the design below.

## 1. Context

Until now, everything the infrastructure computes or holds reached people through two web applications. The data library's products are files in its cache; the domains' algorithms are JavaScript modules imported by pages; the plug-in hosting lived inside the Earth Atlas. A program of the user's own (a model, a notebook, a sensor's firmware, a research code in C++) had no stable way in.

The repository's principle applies: *own the interfaces, borrow the implementations*. The interface is what must stay stable while the implementations behind it change, as the domains' provisional JavaScript will give way to C++ ([ADR 0003](0003-cpp-computational-core.md)).

## 2. Decision

1. **A new edge, `api/`**, beside `data/`, `visualization/` and `observation/`. It reads the data library's cache and registry and the sensors' store, and calls the domains and foundations only through their public interfaces. Nothing in `foundations/`, `domains/`, `data/` or `observation/` depends on it; an application may mount it in its server.
2. **`/api/v1` over HTTP is the owned interface**, described in OpenAPI 3.1 (`api/openapi.json`). Its resources are:
   - fields: description, samples, series, grids, exports in standard formats;
   - observations: of stations, buoys and the user's sensors, read and posted;
   - events, tracks, regions and the scorecard;
   - the domains' functions;
   - plug-ins and their jobs.

   The paths are kept within v1; what answers them may change.
3. **One implementation behind every door.** `science.Local` reads the cache directly; the service answers through it, so that a program on the same computer may use `Local` without a server and get the same answers. Samples follow the page's rule exactly (physis's `grid.js`), and the checks hold them equal to the page's to the last digit.
4. **The compute bridge.** The domains' exported pure functions are called by name through a pool of persistent Node workers, with a time limit per call. Nothing is evaluated; classes and constants are not callable. When a domain gains its C++ core, the same endpoint calls that.
5. **Real time by Server-Sent Events**: a posted reading reaches every listener at once (in about 25 ms here), with heartbeats and resumption by `Last-Event-ID`.
6. **The write policy.** Reads are open to every origin. Writes (posting readings) need the token when one is set, and are refused on a public address without one. Computing (calls into the domains, plug-in runs) is not writing: it is open, bounded in time, and plug-in runs go one at a time, each in its own process. The applications' servers follow the same rule for their data recipes.
7. **Clients** in Python (`api/python/science`, installable, standard library with optional numpy), JavaScript (`api/js/science.mjs`, Node and browsers), and C and C++ (`api/cpp/science.hpp`, header-only, and a C interface in `api/c/`), with runnable examples in each.
8. **The plug-in contract leaves the application.** The Python classes, the reading and writing of field products, and the hosting of server plug-ins move into the package `science` (`api/python/science`); the C++ classes into `api/cpp/plugin.hpp`. The Earth Atlas binds them to its own `plugins/` directory. This amends ADR 0029, whose decisions stand.

## 3. Alternatives Considered

- **`data/api`**, a module of the data library. No new directory; but the API also serves the domains' algorithms and the sensors, which are not the data library's.
- **WebSockets** for real time. Two-way, but a dependency or a hand-written framing on the standard library; Server-Sent Events are one-way, which is what readings need, and every browser and HTTP client reads them.
- **FastAPI or another framework.** Less code, but a dependency for every installation, against the rule that the data library needs only the standard library and numpy.
- **Porting the domains to Python** for the API. A second implementation to keep equal to the first; the bridge keeps one, until C++ replaces it.
- **Plug-in runs in threads.** Simpler, but a thread cannot be stopped when it overruns.

## 4. Consequences

- **Ports:** the standalone service listens on 8800; both atlases serve `/api/v1` beside their pages, under their own write policy and token.
- **The API runs the applications' `plugins/` files**, which sits uneasily with "nothing depends on an application". The files are the plug-ins' home for now; a shared `plugins/` edge may follow if a second application takes plug-ins.
- **Two hosts may run the same server plug-in** (the atlas's `/plugin-host/` and the API's); each runs one at a time within itself, and both write to `data/cache/plugins/<id>/`.
- **The Node workers load the domains once**: a new function appears after the service restarts.
- **Checks:** the API's suite (`api/tests/`) starts its own service and holds the clients equal to `Local`, and `Local` equal to the page.
