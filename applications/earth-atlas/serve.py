"""Serve the Earth Atlas on this computer.

The page is served from this directory; the code it shares, from the
repository's ``domains/``, ``foundations/`` and ``visualization/`` under the
same paths; and the data products from the repository's cache (``data/cache/``)
under ``/data/``. Only 127.0.0.1 is bound unless ``--host`` says otherwise.

Data that cannot be fetched ahead of time, such as the ocean beneath a place
the user has just chosen, are requested on demand through the data library
(data/acquisition/on_demand.py):

    /request?dataset=<id>&<parameter>=<value>...   starts a request, or finds it; returns the job as JSON
    /request/<job id>                                the job's state; when done, the product's path under /data/

The atlas's data recipe (``recipe.json``: the datasets fetched and the products
built, in order) is reported to the page and run on its request:

    GET  /setup                          each step and whether it is done; a run in progress
    POST /setup/run[?steps=a,b&again=1]  run the missing steps, or the chosen ones again

Plug-ins written in Python or C++ run here rather than in the page (py/host.py):

    GET  /plugin-host/describe?file=plugins/x.py   identity, parameters and needs, with the sources that serve each
    POST /plugin-host/run                           {file, chosen, params, context} starts a job; returns {job}
    GET  /plugin-host/job/<id>                      the job's state, progress and results

The API of the low layers (api/README.md) is mounted at /api/v1/, under this
server's write policy and token: fields, observations and their stream, events,
the domains' functions and the plug-ins, for programs in any language.

Usage::

    python applications/earth-atlas/serve.py [--host 127.0.0.1] [--port 8790] [--open]

then open http://127.0.0.1:8790/ in a browser. On a server, ``--host 0.0.0.0``
serves other computers too, and ``--sibling URL`` says where the Celestial
Atlas is served: the page links to it as ``../celestial-atlas/``, and the
server sends the browser on from there.
"""

from __future__ import annotations

import argparse
import email.utils
import gzip
import http.server
import json
import os
import sys
import urllib.parse
import webbrowser
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP.parents[1]
DATA = ROOT / "data" / "cache"
MOUNTS = {
    "/data/": DATA,
    "/domains/": ROOT / "domains",
    "/foundations/": ROOT / "foundations",
    "/visualization/": ROOT / "visualization",
}
JOBS = None          # the data library's jobs on demand, created at the first request
SIBLING = None       # where the Celestial Atlas is served, when it is not beside this page
RECIPE = None        # the job that runs the atlas's data recipe, created at the first request
# Who may start fetches and builds (the data recipe): everyone on this computer;
# on a public address, only the holder of the token, and no one without one.
WRITES = {"allowed": True, "token": None}


def write_refusal(handler) -> str | None:
    """None when this request may write; else the reason it may not."""
    if not WRITES["allowed"]:
        return "this server is read-only: start it with --token to allow fetching with a token"
    token = WRITES["token"]
    if token and handler.headers.get("Authorization", "") != f"Bearer {token}":
        return "fetching on this server needs its token"
    return None


def jobs():
    global JOBS
    if JOBS is None:
        sys.path.insert(0, str(ROOT / "data" / "acquisition"))
        import on_demand
        JOBS = on_demand.Jobs()
    return JOBS


def recipe_job():
    global RECIPE
    if RECIPE is None:
        sys.path.insert(0, str(ROOT / "data" / "acquisition"))
        import recipe
        RECIPE = recipe.Job(APP / "recipe.json", ROOT)
    return RECIPE


API = None           # the API of the low layers (api/service/handler.py), mounted at /api/ by main()
PLUGINS = None       # the host of the plug-ins run here (py/host.py), created at the first request
PLUGIN_JOBS = None   # its jobs


def plugin_host():
    global PLUGINS
    if PLUGINS is None:
        sys.path.insert(0, str(APP / "py"))
        import host
        PLUGINS = host
    return PLUGINS

# Explicit types: on Windows the registry may otherwise map .js to text/plain,
# which browsers refuse for ES modules.
TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".mjs": "text/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json",
    ".bin": "application/octet-stream",
    ".svg": "image/svg+xml",
    ".png": "image/png",
}


# Large products travel compressed when the browser accepts it: a regions raster
# of 13 MB becomes 87 kB. Each file is compressed once per version (its size and
# time of change) and kept in memory, the oldest dropped beyond GZIP_LIMIT bytes.
GZIPPED: dict = {}
GZIP_LIMIT = 512 * 2 ** 20


def gzipped(path: Path) -> bytes | None:
    """The file compressed, or None when that saves too little to be worth it."""
    st = path.stat()
    hit = GZIPPED.get(path)
    if hit and hit[0] == (st.st_mtime_ns, st.st_size):
        return hit[1]
    data = path.read_bytes()
    packed = gzip.compress(data, 6)
    if len(packed) > 0.8 * len(data):
        packed = None
    GZIPPED[path] = ((st.st_mtime_ns, st.st_size), packed)
    while sum(len(v[1] or b"") for v in GZIPPED.values()) > GZIP_LIMIT and len(GZIPPED) > 1:
        GZIPPED.pop(next(iter(GZIPPED)))
    return packed


class Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map, **TYPES}

    def translate_path(self, path: str) -> str:
        path = urllib.parse.unquote(path.split("?", 1)[0].split("#", 1)[0])
        prefix = next((p for p in MOUNTS if path.startswith(p)), None)
        base, rel = (MOUNTS[prefix], path[len(prefix):]) if prefix else (APP, path.lstrip("/"))
        target = (base / rel).resolve()
        if target != base and base not in target.parents:
            return str(base / "__outside__")
        return str(target)

    def do_GET(self) -> None:
        if API is not None and API.handle(self):        # /api/v1/...: the API of the low layers
            return
        url = urllib.parse.urlsplit(self.path)
        if url.path == "/request" or url.path.startswith("/request/"):
            self.on_demand(url)
        elif url.path == "/setup":
            self.reply_json(*self.recipe_request(url))
        elif url.path.startswith("/plugin-host/"):
            self.reply_json(*self.plugin_request(url))
        elif SIBLING and url.path.startswith("/celestial-atlas/"):
            # The page links to the Celestial Atlas as ../celestial-atlas/, keeping
            # the place and instant in the fragment, which the browser carries over
            # a redirect.
            self.send_response(302)
            self.send_header("Location", SIBLING + self.path[len("/celestial-atlas/"):])
            self.end_headers()
        elif not self.send_compressed():
            super().do_GET()

    def send_compressed(self) -> bool:
        """Sends a large .json or .bin product gzipped, if the browser accepts it;
        answers 304 when it has the version already. False when not done here."""
        if "gzip" not in self.headers.get("Accept-Encoding", ""):
            return False
        path = Path(self.translate_path(self.path))
        if path.suffix not in (".json", ".bin") or not path.is_file() or path.stat().st_size < 256 * 1024:
            return False
        st = path.stat()
        since = self.headers.get("If-Modified-Since")
        if since:
            try:
                if email.utils.parsedate_to_datetime(since).timestamp() >= int(st.st_mtime):
                    self.send_response(304)
                    self.end_headers()
                    return True
            except (TypeError, ValueError, OverflowError):
                pass
        packed = gzipped(path)
        if packed is None:
            return False
        self.send_response(200)
        self.send_header("Content-Type", self.guess_type(str(path)))
        self.send_header("Content-Encoding", "gzip")
        self.send_header("Content-Length", str(len(packed)))
        self.send_header("Last-Modified", self.date_time_string(st.st_mtime))
        self.send_header("Vary", "Accept-Encoding")
        self.end_headers()
        self.wfile.write(packed)
        return True

    def do_POST(self) -> None:
        if API is not None and API.handle(self):
            return
        url = urllib.parse.urlsplit(self.path)
        if url.path == "/setup/run":
            self.reply_json(*self.recipe_request(url))
        elif url.path == "/plugin-host/run":
            self.reply_json(*self.plugin_request(url, self.body()))
        else:
            self.send_error(404)

    def do_OPTIONS(self) -> None:
        if API is None or not API.handle(self):           # the API's preflights; nothing else here answers OPTIONS
            self.send_error(404)

    def body(self) -> dict:
        length = int(self.headers.get("Content-Length") or 0)
        return json.loads(self.rfile.read(length).decode("utf-8")) if length else {}

    def plugin_request(self, url, body: dict | None = None) -> tuple[dict, int]:
        """A plug-in run here: its description (GET /plugin-host/describe?file=), a run
        (POST /plugin-host/run) and a job's state (GET /plugin-host/job/<id>)."""
        global PLUGIN_JOBS
        try:
            host = plugin_host()
            if url.path == "/plugin-host/describe":
                return host.describe(dict(urllib.parse.parse_qsl(url.query)).get("file", "")), 200
            if PLUGIN_JOBS is None:
                PLUGIN_JOBS = host.Jobs()
            if url.path == "/plugin-host/run":
                b = body or {}
                return PLUGIN_JOBS.start(b.get("file", ""), b.get("chosen") or {}, b.get("params") or {}, b.get("context") or {}), 200
            if url.path.startswith("/plugin-host/job/"):
                return PLUGIN_JOBS.status(url.path[len("/plugin-host/job/"):]), 200
            return {"error": "unknown"}, 404
        except ImportError as error:
            return {"error": f"plug-ins run here need numpy: {error}"}, 503
        except ValueError as error:          # a file outside plugins/, a bad request
            return {"error": str(error)}, 400
        except RuntimeError as error:        # no compiler, a compilation error, a failed program
            return {"error": str(error)}, 503
        except Exception as error:           # anything else a script did at import or describe
            return {"error": f"{type(error).__name__}: {error}"}, 500

    def on_demand(self, url) -> None:
        try:
            if url.path == "/request":
                query = dict(urllib.parse.parse_qsl(url.query))
                dataset = query.pop("dataset", "")
                body, code = jobs().request(dataset, query), 200
            else:
                body, code = jobs().status(url.path[len("/request/"):]), 200
        except ImportError as error:
            body, code = {"state": "failed", "message": f"requests on demand need numpy: {error}"}, 503
        except ValueError as error:          # RequestError: a refused request
            body, code = {"state": "failed", "message": str(error)}, 400
        self.reply_json(body, code)

    def recipe_request(self, url) -> tuple[dict, int]:
        """The state of the data recipe (GET /setup), or a run of its missing
        steps, or of chosen steps again (POST /setup/run?steps=a,b&again=1)."""
        refusal = write_refusal(self)
        policy = {"writable": refusal is None, "writes": "token" if WRITES["token"] else ("open" if WRITES["allowed"] else "closed"),
                  "refusal": refusal}
        try:
            if url.path == "/setup":
                return {**recipe_job().view(), **policy}, 200
            if refusal:
                return {"message": refusal, **policy}, 403
            query = dict(urllib.parse.parse_qsl(url.query))
            steps = [s for s in query.get("steps", "").split(",") if s] or None
            return {**recipe_job().start(steps, again=query.get("again") == "1"), **policy}, 200
        except RuntimeError as error:          # a run is already in progress
            return {"message": str(error)}, 409

    def reply_json(self, body: dict, code: int) -> None:
        data = json.dumps(body, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def end_headers(self) -> None:
        # Products are rebuilt in place when data are fetched again, so they are
        # revalidated rather than cached for long.
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, format: str, *args) -> None:
        if args and str(args[1]).startswith(("4", "5")):
            super().log_message(format, *args)


class Server(http.server.ThreadingHTTPServer):
    # The page asks for many products at once; the default backlog of 5
    # connections makes Windows refuse the rest.
    request_queue_size = 128
    daemon_threads = True


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Serve the Earth Atlas locally.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8790)
    parser.add_argument("--sibling", metavar="URL",
                        help="where the Celestial Atlas is served when it is not beside this page")
    parser.add_argument("--token", default=os.environ.get("SCIENCE_TOKEN"),
                        help="the token that allows fetching through this server (default: $SCIENCE_TOKEN)")
    parser.add_argument("--read-only", action="store_true", help="refuse every fetch and build requested through the server")
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    args = parser.parse_args(argv)
    global SIBLING, API
    SIBLING = args.sibling.rstrip("/") + "/" if args.sibling else None
    WRITES["token"] = args.token or None
    # On a public address, writes need a token: secure by default.
    WRITES["allowed"] = not args.read_only and (args.host in ("127.0.0.1", "localhost", "::1") or bool(args.token))
    if not WRITES["allowed"] and not args.read_only:
        print("Serving other computers without --token: fetching through the page is refused.", file=sys.stderr)
    if not (DATA / "products" / "geography-ne50m.json").exists():
        print("The geography product is missing. From the repository root, run:\n"
              "  python data/acquisition/fetch.py natural-earth-50m\n"
              "  python data/products/geography.py\n"
              "and see applications/earth-atlas/README.md for the other data.", file=sys.stderr)
    # The API of the low layers at /api/v1/, under this server's write policy and token.
    sys.path.insert(0, str(ROOT / "api" / "service"))
    from handler import Api
    API = Api(ROOT, allow_writes=WRITES["allowed"], token=WRITES["token"],
              read_only_reason=None if WRITES["allowed"] else "this server accepts no writes: start it with --token")
    server = Server((args.host, args.port), Handler)
    url = f"http://127.0.0.1:{args.port}/"
    where = url if args.host == "127.0.0.1" else f"{url} and on {args.host}"
    print(f"Earth Atlas at {where}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
