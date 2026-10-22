"""Serve the Celestial Atlas on this computer.

The page is served from this directory; the code it shares, from the
repository's ``domains/``, ``foundations/`` and ``visualization/`` under the
same paths; and the data from the repository's cache (``data/cache/``)
under ``/data/``. Only 127.0.0.1 is bound unless ``--host`` says otherwise.

The atlas's data recipe (``recipe.json``: the datasets fetched and the products
built, in order) is reported to the page and run on its request:

    GET  /setup                          each step and whether it is done; a run in progress
    POST /setup/run[?steps=a,b&again=1]  run the missing steps, or the chosen ones again

The API of the low layers (api/README.md) is mounted at /api/v1/, under this
server's write policy and token: fields, observations and their stream, events,
the domains' functions and the plug-ins, for programs in any language.

Usage::

    python applications/celestial-atlas/serve.py [--host 127.0.0.1] [--port 8770] [--open]

then open http://127.0.0.1:8770/ in a browser. On a server, ``--host 0.0.0.0``
serves other computers too, and ``--sibling URL`` says where the Earth Atlas
is served: the page links to it as ``../earth-atlas/``, and the server sends
the browser on from there.
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
SIBLING = None       # where the Earth Atlas is served, when it is not beside this page
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


def recipe_job():
    global RECIPE
    if RECIPE is None:
        sys.path.insert(0, str(ROOT / "data" / "acquisition"))
        import recipe
        RECIPE = recipe.Job(APP / "recipe.json", ROOT)
    return RECIPE

# Explicit types: on Windows the registry may otherwise map .js to text/plain,
# which browsers refuse for ES modules.
TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".mjs": "text/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json",
    ".bin": "application/octet-stream",
    ".bsp": "application/octet-stream",
    ".svg": "image/svg+xml",
    ".png": "image/png",
}


class Server(http.server.ThreadingHTTPServer):
    # A page asks for its modules, maps and data at once: a deeper queue than
    # the default of five keeps the connections from being refused.
    request_queue_size = 128
    daemon_threads = True


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


API = None           # the API of the low layers (api/service/handler.py), mounted at /api/ by main()


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
        if url.path == "/setup":
            self.reply_json(*self.recipe_request(url))
        # The page links to the Earth Atlas as ../earth-atlas/, keeping the place
        # and instant in the fragment, which the browser carries over a redirect.
        elif SIBLING and url.path.startswith("/earth-atlas/"):
            self.send_response(302)
            self.send_header("Location", SIBLING + self.path[len("/earth-atlas/"):])
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
        else:
            self.send_error(404)

    def do_OPTIONS(self) -> None:
        if API is None or not API.handle(self):           # the API's preflights; nothing else here answers OPTIONS
            self.send_error(404)

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
        # Fixed downloads (the ephemeris) are kept for a day; products, which
        # change when rebuilt, are checked again each time (a 304 if unchanged).
        if self.path.startswith("/data/") and not self.path.startswith("/data/products/"):
            self.send_header("Cache-Control", "public, max-age=86400")
        else:
            self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, format: str, *args) -> None:
        if args and str(args[1]).startswith(("4", "5")):
            super().log_message(format, *args)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Serve the Celestial Atlas locally.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8770)
    parser.add_argument("--sibling", metavar="URL",
                        help="where the Earth Atlas is served when it is not beside this page")
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
    if not (DATA / "jpl-de440s" / "de440s.bsp").exists():
        print("The data cache is empty. From the repository root, run:\n"
              "  python data/acquisition/fetch.py\n  python data/products/star_catalog.py", file=sys.stderr)
    # The API of the low layers at /api/v1/, under this server's write policy and token.
    sys.path.insert(0, str(ROOT / "api" / "service"))
    from handler import Api
    API = Api(ROOT, allow_writes=WRITES["allowed"], token=WRITES["token"],
              read_only_reason=None if WRITES["allowed"] else "this server accepts no writes: start it with --token")
    server = Server((args.host, args.port), Handler)
    url = f"http://127.0.0.1:{args.port}/"
    where = url if args.host == "127.0.0.1" else f"{url} and on {args.host}"
    print(f"Celestial Atlas at {where}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
