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

Usage::

    python applications/earth-atlas/serve.py [--host 127.0.0.1] [--port 8790] [--open]

then open http://127.0.0.1:8790/ in a browser. On a server, ``--host 0.0.0.0``
serves other computers too, and ``--sibling URL`` says where the Celestial
Atlas is served: the page links to it as ``../celestial-atlas/``, and the
server sends the browser on from there.
"""

from __future__ import annotations

import argparse
import http.server
import json
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
        url = urllib.parse.urlsplit(self.path)
        if url.path == "/request" or url.path.startswith("/request/"):
            self.on_demand(url)
        elif url.path == "/setup":
            self.reply_json(*self.recipe_request(url))
        elif SIBLING and url.path.startswith("/celestial-atlas/"):
            # The page links to the Celestial Atlas as ../celestial-atlas/, keeping
            # the place and instant in the fragment, which the browser carries over
            # a redirect.
            self.send_response(302)
            self.send_header("Location", SIBLING + self.path[len("/celestial-atlas/"):])
            self.end_headers()
        else:
            super().do_GET()

    def do_POST(self) -> None:
        url = urllib.parse.urlsplit(self.path)
        if url.path == "/setup/run":
            self.reply_json(*self.recipe_request(url))
        else:
            self.send_error(404)

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
        try:
            if url.path == "/setup":
                return recipe_job().view(), 200
            query = dict(urllib.parse.parse_qsl(url.query))
            steps = [s for s in query.get("steps", "").split(",") if s] or None
            return recipe_job().start(steps, again=query.get("again") == "1"), 200
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
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    args = parser.parse_args(argv)
    global SIBLING
    SIBLING = args.sibling.rstrip("/") + "/" if args.sibling else None
    if not (DATA / "products" / "geography-ne50m.json").exists():
        print("The geography product is missing. From the repository root, run:\n"
              "  python data/acquisition/fetch.py natural-earth-50m\n"
              "  python data/products/geography.py\n"
              "and see applications/earth-atlas/README.md for the other data.", file=sys.stderr)
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
