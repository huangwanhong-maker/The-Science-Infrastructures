"""Serve the Celestial Atlas on this computer.

The page is served from this directory; the code it shares, from the
repository's ``domains/``, ``foundations/`` and ``visualization/`` under the
same paths; and the data from the repository's cache (``data/cache/``)
under ``/data/``. Only 127.0.0.1 is bound unless ``--host`` says otherwise.

The atlas's data recipe (``recipe.json``: the datasets fetched and the products
built, in order) is reported to the page and run on its request:

    GET  /setup                          each step and whether it is done; a run in progress
    POST /setup/run[?steps=a,b&again=1]  run the missing steps, or the chosen ones again

Usage::

    python applications/celestial-atlas/serve.py [--host 127.0.0.1] [--port 8770] [--open]

then open http://127.0.0.1:8770/ in a browser. On a server, ``--host 0.0.0.0``
serves other computers too, and ``--sibling URL`` says where the Earth Atlas
is served: the page links to it as ``../earth-atlas/``, and the server sends
the browser on from there.
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
SIBLING = None       # where the Earth Atlas is served, when it is not beside this page
RECIPE = None        # the job that runs the atlas's data recipe, created at the first request


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
        if url.path == "/setup":
            self.reply_json(*self.recipe_request(url))
        # The page links to the Earth Atlas as ../earth-atlas/, keeping the place
        # and instant in the fragment, which the browser carries over a redirect.
        elif SIBLING and url.path.startswith("/earth-atlas/"):
            self.send_response(302)
            self.send_header("Location", SIBLING + self.path[len("/earth-atlas/"):])
            self.end_headers()
        else:
            super().do_GET()

    def do_POST(self) -> None:
        url = urllib.parse.urlsplit(self.path)
        if url.path == "/setup/run":
            self.reply_json(*self.recipe_request(url))
        else:
            self.send_error(404)

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
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    args = parser.parse_args(argv)
    global SIBLING
    SIBLING = args.sibling.rstrip("/") + "/" if args.sibling else None
    if not (DATA / "jpl-de440s" / "de440s.bsp").exists():
        print("The data cache is empty. From the repository root, run:\n"
              "  python data/acquisition/fetch.py\n  python data/products/star_catalog.py", file=sys.stderr)
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
