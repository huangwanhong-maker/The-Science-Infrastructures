"""Serve the Celestial Atlas on this computer.

The page is served from this directory, and the data from the repository's
cache (``data/cache/``) under ``/data/``. Only 127.0.0.1 is bound.

Usage::

    python applications/astronomy/celestial-atlas/serve.py [--port 8770] [--open]

then open http://127.0.0.1:8770/ in a browser.
"""

from __future__ import annotations

import argparse
import http.server
import sys
import urllib.parse
import webbrowser
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP.parents[2]
DATA = ROOT / "data" / "cache"

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
        base, rel = (DATA, path[len("/data/"):]) if path.startswith("/data/") else (APP, path.lstrip("/"))
        target = (base / rel).resolve()
        if target != base and base not in target.parents:
            return str(base / "__outside__")
        return str(target)

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
    parser.add_argument("--port", type=int, default=8770)
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    args = parser.parse_args(argv)
    if not (DATA / "jpl-de440s" / "de440s.bsp").exists():
        print("The data cache is empty. From the repository root, run:\n"
              "  python data/acquisition/fetch.py\n  python data/products/star_catalog.py", file=sys.stderr)
    server = Server(("127.0.0.1", args.port), Handler)
    url = f"http://127.0.0.1:{args.port}/"
    print(f"Celestial Atlas at {url}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
